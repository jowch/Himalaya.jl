using Test, HTTP, JSON3, SQLite, DBInterface, Tables, SHA
using HimalayaUI

# ─────────────────────────────────────────────────────────────────────────────
# Phase 2 Task 2.1 — comparisons.jl helper coverage.
#
# Each helper has its own @testset; the route tests live in
# `test_routes_comparisons.jl`. Helpers that touch the dispatcher rely on
# the apply_event! pre-mint pattern from test_events.jl (`_premint_comparison!`).
# ─────────────────────────────────────────────────────────────────────────────

# ── Fixture helpers ─────────────────────────────────────────────────────────

# Pre-mint a comparisons row at a known id so the dispatcher can run UPDATEs
# (matches the route handler's two-step "INSERT placeholder, then dispatcher
# fills in" pattern). Mirrors `_premint_comparison!` in test_events.jl.
# Uses NULL placeholders to mirror the post-#67 route — the dispatcher's
# `COALESCE(col, ?)` then stamps real values on first fold.
function _premint_cmp!(db, id::Int)
    DBInterface.execute(db,
        """INSERT INTO comparisons (id, title, content_hash, created_at, updated_at)
           VALUES (?, NULL, NULL, NULL, NULL)""", [id])
    nothing
end

function _member_payload(; id=nothing, exposure_id=nothing, display_order::Int=0,
                          band_height::Float64=1.0, y_offset::Float64=0.0,
                          normalization::String="none",
                          color_override=nothing, label_override=nothing,
                          q_window_min=nothing, q_window_max=nothing,
                          peak_display=nothing,
                          snapshot=Dict(:effective_peaks => [],
                                        :confirmed_index => nothing,
                                        :analysis_inputs_hash => "sha256:zero"))
    Dict{Symbol,Any}(
        :id             => id,
        :exposure_id    => exposure_id,
        :display_order  => display_order,
        :band_height    => band_height,
        :y_offset       => y_offset,
        :normalization  => normalization,
        :color_override => color_override,
        :label_override => label_override,
        :q_window_min   => q_window_min,
        :q_window_max   => q_window_max,
        :peak_display   => peak_display,
        :snapshot       => snapshot,
    )
end

function _setup_db()
    tmp = mktempdir()
    db = open_prepared_clone(tmp)
    HimalayaUI.bind_db!(db)
    (db = db, tmp = tmp)
end

# Build an analyzed exposure. Returns (db, experiment_id, sample_id, exposure_id, analysis_dir).
function _setup_analyzed_exposure(tmp::String;
                                   datfile::String = "example_tot.dat",
                                   filename::String = "example_tot")
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", datfile),
       joinpath(analysis_dir, datfile))
    db = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename=filename)
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
    (db = db, experiment_id = exp_id, sample_id = s_id,
     exposure_id = e_id, analysis_dir = analysis_dir)
end

@testset "comparisons.jl helpers" begin

    # I3.6 (#177): `current_content_hash` + `is_author` were Compare-route-only
    # helpers, deleted with `routes_comparisons.jl`. Their tests went with them.

    @testset "compute_member_snapshot: empty effective peaks + no confirmed index" begin
        mktempdir() do tmp
            ctx = _setup_db()
            try
                exp_id = HimalayaUI.create_experiment!(ctx.db; path="/x",
                    data_dir="/x", analysis_dir="/x")
                s_id = HimalayaUI.create_sample!(ctx.db; experiment_id=exp_id)
                e_id = HimalayaUI.create_exposure!(ctx.db; sample_id=s_id)
                snap = HimalayaUI.compute_member_snapshot(ctx.db, e_id)
                @test snap[:effective_peaks] == []
                @test snap[:confirmed_index] === nothing
                @test snap[:analysis_inputs_hash] === nothing
            finally
                close(ctx.db)
            end
        end
    end

    @testset "compute_member_snapshot: auto + manual peaks shape" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            # Add a manual peak so the snapshot has at least one source="manual"
            # entry with intensity = nothing.
            DBInterface.execute(ctx.db,
                "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', ?)",
                [ctx.exposure_id, 0.5])
            snap = HimalayaUI.compute_member_snapshot(ctx.db, ctx.exposure_id)
            @test !isempty(snap[:effective_peaks])
            sources = Set(p[:source] for p in snap[:effective_peaks])
            @test "auto"   in sources
            @test "manual" in sources
            for p in snap[:effective_peaks]
                @test haskey(p, :id)
                @test haskey(p, :q)
                @test haskey(p, :intensity)
                @test haskey(p, :sharpness)
                @test haskey(p, :source)
                if p[:source] == "manual"
                    @test p[:intensity] === nothing
                    # Snapshot contract: `sharpness` is never `nothing` so
                    # the JS-side `MemberSnapshotPeak.sharpness: number`
                    # (non-null) holds. Manual peaks default to 0.0 — the
                    # client coerces null to 0 too (`p.sharpness ?? 0`),
                    # so server + client hash the same canonical bytes.
                    @test p[:sharpness] === 0.0
                else
                    @test p[:intensity] isa Real
                    @test p[:sharpness] isa Real
                end
            end
            @test snap[:analysis_inputs_hash] isa AbstractString
        end
    end

    # Regression (Fix 4 / spec auditor): a NULL-sharpness auto peak must
    # surface in the snapshot as `sharpness = 0.0`, not `nothing`. The JS
    # `MemberSnapshotPeak.sharpness: number` type is non-nullable, and the
    # client `computeMemberSnapshot` already coerces null → 0 (snapshot.ts).
    # Without the server-side coercion, a re-hash of a GET-fetched snapshot
    # diverges from the locally-computed hash and `content_hash` parity
    # breaks. The client-vs-server contract for hash inputs lives in
    # contentHash.test.ts (cross-language fixture parity).
    @testset "compute_member_snapshot: NULL auto sharpness coerces to 0.0" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            # Force one auto peak's sharpness to NULL — simulates a legacy
            # row from before sharpness was always populated.
            target = first(Tables.rowtable(DBInterface.execute(ctx.db,
                "SELECT id FROM auto_peaks WHERE exposure_id = ? LIMIT 1",
                [ctx.exposure_id])))
            DBInterface.execute(ctx.db,
                "UPDATE auto_peaks SET sharpness = NULL WHERE id = ?",
                [target.id])
            snap = HimalayaUI.compute_member_snapshot(ctx.db, ctx.exposure_id)
            patched = first(p for p in snap[:effective_peaks]
                            if p[:id] == target.id)
            @test patched[:source] == "auto"
            @test patched[:sharpness] === 0.0
            # And every entry honours the non-null contract — broad sweep.
            for p in snap[:effective_peaks]
                @test p[:sharpness] isa Real
                @test !(p[:sharpness] === nothing)
            end
        end
    end

    # Regression: pin compute_member_snapshot's parallel SQL implementation
    # against the canonical effective_peaks helper from pipeline.jl. Both
    # produce the same q-tolerance union of (auto_peaks − exclude curations
    # ∪ add curations); this test fails loudly if either side drifts.
    @testset "compute_member_snapshot agrees with effective_peaks (auto + exclude + add)" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            # Pick one auto peak to exclude and add a manual peak at a
            # distinct q so the union exercises both curation kinds.
            auto_q = first(Tables.rowtable(DBInterface.execute(ctx.db,
                "SELECT q FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 1",
                [ctx.exposure_id]))).q
            DBInterface.execute(ctx.db,
                "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', ?)",
                [ctx.exposure_id, auto_q])
            DBInterface.execute(ctx.db,
                "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', ?)",
                [ctx.exposure_id, 0.5])

            # Canonical set via the pipeline.jl helper: load the trace once,
            # then call effective_peaks(db, eid, q, I).
            dat_path = joinpath(ctx.analysis_dir, "example_tot.dat")
            q, I, _ = HimalayaUI.load_dat(dat_path)
            canonical = HimalayaUI.effective_peaks(ctx.db, ctx.exposure_id, q, I)

            # Snapshot set via the parallel SQL.
            snap = HimalayaUI.compute_member_snapshot(ctx.db, ctx.exposure_id)
            snap_ids = sort([p[:id] for p in snap[:effective_peaks]])
            snap_qs  = sort([p[:q]  for p in snap[:effective_peaks]])

            @test snap_ids == sort(canonical.peak_id)
            @test snap_qs  ≈ sort(canonical.q)
            # Sanity: the exclude actually dropped one peak, and the add
            # actually contributed one peak, so the diff is nontrivial.
            @test 0.5 in snap_qs
            @test !(auto_q in snap_qs)
        end
    end

    @testset "compute_member_snapshot: confirmed_index sources from the assignment" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp; datfile="cubic_tot.dat",
                                            filename="cubic_tot")
            # D-10: confirmed_index is the highest-scored member of the durable
            # assignment, NOT the legacy active custom group, and there is NO R²
            # gate. A fresh analyze NO LONGER seeds the assignment (auto-grouping
            # is not a durable concept), so an analyzed-but-uncurated exposure
            # has a NULL confirmed_index until a human checks a candidate.
            @test HimalayaUI.compute_member_snapshot(
                ctx.db, ctx.exposure_id)[:confirmed_index] === nothing

            # Check the top-scored auto candidate explicitly → it becomes the
            # confirmed_index (highest-scored member, no R² gate).
            top = first(Tables.rowtable(DBInterface.execute(ctx.db,
                """SELECT id FROM indices WHERE exposure_id = ?
                   ORDER BY score DESC NULLS LAST, id ASC LIMIT 1""",
                [ctx.exposure_id])))
            DBInterface.execute(ctx.db,
                "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)",
                [ctx.exposure_id, Int(top.id)])
            snap = HimalayaUI.compute_member_snapshot(ctx.db, ctx.exposure_id)
            @test snap[:confirmed_index] !== nothing
            ci = snap[:confirmed_index]
            @test ci[:id] == Int(top.id)
            @test ci[:phase] isa AbstractString
            @test ci[:peak_ids] isa AbstractVector

            # A low-R² member stays confirmed — an explicit assignment supersedes
            # the old auto-confirm R² heuristic (no gate post-D-10).
            DBInterface.execute(ctx.db,
                "UPDATE indices SET r_squared = 0.10 WHERE exposure_id = ?",
                [ctx.exposure_id])
            @test HimalayaUI.compute_member_snapshot(
                ctx.db, ctx.exposure_id)[:confirmed_index] !== nothing

            # form_factor / null state → no representative index (state-gated).
            DBInterface.execute(ctx.db,
                """INSERT INTO assignments (exposure_id, state) VALUES (?, 'form_factor')
                   ON CONFLICT(exposure_id) DO UPDATE SET state = 'form_factor'""",
                [ctx.exposure_id])
            @test HimalayaUI.compute_member_snapshot(
                ctx.db, ctx.exposure_id)[:confirmed_index] === nothing

            # An assignment narrowed to one explicit member → that member is the
            # confirmed_index, and confirmed_index tracks assignment_members (not
            # any legacy group).
            DBInterface.execute(ctx.db,
                "UPDATE assignments SET state = 'indexed' WHERE exposure_id = ?",
                [ctx.exposure_id])
            DBInterface.execute(ctx.db,
                "DELETE FROM assignment_members WHERE exposure_id = ?", [ctx.exposure_id])
            DBInterface.execute(ctx.db,
                """INSERT INTO indices (exposure_id, phase, basis, score, r_squared, lattice_d, status, kind)
                   VALUES (?, 'Pn3m', 0.10, 0.42, 0.20, 100.0, 'candidate', 'auto')""",
                [ctx.exposure_id])
            only_id = first(Tables.rowtable(DBInterface.execute(ctx.db,
                "SELECT id FROM indices WHERE exposure_id=? AND phase='Pn3m' ORDER BY id DESC LIMIT 1",
                [ctx.exposure_id]))).id
            DBInterface.execute(ctx.db,
                "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)",
                [ctx.exposure_id, Int(only_id)])
            ci2 = HimalayaUI.compute_member_snapshot(ctx.db, ctx.exposure_id)[:confirmed_index]
            @test ci2 !== nothing
            @test ci2[:id] == Int(only_id)
            @test ci2[:phase] == "Pn3m"
            @test ci2[:lattice_d] == 100.0
        end
    end

    # Plan E (E-4/E-7): the snapshot carries the durable assignment STATE
    # (indexed | form_factor | null) so the Series surface can distinguish a
    # form-factor member from a null member (confirmed_index is null for BOTH),
    # plus `confirmed_phases` — the distinct phases assigned to the member —
    # so coexistence reads/rows/strip cells can self-decode.
    @testset "compute_member_snapshot: assignment_state + confirmed_phases" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp; datfile="cubic_tot.dat",
                                            filename="cubic_tot")
            # Default (no assignments row) → state "indexed".
            snap0 = HimalayaUI.compute_member_snapshot(ctx.db, ctx.exposure_id)
            @test snap0[:assignment_state] == "indexed"
            @test snap0[:confirmed_phases] isa AbstractVector

            # Drive the member to form_factor via the durable assignment table.
            DBInterface.execute(ctx.db,
                """INSERT INTO assignments (exposure_id, state) VALUES (?, 'form_factor')
                   ON CONFLICT(exposure_id) DO UPDATE SET state = 'form_factor'""",
                [ctx.exposure_id])
            snap_ff = HimalayaUI.compute_member_snapshot(ctx.db, ctx.exposure_id)
            @test snap_ff[:assignment_state] == "form_factor"
            # A form_factor member carries NO lattice phases and NO representative
            # index — both confirmed_phases and confirmed_index are state-gated to
            # empty/null even if durable member rows linger (the auto-seeded
            # assignment members are not cleared on a state flip).
            @test isempty(snap_ff[:confirmed_phases])
            @test snap_ff[:confirmed_index] === nothing

            # And to null — distinct from form_factor though both have a null
            # confirmed_index; also carries no lattice phases.
            DBInterface.execute(ctx.db,
                "UPDATE assignments SET state = 'null' WHERE exposure_id = ?",
                [ctx.exposure_id])
            snap_null = HimalayaUI.compute_member_snapshot(ctx.db, ctx.exposure_id)
            @test snap_null[:assignment_state] == "null"
            @test isempty(snap_null[:confirmed_phases])
            @test snap_null[:confirmed_index] === nothing

            # ── Coexistence: two distinct-phase assigned indices → two
            #    confirmed_phases, EACH with its own lattice_d, ordered score-desc.
            DBInterface.execute(ctx.db,
                "UPDATE assignments SET state = 'indexed' WHERE exposure_id = ?",
                [ctx.exposure_id])
            # Clear any migration-backfilled members (from the active auto-group)
            # so this asserts an UNAMBIGUOUS two-phase coexistence, not whatever
            # the fixture's auto group happened to contain.
            DBInterface.execute(ctx.db,
                "DELETE FROM assignment_members WHERE exposure_id = ?", [ctx.exposure_id])
            # Insert two explicit distinct-phase indices with known lattice + score.
            DBInterface.execute(ctx.db,
                """INSERT INTO indices (exposure_id, phase, basis, score, r_squared, lattice_d, status, kind)
                   VALUES (?, 'Pn3m', 0.10, 0.90, 0.99, 100.0, 'candidate', 'auto')""",
                [ctx.exposure_id])
            id_pn = first(Tables.rowtable(DBInterface.execute(ctx.db,
                "SELECT id FROM indices WHERE exposure_id=? AND phase='Pn3m' ORDER BY id DESC LIMIT 1", [ctx.exposure_id]))).id
            DBInterface.execute(ctx.db,
                """INSERT INTO indices (exposure_id, phase, basis, score, r_squared, lattice_d, status, kind)
                   VALUES (?, 'Im3m', 0.10, 0.50, 0.99, 128.0, 'candidate', 'auto')""",
                [ctx.exposure_id])
            id_im = first(Tables.rowtable(DBInterface.execute(ctx.db,
                "SELECT id FROM indices WHERE exposure_id=? AND phase='Im3m' ORDER BY id DESC LIMIT 1", [ctx.exposure_id]))).id
            DBInterface.execute(ctx.db,
                """INSERT OR IGNORE INTO assignment_members (exposure_id, index_id) VALUES (?, ?)""",
                [ctx.exposure_id, Int(id_pn)])
            DBInterface.execute(ctx.db,
                """INSERT OR IGNORE INTO assignment_members (exposure_id, index_id) VALUES (?, ?)""",
                [ctx.exposure_id, Int(id_im)])
            snap_cx = HimalayaUI.compute_member_snapshot(ctx.db, ctx.exposure_id)
            cps = snap_cx[:confirmed_phases]
            @test length(cps) == 2
            # Both phases present, each with its OWN lattice_d (the point of the
            # per-phase shape for coexistence reads/rows).
            by_phase = Dict(cp[:phase] => cp[:lattice_d] for cp in cps)
            @test by_phase["Pn3m"] == 100.0
            @test by_phase["Im3m"] == 128.0
            # Ordered score-desc (Pn3m score 0.90 > Im3m 0.50) so the dominant
            # phase reads first.
            @test cps[1][:phase] == "Pn3m"
        end
    end

    # D-10 retired the confirmed_phases "legacy fallback" (an indexed member with
    # a confirmed_index but no durable assignment_members): confirmed_index now
    # ALSO sources from assignment_members, so with zero members both are empty —
    # there is nothing to fall back FROM. A never-migrated, never-seeded exposure
    # simply reads as indexed/empty until its assignment is seeded.
    @testset "compute_member_snapshot: indexed but no assignment members → empty" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp; datfile="cubic_tot.dat",
                                            filename="cubic_tot")
            DBInterface.execute(ctx.db,
                "DELETE FROM assignment_members WHERE exposure_id = ?", [ctx.exposure_id])
            snap = HimalayaUI.compute_member_snapshot(ctx.db, ctx.exposure_id)
            @test snap[:assignment_state] == "indexed"
            @test snap[:confirmed_index] === nothing
            @test isempty(snap[:confirmed_phases])
        end
    end

    @testset "is_member_stale" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            # Build a fake member_row NamedTuple-style (mirrors what
            # fetch_comparison_with_members feeds into the helper).
            current = HimalayaUI.read_inputs_hash(ctx.db, ctx.exposure_id)
            # Fresh: snapshot hash matches.
            fresh_snap = JSON3.write(Dict(:analysis_inputs_hash => current))
            fresh_row = (exposure_id = ctx.exposure_id, snapshot = fresh_snap)
            @test HimalayaUI.is_member_stale(ctx.db, fresh_row) == false

            # Stale: snapshot hash differs.
            stale_row = (exposure_id = ctx.exposure_id,
                         snapshot = JSON3.write(Dict(:analysis_inputs_hash => "sha256:other")))
            @test HimalayaUI.is_member_stale(ctx.db, stale_row) == true

            # Orphan: exposure_id NULL → not stale (no live exposure).
            orphan_row = (exposure_id = missing, snapshot = stale_row.snapshot)
            @test HimalayaUI.is_member_stale(ctx.db, orphan_row) == false
        end
    end

    # I3.6 (#177): `fetch_comparison_with_members` was a Compare-route-only
    # response builder, deleted with `routes_comparisons.jl`; its tests went
    # with it. `is_member_stale` (above) keeps its own direct test.

    @testset "member_ids_for_comparison" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 33)
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            HimalayaUI.apply_event!(ctx.db, req;
                kind="comparison_created", entity_type="comparison", entity_id=33,
                payload=Dict(:title => "T", :description => nothing,
                             :forked_from_id => nothing, :forked_at_hash => nothing,
                             :members => [
                                _member_payload(exposure_id=ctx.exposure_id, display_order=0),
                                _member_payload(exposure_id=ctx.exposure_id, display_order=1),
                             ]))
            ids = HimalayaUI.member_ids_for_comparison(ctx.db, 33)
            @test length(ids) == 2
            @test all(id -> id isa Int, ids)
        end
    end

    @testset "recently_used_exposures" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            # Add two more exposures and seed members for alice across them.
            e2 = HimalayaUI.create_exposure!(ctx.db; sample_id=ctx.sample_id, filename="e2")
            e3 = HimalayaUI.create_exposure!(ctx.db; sample_id=ctx.sample_id, filename="e3")
            alice = HimalayaUI.get_or_create_user!(ctx.db, "alice")
            _premint_cmp!(ctx.db, 44)
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            HimalayaUI.apply_event!(ctx.db, req;
                kind="comparison_created", entity_type="comparison", entity_id=44,
                payload=Dict(:title => "T", :description => nothing,
                             :forked_from_id => nothing, :forked_at_hash => nothing,
                             :members => [
                                _member_payload(exposure_id=ctx.exposure_id, display_order=0),
                                _member_payload(exposure_id=e2, display_order=1),
                                _member_payload(exposure_id=e3, display_order=2),
                             ]))
            recent = HimalayaUI.recently_used_exposures(ctx.db, alice; limit=20)
            @test sort(recent) == sort([ctx.exposure_id, e2, e3])
            # Limit honored.
            recent2 = HimalayaUI.recently_used_exposures(ctx.db, alice; limit=2)
            @test length(recent2) == 2
        end
    end

    @testset "recently_used_exposures: excludes orphan exposure_id IS NULL" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            alice = HimalayaUI.get_or_create_user!(ctx.db, "alice")
            _premint_cmp!(ctx.db, 55)
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            HimalayaUI.apply_event!(ctx.db, req;
                kind="comparison_created", entity_type="comparison", entity_id=55,
                payload=Dict(:title => "T", :description => nothing,
                             :forked_from_id => nothing, :forked_at_hash => nothing,
                             :members => [_member_payload(exposure_id=ctx.exposure_id)]))
            DBInterface.execute(ctx.db,
                "UPDATE comparison_members SET exposure_id = NULL WHERE comparison_id = ?",
                [55])
            @test HimalayaUI.recently_used_exposures(ctx.db, alice) == Int[]
        end
    end

    # I3.6 (#177): `comparisons_for_experiment` + `forks_of_comparison` were
    # Compare-route-only listing functions, deleted with `routes_comparisons.jl`.
    # Their tests went with them.

    @testset "get_user_id_for_request" begin
        ctx = _setup_db()
        try
            req_anon = HTTP.Request("GET", "/", HTTP.Header[], UInt8[])
            @test HimalayaUI.get_user_id_for_request(ctx.db, req_anon) === nothing
            req_alice = HTTP.Request("GET", "/", ["X-Username" => "alice"], UInt8[])
            uid = HimalayaUI.get_user_id_for_request(ctx.db, req_alice)
            @test uid isa Int
            # Idempotent: same call returns the same id.
            @test HimalayaUI.get_user_id_for_request(ctx.db, req_alice) == uid
        finally
            close(ctx.db)
        end
    end

    @testset "canonical_json: alphabetical key ordering (matches JS)" begin
        # Pin the spec: `canonical_json` sorts keys at every nesting level
        # so cross-language hash parity holds against
        # `frontend/src/lib/comparison/contentHash.ts::canonicalJson`.
        @test HimalayaUI.canonical_json(Dict(:b => 2, :a => 1)) == "{\"a\":1,\"b\":2}"
        @test HimalayaUI.canonical_json(Dict(:z => Dict(:b => 2, :a => 1), :a => 1)) ==
              "{\"a\":1,\"z\":{\"a\":1,\"b\":2}}"
        # Arrays preserve order
        @test HimalayaUI.canonical_json([3, 1, 2]) == "[3,1,2]"
        # nothing/missing → null
        @test HimalayaUI.canonical_json(nothing) == "null"
        @test HimalayaUI.canonical_json(missing) == "null"
        # Float-canonicalization: integer-valued floats drop the .0 to match JS
        @test HimalayaUI.canonical_json(1.0) == "1"
        @test HimalayaUI.canonical_json(-2.0) == "-2"
        @test HimalayaUI.canonical_json(1.5) == "1.5"
        @test HimalayaUI.canonical_json(0.05) == "0.05"
    end

    @testset "cross-language contentHash fixture parity" begin
        # Spec: chat citations like `@comparison:42@<hash8>` resolve to a
        # frozen snapshot. If client and server canonicalize differently,
        # the @-resolver mismatches. The fixture in
        # `frontend/test/fixtures/contentHash.fixture.json` pins both
        # implementations to byte-identical output for the same input.
        # Hashed by both sides; this test is the Julia leg.
        fixture_path = joinpath(@__DIR__, "..", "frontend", "test",
                                "fixtures", "contentHash.fixture.json")
        @test isfile(fixture_path)
        fixture = JSON3.read(read(fixture_path, String))
        @test haskey(fixture, :input)
        @test haskey(fixture, :expected_canonical)
        @test haskey(fixture, :expected_hash)

        canonical = HimalayaUI.canonical_json(fixture.input)
        @test canonical == String(fixture.expected_canonical)

        # The hash also matches — the canonical bytes determine it.
        actual_hash = "sha256:" * bytes2hex(SHA.sha256(canonical))
        @test actual_hash == String(fixture.expected_hash)
    end

    # I3.6 (#177): the Phase 9.6 stale-flip + low-R² snapshot regression
    # testsets observed staleness through `fetch_comparison_with_members`
    # (now deleted). The `is_member_stale` helper keeps its own direct test
    # (above), and confirmed_index sourcing is covered by "compute_member_snapshot:
    # confirmed_index sources from the assignment".
end

# I3.6 (#177): the "Comparisons REST routes" testset block exercised the routes
# in `routes_comparisons.jl` (deleted). The kept `comparison_*` dispatcher
# branches + idempotency-replay are still covered by `test_events.jl` and
# `test_idempotency_replay_invariant.jl`.
@testset "Comparisons REST routes (retired with the Compare routes, #177)" begin
    # All route tests here exercised `routes_comparisons.jl` (deleted): create,
    # get, listing, forks, submit (incl. 403/404/409 author + conflict gates),
    # delete + cascade, messages, and idempotent-replay of the submit op. The
    # kept `comparison_*` dispatcher branches + idempotency replay are still
    # covered by `test_events.jl` and `test_idempotency_replay_invariant.jl`;
    # sample-message routing is covered by `test_routes_messages.jl`.
    @test true  # placeholder so the testset is non-empty
end

@testset "_count_distinct_phases — Compare UX F" begin
    f = HimalayaUI._count_distinct_phases
    @test f("") == 0
    @test f("Pn3m#0") == 1
    @test f("Pn3m#0|Pn3m#1") == 1                       # dups collapse
    @test f("Pn3m#0|Hex#1|Lam#2") == 3
    @test f("Pn3m#0|Im3m#1|Ia3d#2|Hex#3|Lam#4") == 5    # exceeds the top-3 cap
    @test f("Pn3m#0||Hex#1") == 2                        # empty token skipped
    @test f("malformed") == 0                            # no '#' → skipped
end

# I3.6 (#177): the listing-projection / forks_of_comparison / view_* / NULL-hash
# `has_stale_members` testsets exercised `comparisons_listing` /
# `forks_of_comparison` / `fetch_comparison_with_members` (all deleted with the
# Compare routes). `_count_distinct_phases` (kept; reused by series) keeps its
# Compare-UX-F test above; `is_member_stale` keeps its direct test in the
# helpers block.

