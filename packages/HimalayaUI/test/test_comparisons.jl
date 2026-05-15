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
    db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
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
    db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename=filename)
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
    (db = db, experiment_id = exp_id, sample_id = s_id,
     exposure_id = e_id, analysis_dir = analysis_dir)
end

@testset "comparisons.jl helpers" begin

    @testset "current_content_hash" begin
        ctx = _setup_db()
        try
            @test HimalayaUI.current_content_hash(ctx.db, 1) === nothing
            _premint_cmp!(ctx.db, 1)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET content_hash = ? WHERE id = ?",
                ["sha256:abc", 1])
            @test HimalayaUI.current_content_hash(ctx.db, 1) == "sha256:abc"
        finally
            close(ctx.db)
        end
    end

    @testset "is_author: comparison missing → false" begin
        ctx = _setup_db()
        try
            uid = HimalayaUI.get_or_create_user!(ctx.db, "alice")
            @test HimalayaUI.is_author(ctx.db, 999, uid) == false
        finally
            close(ctx.db)
        end
    end

    @testset "is_author: nothing user_id → false" begin
        ctx = _setup_db()
        try
            uid = HimalayaUI.get_or_create_user!(ctx.db, "alice")
            _premint_cmp!(ctx.db, 1)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET created_by = ? WHERE id = ?", [uid, 1])
            @test HimalayaUI.is_author(ctx.db, 1, nothing) == false
        finally
            close(ctx.db)
        end
    end

    @testset "is_author: NULL created_by → false (orphan author)" begin
        ctx = _setup_db()
        try
            uid = HimalayaUI.get_or_create_user!(ctx.db, "alice")
            _premint_cmp!(ctx.db, 1)
            # created_by stays NULL.
            @test HimalayaUI.is_author(ctx.db, 1, uid) == false
        finally
            close(ctx.db)
        end
    end

    @testset "is_author: matching user → true; different user → false" begin
        ctx = _setup_db()
        try
            alice = HimalayaUI.get_or_create_user!(ctx.db, "alice")
            bob   = HimalayaUI.get_or_create_user!(ctx.db, "bob")
            _premint_cmp!(ctx.db, 1)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET created_by = ? WHERE id = ?", [alice, 1])
            @test HimalayaUI.is_author(ctx.db, 1, alice) == true
            @test HimalayaUI.is_author(ctx.db, 1, bob)   == false
        finally
            close(ctx.db)
        end
    end

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

    @testset "compute_member_snapshot: confirmed_index R²-gated" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp; datfile="cubic_tot.dat",
                                            filename="cubic_tot")
            # Force every index below the gate first → confirmed_index = nothing.
            DBInterface.execute(ctx.db,
                "UPDATE indices SET r_squared = 0.5 WHERE exposure_id = ?",
                [ctx.exposure_id])
            # Pick the highest-scored index, confirm it via the standard route
            # path (ensure_custom_group + insert member).
            ix = first(Tables.rowtable(DBInterface.execute(ctx.db,
                """SELECT id FROM indices WHERE exposure_id = ?
                   ORDER BY score DESC LIMIT 1""", [ctx.exposure_id])))
            ix_id = Int(ix.id)
            (custom_id, _) = HimalayaUI.ensure_custom_group!(ctx.db, ctx.exposure_id)
            # Ensure the custom group is the active one for this test (the
            # helper requires g.active = 1). `ensure_custom_group!` clones
            # the auto group's members into the custom group, so the index
            # may already be present — use OR IGNORE.
            DBInterface.execute(ctx.db,
                "UPDATE index_groups SET active = 1 WHERE id = ?", [custom_id])
            DBInterface.execute(ctx.db,
                "INSERT OR IGNORE INTO index_group_members (group_id, index_id) VALUES (?, ?)",
                [custom_id, ix_id])

            snap_below = HimalayaUI.compute_member_snapshot(ctx.db, ctx.exposure_id)
            @test snap_below[:confirmed_index] === nothing  # R²=0.5 below gate

            # Bump that single index above the gate.
            DBInterface.execute(ctx.db,
                "UPDATE indices SET r_squared = 0.99 WHERE id = ?", [ix_id])
            snap = HimalayaUI.compute_member_snapshot(ctx.db, ctx.exposure_id)
            @test snap[:confirmed_index] !== nothing
            ci = snap[:confirmed_index]
            @test ci[:id] == ix_id
            @test ci[:phase] isa AbstractString
            @test ci[:r_squared] >= HimalayaUI.CONFIRMED_INDEX_R2_GATE
            @test ci[:peak_ids] isa AbstractVector
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

    @testset "fetch_comparison_with_members: missing → nothing" begin
        ctx = _setup_db()
        try
            @test HimalayaUI.fetch_comparison_with_members(ctx.db, 999) === nothing
        finally
            close(ctx.db)
        end
    end

    @testset "fetch_comparison_with_members: nested shape + is_stale per member" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            # Build a comparison directly via apply_event! (the live route path
            # used in production).
            _premint_cmp!(ctx.db, 11)
            current_hash = HimalayaUI.read_inputs_hash(ctx.db, ctx.exposure_id)
            fresh_snap = Dict(:effective_peaks => [], :confirmed_index => nothing,
                              :analysis_inputs_hash => current_hash)
            stale_snap = Dict(:effective_peaks => [], :confirmed_index => nothing,
                              :analysis_inputs_hash => "sha256:older")
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            HimalayaUI.apply_event!(ctx.db, req;
                kind="comparison_created", entity_type="comparison", entity_id=11,
                payload=Dict(:title => "T", :description => "d",
                             :forked_from_id => nothing, :forked_at_hash => nothing,
                             :members => [
                                _member_payload(exposure_id=ctx.exposure_id,
                                                display_order=0, snapshot=fresh_snap),
                                _member_payload(exposure_id=ctx.exposure_id,
                                                display_order=1, snapshot=stale_snap),
                             ]))

            out = HimalayaUI.fetch_comparison_with_members(ctx.db, 11)
            @test out !== nothing
            @test out[:id] == 11
            @test out[:title] == "T"
            @test out[:description] == "d"
            @test startswith(out[:content_hash], "sha256:")
            @test out[:created_by] !== nothing
            @test out[:forked_from_id] === nothing
            @test length(out[:members]) == 2
            @test out[:members][1][:is_stale] == false
            @test out[:members][2][:is_stale] == true
            for m in out[:members]
                @test m[:exposure_id] == ctx.exposure_id
                @test m[:snapshot] !== nothing  # parsed JSON
            end
        end
    end

    @testset "fetch_comparison_with_members: orphan member with NULL exposure_id" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 22)
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            HimalayaUI.apply_event!(ctx.db, req;
                kind="comparison_created", entity_type="comparison", entity_id=22,
                payload=Dict(:title => "T", :description => nothing,
                             :forked_from_id => nothing, :forked_at_hash => nothing,
                             :members => [_member_payload(exposure_id=ctx.exposure_id)]))
            # Force the exposure_id to NULL (simulates exposure deletion under
            # ON DELETE SET NULL).
            DBInterface.execute(ctx.db,
                "UPDATE comparison_members SET exposure_id = NULL WHERE comparison_id = ?",
                [22])
            out = HimalayaUI.fetch_comparison_with_members(ctx.db, 22)
            @test out !== nothing
            @test length(out[:members]) == 1
            @test out[:members][1][:exposure_id] === nothing
            @test out[:members][1][:is_stale] == false  # nothing to drift from
        end
    end

    @testset "fetch_comparison_with_members: lineage info (forked)" begin
        ctx = _setup_db()
        try
            _premint_cmp!(ctx.db, 1)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET title = 'parent' WHERE id = ?", [1])
            _premint_cmp!(ctx.db, 2)
            DBInterface.execute(ctx.db,
                """UPDATE comparisons SET forked_from_id = 1, forked_at_hash = 'sha256:p',
                                          title = 'child' WHERE id = ?""", [2])
            out = HimalayaUI.fetch_comparison_with_members(ctx.db, 2)
            @test out[:forked_from_id]    == 1
            @test out[:forked_at_hash]    == "sha256:p"
            @test out[:forked_from_title] == "parent"
        finally
            close(ctx.db)
        end
    end

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

    @testset "comparisons_for_experiment" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            # Build a comparison whose member references the only exposure under
            # this experiment.
            _premint_cmp!(ctx.db, 100)
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            HimalayaUI.apply_event!(ctx.db, req;
                kind="comparison_created", entity_type="comparison", entity_id=100,
                payload=Dict(:title => "First", :description => nothing,
                             :forked_from_id => nothing, :forked_at_hash => nothing,
                             :members => [_member_payload(exposure_id=ctx.exposure_id)]))
            list = HimalayaUI.comparisons_for_experiment(ctx.db, ctx.experiment_id)
            @test length(list) == 1
            @test list[1][:id] == 100
            @test list[1][:title] == "First"

            # An unrelated experiment — should NOT appear.
            other_exp = HimalayaUI.create_experiment!(ctx.db; path="/y",
                data_dir="/y", analysis_dir="/y")
            @test HimalayaUI.comparisons_for_experiment(ctx.db, other_exp) == Dict{Symbol,Any}[]
        end
    end

    @testset "forks_of_comparison" begin
        ctx = _setup_db()
        try
            _premint_cmp!(ctx.db, 1)
            _premint_cmp!(ctx.db, 2)
            _premint_cmp!(ctx.db, 3)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET forked_from_id = 1, title = 'fork1' WHERE id = ?", [2])
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET forked_from_id = 1, title = 'fork2' WHERE id = ?", [3])
            forks = HimalayaUI.forks_of_comparison(ctx.db, 1)
            @test length(forks) == 2
            ids = Set(f[:id] for f in forks)
            @test ids == Set([2, 3])
            @test HimalayaUI.forks_of_comparison(ctx.db, 2) == Dict{Symbol,Any}[]
        finally
            close(ctx.db)
        end
    end

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

    # ─────────────────────────────────────────────────────────────────────
    # Phase 9.6 — backend stale-flip + low-R² snapshot regression tests.
    # ─────────────────────────────────────────────────────────────────────

    @testset "Phase 9.6 stale-flip: change underlying exposure → is_stale flips" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            current_hash = HimalayaUI.read_inputs_hash(ctx.db, ctx.exposure_id)
            @test current_hash isa AbstractString

            # Build a fresh comparison with the live snapshot hash.
            _premint_cmp!(ctx.db, 71)
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            HimalayaUI.apply_event!(ctx.db, req;
                kind="comparison_created", entity_type="comparison", entity_id=71,
                payload=Dict(:title => "Phase 9.6", :description => nothing,
                             :forked_from_id => nothing, :forked_at_hash => nothing,
                             :members => [
                                _member_payload(
                                    exposure_id=ctx.exposure_id,
                                    snapshot=Dict(:effective_peaks => [],
                                                  :confirmed_index => nothing,
                                                  :analysis_inputs_hash => current_hash)),
                             ]))

            # Sanity: fresh fetch reports is_stale=false.
            out = HimalayaUI.fetch_comparison_with_members(ctx.db, 71)
            @test length(out[:members]) == 1
            @test out[:members][1][:is_stale] == false

            # Drift the live exposure's hash (simulates a peak set change
            # under the hood without re-running the full pipeline).
            DBInterface.execute(ctx.db,
                "UPDATE exposures SET analysis_inputs_hash = ? WHERE id = ?",
                ["sha256:after-curation", ctx.exposure_id])

            out2 = HimalayaUI.fetch_comparison_with_members(ctx.db, 71)
            @test out2[:members][1][:is_stale] == true
        end
    end

    @testset "Phase 9.6 stale-flip: only affected member flips, others stay fresh" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            # A second exposure on the same sample, also analyzed. Both
            # members start fresh; we drift exposure A only.
            # Distinct filename — `exposures_unique_filename` UNIQUE INDEX
            # forbids `(sample_id, filename)` collisions; reuse the same .dat
            # via a copy so analyze_exposure! succeeds.
            cp(joinpath(ctx.analysis_dir, "example_tot.dat"),
               joinpath(ctx.analysis_dir, "example_tot_b.dat"))
            e2 = HimalayaUI.create_exposure!(ctx.db; sample_id=ctx.sample_id,
                                              filename="example_tot_b")
            HimalayaUI.analyze_exposure!(ctx.db, e2, ctx.analysis_dir)
            ha = HimalayaUI.read_inputs_hash(ctx.db, ctx.exposure_id)
            hb = HimalayaUI.read_inputs_hash(ctx.db, e2)

            _premint_cmp!(ctx.db, 72)
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            HimalayaUI.apply_event!(ctx.db, req;
                kind="comparison_created", entity_type="comparison", entity_id=72,
                payload=Dict(:title => "AB", :description => nothing,
                             :forked_from_id => nothing, :forked_at_hash => nothing,
                             :members => [
                                _member_payload(
                                    exposure_id=ctx.exposure_id, display_order=0,
                                    snapshot=Dict(:effective_peaks => [],
                                                  :confirmed_index => nothing,
                                                  :analysis_inputs_hash => ha)),
                                _member_payload(
                                    exposure_id=e2, display_order=1,
                                    snapshot=Dict(:effective_peaks => [],
                                                  :confirmed_index => nothing,
                                                  :analysis_inputs_hash => hb)),
                             ]))
            # Drift only A's hash.
            DBInterface.execute(ctx.db,
                "UPDATE exposures SET analysis_inputs_hash = ? WHERE id = ?",
                ["sha256:drifted-A", ctx.exposure_id])
            out = HimalayaUI.fetch_comparison_with_members(ctx.db, 72)
            @test length(out[:members]) == 2
            @test out[:members][1][:is_stale] == true   # A drifted
            @test out[:members][2][:is_stale] == false  # B unchanged
        end
    end

    @testset "Phase 9.6 low-R² confirmed_index snapshot path: below-gate → null" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp; datfile="cubic_tot.dat",
                                            filename="cubic_tot")
            # Force every index below the gate AND confirm one via a custom
            # active group, simulating a user who clicked Confirm on an
            # imperfect fit. The snapshot must NOT carry the index.
            DBInterface.execute(ctx.db,
                "UPDATE indices SET r_squared = 0.95 WHERE exposure_id = ?",
                [ctx.exposure_id])
            ix_id = first(Tables.rowtable(DBInterface.execute(ctx.db,
                "SELECT id FROM indices WHERE exposure_id = ? ORDER BY score DESC LIMIT 1",
                [ctx.exposure_id]))).id |> Int
            (custom_id, _) = HimalayaUI.ensure_custom_group!(ctx.db, ctx.exposure_id)
            DBInterface.execute(ctx.db,
                "UPDATE index_groups SET active = 1 WHERE id = ?", [custom_id])
            DBInterface.execute(ctx.db,
                "INSERT OR IGNORE INTO index_group_members (group_id, index_id) VALUES (?, ?)",
                [custom_id, ix_id])

            # The snapshot computed at submit time MUST reflect the gate
            # (R² >= 0.98). 0.95 is below → null.
            snap = HimalayaUI.compute_member_snapshot(ctx.db, ctx.exposure_id)
            @test snap[:confirmed_index] === nothing

            # Submit the comparison with this snapshot; the fetch round-trip
            # surfaces null in the JSON shape too.
            _premint_cmp!(ctx.db, 73)
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            HimalayaUI.apply_event!(ctx.db, req;
                kind="comparison_created", entity_type="comparison", entity_id=73,
                payload=Dict(:title => "low R²", :description => nothing,
                             :forked_from_id => nothing, :forked_at_hash => nothing,
                             :members => [
                                _member_payload(
                                    exposure_id=ctx.exposure_id,
                                    snapshot=snap),
                             ]))
            out = HimalayaUI.fetch_comparison_with_members(ctx.db, 73)
            # `snapshot` JSON has `confirmed_index === nothing` — round-tripped.
            ms = out[:members][1][:snapshot]
            @test ms isa AbstractDict
            @test ms[:confirmed_index] === nothing
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 2 Task 2.2 — REST routes (`routes_comparisons.jl`).
#
# Covers happy + sad paths for every route per spec §REST API + the
# 409-on-stale-hash contract per spec lines 273-310. Idempotency-replay
# rows for the comparison kinds live in `test_idempotency_replay_invariant.jl`.
# ─────────────────────────────────────────────────────────────────────────────

# Build a request body for POST /api/comparisons.
function _create_body(exposure_id::Int; title::String = "T", members_extra...)
    snap = Dict(:effective_peaks => [], :confirmed_index => nothing,
                :analysis_inputs_hash => "sha256:zero")
    members = [Dict{Symbol,Any}(:exposure_id => exposure_id,
                                :display_order => 0,
                                :snapshot => snap)]
    Dict{Symbol, Any}(:title => title, :members => members, members_extra...)
end

@testset "Comparisons REST routes" begin

    @testset "POST /api/comparisons: 400 on zero members" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(Dict(:title => "T", :members => [])),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"],
                    status_exception = false)
                @test r.status == 400
            end
        end
    end

    @testset "POST /api/comparisons: 201 + canonical body shape" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                # Use the live analysis hash in the snapshot so is_stale is false.
                cur_hash = HimalayaUI.read_inputs_hash(ctx.db, ctx.exposure_id)
                fresh_snap = Dict(:effective_peaks => [], :confirmed_index => nothing,
                                  :analysis_inputs_hash => cur_hash)
                body_in = Dict{Symbol,Any}(:title => "T",
                    :members => [Dict(:exposure_id => ctx.exposure_id,
                                      :display_order => 0,
                                      :snapshot => fresh_snap)])
                r = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(body_in),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r.status == 201
                body = JSON3.read(String(r.body))
                @test body.id isa Integer
                @test body.title == "T"
                @test startswith(body.content_hash, "sha256:")
                @test body.created_by !== nothing
                @test body.forked_from_id === nothing
                @test length(body.members) == 1
                @test body.members[1].exposure_id == ctx.exposure_id
                @test body.members[1].is_stale == false
                @test body.members[1].snapshot !== nothing
            end
        end
    end

    @testset "POST /api/comparisons: dispatcher fallback fills missing snapshot from server" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                # Send a member WITHOUT a snapshot — the route should fall
                # back to compute_member_snapshot(db, exposure_id) before
                # apply_event! (the dispatcher errors on missing snapshot).
                no_snap_body = Dict{Symbol,Any}(:title => "T",
                    :members => [Dict(:exposure_id => ctx.exposure_id,
                                      :display_order => 0)])
                r = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(no_snap_body),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r.status == 201
                body = JSON3.read(String(r.body))
                snap = body.members[1].snapshot
                @test haskey(snap, :effective_peaks)
                @test haskey(snap, :confirmed_index)
                @test haskey(snap, :analysis_inputs_hash)
            end
        end
    end

    @testset "POST /api/comparisons: fork payload echoes forked_from_id + forked_at_hash" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                # First create a parent.
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id; title="parent")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                parent = JSON3.read(String(r1.body))
                # Now fork it.
                fork_body = _create_body(ctx.exposure_id; title="fork-of-parent",
                    forked_from_id = parent.id,
                    forked_at_hash = parent.content_hash)
                r2 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(fork_body),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "bob"])
                @test r2.status == 201
                fork = JSON3.read(String(r2.body))
                @test fork.forked_from_id    == parent.id
                @test fork.forked_at_hash    == parent.content_hash
                @test fork.forked_from_title == "parent"
            end
        end
    end

    @testset "GET /api/comparisons/:id: 404 when missing, 200 when present" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r404 = HTTP.get("$base/api/comparisons/9999";
                    status_exception = false)
                @test r404.status == 404

                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                cmp = JSON3.read(String(r1.body))
                r2 = HTTP.get("$base/api/comparisons/$(cmp.id)")
                @test r2.status == 200
                body = JSON3.read(String(r2.body))
                @test body.id == cmp.id
                @test length(body.members) == 1
            end
        end
    end

    @testset "GET /api/comparisons (global listing) and per-experiment listing" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                r = HTTP.get("$base/api/comparisons")
                @test r.status == 200
                @test length(JSON3.read(String(r.body))) == 1

                r2 = HTTP.get("$base/api/experiments/$(ctx.experiment_id)/comparisons")
                @test r2.status == 200
                @test length(JSON3.read(String(r2.body))) == 1
            end
        end
    end

    @testset "GET /api/comparisons/:id/forks lists forks (and is empty when none)" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id; title="parent")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                parent = JSON3.read(String(r1.body))

                forks_empty = JSON3.read(String(HTTP.get(
                    "$base/api/comparisons/$(parent.id)/forks").body))
                @test isempty(forks_empty)

                HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id; title="fork",
                        forked_from_id = parent.id, forked_at_hash = parent.content_hash)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "bob"])
                forks = JSON3.read(String(HTTP.get(
                    "$base/api/comparisons/$(parent.id)/forks").body))
                @test length(forks) == 1
            end
        end
    end

    @testset "POST /submit: 403 for non-author" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                cmp = JSON3.read(String(r1.body))

                # Bob tries to submit Alice's comparison → 403.
                r2 = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                    body = JSON3.write(Dict(:title => "T", :members => [],
                                            :expected_content_hash => cmp.content_hash)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "bob"],
                    status_exception = false)
                @test r2.status == 403
            end
        end
    end

    @testset "POST /submit: 403 for orphaned author (NULL created_by)" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                cmp = JSON3.read(String(r1.body))

                # Force created_by to NULL (simulates the user row being deleted).
                DBInterface.execute(ctx.db,
                    "UPDATE comparisons SET created_by = NULL WHERE id = ?", [cmp.id])

                # Even Alice can't submit now.
                r2 = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                    body = JSON3.write(Dict(:title => "T", :members => [],
                                            :expected_content_hash => cmp.content_hash)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"],
                    status_exception = false)
                @test r2.status == 403
                # Spec §Authorship: "no author" detail in the error body.
                # Body shape is just {error}; the orphan path simply 403s.
                err = JSON3.read(String(r2.body))
                @test haskey(err, :error)
            end
        end
    end

    @testset "POST /submit: 409 on expected_content_hash mismatch (full body shape)" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                cmp = JSON3.read(String(r1.body))
                m_id = cmp.members[1].id

                # Submit with a stale expected_content_hash → 409.
                r2 = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                    body = JSON3.write(Dict(:title => "T", :description => nothing,
                        :expected_content_hash => "sha256:oldgarbage",
                        :members => [Dict(:id => m_id, :exposure_id => ctx.exposure_id,
                                          :display_order => 0,
                                          :snapshot => Dict(:effective_peaks => [],
                                                            :confirmed_index => nothing,
                                                            :analysis_inputs_hash => "sha256:zero"))])),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"],
                    status_exception = false)
                @test r2.status == 409
                body = JSON3.read(String(r2.body))
                @test body.error == "conflict"
                @test body.current_hash == cmp.content_hash
                # current_state mirrors GET /api/comparisons/:id shape.
                @test body.current_state.id == cmp.id
                @test body.current_state.title == cmp.title
                @test haskey(body.current_state, :members)
                @test length(body.current_state.members) == 1
            end
        end
    end

    @testset "POST /submit: 200 + new state with correct expected_content_hash" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id; title="orig")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                cmp = JSON3.read(String(r1.body))
                m_id = cmp.members[1].id

                r2 = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                    body = JSON3.write(Dict(:title => "renamed", :description => "d2",
                        :expected_content_hash => cmp.content_hash,
                        :members => [Dict(:id => m_id, :exposure_id => ctx.exposure_id,
                                          :display_order => 0,
                                          :snapshot => Dict(:effective_peaks => [],
                                                            :confirmed_index => nothing,
                                                            :analysis_inputs_hash => "sha256:zero"))])),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r2.status == 200
                body = JSON3.read(String(r2.body))
                @test body.title == "renamed"
                @test body.description == "d2"
                @test body.content_hash != cmp.content_hash
            end
        end
    end

    @testset "POST /submit: empty analysis_inputs_hash in overwrite payload → is_stale" begin
        # Issue #94 — server-side pin for the bug class fixed by #74. The
        # frontend's ConflictModal::buildOverwritePayload prefetches cache
        # entries so a cold cache no longer lands `analysis_inputs_hash = ""`
        # on outgoing snapshots. This test fixes the *server-side*
        # consequence so the prefetch fix stays load-bearing: if anyone
        # ever relaxes `is_member_stale` to treat empty hash as "unknown,
        # preserve previous", the frontend prefetch quietly stops mattering
        # and bug #74 reappears without any frontend test catching it.
        #
        # Six-layer rule (docs/contract-testing.md): frontend layers (cache
        # shape, onMutate payload) are already pinned in Vitest. This row
        # covers the two missing server layers — submit-route response
        # body and the view-fold staleness derivation on follow-up GET.
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                cur_hash = HimalayaUI.read_inputs_hash(ctx.db, ctx.exposure_id)
                @test cur_hash isa AbstractString && !isempty(cur_hash)

                # Create with a fresh snapshot so the initial state is not stale.
                fresh_snap = Dict(:effective_peaks => [], :confirmed_index => nothing,
                                  :analysis_inputs_hash => cur_hash)
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(Dict(:title => "T",
                        :members => [Dict(:exposure_id => ctx.exposure_id,
                                          :display_order => 0,
                                          :snapshot => fresh_snap)])),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r1.status == 201
                cmp = JSON3.read(String(r1.body))
                @test cmp.members[1].is_stale == false
                m_id = cmp.members[1].id

                # Overwrite via /submit with analysis_inputs_hash = "". The
                # canonical staleness predicate `snapshot_hash != current_hash`
                # MUST treat "" as drift from the live, non-empty cur_hash.
                empty_snap = Dict(:effective_peaks => [], :confirmed_index => nothing,
                                  :analysis_inputs_hash => "")
                r2 = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                    body = JSON3.write(Dict(:title => "T", :description => nothing,
                        :expected_content_hash => cmp.content_hash,
                        :members => [Dict(:id => m_id, :exposure_id => ctx.exposure_id,
                                          :display_order => 0,
                                          :snapshot => empty_snap)])),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r2.status == 200
                submit_body = JSON3.read(String(r2.body))
                @test submit_body.members[1].is_stale == true

                # Follow-up GET pins the view-fold derivation specifically:
                # the empty-hash snapshot was persisted and `is_member_stale`
                # re-derives drift on every read.
                r3 = HTTP.get("$base/api/comparisons/$(cmp.id)")
                @test r3.status == 200
                got = JSON3.read(String(r3.body))
                @test got.members[1].is_stale == true
                # Pin literal empty-string persistence — guards against a
                # different failure mode than the staleness predicate above:
                # a future serializer that drops empty fields, or a SQLite
                # TEXT round-trip that coerces "" → NULL, would still leave
                # is_stale == true (NULL ≠ live hash) but break this assert.
                @test got.members[1].snapshot.analysis_inputs_hash == ""
            end
        end
    end

    @testset "DELETE /:id: 403 for non-author, 200 + cascade for author" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                cmp = JSON3.read(String(r1.body))

                # Post a chat message so the cascade can be observed.
                HTTP.post("$base/api/comparisons/$(cmp.id)/messages";
                    body = JSON3.write(Dict(:body => "hello")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])

                # Bob can't delete.
                r403 = HTTP.delete("$base/api/comparisons/$(cmp.id)";
                    headers = ["X-Username" => "bob"],
                    status_exception = false)
                @test r403.status == 403

                # Alice can.
                r2 = HTTP.delete("$base/api/comparisons/$(cmp.id)";
                    headers = ["X-Username" => "alice"])
                @test r2.status == 200
                body = JSON3.read(String(r2.body))
                @test body.deleted == true

                # Cascade: comparisons + members + messages all gone.
                @test HimalayaUI.fetch_comparison_with_members(ctx.db, cmp.id) === nothing
                n_members = first(Tables.rowtable(DBInterface.execute(ctx.db,
                    "SELECT COUNT(*) AS c FROM comparison_members WHERE comparison_id = ?",
                    [cmp.id]))).c
                @test n_members == 0
                n_msgs = first(Tables.rowtable(DBInterface.execute(ctx.db,
                    "SELECT COUNT(*) AS c FROM comparison_messages WHERE comparison_id = ?",
                    [cmp.id]))).c
                @test n_msgs == 0
            end
        end
    end

    @testset "DELETE /:id: 403 for orphaned author (NULL created_by)" begin
        # Spec §Authorship: when `created_by IS NULL` (the original author's
        # user row was deleted), the comparison enters a "fork-only" state.
        # No user matches `is_author`, so even the original actor — and any
        # other user — gets 403 on author-gated mutations. The frontend
        # disambiguates orphan vs non-author via `created_by === null` from
        # the GET payload (Phase 2 design note); the unified 403 keeps the
        # backend simple.
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                cmp = JSON3.read(String(r1.body))

                # Force created_by to NULL (simulates ON DELETE SET NULL).
                DBInterface.execute(ctx.db,
                    "UPDATE comparisons SET created_by = NULL WHERE id = ?", [cmp.id])

                # Even Alice (the original author) can't delete now.
                r403_alice = HTTP.delete("$base/api/comparisons/$(cmp.id)";
                    headers = ["X-Username" => "alice"],
                    status_exception = false)
                @test r403_alice.status == 403

                # And Bob certainly can't either.
                r403_bob = HTTP.delete("$base/api/comparisons/$(cmp.id)";
                    headers = ["X-Username" => "bob"],
                    status_exception = false)
                @test r403_bob.status == 403
            end
        end
    end

    @testset "GET /:id: orphaned-author comparison still readable by anyone" begin
        # "Fork-only" status from spec §Authorship: orphan-author comparisons
        # are not hidden — anyone can GET them so they can still be forked.
        # The frontend uses `created_by === null` from the GET payload to
        # disambiguate orphan from non-author and shows Fork (not Edit).
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                cmp = JSON3.read(String(r1.body))

                DBInterface.execute(ctx.db,
                    "UPDATE comparisons SET created_by = NULL WHERE id = ?", [cmp.id])

                # Bob can GET the orphan-author comparison.
                r_get = HTTP.get("$base/api/comparisons/$(cmp.id)";
                    headers = ["X-Username" => "bob"])
                @test r_get.status == 200
                body = JSON3.read(String(r_get.body))
                @test body.id == cmp.id
                @test body.created_by === nothing  # JSON null → nothing
            end
        end
    end

    # HTTP semantics: 404 must precede 403 — a request for a resource that
    # does not exist cannot be "forbidden" because there is nothing to
    # forbid. `is_author` returns false for missing comparisons, which
    # previously made author-gated routes return 403 in the missing case.
    # The fix probes existence via `current_content_hash === nothing` first.
    @testset "POST /submit: 404 (not 403) for never-existed comparison id" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.post("$base/api/comparisons/9999/submit";
                    body = JSON3.write(Dict(:title => "T", :members => [])),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"],
                    status_exception = false)
                @test r.status == 404
                err = JSON3.read(String(r.body))
                @test err.error == "comparison not found"
            end
        end
    end

    @testset "DELETE /:id: 404 (not 403) for never-existed comparison id" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.delete("$base/api/comparisons/9999";
                    headers = ["X-Username" => "alice"],
                    status_exception = false)
                @test r.status == 404
                err = JSON3.read(String(r.body))
                @test err.error == "comparison not found"
            end
        end
    end

    # Regression: 403 still fires when the comparison exists but the actor
    # is not the author (the most common author-gate path).
    @testset "POST /submit: 403 for existing comparison + non-author (regression)" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                cmp = JSON3.read(String(r1.body))

                # Bob tries to submit Alice's existing comparison → 403,
                # not 404 (the comparison exists, just not bob's).
                r2 = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                    body = JSON3.write(Dict(:title => "T", :members => [],
                                            :expected_content_hash => cmp.content_hash)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "bob"],
                    status_exception = false)
                @test r2.status == 403
            end
        end
    end

    @testset "GET / POST /api/comparisons/:id/messages" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                cmp = JSON3.read(String(r1.body))

                # Empty initially.
                r_empty = HTTP.get("$base/api/comparisons/$(cmp.id)/messages")
                @test r_empty.status == 200
                @test isempty(JSON3.read(String(r_empty.body)))

                # 401 without X-Username.
                r401 = HTTP.post("$base/api/comparisons/$(cmp.id)/messages";
                    body = JSON3.write(Dict(:body => "anon")),
                    headers = ["Content-Type" => "application/json"],
                    status_exception = false)
                @test r401.status == 401

                # 400 on empty body.
                r400 = HTTP.post("$base/api/comparisons/$(cmp.id)/messages";
                    body = JSON3.write(Dict(:body => "   ")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"],
                    status_exception = false)
                @test r400.status == 400

                # 201 + chat thread populated.
                r2 = HTTP.post("$base/api/comparisons/$(cmp.id)/messages";
                    body = JSON3.write(Dict(:body => "looks good")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "bob"])
                @test r2.status == 201
                msg = JSON3.read(String(r2.body))
                @test msg.author == "bob"
                @test msg.body == "looks good"
                @test msg.comparison_id == cmp.id

                r_list = HTTP.get("$base/api/comparisons/$(cmp.id)/messages")
                @test length(JSON3.read(String(r_list.body))) == 1
            end
        end
    end

    @testset "Idempotent retry of successful submit (200) → cached, no duplicate event" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id; title="orig")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                cmp = JSON3.read(String(r1.body))
                m_id = cmp.members[1].id

                op_id = "submit-replay-$(rand(UInt32))"
                payload_body = JSON3.write(Dict(:title => "renamed",
                    :description => nothing,
                    :expected_content_hash => cmp.content_hash,
                    :members => [Dict(:id => m_id, :exposure_id => ctx.exposure_id,
                                      :display_order => 0,
                                      :snapshot => Dict(:effective_peaks => [],
                                                        :confirmed_index => nothing,
                                                        :analysis_inputs_hash => "sha256:zero"))]))
                hdrs = ["Content-Type" => "application/json",
                        "X-Username" => "alice",
                        "X-Client-Op-Id" => op_id]

                pre_count = first(Tables.rowtable(DBInterface.execute(ctx.db,
                    "SELECT COUNT(*) AS c FROM user_actions WHERE action = 'comparison_submitted'"))).c

                r2 = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                    body = payload_body, headers = hdrs)
                @test r2.status == 200
                body1 = String(copy(r2.body))

                # Replay → same op_id → cached, byte-identical, no extra event row.
                r3 = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                    body = payload_body, headers = hdrs)
                @test r3.status == 200
                @test String(copy(r3.body)) == body1

                post_count = first(Tables.rowtable(DBInterface.execute(ctx.db,
                    "SELECT COUNT(*) AS c FROM user_actions WHERE action = 'comparison_submitted'"))).c
                @test post_count - pre_count == 1
            end
        end
    end

    # Phase 2 Task 2.3 — `post_message` event routing by entity_type.
    # Verifies the message routes write to distinct tables and that an
    # unknown entity_type at the apply_event! layer doesn't silently
    # corrupt either table. The route-level dispatch is what makes this
    # work — `update_view_for_event!` still no-ops on `post_message`
    # (sample_messages and comparison_messages are written by the route
    # handlers, not the dispatcher).
    @testset "post_message routing: comparison_message vs sample_message tables" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                cmp = JSON3.read(String(r1.body))

                # Sample message → sample_messages, NOT comparison_messages.
                HTTP.post("$base/api/samples/$(ctx.sample_id)/messages";
                    body = JSON3.write(Dict(:body => "sample-only")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                # Comparison message → comparison_messages, NOT sample_messages.
                HTTP.post("$base/api/comparisons/$(cmp.id)/messages";
                    body = JSON3.write(Dict(:body => "comparison-only")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])

                n_sample = first(Tables.rowtable(DBInterface.execute(ctx.db,
                    "SELECT COUNT(*) AS c FROM sample_messages WHERE sample_id = ?",
                    [ctx.sample_id]))).c
                n_compare = first(Tables.rowtable(DBInterface.execute(ctx.db,
                    "SELECT COUNT(*) AS c FROM comparison_messages WHERE comparison_id = ?",
                    [cmp.id]))).c
                @test n_sample == 1
                @test n_compare == 1

                # No comparison_message rows landed in sample_messages and
                # vice versa.
                bodies_sample = [String(r.body) for r in Tables.rowtable(DBInterface.execute(ctx.db,
                    "SELECT body FROM sample_messages WHERE sample_id = ?",
                    [ctx.sample_id]))]
                @test bodies_sample == ["sample-only"]
                bodies_compare = [String(r.body) for r in Tables.rowtable(DBInterface.execute(ctx.db,
                    "SELECT body FROM comparison_messages WHERE comparison_id = ?",
                    [cmp.id]))]
                @test bodies_compare == ["comparison-only"]

                # The user_actions log carries entity_type so replay can
                # reconstruct routing. Both kinds are 'post_message'.
                events = Tables.rowtable(DBInterface.execute(ctx.db,
                    "SELECT entity_type, entity_id FROM user_actions WHERE action = 'post_message'"))
                etypes = Set(String(e.entity_type) for e in events)
                @test "sample_message"     in etypes
                @test "comparison_message" in etypes
            end
        end
    end

    @testset "409 retry contract: status >= 400 NOT cached → conflict re-evaluates" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body(ctx.exposure_id; title="orig")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                cmp = JSON3.read(String(r1.body))
                m_id = cmp.members[1].id

                op_id = "conflict-retry-$(rand(UInt32))"
                stale_payload = JSON3.write(Dict(:title => "x",
                    :description => nothing,
                    :expected_content_hash => "sha256:stale",
                    :members => [Dict(:id => m_id, :exposure_id => ctx.exposure_id,
                                      :display_order => 0,
                                      :snapshot => Dict(:effective_peaks => [],
                                                        :confirmed_index => nothing,
                                                        :analysis_inputs_hash => "sha256:zero"))]))
                hdrs = ["Content-Type" => "application/json",
                        "X-Username" => "alice",
                        "X-Client-Op-Id" => op_id]

                # First call → 409 (stale hash). 4xx not cached.
                r_a = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                    body = stale_payload, headers = hdrs,
                    status_exception = false)
                @test r_a.status == 409

                # Second call → still 409, conflict re-evaluated (not a
                # cached replay of the prior 409 — the OP_LOCKS entry is
                # cleaned up after a 4xx, and the cache is empty).
                r_b = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                    body = stale_payload, headers = hdrs,
                    status_exception = false)
                @test r_b.status == 409

                # Now resolve the conflict by submitting with the correct
                # hash but a fresh op_id (the spec's prescribed escape valve).
                resolve_payload = JSON3.write(Dict(:title => "renamed",
                    :description => nothing,
                    :expected_content_hash => cmp.content_hash,
                    :members => [Dict(:id => m_id, :exposure_id => ctx.exposure_id,
                                      :display_order => 0,
                                      :snapshot => Dict(:effective_peaks => [],
                                                        :confirmed_index => nothing,
                                                        :analysis_inputs_hash => "sha256:zero"))]))
                r_c = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                    body = resolve_payload,
                    headers = ["Content-Type" => "application/json",
                               "X-Username" => "alice",
                               "X-Client-Op-Id" => "fresh-$(rand(UInt32))"])
                @test r_c.status == 200
            end
        end
    end
end

@testset "listing projection — Compare UX A-3" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))

        # Minimal fixture: user, experiment, sample, exposure, comparison
        # with one member whose snapshot pins a Pn3m index.
        DBInterface.execute(db, "INSERT INTO users (id, username) VALUES (1, 'alice')")
        DBInterface.execute(db, """INSERT INTO experiments
                                   (id, name, path, data_dir, analysis_dir)
                                   VALUES (10, 'exp', '/x', '/x/d', '/x/a')""")
        DBInterface.execute(db, """INSERT INTO samples (id, experiment_id, name)
                                   VALUES (100, 10, 'sA')""")
        DBInterface.execute(db, """INSERT INTO exposures
                                   (id, sample_id, filename, selected)
                                   VALUES (1000, 100, 'JC001', 1)""")
        DBInterface.execute(db, """INSERT INTO comparisons
                                   (id, title, content_hash, created_by, updated_at)
                                   VALUES (1, 'Cubic vs Hex', 'h', 1, '2026-05-14T10:00:00Z')""")
        snap = JSON3.write(Dict(
            :confirmed_index => Dict(:id => 1, :phase => "Pn3m"),
            :analysis_inputs_hash => "ih1",
        ))
        DBInterface.execute(db, """INSERT INTO comparison_members
                                   (comparison_id, exposure_id, display_order,
                                    band_height, y_offset, normalization, snapshot, created_at)
                                   VALUES (1, 1000, 0, 1.0, 0.0, 'max', ?,
                                           '2026-05-14T10:00:00Z')""",
                            [snap])

        rows = HimalayaUI.comparisons_listing(db)
        @test length(rows) == 1
        r = rows[1]
        @test r[:author_username]    == "alice"
        @test r[:member_count]       == 1
        @test r[:member_phases]      == ["Pn3m"]
        @test r[:has_stale_members]  == false
        @test r[:last_event_at] isa Union{String, Nothing}
        close(db)
    end
end

@testset "forks_of_comparison projects new aggregates (Compare UX A-3)" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
        DBInterface.execute(db, "INSERT INTO users (id, username) VALUES (1, 'alice')")
        DBInterface.execute(db, """INSERT INTO experiments
                                   (id, name, path, data_dir, analysis_dir)
                                   VALUES (10, 'exp', '/x', '/x/d', '/x/a')""")
        DBInterface.execute(db, """INSERT INTO samples (id, experiment_id, name)
                                   VALUES (100, 10, 'sA')""")
        DBInterface.execute(db, """INSERT INTO exposures
                                   (id, sample_id, filename, selected)
                                   VALUES (1000, 100, 'JC001', 1)""")
        # Parent comparison + one Pn3m member.
        DBInterface.execute(db, """INSERT INTO comparisons
                                   (id, title, content_hash, created_by, updated_at)
                                   VALUES (1, 'Parent', 'h1', 1, '2026-05-14T10:00:00Z')""")
        parent_snap = JSON3.write(Dict(
            :confirmed_index => Dict(:id => 1, :phase => "Pn3m"),
            :analysis_inputs_hash => "ih1"))
        DBInterface.execute(db, """INSERT INTO comparison_members
                                   (comparison_id, exposure_id, display_order,
                                    band_height, y_offset, normalization, snapshot, created_at)
                                   VALUES (1, 1000, 0, 1.0, 0.0, 'max', ?,
                                           '2026-05-14T10:00:00Z')""",
                            [parent_snap])
        # Fork (forked_from_id = 1) + one Hex member.
        DBInterface.execute(db, """INSERT INTO comparisons
                                   (id, title, content_hash, created_by, updated_at,
                                    forked_from_id)
                                   VALUES (2, 'Fork', 'h2', 1, '2026-05-14T11:00:00Z', 1)""")
        fork_snap = JSON3.write(Dict(
            :confirmed_index => Dict(:id => 2, :phase => "Hex"),
            :analysis_inputs_hash => "ih2"))
        DBInterface.execute(db, """INSERT INTO comparison_members
                                   (comparison_id, exposure_id, display_order,
                                    band_height, y_offset, normalization, snapshot, created_at)
                                   VALUES (2, 1000, 0, 1.0, 0.0, 'max', ?,
                                           '2026-05-14T11:00:00Z')""",
                            [fork_snap])

        forks = HimalayaUI.forks_of_comparison(db, 1)
        @test length(forks) == 1
        f = forks[1]
        @test f[:author_username]   == "alice"
        @test f[:member_count]      == 1
        @test f[:member_phases]     == ["Hex"]
        @test f[:has_stale_members] == false
        @test haskey(f, :view_grouping_mode)
        @test haskey(f, :last_event_at)
        close(db)
    end
end

