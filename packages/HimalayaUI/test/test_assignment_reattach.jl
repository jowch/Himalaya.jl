# F-WIPE W1: reanalysis must not silently wipe assignment members.
#
# `assignment_members.index_id` carries ON DELETE CASCADE, and persist_analysis!
# deletes every kind='auto' index on each reanalysis — so without a semantic
# snapshot/re-attach (mirroring the custom-group block), every auto assignment
# member silently vanished and seed_assignment_if_absent! then REPLACED a
# curated assignment with the fresh machine picks. These tests pin the fixed
# contract:
#   a) a curated all-auto assignment re-attaches by phase + nearest basis;
#   b) a mixed speculative+auto assignment keeps the speculative member
#      untouched, re-attaches the auto member, and is never re-seeded;
#   c) a member whose identity fails to re-attach is DROPPED — the assignment
#      stays EMPTY (announced via dropped_assignment_phases, no re-seed);
#   d) a fresh exposure with no members still seeds from the auto selection;
#   e) peak-route SSE frames carry post_state.assignment ({state, members},
#      same shape as the assignment_* frames) + assignment_dropped when the
#      drop case fires.
using Test
using HimalayaUI
using Himalaya
using SQLite, DBInterface, Tables
using SparseArrays
using JSON3, HTTP

# Standalone-run support: under runtests.jl the shared harness is already
# loaded (test_http.jl is included first); on a standalone run pull it in.
if !isdefined(@__MODULE__, :with_test_server)
    include("test_http.jl")
end

const _FWIPE_DAT = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")

"""Analyzed-exposure fixture: real trace, one analyze run (seeds the assignment)."""
function _fwipe_setup(tmp)
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(_FWIPE_DAT, joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp, "data"), analysis_dir=analysis_dir)
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
    (db = db, e_id = e_id, analysis_dir = analysis_dir)
end

_fwipe_members(db, e_id) = Tables.rowtable(DBInterface.execute(db, """
    SELECT m.index_id, i.phase, i.basis, i.kind
    FROM assignment_members m JOIN indices i ON i.id = m.index_id
    WHERE m.exposure_id = ? ORDER BY m.index_id""", [e_id]))

# In-process SSE capture (same pattern as test_idempotency_replay_invariant.jl's
# _capture_sse_during; distinct name so both files can coexist under runtests).
function _fwipe_capture_sse(f::Function, kind_filter::String)
    pending = Channel{String}(64)
    sub = (pending = pending,)
    lock(HimalayaUI.SSE_LOCK) do
        push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
    end
    try
        f()
        # Allow the post-commit broadcast queue to flush.
        sleep(0.3)
    finally
        lock(HimalayaUI.SSE_LOCK) do
            filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
        end
        close(pending)
    end
    frames = String[]
    for frame in pending
        startswith(frame, ":") && continue
        occursin("\"kind\":\"$kind_filter\"", frame) && push!(frames, frame)
    end
    frames
end

# "event: curation\ndata: {json}\n\n" → parsed JSON3 object.
_fwipe_frame_json(frame::String) = JSON3.read(split(frame, "data: "; limit=2)[2])

"""Hand-built 2-peak Pn3m candidate at `basis` (normalized-q₁-slope convention)."""
function _fwipe_pn3m_candidate(basis::Float64)
    rn = Himalaya.phaseratios(Himalaya.Pn3m; normalize = true)
    n  = length(rn)
    Himalaya.Index{Himalaya.Pn3m}(basis,
        SparseArrays.SparseVector{Float64, Int}(n, [1, 2], [basis * rn[1], basis * rn[2]]),
        SparseArrays.SparseVector{Float64, Int}(n, [1, 2], [1.0, 1.0]))
end

# Minimal direct-persist fixture: exposure + N curated same-phase auto members
# (no peaks/eff needed — the re-attach path never reads them).
function _fwipe_merge_setup(tmp, member_bases::Vector{Float64})
    db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))   # FK ON (cascade live)
    exp_id = HimalayaUI.create_experiment!(db; path=tmp,
        data_dir=joinpath(tmp, "data"), analysis_dir=joinpath(tmp, "analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    DBInterface.execute(db,
        "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [e_id])
    for b in member_bases
        res = DBInterface.execute(db,
            "INSERT INTO indices (exposure_id, phase, basis) VALUES (?, 'Pn3m', ?)",
            [e_id, b])
        DBInterface.execute(db,
            "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)",
            [e_id, Int(DBInterface.lastrowid(res))])
    end
    (db = db, e_id = e_id)
end

_fwipe_empty_eff() = (q = Float64[], sharpness = Float64[],
                      peak_id = Int[], peak_kind = Symbol[])
_fwipe_empty_pr()  = (q = Float64[], indices = Int[],
                      prominence = Float64[], sharpness = Float64[])

@testset "F-WIPE W1: assignment survives reanalysis" begin

    @testset "F-WIPE a: curated all-auto assignment re-attaches across reanalysis" begin
        # The example_tot fixture indexes to a SINGLE candidate, so a re-seed and
        # a re-attach are indistinguishable at the analyze_exposure! level. Pin
        # persist_analysis! directly with a controlled two-candidate set (the
        # real Pn3m + a hand-built Lamellar that claims one of Pn3m's peaks):
        # auto_group then selects only one of them, so the curated assignment
        # {both} differs from the machine selection — a silent re-seed (the
        # pre-fix behavior) shrinks the assignment back to the selection.
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))   # FK ON (cascade live)
            exp_id = HimalayaUI.create_experiment!(db; path=tmp,
                data_dir=joinpath(tmp, "data"), analysis_dir=joinpath(tmp, "analysis"))
            s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
            e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")

            q, I, σ = HimalayaUI.load_dat(_FWIPE_DAT)
            pr = Himalaya.findpeaks(q, I, σ)
            HimalayaUI.diff_update_auto_peaks!(db, e_id, pr, I)
            eff = HimalayaUI.effective_peaks(db, e_id, q, I)
            candidates = Himalaya.indexpeaks(eff.q, eff.sharpness)
            @test !isempty(candidates)

            # Hand-built Lamellar anchored on a peak the first candidate claims,
            # so the two candidates conflict and auto_group keeps only one.
            q1 = first(SparseArrays.findnz(candidates[1].peaks)[2])
            rl = Himalaya.phaseratios(Himalaya.Lamellar)
            nl = length(rl)
            lam = Himalaya.Index{Himalaya.Lamellar}(
                rl[1:2] \ [q1, 2q1],
                SparseArrays.SparseVector{Float64, Int}(nl, [1, 2], [q1, 2q1]),
                SparseArrays.SparseVector{Float64, Int}(nl, [1, 2], [1.0, 1.0]))
            push!(candidates, lam)
            group = HimalayaUI.auto_group(candidates, eff)
            @test length(group) < length(candidates)   # selection is a strict subset

            HimalayaUI.persist_analysis!(db, e_id, q, I, pr, candidates, group, eff)

            # Seeded from the auto selection; curate by assigning EVERY candidate
            # (the realistic assignment_add flow, by direct insert).
            for r in Tables.rowtable(DBInterface.execute(db,
                    "SELECT id FROM indices WHERE exposure_id = ? AND kind = 'auto'", [e_id]))
                DBInterface.execute(db,
                    "INSERT OR IGNORE INTO assignment_members (exposure_id, index_id) VALUES (?, ?)",
                    [e_id, Int(r.id)])
            end
            kept        = _fwipe_members(db, e_id)
            kept_phases = sort(String[String(r.phase) for r in kept])
            kept_ids    = Set(Int(r.index_id) for r in kept)
            @test length(kept) == length(candidates)
            @test length(kept) > length(group)          # curated ≠ machine selection

            # Reanalysis re-persists the same candidate set with fresh index ids.
            HimalayaUI.persist_analysis!(db, e_id, q, I, pr, candidates, group, eff)

            after        = _fwipe_members(db, e_id)
            after_phases = sort(String[String(r.phase) for r in after])
            # Same curated membership by semantic identity (phase + basis)…
            @test after_phases == kept_phases
            # …carried by FRESH index ids (auto indices were rebuilt; AUTOINCREMENT
            # never reuses the cascaded ids).
            @test isempty(intersect(Set(Int(r.index_id) for r in after), kept_ids))
        end
    end

    @testset "F-WIPE merge (i): two same-phase members collapsing onto one candidate shrink AND report" begin
        # Two curated Pn3m members (bases 0.100 / 0.103) both sit within
        # MEMBER_REATTACH_RELTOL of the single new Pn3m candidate at 0.101.
        # Without claim tracking, the second INSERT OR IGNORE is a silent PK
        # no-op: membership shrinks 2 → 1 with an EMPTY dropped list — the same
        # species of silent loss W1 exists to kill. The merge must be reported:
        # one per-member drop entry for the losing identity.
        mktempdir() do tmp
            ctx = _fwipe_merge_setup(tmp, [0.100, 0.103])
            db, e_id = ctx.db, ctx.e_id

            cand = _fwipe_pn3m_candidate(0.101)
            dropped = HimalayaUI.persist_analysis!(db, e_id, Float64[], Float64[],
                _fwipe_empty_pr(), [cand], [cand], _fwipe_empty_eff())

            after = _fwipe_members(db, e_id)
            @test length(after) == 1
            @test Float64(after[1].basis) ≈ 0.101
            # Per-member reporting: the merged-away identity counts as dropped,
            # exactly once (phases may repeat when several members of one phase drop).
            @test dropped == ["Pn3m"]
        end
    end

    @testset "F-WIPE merge (ii): claiming is deterministic — first identity in (phase, basis) order wins, no re-pairing" begin
        # Members 0.100 / 0.103; candidates 0.101 / 0.108. Snapshot order is
        # ORDER BY phase, basis, so 0.100 claims first: its nearest candidate is
        # 0.101 (Δ=0.001). 0.103's nearest is ALSO 0.101 (Δ=0.002 < Δ=0.005 to
        # 0.108) — already claimed, so 0.103 DROPS. There is deliberately no
        # fallback re-pairing onto the next-nearest 0.108: re-attach is an
        # identity match, not an assignment problem, and 0.108 within tol of
        # nothing it actually was must stay out of the membership.
        mktempdir() do tmp
            ctx = _fwipe_merge_setup(tmp, [0.100, 0.103])
            db, e_id = ctx.db, ctx.e_id

            cands = [_fwipe_pn3m_candidate(0.101), _fwipe_pn3m_candidate(0.108)]
            dropped = HimalayaUI.persist_analysis!(db, e_id, Float64[], Float64[],
                _fwipe_empty_pr(), cands, cands, _fwipe_empty_eff())

            after = _fwipe_members(db, e_id)
            @test length(after) == 1
            # The survivor is the first-sorted claimant's nearest candidate…
            @test Float64(after[1].basis) ≈ 0.101
            # …and the loser was NOT silently re-paired onto 0.108.
            @test all(!isapprox(Float64(r.basis), 0.108) for r in after)
            @test dropped == ["Pn3m"]
        end
    end

    @testset "F-WIPE b: mixed speculative+auto — speculative untouched, auto re-attaches, no re-seed" begin
        mktempdir() do tmp
            ctx = _fwipe_setup(tmp)
            db, e_id, dir = ctx.db, ctx.e_id, ctx.analysis_dir

            before = _fwipe_members(db, e_id)
            @test !isempty(before)
            # Curate down to exactly ONE auto member.
            keep = first(before)
            for r in before[2:end]
                DBInterface.execute(db,
                    "DELETE FROM assignment_members WHERE exposure_id = ? AND index_id = ?",
                    [e_id, Int(r.index_id)])
            end
            # Add a custom-committed (kind='speculative', stable-id) second member.
            P = Himalaya.Pn3m
            spec_basis = 2π / 150.0 * first(Himalaya.phaseratios(P))
            spec_id = HimalayaUI.insert_custom_index!(db, e_id, P, spec_basis)
            DBInterface.execute(db,
                "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)",
                [e_id, spec_id])

            DBInterface.execute(db,
                "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.123)",
                [e_id])
            HimalayaUI.analyze_exposure!(db, e_id, dir)

            after     = _fwipe_members(db, e_id)
            after_ids = Set(Int(r.index_id) for r in after)
            @test spec_id in after_ids        # speculative member survives untouched
            @test length(after) == 2          # spec + one re-attached auto; NOT the re-seeded full selection
            auto_after = [r for r in after if String(r.kind) == "auto"]
            @test length(auto_after) == 1
            @test String(auto_after[1].phase) == String(keep.phase)
        end
    end

    @testset "F-WIPE c: dropped member leaves assignment EMPTY (no re-seed) and reports the phase" begin
        mktempdir() do tmp
            ctx = _fwipe_setup(tmp)
            db, e_id, dir = ctx.db, ctx.e_id, ctx.analysis_dir

            before = _fwipe_members(db, e_id)
            @test !isempty(before)
            before_phases = sort(String[String(r.phase) for r in before])

            # Exclude every auto peak except one — `minpeaks` ≥ 2 for every phase,
            # so a single effective peak yields ZERO candidates and every assigned
            # identity fails to re-attach.
            auto_qs = Tables.rowtable(DBInterface.execute(db,
                "SELECT q FROM auto_peaks WHERE exposure_id = ? ORDER BY q", [e_id]))
            for r in auto_qs[2:end]
                DBInterface.execute(db,
                    "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', ?)",
                    [e_id, Float64(r.q)])
            end
            res = HimalayaUI.analyze_exposure!(db, e_id, dir)

            # Assignment goes EMPTY — announced, never silently re-seeded.
            @test isempty(_fwipe_members(db, e_id))
            # The dropped identities are reported up through analyze_exposure!.
            @test res isa NamedTuple
            @test sort(res.dropped_assignment_phases) == before_phases
        end
    end

    @testset "F-WIPE d: fresh exposure with no members — seeding behavior unchanged" begin
        mktempdir() do tmp
            ctx = _fwipe_setup(tmp)
            db, e_id = ctx.db, ctx.e_id

            # First analyze of a member-less exposure seeds from the auto selection.
            gid = Int(first(Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM index_groups WHERE exposure_id = ? AND kind = 'auto'",
                [e_id]))).id)
            auto_sel = sort(Int[Int(r.index_id) for r in Tables.rowtable(DBInterface.execute(db,
                "SELECT index_id FROM index_group_members WHERE group_id = ?", [gid]))])
            seeded = [Int(r.index_id) for r in _fwipe_members(db, e_id)]
            @test !isempty(seeded)
            @test seeded == auto_sel
            @test HimalayaUI._assignment_body(db, e_id)[:state] == "indexed"
        end
    end

    if isdefined(@__MODULE__, :with_test_server)
        @testset "F-WIPE e: peak-route SSE post_state carries assignment {state, members}" begin
            mktempdir() do tmp
                ctx = _fwipe_setup(tmp)
                db, e_id = ctx.db, ctx.e_id

                # A manual peak the route can DELETE (direct insert is fine in tests;
                # the route's reanalyze recomputes from the curation tables).
                res = DBInterface.execute(db,
                    "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.123)",
                    [e_id])
                cur_id = Int(DBInterface.lastrowid(res))

                with_test_server(db) do port, base
                    frames = _fwipe_capture_sse("peak_removed") do
                        r = HTTP.delete("$base/api/peaks/$cur_id",
                            ["X-Username" => "alice"])
                        @test r.status == 200
                    end
                    @test length(frames) == 1
                    evt = _fwipe_frame_json(frames[1])
                    ps  = evt.post_state
                    # Existing curation post_state contract still holds…
                    @test haskey(ps, :analysis_inputs_hash)
                    @test haskey(ps, :indices)
                    # …now ALSO the assignment, in the assignment_* frame shape.
                    @test haskey(ps, :assignment)
                    @test ps.assignment.state == "indexed"
                    @test sort(collect(Int, ps.assignment.members)) ==
                          HimalayaUI._assignment_body(db, e_id)[:members]
                    @test !isempty(ps.assignment.members)
                    # No drop fired → key omitted entirely.
                    @test !haskey(ps, :assignment_dropped)
                end
            end
        end

        @testset "F-WIPE e: drop case fires assignment_dropped on the peak_excluded frame" begin
            mktempdir() do tmp
                ctx = _fwipe_setup(tmp)
                db, e_id = ctx.db, ctx.e_id

                before_phases = sort(String[String(r.phase) for r in _fwipe_members(db, e_id)])
                @test !isempty(before_phases)

                # Exclude all auto peaks except the last two by direct insert, then
                # exclude one of the remaining two THROUGH the route — its reanalyze
                # sees a single effective peak, zero candidates, full member drop.
                auto_rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, q FROM auto_peaks WHERE exposure_id = ? ORDER BY q", [e_id]))
                @test length(auto_rows) >= 2
                for r in auto_rows[1:end-2]
                    DBInterface.execute(db,
                        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', ?)",
                        [e_id, Float64(r.q)])
                end
                patch_id = Int(auto_rows[end-1].id)

                with_test_server(db) do port, base
                    frames = _fwipe_capture_sse("peak_excluded") do
                        r = HTTP.patch("$base/api/peaks/$patch_id",
                            ["Content-Type" => "application/json", "X-Username" => "alice"],
                            JSON3.write(Dict(:excluded => true)))
                        @test r.status == 200
                    end
                    @test length(frames) == 1
                    ps = _fwipe_frame_json(frames[1]).post_state
                    @test haskey(ps, :assignment)
                    @test isempty(ps.assignment.members)
                    @test haskey(ps, :assignment_dropped)
                    @test sort(String[String(p) for p in ps.assignment_dropped]) == before_phases
                end
            end
        end

        @testset "F-WIPE f: manual analyze SSE frame carries assignment {state, members}" begin
            # The wipe/re-attach runs on the manual Re-analyze path too. With the
            # frontend's envelope-absent-means-do-not-touch rule, an analyze_run
            # frame without the assignment envelope strands every tab's cached
            # members on dead (cascaded) index ids until a hard refetch.
            mktempdir() do tmp
                ctx = _fwipe_setup(tmp)
                db, e_id = ctx.db, ctx.e_id

                # Change the inputs hash so the reanalyze is NOT a no-op.
                DBInterface.execute(db,
                    "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.123)",
                    [e_id])

                with_test_server(db) do port, base
                    frames = _fwipe_capture_sse("analyze_run") do
                        r = HTTP.post("$base/api/exposures/$e_id/analyze",
                            ["X-Username" => "alice",
                             "X-Client-Op-Id" => "op-fwipe-analyze-1"])
                        @test r.status == 200
                    end
                    @test length(frames) == 1
                    ps = _fwipe_frame_json(frames[1]).post_state
                    # Existing analyze_run post_state contract still holds…
                    @test haskey(ps, :analysis_inputs_hash)
                    @test haskey(ps, :indices)
                    # …now ALSO the assignment, mirroring _enrich_curation_post_state.
                    @test haskey(ps, :assignment)
                    @test ps.assignment.state == "indexed"
                    @test sort(collect(Int, ps.assignment.members)) ==
                          HimalayaUI._assignment_body(db, e_id)[:members]
                    @test !isempty(ps.assignment.members)
                    @test !haskey(ps, :assignment_dropped)
                end
            end
        end

        @testset "F-WIPE f: manual analyze drop case carries assignment_dropped" begin
            mktempdir() do tmp
                ctx = _fwipe_setup(tmp)
                db, e_id = ctx.db, ctx.e_id

                before_phases = sort(String[String(r.phase) for r in _fwipe_members(db, e_id)])
                @test !isempty(before_phases)

                # Exclude every auto peak except one — a single effective peak
                # yields zero candidates, so every assigned identity drops.
                auto_qs = Tables.rowtable(DBInterface.execute(db,
                    "SELECT q FROM auto_peaks WHERE exposure_id = ? ORDER BY q", [e_id]))
                for r in auto_qs[2:end]
                    DBInterface.execute(db,
                        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', ?)",
                        [e_id, Float64(r.q)])
                end

                with_test_server(db) do port, base
                    frames = _fwipe_capture_sse("analyze_run") do
                        r = HTTP.post("$base/api/exposures/$e_id/analyze",
                            ["X-Username" => "alice",
                             "X-Client-Op-Id" => "op-fwipe-analyze-2"])
                        @test r.status == 200
                    end
                    @test length(frames) == 1
                    ps = _fwipe_frame_json(frames[1]).post_state
                    @test haskey(ps, :assignment)
                    @test isempty(ps.assignment.members)
                    @test haskey(ps, :assignment_dropped)
                    @test sort(String[String(p) for p in ps.assignment_dropped]) == before_phases
                end
            end
        end

        @testset "F-WIPE f: no-op manual analyze emits no analyze_run frame" begin
            # M0.4 contract: a reanalyze where both findpeaks and indexpeaks are
            # skipped announces nothing — the hashes prove no-op-ness. The manual
            # route's hand-enqueued broadcast must honor the same suppression the
            # _maybe_broadcast_event! path applies (and the assignment-envelope
            # extension must not un-suppress it).
            mktempdir() do tmp
                ctx = _fwipe_setup(tmp)
                db, e_id = ctx.db, ctx.e_id

                with_test_server(db) do port, base
                    frames = _fwipe_capture_sse("analyze_run") do
                        r = HTTP.post("$base/api/exposures/$e_id/analyze",
                            ["X-Username" => "alice",
                             "X-Client-Op-Id" => "op-fwipe-analyze-3"])
                        @test r.status == 200   # response unchanged; only the frame is suppressed
                    end
                    @test isempty(frames)
                end
            end
        end
    end

    @testset "F-WIPE f: non-deferred analyze_exposure! frame carries the assignment envelope" begin
        # The CLI / POST /api/experiments/{id}/analyze paths call
        # analyze_exposure! WITHOUT defer_broadcast — the frame's post_state is
        # built inside the pipeline, not by the route. It must carry the same
        # assignment envelope or a bulk reanalyze strands every tab's cache.
        mktempdir() do tmp
            ctx = _fwipe_setup(tmp)
            db, e_id, dir = ctx.db, ctx.e_id, ctx.analysis_dir

            DBInterface.execute(db,
                "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.123)",
                [e_id])
            # (broadcast_event! never touches current_db() here — the system
            # request has no user, so the actor lookup short-circuits.)
            frames = _fwipe_capture_sse("analyze_run") do
                HimalayaUI.analyze_exposure!(db, e_id, dir)
            end
            @test length(frames) == 1
            ps = _fwipe_frame_json(frames[1]).post_state
            @test haskey(ps, :assignment)
            @test ps.assignment.state == "indexed"
            @test sort(collect(Int, ps.assignment.members)) ==
                  HimalayaUI._assignment_body(db, e_id)[:members]
            @test !haskey(ps, :assignment_dropped)
        end
    end

end
