using Test, HTTP, JSON3, SQLite, DBInterface, Tables
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
function _premint_cmp!(db, id::Int)
    DBInterface.execute(db,
        """INSERT INTO comparisons (id, title, content_hash, created_at, updated_at)
           VALUES (?, '', '', '', '')""", [id])
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
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
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
                    @test p[:sharpness] === nothing
                else
                    @test p[:intensity] isa Real
                end
            end
            @test snap[:analysis_inputs_hash] isa AbstractString
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
end
