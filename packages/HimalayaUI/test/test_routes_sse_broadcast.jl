using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# Regression test for PR review issue #5: every migrated route that calls
# `apply_event!(InTransaction(), ...)` must explicitly enqueue an SSE broadcast
# (the InTransaction variant does NOT broadcast — caller is responsible).
# Without this, multiplayer is silently broken for these event kinds.
#
# Each subtest:
# 1. Hooks a fake subscriber into `SSE_SUBSCRIBERS`.
# 2. Hits the route via `with_test_server`.
# 3. Asserts the channel saw a frame whose `kind` matches the route's event kind.

function _drain_frames(ch::Channel{String}; timeout_s = 1.0)
    frames = String[]
    deadline = time() + timeout_s
    while time() < deadline
        if isready(ch)
            push!(frames, take!(ch))
        else
            sleep(0.02)
        end
    end
    frames
end

function _frame_kinds(frames::Vector{String})
    out = String[]
    for f in frames
        m = match(r"\"kind\":\"([a-z_]+)\"", f)
        m === nothing || push!(out, m.captures[1])
    end
    out
end

function _with_sub(f::Function)
    pending = Channel{String}(64)
    sub = (pending = pending,)
    lock(HimalayaUI.SSE_LOCK) do
        push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
    end
    try
        f(pending)
    finally
        lock(HimalayaUI.SSE_LOCK) do
            filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
        end
        close(pending)
        HimalayaUI.SSE_SUBSCRIBERS[] = []
    end
end

@testset "Migrated routes broadcast SSE frames (issue #5)" begin
    @testset "POST /api/samples/:id/messages → post_message frame" begin
        tmp    = mktempdir()
        db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path=tmp,
            data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")

        _with_sub() do pending
            with_test_server(db) do port, base
                r = HTTP.post("$base/api/samples/$s_id/messages";
                    body = JSON3.write(Dict(:body => "hello")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r.status == 201
            end
            kinds = _frame_kinds(_drain_frames(pending))
            @test "post_message" in kinds
        end
    end

    @testset "PATCH /api/samples/:id → update_sample frame" begin
        tmp    = mktempdir()
        db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path=tmp,
            data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")

        _with_sub() do pending
            with_test_server(db) do port, base
                r = HTTP.patch("$base/api/samples/$s_id";
                    body = JSON3.write(Dict(:notes => "n")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r.status == 200
            end
            kinds = _frame_kinds(_drain_frames(pending))
            @test "update_sample" in kinds
        end
    end

    @testset "POST/DELETE /api/samples/:id/tags → add_tag/remove_tag frames" begin
        tmp    = mktempdir()
        db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path=tmp,
            data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")

        _with_sub() do pending
            with_test_server(db) do port, base
                r = HTTP.post("$base/api/samples/$s_id/tags";
                    body = JSON3.write(Dict(:key => "k", :value => "v")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r.status == 201
                tag_id = JSON3.read(String(r.body)).id

                r2 = HTTP.delete("$base/api/samples/$s_id/tags/$tag_id";
                    headers = ["X-Username" => "alice"])
                @test r2.status == 204
            end
            kinds = _frame_kinds(_drain_frames(pending))
            @test "add_tag" in kinds
            @test "remove_tag" in kinds
        end
    end

    @testset "PATCH /api/exposures/:id/status → set_exposure_status frame" begin
        tmp    = mktempdir()
        db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path=tmp,
            data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="x")

        _with_sub() do pending
            with_test_server(db) do port, base
                r = HTTP.patch("$base/api/exposures/$e_id/status";
                    body = JSON3.write(Dict(:status => "accepted")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r.status == 200
            end
            kinds = _frame_kinds(_drain_frames(pending))
            @test "set_exposure_status" in kinds
        end
    end

    @testset "PATCH /api/exposures/:id/select → select_exposure frame" begin
        tmp    = mktempdir()
        db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path=tmp,
            data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="x")

        _with_sub() do pending
            with_test_server(db) do port, base
                r = HTTP.patch("$base/api/exposures/$e_id/select";
                    headers = ["X-Username" => "alice"])
                @test r.status == 200
            end
            kinds = _frame_kinds(_drain_frames(pending))
            @test "select_exposure" in kinds
        end
    end

    @testset "DELETE /api/groups/:id/members/:idx → index_unconfirmed frame" begin
        tmp = mktempdir()
        analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
        mkpath(analysis_dir)
        cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
           joinpath(analysis_dir, "example_tot.dat"))
        db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        exp_id = HimalayaUI.init_experiment!(db; path=tmp,
            data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
        HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

        _with_sub() do pending
            with_test_server(db) do port, base
                r = HTTP.get("$base/api/exposures/$e_id/groups")
                groups = JSON3.read(String(r.body))
                gid = groups[1].id
                isempty(groups[1].members) && return
                idx_id = first(groups[1].members)
                r2 = HTTP.delete("$base/api/groups/$gid/members/$idx_id";
                    headers = ["X-Username" => "alice"])
                @test r2.status == 200
            end
            kinds = _frame_kinds(_drain_frames(pending))
            @test "index_unconfirmed" in kinds
        end
    end

    @testset "POST /api/exposures/:id/speculative → speculative_created frame" begin
        tmp = mktempdir()
        analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
        mkpath(analysis_dir)
        cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
           joinpath(analysis_dir, "example_tot.dat"))
        db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        exp_id = HimalayaUI.init_experiment!(db; path=tmp,
            data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
        HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

        peaks = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 2", [e_id]))
        p1 = Int(peaks[1].id); p2 = Int(peaks[2].id)

        _with_sub() do pending
            with_test_server(db) do port, base
                body = Dict(:phase => "Lamellar",
                            :anchor_peak_id => p1, :anchor_ratio => 1,
                            :additional => [Dict(:ratio_position => 2, :peak_id => p2)])
                r = HTTP.post("$base/api/exposures/$e_id/speculative";
                    body = JSON3.write(body),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r.status == 200
                idx_id = JSON3.read(String(r.body)).id

                # Also exercise the DELETE path while we have a fresh index
                r2 = HTTP.delete("$base/api/indices/$idx_id";
                    headers = ["X-Username" => "alice"])
                @test r2.status == 200
            end
            kinds = _frame_kinds(_drain_frames(pending))
            @test "speculative_created" in kinds
            @test "speculative_deleted" in kinds
        end
    end

    @testset "POST /api/exposures/:id/analyze → analyze_run frame with post_state" begin
        tmp = mktempdir()
        analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
        mkpath(analysis_dir)
        cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
           joinpath(analysis_dir, "example_tot.dat"))
        db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        exp_id = HimalayaUI.init_experiment!(db; path=tmp,
            data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
        HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

        # Force a non-no-op reanalyze: drop the inputs hash so it falls through
        # to the slow path and emits a fresh analyze_run row.
        DBInterface.execute(db,
            "UPDATE exposures SET analysis_inputs_hash = NULL WHERE id = ?", [e_id])

        _with_sub() do pending
            with_test_server(db) do port, base
                r = HTTP.post("$base/api/exposures/$e_id/analyze";
                    headers = ["X-Username"     => "alice",
                               "X-Client-Op-Id" => "uuid-analyze-broadcast"])
                @test r.status == 200
            end
            frames = _drain_frames(pending)
            kinds  = _frame_kinds(frames)
            @test "analyze_run" in kinds
            # Spec contract: analyze_run frames carry post_state
            # (analysis_inputs_hash + indices). Verify presence on the frame.
            analyze_frame = first(filter(f -> occursin("\"kind\":\"analyze_run\"", f), frames))
            @test occursin("\"post_state\"", analyze_frame)
            @test occursin("\"analysis_inputs_hash\"", analyze_frame)
            @test occursin("\"indices\"", analyze_frame)
        end
    end

    @testset "POST /api/groups/:id/members → index_confirmed frame" begin
        tmp = mktempdir()
        analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
        mkpath(analysis_dir)
        cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
           joinpath(analysis_dir, "example_tot.dat"))
        db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        exp_id = HimalayaUI.init_experiment!(db; path=tmp,
            data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
        HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

        _with_sub() do pending
            with_test_server(db) do port, base
                r = HTTP.get("$base/api/exposures/$e_id/groups")
                @test r.status == 200
                groups = JSON3.read(String(r.body))
                @test length(groups) >= 1
                gid = groups[1].id

                r2 = HTTP.get("$base/api/exposures/$e_id/indices")
                indices = JSON3.read(String(r2.body))
                # The route adds to the custom group, independent of auto-group
                # membership, so any candidate index exercises index_confirmed.
                idx = first(indices)

                r3 = HTTP.post("$base/api/groups/$gid/members";
                    body = JSON3.write(Dict(:index_id => idx.id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r3.status == 200
            end
            kinds = _frame_kinds(_drain_frames(pending))
            @test "index_confirmed" in kinds
        end
    end

    @testset "POST/DELETE /api/exposures/:id/tags → add_tag/remove_tag frames" begin
        tmp    = mktempdir()
        db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path=tmp,
            data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="x")

        _with_sub() do pending
            with_test_server(db) do port, base
                r = HTTP.post("$base/api/exposures/$e_id/tags";
                    body = JSON3.write(Dict(:key => "k", :value => "v")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r.status == 201
                tag_id = JSON3.read(String(r.body)).id

                r2 = HTTP.delete("$base/api/exposures/$e_id/tags/$tag_id";
                    headers = ["X-Username" => "alice"])
                @test r2.status == 204
            end
            kinds = _frame_kinds(_drain_frames(pending))
            @test "add_tag" in kinds
            @test "remove_tag" in kinds
        end
    end
end
