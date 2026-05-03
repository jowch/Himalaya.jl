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
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")

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
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")

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
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")

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
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
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
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
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

    @testset "POST /api/groups/:id/members → index_confirmed frame" begin
        tmp = mktempdir()
        analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
        mkpath(analysis_dir)
        cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
           joinpath(analysis_dir, "example_tot.dat"))
        db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        exp_id = HimalayaUI.init_experiment!(db; path=tmp,
            data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
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
                idx = first(filter(i -> !(i.id in groups[1].members), indices))

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
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
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
