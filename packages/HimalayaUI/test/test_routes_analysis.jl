using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "indices + groups routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    with_test_server(db) do port, base
        # Indices
        r = HTTP.get("$base/api/exposures/$e_id/indices")
        @test r.status == 200
        indices = JSON3.read(String(r.body))
        @test length(indices) >= 1
        @test haskey(indices[1], :peaks)

        # Groups — auto only, active
        r = HTTP.get("$base/api/exposures/$e_id/groups")
        @test r.status == 200
        groups = JSON3.read(String(r.body))
        @test length(groups) == 1
        @test groups[1].kind   == "auto"
        @test groups[1].active === true
        auto_gid  = groups[1].id
        auto_mems = groups[1].members

        # Find a candidate index NOT in the auto group
        extra_candidate = nothing
        for ix in indices
            ix.id in auto_mems && continue
            extra_candidate = ix.id
            break
        end

        if extra_candidate !== nothing
            r = HTTP.post("$base/api/groups/$auto_gid/members";
                body = JSON3.write(Dict(:index_id => extra_candidate)),
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice"])
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test body.kind   == "custom"
            @test body.active === true
            @test extra_candidate in body.members

            r = HTTP.get("$base/api/exposures/$e_id/groups")
            groups = JSON3.read(String(r.body))
            @test length(groups) == 2
            auto_g  = first(filter(g -> g.kind == "auto",   groups))
            cust_g  = first(filter(g -> g.kind == "custom", groups))
            @test auto_g.active === false
            @test cust_g.active === true
            @test extra_candidate in cust_g.members
        end

        # Reset groups and re-run analyze for a clean DELETE test
        DBInterface.execute(db,
            "DELETE FROM index_group_members WHERE group_id IN
             (SELECT id FROM index_groups WHERE exposure_id = ?)", [e_id])
        DBInterface.execute(db,
            "DELETE FROM index_groups WHERE exposure_id = ?", [e_id])
        HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

        r = HTTP.get("$base/api/exposures/$e_id/groups")
        groups = JSON3.read(String(r.body))
        @test length(groups) == 1
        auto_gid  = groups[1].id
        auto_mems = groups[1].members
        if !isempty(auto_mems)
            removed = first(auto_mems)
            r = HTTP.delete("$base/api/groups/$auto_gid/members/$removed";
                headers = ["X-Username" => "alice"])
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test body.kind == "custom"
            @test !(removed in body.members)
        end

        # 404 on unknown group
        r = HTTP.post("$base/api/groups/99999/members";
            body = JSON3.write(Dict(:index_id => 1)),
            headers = ["Content-Type" => "application/json"],
            status_exception = false)
        @test r.status == 404
    end
end
