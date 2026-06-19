using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "exposures routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))

    db     = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id,
        name="D1", display_name="UX1")
    e1 = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    e2 = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="other")

    with_inproc_routes(db) do call
        # List
        r = call("GET", "/api/samples/$s_id/exposures")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test length(list) == 2
        @test list[1].filename == "example_tot"
        @test list[1].selected === false
        @test list[1].tags     == []
        @test list[1].sources  == []

        # Select e2
        r = call("PATCH", "/api/exposures/$e2/select";
            headers = ["X-Username" => "alice"])
        @test r.status == 200

        r    = call("GET", "/api/samples/$s_id/exposures")
        list = JSON3.read(String(r.body))
        sel  = findfirst(x -> x.id == e2, list)
        @test list[sel].selected === true
        other = findfirst(x -> x.id == e1, list)
        @test list[other].selected === false

        # Add tag to e1
        r = call("POST", "/api/exposures/$e1/tags";
            headers = ["Content-Type"=>"application/json", "X-Username"=>"alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:key=>"flag", :value=>"good"))))
        @test r.status == 201
        tag_id = JSON3.read(String(r.body)).id

        r    = call("GET", "/api/samples/$s_id/exposures")
        list = JSON3.read(String(r.body))
        other = findfirst(x -> x.id == e1, list)
        @test length(list[other].tags) == 1

        # Delete tag
        r = call("DELETE", "/api/exposures/$e1/tags/$tag_id";
            headers = ["X-Username" => "alice"])
        @test r.status == 204

        # Analyze single exposure
        r = call("POST", "/api/exposures/$e1/analyze";
            headers = ["X-Username" => "alice"])
        @test r.status == 200
        @test JSON3.read(String(r.body)).analyzed === true

        # 404 for nonexistent
        r = call("POST", "/api/exposures/9999/analyze";
            headers = ["X-Username" => "alice"])
        @test r.status == 404
    end
end

@testset "PATCH /api/exposures/:id/select is idempotent under retry" begin
    tmp = mktempdir()
    db     = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e1 = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="A")
    e2 = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="B")

    with_inproc_routes(db) do call
        op_id = "uuid-select-retry-1"
        hdrs  = ["X-Username"     => "alice",
                 "X-Client-Op-Id" => op_id]

        r1 = call("PATCH", "/api/exposures/$e2/select"; headers=hdrs)
        @test r1.status == 200
        body1 = String(copy(r1.body))

        r2 = call("PATCH", "/api/exposures/$e2/select"; headers=hdrs)
        @test r2.status == 200
        @test String(copy(r2.body)) == body1

        # Exactly one selected, and it's e2
        n_sel = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM exposures WHERE sample_id = ? AND selected = 1", [s_id]))).c
        @test n_sel == 1

        n_events = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions WHERE action = 'select_exposure' AND entity_id = ?",
            [e2]))).c
        @test n_events == 1
    end
end
