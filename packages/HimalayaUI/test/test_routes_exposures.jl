using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "exposures routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))

    db     = HimalayaUI.open_db(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id,
        label="D1", name="UX1")
    e1 = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    e2 = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="other")

    with_test_server(db) do port, base
        # List
        r = HTTP.get("$base/api/samples/$s_id/exposures")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test length(list) == 2
        @test list[1].filename == "example_tot"
        @test list[1].selected === false
        @test list[1].tags     == []
        @test list[1].sources  == []

        # Select e2
        r = HTTP.patch("$base/api/exposures/$e2/select";
            headers = ["X-Username" => "alice"])
        @test r.status == 200

        r    = HTTP.get("$base/api/samples/$s_id/exposures")
        list = JSON3.read(String(r.body))
        sel  = findfirst(x -> x.id == e2, list)
        @test list[sel].selected === true
        other = findfirst(x -> x.id == e1, list)
        @test list[other].selected === false

        # Add tag to e1
        r = HTTP.post("$base/api/exposures/$e1/tags";
            body = JSON3.write(Dict(:key=>"flag", :value=>"good")),
            headers = ["Content-Type"=>"application/json", "X-Username"=>"alice"])
        @test r.status == 201
        tag_id = JSON3.read(String(r.body)).id

        r    = HTTP.get("$base/api/samples/$s_id/exposures")
        list = JSON3.read(String(r.body))
        other = findfirst(x -> x.id == e1, list)
        @test length(list[other].tags) == 1

        # Delete tag
        r = HTTP.delete("$base/api/exposures/$e1/tags/$tag_id";
            headers = ["X-Username" => "alice"])
        @test r.status == 204

        # Analyze single exposure
        r = HTTP.post("$base/api/exposures/$e1/analyze";
            headers = ["X-Username" => "alice"])
        @test r.status == 200
        @test JSON3.read(String(r.body)).analyzed === true

        # 404 for nonexistent
        r = HTTP.post("$base/api/exposures/9999/analyze";
            headers = ["X-Username" => "alice"], status_exception = false)
        @test r.status == 404
    end
end
