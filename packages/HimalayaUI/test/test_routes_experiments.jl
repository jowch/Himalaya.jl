using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "experiments routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))

    db = HimalayaUI.open_db(tmp)
    exp_id = HimalayaUI.init_experiment!(db;
        name = "E1", path = tmp,
        data_dir = joinpath(tmp, "data"),
        analysis_dir = analysis_dir)
    s_id = HimalayaUI.create_sample!(db; experiment_id = exp_id,
        label = "D1", name = "UX1")
    HimalayaUI.create_exposure!(db; sample_id = s_id, filename = "example_tot")

    with_test_server(db) do port, base
        # GET
        r = HTTP.get("$base/api/experiments/$exp_id")
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test body.id == exp_id
        @test body.name == "E1"

        # PATCH name
        r = HTTP.patch("$base/api/experiments/$exp_id";
            body = JSON3.write(Dict(:name => "E1-renamed")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 200
        @test JSON3.read(String(r.body)).name == "E1-renamed"

        # POST analyze
        r = HTTP.post("$base/api/experiments/$exp_id/analyze";
            headers = ["X-Username" => "alice"])
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test body.analyzed == 1

        peak_count = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM peaks"))[1].c
        @test peak_count > 0

        # 404
        r = HTTP.get("$base/api/experiments/999"; status_exception = false)
        @test r.status == 404
    end
end
