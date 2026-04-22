using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "peaks routes" begin
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
        # GET peaks
        r = HTTP.get("$base/api/exposures/$e_id/peaks")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test length(list) > 0
        @test all(p -> p.source == "auto", list)

        initial_indices = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM indices WHERE exposure_id = ?", [e_id]))
        @test !isempty(initial_indices)

        # POST manual peak
        r = HTTP.post("$base/api/exposures/$e_id/peaks";
            body = JSON3.write(Dict(:q => 0.5)),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 201
        body = JSON3.read(String(r.body))
        @test body.q == 0.5
        @test body.source == "manual"
        peak_id = body.id
        @test body.stale_indices == length(initial_indices)

        stale = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM indices
             WHERE exposure_id = ? AND status = 'stale'", [e_id]))
        @test stale[1].c == length(initial_indices)

        r = HTTP.get("$base/api/exposures/$e_id/peaks")
        list = JSON3.read(String(r.body))
        manual = filter(p -> p.source == "manual", list)
        @test length(manual) == 1

        # DELETE peak
        r = HTTP.delete("$base/api/peaks/$peak_id";
            headers = ["X-Username" => "alice"])
        @test r.status == 204

        r = HTTP.get("$base/api/exposures/$e_id/peaks")
        list = JSON3.read(String(r.body))
        @test !any(p -> p.id == peak_id, list)

        # 404 delete
        r = HTTP.delete("$base/api/peaks/99999";
            headers = ["X-Username" => "alice"], status_exception = false)
        @test r.status == 404
    end
end
