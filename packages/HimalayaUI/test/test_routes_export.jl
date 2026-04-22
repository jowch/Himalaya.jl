using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "export routes" begin
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
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    with_test_server(db) do port, base
        # JSON
        r = HTTP.get("$base/api/experiments/$exp_id/export?format=json")
        @test r.status == 200
        @test occursin("application/json", HTTP.header(r, "Content-Type"))
        body = JSON3.read(String(r.body))
        @test length(body) == 1
        s = body[1]
        @test s.label == "D1"
        @test length(s.exposures) == 1
        @test s.exposures[1].filename == "example_tot"
        @test length(s.exposures[1].indices) >= 1

        # CSV
        r = HTTP.get("$base/api/experiments/$exp_id/export?format=csv")
        @test r.status == 200
        @test occursin("text/csv", HTTP.header(r, "Content-Type"))
        csv_body = String(r.body)
        @test startswith(csv_body, "sample_label,sample_name,exposure_filename,phases")
        @test occursin("D1,UX1,example_tot", csv_body)

        # Default = json
        r = HTTP.get("$base/api/experiments/$exp_id/export")
        @test occursin("application/json", HTTP.header(r, "Content-Type"))

        # Invalid format
        r = HTTP.get("$base/api/experiments/$exp_id/export?format=xml";
                     status_exception = false)
        @test r.status == 400
    end
end
