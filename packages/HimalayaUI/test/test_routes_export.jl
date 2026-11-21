using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "export routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id,
        name="D1")
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
    for (label, state) in (("Buffer", "confirmed"), ("maybe", "suggested"), ("37C", "confirmed"), ("gone", "removed"))
        DBInterface.execute(db, "INSERT INTO sample_labels (sample_id, label, state) VALUES (?, ?, ?)", [s_id, label, state])
    end

    with_inproc_routes(db) do call
        # JSON
        r = call("GET", "/api/experiments/$exp_id/export?format=json")
        @test r.status == 200
        @test occursin("application/json", HTTP.header(r, "Content-Type"))
        body = JSON3.read(String(r.body))
        @test length(body) == 1
        s = body[1]
        @test s.name == "D1"
        @test length(s.exposures) == 1
        @test s.exposures[1].filename == "example_tot"
        @test length(s.exposures[1].indices) >= 1
        @test collect(s.labels) == ["Buffer", "37C"]
        @test !haskey(s, :tags)

        # CSV
        r = call("GET", "/api/experiments/$exp_id/export?format=csv")
        @test r.status == 200
        @test occursin("text/csv", HTTP.header(r, "Content-Type"))
        csv_body = String(r.body)
        @test startswith(csv_body, "sample_name,exposure_filename,phases")
        @test occursin("D1,example_tot", csv_body)
        @test strip(last(split(first(split(csv_body, "\n")), ","))) == "labels"
        @test endswith(strip(split(csv_body, "\n")[2]), ",Buffer;37C")

        # Default = json
        r = call("GET", "/api/experiments/$exp_id/export")
        @test occursin("application/json", HTTP.header(r, "Content-Type"))

        # Invalid format
        r = call("GET", "/api/experiments/$exp_id/export?format=xml")
        @test r.status == 400
    end
end
