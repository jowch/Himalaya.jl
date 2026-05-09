using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "experiments routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))

    db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db;
        name = "E1", path = tmp,
        data_dir = joinpath(tmp, "data"),
        analysis_dir = analysis_dir)
    s_id = HimalayaUI.create_sample!(db; experiment_id = exp_id,
        name = "D1", display_name = "UX1")
    HimalayaUI.create_exposure!(db; sample_id = s_id, filename = "example_tot")

    with_test_server(db) do port, base
        # GET
        r = HTTP.get("$base/api/experiments/$exp_id")
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test body.id == exp_id
        @test body.name == "E1"

        # PATCH name is no longer allowed — experiments.name is derived from
        # experiment.toml and must change via reingest.
        r = HTTP.patch("$base/api/experiments/$exp_id";
            body = JSON3.write(Dict(:name => "E1-renamed")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            status_exception = false)
        @test r.status == 400

        # POST analyze
        r = HTTP.post("$base/api/experiments/$exp_id/analyze";
            headers = ["X-Username" => "alice"])
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test body.analyzed == 1

        peak_count = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM auto_peaks"))[1].c
        @test peak_count > 0

        # q_units present in response (default A-1 when not set in config; UI prettifies)
        r2 = HTTP.get("$base/api/experiments/$exp_id")
        body2 = JSON3.read(String(r2.body))
        @test haskey(body2, :q_units)
        @test body2.q_units == "A-1"

        # q_units in list response too
        r3 = HTTP.get("$base/api/experiments")
        list = JSON3.read(String(r3.body))
        @test length(list) >= 1
        @test haskey(list[1], :q_units)

        # Malformed config blob must not 500 the route — should fall back to default.
        # This guards against TOML.parse exceptions taking down GET /api/experiments.
        DBInterface.execute(db,
            "UPDATE experiments SET config = ? WHERE id = ?",
            ["[beamline\nq_units = \"x", exp_id])  # unbalanced bracket → invalid TOML
        r4 = HTTP.get("$base/api/experiments/$exp_id")
        @test r4.status == 200
        body4 = JSON3.read(String(r4.body))
        @test body4.q_units == "A-1"
        r5 = HTTP.get("$base/api/experiments")
        @test r5.status == 200

        # 404
        r = HTTP.get("$base/api/experiments/999"; status_exception = false)
        @test r.status == 404

        # PATCH must not touch path fields — those go through reingest.
        original_data_dir = Tables.rowtable(DBInterface.execute(db,
            "SELECT data_dir FROM experiments WHERE id = ?", [exp_id]))[1].data_dir
        r = HTTP.patch("$base/api/experiments/$exp_id";
            body = JSON3.write(Dict(:data_dir => "/somewhere/else",
                                    :analysis_dir => "/elsewhere",
                                    :manifest_path => "/no.csv")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            status_exception = false)
        @test r.status == 400
        # Row must be unchanged.
        @test Tables.rowtable(DBInterface.execute(db,
            "SELECT data_dir FROM experiments WHERE id = ?", [exp_id]))[1].data_dir ==
              original_data_dir
    end
end

@testset "POST /api/experiments/:id/reingest" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    # Drop one .dat file under analysis_dir so reingest discovers an exposure.
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "S001.dat"))
    # Minimal experiment.toml — uses simple.toml defaults; only [experiment].name set.
    write(joinpath(tmp, "experiment.toml"), """
    [experiment]
    name        = "ReE"
    description = ""
    manifest    = "manifest.csv"
    """)
    # Manifest: tab-delimited, 1 header row, columns 1/2/3/9/10/11.
    write(joinpath(tmp, "manifest.csv"),
        "id\tlabel\tname\tc4\tc5\tc6\tc7\tc8\tfile\tnsamp\tnexp\n" *
        "1\tS1\tSample-1\t.\t.\t.\t.\t.\tS001\t\t\n")

    db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db;
        name = "ReE", path = tmp,
        data_dir = joinpath(tmp, "data"),
        analysis_dir = analysis_dir)

    with_test_server(db) do port, base
        # Happy path
        r = HTTP.post("$base/api/experiments/$exp_id/reingest";
            headers = ["X-Username" => "alice"])
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test body.status == "ok"
        @test body.added_samples   >= 1
        @test body.added_exposures >= 1

        n_exp = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM exposures"))[1].c
        @test n_exp >= 1

        # No-manifest path: delete manifest, reingest again — should report
        # status=no_manifest with 0 counts but still succeed (config updated).
        rm(joinpath(tmp, "manifest.csv"))
        r = HTTP.post("$base/api/experiments/$exp_id/reingest";
            headers = ["X-Username" => "alice"])
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test body.status == "no_manifest"
        @test body.added_samples   == 0
        @test body.added_exposures == 0

        # 404 for unknown experiment
        r = HTTP.post("$base/api/experiments/9999/reingest";
            headers = ["X-Username" => "alice"], status_exception = false)
        @test r.status == 404
    end
end

@testset "PATCH /api/experiments/:id no longer accepts name" begin
    tmp = mktempdir()
    db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    eid = HimalayaUI.init_experiment!(db;
        name = "PatchTest", path = tmp,
        data_dir = joinpath(tmp, "data"),
        analysis_dir = joinpath(tmp, "analysis"))

    with_test_server(db) do port, base
        r = HTTP.request("PATCH", "$base/api/experiments/$eid",
            ["Content-Type" => "application/json",
             "X-Username"   => "alice"],
            JSON3.write(Dict(:name => "newname")); status_exception=false)
        @test r.status == 400
    end
end

@testset "POST /api/experiments/:id/reingest returns 400 on validation error" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    # experiment.toml — simple layout.
    write(joinpath(tmp, "experiment.toml"), """
    [experiment]
    name        = "DupTest"
    description = ""
    manifest    = "manifest.csv"
    """)
    # Manifest with duplicate sample names in column 2 (the `name` / stable-id
    # column per simple.toml defaults) — triggers :duplicate_name violation.
    # Column layout: col1=skip_id  col2=name  col3=display_name  col9=filenames
    write(joinpath(tmp, "manifest.csv"),
        "skip_header_row\n" *
        "1\tDupSample\tLabel-A\t.\t.\t.\t.\t.\tS001\t\t\n" *
        "2\tDupSample\tLabel-B\t.\t.\t.\t.\t.\tS002\t\t\n")

    db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    eid = HimalayaUI.init_experiment!(db;
        name = "DupTest", path = tmp,
        data_dir = joinpath(tmp, "data"),
        analysis_dir = analysis_dir)

    with_test_server(db) do port, base
        r = HTTP.post("$base/api/experiments/$eid/reingest";
            headers = ["X-Username" => "alice"], status_exception = false)
        @test r.status == 400
        body = JSON3.read(String(r.body))
        @test body[:error] == "manifest_invalid"
        @test !isempty(body[:violations])
    end
end
