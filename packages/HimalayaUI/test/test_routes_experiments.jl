using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "experiments routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))

    db = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db;
        name = "E1", path = tmp,
        data_dir = joinpath(tmp, "data"),
        analysis_dir = analysis_dir)
    s_id = HimalayaUI.create_sample!(db; experiment_id = exp_id,
        name = "D1")
    HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = s_id, filename = "example_tot")

    with_inproc_routes(db) do call
        # GET
        r = call("GET", "/api/experiments/$exp_id")
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test body.id == exp_id
        @test body.name == "E1"

        # PATCH name is no longer allowed — experiments.name is derived from
        # experiment.toml and must change via reingest.
        r = call("PATCH", "/api/experiments/$exp_id";
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:name => "E1-renamed"))))
        @test r.status == 400

        # POST analyze
        r = call("POST", "/api/experiments/$exp_id/analyze";
            headers = ["X-Username" => "alice"])
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test body.analyzed == 1

        peak_count = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM auto_peaks"))[1].c
        @test peak_count > 0

        # q_units present in response (default A-1 when not set in config; UI prettifies)
        r2 = call("GET", "/api/experiments/$exp_id")
        body2 = JSON3.read(String(r2.body))
        @test haskey(body2, :q_units)
        @test body2.q_units == "A-1"

        # q_units in list response too
        r3 = call("GET", "/api/experiments")
        list = JSON3.read(String(r3.body))
        @test length(list) >= 1
        @test haskey(list[1], :q_units)

        # Malformed config blob must not 500 the route — should fall back to default.
        # This guards against TOML.parse exceptions taking down GET /api/experiments.
        DBInterface.execute(db,
            "UPDATE experiments SET config = ? WHERE id = ?",
            ["[beamline\nq_units = \"x", exp_id])  # unbalanced bracket → invalid TOML
        r4 = call("GET", "/api/experiments/$exp_id")
        @test r4.status == 200
        body4 = JSON3.read(String(r4.body))
        @test body4.q_units == "A-1"
        r5 = call("GET", "/api/experiments")
        @test r5.status == 200

        # 404
        r = call("GET", "/api/experiments/999")
        @test r.status == 404

        # PATCH must not touch path fields — those go through reingest.
        original_data_dir = Tables.rowtable(DBInterface.execute(db,
            "SELECT data_dir FROM experiments WHERE id = ?", [exp_id]))[1].data_dir
        r = call("PATCH", "/api/experiments/$exp_id";
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:data_dir => "/somewhere/else",
                                    :analysis_dir => "/elsewhere",
                                    :manifest_path => "/no.csv"))))
        @test r.status == 400
        # Row must be unchanged.
        @test Tables.rowtable(DBInterface.execute(db,
            "SELECT data_dir FROM experiments WHERE id = ?", [exp_id]))[1].data_dir ==
              original_data_dir
    end
end

@testset "PATCH /api/experiments/:id no longer accepts name" begin
    tmp = mktempdir()
    db = open_prepared_clone(tmp)
    eid = HimalayaUI.init_experiment!(db;
        name = "PatchTest", path = tmp,
        data_dir = joinpath(tmp, "data"),
        analysis_dir = joinpath(tmp, "analysis"))

    with_inproc_routes(db) do call
        r = call("PATCH", "/api/experiments/$eid";
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:name => "newname"))))
        @test r.status == 400
    end
end

@testset "experiment route surfaces beam center + pixel size" begin
    tmp = mktempdir()
    db = open_prepared_clone(tmp)
    exp_id = HimalayaUI.create_experiment!(db; path=tmp, data_dir="data", analysis_dir="analysis")

    DBInterface.execute(db, "UPDATE experiments SET config = ? WHERE id = ?", [
        """
        [beamline]
        beam_center_x = 420.791
        beam_center_y = 838.83
        pixel_size_um = 172
        """, exp_id])

    with_inproc_routes(db) do call
        body = JSON3.read(String(call("GET", "/api/experiments/$exp_id").body))
        @test body.beam_center_x == 420.791
        @test body.beam_center_y == 838.83
        @test body.pixel_size_um == 172.0   # bare int coerced to Float64

        # Garbage numeric value must not 500 — on the single route OR the list
        # route (the spec's primary protected case: one bad config in a loop
        # must not break the whole list).
        DBInterface.execute(db, "UPDATE experiments SET config = ? WHERE id = ?",
            ["[beamline]\nbeam_center_x = \"oops\"\n", exp_id])
        r = call("GET", "/api/experiments/$exp_id")
        @test r.status == 200
        @test JSON3.read(String(r.body)).beam_center_x === nothing

        rl = call("GET", "/api/experiments")
        @test rl.status == 200
        @test JSON3.read(String(rl.body))[1].beam_center_x === nothing

        # A non-string q_units must not 500 the list route either (the shim's
        # ::String contract is upheld by coercing to the default).
        DBInterface.execute(db, "UPDATE experiments SET config = ? WHERE id = ?",
            ["[beamline]\nq_units = 5\n", exp_id])
        rq = call("GET", "/api/experiments")
        @test rq.status == 200
        @test JSON3.read(String(rq.body))[1].q_units == "A-1"
    end
end
