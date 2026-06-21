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

        # PATCH name is now allowed (Phase E1: name became editable).
        r = call("PATCH", "/api/experiments/$exp_id";
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:name => "E1-renamed"))))
        @test r.status == 200

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

@testset "PATCH /api/experiments/:id accepts name (Phase E1)" begin
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
        @test r.status == 200
        rb = call("GET", "/api/experiments/$eid")
        @test JSON3.read(String(rb.body)).name == "newname"
    end
end

@testset "experiment stats: span_hours + derived sessions + list parity" begin
    tmp = mktempdir()
    db = open_prepared_clone(tmp)
    eid = HimalayaUI.init_experiment!(db;
        name = "Stats", path = tmp,
        data_dir = joinpath(tmp, "data"),
        analysis_dir = joinpath(tmp, "analysis"))
    # Two loads 3h apart → 2 macro-sessions. session_id persisted at ingest.
    l1 = HimalayaUI.create_load!(db; experiment_id = eid, load_index = 1,
        start_time = "2026-04-25T10:00:00", session_id = 1)
    l2 = HimalayaUI.create_load!(db; experiment_id = eid, load_index = 2,
        start_time = "2026-04-25T13:00:00", session_id = 2)
    s = HimalayaUI.create_sample!(db; experiment_id = eid, name = "D1")
    HimalayaUI.create_exposure!(db; experiment_id = eid, sample_id = s, load_id = l1,
        filename = "a", timestamp = "2026-04-25T10:00:00")
    HimalayaUI.create_exposure!(db; experiment_id = eid, sample_id = s, load_id = l2,
        filename = "b", timestamp = "2026-04-25T13:00:00")

    with_inproc_routes(db) do call
        # Detail endpoint: stats carry span_hours + derived sessions.
        body = JSON3.read(String(call("GET", "/api/experiments/$eid").body))
        @test body.stats.loads == 2
        @test body.stats.exposures == 2
        @test body.stats.sessions == 2                       # 3h gap > 2h threshold
        @test isapprox(Float64(body.stats.span_hours), 3.0; atol = 0.05)

        # List endpoint now carries the same stats (was the gap: list omitted stats).
        list = JSON3.read(String(call("GET", "/api/experiments").body))
        row = first(filter(e -> e.id == eid, list))
        @test haskey(row, :stats)
        @test row.stats.loads == 2
        @test row.stats.sessions == 2
        @test isapprox(Float64(row.stats.span_hours), 3.0; atol = 0.05)
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

@testset "GET /experiments/:id stats includes sessions count" begin
    tmp = mktempdir()
    db = open_prepared_clone(tmp)
    exp_id = HimalayaUI.create_experiment!(db; path=tmp, data_dir="data", analysis_dir="analysis")

    # Insert 3 loads across 2 distinct sessions (one load has session_id=nothing).
    HimalayaUI.create_load!(db; experiment_id=exp_id, load_index=1, session_id=10)
    HimalayaUI.create_load!(db; experiment_id=exp_id, load_index=2, session_id=10)
    HimalayaUI.create_load!(db; experiment_id=exp_id, load_index=3, session_id=20)
    HimalayaUI.create_load!(db; experiment_id=exp_id, load_index=4, session_id=nothing)

    with_inproc_routes(db) do call
        r = call("GET", "/api/experiments/$exp_id")
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test haskey(body, :stats)
        @test body.stats.loads == 4
        # sessions = COUNT(DISTINCT session_id) WHERE session_id IS NOT NULL = 2 (10 and 20)
        @test body.stats.sessions == 2
    end
end

@testset "experiment stats: started_at = min exposure timestamp" begin
    tmp = mktempdir()
    db = open_prepared_clone(tmp)
    eid = HimalayaUI.init_experiment!(db;
        name = "StartedAt", path = tmp,
        data_dir = joinpath(tmp, "data"),
        analysis_dir = joinpath(tmp, "analysis"))
    l = HimalayaUI.create_load!(db; experiment_id = eid, load_index = 1)
    s = HimalayaUI.create_sample!(db; experiment_id = eid, name = "S1", load_id = l)
    # Two exposures; the earlier timestamp should be the started_at.
    HimalayaUI.create_exposure!(db; experiment_id = eid, sample_id = s, load_id = l,
        filename = "a", timestamp = "2026-03-01T08:00:00", frame_no = 1)
    HimalayaUI.create_exposure!(db; experiment_id = eid, sample_id = s, load_id = l,
        filename = "b", timestamp = "2026-04-15T12:00:00", frame_no = 2)

    with_inproc_routes(db) do call
        # Detail endpoint: stats carries started_at.
        body = JSON3.read(String(call("GET", "/api/experiments/$eid").body))
        @test haskey(body.stats, :started_at)
        @test startswith(String(body.stats.started_at), "2026-03-01")

        # List endpoint: same started_at propagates.
        list = JSON3.read(String(call("GET", "/api/experiments").body))
        row = first(filter(e -> e.id == eid, list))
        @test haskey(row.stats, :started_at)
        @test startswith(String(row.stats.started_at), "2026-03-01")
    end
end

@testset "experiment list: review_count" begin
    tmp = mktempdir()
    db = open_prepared_clone(tmp)
    eid = HimalayaUI.init_experiment!(db;
        name = "ReviewCount", path = tmp,
        data_dir = joinpath(tmp, "data"),
        analysis_dir = joinpath(tmp, "analysis"))

    with_inproc_routes(db) do call
        # CLEAN experiment: one load, one sample, two exposures at the SAME
        # horizontal_position → no split flag → review_count == 0.
        l = HimalayaUI.create_load!(db; experiment_id = eid, load_index = 1)
        s = HimalayaUI.create_sample!(db; experiment_id = eid, name = "Clean",
            load_id = l, slot_index = 1)
        HimalayaUI.create_exposure!(db; experiment_id = eid, sample_id = s, load_id = l,
            filename = "c1", horizontal_position = 10.0, frame_no = 1)
        HimalayaUI.create_exposure!(db; experiment_id = eid, sample_id = s, load_id = l,
            filename = "c2", horizontal_position = 10.1, frame_no = 2)  # ~0.1 mm jitter, no flag

        list = JSON3.read(String(call("GET", "/api/experiments").body))
        row = first(filter(e -> e.id == eid, list))
        @test haskey(row, :review_count)
        @test row.review_count == 0

        # FLAGGED experiment: second experiment, one sample, two exposures with
        # horizontal_position jump > 0.5 mm → SplitFlag → review_count >= 1.
        eid2 = HimalayaUI.init_experiment!(db;
            name = "Flagged", path = mktempdir(),
            data_dir = joinpath(mktempdir(), "data"),
            analysis_dir = joinpath(mktempdir(), "analysis"))
        l2 = HimalayaUI.create_load!(db; experiment_id = eid2, load_index = 1)
        s2 = HimalayaUI.create_sample!(db; experiment_id = eid2, name = "Split",
            load_id = l2, slot_index = 1)
        HimalayaUI.create_exposure!(db; experiment_id = eid2, sample_id = s2, load_id = l2,
            filename = "f1", horizontal_position = 10.0, frame_no = 1)
        HimalayaUI.create_exposure!(db; experiment_id = eid2, sample_id = s2, load_id = l2,
            filename = "f2", horizontal_position = 20.0, frame_no = 2)  # 10 mm jump → SplitFlag

        list2 = JSON3.read(String(call("GET", "/api/experiments").body))
        row2 = first(filter(e -> e.id == eid2, list2))
        @test haskey(row2, :review_count)
        @test row2.review_count >= 1
    end
end

@testset "POST /api/experiments rejects a duplicate data_dir (409)" begin
    tmp = mktempdir()
    data_dir = joinpath(tmp, "data"); mkpath(data_dir)
    db = open_prepared_clone(tmp)
    # An experiment already uses this directory.
    HimalayaUI.init_experiment!(db; name = "E1", path = tmp,
        data_dir = data_dir, analysis_dir = data_dir)

    with_inproc_routes(db) do call
        r = call("POST", "/api/experiments";
            headers = ["Content-Type" => "application/json"],
            body = Vector{UInt8}(JSON3.write(Dict(:path => data_dir))))
        @test r.status == 409
        body = JSON3.read(String(r.body))
        @test occursin("director", lowercase(String(body.error)))
        # No second experiment was created for that directory.
        n = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM experiments WHERE data_dir = ?", [data_dir]))[1].c
        @test n == 1
    end
end
