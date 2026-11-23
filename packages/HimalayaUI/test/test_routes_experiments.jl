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

        for sm in HimalayaUI.get_samples(db, exp_id), ex in HimalayaUI.get_exposures(db, Int(sm.id))
            HimalayaUI.analyze_exposure!(db, Int(ex.id))
        end

        peak_count = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM auto_peaks"))[1].c
        @test peak_count > 0

        # q_units present in response (default A-1 when not set in config; UI prettifies)
        r2 = call("GET", "/api/experiments/$exp_id")
        body2 = JSON3.read(String(r2.body))
        @test length(body2.setups) == 1
        @test body2.setups[1].q_units == "A-1"
        @test body2.setups[1].name == "$(basename(tmp))"
        @test body2.setups[1].included == true

        # q_units in list response too
        r3 = call("GET", "/api/experiments")
        list = JSON3.read(String(r3.body))
        @test length(list) >= 1
        @test haskey(list[1].setups[1], :q_units)

        # Typed-geometry row with a NULL q_units (a migrated/fresh experiment that
        # carries real beam center etc.) must STILL default q_units to A-1. The
        # typed path used to serve the raw NULL column, blanking the geometry ledger.
        DBInterface.execute(db,
            "UPDATE setups SET beam_center_x = 421.4, beam_center_x_source = 'setup', q_units = NULL WHERE experiment_id = ?",
            [exp_id])
        r4 = call("GET", "/api/experiments/$exp_id")
        body4 = JSON3.read(String(r4.body))
        @test body4.setups[1].beam_center_x == 421.4
        @test body4.setups[1].q_units == "A-1"
        @test body4.setups[1].q_units_source == "default"

        r5 = call("GET", "/api/experiments")
        @test r5.status == 200

        # 404
        r = call("GET", "/api/experiments/999")
        @test r.status == 404

        # PATCH must not touch path fields — those go through reingest.
        original_data_dir = Tables.rowtable(DBInterface.execute(db,
            "SELECT data_dir FROM setups WHERE experiment_id = ?", [exp_id]))[1].data_dir
        r = call("PATCH", "/api/experiments/$exp_id";
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:data_dir => "/somewhere/else",
                                    :analysis_dir => "/elsewhere",
                                    :manifest_path => "/no.csv"))))
        @test r.status == 400
        # Row must be unchanged.
        @test Tables.rowtable(DBInterface.execute(db,
            "SELECT data_dir FROM setups WHERE experiment_id = ?", [exp_id]))[1].data_dir ==
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

@testset "experiment stats carry a collection timeline of samples in measured order" begin
    tmp = mktempdir()
    db = open_prepared_clone(tmp)
    eid = HimalayaUI.init_experiment!(db;
        name = "Timeline", path = tmp,
        data_dir = joinpath(tmp, "data"),
        analysis_dir = joinpath(tmp, "analysis"))
    later = HimalayaUI.create_sample!(db; experiment_id = eid, name = "later")
    first_measured = HimalayaUI.create_sample!(db; experiment_id = eid, name = "first")
    HimalayaUI.create_exposure!(db; experiment_id = eid, sample_id = first_measured,
        filename = "f1", timestamp = "2026-04-25T10:00:15", exposure_time = 15.0)
    HimalayaUI.create_exposure!(db; experiment_id = eid, sample_id = first_measured,
        filename = "f2", timestamp = "2026-04-25T10:02:15", exposure_time = 15.0)
    off = HimalayaUI.create_exposure!(db; experiment_id = eid, sample_id = later,
        filename = "l1", timestamp = "2026-04-25T11:00:00")
    HimalayaUI.create_exposure!(db; experiment_id = eid, sample_id = later,
        filename = "untimed")
    HimalayaUI.create_exposure!(db; experiment_id = eid, sample_id = later,
        kind = "combined", timestamp = "2026-04-25T12:00:00")
    DBInterface.execute(db, "UPDATE exposures SET trace_median = 0.2 WHERE id = ?", [off])
    DBInterface.execute(db, "INSERT INTO sample_calls (sample_id, state) VALUES (?, 'form_factor')", [first_measured])
    DBInterface.execute(db, "INSERT INTO sample_calls (sample_id, state) VALUES (?, 'null')", [later])

    hidden = HimalayaUI.create_setup!(db; experiment_id = eid, folder = "hidden",
        root_dir = tmp, data_dir = tmp, analysis_dir = tmp, included = false)
    h = HimalayaUI.create_sample!(db; experiment_id = eid, name = "h", setup_id = hidden)
    HimalayaUI.create_exposure!(db; experiment_id = eid, sample_id = h, setup_id = hidden,
        filename = "h1", timestamp = "2026-04-25T09:00:00")

    with_inproc_routes(db) do call
        tl = JSON3.read(String(call("GET", "/api/experiments/$eid").body)).stats.timeline
        @test tl.start == "2026-04-25T10:00:00"
        @test [collect(Int, s) for s in tl.samples] == [[0, 135, 1], [3600, 3600, 0]]
        @test collect(Int, tl.beam_off) == [3600]

        listed = first(filter(e -> e.id == eid, JSON3.read(String(call("GET", "/api/experiments").body))))
        @test listed.stats.timeline == tl
    end
end

@testset "an experiment with no timed images has no timeline" begin
    tmp = mktempdir()
    db = open_prepared_clone(tmp)
    eid = HimalayaUI.create_experiment!(db; path = tmp, data_dir = "data", analysis_dir = "analysis")
    with_inproc_routes(db) do call
        @test JSON3.read(String(call("GET", "/api/experiments/$eid").body)).stats.timeline === nothing
    end
end

@testset "an experiment sends its configuration only under setups" begin
    tmp = mktempdir()
    db = open_prepared_clone(tmp)
    exp_id = HimalayaUI.create_experiment!(db; path=tmp, data_dir="data", analysis_dir="analysis")

    with_inproc_routes(db) do call
        body = JSON3.read(String(call("GET", "/api/experiments/$exp_id").body))
        for k in (:data_dir, :analysis_dir, :beam_center_x, :q_units, :config, :integration_pattern)
            @test !haskey(body, k)
        end
        @test body.setups[1].beam_center_x === nothing
        @test call("GET", "/api/experiments").status == 200
    end
end

@testset "typed q_units is honored even with no other geometry (6e)" begin
    tmp = mktempdir()
    db = open_prepared_clone(tmp)
    exp_id = HimalayaUI.create_experiment!(db; path=tmp, data_dir="data",
        analysis_dir="analysis", q_units="nm^-1", q_units_source="user")

    with_inproc_routes(db) do call
        body = JSON3.read(String(call("GET", "/api/experiments/$exp_id").body))
        @test body.setups[1].q_units == "nm^-1"
        @test body.setups[1].q_units_source == "user"
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

@testset "POST /api/experiments persists first-run geometry overrides as source='user'" begin
    tmp = mktempdir()
    data_dir = joinpath(tmp, "data"); mkpath(data_dir)   # empty → scan derives nothing
    db = open_prepared_clone(tmp)

    with_inproc_routes(db) do call
        r = call("POST", "/api/experiments";
            headers = ["Content-Type" => "application/json"],
            body = Vector{UInt8}(JSON3.write(Dict(
                :data_dir => data_dir,
                :name     => "Geo",
                :geometry => Dict(:beam_center_x => 421.3, :beam_center_y => 836.7,
                                  :energy_kev => 11.2)))))
        @test r.status == 202
        eid = Int(JSON3.read(String(r.body)).id)

        row = Tables.rowtable(DBInterface.execute(db,
            "SELECT beam_center_x, beam_center_x_source, beam_center_y, beam_center_y_source,
                    energy_kev, energy_kev_source, flight_path_m, flight_path_m_source
             FROM setups WHERE experiment_id = ?", [eid]))[1]
        # Provided fields are persisted and stamped 'user' (the scan never clobbers them).
        @test row.beam_center_x == 421.3
        @test row.beam_center_x_source == "user"
        @test row.beam_center_y == 836.7
        @test row.beam_center_y_source == "user"
        @test row.energy_kev == 11.2
        @test row.energy_kev_source == "user"
        # An un-provided field stays 'default' (left for the scan to derive).
        @test row.flight_path_m_source == "default"
    end
end

@testset "POST /api/experiments persists committed geometry with its honest source" begin
    tmp = mktempdir()
    data_dir = joinpath(tmp, "data"); mkpath(data_dir)
    db = open_prepared_clone(tmp)

    with_inproc_routes(db) do call
        # The funnel commits the WHOLE preview geometry: each value carries its
        # real source ('setup'/'prp'), and only an edited field is 'user'.
        r = call("POST", "/api/experiments";
            headers = ["Content-Type" => "application/json"],
            body = Vector{UInt8}(JSON3.write(Dict(
                :data_dir => data_dir, :name => "Geo",
                :geometry => Dict(
                    :beam_center_x => 421.4, :beam_center_x_source => "setup",
                    :flight_path_m => 1.8095, :flight_path_m_source => "setup",
                    :pixel_size_um => 172.0, :pixel_size_um_source => "prp",
                    :energy_kev    => 11.2,  :energy_kev_source    => "user")))))
        @test r.status == 202
        eid = Int(JSON3.read(String(r.body)).id)
        row = Tables.rowtable(DBInterface.execute(db,
            "SELECT beam_center_x_source, flight_path_m_source, pixel_size_um_source,
                    energy_kev_source, flight_path_m FROM setups WHERE experiment_id = ?", [eid]))[1]
        @test row.beam_center_x_source == "setup"   # honest provenance, not 'user'
        @test row.flight_path_m_source == "setup"
        @test row.flight_path_m == 1.8095           # the setup value, not a PRP re-derive
        @test row.pixel_size_um_source == "prp"
        @test row.energy_kev_source == "user"       # the one edited field
    end
end

@testset "_fill_unset_geometry! fills only 'default' fields, never an established one" begin
    tmp = mktempdir()
    db = open_prepared_clone(tmp)
    eid = HimalayaUI.create_experiment!(db; path = tmp, data_dir = joinpath(tmp, "d"),
        analysis_dir = joinpath(tmp, "a"),
        flight_path_m = 1.8095, flight_path_m_source = "setup",      # established
        beam_center_x = 421.4, beam_center_x_source = "setup")        # established
    # A later scan re-derives DIFFERENT values (e.g. from PRP only, no setup file).
    geo = (energy_kev = 9.0, energy_kev_source = "prp",
           flight_path_m = 1.7, flight_path_m_source = "prp",          # would clobber
           beam_center_x = missing, beam_center_x_source = "prp",
           beam_center_y = missing, beam_center_y_source = "prp",
           pixel_size_um = 172.0, pixel_size_um_source = "prp")
    HimalayaUI._fill_unset_geometry!(db, HimalayaUI.default_setup_id(db, eid), geo)
    row = Tables.rowtable(DBInterface.execute(db,
        "SELECT flight_path_m, flight_path_m_source, energy_kev, energy_kev_source,
                beam_center_x, pixel_size_um FROM setups WHERE experiment_id = ?", [eid]))[1]
    @test row.flight_path_m == 1.8095                # established 'setup' value PRESERVED
    @test row.flight_path_m_source == "setup"
    @test row.energy_kev == 9.0                      # was 'default' → filled
    @test row.energy_kev_source == "prp"
    @test row.beam_center_x == 421.4                 # established, untouched
    @test row.pixel_size_um == 172.0                 # was 'default' → filled
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
            "SELECT COUNT(*) AS c FROM setups WHERE data_dir = ?", [data_dir]))[1].c
        @test n == 1
    end
end
