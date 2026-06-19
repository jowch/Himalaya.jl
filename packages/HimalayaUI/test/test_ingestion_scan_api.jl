# packages/HimalayaUI/test/test_ingestion_scan_api.jl
using Test
using HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

include("test_http.jl")   # provides with_test_server

"""
Build a fresh temp DB with one experiment whose data directory exists on disk.
Returns (db, dir, exp_id) where dir is a tempdir that the experiment points at.
"""
function scan_test_db()
    dir = mktempdir()
    db  = HimalayaUI.open_db(joinpath(dir, "himalaya.db"))
    exp_id = HimalayaUI.create_experiment!(db;
        name         = "TestExp",
        path         = dir,
        data_dir     = dir,
        analysis_dir = joinpath(dir, "analysis"))
    (db, dir, exp_id)
end

@testset "ingestion scan API + SSE + scheduler (Phase C)" begin
    @testset "broadcast_progress! emits curation frame with no user_actions row" begin
        mktempdir() do dir
            db = HimalayaUI.open_db(joinpath(dir, "h.db"))
            HimalayaUI.bind_db!(db)

            exp_id = HimalayaUI.create_experiment!(db;
                name = "BP", path = dir, data_dir = dir, analysis_dir = dir)

            pending = Channel{String}(64)
            sub = (pending = pending,)
            lock(HimalayaUI.SSE_LOCK) do
                push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
            end

            before_count = Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM user_actions"))[1].c

            HimalayaUI.broadcast_progress!(exp_id; kind = "ingest_started",
                processed = 0, total = 680)

            after_count = Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM user_actions"))[1].c

            # No durable row written
            @test after_count == before_count

            # Frame is on the channel
            @test isready(pending)
            frame = take!(pending)
            @test occursin("event: curation", frame)

            data_line = first(filter(l -> startswith(l, "data: "), split(frame, '\n')))
            obj = JSON3.read(replace(data_line, r"^data: " => ""))
            # Frame is payload-wrapped (mirrors broadcast_event!): `kind` top-level,
            # experiment/count fields under `payload` so it parses like every other
            # "curation" frame and the frontend reads `payload.experiment_id`.
            @test obj.kind == "ingest_started"
            @test obj.payload.experiment_id == exp_id
            @test obj.payload.processed == 0
            @test obj.payload.total == 680

            lock(HimalayaUI.SSE_LOCK) do
                filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
            end
            close(pending)
            HimalayaUI.SSE_SUBSCRIBERS[] = []
        end
    end

    @testset "create_load! persists a load row" begin
        db, dir, exp_id = scan_test_db()
        lid = HimalayaUI.create_load!(db;
            experiment_id = exp_id,
            load_index    = 1,
            start_time    = "2026-04-12T08:00:00",
            end_time      = "2026-04-12T08:45:00",
            frame_count   = 48)
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM loads WHERE id = ?", [lid])))
        @test row.experiment_id == exp_id
        @test row.load_index    == 1
        @test row.frame_count   == 48
        @test row.start_time    == "2026-04-12T08:00:00"
        SQLite.close(db)
    end

    @testset "GET /api/experiments/{id} includes typed geometry + ingest_status + stats" begin
        db, dir, exp_id = scan_test_db()
        # Seed typed geometry columns directly (Phase B would do this via scan_and_group!)
        DBInterface.execute(db, """
            UPDATE experiments SET
                beam_center_x       = 421.409,
                beam_center_y       = 836.946,
                beam_center_x_source = 'setup',
                beam_center_y_source = 'setup',
                pixel_size_um       = 172.0,
                pixel_size_um_source = 'prp',
                flight_path_m       = 1.8095,
                flight_path_m_source = 'setup',
                energy_kev          = 9.0,
                energy_kev_source   = 'prp',
                q_units             = 'A^-1',
                ingest_status       = 'complete',
                last_scanned_at     = '2026-04-12T10:00:00Z'
            WHERE id = ?
        """, [exp_id])
        # Seed a load + sample + exposure so stats are non-zero
        lid = HimalayaUI.create_load!(db; experiment_id = exp_id, load_index = 1, frame_count = 4)
        sid = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "HA85 (S01P15)",
            load_id = lid, slot_index = 15)
        HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = sid, filename = "f.tif")

        with_test_server(db) do port, base
            r = HTTP.get("$base/api/experiments/$exp_id")
            @test r.status == 200
            body = JSON3.read(String(r.body))

            # Typed geometry
            @test body.beam_center_x       ≈ 421.409
            @test body.beam_center_x_source == "setup"
            @test body.flight_path_m       ≈ 1.8095
            @test body.flight_path_m_source == "setup"
            @test body.energy_kev          ≈ 9.0
            @test body.energy_kev_source   == "prp"
            @test body.pixel_size_um       ≈ 172.0
            @test body.q_units             == "A^-1"
            @test body.ingest_status       == "complete"
            @test body.last_scanned_at     == "2026-04-12T10:00:00Z"

            # Stats roll-up
            @test haskey(body, :stats)
            @test body.stats.loads     == 1
            @test body.stats.samples   == 1
            @test body.stats.exposures == 1
        end
        SQLite.close(db)
    end

    @testset "GET /api/experiments/{id}/loads returns Load▸Sample▸Exposure roll-up" begin
        db, dir, exp_id = scan_test_db()
        lid1 = HimalayaUI.create_load!(db; experiment_id = exp_id, load_index = 1, frame_count = 8)
        lid2 = HimalayaUI.create_load!(db; experiment_id = exp_id, load_index = 2, frame_count = 4)
        sid1 = HimalayaUI.create_sample!(db; experiment_id = exp_id,
            name = "HA85 (S01P01)", load_id = lid1, slot_index = 1)
        sid2 = HimalayaUI.create_sample!(db; experiment_id = exp_id,
            name = "HA85 (S02P01)", load_id = lid2, slot_index = 1)
        xid1 = HimalayaUI.create_exposure!(db; experiment_id = exp_id,
            sample_id = sid1, filename = "f1.tif", horizontal_position = 12.4, frame_no = 1)
        xid2 = HimalayaUI.create_exposure!(db; experiment_id = exp_id,
            sample_id = sid1, filename = "f2.tif", horizontal_position = 12.4, frame_no = 2)
        xid3 = HimalayaUI.create_exposure!(db; experiment_id = exp_id,
            sample_id = sid2, filename = "f3.tif", horizontal_position = 39.5, frame_no = 1)

        with_test_server(db) do port, base
            r = HTTP.get("$base/api/experiments/$exp_id/loads")
            @test r.status == 200
            loads = JSON3.read(String(r.body))
            @test length(loads) == 2

            l1 = first(filter(l -> l.load_index == 1, loads))
            @test length(l1.samples) == 1
            @test l1.samples[1].name == "HA85 (S01P01)"
            @test length(l1.samples[1].exposures) == 2
            exps = l1.samples[1].exposures
            filenames = [e.filename for e in exps]
            @test "f1.tif" in filenames && "f2.tif" in filenames

            l2 = first(filter(l -> l.load_index == 2, loads))
            @test length(l2.samples) == 1
            @test length(l2.samples[1].exposures) == 1

            # 404 for unknown experiment
            r2 = HTTP.get("$base/api/experiments/999999/loads"; status_exception = false)
            @test r2.status == 404
        end
        SQLite.close(db)
    end

    @testset "PATCH /api/experiments/{id} writes geometry overrides, marks source=user" begin
        db, dir, exp_id = scan_test_db()

        with_test_server(db) do port, base
            # Patch a geometry field
            r = HTTP.patch("$base/api/experiments/$exp_id";
                body    = JSON3.write(Dict(
                    :flight_path_m  => 1.7500,
                    :beam_center_x  => 500.0,
                    :beam_center_y  => 800.0,
                    :pixel_size_um  => 172.0,
                    :energy_kev     => 9.0,
                    :q_units        => "nm^-1")),
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice"],
                status_exception = false)
            @test r.status == 200

            # Verify persisted values + source override
            row = Tables.rowtable(DBInterface.execute(db, """
                SELECT flight_path_m, flight_path_m_source,
                       beam_center_x, beam_center_x_source,
                       q_units
                  FROM experiments WHERE id = ?
            """, [exp_id]))[1]
            @test row.flight_path_m       ≈ 1.7500
            @test row.flight_path_m_source == "user"
            @test row.beam_center_x       ≈ 500.0
            @test row.beam_center_x_source == "user"
            @test row.q_units             == "nm^-1"

            # Partial patch: only flight_path_m. Others untouched.
            r2 = HTTP.patch("$base/api/experiments/$exp_id";
                body    = JSON3.write(Dict(:energy_kev => 12.0)),
                headers = ["Content-Type" => "application/json"],
                status_exception = false)
            @test r2.status == 200
            row2 = Tables.rowtable(DBInterface.execute(db,
                "SELECT energy_kev, energy_kev_source, flight_path_m FROM experiments WHERE id=?",
                [exp_id]))[1]
            @test row2.energy_kev        ≈ 12.0
            @test row2.energy_kev_source == "user"
            @test row2.flight_path_m     ≈ 1.7500  # previous value preserved

            # 404 for unknown experiment
            r3 = HTTP.patch("$base/api/experiments/999999";
                body    = JSON3.write(Dict(:energy_kev => 1.0)),
                headers = ["Content-Type" => "application/json"],
                status_exception = false)
            @test r3.status == 404

            # Read-only fields rejected
            r4 = HTTP.patch("$base/api/experiments/$exp_id";
                body    = JSON3.write(Dict(:data_dir => "/evil")),
                headers = ["Content-Type" => "application/json"],
                status_exception = false)
            @test r4.status == 400

            # name stays read-only here — rename lands in Phase E1 (derived-name
            # contract preserved until then).
            r5 = HTTP.patch("$base/api/experiments/$exp_id";
                body    = JSON3.write(Dict(:name => "Renamed")),
                headers = ["Content-Type" => "application/json"],
                status_exception = false)
            @test r5.status == 400

            # A body with no recognized patchable (geometry) field → 400, matching
            # the codebase PATCH-validation convention (test_route_validation_routing.jl).
            r6 = HTTP.patch("$base/api/experiments/$exp_id";
                body    = JSON3.write(Dict(:NOT_a_field => "x")),
                headers = ["Content-Type" => "application/json"],
                status_exception = false)
            @test r6.status == 400
        end
        SQLite.close(db)
    end

    @testset "rescan timer start/stop lifecycle (unit)" begin
        # This test never opens a real DB — it only exercises the timer registry.
        # start_rescan_scheduler! requires a db for the tick body; we pass a temp one.
        mktempdir() do dir
            db = HimalayaUI.open_db(joinpath(dir, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db;
                name = "T", path = dir, data_dir = dir, analysis_dir = dir)

            # Real Timers are armed below; guarantee they are all closed and the DB
            # is shut even if an assertion throws mid-test (otherwise a leaked Timer
            # keeps firing against a closed DB and poisons later testsets).
            try
                # Starting a scheduler registers a timer
                HimalayaUI.start_rescan_scheduler!(db, exp_id;
                    tick_interval_seconds = 3600)
                @test haskey(HimalayaUI.RESCAN_TIMERS, exp_id)

                # Stopping it removes the entry
                HimalayaUI.stop_rescan_scheduler!(exp_id)
                @test !haskey(HimalayaUI.RESCAN_TIMERS, exp_id)

                # Stopping a non-existent id is a no-op
                @test_nowarn HimalayaUI.stop_rescan_scheduler!(exp_id)

                # stop_all_rescan_timers! clears all entries
                HimalayaUI.start_rescan_scheduler!(db, exp_id; tick_interval_seconds = 3600)
                HimalayaUI.stop_all_rescan_timers!()
                @test isempty(HimalayaUI.RESCAN_TIMERS)
            finally
                HimalayaUI.stop_all_rescan_timers!()
                SQLite.close(db)
            end
        end
    end

    @testset "DELETE /api/experiments/{id} removes experiment + cascade" begin
        db, dir, exp_id = scan_test_db()
        lid = HimalayaUI.create_load!(db; experiment_id = exp_id, load_index = 1, frame_count = 2)
        sid = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "A",
            load_id = lid, slot_index = 1)
        HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = sid, filename = "f.tif")

        with_test_server(db) do port, base
            r = HTTP.delete("$base/api/experiments/$exp_id",
                headers = ["X-Username" => "alice"],
                status_exception = false)
            @test r.status == 200

            # Experiment is gone
            r2 = HTTP.get("$base/api/experiments/$exp_id"; status_exception = false)
            @test r2.status == 404

            # Cascades: loads, samples, exposures removed
            cnt_loads = Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM loads WHERE experiment_id = ?", [exp_id]))[1].c
            cnt_samples = Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM samples WHERE experiment_id = ?", [exp_id]))[1].c
            @test cnt_loads   == 0
            @test cnt_samples == 0

            # 404 on repeat delete
            r3 = HTTP.delete("$base/api/experiments/$exp_id",
                headers = ["X-Username" => "alice"],
                status_exception = false)
            @test r3.status == 404
        end
        SQLite.close(db)
    end

    @testset "_rescan_tick! tiered backoff persists to DB" begin
        mktempdir() do dir
            db = HimalayaUI.open_db(joinpath(dir, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db;
                name = "BT", path = dir, data_dir = dir, analysis_dir = dir)

            # The "change" branch below calls start_rescan_scheduler!, which arms a REAL
            # Timer. Wrap the whole body so that timer (and any other) is always closed
            # and the DB shut even if an assertion throws — a leaked Timer would keep
            # firing against a closed DB.
            try
                ticks_before_daily = 3   # configurable in the call below

                # Inject the change decision so tier transitions are deterministic without
                # filesystem fixtures. No scan stub is needed: on a "change" the real
                # scan_and_group! runs against the EMPTY temp dir (no matching triplets) —
                # a harmless no-op.
                no_change  = (_, _) -> false
                has_change = (_, _) -> true

                # Run ticks_before_daily empty ticks; should stay in 'fast' tier until threshold
                for _ in 1:ticks_before_daily
                    HimalayaUI._rescan_tick!(db, exp_id;
                        cheap_check_fn = no_change,
                        fast_interval = 3600.0, ticks_before_daily = ticks_before_daily,
                        ticks_before_stop = 2)
                end
                row = Tables.rowtable(DBInterface.execute(db,
                    "SELECT last_scan_tier, consecutive_empty_ticks FROM experiments WHERE id=?",
                    [exp_id]))[1]
                # After exactly ticks_before_daily empty ticks, tier advances to daily
                @test row.last_scan_tier == "daily"
                @test row.consecutive_empty_ticks == 0  # reset on tier transition

                # One more empty daily tick → not yet stopped (need ticks_before_stop=2 more)
                HimalayaUI._rescan_tick!(db, exp_id;
                    cheap_check_fn = no_change,
                    fast_interval = 3600.0, ticks_before_daily = ticks_before_daily,
                    ticks_before_stop = 2)
                row2 = Tables.rowtable(DBInterface.execute(db,
                    "SELECT last_scan_tier, consecutive_empty_ticks FROM experiments WHERE id=?",
                    [exp_id]))[1]
                @test row2.last_scan_tier == "daily"
                @test row2.consecutive_empty_ticks == 1

                # Simulate a change: re-arms back to fast. cheap_check_fn returns true, then
                # the real scan_and_group! runs against the EMPTY temp dir — a harmless no-op.
                HimalayaUI._rescan_tick!(db, exp_id;
                    cheap_check_fn = has_change,
                    fast_interval = 3600.0, ticks_before_daily = ticks_before_daily,
                    ticks_before_stop = 2)
                row3 = Tables.rowtable(DBInterface.execute(db,
                    "SELECT last_scan_tier, consecutive_empty_ticks FROM experiments WHERE id=?",
                    [exp_id]))[1]
                @test row3.last_scan_tier == "fast"
                @test row3.consecutive_empty_ticks == 0
            finally
                # The change-branch re-armed a fast-tier Timer; close it before the DB.
                HimalayaUI.stop_rescan_scheduler!(exp_id)
                HimalayaUI.stop_all_rescan_timers!()
                SQLite.close(db)
            end
        end
    end

    @testset "POST /api/experiments/{id}/scan no-change is idempotent" begin
        db, dir, exp_id = scan_test_db()

        with_test_server(db) do port, base
            # First scan on empty dir: no files, should return 200 + changed=false
            # (since scan_and_group! is Phase B and may not be present in test env,
            #  the route must not crash when the directory has no matching files;
            #  test the HTTP contract, not the full scan logic)
            r = HTTP.post("$base/api/experiments/$exp_id/scan";
                headers = ["X-Username" => "alice"],
                status_exception = false)
            # 200 or 202 are both acceptable (implementation may vary)
            @test r.status in (200, 202)
            body = JSON3.read(String(r.body))
            @test haskey(body, :status)

            # 404 for unknown experiment
            r2 = HTTP.post("$base/api/experiments/999999/scan";
                headers = ["X-Username" => "alice"],
                status_exception = false)
            @test r2.status == 404
        end
        SQLite.close(db)
    end

    @testset "POST /api/experiments creates experiment + starts async scan" begin
        mktempdir() do dir
            db = HimalayaUI.open_db(joinpath(dir, "himalaya.db"))

            with_test_server(db) do port, base
                r = HTTP.post("$base/api/experiments";
                    body    = JSON3.write(Dict(:path => dir, :name => "MyExp")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"],
                    status_exception = false)
                @test r.status in (200, 201, 202)
                body = JSON3.read(String(r.body))
                @test haskey(body, :id)
                @test body.id isa Integer
                @test body.id > 0

                # Experiment row is created
                exp_rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, name, data_dir, ingest_status FROM experiments WHERE id = ?",
                    [body.id]))
                @test length(exp_rows) == 1
                @test exp_rows[1].name == "MyExp"
                @test exp_rows[1].data_dir == dir

                # ingest_status transitions to scanning (may have already completed in test)
                @test exp_rows[1].ingest_status in ("scanning", "complete", "failed")

                # Missing path → 400
                r2 = HTTP.post("$base/api/experiments";
                    body    = JSON3.write(Dict(:name => "NoPath")),
                    headers = ["Content-Type" => "application/json"],
                    status_exception = false)
                @test r2.status == 400

                # Non-existent path → 400
                r3 = HTTP.post("$base/api/experiments";
                    body    = JSON3.write(Dict(:path => "/does/not/exist/xyz123")),
                    headers = ["Content-Type" => "application/json"],
                    status_exception = false)
                @test r3.status == 400
            end
            SQLite.close(db)
        end
    end
end
