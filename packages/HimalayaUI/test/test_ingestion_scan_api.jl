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
end
