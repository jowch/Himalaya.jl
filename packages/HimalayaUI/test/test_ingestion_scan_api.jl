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

@testset "ingestion scan API + SSE" begin
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

    @testset "broadcast_progress! carries a phase discriminator in the payload" begin
        # The rescan route tags its frames phase="rescan" so the frontend can map
        # them to the "analyzing" (re-analysis) surface vs "scanning" (initial).
        # broadcast_progress! merges kwargs into the payload, so phase rides through.
        mktempdir() do dir
            db = HimalayaUI.open_db(joinpath(dir, "h.db"))
            HimalayaUI.bind_db!(db)
            exp_id = HimalayaUI.create_experiment!(db;
                name = "RP", path = dir, data_dir = dir, analysis_dir = dir)

            pending = Channel{String}(64)
            sub = (pending = pending,)
            lock(HimalayaUI.SSE_LOCK) do
                push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
            end

            HimalayaUI.broadcast_progress!(exp_id; kind = "ingest_progress",
                processed = 12, total = 40, phase = "rescan")

            @test isready(pending)
            frame = take!(pending)
            data_line = first(filter(l -> startswith(l, "data: "), split(frame, '\n')))
            obj = JSON3.read(replace(data_line, r"^data: " => ""))
            @test obj.kind == "ingest_progress"
            @test obj.payload.phase == "rescan"
            @test obj.payload.processed == 12

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

    @testset "GET /api/experiments/{id} includes setup geometry + ingest_status + stats" begin
        db, dir, exp_id = scan_test_db()
        # Seed typed geometry columns directly (Phase B would do this via scan_and_group!)
        DBInterface.execute(db, """
            UPDATE setups SET
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
                q_units             = 'A^-1'
            WHERE experiment_id = ?
        """, [exp_id])
        DBInterface.execute(db,
            "UPDATE experiments SET ingest_status = 'complete', last_scanned_at = '2026-04-12T10:00:00Z' WHERE id = ?",
            [exp_id])
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
            st = only(body.setups)
            @test st.beam_center_x       ≈ 421.409
            @test st.beam_center_x_source == "setup"
            @test st.flight_path_m       ≈ 1.8095
            @test st.flight_path_m_source == "setup"
            @test st.energy_kev          ≈ 9.0
            @test st.energy_kev_source   == "prp"
            @test st.pixel_size_um       ≈ 172.0
            @test st.q_units             == "A^-1"
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

    @testset "DELETE /api/experiments/{id} refuses while an ingest runs" begin
        db, dir, exp_id = scan_test_db()
        DBInterface.execute(db, "UPDATE experiments SET ingest_status = 'scanning' WHERE id = ?", [exp_id])
        with_inproc_routes(db) do call
            @test call("DELETE", "/api/experiments/$exp_id").status == 409
            DBInterface.execute(db, "UPDATE experiments SET ingest_status = 'failed' WHERE id = ?", [exp_id])
            @test call("DELETE", "/api/experiments/$exp_id").status == 200
        end
        SQLite.close(db)
    end

    @testset "DELETE /api/experiments/{id} tells other clients the experiment is gone" begin
        db, dir, exp_id = scan_test_db()
        pending = Channel{String}(64)
        sub = (pending = pending,)
        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end
        try
            with_inproc_routes(db) do call
                @test call("DELETE", "/api/experiments/$exp_id").status == 200
            end
            frames = String[]
            while isready(pending)
                push!(frames, take!(pending))
            end
            objs = [JSON3.read(replace(first(filter(l -> startswith(l, "data: "), split(f, '\n'))), r"^data: " => "")) for f in frames]
            deleted = filter(o -> o.kind == "experiment_deleted", objs)
            @test length(deleted) == 1
            @test only(deleted).entity_type == "experiment"
            @test only(deleted).entity_id == exp_id
        finally
            lock(HimalayaUI.SSE_LOCK) do
                filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
            end
            close(pending)
        end
        SQLite.close(db)
    end

    @testset "a rescan pass emits rescan-phase ingest frames and a terminal frame" begin
        mktempdir() do dir
            db = HimalayaUI.open_db(joinpath(dir, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db;
                name = "RT", path = dir, data_dir = dir, analysis_dir = dir)
            pending = Channel{String}(64)
            sub = (pending = pending,)
            lock(HimalayaUI.SSE_LOCK) do
                push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
            end
            try
                @test HimalayaUI.run_ingest!(db, exp_id; phase = "rescan") === :ok
                HimalayaUI.wait_prewarm(exp_id)

                kinds  = String[]
                phases = String[]
                while isready(pending)
                    frame = take!(pending)
                    data_line = first(filter(l -> startswith(l, "data: "), split(frame, '\n')))
                    obj = JSON3.read(replace(data_line, r"^data: " => ""))
                    push!(kinds, String(obj.kind))
                    haskey(obj.payload, :phase) && push!(phases, String(obj.payload.phase))
                end
                @test "ingest_started" in kinds
                @test "ingest_complete" in kinds
                @test "rescan" in phases
            finally
                lock(HimalayaUI.SSE_LOCK) do
                    filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
                end
                close(pending)
                HimalayaUI.SSE_SUBSCRIBERS[] = []
                SQLite.close(db)
            end
        end
    end

    @testset "POST /api/experiments/{id}/scan answers each case of its contract" begin
        db, dir, exp_id = scan_test_db()
        pending = Channel{String}(256)
        sub = (pending = pending,)
        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end
        kinds() = begin
            out = String[]
            while isready(pending)
                frame = take!(pending)
                lines = filter(l -> startswith(l, "data: "), split(frame, '\n'))
                isempty(lines) && continue
                obj = JSON3.read(replace(first(lines), r"^data: " => ""))
                haskey(obj, :kind) && push!(out, String(obj.kind))
            end
            out
        end
        wait_until(f; timeout = 30) = begin
            t0 = time()
            while !f() && time() - t0 < timeout
                sleep(0.05)
            end
            f()
        end
        post(body = nothing) = HTTP.post(base_url[] * "/api/experiments/$exp_id/scan";
            body = body === nothing ? UInt8[] : JSON3.write(body),
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
            status_exception = false)
        status_of() = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT ingest_status FROM experiments WHERE id = ?", [exp_id]))).ingest_status
        base_url = Ref("")
        try
            with_test_server(db) do port, base
                base_url[] = base
                DBInterface.execute(db, "UPDATE experiments SET ingest_status = 'complete' WHERE id = ?", [exp_id])
                kinds()

                r = post()
                @test r.status == 200
                body = JSON3.read(String(r.body))
                @test body.status == "complete"
                @test body.changed == false
                sleep(0.5)
                @test !("ingest_started" in kinds())

                DBInterface.execute(db, "UPDATE experiments SET ingest_status = 'failed' WHERE id = ?", [exp_id])
                r = post()
                @test r.status == 200
                @test JSON3.read(String(r.body)).status == "failed"
                sleep(0.5)
                @test !("ingest_started" in kinds())

                write(joinpath(dir, "NEW_S0001_0_001.tif"), "fake tif")
                r = post()
                @test r.status == 202
                @test JSON3.read(String(r.body)).status == "scanning"
                @test wait_until(() -> status_of() == "complete" && !HimalayaUI.ingest_running(exp_id))
                HimalayaUI.wait_prewarm(exp_id)
                seen = kinds()
                @test "ingest_started" in seen
                @test "ingest_complete" in seen

                r = post()
                @test r.status == 200
                sleep(0.5)
                @test !("ingest_started" in kinds())

                r = post(Dict(:force => true))
                @test r.status == 202
                @test JSON3.read(String(r.body)).status == "scanning"
                @test wait_until(() -> status_of() == "complete" && !HimalayaUI.ingest_running(exp_id))
                HimalayaUI.wait_prewarm(exp_id)
                @test "ingest_started" in kinds()

                r = post(Dict(:force => false))
                @test r.status == 200

                r404 = HTTP.post("$base/api/experiments/999999/scan";
                    headers = ["X-Username" => "alice"], status_exception = false)
                @test r404.status == 404
            end
        finally
            lock(HimalayaUI.SSE_LOCK) do
                filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
            end
            close(pending)
            SQLite.close(db)
        end
    end

    @testset "a plain rescan with no new files still backfills missing beam monitors" begin
        db, dir, exp_id = scan_test_db()
        stems = ["MON_S$(lpad(i, 4, '0'))_0_001" for i in 1:2]
        for (i, stem) in enumerate(stems)
            write(joinpath(dir, "$stem.tif"), "fake tif")
            write(joinpath(dir, "$stem.prp"),
                "I_signal=$(7000.0 + i) (corrected for dark counts)\nI_0(daqIon3.0)=1.5\n" *
                "Exposure time=15\nCounting time=15.8\n")
            HimalayaUI.create_exposure!(db; experiment_id = exp_id, filename = stem,
                prp_path = joinpath(dir, "$stem.prp"), monitors_read = false)
        end
        monitors() = Tables.rowtable(DBInterface.execute(db,
            "SELECT transmitted_counts, counting_time FROM exposures WHERE experiment_id = ? ORDER BY filename",
            [exp_id]))
        status() = Tables.rowtable(DBInterface.execute(db,
            "SELECT ingest_status FROM experiments WHERE id = ?", [exp_id]))[1].ingest_status
        try
            with_test_server(db) do port, base
                r = HTTP.post("$base/api/experiments/$exp_id/scan";
                    body = JSON3.write(Dict(:force => false)),
                    headers = ["Content-Type" => "application/json"],
                    status_exception = false)
                @test r.status in (200, 202)
                @test timedwait(() -> all(!ismissing(m.counting_time) for m in monitors()) &&
                                      status() != "scanning", 30.0) === :ok
            end
            @test [m.transmitted_counts for m in monitors()] ≈ [7001.0, 7002.0]
            @test all(m.counting_time ≈ 15.8 for m in monitors())
        finally
            SQLite.close(db)
        end
    end

    @testset "an ingest pass backfills missing beam monitors" begin
        db, dir, exp_id = scan_test_db()
        stem = "TICK_S0001_0_001"
        write(joinpath(dir, "$stem.prp"), "I_signal=6500.0\nExposure time=15\nCounting time=15.8\n")
        HimalayaUI.create_exposure!(db; experiment_id = exp_id, filename = stem,
            prp_path = joinpath(dir, "$stem.prp"), monitors_read = false)
        try
            HimalayaUI.run_ingest!(db, exp_id)
            HimalayaUI.wait_prewarm(exp_id)
            r = Tables.rowtable(DBInterface.execute(db,
                "SELECT transmitted_counts FROM exposures WHERE experiment_id = ?", [exp_id]))[1]
            @test r.transmitted_counts ≈ 6500.0
        finally
            SQLite.close(db)
        end
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
                    "SELECT x.id, x.name, st.data_dir, x.ingest_status FROM experiments x JOIN setups st ON st.experiment_id = x.id WHERE x.id = ?",
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

    @testset "POST /api/experiments persists patterns from the nested body" begin
        db = HimalayaUI.open_db(joinpath(mktempdir(), "h.db"))
        with_inproc_routes(db) do call
            dir = mktempdir(); touch(joinpath(dir, "s.tif"))
            body = Vector{UInt8}(JSON3.write(Dict(
                "path"     => dir,
                "patterns" => Dict("image" => "{name}.tiff"))))
            resp = call("POST", "/api/experiments";
                headers = ["Content-Type" => "application/json"], body = body)
            @test resp.status == 202
            id  = JSON3.read(resp.body).id
            row = only(Tables.rowtable(DBInterface.execute(db,
                "SELECT image_pattern FROM setups WHERE experiment_id = ?", [id])))
            @test row.image_pattern == "{name}.tiff"
            # Poll for the async scan to finish before with_inproc_routes closes the DB.
            deadline = time() + 5.0
            while time() < deadline
                row2 = only(Tables.rowtable(DBInterface.execute(db,
                    "SELECT ingest_status FROM experiments WHERE id = ?", [id])))
                row2.ingest_status in ("complete", "failed") && break
                sleep(0.05)
            end
            final_row = only(Tables.rowtable(DBInterface.execute(db,
                "SELECT ingest_status FROM experiments WHERE id = ?", [id])))
            @test final_row.ingest_status == "complete"
        end
        SQLite.close(db)
    end
end

@testset "ingest_progress carries `stage` on the wire, distinct from `phase`" begin
    # broadcast_progress! merges kwargs... into the payload UNTYPED (events.jl), so a
    # renamed or typo'd kwarg would ship silently — nothing else asserts that `stage`
    # actually reaches the frame. docs/contract-testing.md wants a row at every
    # layer; this is the wire one.
    pending = Channel{String}(64)
    sub = (pending = pending,)
    lock(HimalayaUI.SSE_LOCK) do
        empty!(HimalayaUI.SSE_SUBSCRIBERS[])
        push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
    end
    try
        HimalayaUI.broadcast_progress!(9; kind = "ingest_progress",
            processed = 780, total = 1100, phase = "rescan", stage = "discovery")
        @test isready(pending)
        frame = take!(pending)
        data_line = first(filter(l -> startswith(l, "data: "), split(frame, '\n')))
        obj = JSON3.read(replace(data_line, r"^data: " => ""))
        # BOTH discriminators ride the payload independently: `phase` selects the
        # frontend surface, `stage` selects the progress-bar segment.
        @test obj.payload.stage == "discovery"
        @test obj.payload.phase == "rescan"
        @test obj.payload.processed == 780
        @test obj.payload.total == 1100

        # A frame with no stage must OMIT the key (the frontend falls back to the
        # plain single-track bar), not send null or an empty string.
        HimalayaUI.broadcast_progress!(9; kind = "ingest_started", processed = 0, total = 0)
        frame2 = take!(pending)
        line2 = first(filter(l -> startswith(l, "data: "), split(frame2, '\n')))
        obj2 = JSON3.read(replace(line2, r"^data: " => ""))
        @test !haskey(obj2.payload, :stage)
    finally
        lock(HimalayaUI.SSE_LOCK) do
            empty!(HimalayaUI.SSE_SUBSCRIBERS[])
        end
    end
end

@testset "an ingest pass ingests a NON-empty directory (3-arg on_progress)" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        data_dir = mktempdir()
        touch(joinpath(data_dir, "e1.tif"))
        write(joinpath(data_dir, "e1.prp"), "Detector: Pilatus 1M\n")

        id = HimalayaUI.create_experiment!(db; name = "rescan-nonempty", path = data_dir,
                                          data_dir = data_dir, analysis_dir = data_dir,
                                          ingest_status = "complete")
        HimalayaUI.run_ingest!(db, id; phase = "rescan")
        HimalayaUI.wait_prewarm(id)

        n = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) c FROM exposures WHERE experiment_id = ?", [id]))).c
        @test n == 1
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT ingest_status FROM experiments WHERE id = ?", [id])))
        @test row.ingest_status != "failed"
        SQLite.close(db)
    end
end
