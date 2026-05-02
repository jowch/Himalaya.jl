using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

@testset "apply_event! writes to user_actions and returns event_id" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x",
            ["X-Username" => "alice"], UInt8[])

        result = HimalayaUI.apply_event!(db, req;
            kind = "test_kind", entity_type = "exposure", entity_id = 42,
            payload = Dict(:foo => "bar"))
        eid = result.event_id
        @test eid > 0

        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM user_actions WHERE id = ?", [eid])))
        @test String(row.action) == "test_kind"
        @test String(row.entity_type) == "exposure"
        @test Int(row.entity_id) == 42
        @test JSON3.read(row.payload).foo == "bar"
        # user_id should resolve to alice's id.
        @test !ismissing(row.user_id)
    end
end

@testset "apply_event! with no payload writes NULL payload and skips dispatcher" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "bob"], UInt8[])
        result = HimalayaUI.apply_event!(db, req;
            kind = "no_payload_event", entity_type = "exposure", entity_id = 1)
        eid = result.event_id
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT payload FROM user_actions WHERE id = ?", [eid])))
        @test ismissing(row.payload)
        @test result.view_row_id === nothing
    end
end

@testset "apply_event! with no X-Username writes user_id = NULL" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", Pair{String,String}[], UInt8[])
        result = HimalayaUI.apply_event!(db, req;
            kind = "anon_event", entity_type = "exposure", entity_id = 1,
            payload = Dict(:k => "v"))
        eid = result.event_id
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT user_id FROM user_actions WHERE id = ?", [eid])))
        @test ismissing(row.user_id)
    end
end

@testset "log_action! still works (legacy wrapper)" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        # Existing call sites pass these args; should write a row with the
        # legacy column shape and a payload of {:note => "..."} when given.
        HimalayaUI.log_action!(db, req;
            action = "set_status", entity_type = "exposure", entity_id = 7,
            note = "rejected")
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT action, entity_type, entity_id, payload FROM user_actions
             ORDER BY id DESC LIMIT 1")))
        @test String(row.action) == "set_status"
        @test JSON3.read(row.payload).note == "rejected"
    end
end

@testset "idx_events_by_exposure index exists" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='index' AND name='idx_events_by_exposure'"))
        @test !isempty(rows)
    end
end

@testset "user_actions has payload and undoes_event_id columns after migration" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        cols = [String(c.name) for c in Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(user_actions)"))]
        @test "payload" in cols
        @test "undoes_event_id" in cols
    end
end

@testset "rebuild_views_from_log! reproduces incremental view state" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])

        # Set up an exposure with one auto peak (so peak_excluded has a target).
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, 0.10, 1.0)", [e_id])

        # Apply a sequence of view-producing events.
        r1 = HimalayaUI.apply_event!(db, req;
            kind = "peak_added", entity_type = "exposure", entity_id = e_id,
            payload = Dict(:q => 0.20))
        @test r1.event_id > 0
        @test r1.view_row_id isa Int
        r2 = HimalayaUI.apply_event!(db, req;
            kind = "peak_excluded", entity_type = "exposure", entity_id = e_id,
            payload = Dict(:q => 0.10, :auto_peak_id => 1))
        @test r2.event_id > r1.event_id
        @test r2.view_row_id isa Int

        # Snapshot view state (sort by q for determinism).
        curations_before = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, kind, q FROM peak_curations WHERE exposure_id = ? ORDER BY q",
            [e_id]))

        # Wipe the views (simulating disaster recovery / migration).
        DBInterface.execute(db, "DELETE FROM peak_curations WHERE exposure_id = ?", [e_id])

        # Sanity: wiped.
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM peak_curations WHERE exposure_id = ?", [e_id])))

        # Rebuild from the log.
        HimalayaUI.rebuild_views_from_log!(db, e_id)

        curations_after = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, kind, q FROM peak_curations WHERE exposure_id = ? ORDER BY q",
            [e_id]))

        # Assert view state matches: same number of rows and same (kind, q) tuples.
        @test length(curations_before) == length(curations_after) == 2
        @test Set([(String(c.kind), Float64(c.q)) for c in curations_before]) ==
              Set([(String(c.kind), Float64(c.q)) for c in curations_after])
    end
end

@testset "migrate_r4_rebase_entity_type! rewrites legacy peak/index entity_type rows" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])

        # Set up minimal data.
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, 0.10, 1.0)", [e_id])
        ap_id = Int(DBInterface.lastrowid(DBInterface.execute(db, "SELECT last_insert_rowid()")))

        # Simulate a legacy user_actions row with entity_type='peak'.
        DBInterface.execute(db,
            """INSERT INTO user_actions (action, entity_type, entity_id)
               VALUES ('add_peak', 'peak', ?)""", [ap_id])

        # Run the migration.
        HimalayaUI.migrate_r4_rebase_entity_type!(db)

        # The row should now have entity_type='exposure'.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT entity_type, entity_id FROM user_actions WHERE action = 'add_peak'"))
        @test !isempty(rows)
        @test String(rows[1].entity_type) == "exposure"
        @test Int(rows[1].entity_id) == e_id

        # Idempotent: running again does nothing.
        HimalayaUI.migrate_r4_rebase_entity_type!(db)
        rows2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT entity_type FROM user_actions WHERE action = 'add_peak'"))
        @test String(rows2[1].entity_type) == "exposure"
    end
end
