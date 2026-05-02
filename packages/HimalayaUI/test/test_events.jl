using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

@testset "apply_event! writes to user_actions and returns event_id" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x",
            ["X-Username" => "alice"], UInt8[])

        eid = HimalayaUI.apply_event!(db, req;
            kind = "test_kind", entity_type = "exposure", entity_id = 42,
            payload = Dict(:foo => "bar"))
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
        eid = HimalayaUI.apply_event!(db, req;
            kind = "no_payload_event", entity_type = "exposure", entity_id = 1)
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT payload FROM user_actions WHERE id = ?", [eid])))
        @test ismissing(row.payload)
    end
end

@testset "apply_event! with no X-Username writes user_id = NULL" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", Pair{String,String}[], UInt8[])
        eid = HimalayaUI.apply_event!(db, req;
            kind = "anon_event", entity_type = "exposure", entity_id = 1,
            payload = Dict(:k => "v"))
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
