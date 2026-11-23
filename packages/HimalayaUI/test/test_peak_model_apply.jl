using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI
import HimalayaPeakModel as HPM
isdefined(@__MODULE__, :_pt_db) || include("peak_training_fixture.jl")

_pma_req() = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
_pma_one(db, sql, args = []) = only(Tables.rowtable(DBInterface.execute(db, sql, args)))

@testset "peak_model_version: column, migration, wire" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        @test "peak_model_version" in HimalayaUI.cols_of(db, "exposures")
        DBInterface.execute(db, "ALTER TABLE exposures DROP COLUMN peak_model_version")
        HimalayaUI.migrate_schema!(db)
        HimalayaUI.migrate_schema!(db)
        @test "peak_model_version" in HimalayaUI.cols_of(db, "exposures")
        xid = HimalayaUI.create_experiment!(db; path = "/x", data_dir = "/x", analysis_dir = "/x")
        sid = HimalayaUI.create_sample!(db; experiment_id = xid, name = "S")
        eid = HimalayaUI.create_exposure!(db; experiment_id = xid, sample_id = sid)
        DBInterface.execute(db, "UPDATE exposures SET peak_model_version = 7 WHERE id = ?", [eid])
        with_inproc_routes(db) do call
            r = call("GET", "/api/samples/$sid/exposures")
            body = JSON3.read(String(r.body))
            @test length(body) == 1 && !haskey(body[1], :peak_model_version)
        end
        DBInterface.close!(db)
    end
end

@testset "review snapshots record the model on screen" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        xid = HimalayaUI.create_experiment!(db; path = "/x", data_dir = "/x", analysis_dir = "/x")
        sid = HimalayaUI.create_sample!(db; experiment_id = xid, name = "S")
        a = HimalayaUI.create_exposure!(db; experiment_id = xid, sample_id = sid)
        b = HimalayaUI.create_exposure!(db; experiment_id = xid, sample_id = sid)
        DBInterface.execute(db, "UPDATE exposures SET peak_model_version = 7 WHERE id = ?", [a])
        for e in (a, b)
            HimalayaUI.apply_event!(db, _pma_req(); kind = "peak_added", entity_type = "exposure",
                                    entity_id = e, payload = Dict(:q => 0.15), defer_broadcast = true)
        end
        shown(e) = _pma_one(db, "SELECT shown_model_id FROM peak_reviews WHERE exposure_id = ?", [e]).shown_model_id
        @test shown(a) == 7
        @test ismissing(shown(b))
        DBInterface.close!(db)
    end
end
