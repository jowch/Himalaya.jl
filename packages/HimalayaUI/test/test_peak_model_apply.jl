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

const _PMA_TRAINED = Ref{Any}(nothing)
function _pma_trained()
    _PMA_TRAINED[] === nothing || return _PMA_TRAINED[]
    mktempdir() do dir
        b = _pt_db(dir)
        DBInterface.close!(b.db)
        ro = HimalayaUI.open_readonly_db(b.path)
        r = HimalayaUI.run_peak_training(ro; recipe = HPM.default_recipe(n_trees = 15), io = devnull,
                                         peaks_slack = 100, ff_slack = 100)
        DBInterface.close!(ro)
        _PMA_TRAINED[] = JSON3.read(JSON3.write(r))
    end
    _PMA_TRAINED[]
end

function _pma_result(; curations = nothing, watermark = nothing)
    r = Dict{String, Any}(String(k) => v for (k, v) in pairs(_pma_trained()))
    r["status"] = "passed"
    curations === nothing || (r["counts"] = merge(Dict{String, Any}(String(k) => v for (k, v) in pairs(r["counts"])),
                                                  Dict("curations" => curations)))
    watermark === nothing || (r["watermark"] = watermark)
    r
end

function _pma_add_model!(db; status = "passed", kw...)
    id = HimalayaUI.insert_peak_model!(db, _pma_result(; kw...))
    status == "passed" || DBInterface.execute(db, "UPDATE peak_models SET status = ? WHERE id = ?", [status, id])
    id
end

_pma_status(db, id) = _pma_one(db, "SELECT status, fitted_json IS NOT NULL AS has_fitted FROM peak_models WHERE id = ?", [id])

@testset "the resolver returns the active model, caches it and falls back to the hybrid" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        @test _pma_trained().status == "passed"
        @test HimalayaUI.production_detector(db) === nothing
        id = _pma_add_model!(db)
        @test HimalayaUI.production_detector(db) === nothing
        DBInterface.execute(db, "UPDATE peak_models SET status = 'active' WHERE id = ?", [id])
        m = HimalayaUI.production_detector(db)
        @test m.id == id && m.model isa HPM.TrainedModel
        @test HimalayaUI.production_detector(db) === m
        HimalayaUI.invalidate_peak_model_cache!()
        m2 = HimalayaUI.production_detector(db)
        @test m2 !== m && m2.id == id
        id2 = _pma_add_model!(db; status = "active")
        DBInterface.execute(db, "UPDATE peak_models SET status = 'retired' WHERE id = ?", [id])
        @test HimalayaUI.production_detector(db).id == id2
        DBInterface.execute(db, "UPDATE peak_models SET status = 'retired'")
        @test HimalayaUI.production_detector(db) === nothing
        DBInterface.close!(db)
    end
end

@testset "a model with another feature version, or one that cannot load, is not used" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        id = _pma_add_model!(db; status = "active")
        @test HimalayaUI.production_detector(db) !== nothing
        HimalayaUI.invalidate_peak_model_cache!()
        DBInterface.execute(db, "UPDATE peak_models SET feature_version = 'default-1;old' WHERE id = ?", [id])
        @test (@test_logs (:warn, r"peak model") HimalayaUI.production_detector(db)) === nothing
        DBInterface.execute(db, "UPDATE peak_models SET feature_version = ?, fitted_json = '{\"kind\":\"nope\"}' WHERE id = ?",
                            [_pma_trained().feature_version, id])
        HimalayaUI.invalidate_peak_model_cache!()
        @test (@test_logs (:warn, r"peak model") HimalayaUI.production_detector(db)) === nothing
        DBInterface.close!(db)
    end
end

@testset "boot retires models whose feature version is out of date" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        good = _pma_add_model!(db; status = "active")
        @test isempty(HimalayaUI.retire_stale_peak_models!(db))
        @test _pma_status(db, good).status == "active"
        stale = _pma_add_model!(db; status = "active")
        DBInterface.execute(db, "UPDATE peak_models SET feature_version = 'default-1;old' WHERE id = ?", [stale])
        @test HimalayaUI.retire_stale_peak_models!(db) == [stale]
        r = _pma_status(db, stale)
        @test r.status == "retired" && r.has_fitted == 0
        @test _pma_status(db, good).status == "active"
        DBInterface.close!(db)
    end
end

@testset "a stale model makes a retrain due without new labels, once per process" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        due() = HimalayaUI.peak_training_due(db; min_events = 20, quiet_seconds = 600.0)
        @test !due()
        id = _pma_add_model!(db; status = "retired")
        @test !due()
        for _ in 1:3
            DBInterface.execute(db,
                "INSERT INTO user_actions (action, entity_type, entity_id, timestamp) VALUES ('peak_added', 'exposure', 1, datetime('now', '-1 hour'))")
        end
        @test !due()
        DBInterface.execute(db, "UPDATE peak_models SET feature_version = 'default-1;old' WHERE id = ?", [id])
        @test due()
        DBInterface.execute(db, "INSERT INTO peak_models (created_at, status, label_watermark, core_version) VALUES ('t', 'failed', 0, 'v')")
        @test due()
        DBInterface.close!(db)
    end
end
