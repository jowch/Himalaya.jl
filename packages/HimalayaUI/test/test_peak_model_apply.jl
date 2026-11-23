using Test, HTTP, JSON3, Random, SQLite, DBInterface, Tables
using HimalayaUI
import HimalayaPeakModel as HPM
isdefined(@__MODULE__, :_pt_db) || include("peak_training_fixture.jl")
isdefined(@__MODULE__, :_pma_db) || include("peak_model_apply_fixture.jl")


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

@testset "an unreviewed exposure gets the model's peaks with provenance" begin
    mktempdir() do dir
        b = _pma_db(dir)
        model_id = _pma_add_model!(b.db; status = "active")
        pm = HimalayaUI.production_detector(b.db)
        e = b.eids[1]
        HimalayaUI.analyze_exposure!(b.db, e, b.adir)
        q, I, σ = b.traces[e]
        expected = HPM.predict_peaks(pm.model, q, I, σ)
        rows = _pma_peaks(b.db, e)
        @test !isempty(expected)
        @test [r.q for r in rows] == expected
        hyb = Himalaya.findpeaks(q, I, σ; method = :hybrid)
        step = HimalayaUI.grid_step(q)
        for r in rows
            near = any(h -> abs(r.q - h) <= HimalayaUI.peak_match_tol(r.q, step), hyb.q)
            @test r.found_by == (near ? "both" : "model")
            @test r.model_id == model_id
            @test r.sharpness > 0 && r.prominence >= 0
            @test q[r.findpeaks_index] == r.q
        end
        v = _pma_version(b.db, e)
        @test v.v == model_id && v.analysis_version == HimalayaUI.ANALYSIS_VERSION
        DBInterface.close!(b.db)
    end
end

@testset "reviewed exposures, and databases without a model, get the hybrid's peaks" begin
    mktempdir() do dir
        b = _pma_db(dir)
        a, r = b.eids[1], b.eids[2]
        HimalayaUI.analyze_exposure!(b.db, a, b.adir)
        @test all(p -> p.found_by == "hybrid" && ismissing(p.model_id), _pma_peaks(b.db, a))
        @test ismissing(_pma_version(b.db, a).v)
        _pma_add_model!(b.db; status = "active")
        _pma_review!(b.db, r)
        HimalayaUI.analyze_exposure!(b.db, r, b.adir)
        q, I, σ = b.traces[r]
        hyb = Himalaya.findpeaks(q, I, σ; method = :hybrid)
        @test [p.q for p in _pma_peaks(b.db, r)] == hyb.q
        @test all(p -> p.found_by == "hybrid" && ismissing(p.model_id), _pma_peaks(b.db, r))
        @test ismissing(_pma_version(b.db, r).v)
        DBInterface.close!(b.db)
    end
end

@testset "a reviewed exposure re-detected after a version bump goes back to the hybrid's provenance" begin
    mktempdir() do dir
        b = _pma_db(dir)
        e = b.eids[1]
        id = _pma_add_model!(b.db; status = "active")
        HimalayaUI.analyze_exposure!(b.db, e, b.adir)
        @test _pma_version(b.db, e).v == id
        _pma_review!(b.db, e)
        DBInterface.execute(b.db, "UPDATE exposures SET analysis_version = 'older' WHERE id = ?", [e])
        HimalayaUI.analyze_exposure!(b.db, e, b.adir)
        @test all(p -> p.found_by == "hybrid" && ismissing(p.model_id), _pma_peaks(b.db, e))
        @test ismissing(_pma_version(b.db, e).v)
        DBInterface.close!(b.db)
    end
end

@testset "a failing prediction, or a trace too short for the model, keeps the hybrid's peaks" begin
    mktempdir() do dir
        b = _pma_db(dir)
        _pma_add_model!(b.db; status = "active")
        old = HimalayaUI.PEAK_MODEL_PREDICT[]
        HimalayaUI.PEAK_MODEL_PREDICT[] = (args...; kw...) -> error("boom")
        try
            e = b.eids[1]
            @test_logs (:warn, r"peak model") HimalayaUI.analyze_exposure!(b.db, e, b.adir)
            q, I, σ = b.traces[e]
            @test [p.q for p in _pma_peaks(b.db, e)] == Himalaya.findpeaks(q, I, σ; method = :hybrid).q
            @test all(p -> p.found_by == "hybrid" && ismissing(p.model_id), _pma_peaks(b.db, e))
            @test ismissing(_pma_version(b.db, e).v)
            @test ismissing(_pma_one(b.db, "SELECT analysis_error AS e FROM exposures WHERE id = ?", [e]).e)
        finally
            HimalayaUI.PEAK_MODEL_PREDICT[] = old
        end
        DBInterface.close!(b.db)
    end
    mktempdir() do dir
        b = _pma_db(dir; n = 1, points = 100)
        _pma_add_model!(b.db; status = "active")
        e = b.eids[1]
        HimalayaUI.analyze_exposure!(b.db, e, b.adir)
        @test all(p -> p.found_by == "hybrid" && ismissing(p.model_id), _pma_peaks(b.db, e))
        @test ismissing(_pma_version(b.db, e).v)
        DBInterface.close!(b.db)
    end
end

@testset "force_detect re-detects despite matching hashes; only_unreviewed leaves reviewed exposures alone" begin
    mktempdir() do dir
        b = _pma_db(dir)
        a, r = b.eids[1], b.eids[2]
        for e in (a, r)
            HimalayaUI.analyze_exposure!(b.db, e, b.adir)
        end
        before = _pma_peaks(b.db, a)
        id = _pma_add_model!(b.db; status = "active")
        HimalayaUI.analyze_exposure!(b.db, a, b.adir)
        @test isequal(_pma_peaks(b.db, a), before) && ismissing(_pma_version(b.db, a).v)
        _pma_review!(b.db, r)
        before_r = _pma_peaks(b.db, r)
        HimalayaUI.analyze_exposure!(b.db, r, b.adir; force_detect = true, only_unreviewed = true)
        @test isequal(_pma_peaks(b.db, r), before_r) && ismissing(_pma_version(b.db, r).v)
        HimalayaUI.analyze_exposure!(b.db, a, b.adir; force_detect = true, only_unreviewed = true)
        @test _pma_version(b.db, a).v == id
        @test all(p -> p.model_id == id, _pma_peaks(b.db, a))
        DBInterface.close!(b.db)
    end
end

@testset "a passed run becomes active when nothing is active" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        id = _pma_add_model!(db; curations = 100)
        @test HimalayaUI.promote_peak_model!(db, id)
        @test _pma_status(db, id) == (status = "active", has_fitted = 1)
        @test HimalayaUI.production_detector(db).id == id
        @test !HimalayaUI.promote_peak_model!(db, id)
        DBInterface.close!(db)
    end
end

@testset "a later run replaces the active model only after 25% more curations" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        first = _pma_add_model!(db; curations = 100)
        @test HimalayaUI.promote_peak_model!(db, first)
        cached = HimalayaUI.production_detector(db)
        small = _pma_add_model!(db; curations = 124)
        @test !HimalayaUI.promote_peak_model!(db, small)
        @test _pma_status(db, small) == (status = "passed", has_fitted = 1)
        @test _pma_status(db, first).status == "active"
        @test HimalayaUI.production_detector(db) === cached
        grown = _pma_add_model!(db; curations = 125)
        @test _pma_status(db, small) == (status = "passed", has_fitted = 0)
        @test HimalayaUI.promote_peak_model!(db, grown)
        @test _pma_status(db, grown) == (status = "active", has_fitted = 1)
        @test _pma_status(db, first) == (status = "retired", has_fitted = 0)
        @test HimalayaUI.production_detector(db).id == grown
        DBInterface.close!(db)
    end
end

@testset "growth is measured from the active model, and an unusable active model does not block" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        a = _pma_add_model!(db; curations = 100)
        HimalayaUI.promote_peak_model!(db, a)
        b = _pma_add_model!(db; curations = 130)
        @test HimalayaUI.promote_peak_model!(db, b)
        c = _pma_add_model!(db; curations = 140)
        @test !HimalayaUI.promote_peak_model!(db, c)
        DBInterface.execute(db, "UPDATE peak_models SET feature_version = 'default-1;old' WHERE id = ?", [b])
        HimalayaUI.invalidate_peak_model_cache!()
        @test @test_logs (:warn, r"peak model") HimalayaUI.promote_peak_model!(db, c)
        @test _pma_status(db, b).status == "retired" && _pma_status(db, c).status == "active"
        DBInterface.close!(db)
    end
end

@testset "a run whose model cannot be loaded stays passed, as does a rejected one" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        id = HimalayaUI.insert_peak_model!(db, Dict{String, Any}("status" => "passed", "watermark" => 1,
            "counts" => Dict("curations" => 1), "recipe" => Dict("version" => 1), "recipe_id" => "r",
            "feature_version" => "fv", "core_version" => "0", "threshold" => 0.5, "nms" => 8,
            "metrics" => Dict("gate" => Dict("passed" => true)), "fitted" => Dict("kind" => "fake")))
        @test (@test_logs (:warn, r"peak model") HimalayaUI.promote_peak_model!(db, id)) == false
        @test _pma_status(db, id).status == "passed"
        rej = HimalayaUI.insert_peak_model!(db, merge(_pma_result(), Dict{String, Any}("status" => "rejected")))
        @test !HimalayaUI.promote_peak_model!(db, rej)
        @test HimalayaUI.production_detector(db) === nothing
        DBInterface.close!(db)
    end
end

@testset "labels from a model-era snapshot: model-only peaks are dropped, shared ones count as kept" begin
    mktempdir() do dir
        b = _pma_db(dir; n = 2)
        _pma_add_model!(b.db; status = "active")
        e, hybrid_only = b.eids
        HimalayaUI.analyze_exposure!(b.db, e, b.adir)
        HimalayaUI.analyze_exposure!(b.db, hybrid_only, b.adir)
        DBInterface.execute(b.db, "UPDATE auto_peaks SET found_by = 'hybrid', model_id = NULL WHERE exposure_id = ?", [hybrid_only])
        DBInterface.execute(b.db, "UPDATE exposures SET peak_model_version = NULL WHERE id = ?", [hybrid_only])
        peaks = _pma_peaks(b.db, e)
        DBInterface.execute(b.db, "UPDATE auto_peaks SET found_by = 'model' WHERE exposure_id = ? AND q = ?", [e, peaks[1].q])
        DBInterface.execute(b.db, "UPDATE auto_peaks SET found_by = 'both' WHERE exposure_id = ? AND q != ?", [e, peaks[1].q])
        for x in (e, hybrid_only)
            HimalayaUI.apply_event!(b.db, _pma_req(); kind = "peak_added", entity_type = "exposure",
                                    entity_id = x, payload = Dict(:q => 0.321), defer_broadcast = true)
        end
        exps = Dict(x.id => x for x in HimalayaUI.labeled_exposures(b.db; io = devnull))
        @test exps[e].kept_model == [peaks[1].q]
        @test sort(exps[e].kept) == [p.q for p in peaks[2:end]]
        @test isempty(exps[hybrid_only].kept_model) && !isempty(exps[hybrid_only].kept)
        @test HimalayaUI.no_model_shown_ids(b.db) == Set([hybrid_only])
        DBInterface.close!(b.db)
    end
end
