using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

_pp_req() = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])

function _pp_fixture(dir; exposures = 1)
    db = open_prepared_clone(dir)
    exp_id = HimalayaUI.create_experiment!(db; path = "/x", data_dir = "/x", analysis_dir = "/x")
    s_id = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "S")
    other = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "O")
    eids = [HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = s_id) for _ in 1:exposures]
    oid = HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = other)
    for e in (eids..., oid)
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, found_by, model_id) VALUES (?, 0.1, 'hybrid', NULL)", [e])
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, found_by, model_id) VALUES (?, 0.2, 'model', 3)", [e])
    end
    ix = Int(DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO indices (exposure_id, phase, basis) VALUES (?, 'Pn3m', 0.1)", [eids[1]])))
    (; db, exp_id, s_id, eids, oid, ix)
end

_reviewed(db, id) = !ismissing(only(Tables.rowtable(DBInterface.execute(db,
    "SELECT reviewed_at FROM exposures WHERE id = ?", [id]))).reviewed_at)
_reviews(db, id) = Tables.rowtable(DBInterface.execute(db,
    "SELECT * FROM peak_reviews WHERE exposure_id = ? ORDER BY id", [id]))
_n_events(db) = only(Tables.rowtable(DBInterface.execute(db, "SELECT COUNT(*) AS n FROM user_actions"))).n

@testset "reviewing events pin the exposure and snapshot its auto peaks" begin
    cases = [
        ("peak_added",           "exposure", fx -> Dict(:q => 0.15)),
        ("peak_excluded",        "exposure", fx -> Dict(:q => 0.1)),
        ("peak_unexcluded",      "exposure", fx -> Dict(:q => 0.1, :curation_qs => Float64[])),
        ("peak_removed",         "exposure", fx -> Dict(:q => 0.15)),
        ("assignment_add",       "exposure", fx -> Dict(:index_id => fx.ix)),
        ("assignment_remove",    "exposure", fx -> Dict(:index_id => fx.ix)),
        ("assignment_set_state", "exposure", fx -> Dict(:state => "null")),
        ("set_exposure_status",  "exposure", fx -> Dict(:status => "rejected")),
    ]
    for (kind, etype, payload) in cases
        mktempdir() do dir
            fx = _pp_fixture(dir)
            e = fx.eids[1]
            @test !_reviewed(fx.db, e)
            r = HimalayaUI.apply_event!(fx.db, _pp_req(); kind, entity_type = etype,
                                        entity_id = e, payload = payload(fx))
            @test _reviewed(fx.db, e)
            @test !_reviewed(fx.db, fx.oid)
            snap = only(_reviews(fx.db, e))
            @test snap.event_id == r.event_id
            @test snap.backfilled == 0
            @test snap.direct == (kind == "set_exposure_status" ? 0 : 1)
            @test snap.shown_detector == HimalayaUI.ANALYSIS_VERSION
            @test ismissing(snap.shown_model_id)
            peaks = JSON3.read(snap.peaks_json)
            @test [p.q for p in peaks] == [0.1, 0.2]
            @test [p.found_by for p in peaks] == ["hybrid", "model"]
            @test [p.model_id for p in peaks] == [nothing, 3]
        end
    end
end

@testset "set_exposure_status without a payload still pins" begin
    mktempdir() do dir
        fx = _pp_fixture(dir)
        HimalayaUI.log_action!(fx.db, _pp_req(); action = "set_exposure_status",
                               entity_type = "exposure", entity_id = fx.eids[1])
        @test _reviewed(fx.db, fx.eids[1])
        @test length(_reviews(fx.db, fx.eids[1])) == 1
    end
end

@testset "other events do not pin" begin
    mktempdir() do dir
        fx = _pp_fixture(dir)
        for (kind, etype, eid, payload) in [
                ("analyze_run", "exposure", fx.eids[1], Dict(:findpeaks_skipped => false)),
                ("label_added", "sample", fx.s_id, Dict(:label => "x", :experiment_id => fx.exp_id)),
                ("update_sample", "sample", fx.s_id, Dict(:name => "n"))]
            HimalayaUI.apply_event!(fx.db, _pp_req(); kind, entity_type = etype, entity_id = eid, payload)
        end
        @test !_reviewed(fx.db, fx.eids[1])
        @test isempty(_reviews(fx.db, fx.eids[1]))
    end
end

@testset "reviewed_at keeps the first event's time; every event adds a snapshot" begin
    mktempdir() do dir
        fx = _pp_fixture(dir)
        e = fx.eids[1]
        HimalayaUI.apply_event!(fx.db, _pp_req(); kind = "peak_added", entity_type = "exposure",
                                entity_id = e, payload = Dict(:q => 0.15))
        DBInterface.execute(fx.db, "UPDATE exposures SET reviewed_at = 'first' WHERE id = ?", [e])
        DBInterface.execute(fx.db, "UPDATE auto_peaks SET found_by = 'both' WHERE exposure_id = ? AND q = 0.2", [e])
        HimalayaUI.apply_event!(fx.db, _pp_req(); kind = "peak_excluded", entity_type = "exposure",
                                entity_id = e, payload = Dict(:q => 0.1))
        @test only(Tables.rowtable(DBInterface.execute(fx.db,
            "SELECT reviewed_at FROM exposures WHERE id = ?", [e]))).reviewed_at == "first"
        snaps = _reviews(fx.db, e)
        @test length(snaps) == 2
        @test [p.found_by for p in JSON3.read(snaps[2].peaks_json)] == ["hybrid", "both"]
    end
end

@testset "a failing event writes neither the pin nor the snapshot" begin
    mktempdir() do dir
        fx = _pp_fixture(dir)
        e = fx.eids[1]
        n = _n_events(fx.db)
        @test_throws Exception HimalayaUI.apply_event!(fx.db, _pp_req(); kind = "assignment_add",
            entity_type = "exposure", entity_id = e, payload = Dict(:nope => 1))
        @test _n_events(fx.db) == n
        @test !_reviewed(fx.db, e)
        @test isempty(_reviews(fx.db, e))
    end
end

@testset "an idempotent retry does not snapshot again" begin
    mktempdir() do dir
        fx = _pp_fixture(dir)
        e = fx.eids[1]
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice", "X-Client-Op-Id" => "op-1"], UInt8[])
        for _ in 1:2
            HimalayaUI.apply_event!(fx.db, req; kind = "peak_excluded", entity_type = "exposure",
                                    entity_id = e, payload = Dict(:q => 0.1))
        end
        @test length(_reviews(fx.db, e)) == 1
    end
end

@testset "a sample call pins every exposure of the sample and no other" begin
    mktempdir() do dir
        fx = _pp_fixture(dir; exposures = 3)
        r = HimalayaUI.apply_event!(fx.db, _pp_req(); kind = "sample_call_set_state",
            entity_type = "sample", entity_id = fx.s_id, payload = Dict(:state => "form_factor"))
        for e in fx.eids
            @test _reviewed(fx.db, e)
            @test only(_reviews(fx.db, e)).event_id == r.event_id
            @test only(_reviews(fx.db, e)).direct == 0
        end
        @test !_reviewed(fx.db, fx.oid)
        @test isempty(_reviews(fx.db, fx.oid))
    end
end

@testset "rebuild_views_from_log! leaves reviewed_at and peak_reviews alone" begin
    mktempdir() do dir
        fx = _pp_fixture(dir)
        e = fx.eids[1]
        HimalayaUI.apply_event!(fx.db, _pp_req(); kind = "peak_excluded", entity_type = "exposure",
                                entity_id = e, payload = Dict(:q => 0.1))
        before = (_reviews(fx.db, e), _reviewed(fx.db, e))
        DBInterface.execute(fx.db, "DELETE FROM peak_curations WHERE exposure_id = ?", [e])
        HimalayaUI.rebuild_views_from_log!(fx.db, e)
        @test isequal((_reviews(fx.db, e), _reviewed(fx.db, e)), before)
        @test only(Tables.rowtable(DBInterface.execute(fx.db,
            "SELECT COUNT(*) AS n FROM peak_curations WHERE exposure_id = ?", [e]))).n == 1
    end
end

@testset "deleting an experiment removes its review snapshots" begin
    mktempdir() do dir
        fx = _pp_fixture(dir)
        HimalayaUI.apply_event!(fx.db, _pp_req(); kind = "peak_excluded", entity_type = "exposure",
                                entity_id = fx.eids[1], payload = Dict(:q => 0.1))
        HimalayaUI._delete_experiment_rows!(fx.db, fx.exp_id)
        @test only(Tables.rowtable(DBInterface.execute(fx.db, "SELECT COUNT(*) AS n FROM peak_reviews"))).n == 0
    end
end

@testset "diff_update_auto_peaks! writes found_by and model_id" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        exp_id = HimalayaUI.create_experiment!(db; path = "/x", data_dir = "/x", analysis_dir = "/x")
        s = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "S")
        e = HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = s)
        I = fill(10.0, 100)
        pr(qs) = (q = qs, indices = [round(Int, 100q) for q in qs], prominence = fill(1.0, length(qs)),
                  sharpness = fill(1.0, length(qs)))
        rows() = [(r.q, r.found_by, r.model_id) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT q, found_by, model_id FROM auto_peaks WHERE exposure_id = ? ORDER BY q", [e]))]
        HimalayaUI.diff_update_auto_peaks!(db, e, pr([0.3, 0.6]), I; q_step = 0.01)
        @test isequal(rows(), [(0.3, "hybrid", missing), (0.6, "hybrid", missing)])
        HimalayaUI.diff_update_auto_peaks!(db, e, pr([0.3, 0.6]), I; q_step = 0.01,
                                           found_by = ["both", "model"], model_id = [5, 5])
        @test isequal(rows(), [(0.3, "both", 5), (0.6, "model", 5)])
        HimalayaUI.diff_update_auto_peaks!(db, e, pr([0.3, 0.8]), I; q_step = 0.01, found_by = "model", model_id = 9)
        @test isequal(rows(), [(0.3, "model", 9), (0.8, "model", 9)])
        HimalayaUI.diff_update_auto_peaks!(db, e, pr([0.3, 0.8]), I; q_step = 0.01)
        @test isequal(rows(), [(0.3, "hybrid", missing), (0.8, "hybrid", missing)])
        @test_throws Exception HimalayaUI.diff_update_auto_peaks!(db, e, pr([0.3]), I; q_step = 0.01, found_by = "bogus")
    end
end

function _pp_legacy(dir)
    path = joinpath(dir, "legacy.db")
    db = HimalayaUI.open_db(path)
    exp_id = HimalayaUI.create_experiment!(db; path = "/x", data_dir = "/x", analysis_dir = "/x")
    mk(name) = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = name)
    s_cur, s_asg, s_call, s_stat, s_none = mk("cur"), mk("asg"), mk("call"), mk("stat"), mk("none")
    ex(s) = HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = s)
    e = (cur = ex(s_cur), asg = ex(s_asg), call1 = ex(s_call), call2 = ex(s_call),
         stat = ex(s_stat), none = ex(s_none))
    for id in values(e)
        DBInterface.execute(db, "INSERT INTO auto_peaks (exposure_id, q) VALUES (?, 0.1)", [id])
    end
    DBInterface.execute(db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.3)", [e.cur])
    DBInterface.execute(db,
        "INSERT INTO user_actions (action, entity_type, entity_id, timestamp) VALUES ('peak_added', 'exposure', ?, '2026-01-02 03:04:05')", [e.cur])
    ix = Int(DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO indices (exposure_id, phase, basis) VALUES (?, 'Pn3m', 0.1)", [e.asg])))
    DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [e.asg, ix])
    DBInterface.execute(db, "INSERT INTO sample_calls (sample_id, state) VALUES (?, 'form_factor')", [s_call])
    DBInterface.execute(db, "UPDATE exposures SET status = 'rejected' WHERE id = ?", [e.stat])
    DBInterface.execute(db, "ALTER TABLE auto_peaks DROP COLUMN found_by")
    DBInterface.execute(db, "ALTER TABLE auto_peaks DROP COLUMN model_id")
    DBInterface.execute(db, "ALTER TABLE exposures DROP COLUMN reviewed_at")
    DBInterface.execute(db, "DROP TABLE peak_reviews")
    DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_PEAK_PROVENANCE])
    close(db)
    (; path, e)
end

@testset "migration backfills provenance and review snapshots, once" begin
    mktempdir() do dir
        lg = _pp_legacy(dir)
        db = HimalayaUI.open_db(lg.path)
        e = lg.e
        for id in (e.cur, e.asg, e.call1, e.call2, e.stat)
            @test _reviewed(db, id)
            snap = only(_reviews(db, id))
            @test snap.backfilled == 1 && ismissing(snap.event_id)
            @test snap.direct == (id in (e.cur, e.asg) ? 1 : 0)
            @test snap.shown_detector == HimalayaUI.ANALYSIS_VERSION
            @test [p.found_by for p in JSON3.read(snap.peaks_json)] == ["hybrid"]
        end
        @test !_reviewed(db, e.none) && isempty(_reviews(db, e.none))
        @test only(Tables.rowtable(DBInterface.execute(db,
            "SELECT reviewed_at FROM exposures WHERE id = ?", [e.cur]))).reviewed_at == "2026-01-02T03:04:05.000Z"
        @test isequal(only(Tables.rowtable(DBInterface.execute(db,
            "SELECT found_by, model_id FROM auto_peaks WHERE exposure_id = ?", [e.cur]))),
            (found_by = "hybrid", model_id = missing))
        n = only(Tables.rowtable(DBInterface.execute(db, "SELECT COUNT(*) AS n FROM peak_reviews"))).n
        close(db)
        db = HimalayaUI.open_db(lg.path)
        @test only(Tables.rowtable(DBInterface.execute(db, "SELECT COUNT(*) AS n FROM peak_reviews"))).n == n == 5
        @test only(Tables.rowtable(DBInterface.execute(db, "PRAGMA integrity_check"))).integrity_check == "ok"
        close(db)
    end
end

@testset "migration is safe when re-run after a partial apply" begin
    mktempdir() do dir
        lg = _pp_legacy(dir)
        db = SQLite.DB(lg.path)
        HimalayaUI.create_schema!(db)
        DBInterface.execute(db, "ALTER TABLE exposures ADD COLUMN reviewed_at TEXT")
        HimalayaUI.migrate_peak_provenance!(db)
        @test length(_reviews(db, lg.e.cur)) == 1
        @test HimalayaUI._migrated(db, HimalayaUI.MIGRATION_PEAK_PROVENANCE)
        close(db)
    end
end
