using Test, SQLite, DBInterface, Tables
import JSON
using HimalayaUI

_pm_rows(db) = Tables.rowtable(DBInterface.execute(db, "SELECT * FROM peak_models ORDER BY id"))

function _pm_event!(db, kind; age = "-1 hour", entity = "exposure")
    DBInterface.execute(db,
        "INSERT INTO user_actions (action, entity_type, entity_id, timestamp) VALUES (?, ?, 1, datetime('now', ?))",
        [kind, entity, age])
    Int(DBInterface.lastrowid(DBInterface.execute(db, "SELECT last_insert_rowid()")))
end

_pm_result(; status = "passed", watermark = 5, fitted = Dict("kind" => "rf", "trees" => [1, 2]),
           reasons = String[]) = Dict{String,Any}(
    "status" => status, "watermark" => watermark,
    "counts" => Dict("curations" => 120),
    "recipe" => Dict("version" => 1), "recipe_id" => "abc123", "feature_version" => HimalayaUI.current_peak_feature_version(),
    "core_version" => "0.5.3", "threshold" => 0.55, "nms" => 8,
    "metrics" => Dict("gate" => Dict("passed" => status == "passed", "reasons" => reasons)),
    "fitted" => fitted)

@testset "peak_models table" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        cols = HimalayaUI.cols_of(db, "peak_models")
        for c in ("id", "created_at", "status", "label_watermark", "counts_json", "recipe_json", "recipe_id",
                  "feature_version", "core_version", "threshold", "nms", "metrics_json", "fitted_json")
            @test c in cols
        end
        HimalayaUI.migrate_peak_models!(db)
        HimalayaUI.migrate_peak_models!(db)
        @test isempty(_pm_rows(db))
        @test_throws SQLite.SQLiteException DBInterface.execute(db,
            "INSERT INTO peak_models (created_at, status, label_watermark, core_version) VALUES ('t', 'bogus', 0, 'v')")
        for s in ("passed", "rejected", "active", "retired", "failed")
            DBInterface.execute(db,
                "INSERT INTO peak_models (created_at, status, label_watermark, core_version) VALUES ('t', ?, 0, 'v')", [s])
        end
        @test length(_pm_rows(db)) == 5

        DBInterface.execute(db, "DROP TABLE peak_models")
        HimalayaUI.migrate_schema!(db)
        @test HimalayaUI._table_exists(db, "peak_models")
        DBInterface.close!(db)
    end
end

@testset "label events and dueness" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        due(; kw...) = HimalayaUI.peak_training_due(db; min_events = 3, quiet_seconds = 60.0, kw...)
        @test HimalayaUI.label_watermark(db) == 0
        @test !due()

        _pm_event!(db, "peak_added"); _pm_event!(db, "peak_excluded")
        @test !due()
        for k in ("set_exposure_status", "sample_call_set_state", "label_added", "analyze_run")
            _pm_event!(db, k)
        end
        _pm_event!(db, "sample_call_set_state"; entity = "sample")
        @test !due()
        @test HimalayaUI.label_watermark(db) == 2
        last = _pm_event!(db, "assignment_set_state")
        @test due()
        @test HimalayaUI.label_watermark(db) == last

        _pm_event!(db, "peak_removed"; age = "-5 seconds")
        @test !due()
        @test due(quiet_seconds = 1.0)

        w = HimalayaUI.label_watermark(db)
        HimalayaUI.insert_peak_model!(db, _pm_result(status = "rejected", watermark = w, fitted = nothing))
        @test HimalayaUI.last_training_watermark(db) == w
        @test !due(quiet_seconds = 1.0)
        for _ in 1:2; _pm_event!(db, "peak_unexcluded"); end
        @test !due(quiet_seconds = 1.0)
        _pm_event!(db, "peak_unexcluded")
        @test due(quiet_seconds = 1.0)

        w2 = HimalayaUI.label_watermark(db)
        HimalayaUI.insert_failed_peak_model!(db; watermark = w2, error = "boom")
        @test HimalayaUI.last_training_watermark(db) == w2
        @test !due(quiet_seconds = 1.0)
        DBInterface.close!(db)
    end
end

@testset "insert_peak_model! stores the fitted model only where the design says" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        fitted(id) = only(filter(r -> r.id == id, _pm_rows(db))).fitted_json
        a = HimalayaUI.insert_peak_model!(db, _pm_result())
        row = only(_pm_rows(db))
        @test row.status == "passed" && row.label_watermark == 5 && row.recipe_id == "abc123"
        @test row.threshold == 0.55 && row.nms == 8 && row.core_version == "0.5.3"
        @test JSON.parse(row.fitted_json).kind == "rf"
        @test JSON.parse(row.counts_json).curations == 120
        @test JSON.parse(row.metrics_json).gate.passed
        @test JSON.parse(row.recipe_json).version == 1

        b = HimalayaUI.insert_peak_model!(db, _pm_result(watermark = 9))
        @test ismissing(fitted(a)) && !ismissing(fitted(b))

        DBInterface.execute(db, "UPDATE peak_models SET status = 'active' WHERE id = ?", [b])
        c = HimalayaUI.insert_peak_model!(db, _pm_result(watermark = 12))
        @test !ismissing(fitted(b)) && !ismissing(fitted(c))

        r = HimalayaUI.insert_peak_model!(db, _pm_result(status = "rejected", watermark = 14,
            reasons = ["peaks per exposure too high"]))
        @test ismissing(fitted(r)) && !ismissing(fitted(c))
        @test JSON.parse(only(filter(x -> x.id == r, _pm_rows(db))).metrics_json).gate.reasons == ["peaks per exposure too high"]

        f = HimalayaUI.insert_failed_peak_model!(db; watermark = 20, error = "x" ^ 5000)
        row = only(filter(x -> x.id == f, _pm_rows(db)))
        @test row.status == "failed" && row.label_watermark == 20 && ismissing(row.fitted_json)
        @test length(JSON.parse(row.metrics_json).error) == 2000
        @test !ismissing(fitted(c))
        @test_throws ArgumentError HimalayaUI.insert_peak_model!(db, _pm_result(status = "active"))
        DBInterface.close!(db)
    end
end
