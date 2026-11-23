using JSON3

const PEAK_TRAINING_MIN_EVENTS = 20
const PEAK_TRAINING_QUIET_SECONDS = 600.0
const PEAK_MODEL_STATUSES = ("passed", "rejected", "active", "retired", "failed")
const _LABEL_EVENT_SQL = "action IN ($(_sql_list(DIRECT_REVIEW_KINDS))) AND entity_type = 'exposure'"

function migrate_peak_models!(db::SQLite.DB)
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS peak_models (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            created_at      TEXT NOT NULL,
            status          TEXT NOT NULL CHECK (status IN ($(join(("'$s'" for s in PEAK_MODEL_STATUSES), ", ")))),
            label_watermark INTEGER NOT NULL,
            counts_json     TEXT,
            recipe_json     TEXT,
            recipe_id       TEXT,
            feature_version TEXT,
            core_version    TEXT NOT NULL,
            threshold       REAL,
            nms             INTEGER,
            metrics_json    TEXT,
            fitted_json     TEXT
        )""")
    nothing
end

label_watermark(db::SQLite.DB) = Int(coalesce(only(sqlrows(DBInterface.execute(db,
    "SELECT MAX(id) AS w FROM user_actions WHERE $_LABEL_EVENT_SQL"))).w, 0))

last_training_watermark(db::SQLite.DB) = Int(coalesce(only(sqlrows(DBInterface.execute(db,
    "SELECT MAX(label_watermark) AS w FROM peak_models"))).w, 0))

function peak_training_due(db::SQLite.DB; min_events::Integer = PEAK_TRAINING_MIN_EVENTS,
                           quiet_seconds::Real = PEAK_TRAINING_QUIET_SECONDS,
                           after::Integer = 0)
    r = only(sqlrows(DBInterface.execute(db, """
        SELECT COUNT(*) AS n,
               (julianday('now') - julianday(MAX(timestamp))) * 86400.0 AS idle
          FROM user_actions WHERE $_LABEL_EVENT_SQL AND id > ?""", [max(last_training_watermark(db), Int(after))])))
    r.n >= min_events && r.n > 0 && r.idle >= quiet_seconds
end

_json_or_null(x) = x === nothing ? nothing : JSON3.write(x)

function _insert_peak_model_row!(db::SQLite.DB; status, watermark, counts = nothing, recipe = nothing,
                                 recipe_id = nothing, feature_version = nothing,
                                 core_version = string(pkgversion(Himalaya)), threshold = nothing,
                                 nms = nothing, metrics = nothing, fitted = nothing)
    result = DBInterface.execute(db, """
        INSERT INTO peak_models (created_at, status, label_watermark, counts_json, recipe_json, recipe_id,
                                 feature_version, core_version, threshold, nms, metrics_json, fitted_json)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        [now_iso(), status, Int(watermark), _json_or_null(counts), _json_or_null(recipe), recipe_id,
         feature_version, core_version, threshold, nms, _json_or_null(metrics), _json_or_null(fitted)])
    Int(DBInterface.lastrowid(result))
end

function insert_peak_model!(db::SQLite.DB, r::AbstractDict)
    status = String(r["status"])
    status in ("passed", "rejected") ||
        throw(ArgumentError("a run is recorded as passed or rejected, not \"$status\""))
    get(r, "metrics", nothing) !== nothing || throw(ArgumentError("a run needs metrics"))
    get(r, "core_version", nothing) !== nothing || throw(ArgumentError("a run needs its core version"))
    SQLite.transaction(db) do
        status == "passed" && DBInterface.execute(db,
            "UPDATE peak_models SET fitted_json = NULL WHERE status = 'passed'")
        _insert_peak_model_row!(db; status, watermark = r["watermark"], counts = get(r, "counts", nothing),
            recipe = get(r, "recipe", nothing), recipe_id = get(r, "recipe_id", nothing),
            feature_version = get(r, "feature_version", nothing), core_version = String(r["core_version"]),
            threshold = get(r, "threshold", nothing), nms = get(r, "nms", nothing),
            metrics = r["metrics"], fitted = status == "passed" ? get(r, "fitted", nothing) : nothing)
    end
end

function insert_failed_peak_model!(db::SQLite.DB; watermark::Integer, error::AbstractString)
    _insert_peak_model_row!(db; status = "failed", watermark,
        metrics = Dict("error" => first(String(error), 2000)))
end
