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

current_peak_feature_version() =
    HimalayaPeakModel.feature_version(HimalayaPeakModel.default_recipe().features)

function peak_training_stale(db::SQLite.DB)
    rows = sqlrows(DBInterface.execute(db,
        "SELECT feature_version FROM peak_models WHERE status != 'failed' ORDER BY id DESC LIMIT 1"))
    !isempty(rows) && coalesce(rows[1].feature_version, "") != current_peak_feature_version() &&
        label_watermark(db) > 0
end

function peak_training_due(db::SQLite.DB; min_events::Integer = PEAK_TRAINING_MIN_EVENTS,
                           quiet_seconds::Real = PEAK_TRAINING_QUIET_SECONDS,
                           after::Integer = 0, stale_attempted::Bool = false)
    stale_attempted || !peak_training_stale(db) || return true
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

struct ProductionModel
    id::Int
    model::HimalayaPeakModel.TrainedModel
end

const PEAK_MODEL_CACHE = Dict{String, Tuple{Int, Union{Nothing, ProductionModel}}}()
const PEAK_MODEL_CACHE_LOCK = ReentrantLock()

invalidate_peak_model_cache!() = lock(() -> empty!(PEAK_MODEL_CACHE), PEAK_MODEL_CACHE_LOCK)

function trained_model_from_row(row)::HimalayaPeakModel.TrainedModel
    recipe = HimalayaPeakModel.from_dict(HimalayaPeakModel.Recipe, JSON3.read(String(row.recipe_json)))
    row.feature_version == HimalayaPeakModel.feature_version(recipe.features) ||
        error("the model was trained with feature version $(row.feature_version), not the current one")
    fitted = HimalayaPeakModel.from_dict(HimalayaPeakModel.Fitted, JSON3.read(String(row.fitted_json)))
    counts = ismissing(row.counts_json) ? nothing : JSON3.read(String(row.counts_json))
    trained_on = (positive = counts === nothing ? 0 : Int(get(counts, :positive, 0)),
                  negative = counts === nothing ? 0 : Int(get(counts, :negative, 0)))
    HimalayaPeakModel.TrainedModel(recipe, fitted, recipe.threshold, recipe.nms,
        HimalayaPeakModel.feature_names(recipe.features), trained_on, String(row.feature_version))
end

function _load_production_model(db::SQLite.DB, id::Int)::Union{Nothing, ProductionModel}
    row = only(sqlrows(DBInterface.execute(db,
        "SELECT recipe_json, feature_version, fitted_json, counts_json FROM peak_models WHERE id = ?", [id])))
    (ismissing(row.fitted_json) || ismissing(row.recipe_json)) && return nothing
    try
        ProductionModel(id, trained_model_from_row(row))
    catch err
        @warn "peak model $id cannot be used; peaks come from the hybrid detector" exception = err
        nothing
    end
end

function production_detector(db::SQLite.DB)::Union{Nothing, ProductionModel}
    rows = sqlrows(DBInterface.execute(db,
        "SELECT id FROM peak_models WHERE status = 'active' ORDER BY id DESC LIMIT 1"))
    isempty(rows) && return nothing
    id = Int(rows[1].id)
    _poolable(db.file) || return _load_production_model(db, id)
    key = String(db.file)
    lock(PEAK_MODEL_CACHE_LOCK) do
        hit = get(PEAK_MODEL_CACHE, key, nothing)
        hit !== nothing && hit[1] == id && return hit[2]
        loaded = _load_production_model(db, id)
        PEAK_MODEL_CACHE[key] = (id, loaded)
        loaded
    end
end

function retire_stale_peak_models!(db::SQLite.DB)
    lock(_DB_WRITE_LOCK) do
        retired = Int[]
        for r in sqlrows(DBInterface.execute(db,
                "SELECT id, recipe_json, feature_version FROM peak_models WHERE status = 'active' ORDER BY id"))
            usable = try
                !ismissing(r.recipe_json) && r.feature_version == HimalayaPeakModel.feature_version(
                    HimalayaPeakModel.from_dict(HimalayaPeakModel.Recipe, JSON3.read(String(r.recipe_json))).features)
            catch
                false
            end
            usable && continue
            DBInterface.execute(db, "UPDATE peak_models SET status = 'retired', fitted_json = NULL WHERE id = ?", [Int(r.id)])
            push!(retired, Int(r.id))
        end
        isempty(retired) || invalidate_peak_model_cache!()
        retired
    end
end

function insert_failed_peak_model!(db::SQLite.DB; watermark::Integer, error::AbstractString)
    _insert_peak_model_row!(db; status = "failed", watermark,
        metrics = Dict("error" => first(String(error), 2000)))
end
