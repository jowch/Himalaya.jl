const REVIEW_EXPOSURE_KINDS = ("peak_added", "peak_removed", "peak_excluded", "peak_unexcluded",
    "assignment_add", "assignment_remove", "assignment_set_state", "set_exposure_status")
const DIRECT_REVIEW_KINDS = ("peak_added", "peak_removed", "peak_excluded", "peak_unexcluded",
    "assignment_add", "assignment_remove", "assignment_set_state")
const REVIEW_SAMPLE_KINDS = ("sample_call_set_state",)

function review_exposure_ids(db::SQLite.DB, kind::AbstractString,
                             entity_type::AbstractString, entity_id::Integer)
    if kind in REVIEW_EXPOSURE_KINDS && entity_type == "exposure"
        return [Int(entity_id)]
    elseif kind in REVIEW_SAMPLE_KINDS && entity_type == "sample"
        return Int[r.id for r in sqlrows(DBInterface.execute(db,
            "SELECT id FROM exposures WHERE sample_id = ? ORDER BY id", [Int(entity_id)]))]
    end
    Int[]
end

function _snapshot_peaks_json(db::SQLite.DB, exposure_id::Integer)
    rows = sqlrows(DBInterface.execute(db,
        "SELECT q, found_by, model_id FROM auto_peaks WHERE exposure_id = ? ORDER BY id",
        [Int(exposure_id)]))
    JSON3.write([Dict(:q => Float64(r.q), :found_by => String(r.found_by),
                      :model_id => ismissing(r.model_id) ? nothing : Int(r.model_id))
                 for r in rows])
end

function record_review!(db::SQLite.DB, exposure_id::Integer;
                        event_id::Union{Int,Nothing} = nothing,
                        at::AbstractString = now_iso(),
                        backfilled::Bool = false,
                        direct::Bool = false)
    DBInterface.execute(db,
        "UPDATE exposures SET reviewed_at = ? WHERE id = ? AND reviewed_at IS NULL",
        [String(at), Int(exposure_id)])
    DBInterface.execute(db,
        """INSERT INTO peak_reviews
               (exposure_id, event_id, created_at, shown_detector, shown_model_id, backfilled, direct, peaks_json)
           VALUES (?, ?, ?, ?, NULL, ?, ?, ?)""",
        [Int(exposure_id), event_id, String(at), ANALYSIS_VERSION, backfilled ? 1 : 0,
         direct ? 1 : 0, _snapshot_peaks_json(db, exposure_id)])
    nothing
end

const NO_REVIEW_PINS_KEY = :himalaya_no_review_pins

function without_review_pins(f)
    task_local_storage(f, NO_REVIEW_PINS_KEY, true)
end

function record_event_reviews!(db::SQLite.DB, kind::AbstractString,
                               entity_type::AbstractString, entity_id::Integer, event_id::Integer;
                               direct::Union{Bool, Nothing} = nothing)
    get(task_local_storage(), NO_REVIEW_PINS_KEY, false) === true && return nothing
    is_direct = direct === nothing ? (kind in DIRECT_REVIEW_KINDS && entity_type == "exposure") : direct
    for id in review_exposure_ids(db, kind, entity_type, entity_id)
        record_review!(db, id; event_id = Int(event_id), direct = is_direct)
    end
    nothing
end

const LEGACY_DIRECT_KINDS = ("index_confirmed", "index_unconfirmed")
_sql_list(kinds) = join(("'$k'" for k in kinds), ", ")
const _REVIEW_EVENT_KIND_SQL = _sql_list((REVIEW_EXPOSURE_KINDS..., REVIEW_SAMPLE_KINDS...))
const _DIRECT_KIND_SQL = _sql_list((DIRECT_REVIEW_KINDS..., LEGACY_DIRECT_KINDS...))

function migrate_peak_provenance_schema!(db::SQLite.DB)
    "found_by" in cols_of(db, "auto_peaks") || DBInterface.execute(db,
        """ALTER TABLE auto_peaks ADD COLUMN found_by TEXT NOT NULL DEFAULT 'hybrid'
           CHECK (found_by IN ('hybrid', 'model', 'both'))""")
    "model_id" in cols_of(db, "auto_peaks") ||
        DBInterface.execute(db, "ALTER TABLE auto_peaks ADD COLUMN model_id INTEGER")
    "reviewed_at" in cols_of(db, "exposures") ||
        DBInterface.execute(db, "ALTER TABLE exposures ADD COLUMN reviewed_at TEXT")
    "direct" in cols_of(db, "peak_reviews") ||
        DBInterface.execute(db, "ALTER TABLE peak_reviews ADD COLUMN direct INTEGER NOT NULL DEFAULT 0")
    nothing
end

function migrate_peak_provenance!(db::SQLite.DB)
    _migrated(db, MIGRATION_PEAK_PROVENANCE) && return nothing
    SQLite.transaction(db) do
        migrate_peak_provenance_schema!(db)
        migrated_at = now_iso()
        rows = sqlrows(DBInterface.execute(db, """
            SELECT e.id,
                   (EXISTS (SELECT 1 FROM peak_curations c WHERE c.exposure_id = e.id)
                    OR EXISTS (SELECT 1 FROM user_actions ua
                                WHERE ua.entity_type = 'exposure' AND ua.entity_id = e.id
                                  AND ua.action IN ($_DIRECT_KIND_SQL))) AS direct,
                   (SELECT strftime('%Y-%m-%dT%H:%M:%fZ', MIN(ua.timestamp)) FROM user_actions ua
                     WHERE ua.action IN ($_REVIEW_EVENT_KIND_SQL, 'index_confirmed', 'index_unconfirmed')
                       AND ((ua.entity_type = 'exposure' AND ua.entity_id = e.id)
                         OR (ua.entity_type = 'sample' AND ua.entity_id = e.sample_id))) AS first_event
              FROM exposures e
             WHERE e.reviewed_at IS NULL
               AND (EXISTS (SELECT 1 FROM peak_curations c WHERE c.exposure_id = e.id)
                 OR EXISTS (SELECT 1 FROM assignment_members m WHERE m.exposure_id = e.id)
                 OR EXISTS (SELECT 1 FROM user_actions ua
                             WHERE ua.entity_type = 'exposure' AND ua.entity_id = e.id
                               AND ua.action IN ($_DIRECT_KIND_SQL))
                 OR e.status IS NOT NULL
                 OR EXISTS (SELECT 1 FROM sample_calls sc WHERE sc.sample_id = e.sample_id))
             ORDER BY e.id"""))
        for r in rows
            record_review!(db, Int(r.id);
                           at = ismissing(r.first_event) ? migrated_at : String(r.first_event),
                           backfilled = true, direct = r.direct == 1)
        end
        _record_migration!(db, MIGRATION_PEAK_PROVENANCE)
    end
    nothing
end
