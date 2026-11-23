using DBInterface, SQLite

function trace_median(I::AbstractVector{<:Real})
    isempty(I) && return nothing
    s = sort(I)
    n = length(s)
    isodd(n) ? Float64(s[(n + 1) ÷ 2]) : (Float64(s[n ÷ 2]) + Float64(s[n ÷ 2 + 1])) / 2
end

_known(x) = !(x === nothing || x === missing)

# physics: a _tot.dat holds raw detector counts per pixel. A beam-off image records
# physics: only the dark floor, so its median is ~0, while even a dim real image keeps
# physics: ~10 counts. The detector decides when there is a trace, whatever I₀ reads,
# physics: because a frame with a negative I₀ can still hold real scattering. Without a
# physics: trace, an I₀ at or below zero (the electrometer's offset with the shutter
# physics: closed) is the only evidence left, so there it decides alone.
const DARK_MEDIAN_COUNTS = 1.0

function is_beam_off(median, analysis_error, incident_monitor)::Bool
    _known(median) && return Float64(median) < DARK_MEDIAN_COUNTS
    _known(analysis_error) && String(analysis_error) == NO_TRACE &&
        _known(incident_monitor) && Float64(incident_monitor) <= 0
end

function backfill_trace_medians!(db::SQLite.DB)::Int
    rows = sqlrows(DBInterface.execute(db,
        "SELECT id, filename FROM exposures
          WHERE trace_median IS NULL AND trace_hash IS NOT NULL AND filename IS NOT NULL
          ORDER BY id"))
    filled = 0
    for r in rows
        id = Int(r.id)
        try
            dir = _resolve_analysis_dir(db, id)
            dir === nothing && continue
            path = resolve_trace_path(dir, integration_pattern_for_exposure(db, id), String(r.filename))
            path === nothing && continue
            m = trace_median(load_dat(path)[2])
            m === nothing && continue
            lock(_DB_WRITE_LOCK) do
                DBInterface.execute(db, "UPDATE exposures SET trace_median = ? WHERE id = ?", [m, id])
            end
            filled += 1
        catch err
            @warn "trace median backfill failed" exposure_id = id exception = err
        end
    end
    filled
end

function deselect_beam_off!(db::SQLite.DB, ids::AbstractVector{<:Integer})::Int
    isempty(ids) && return 0
    ph = join(fill("?", length(ids)), ",")
    dark = [Int(r.id) for r in sqlrows(DBInterface.execute(db, """
        SELECT id, trace_median, analysis_error, incident_monitor FROM exposures
         WHERE id IN ($ph) AND kind = 'file' AND deselected = 0
           AND (status IS NULL OR status != 'rejected')""", collect(Int, ids)))
            if is_beam_off(r.trace_median, r.analysis_error, r.incident_monitor)]
    isempty(dark) && return 0
    DBInterface.execute(db, "UPDATE exposures SET deselected = 1 WHERE id IN ($(join(fill("?", length(dark)), ",")))", dark)
    length(dark)
end

function deselect_beam_off_once!(db::SQLite.DB)::Int
    _migrated(db, MIGRATION_BEAM_OFF_DESELECTED) && return 0
    untouched = [Int(r.id) for r in sqlrows(DBInterface.execute(db, """
        SELECT e.id FROM exposures e
         WHERE e.kind = 'file' AND e.sample_id IS NOT NULL
           AND NOT EXISTS (SELECT 1 FROM exposures o
                            WHERE o.sample_id = e.sample_id AND o.kind = 'file' AND o.deselected = 1
                              AND (o.status IS NULL OR o.status != 'rejected'))"""))]
    n = lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            n = deselect_beam_off!(db, untouched)
            _record_migration!(db, MIGRATION_BEAM_OFF_DESELECTED)
            n
        end
    end
    n
end
