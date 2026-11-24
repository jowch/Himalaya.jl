using DBInterface, SQLite, Dates

const _TIMELINE_FORMAT = dateformat"yyyy-mm-ddTHH:MM:SS"

_timeline_time(ts) = ismissing(ts) ? nothing : tryparse(DateTime, String(ts), _TIMELINE_FORMAT)

_seconds_between(a::DateTime, b::DateTime) = Int(fld(Dates.value(b - a), 1000))

function _timed_images(db::SQLite.DB, exp_id::Integer)
    rows = sqlrows(DBInterface.execute(db, """
        SELECT e.sample_id, e.timestamp, e.exposure_time, e.trace_median,
               e.analysis_error, e.incident_monitor
          FROM exposures e
          JOIN samples s ON s.id = e.sample_id
          LEFT JOIN setups st ON st.id = s.setup_id
         WHERE e.experiment_id = ? AND e.kind = 'file' AND e.timestamp IS NOT NULL
           AND (st.id IS NULL OR st.included = 1)""", [exp_id]))
    images = NamedTuple{(:sample_id, :ends, :starts, :beam_off), Tuple{Int, DateTime, DateTime, Bool}}[]
    for r in rows
        t = _timeline_time(r.timestamp)
        t === nothing && continue
        exposure = ismissing(r.exposure_time) ? 0.0 : Float64(r.exposure_time)
        push!(images, (sample_id = Int(r.sample_id), ends = t,
            starts = t - Millisecond(round(Int, 1000 * exposure)),
            beam_off = is_beam_off(r.trace_median, r.analysis_error, r.incident_monitor)))
    end
    images
end

function reviewed_samples(db::SQLite.DB, sample_ids::Vector{Int})
    states = sample_call_states(db, sample_ids)
    members = sample_call_members(db, sample_ids)
    Set(sid for sid in sample_ids
        if export_call_state(get(states, sid, "indexed"), get(members, sid, [])) != "uncalled")
end

function shown_sample_ids(db::SQLite.DB, exp_id::Integer)
    Int[Int(r.id) for r in sqlrows(DBInterface.execute(db, """
        SELECT s.id FROM samples s LEFT JOIN setups st ON st.id = s.setup_id
         WHERE s.experiment_id = ? AND (st.id IS NULL OR st.included = 1)""", [exp_id]))]
end

function review_progress(db::SQLite.DB, exp_id::Integer)
    shown = shown_sample_ids(db, exp_id)
    (reviewed = length(reviewed_samples(db, shown)), total = length(shown))
end

function collection_timeline(db::SQLite.DB, exp_id::Integer)
    images = _timed_images(db, exp_id)
    isempty(images) && return nothing
    origin = minimum(im.starts for im in images)
    spans = Dict{Int, Tuple{DateTime, DateTime}}()
    for im in images
        a, b = get(spans, im.sample_id, (im.starts, im.ends))
        spans[im.sample_id] = (min(a, im.starts), max(b, im.ends))
    end
    reviewed = reviewed_samples(db, collect(keys(spans)))
    ordered = sort(collect(spans); by = kv -> (kv[2][1], kv[1]))
    (start = Dates.format(origin, _TIMELINE_FORMAT),
     samples = [[_seconds_between(origin, a), _seconds_between(origin, b), Int(sid in reviewed)]
                for (sid, (a, b)) in ordered],
     beam_off = sort([_seconds_between(origin, im.ends) for im in images if im.beam_off]))
end
