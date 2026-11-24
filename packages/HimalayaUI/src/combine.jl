using SHA

const COMBINE_METHOD = "time_weighted_sum"
const MIN_COMBINED_SOURCES = 2
const COMBINED_KIND = "combined"
const Q_GRID_RTOL = 1e-6

# physics: the average is the unweighted mean ΣIᵢ / N. SASTool usually divides each
# physics: trace by its own I_2 (the setup's trace_normalisation), so each image is
# physics: already intensity per unit of transmitted beam and the mean is unbiased,
# physics: though not the lowest-noise estimate: that weights each image by its I_2,
# physics: Σ(Iᵢ·I₂ᵢ) / ΣI₂ᵢ. On a raw-count setup ΣIᵢ / N is the pooled sum put back on
# physics: one image's scale, which stays unbiased when the beam decays between
# physics: images. Errors add in quadrature either way: σ = √(Σσᵢ²) / N.
function combine_traces(traces::AbstractVector)
    isempty(traces) && throw(ArgumentError("nothing to combine"))
    q = traces[1][1]
    for t in traces
        length(t[1]) == length(q) || throw(ArgumentError("images have different q grids"))
        all(isapprox.(t[1], q; rtol = Q_GRID_RTOL)) || throw(ArgumentError("images have different q grids"))
    end
    n = length(traces)
    I = sum(t[2] for t in traces) ./ n
    σ = sqrt.(sum(t[3] .^ 2 for t in traces)) ./ n
    copy(q), I, σ
end

is_combined_kind(kind) = kind isa AbstractString && kind == COMBINED_KIND

function selected_image_ids(db::SQLite.DB, sample_id::Int)::Vector{Int}
    Int[Int(r.id) for r in sqlrows(DBInterface.execute(db, """
        SELECT id FROM exposures
         WHERE sample_id = ? AND kind = 'file'
           AND (status IS NULL OR status != 'rejected')
           AND deselected = 0
         ORDER BY id""", [sample_id]))]
end

function _all_traced(db::SQLite.DB, ids::Vector{Int})::Bool
    isempty(ids) && return false
    ph = join(fill("?", length(ids)), ",")
    n = first(sqlrows(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM exposures WHERE id IN ($ph) AND trace_hash IS NOT NULL", ids))).n
    Int(n) == length(ids)
end

function combined_exposure_ids(db::SQLite.DB, sample_id::Int)::Vector{Int}
    Int[Int(r.id) for r in sqlrows(DBInterface.execute(db,
        "SELECT id FROM exposures WHERE sample_id = ? AND kind = 'combined' ORDER BY id", [sample_id]))]
end

function combined_source_ids(db::SQLite.DB, combined_id::Int)::Vector{Int}
    Int[Int(r.source_exposure_id) for r in sqlrows(DBInterface.execute(db,
        "SELECT source_exposure_id FROM exposure_sources WHERE averaged_exposure_id = ? ORDER BY source_exposure_id",
        [combined_id]))]
end

function combined_for_set(db::SQLite.DB, sample_id::Int, ids::Vector{Int})::Union{Int,Nothing}
    want = sort(ids)
    for cid in combined_exposure_ids(db, sample_id)
        combined_source_ids(db, cid) == want && return cid
    end
    nothing
end

function analyzed_exposure_id(db::SQLite.DB, sample_id::Int)::Union{Int,Nothing}
    sel = selected_image_ids(db, sample_id)
    isempty(sel) && return nothing
    length(sel) == 1 && return sel[1]
    combined_for_set(db, sample_id, sel)
end

function ensure_combined!(db::SQLite.DB, sample_id::Int, ids::Vector{Int})::Union{Int,Nothing}
    length(ids) < MIN_COMBINED_SOURCES && return nothing
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            existing = combined_for_set(db, sample_id, ids)
            existing === nothing || return existing
            _create_combined!(db, sample_id, sort(ids))
        end
    end
end

function _create_combined!(db::SQLite.DB, sample_id::Int, sources::Vector{Int})::Int
    sample = first(sqlrows(DBInterface.execute(db,
        "SELECT experiment_id, setup_id, load_id FROM samples WHERE id = ?", [sample_id])))
    exp_id = Int(sample.experiment_id)
    cid = create_exposure!(db;
        experiment_id = exp_id,
        setup_id      = ismissing(sample.setup_id) ? default_setup_id(db, exp_id) : Int(sample.setup_id),
        sample_id     = sample_id,
        kind          = COMBINED_KIND,
        load_id       = ismissing(sample.load_id) ? nothing : Int(sample.load_id),
        monitors_read = true)
    for sid in sources
        DBInterface.execute(db,
            "INSERT INTO exposure_sources (averaged_exposure_id, source_exposure_id, role) VALUES (?, ?, 'signal')",
            [cid, sid])
    end
    _refresh_combined_metadata!(db, cid, sources)
    cid
end

function _refresh_combined_metadata!(db::SQLite.DB, cid::Int, sources::Vector{Int})
    if isempty(sources)
        DBInterface.execute(db,
            "UPDATE exposures SET exposure_time = NULL, timestamp = NULL WHERE id = ?", [cid])
        return nothing
    end
    ph = join(fill("?", length(sources)), ",")
    agg = first(sqlrows(DBInterface.execute(db, """
        SELECT SUM(exposure_time) AS total_time, COUNT(exposure_time) AS timed, MAX(timestamp) AS latest
          FROM exposures WHERE id IN ($ph)""", sources)))
    total = Int(agg.timed) == length(sources) ? Float64(agg.total_time) : nothing
    DBInterface.execute(db,
        "UPDATE exposures SET exposure_time = ?, timestamp = ? WHERE id = ?",
        [total, coalesce(agg.latest, nothing), cid])
    nothing
end

function combined_trace_hash(db::SQLite.DB, cid::Int)::Union{String,Nothing}
    rows = sqlrows(DBInterface.execute(db, """
        SELECT e.id, e.trace_hash FROM exposure_sources s JOIN exposures e ON e.id = s.source_exposure_id
         WHERE s.averaged_exposure_id = ? ORDER BY e.id""", [cid]))
    length(rows) < MIN_COMBINED_SOURCES && return nothing
    any(r -> ismissing(r.trace_hash), rows) && return nothing
    key = join([string(Int(r.id), ":", String(r.trace_hash)) for r in rows], ",")
    bytes2hex(SHA.sha256(string(COMBINE_METHOD, "|", key)))
end

function load_combined_trace(db::SQLite.DB, cid::Int)
    sources = combined_source_ids(db, cid)
    length(sources) < MIN_COMBINED_SOURCES &&
        throw(ArgumentError("fewer than $MIN_COMBINED_SOURCES images to combine"))
    traces = map(sources) do sid
        path = _source_trace_path(db, sid)
        path === nothing && throw(ArgumentError("no trace for image $sid"))
        load_dat(path)
    end
    combine_traces(traces)
end

function _source_trace_path(db::SQLite.DB, exposure_id::Int)
    rows = sqlrows(DBInterface.execute(db, """
        SELECT e.filename, st.analysis_dir FROM exposures e JOIN setups st ON st.id = e.setup_id
         WHERE e.id = ?""", [exposure_id]))
    isempty(rows) && return nothing
    r = rows[1]
    (ismissing(r.filename) || ismissing(r.analysis_dir)) && return nothing
    resolve_trace_path(String(r.analysis_dir), integration_pattern_for_exposure(db, exposure_id), String(r.filename))
end

function _combined_normalisation(db::SQLite.DB, cid::Int)::String
    rows = sqlrows(DBInterface.execute(db,
        "SELECT st.trace_normalisation FROM exposures e JOIN setups st ON st.id = e.setup_id WHERE e.id = ?", [cid]))
    isempty(rows) ? NORMALISATION_UNKNOWN : coalesce(rows[1].trace_normalisation, NORMALISATION_UNKNOWN)
end

function combined_monitors(db::SQLite.DB, cid::Int)
    sources = combined_source_ids(db, cid)
    normalised_by = _combined_normalisation(db, cid)
    none = Dict{Symbol,Any}(:transmitted_counts => nothing, :exposure_time => nothing, :counting_time => nothing,
                            :normalised_by => normalised_by, :beam_off => false)
    isempty(sources) && return none
    ph = join(fill("?", length(sources)), ",")
    rows = sqlrows(DBInterface.execute(db,
        "SELECT transmitted_counts, exposure_time, counting_time, trace_median, analysis_error, incident_monitor
           FROM exposures WHERE id IN ($ph)", sources))
    beam_off = any(r -> is_beam_off(r.trace_median, r.analysis_error, r.incident_monitor), rows)
    usable(v) = !ismissing(v) && isfinite(Float64(v)) && Float64(v) > 0
    all(r -> usable(r.transmitted_counts) && usable(r.exposure_time) && usable(r.counting_time), rows) ||
        return merge(none, Dict{Symbol,Any}(:beam_off => beam_off))
    flux = sum(Float64(r.transmitted_counts) * Float64(r.exposure_time) / Float64(r.counting_time) for r in rows) / length(rows)
    Dict{Symbol,Any}(:transmitted_counts => flux, :exposure_time => 1.0, :counting_time => 1.0,
                     :normalised_by => normalised_by, :beam_off => beam_off)
end

function combined_due(db::SQLite.DB, cid::Int)::Bool
    row = first(sqlrows(DBInterface.execute(db,
        "SELECT trace_hash, analysis_version, analysis_error, analysis_trace_stat FROM exposures WHERE id = ?", [cid])))
    current = combined_trace_hash(db, cid)
    current === nothing && return false
    ismissing(row.analysis_error) || return coalesce(row.analysis_version, "") != ANALYSIS_VERSION ||
        coalesce(row.analysis_trace_stat, "") != current
    ismissing(row.trace_hash) && return true
    String(row.trace_hash) != current && return true
    coalesce(row.analysis_version, "") != ANALYSIS_VERSION && !peak_reviewed(db, cid)
end

function refresh_selection!(db::SQLite.DB, sample_id::Int)::Union{Int,Nothing}
    is_autosampler_sample(db, sample_id) && return analyzed_exposure_id(db, sample_id)
    sel = selected_image_ids(db, sample_id)
    length(sel) < MIN_COMBINED_SOURCES && return isempty(sel) ? nothing : sel[1]
    _all_traced(db, sel) || return nothing
    cid = ensure_combined!(db, sample_id, sel)
    combined_due(db, cid) || return cid
    try
        analyze_exposure!(db, cid)
    catch e
        @warn "analyze combined image failed" sample_id exception = e
        lock(_DB_WRITE_LOCK) do
            DBInterface.execute(db,
                "UPDATE exposures SET analysis_error = ?, analysis_version = ?, analysis_trace_stat = ? WHERE id = ?",
                [first(sprint(showerror, e), 500), ANALYSIS_VERSION, combined_trace_hash(db, cid), cid])
        end
    end
    cid
end

function backfill_averages!(db::SQLite.DB)
    samples = [Int(r.sample_id) for r in sqlrows(DBInterface.execute(db, """
        SELECT sample_id FROM exposures
         WHERE kind = 'file' AND sample_id IS NOT NULL
           AND (status IS NULL OR status != 'rejected') AND deselected = 0
         GROUP BY sample_id HAVING COUNT(*) >= ?""", [MIN_COMBINED_SOURCES]))]
    foreach(sid -> refresh_selection!(db, sid), samples)
    length(samples)
end
