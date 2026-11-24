"""
    experiments(db) -> Vector{<:NamedTuple}

All experiments (one per beamtime), ordered by id. Paths, file patterns and
geometry live on each setup: see [`setups`](@ref).
"""
experiments(db::SQLite.DB) = Tables.rowtable(DBInterface.execute(db,
    "SELECT id, name, path, description, created_at, last_scanned_at, ingest_status FROM experiments ORDER BY id"))

const _SETUP_COLS = """id, experiment_id, folder, display_name, position, included,
    root_dir, data_dir, analysis_dir, image_pattern, metadata_pattern, integration_pattern,
    energy_kev, flight_path_m, beam_center_x, beam_center_y, pixel_size_um, q_units,
    trace_normalisation, normalisation_base"""

"""
    setups(db; experiment=nothing) -> Vector{<:NamedTuple}

Setups (one per beamtime subfolder), optionally for one experiment id, in UI order.
Rows carry directories, file patterns, geometry (`energy_kev`, `flight_path_m`, beam
centre, `pixel_size_um`, `q_units`) and `trace_normalisation` (`"i2"`: traces divided
by I_2, `"none"`: raw counts) with `normalisation_base`. Join by `setup_id`:
setups of one beamtime can differ in flight path and energy.
"""
function setups(db::SQLite.DB; experiment::Union{Integer,Nothing}=nothing)
    if experiment === nothing
        Tables.rowtable(DBInterface.execute(db,
            "SELECT $_SETUP_COLS FROM setups ORDER BY experiment_id, position, id"))
    else
        Tables.rowtable(DBInterface.execute(db,
            "SELECT $_SETUP_COLS FROM setups WHERE experiment_id = ? ORDER BY position, id",
            [Int(experiment)]))
    end
end

const _SAMPLE_COLS = "id, experiment_id, setup_id, name, notes"

"""
    samples(db; experiment=nothing) -> Vector{<:NamedTuple}

Samples, optionally filtered to one experiment id.
"""
function samples(db::SQLite.DB; experiment::Union{Integer,Nothing}=nothing)
    if experiment === nothing
        Tables.rowtable(DBInterface.execute(db,
            "SELECT $_SAMPLE_COLS FROM samples ORDER BY id"))
    else
        Tables.rowtable(DBInterface.execute(db,
            "SELECT $_SAMPLE_COLS FROM samples WHERE experiment_id = ? ORDER BY id",
            [Int(experiment)]))
    end
end

const _EXPOSURE_COLS = """id, experiment_id, setup_id, sample_id, filename, kind,
    status, image_path, timestamp, exposure_time, frame_no,
    transmitted_counts, incident_monitor"""

"""
    exposures(db; sample=nothing) -> Vector{<:NamedTuple}

Exposures (HimalayaUI calls them frames), optionally filtered to one sample id.
`transmitted_counts` (beamstop `I_2`) and `incident_monitor` (`I_0`) are the beam
monitors read from the exposure's `.prp` file (`missing` when none was read).
"""
function exposures(db::SQLite.DB; sample::Union{Integer,Nothing}=nothing)
    if sample === nothing
        Tables.rowtable(DBInterface.execute(db,
            "SELECT $_EXPOSURE_COLS FROM exposures ORDER BY id"))
    else
        Tables.rowtable(DBInterface.execute(db,
            "SELECT $_EXPOSURE_COLS FROM exposures WHERE sample_id = ? ORDER BY id",
            [Int(sample)]))
    end
end

"""
    curated_peaks(db, exposure_id) -> Vector{<:NamedTuple}

The effective (curated) peak set for an exposure: `auto_peaks ∪ adds − excludes`.
Each row is tagged `source ∈ {"auto","manual"}` and `excluded ∈ {0,1}` — excluded
peaks are returned tagged, not removed from the result, so the truly-effective
set is `filter(p -> p.excluded == 0, peaks)`.
Excludes are matched to auto peaks by q within `MAX(1e-6, ABS(q)*0.001)`.

Mirrors `HimalayaUI.get_peaks_for_exposure`; the contract test guards drift.
"""
curated_peaks(db::SQLite.DB, exposure_id::Integer) = Tables.rowtable(DBInterface.execute(db, """
    SELECT a.id, a.exposure_id, a.q, a.intensity, a.prominence, a.sharpness,
           'auto' AS source,
           CASE WHEN c.q IS NOT NULL THEN 1 ELSE 0 END AS excluded
    FROM auto_peaks a
    LEFT JOIN peak_curations c
        ON c.exposure_id = a.exposure_id
       AND c.kind = 'exclude'
       AND ABS(c.q - a.q) <= MAX(1e-6, ABS(a.q) * 0.001)
    WHERE a.exposure_id = ?
    UNION ALL
    SELECT id, exposure_id, q, NULL AS intensity, NULL AS prominence, NULL AS sharpness,
           'manual' AS source, 0 AS excluded
    FROM peak_curations
    WHERE exposure_id = ? AND kind = 'add'
    ORDER BY q
""", [Int(exposure_id), Int(exposure_id)]))

const _INDEX_COLS = "id, exposure_id, phase, basis, score, r_squared, lattice_d, status, kind, inputs_hash"

"""
    index_candidates(db, exposure_id) -> Vector{<:NamedTuple}

All candidate index-choices for an exposure (auto + speculative), sorted by score.
Mirrors `HimalayaUI.get_indices_for_exposure`.
"""
index_candidates(db::SQLite.DB, exposure_id::Integer) = Tables.rowtable(DBInterface.execute(db,
    "SELECT $_INDEX_COLS FROM indices WHERE exposure_id = ? ORDER BY score DESC",
    [Int(exposure_id)]))

"""
    confirmed_indices(db, exposure_id) -> Vector{<:NamedTuple}

The exposure's durable *indexed assignment*: the index rows in `assignment_members`
when the exposure's `assignments.state` is `'indexed'` (the state defaults to
`'indexed'` when no assignments row exists). This is the same "confirmed index"
HimalayaUI reads. Empty when the state is `form_factor`/`null`. Sorted by score.

The phase call HimalayaUI shows belongs to the sample, not the exposure: use
[`phase_call`](@ref) for it.

A non-empty result is NOT by itself proof of an explicit human decision: HimalayaUI's
D-10 migration backfills the auto group's pick into `assignment_members` at
`state='indexed'` for already-analyzed exposures, so on upgraded databases this can
return an auto-seeded index. To tell an auto seed from a curator's choice, read the
assignment state.

Mirrors HimalayaUI's confirmed-index read: `assignment_members` JOIN `indices`, gated
by `assignments.state = 'indexed'` (the `COALESCE` keeps the no-row default faithful).
The legacy `index_groups`/`index_group_members` `kind='custom'` path this used to read
was retired by the D-10 redesign; those tables persist only as historical data.
"""
confirmed_indices(db::SQLite.DB, exposure_id::Integer) = Tables.rowtable(DBInterface.execute(db, """
    SELECT i.id, i.exposure_id, i.phase, i.basis, i.score, i.r_squared,
           i.lattice_d, i.status, i.kind, i.inputs_hash
    FROM assignment_members m
    JOIN indices i          ON i.id = m.index_id
    LEFT JOIN assignments a  ON a.exposure_id = m.exposure_id
    WHERE m.exposure_id = ? AND COALESCE(a.state, 'indexed') = 'indexed'
    ORDER BY i.score DESC NULLS LAST, i.id ASC
""", [Int(exposure_id)]))
