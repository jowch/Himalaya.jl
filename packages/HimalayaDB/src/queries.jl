"""
    experiments(db) -> Vector{<:NamedTuple}

All experiments, ordered by id.
"""
experiments(db::SQLite.DB) = Tables.rowtable(DBInterface.execute(db,
    "SELECT id, name, path, data_dir, analysis_dir, experiment_type, energy_kev, flight_path_m, created_at FROM experiments ORDER BY id"))

"""
    samples(db; experiment=nothing) -> Vector{<:NamedTuple}

Samples, optionally filtered to one experiment id.
"""
function samples(db::SQLite.DB; experiment::Union{Integer,Nothing}=nothing)
    if experiment === nothing
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, experiment_id, name, notes FROM samples ORDER BY id"))
    else
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, experiment_id, name, notes FROM samples WHERE experiment_id = ? ORDER BY id",
            [Int(experiment)]))
    end
end

"""
    exposures(db; sample=nothing) -> Vector{<:NamedTuple}

Exposures, optionally filtered to one sample id.
"""
function exposures(db::SQLite.DB; sample::Union{Integer,Nothing}=nothing)
    if sample === nothing
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, sample_id, filename, kind, selected, status, image_path FROM exposures ORDER BY id"))
    else
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, sample_id, filename, kind, selected, status, image_path FROM exposures WHERE sample_id = ? ORDER BY id",
            [Int(sample)]))
    end
end

"""
    curated_peaks(db, exposure_id) -> Vector{<:NamedTuple}

The effective (curated) peak set for an exposure: `auto_peaks ∪ adds − excludes`.
Each row is tagged `source ∈ {"auto","manual"}` and `excluded ∈ {0,1}`.
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

The human-confirmed indices: members of the exposure's `kind='custom'` index
group — what the curator committed to. Sorted by score. Empty when the curator
has never touched the exposure (no custom group exists yet).

Verified against HimalayaUI's write/read paths: human curation always lands in
the on-demand `kind='custom'` group (routes_analysis.jl `ensure_custom_group!`),
and `kind='custom' ⟹ active=1`, so filtering on `active=1` alone would wrongly
include the pre-curation auto group. No single HimalayaUI getter returns this
shape, so the members join is composed here.
"""
confirmed_indices(db::SQLite.DB, exposure_id::Integer) = Tables.rowtable(DBInterface.execute(db, """
    SELECT i.id, i.exposure_id, i.phase, i.basis, i.score, i.r_squared,
           i.lattice_d, i.status, i.kind, i.inputs_hash
    FROM indices i
    JOIN index_group_members m ON m.index_id = i.id
    JOIN index_groups g        ON g.id = m.group_id
    WHERE g.exposure_id = ? AND g.kind = 'custom'
    ORDER BY i.score DESC
""", [Int(exposure_id)]))
