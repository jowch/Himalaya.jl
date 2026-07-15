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
