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
