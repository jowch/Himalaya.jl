using TOML

const DEFAULT_INTEGRATION_PATTERN = "{name}.dat"

_name_pattern(v) = (v isa AbstractString && occursin("{name}", v)) ? String(v) : nothing

function _legacy_blob_integration(blob)
    blob isa AbstractString || return nothing
    files = try
        get(TOML.parse(blob), "files", Dict())
    catch
        return nothing
    end
    files isa AbstractDict ? _name_pattern(get(files, "integration", nothing)) : nothing
end

"""
    integration_pattern_for(db, experiment_id) -> String

The integration-file pattern for an experiment: its `integration_pattern`
column, else `[files].integration` from a legacy `experiments.config` blob,
else `$(DEFAULT_INTEGRATION_PATTERN)`. Values without `{name}` are ignored.
"""
function integration_pattern_for(db::SQLite.DB, experiment_id::Integer)::String
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT config, integration_pattern FROM experiments WHERE id = ?",
        [Int(experiment_id)]))
    isempty(rows) && error("Experiment $experiment_id not found")
    r = rows[1]
    something(_name_pattern(r.integration_pattern),
              _legacy_blob_integration(r.config),
              DEFAULT_INTEGRATION_PATTERN)
end
