const DEFAULT_INTEGRATION_PATTERN = "{name}.dat"

_name_pattern(v) = (v isa AbstractString && occursin("{name}", v)) ? String(v) : nothing

"""
    integration_pattern_for(db, experiment_id) -> String

The integration-file pattern for an experiment: its `integration_pattern`
column, else `$(DEFAULT_INTEGRATION_PATTERN)`. Values without `{name}` are ignored.
"""
function integration_pattern_for(db::SQLite.DB, experiment_id::Integer)::String
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT integration_pattern FROM experiments WHERE id = ?",
        [Int(experiment_id)]))
    isempty(rows) && error("Experiment $experiment_id not found")
    something(_name_pattern(rows[1].integration_pattern), DEFAULT_INTEGRATION_PATTERN)
end
