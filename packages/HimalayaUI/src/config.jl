const DEFAULT_INTEGRATION_PATTERN = "{name}.dat"

_name_pattern(v) = (v isa AbstractString && occursin("{name}", v)) ? String(v) : nothing

integration_pattern_for(setup) =
    something(_name_pattern(coalesce(setup.integration_pattern, nothing)), DEFAULT_INTEGRATION_PATTERN)

function integration_pattern_for_exposure(db::SQLite.DB, exposure_id::Integer)::String
    setup = exposure_setup(db, exposure_id)
    setup === nothing && error("exposure $exposure_id has no setup")
    integration_pattern_for(setup)
end
