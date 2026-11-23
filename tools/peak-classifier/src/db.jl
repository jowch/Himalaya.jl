using SQLite, DBInterface

"""
    open_readonly(db_path) -> SQLite.DB

Opens with the URI `mode=ro`, so a missing file is an error rather than a new
empty database, and closing the connection never checkpoints a `-wal` file.
`PRAGMA query_only` is a second guard. A WAL database may still get empty
`-shm` and `-wal` side files from a read-only open. This is
`HimalayaUI.open_readonly_db`.
"""
open_readonly(db_path::AbstractString) = HimalayaUI.open_readonly_db(db_path)

"""
    load_labeled_exposures(db_path; group_by = :experiment, include_form_factor = true,
                           io = stderr) -> Vector{LabeledExposure}

Opens the database read-only (`open_readonly`, no migration) and reads it with
`HimalayaUI.labeled_exposures`.
"""
function load_labeled_exposures(db_path::AbstractString; kwargs...)
    db = open_readonly(db_path)
    try
        HimalayaUI.labeled_exposures(db; kwargs...)
    finally
        close(db)
    end
end
