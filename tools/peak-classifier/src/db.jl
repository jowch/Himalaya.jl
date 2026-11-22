using SQLite, DBInterface

# Percent-encode the characters that are special in a SQLite URI path.
_uri_path(p) = replace(abspath(p), "%" => "%25", "?" => "%3f", "#" => "%23")

"""
    open_readonly(db_path) -> SQLite.DB

Opens with the URI `mode=ro`, so a missing file is an error rather than a new
empty database, and closing the connection never checkpoints a `-wal` file.
`PRAGMA query_only` is a second guard. A WAL database may still get empty
`-shm` and `-wal` side files from a read-only open.
"""
function open_readonly(db_path::AbstractString)
    db = SQLite.DB("file:" * _uri_path(db_path) * "?mode=ro")
    try
        DBInterface.execute(db, "PRAGMA query_only = ON")
    catch
        close(db)
        rethrow()
    end
    db
end

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
