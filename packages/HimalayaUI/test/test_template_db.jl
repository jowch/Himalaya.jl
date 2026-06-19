# M2: build a fully-migrated, zero-row template DB ONCE per process, then clone
# it per test instead of re-running create_schema!+migrate_schema! each time.
# Scoped to the FK-ON open_db family only. VACUUM INTO produces a compacted,
# WAL-free single file; the clone re-applies the connection-scoped setup that
# open_db does (foreign_keys=ON, WAL, finalize_statements!) — none of which
# survive a file copy.
using SQLite, DBInterface, Tables
using HimalayaUI

const _TEMPLATE_REF = Ref{Union{String,Nothing}}(nothing)

"Build (once) and return the path to the migrated, checkpointed template DB."
function prepared_template_path()
    _TEMPLATE_REF[] === nothing || return _TEMPLATE_REF[]
    tmpdir = mktempdir(; cleanup = false)         # lives for the process
    src = joinpath(tmpdir, "source.db")
    db = HimalayaUI.open_db(src)                  # full create_schema! + migrate_schema!
    # Checkpoint the WAL into the main file and close, so VACUUM INTO copies a
    # quiescent single file with no -wal sidecar. Drain the result iterator
    # (same pattern as open_db's journal_mode=WAL) so no statement is in-flight
    # when VACUUM INTO runs — SQLite rejects VACUUM while any statement is open.
    Tables.rowtable(DBInterface.execute(db, "PRAGMA wal_checkpoint(TRUNCATE)"))
    SQLite.finalize_statements!(db)
    tmpl = joinpath(tmpdir, "template.db")
    DBInterface.execute(db, "VACUUM INTO '$(tmpl)'")   # target must NOT pre-exist
    SQLite.close(db)
    _TEMPLATE_REF[] = tmpl
    return tmpl
end

"Clone the template into `dir` and return a connection equivalent to open_db's."
function open_prepared_clone(dir::AbstractString)
    tmpl = prepared_template_path()
    dest = joinpath(dir, "h.db")
    cp(tmpl, dest; force = false)                 # copy ONLY the main file
    db = SQLite.DB(dest)
    SQLite.finalize_statements!(db)
    DBInterface.execute(db, "PRAGMA foreign_keys = ON")
    Tables.rowtable(DBInterface.execute(db, "PRAGMA journal_mode = WAL"))
    SQLite.finalize_statements!(db)
    return db
end
