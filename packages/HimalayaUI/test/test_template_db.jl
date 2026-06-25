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
    # Checkpoint the WAL into the main file so VACUUM INTO copies a quiescent single
    # file with no -wal sidecar. wal_checkpoint(TRUNCATE) returns SQLITE_LOCKED
    # ("database table is locked") while ANY statement holds a read mark on the
    # connection. Two sources of such a mark here: (a) a cached prepared statement
    # left by open_db's migration chain, cleared by finalize_statements!; (b) under
    # the multithreaded suite (nthreads≫1) a background GC finalizer freeing a prior
    # test's SQLite.Stmt can transiently grab the lock in the window before the
    # checkpoint. So: flush pending finalizers, finalize, then bounded-retry. The
    # result iterator is drained (rowtable) so no statement is in-flight for the
    # subsequent VACUUM INTO. ponytail: 5 attempts is the ceiling — one retry has
    # always sufficed in practice; raise it only if a checkpoint ever stays locked.
    local checkpointed = false
    for _ in 1:5
        GC.gc()                              # run pending SQLite.Stmt finalizers first
        SQLite.finalize_statements!(db)
        try
            Tables.rowtable(DBInterface.execute(db, "PRAGMA wal_checkpoint(TRUNCATE)"))
            checkpointed = true
            break
        catch
            # transient lock — loop to flush finalizers and retry
        end
    end
    checkpointed || error("prepared_template_path: WAL checkpoint stayed locked after 5 retries")
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
