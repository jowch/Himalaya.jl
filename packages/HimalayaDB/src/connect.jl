"""
    default_himalaya_db_path() -> String

Resolve the DB path from `HIMALAYA_DB_PATH`, else `~/.himalaya/himalaya.db`.
Never creates directories (unlike HimalayaUI's `default_db_path`).
"""
default_himalaya_db_path() =
    get(ENV, "HIMALAYA_DB_PATH", joinpath(homedir(), ".himalaya", "himalaya.db"))

"""
    connect(path = default_himalaya_db_path()) -> SQLite.DB

Open the HimalayaUI database read-only. Errors if the file is missing.
Never runs migrations and never chmods the file.
"""
function connect(path::AbstractString = default_himalaya_db_path())
    isfile(path) || throw(ArgumentError(
        "HimalayaDB.connect: no database at $path (set HIMALAYA_DB_PATH?)"))
    db = SQLite.DB(path)
    # ponytail: SQLite.jl 1.8 opens RW at the OS layer (no readonly flag on DB());
    # query_only=ON makes every write/DDL fail loudly, which is all a reader needs.
    # Upgrade path if a truly read-only filesystem is ever required: drop to
    # SQLite.C.sqlite3_open_v2 with SQLITE_OPEN_READONLY.
    DBInterface.execute(db, "PRAGMA query_only = ON")
    _warn_if_unmigrated(db, path)
    return db
end

# #304 renumbered every Hexagonal `ratio_position` past 5 when √11 left
# phaseratios. `reconstruct_index` indexes straight into that series, so on a DB
# predating the migration — an old backup, or a deploy not yet restarted — every
# Hexagonal position 6+ reads one radicand high and `fit` returns a wrong lattice
# with no error. This package opens `query_only` and never migrates, so it cannot
# heal that itself; warning at connect is the only honest signal. Warn, don't
# throw: reading an old backup is legitimate, and every non-Hexagonal phase is
# unaffected.
function _warn_if_unmigrated(db::SQLite.DB, path::AbstractString)
    applied = try
        !isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM schema_migrations WHERE name = 'hex_sqrt11_removal_v1'")))
    catch
        # No schema_migrations table at all — far older than #304, same verdict.
        false
    end
    applied && return nothing
    hex = try
        first(Tables.rowtable(DBInterface.execute(db,
            "SELECT count(*) AS n FROM indices WHERE phase IN ('Hexagonal', 'Himalaya.Hexagonal')"))).n
    catch
        0
    end
    hex == 0 && return nothing   # nothing to mis-read
    @warn """
          $path predates the #304 Hexagonal renumbering (migration sentinel \
          'hex_sqrt11_removal_v1' absent). Hexagonal indices claiming ratio \
          position 6 or higher will reconstruct one reflection too high, so \
          `reconstruct_index` and `Himalaya.fit` return wrong lattice constants \
          for them. Open the database once with HimalayaUI to migrate it, then \
          re-run `analyze`. Other phases are unaffected.
          """ hexagonal_indices = hex
    nothing
end
