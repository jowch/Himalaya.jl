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
    return db
end
