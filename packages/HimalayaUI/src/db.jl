using SQLite, DBInterface, Tables

const SCHEMA = """
CREATE TABLE IF NOT EXISTS users (
    id         INTEGER PRIMARY KEY,
    username   TEXT UNIQUE NOT NULL,
    first_name TEXT,
    last_name  TEXT
);

CREATE TABLE IF NOT EXISTS experiments (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    name            TEXT,
    path            TEXT NOT NULL,
    data_dir        TEXT NOT NULL,
    analysis_dir    TEXT NOT NULL,
    manifest_path   TEXT,
    config          TEXT,
    experiment_type TEXT,
    energy_kev      REAL,
    flight_path_m   REAL,
    created_at      DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS samples (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id INTEGER REFERENCES experiments(id),
    label         TEXT,
    name          TEXT,
    notes         TEXT
);

CREATE TABLE IF NOT EXISTS sample_tags (
    id        INTEGER PRIMARY KEY,
    sample_id INTEGER REFERENCES samples(id),
    key       TEXT NOT NULL,
    value     TEXT NOT NULL,
    source    TEXT DEFAULT 'manual'
);

CREATE TABLE IF NOT EXISTS exposures (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id  INTEGER REFERENCES samples(id),
    filename   TEXT,
    kind       TEXT DEFAULT 'file',
    selected   BOOLEAN DEFAULT FALSE,
    status     TEXT CHECK (status IN ('accepted', 'rejected')),
    image_path TEXT
);

CREATE TABLE IF NOT EXISTS exposure_sources (
    averaged_exposure_id INTEGER REFERENCES exposures(id),
    source_exposure_id   INTEGER REFERENCES exposures(id),
    role                 TEXT DEFAULT 'signal',
    PRIMARY KEY (averaged_exposure_id, source_exposure_id)
);

CREATE TABLE IF NOT EXISTS exposure_tags (
    id          INTEGER PRIMARY KEY,
    exposure_id INTEGER REFERENCES exposures(id),
    key         TEXT NOT NULL,
    value       TEXT NOT NULL,
    source      TEXT DEFAULT 'manual'
);

CREATE TABLE IF NOT EXISTS peaks (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    exposure_id INTEGER REFERENCES exposures(id),
    q           REAL NOT NULL,
    intensity   REAL,
    prominence  REAL,
    sharpness   REAL,
    source      TEXT DEFAULT 'auto',
    excluded    INTEGER DEFAULT 0
);

CREATE TABLE IF NOT EXISTS indices (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    exposure_id INTEGER REFERENCES exposures(id),
    phase       TEXT NOT NULL,
    basis       REAL NOT NULL,
    score       REAL,
    r_squared   REAL,
    lattice_d   REAL,
    status      TEXT DEFAULT 'candidate',
    kind        TEXT NOT NULL DEFAULT 'auto'
);

CREATE TABLE IF NOT EXISTS index_peaks (
    index_id       INTEGER REFERENCES indices(id),
    peak_id        INTEGER REFERENCES peaks(id),
    ratio_position INTEGER,
    residual       REAL,
    PRIMARY KEY (index_id, peak_id)
);

CREATE TABLE IF NOT EXISTS index_groups (
    id          INTEGER PRIMARY KEY,
    exposure_id INTEGER REFERENCES exposures(id),
    kind        TEXT NOT NULL DEFAULT 'auto',
    active      BOOLEAN DEFAULT FALSE,
    created_by  INTEGER REFERENCES users(id) ON DELETE SET NULL,
    created_at  DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS index_group_members (
    group_id  INTEGER REFERENCES index_groups(id),
    index_id  INTEGER REFERENCES indices(id),
    PRIMARY KEY (group_id, index_id)
);

CREATE TABLE IF NOT EXISTS sample_messages (
    id         INTEGER PRIMARY KEY,
    sample_id  INTEGER REFERENCES samples(id),
    author_id  INTEGER REFERENCES users(id) ON DELETE SET NULL,
    body       TEXT NOT NULL,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_sample_messages_sample
    ON sample_messages(sample_id, created_at);

CREATE TABLE IF NOT EXISTS user_actions (
    id          INTEGER PRIMARY KEY,
    user_id     INTEGER REFERENCES users(id) ON DELETE SET NULL,
    timestamp   DATETIME DEFAULT CURRENT_TIMESTAMP,
    action      TEXT,
    entity_type TEXT,
    entity_id   INTEGER,
    note        TEXT
);
"""

function create_schema!(db::SQLite.DB)
    for stmt in split(SCHEMA, ";")
        s = strip(stmt)
        isempty(s) && continue
        DBInterface.execute(db, s)
    end
end

function migrate_schema!(db::SQLite.DB)
    stmts = [
        "ALTER TABLE exposures ADD COLUMN status TEXT CHECK (status IN ('accepted', 'rejected'))",
        "ALTER TABLE exposures ADD COLUMN image_path TEXT",
        "ALTER TABLE experiments ADD COLUMN config TEXT",
        "ALTER TABLE experiments ADD COLUMN experiment_type TEXT",
        "ALTER TABLE experiments ADD COLUMN energy_kev REAL",
        "ALTER TABLE experiments ADD COLUMN flight_path_m REAL",
        "ALTER TABLE users ADD COLUMN first_name TEXT",
        "ALTER TABLE users ADD COLUMN last_name TEXT",
        "ALTER TABLE indices ADD COLUMN kind TEXT NOT NULL DEFAULT 'auto'",
    ]
    for stmt in stmts
        try
            DBInterface.execute(db, stmt)
        catch
            # column already exists — safe to ignore
        end
    end
    migrate_pk_to_autoincrement!(db)
end

"""
    migrate_pk_to_autoincrement!(db)

Rebuild the five entity tables that participate in chat `@`-mentions
(experiments, samples, exposures, peaks, indices) so their primary keys
are `INTEGER PRIMARY KEY AUTOINCREMENT`. SQLite's plain `INTEGER PRIMARY
KEY` is rowid-aliased and **reuses ids on deletion** — so a stale mention
of a deleted index can silently rebind to a new index that takes its id.
AUTOINCREMENT keeps a monotonically-increasing counter in `sqlite_sequence`
so freed ids are never reused.

No-op on fresh DBs (the schema already declares AUTOINCREMENT) and on
DBs that have already been migrated.
"""
function migrate_pk_to_autoincrement!(db::SQLite.DB)
    tables = ["experiments", "samples", "exposures", "peaks", "indices"]

    # Sentinel: skip iff every table in `tables` already has AUTOINCREMENT.
    # Checking just one would let a future addition to `tables` (a 6th
    # mention-target) get silently skipped on already-migrated DBs.
    # Skip the migration entirely if any table is missing — the migration
    # loop below assumes all five tables exist (it ALTER-renames each one),
    # and partial-schema fixtures aren't real production DBs anyway.
    needs_migration = false
    for t in tables
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='table' AND name=?", [t]))
        isempty(rows) && return  # one or more tables missing → leave alone
        if !occursin("AUTOINCREMENT", String(rows[1].sql))
            needs_migration = true
        end
    end
    needs_migration || return

    # FK enforcement must be disabled OUTSIDE a transaction (SQLite docs).
    DBInterface.execute(db, "PRAGMA foreign_keys = OFF")
    try
        # Rename old tables, create the fresh schema (the renamed tables no
        # longer exist so `CREATE TABLE IF NOT EXISTS` fires for them and is
        # a no-op for everyone else), copy rows, drop the renamed originals.
        SQLite.transaction(db) do
            for t in tables
                DBInterface.execute(db, "ALTER TABLE $t RENAME TO _migrate_old_$t")
            end
            create_schema!(db)
            for t in tables
                # Copy only the columns that exist in BOTH the renamed source
                # and the freshly-created destination — the old table may be
                # missing columns added by `migrate_schema!`'s ALTER TABLE pass
                # (those columns will land NULL, which the schema permits).
                new_cols = Set(String[String(c.name) for c in Tables.rowtable(
                    DBInterface.execute(db, "PRAGMA table_info($t)"))])
                old_cols = String[String(c.name) for c in Tables.rowtable(
                    DBInterface.execute(db, "PRAGMA table_info(_migrate_old_$t)"))]
                shared = join(String[c for c in old_cols if c in new_cols], ", ")
                DBInterface.execute(db,
                    "INSERT INTO $t ($shared) SELECT $shared FROM _migrate_old_$t")
                DBInterface.execute(db, "DROP TABLE _migrate_old_$t")
            end
        end
    finally
        DBInterface.execute(db, "PRAGMA foreign_keys = ON")
    end
end

function create_experiment!(db::SQLite.DB;
        name::Union{String,Nothing} = nothing,
        path::String,
        data_dir::String,
        analysis_dir::String,
        manifest_path::Union{String,Nothing} = nothing,
        config::Union{String,Nothing} = nothing,
        experiment_type::Union{String,Nothing} = nothing,
        energy_kev::Union{Float64,Nothing} = nothing,
        flight_path_m::Union{Float64,Nothing} = nothing)
    result = DBInterface.execute(db,
        """INSERT INTO experiments
             (name, path, data_dir, analysis_dir, manifest_path,
              config, experiment_type, energy_kev, flight_path_m)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        [name, path, data_dir, analysis_dir, manifest_path,
         config, experiment_type, energy_kev, flight_path_m])
    Int(DBInterface.lastrowid(result))
end

function create_sample!(db::SQLite.DB;
        experiment_id::Int,
        label::Union{String,Nothing} = nothing,
        name::Union{String,Nothing}  = nothing,
        notes::Union{String,Nothing} = nothing)
    result = DBInterface.execute(db,
        "INSERT INTO samples (experiment_id, label, name, notes) VALUES (?, ?, ?, ?)",
        [experiment_id, label, name, notes])
    Int(DBInterface.lastrowid(result))
end

function create_exposure!(db::SQLite.DB;
        sample_id::Int,
        filename::Union{String,Nothing}  = nothing,
        kind::String                     = "file",
        selected::Bool                   = false,
        status::Union{String,Nothing}    = nothing,
        image_path::Union{String,Nothing} = nothing)
    result = DBInterface.execute(db,
        "INSERT INTO exposures (sample_id, filename, kind, selected, status, image_path)
         VALUES (?, ?, ?, ?, ?, ?)",
        [sample_id, filename, kind, Int(selected), status, image_path])
    Int(DBInterface.lastrowid(result))
end

function get_experiment(db::SQLite.DB, id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM experiments WHERE id = ?", [id]))
    isempty(rows) && error("experiment $id not found")
    first(rows)
end

function get_samples(db::SQLite.DB, experiment_id::Int)
    Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM samples WHERE experiment_id = ? ORDER BY id", [experiment_id]))
end

function get_exposures(db::SQLite.DB, sample_id::Int)
    Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM exposures WHERE sample_id = ? ORDER BY id", [sample_id]))
end

"""
    default_db_path() -> String

Resolve the canonical Himalaya DB path. Reads `HIMALAYA_DB_PATH` from
the environment when set; otherwise falls back to `~/.himalaya/himalaya.db`
(creating the parent directory on first call).
"""
function default_db_path()::String
    haskey(ENV, "HIMALAYA_DB_PATH") && return ENV["HIMALAYA_DB_PATH"]
    dir = joinpath(homedir(), ".himalaya")
    isdir(dir) || mkpath(dir)
    joinpath(dir, "himalaya.db")
end

"""
    open_db(db_path = default_db_path()) -> SQLite.DB

Open the SQLite database at `db_path`, creating the file (and any missing
parent directories) if necessary. Applies schema migrations and enables
foreign-key enforcement on every connection.

Pass an explicit path for tests or alternate deployments; omit the
argument to use [`default_db_path`](@ref).
"""
function open_db(db_path::AbstractString = default_db_path())::SQLite.DB
    parent = dirname(db_path)
    !isempty(parent) && !isdir(parent) && mkpath(parent)
    db = SQLite.DB(db_path)
    create_schema!(db)
    migrate_schema!(db)
    DBInterface.execute(db, "PRAGMA foreign_keys = ON")

    # SQLite hardcodes O_CREAT mode 0644 in os_unix.c — process umask only
    # masks bits OUT, so umask 0002 can't promote 0644 to 0664. For
    # multi-user deploys (curators in a shared group writing the same DB),
    # we need group-write on the file. chmod is idempotent; if we don't
    # own the file (e.g. another user created it), this is a no-op error
    # we swallow rather than failing the whole open.
    if isfile(db_path)
        try
            chmod(db_path, 0o664)
        catch e
            # Swallow only the expected "not our file" / FS-permission errors;
            # let unexpected failures (InterruptException, oddities) propagate
            # so they don't get masked by this best-effort fix-up.
            e isa Base.IOError || e isa SystemError || rethrow()
        end
    end
    db
end
