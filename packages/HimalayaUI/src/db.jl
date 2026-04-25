using SQLite, DBInterface, Tables

const SCHEMA = """
CREATE TABLE IF NOT EXISTS users (
    id       INTEGER PRIMARY KEY,
    username TEXT UNIQUE NOT NULL
);

CREATE TABLE IF NOT EXISTS experiments (
    id             INTEGER PRIMARY KEY,
    name           TEXT,
    path           TEXT NOT NULL,
    data_dir       TEXT NOT NULL,
    analysis_dir   TEXT NOT NULL,
    manifest_path  TEXT,
    created_at     DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS samples (
    id            INTEGER PRIMARY KEY,
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
    id        INTEGER PRIMARY KEY,
    sample_id INTEGER REFERENCES samples(id),
    filename  TEXT,
    kind      TEXT DEFAULT 'file',
    selected  BOOLEAN DEFAULT FALSE
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
    id          INTEGER PRIMARY KEY,
    exposure_id INTEGER REFERENCES exposures(id),
    q           REAL NOT NULL,
    intensity   REAL,
    prominence  REAL,
    sharpness   REAL,
    source      TEXT DEFAULT 'auto',
    excluded    INTEGER DEFAULT 0
);

CREATE TABLE IF NOT EXISTS indices (
    id          INTEGER PRIMARY KEY,
    exposure_id INTEGER REFERENCES exposures(id),
    phase       TEXT NOT NULL,
    basis       REAL NOT NULL,
    score       REAL,
    r_squared   REAL,
    lattice_d   REAL,
    status      TEXT DEFAULT 'candidate'
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
    created_by  INTEGER REFERENCES users(id),
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
    user_id     INTEGER REFERENCES users(id),
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

function create_experiment!(db::SQLite.DB;
        name::Union{String,Nothing} = nothing,
        path::String,
        data_dir::String,
        analysis_dir::String,
        manifest_path::Union{String,Nothing} = nothing)
    result = DBInterface.execute(db,
        "INSERT INTO experiments (name, path, data_dir, analysis_dir, manifest_path)
         VALUES (?, ?, ?, ?, ?)",
        [name, path, data_dir, analysis_dir, manifest_path])
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
        filename::Union{String,Nothing} = nothing,
        kind::String                    = "file",
        selected::Bool                  = false)
    result = DBInterface.execute(db,
        "INSERT INTO exposures (sample_id, filename, kind, selected) VALUES (?, ?, ?, ?)",
        [sample_id, filename, kind, Int(selected)])
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

function open_db(experiment_path::String)::SQLite.DB
    db_path = joinpath(experiment_path, "himalaya.db")
    db = SQLite.DB(db_path)
    create_schema!(db)
    db
end
