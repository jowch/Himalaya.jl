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

-- index_peaks: peak_id references auto_peaks OR peak_curations (peak_kind disambiguates).
-- Existing rows are all 'auto' (manual-peak refs get repointed during migration).
CREATE TABLE IF NOT EXISTS index_peaks (
    index_id       INTEGER REFERENCES indices(id),
    peak_id        INTEGER NOT NULL,
    peak_kind      TEXT NOT NULL DEFAULT 'auto'
                   CHECK (peak_kind IN ('auto', 'curation')),
    ratio_position INTEGER,
    residual       REAL,
    PRIMARY KEY (index_id, peak_id, peak_kind)
);

CREATE TABLE IF NOT EXISTS auto_peaks (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    exposure_id     INTEGER REFERENCES exposures(id),
    q               REAL NOT NULL,
    intensity       REAL,
    prominence      REAL,
    sharpness       REAL,
    findpeaks_index INTEGER
);

CREATE TABLE IF NOT EXISTS peak_curations (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    exposure_id INTEGER REFERENCES exposures(id),
    kind        TEXT NOT NULL CHECK (kind IN ('exclude', 'add')),
    q           REAL NOT NULL,
    created_by  INTEGER REFERENCES users(id) ON DELETE SET NULL,
    created_at  DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_auto_peaks_exposure
    ON auto_peaks(exposure_id);
CREATE INDEX IF NOT EXISTS idx_peak_curations_exposure
    ON peak_curations(exposure_id);

CREATE TABLE IF NOT EXISTS index_groups (
    id          INTEGER PRIMARY KEY,
    exposure_id INTEGER REFERENCES exposures(id),
    kind        TEXT NOT NULL DEFAULT 'auto',
    active      BOOLEAN DEFAULT FALSE,
    created_by  INTEGER REFERENCES users(id) ON DELETE SET NULL,
    created_at  DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE UNIQUE INDEX IF NOT EXISTS idx_one_custom_group_per_exposure
    ON index_groups(exposure_id) WHERE kind = 'custom';

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

"""
    preflight_index_groups_uniqueness!(db)

If a multiplayer-era duplicate-custom-group row exists in a pre-R0.1 DB,
fail loudly with a useful message rather than letting `CREATE UNIQUE INDEX`
produce SQLite's terse "UNIQUE constraint failed" error — operators wouldn't
know that "merge the duplicate custom groups" is the right next step.
No-op on truly-fresh DBs (the `index_groups` table doesn't exist yet).
"""
function preflight_index_groups_uniqueness!(db::SQLite.DB)
    has_table = !isempty(Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='index_groups'")))
    has_table || return
    dups = Tables.rowtable(DBInterface.execute(db, """
        SELECT exposure_id, COUNT(*) AS n FROM index_groups
        WHERE kind = 'custom' GROUP BY exposure_id HAVING n > 1
    """))
    if !isempty(dups)
        error("DB has duplicate 'custom' index_groups for exposures " *
              join([string(d.exposure_id) for d in dups], ", ") *
              " — manual merge required before idx_one_custom_group_per_exposure can be enforced")
    end
end

function create_schema!(db::SQLite.DB)
    for stmt in split(SCHEMA, ";")
        s = strip(stmt)
        isempty(s) && continue
        # Defensive: skip fragments that are purely SQL comments. Dead on the
        # current SCHEMA (no `;` inside any `--` comment), but a future SCHEMA
        # edit that puts `;` inside a comment would otherwise leave a
        # pure-comment fragment that DBInterface.execute rejects.
        all(l -> isempty(strip(l)) || startswith(strip(l), "--"), split(s, "\n")) && continue
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
    migrate_r2_widen_index_peaks_pk!(db)  # rebuild with widened PK first
    migrate_r2_split_peaks!(db)            # then repoint manual-peak refs
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

    # Sentinel: skip iff every table in `tables` that EXISTS already has
    # AUTOINCREMENT. "peaks" may no longer exist in R2.2+ DBs (removed from
    # SCHEMA); skip it if absent so migrate_pk_to_autoincrement! stays
    # idempotent on fresh DBs.
    needs_migration = false
    for t in tables
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='table' AND name=?", [t]))
        isempty(rows) && continue  # table absent (e.g. peaks on R2.2+ DBs) — skip
        if !occursin("AUTOINCREMENT", String(rows[1].sql))
            needs_migration = true
        end
    end
    # If ALL tables are either absent or already have AUTOINCREMENT, no-op.
    needs_migration || return

    # Only migrate tables that actually exist AND are still in SCHEMA (can be
    # recreated by create_schema!). `peaks` was removed from SCHEMA in R2.2;
    # legacy `peaks` rows are handled by migrate_r2_split_peaks! instead.
    schema_tables = let db_tmp = SQLite.DB()
        create_schema!(db_tmp)
        Set(String[String(r.name) for r in Tables.rowtable(DBInterface.execute(db_tmp,
            "SELECT name FROM sqlite_master WHERE type='table'"))])
    end
    tables_to_migrate = filter(t ->
        !isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", [t]))) &&
        t in schema_tables,
        tables)

    # FK enforcement must be disabled OUTSIDE a transaction (SQLite docs).
    DBInterface.execute(db, "PRAGMA foreign_keys = OFF")
    try
        SQLite.transaction(db) do
            # Rename old tables and recreate them with AUTOINCREMENT.
            # IMPORTANT: we call create_schema! here ONLY to recreate the five
            # entity tables. We do NOT want create_schema! to also create tables
            # that have FK references to `exposures` within this transaction,
            # because SQLite's ALTER TABLE RENAME tracking would corrupt those
            # FK references (storing "_migrate_old_exposures" instead of
            # "exposures" in the new table's schema). The deferred create_schema!
            # call in open_db (which already ran before migrate_schema!) handles
            # non-entity tables — they exist or will be created after this
            # transaction commits when the FK tracking state is clean.
            for t in tables_to_migrate
                DBInterface.execute(db, "ALTER TABLE $t RENAME TO _migrate_old_$t")
            end
            create_schema!(db)
            for t in tables_to_migrate
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
        # After the transaction commits, fix any FK references that were corrupted
        # by the ALTER TABLE RENAME tracking (SQLite updates all FK refs to the
        # renamed table name, including tables created INSIDE the transaction that
        # referenced `exposures` by name AFTER the rename — they end up stored as
        # REFERENCES "_migrate_old_exposures" which no longer exists).
        # Drop and recreate affected tables now that the rename transaction is done.
        _fix_fk_references_after_autoincrement_migration!(db)
    finally
        DBInterface.execute(db, "PRAGMA foreign_keys = ON")
    end
end

"""
    _fix_fk_references_after_autoincrement_migration!(db)

After `migrate_pk_to_autoincrement!` runs, SQLite's ALTER TABLE RENAME tracking
may have stored corrupted FK references (pointing to `_migrate_old_exposures`
instead of `exposures`) in tables created by `create_schema!` INSIDE the rename
transaction. Fix by dropping and recreating those tables. This is safe because
these tables were freshly created (empty) by create_schema! in the migration.
"""
function _fix_fk_references_after_autoincrement_migration!(db::SQLite.DB)
    # Find tables with FK references to non-existent tables (the renamed ones).
    # Only check tables we created in this run (non-entity new R2.1 tables).
    candidates = ["auto_peaks", "peak_curations"]
    for t in candidates
        # Check if the table exists at all.
        exists = !isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", [t])))
        exists || continue
        # Check if any FK points to a non-existent table.
        fks = Tables.rowtable(DBInterface.execute(db, "PRAGMA foreign_key_list($t)"))
        needs_fix = any(fk -> begin
            ref_table = String(fk.table)
            isempty(Tables.rowtable(DBInterface.execute(db,
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", [ref_table])))
        end, fks)
        needs_fix || continue
        # Drop and recreate: the table was just created (empty) so no data is lost.
        @debug "_fix_fk_references_after_autoincrement_migration!: rebuilding $t"
        DBInterface.execute(db, "DROP TABLE $t")
    end
    # Recreate all affected tables via create_schema! (IF NOT EXISTS is safe).
    create_schema!(db)
end

"""
    migrate_r2_widen_index_peaks_pk!(db)

Rebuild `index_peaks` to widen its PRIMARY KEY from `(index_id, peak_id)` to
`(index_id, peak_id, peak_kind)` and add the `peak_kind` discriminator column.
Idempotent: no-op if the column is already part of the PK or if the table
doesn't exist yet (fresh DBs created from the new SCHEMA already have the
right shape).
"""
function migrate_r2_widen_index_peaks_pk!(db::SQLite.DB)
    # Sentinel: skip if peak_kind already in PK (rebuilt previously).
    info = Tables.rowtable(DBInterface.execute(db,
        "SELECT name, pk FROM pragma_table_info('index_peaks')"))
    already_widened = any(r -> String(r.name) == "peak_kind" && Int(r.pk) > 0, info)
    already_widened && return
    # Skip on fresh DBs that already have the new shape via CREATE.
    isempty(info) && return

    DBInterface.execute(db, "PRAGMA foreign_keys = OFF")
    try
        SQLite.transaction(db) do
            DBInterface.execute(db, "ALTER TABLE index_peaks RENAME TO _index_peaks_old")
            # Re-create with the new shape (matches SCHEMA above).
            DBInterface.execute(db, """
                CREATE TABLE index_peaks (
                    index_id       INTEGER REFERENCES indices(id),
                    peak_id        INTEGER NOT NULL,
                    peak_kind      TEXT NOT NULL DEFAULT 'auto'
                                   CHECK (peak_kind IN ('auto', 'curation')),
                    ratio_position INTEGER,
                    residual       REAL,
                    PRIMARY KEY (index_id, peak_id, peak_kind)
                )
            """)
            # Old rows are all 'auto' (this runs before migrate_r2_split_peaks!,
            # which is responsible for repointing manual-peak refs).
            DBInterface.execute(db, """
                INSERT INTO index_peaks (index_id, peak_id, peak_kind, ratio_position, residual)
                SELECT index_id, peak_id, 'auto', ratio_position, residual
                FROM _index_peaks_old
            """)
            DBInterface.execute(db, "DROP TABLE _index_peaks_old")
        end
    finally
        DBInterface.execute(db, "PRAGMA foreign_keys = ON")
    end
end

"""
    migrate_r2_split_peaks!(db)

Backfill `auto_peaks` and `peak_curations` from the legacy `peaks` table,
repointing `index_peaks.peak_id` for manual-peak references so user-built
speculatives survive the migration. Idempotent: returns early if `peaks`
no longer exists.
"""
function migrate_r2_split_peaks!(db::SQLite.DB)
    # Sentinel: if `peaks` table is gone, migration already ran (or was never
    # needed). This is the normal state for all R2.1+ DBs after first run.
    peaks_exists = !isempty(Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='peaks'")))
    peaks_exists || return

    # Sentinel: distinguish fresh R2.1 DBs (peaks exists but was never written
    # to by the pipeline) from legacy DBs (peaks has data from pre-R2.1 use).
    # Two checks:
    # (a) sqlite_sequence: AUTOINCREMENT tables get an entry here after first
    #     INSERT — catches R2.1 DBs where peaks had AUTOINCREMENT.
    # (b) Direct row count: pre-R2.1 DBs used plain INTEGER PRIMARY KEY (no
    #     AUTOINCREMENT), so sqlite_sequence has no entry even if rows exist —
    #     fall back to a COUNT(*) check.
    # Fresh DBs need no migration; peaks will be removed in R2.2.
    peaks_ever_written = !isempty(Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM sqlite_sequence WHERE name = 'peaks'"))) ||
        (first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM peaks"))).n > 0)
    peaks_ever_written || return

    # Sentinel: if auto_peaks already has rows, we're partway through —
    # the only safe action is to bail and require operator intervention.
    auto_peaks_exists = !isempty(Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='auto_peaks'")))
    if auto_peaks_exists
        auto_count = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM auto_peaks"))).n
        if Int(auto_count) > 0
            error("migrate_r2_split_peaks!: auto_peaks already has $auto_count rows " *
                  "but peaks table still exists — operator intervention required " *
                  "(restore from backup or manually reconcile)")
        end
    end

    # At this point: `peaks` exists and has been written to (legacy DB).
    # auto_peaks is either absent (true pre-R2.1) or empty (post-create_schema!).
    # Create destination tables if they don't exist yet.

    DBInterface.execute(db, "PRAGMA foreign_keys = OFF")
    try
        SQLite.transaction(db) do
            # Ensure destination tables exist (true pre-R2.1 DBs won't have them).
            DBInterface.execute(db, """
                CREATE TABLE IF NOT EXISTS auto_peaks (
                    id              INTEGER PRIMARY KEY AUTOINCREMENT,
                    exposure_id     INTEGER REFERENCES exposures(id),
                    q               REAL NOT NULL,
                    intensity       REAL,
                    prominence      REAL,
                    sharpness       REAL,
                    findpeaks_index INTEGER
                )
            """)
            DBInterface.execute(db, """
                CREATE TABLE IF NOT EXISTS peak_curations (
                    id          INTEGER PRIMARY KEY AUTOINCREMENT,
                    exposure_id INTEGER REFERENCES exposures(id),
                    kind        TEXT NOT NULL CHECK (kind IN ('exclude', 'add')),
                    q           REAL NOT NULL,
                    created_by  INTEGER REFERENCES users(id) ON DELETE SET NULL,
                    created_at  DATETIME DEFAULT CURRENT_TIMESTAMP
                )
            """)

            # Introspect the legacy `peaks` table to handle minimal/partial schemas
            # (e.g. test DBs that only have id, exposure_id, q).
            peaks_cols = Set(String.(
                Tables.rowtable(DBInterface.execute(db, "PRAGMA table_info(peaks)")) .|>
                r -> r.name))
            has_source     = "source"     ∈ peaks_cols
            has_intensity  = "intensity"  ∈ peaks_cols
            has_prominence = "prominence" ∈ peaks_cols
            has_sharpness  = "sharpness"  ∈ peaks_cols
            has_excluded   = "excluded"   ∈ peaks_cols

            intensity_sel  = has_intensity  ? "intensity"  : "NULL"
            prominence_sel = has_prominence ? "prominence" : "NULL"
            sharpness_sel  = has_sharpness  ? "sharpness"  : "NULL"
            auto_where     = has_source     ? "WHERE source = 'auto'" : ""
            excl_where     = has_source && has_excluded ?
                "WHERE source = 'auto' AND excluded = 1" : "WHERE 1=0"
            manual_where   = has_source     ? "WHERE source = 'manual'" : "WHERE 1=0"

            # 1. Auto peaks: id preserved (peaks PK was AUTOINCREMENT).
            # findpeaks_index left NULL for legacy rows — synthesize_peaks_result
            # falls back to argmin lookup when NULL; the next analyze run that
            # invokes diff_update_auto_peaks! will populate it.
            DBInterface.execute(db, """
                INSERT INTO auto_peaks (id, exposure_id, q, intensity, prominence, sharpness, findpeaks_index)
                SELECT id, exposure_id, q, $intensity_sel, $prominence_sel, $sharpness_sel, NULL
                FROM peaks $auto_where
            """)

            # 2. Exclusion curations: q-value is the binding key.
            DBInterface.execute(db, """
                INSERT INTO peak_curations (exposure_id, kind, q, created_by)
                SELECT exposure_id, 'exclude', q, NULL
                FROM peaks $excl_where
            """)

            # 3. Addition curations: row-by-row to capture old→new id mapping.
            manual_rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, exposure_id, q FROM peaks $manual_where"))
            old_to_new = Dict{Int, Int}()
            for r in manual_rows
                res = DBInterface.execute(db,
                    """INSERT INTO peak_curations (exposure_id, kind, q, created_by)
                       VALUES (?, 'add', ?, NULL)""",
                    [Int(r.exposure_id), Float64(r.q)])
                new_id = Int(DBInterface.lastrowid(res))
                old_to_new[Int(r.id)] = new_id
            end

            # 4. Repoint index_peaks rows whose peak_id was a manual peak.
            for (old_id, new_id) in old_to_new
                DBInterface.execute(db,
                    """UPDATE index_peaks SET peak_id = ?, peak_kind = 'curation'
                       WHERE peak_id = ?""",
                    [new_id, old_id])
            end

            # 5. Verify no orphan index_peaks rows remain (auto refs survive,
            #    manual refs were just repointed). Any remaining row whose
            #    peak_id doesn't resolve in either table is a bug.
            orphans = Tables.rowtable(DBInterface.execute(db, """
                SELECT ip.peak_id, ip.peak_kind, ip.index_id
                FROM index_peaks ip
                WHERE (ip.peak_kind = 'auto'     AND ip.peak_id NOT IN (SELECT id FROM auto_peaks))
                   OR (ip.peak_kind = 'curation' AND ip.peak_id NOT IN (SELECT id FROM peak_curations))
            """))
            if !isempty(orphans)
                error("migrate_r2_split_peaks!: $(length(orphans)) orphaned index_peaks " *
                      "rows after repoint — operator intervention required")
            end

            # 6. Mark all indices stale so the next analyze recomputes basis/score
            #    under the new effective_peaks model.
            DBInterface.execute(db, "UPDATE indices SET status = 'stale'")

            # 7. Drop the old peaks table — fully decomposed and repointed.
            DBInterface.execute(db, "DROP TABLE peaks")
        end
        @info "migrate_r2_split_peaks! complete"
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
    preflight_index_groups_uniqueness!(db)
    create_schema!(db)
    migrate_schema!(db)
    # Flush any SQLite.jl statement cache entries that became stale due to
    # DDL operations (table renames/drops) in the migration functions.
    # Without this, the first user query after migration can fail with
    # "no such table: _migrate_old_*" because a cached prepared statement
    # from inside the migration transaction references a dropped table.
    SQLite.finalize_statements!(db)
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
