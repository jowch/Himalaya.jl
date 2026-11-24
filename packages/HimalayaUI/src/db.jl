using SQLite, DBInterface, Tables
using Dates: now, UTC, format, @dateformat_str

# Sentinel marker name for the speculative peak durability migration (2026-07-14).
const MIGRATION_SPECULATIVE_PEAK_DURABILITY = "speculative_peak_durability"

# Sentinel marker name for the Hexagonal √11 removal (#304).
const MIGRATION_HEX_SQRT11 = "hex_sqrt11_removal_v1"

const MIGRATION_RATIO_POSITION_INT = "index_peaks_ratio_position_int_v1"

# Sentinel marker name for the Plotting redesign Plan A durable-assignment backfill.
const MIGRATION_ASSIGNMENTS = "assignments_v1"

const MIGRATION_SAMPLE_CALLS = "sample_calls_v1"

# Sentinel marker names for the ingestion-redesign Phase A schema migrations.
const MIGRATION_LOADS_TABLE             = "loads_table_v1"
const MIGRATION_EXPOSURES_EXPERIMENT_ID = "exposures_experiment_id_v1"
const MIGRATION_EXPERIMENTS_GEOMETRY    = "experiments_geometry_v1"
const MIGRATION_SAMPLES_NAME_COLLAPSE   = "samples_name_collapse_v1"

# Sentinel marker name for the ingestion-redesign Phase E1 editable config columns.
const MIGRATION_EXPERIMENTS_CONFIG_COLS = "experiments_config_cols_v1"

# Sentinel marker name for the session_id backfill on existing loads.
const MIGRATION_LOAD_SESSIONS_BACKFILL = "load_sessions_backfill_v1"

const MIGRATION_DROP_SCAN_SIGNATURE = "drop_scan_signature_v1"

const MIGRATION_EXPOSURES_BEAM_MONITORS = "exposures_beam_monitors_v1"
const MIGRATION_EXPOSURES_PIPE_LENGTH = "exposures_pipe_length_v1"
const MIGRATION_EXPOSURES_Q_STEP = "exposures_q_step_v1"
const MIGRATION_EXPOSURES_DESELECTED = "exposures_deselected_v1"
const MIGRATION_EXPOSURES_TRACE_MEDIAN = "exposures_trace_median_v1"
const MIGRATION_BEAM_OFF_DESELECTED = "beam_off_deselected_v1"

const MIGRATION_EXPOSURES_INGEST_STATE = "exposures_ingest_state_v1"

const MIGRATION_EXPERIMENTS_INGEST_ERROR = "experiments_ingest_error_v1"
const MIGRATION_LABEL_SUGGESTIONS = "label_suggestions_v1"
const MIGRATION_BRIDGED_LABELS = "bridged_labels_shared_v1"
const MIGRATION_BRIDGED_CALLS = "bridged_calls_shared_v1"
const MIGRATION_PEAK_PROVENANCE = "peak_provenance_v1"

# Sentinel helpers (ingestion redesign): gate-read and marker-write for the
# `schema_migrations` table.
_migrated(db::SQLite.DB, name::AbstractString) = !isempty(sqlrows(DBInterface.execute(db,
    "SELECT 1 FROM schema_migrations WHERE name = ?", [name])))
_table_exists(db::SQLite.DB, name::AbstractString) = !isempty(sqlrows(DBInterface.execute(db,
    "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?", [name])))

_record_migration!(db::SQLite.DB, name::AbstractString) = DBInterface.execute(db,
    "INSERT INTO schema_migrations (name, applied_at) VALUES (?, ?)", [name, now_iso()])

# Column names of a table — lets ADD COLUMN migrations stay idempotent on a
# partially-migrated DB by skipping columns that already exist.
cols_of(db, table) = String.(getproperty.(sqlrows(
    DBInterface.execute(db, "PRAGMA table_info($table)")), :name))

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
    created_at      DATETIME DEFAULT CURRENT_TIMESTAMP,
    description     TEXT,
    last_scanned_at TEXT,
    ingest_status   TEXT DEFAULT 'idle',
    ingest_error    TEXT
);

CREATE TABLE IF NOT EXISTS setups (
    id                   INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id        INTEGER NOT NULL REFERENCES experiments(id) ON DELETE CASCADE,
    folder               TEXT NOT NULL,
    display_name         TEXT,
    position             INTEGER NOT NULL DEFAULT 0,
    included             INTEGER NOT NULL DEFAULT 1,
    root_dir             TEXT NOT NULL,
    data_dir             TEXT NOT NULL,
    analysis_dir         TEXT NOT NULL,
    image_pattern        TEXT,
    metadata_pattern     TEXT,
    integration_pattern  TEXT,
    energy_kev           REAL,
    energy_kev_source    TEXT DEFAULT 'default',
    flight_path_m        REAL,
    flight_path_m_source TEXT DEFAULT 'default',
    beam_center_x        REAL,
    beam_center_x_source TEXT DEFAULT 'default',
    beam_center_y        REAL,
    beam_center_y_source TEXT DEFAULT 'default',
    pixel_size_um        REAL,
    pixel_size_um_source TEXT DEFAULT 'default',
    q_units              TEXT,
    q_units_source       TEXT DEFAULT 'default',
    created_at           DATETIME DEFAULT CURRENT_TIMESTAMP,
    mode                 TEXT,
    mode_source          TEXT DEFAULT 'default'
);

CREATE TABLE IF NOT EXISTS loads (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id INTEGER NOT NULL REFERENCES experiments(id) ON DELETE CASCADE,
    load_index    INTEGER NOT NULL,
    session_id    INTEGER,
    start_time    TEXT,
    end_time      TEXT,
    frame_count   INTEGER NOT NULL DEFAULT 0,
    note          TEXT,
    setup_id      INTEGER REFERENCES setups(id)
);

CREATE TABLE IF NOT EXISTS samples (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id INTEGER REFERENCES experiments(id),
    name          TEXT,
    display_name  TEXT,
    notes         TEXT,
    load_id          INTEGER REFERENCES loads(id),
    slot_index       INTEGER,
    grouping_source  TEXT DEFAULT 'auto_position',
    name_source      TEXT DEFAULT 'auto',
    setup_id         INTEGER REFERENCES setups(id)
);

CREATE TABLE IF NOT EXISTS exposures (
    id                   INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id            INTEGER REFERENCES samples(id),
    filename             TEXT,
    kind                 TEXT DEFAULT 'file',
    status               TEXT CHECK (status IN ('accepted', 'rejected')),
    image_path           TEXT,
    trace_hash           TEXT,
    analysis_inputs_hash TEXT,
    experiment_id        INTEGER REFERENCES experiments(id) ON DELETE CASCADE,
    prp_path             TEXT,
    timestamp            TEXT,
    exposure_time        REAL,
    horizontal_position  REAL,
    scan_id              INTEGER,
    frame_no             INTEGER,
    load_id              INTEGER REFERENCES loads(id),
    content_fingerprint  TEXT,
    transmitted_counts   REAL,
    incident_monitor     REAL,
    counting_time        REAL,
    monitors_read        INTEGER NOT NULL DEFAULT 0,
    pipe_length_m        REAL,
    q_step               REAL,
    trace_median         REAL,
    analysis_version     TEXT,
    analysis_error       TEXT,
    analysis_trace_stat  TEXT,
    setup_id             INTEGER REFERENCES setups(id),
    reviewed_at          TEXT,
    peak_model_version   INTEGER,
    deselected           INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE IF NOT EXISTS exposure_sources (
    averaged_exposure_id INTEGER REFERENCES exposures(id),
    source_exposure_id   INTEGER REFERENCES exposures(id),
    role                 TEXT DEFAULT 'signal',
    PRIMARY KEY (averaged_exposure_id, source_exposure_id)
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
    kind        TEXT NOT NULL DEFAULT 'auto',
    inputs_hash TEXT,
    -- 1 = `basis`/`lattice_d` are a lattice the USER chose (the custom-index
    -- modal), not a fit. `_persist_analysis_inner!` re-resolves a locked
    -- index's peaks like any other speculative but skips the least-squares
    -- refit, so no reanalysis can move the lattice out from under the user.
    -- This is what makes insert_custom_index!'s scan-derived intents safe:
    -- they resolve WHICH peaks are claimed, never WHERE the comb sits.
    basis_locked INTEGER NOT NULL DEFAULT 0
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
    findpeaks_index INTEGER,
    found_by        TEXT NOT NULL DEFAULT 'hybrid' CHECK (found_by IN ('hybrid', 'model', 'both')),
    model_id        INTEGER
);

CREATE TABLE IF NOT EXISTS peak_curations (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    exposure_id INTEGER REFERENCES exposures(id),
    kind        TEXT NOT NULL CHECK (kind IN ('exclude', 'add')),
    q           REAL NOT NULL,
    created_by  INTEGER REFERENCES users(id) ON DELETE SET NULL,
    created_at  DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS peak_reviews (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    exposure_id    INTEGER NOT NULL REFERENCES exposures(id),
    event_id       INTEGER REFERENCES user_actions(id),
    created_at     TEXT NOT NULL,
    shown_detector TEXT NOT NULL,
    shown_model_id INTEGER,
    backfilled     INTEGER NOT NULL DEFAULT 0,
    direct         INTEGER NOT NULL DEFAULT 0,
    peaks_json     TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_peak_reviews_exposure
    ON peak_reviews(exposure_id);

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

CREATE TABLE IF NOT EXISTS assignments (
    exposure_id INTEGER PRIMARY KEY REFERENCES exposures(id),
    state       TEXT NOT NULL DEFAULT 'indexed'
                CHECK (state IN ('indexed', 'form_factor', 'null'))
);

CREATE TABLE IF NOT EXISTS assignment_members (
    exposure_id INTEGER NOT NULL REFERENCES exposures(id),
    index_id    INTEGER NOT NULL REFERENCES indices(id) ON DELETE CASCADE,
    PRIMARY KEY (exposure_id, index_id)
);

CREATE TABLE IF NOT EXISTS sample_calls (
    sample_id INTEGER PRIMARY KEY REFERENCES samples(id) ON DELETE CASCADE,
    state     TEXT NOT NULL DEFAULT 'indexed'
              CHECK (state IN ('indexed', 'form_factor', 'null'))
);

CREATE TABLE IF NOT EXISTS user_actions (
    id              INTEGER PRIMARY KEY,
    user_id         INTEGER REFERENCES users(id) ON DELETE SET NULL,
    timestamp       DATETIME DEFAULT CURRENT_TIMESTAMP,
    action          TEXT,
    entity_type     TEXT,
    entity_id       INTEGER,
    note            TEXT,
    payload         TEXT,
    undoes_event_id INTEGER REFERENCES user_actions(id),
    client_id       TEXT,
    client_op_id    TEXT
);

CREATE INDEX IF NOT EXISTS idx_events_by_exposure
    ON user_actions(entity_type, entity_id, id);

-- I2 partial unique index — installed by migrate_schema! after the legacy
-- ALTER TABLE pass that adds client_op_id, so it works on fresh and
-- legacy DBs alike. The SQL lives in migrate_schema! and reads
-- UNIQUE on the columns (client_op_id, action, entity_id), with a partial
-- WHERE client_op_id IS NOT NULL filter so legacy NULL-op_id rows are
-- excluded from the constraint.
--
-- Earlier drafts of the spec described this index as NOT UNIQUE, on the
-- premise that one request might emit multiple events under one
-- client_op_id. The implementation took the opposite path: each
-- with_idempotency-wrapped route is constrained to emit at most one
-- event row PER (action, entity_id) per request, which lets the unique
-- index serve as the idempotency-retry guard at the DB layer (a retry
-- with the same op_id trips the constraint and apply_event! short-
-- circuits to the prior event_id). The trade-off is documented as a
-- precondition for new routes — see CLAUDE.md Mutation queue section.

CREATE TABLE IF NOT EXISTS idempotent_responses (
    client_op_id  TEXT PRIMARY KEY,
    status_code   INTEGER NOT NULL,
    body          TEXT NOT NULL,
    created_at    DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_idempotent_responses_created
    ON idempotent_responses(created_at);

CREATE TABLE IF NOT EXISTS schema_migrations (
    name        TEXT PRIMARY KEY,
    applied_at  TEXT
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
    has_table = !isempty(sqlrows(DBInterface.execute(db,
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='index_groups'")))
    has_table || return
    dups = sqlrows(DBInterface.execute(db, """
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
    # Each ALTER TABLE adds a column that may already exist on legacy DBs;
    # SQLite has no `ADD COLUMN IF NOT EXISTS`, so we tolerate "duplicate column"
    # specifically. CREATE INDEX/TABLE statements below are `IF NOT EXISTS`-
    # guarded and any unexpected failure must surface — never swallow them
    # under a bare catch (issue #6 from PR review).
    alter_stmts = [
        "ALTER TABLE exposures ADD COLUMN status TEXT CHECK (status IN ('accepted', 'rejected'))",
        "ALTER TABLE exposures ADD COLUMN image_path TEXT",
        "ALTER TABLE users ADD COLUMN first_name TEXT",
        "ALTER TABLE users ADD COLUMN last_name TEXT",
        "ALTER TABLE indices ADD COLUMN kind TEXT NOT NULL DEFAULT 'auto'",
        "ALTER TABLE exposures ADD COLUMN trace_hash TEXT",
        "ALTER TABLE exposures ADD COLUMN analysis_inputs_hash TEXT",
        "ALTER TABLE indices ADD COLUMN inputs_hash TEXT",
        # Legacy rows default to 0 (unlocked) — correct: every index that
        # predates the column was either auto or anchor-snap speculative, both
        # of which are meant to refit. Custom commits set it explicitly.
        "ALTER TABLE indices ADD COLUMN basis_locked INTEGER NOT NULL DEFAULT 0",
        "ALTER TABLE user_actions ADD COLUMN payload TEXT",
        "ALTER TABLE user_actions ADD COLUMN undoes_event_id INTEGER REFERENCES user_actions(id)",
        "ALTER TABLE user_actions ADD COLUMN client_id TEXT",
        "ALTER TABLE user_actions ADD COLUMN client_op_id TEXT",
    ]
    for stmt in alter_stmts
        try
            DBInterface.execute(db, stmt)
        catch err
            # Lowercase the message before matching so a future SQLite or
            # SQLite.jl change in casing/prefix doesn't silently flip a
            # tolerated error into a propagated one (PR review suggestion #9).
            msg = lowercase(sprint(showerror, err))
            # Two errors are expected during incremental upgrades:
            # - "duplicate column name": the column already exists on a DB
            #   that's been migrated before. Idempotent — ignore.
            # - "no such table": a partial legacy DB or fresh-from-create_schema!
            #   that hasn't yet seen this table (e.g. older test fixtures
            #   pre-dating user_actions). Idempotent — ignore.
            # Anything else (typo, FK clash, constraint violation) must
            # propagate so a half-migrated DB isn't silently masked (issue #6).
            (occursin("duplicate column name", msg) || occursin("no such table", msg)) ||
                rethrow()
        end
    end

    migrate_setups!(db; finish = false)

    # The below statements are all IF NOT EXISTS-guarded. The CREATE INDEX
    # variants depend on `user_actions` existing; partial-legacy DBs (e.g.
    # the test fixture in test_db.jl that only has `experiments`) tolerate
    # "no such table" — every other failure must propagate (issue #6).
    function _create_safely(stmt::String)
        try
            DBInterface.execute(db, stmt)
        catch err
            occursin("no such table", lowercase(sprint(showerror, err))) || rethrow()
        end
    end
    _create_safely(
        "CREATE INDEX IF NOT EXISTS idx_events_by_client_op_id ON user_actions(client_op_id) WHERE client_op_id IS NOT NULL")
    DBInterface.execute(db,
        """CREATE TABLE IF NOT EXISTS idempotent_responses (
            client_op_id  TEXT PRIMARY KEY,
            status_code   INTEGER NOT NULL,
            body          TEXT NOT NULL,
            created_at    DATETIME DEFAULT CURRENT_TIMESTAMP
        )""")
    DBInterface.execute(db,
        "CREATE INDEX IF NOT EXISTS idx_idempotent_responses_created ON idempotent_responses(created_at)")
    # Precondition: legacy DBs that already populated `client_op_id` without
    # also installing this UNIQUE index AND that wrote duplicate rows for the
    # same (client_op_id, action, entity_id) tuple would fail this CREATE
    # with `UNIQUE constraint failed`. _create_safely only tolerates "no such
    # table" — a duplicate-constraint failure here propagates intentionally,
    # so a corrupt DB can't silently mask the inconsistency. There is no
    # released version of that partial-deploy state; should one ever exist,
    # operators must dedupe `user_actions` by (client_op_id, action,
    # entity_id) before re-running open_db (suggestion #13).
    _create_safely(
        """CREATE UNIQUE INDEX IF NOT EXISTS idx_events_unique_op
            ON user_actions(client_op_id, action, entity_id)
            WHERE client_op_id IS NOT NULL""")
    # (Issue #88) Run BEFORE the AUTOINCREMENT rebuild so the rebuild's column-copy
    # step sees the canonical (name, display_name) shape on both sides — preserving
    # stable identifiers on pre-Plan-7 DBs that lack AUTOINCREMENT.
    migrate_samples_naming!(db)
    migrate_pk_to_autoincrement!(db)
    # Heal corrupted FKs from prior runs: a DB migrated before the helper
    # was generalized may still have FK refs to `_migrate_old_*` on tables
    # like exposure_sources/index_groups/etc. Idempotent — no-op when clean.
    _fix_fk_references_after_autoincrement_migration!(db)
    # Run AFTER `_fix_fk_references_after_autoincrement_migration!` so the
    # AUTOINCREMENT rebuild has settled and the index attaches to the rebuilt
    # `exposures` table — placing it earlier would have it dropped along with
    # `_migrate_old_exposures` during `migrate_pk_to_autoincrement!`.
    migrate_exposures_unique_filename!(db)
    migrate_r2_widen_index_peaks_pk!(db)  # rebuild with widened PK first
    migrate_r2_split_peaks!(db)            # then repoint manual-peak refs

    # R2.3 sentinel: belt-and-suspenders verifier. `migrate_r2_split_peaks!`
    # is responsible for dropping the legacy `peaks` table; if it's still
    # around at this point, something went wrong in the migration path
    # (e.g. partial state from an aborted run). Warn loudly so operators
    # see it in the daemon log.
    legacy = sqlrows(DBInterface.execute(db,
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='peaks'"))
    if !isempty(legacy)
        @warn "Legacy 'peaks' table still present after R2 migration — investigate"
    end

    # R3.2: Once analysis_inputs_hash is the source of truth for staleness,
    # the old 'stale' enum value is dead. Any existing 'stale' rows came from
    # R2's transitional UPDATE; normalize them to 'candidate' (the next analyze
    # run will set inputs_hash; until then, hash mismatch with NULL on the index
    # renders them as stale to the UI). Unconditional — idempotent (no-op when
    # no 'stale' rows exist).
    DBInterface.execute(db, "UPDATE indices SET status = 'candidate' WHERE status = 'stale'")

    # R4.2: Rewrite legacy entity_type='peak'/'index' rows to 'exposure' so the
    # idx_events_by_exposure index is useful for fold-by-exposure queries.
    # Must run AFTER R2 migrations (which create auto_peaks and drop legacy peaks).
    migrate_r4_rebase_entity_type!(db)

    # Plotting redesign Plan A: durable per-exposure assignment. MUST run AFTER
    # migrate_pk_to_autoincrement! (db.jl, which rebuilds exposures/indices — the
    # assignments/assignment_members FKs must point at the rebuilt tables) and
    # after create_schema! (the assignment tables + the
    # schema_migrations sentinel exist). Backfills from the legacy active group;
    # sentinel-gated; own transaction.
    migrate_assignments!(db)

    # Ingestion redesign (Phase A): additive schema scaffolding. Sentinel-gated
    # stubs for now — later tasks fill in the real DDL. Run AFTER
    # migrate_samples_naming!, in dependency order:
    # loads → exposures denorm → experiments geometry → samples collapse.
    migrate_loads_table!(db)
    migrate_exposures_experiment_id!(db)
    migrate_experiments_geometry!(db)
    migrate_samples_name_collapse!(db)

    migrate_drop_merged_into!(db)

    # Ingestion redesign (Phase E1): editable config columns on experiments
    # (description + 3 pattern fields). Additive; no dependency on Phase D.
    migrate_experiments_config_cols!(db)

    # Session-id backfill: assign macro-session indices to all existing loads whose
    # session_id is NULL. Run AFTER migrate_loads_table! (the loads table must exist).
    backfill_load_sessions!(db)

    migrate_drop_scan_signature!(db)
    migrate_exposures_beam_monitors!(db)
    migrate_exposures_ingest_state!(db)
    migrate_experiments_ingest_error!(db)

    # Speculative peak durability (2026-07-14 spec): intents table + legacy
    # basis rescale + heal trigger. Must run last — see the function docstring
    # for the three ordering constraints (after the PK/AUTOINCREMENT rebuild +
    # FK-heal, after migrate_r2_split_peaks!).
    migrate_speculative_peak_durability!(db)

    # #304: √11 removed from phaseratios(Hexagonal). Renumbers the durable
    # ratio_position store and reopens the indexpeaks memoization gate. Must run
    # after migrate_speculative_peak_durability! (needs speculative_peak_intents).
    migrate_hex_sqrt11!(db)

    migrate_index_peaks_ratio_int!(db)
    migrate_peak_provenance_schema!(db)
    migrate_sample_bridges!(db)
    migrate_sample_compares!(db)
    migrate_sample_labels!(db)
    migrate_label_colors!(db)
    migrate_sample_calls!(db)
    migrate_drop_exposure_selected!(db)
    migrate_drop_retired_tables!(db)
    migrate_label_suggestions!(db)
    migrate_share_bridged_labels!(db)
    migrate_share_bridged_calls!(db)
    migrate_setups!(db)
    migrate_exposures_pipe_length!(db)
    migrate_exposures_q_step!(db)
    migrate_exposures_trace_median!(db)
    migrate_peak_provenance!(db)
    migrate_peak_models!(db)
    migrate_exposures_deselected!(db)
end

function migrate_exposures_deselected!(db::SQLite.DB)
    _migrated(db, MIGRATION_EXPOSURES_DESELECTED) && return nothing
    SQLite.transaction(db) do
        "deselected" in cols_of(db, "exposures") ||
            DBInterface.execute(db, "ALTER TABLE exposures ADD COLUMN deselected INTEGER NOT NULL DEFAULT 0")
        _select_worked_images!(db)
        _record_migration!(db, MIGRATION_EXPOSURES_DESELECTED)
    end
    nothing
end

function _select_worked_images!(db::SQLite.DB)
    kinds = join(("'" * k * "'" for k in DIRECT_REVIEW_KINDS), ",")
    actions = _table_exists(db, "user_actions") ? """
        OR EXISTS (SELECT 1 FROM user_actions a WHERE a.entity_type = 'exposure'
                    AND a.entity_id = e.id AND a.action IN ($kinds))""" : ""
    worked = sqlrows(DBInterface.execute(db, """
        SELECT e.sample_id AS sample_id,
               COALESCE(MAX(CASE WHEN EXISTS (SELECT 1 FROM assignment_members m WHERE m.exposure_id = e.id)
                                 THEN e.id END), MAX(e.id)) AS id
          FROM exposures e
         WHERE e.kind = 'file' AND e.sample_id IS NOT NULL
           AND (e.status IS NULL OR e.status != 'rejected')
           AND (EXISTS (SELECT 1 FROM assignment_members m WHERE m.exposure_id = e.id) $actions)
         GROUP BY e.sample_id"""))
    for w in worked
        DBInterface.execute(db,
            "UPDATE exposures SET deselected = CASE WHEN id = ? THEN 0 ELSE 1 END WHERE sample_id = ? AND kind = 'file'",
            [Int(w.id), Int(w.sample_id)])
    end
    nothing
end

function migrate_label_colors!(db::SQLite.DB)
    has(t) = !isempty(sqlrows(DBInterface.execute(db,
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", [t])))
    if has("tag_colors") && !has("label_colors")
        DBInterface.execute(db, "ALTER TABLE tag_colors RENAME TO label_colors")
    end
    DBInterface.execute(db, "DROP INDEX IF EXISTS idx_tag_colors_label")
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS label_colors (
            label TEXT NOT NULL,
            color TEXT NOT NULL
        )""")
    DBInterface.execute(db,
        "CREATE UNIQUE INDEX IF NOT EXISTS idx_label_colors_label ON label_colors(label COLLATE NOCASE)")
    nothing
end

function migrate_sample_labels!(db::SQLite.DB)
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS sample_labels (
            id        INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id INTEGER NOT NULL REFERENCES samples(id) ON DELETE CASCADE,
            label     TEXT NOT NULL,
            state     TEXT NOT NULL CHECK (state IN ('suggested', 'confirmed', 'removed')),
            source    TEXT NOT NULL DEFAULT 'user'
        )""")
    DBInterface.execute(db,
        "CREATE UNIQUE INDEX IF NOT EXISTS idx_sample_labels_unique ON sample_labels(sample_id, label COLLATE NOCASE)")
    nothing
end

const RETIRED_TABLES = (
    "series_pins", "series_messages", "series_samples", "series_members", "series",
    "comparison_pins", "comparison_messages", "comparison_members", "comparisons",
    "sample_messages", "sample_tags", "exposure_tags",
)

function migrate_drop_retired_tables!(db::SQLite.DB)
    SQLite.transaction(db) do
        for table in RETIRED_TABLES
            DBInterface.execute(db, "DROP TABLE IF EXISTS $table")
        end
    end
    nothing
end

now_iso() = format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ")

function migrate_sample_bridges!(db::SQLite.DB)
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS sample_bridges (
            sample_id INTEGER PRIMARY KEY REFERENCES samples(id) ON DELETE CASCADE,
            group_id  INTEGER NOT NULL
        )""")
    DBInterface.execute(db,
        "CREATE INDEX IF NOT EXISTS idx_sample_bridges_group ON sample_bridges(group_id)")
    nothing
end

function migrate_sample_compares!(db::SQLite.DB)
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS sample_compares (
            id            INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id     INTEGER NOT NULL REFERENCES samples(id) ON DELETE CASCADE,
            ref_sample_id INTEGER NOT NULL REFERENCES samples(id) ON DELETE CASCADE,
            added_by      INTEGER REFERENCES users(id) ON DELETE SET NULL
        )""")
    DBInterface.execute(db,
        "CREATE UNIQUE INDEX IF NOT EXISTS idx_sample_compares_pair ON sample_compares(sample_id, ref_sample_id)")
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS sample_subtractions (
            sample_id     INTEGER PRIMARY KEY REFERENCES samples(id) ON DELETE CASCADE,
            ref_sample_id INTEGER NOT NULL REFERENCES samples(id) ON DELETE CASCADE
        )""")
    nothing
end

# ── Ingestion redesign Phase A migrations ───────────────────────────────────
# Sentinel-gated stubs (Task 1). Later tasks replace the bodies with real DDL;
# each is idempotent via its `schema_migrations` sentinel.

function migrate_loads_table!(db::SQLite.DB)
    _migrated(db, MIGRATION_LOADS_TABLE) && return nothing
    SQLite.transaction(db) do
        DBInterface.execute(db, """
            CREATE TABLE IF NOT EXISTS loads (
                id           INTEGER PRIMARY KEY AUTOINCREMENT,
                experiment_id INTEGER NOT NULL REFERENCES experiments(id) ON DELETE CASCADE,
                load_index   INTEGER NOT NULL,
                session_id   INTEGER,
                start_time   TEXT,
                end_time     TEXT,
                frame_count  INTEGER NOT NULL DEFAULT 0,
                note         TEXT
            )
        """)
        DBInterface.execute(db,
            "CREATE INDEX IF NOT EXISTS loads_experiment_idx ON loads(experiment_id)")
        _record_migration!(db, MIGRATION_LOADS_TABLE)
    end
end

function migrate_exposures_experiment_id!(db::SQLite.DB)
    _migrated(db, MIGRATION_EXPOSURES_EXPERIMENT_ID) && return nothing
    SQLite.transaction(db) do
        existing = cols_of(db, "exposures")
        adds = [
            ("experiment_id", "INTEGER REFERENCES experiments(id) ON DELETE CASCADE"),
            ("prp_path", "TEXT"), ("timestamp", "TEXT"), ("exposure_time", "REAL"),
            ("horizontal_position", "REAL"), ("scan_id", "INTEGER"),
            ("frame_no", "INTEGER"), ("load_id", "INTEGER REFERENCES loads(id)"),
            ("content_fingerprint", "TEXT"),
        ]
        for (name, decl) in adds
            name in existing || DBInterface.execute(db,
                "ALTER TABLE exposures ADD COLUMN $name $decl")
        end

        # backfill experiment_id from the samples JOIN (sample_id may be NULL on some rows)
        DBInterface.execute(db, """
            UPDATE exposures
               SET experiment_id = (SELECT s.experiment_id FROM samples s WHERE s.id = exposures.sample_id)
             WHERE experiment_id IS NULL AND sample_id IS NOT NULL
        """)

        # FAIL-FAST: any exposure with no derivable experiment_id is a data error (spec §10 / P1-5)
        orphans = sqlrows(DBInterface.execute(db,
            "SELECT id, filename FROM exposures WHERE experiment_id IS NULL"))
        isempty(orphans) || error(
            "migrate_exposures_experiment_id!: $(length(orphans)) exposures have no derivable " *
            "experiment_id (orphaned sample_id). Resolve manually before migrating. ids: " *
            join(getproperty.(orphans, :id), ", "))

        # PREFLIGHT (dedupe-then-enforce, spec §10): the new key is (experiment_id, filename).
        # Collisions that were legal under the old (sample_id, filename) key would otherwise
        # make the CREATE UNIQUE INDEX below throw SQLite's terse "UNIQUE constraint failed"
        # on every open_db. Resolve them in two ways:
        #   * A group whose rows share EXACTLY ONE non-NULL image_path (np == 1) is a proven
        #     redundant duplicate (the real-data case: same physical frame re-attached under
        #     two samples). Keep one survivor and DELETE the rest (+ their FK children) — safe.
        #   * Any other group cannot be proven redundant, so we fail loudly rather than risk
        #     deleting real data. Two cases (both absent in real data, both fail the proof):
        #       - np > 1: rows span DISTINCT image_paths — genuinely-distinct files.
        #       - np == 0: ALL image_paths NULL — no evidence the rows are the same frame, so
        #         a silent delete could drop a distinct exposure. (np counts only non-NULLs.)
        # filename IS NOT NULL: a UNIQUE index treats multiple NULLs as distinct, so NULL
        # filenames are not collisions and must be excluded from the GROUP BY.
        dupes = sqlrows(DBInterface.execute(db, """
            SELECT experiment_id, filename,
                   COUNT(*) AS n, COUNT(DISTINCT image_path) AS np
              FROM exposures
             WHERE filename IS NOT NULL
             GROUP BY experiment_id, filename HAVING n > 1
        """))

        # Dedup only groups proven redundant by a single shared non-NULL image_path (np == 1).
        unproven = filter(d -> d.np != 1, dupes)
        isempty(unproven) || error(
            "migrate_exposures_experiment_id!: $(length(unproven)) (experiment_id, filename) groups " *
            "have duplicate exposures that cannot be proven redundant (DISTINCT image_paths, or " *
            "all image_paths NULL) — legal under the old (sample_id, filename) key. Manual merge " *
            "required before the unique index can be enforced. groups: " *
            join(["(exp=$(d.experiment_id), $(d.filename))" for d in unproven], ", "))

        # Same-image_path collisions: pick a survivor per group and delete the rest.
        # Survivor = analyzed-first (trace_hash present), lowest-id tiebreak.
        nonsurvivors = Int[]
        for d in dupes  # only np==1 (proven-redundant) groups remain (the rest errored above)
            rows = sqlrows(DBInterface.execute(db, """
                SELECT id FROM exposures
                 WHERE experiment_id IS ? AND filename = ?
                 ORDER BY (trace_hash IS NULL), id
            """, [d.experiment_id, d.filename]))
            append!(nonsurvivors, Int.(getproperty.(rows[2:end], :id)))  # all but the survivor
        end

        if !isempty(nonsurvivors)
            # Delete the non-survivors' FK children in strict child→parent order, scoped to
            # the non-survivor id set, BEFORE the exposure rows. FK enforcement is ON inside
            # this transaction (earlier migrations leave PRAGMA foreign_keys=ON and the pragma
            # cannot change inside an open transaction), so this order is load-bearing.
            # Mirrors the canonical delete order in routes_experiments.jl.
            idset = join(nonsurvivors, ", ")  # internal integer ids — safe to interpolate
            # grandchildren first: index_peaks / index_group_members via indices JOIN exposures
            DBInterface.execute(db, """
                DELETE FROM index_peaks WHERE index_id IN
                  (SELECT id FROM indices WHERE exposure_id IN ($idset))""")
            DBInterface.execute(db, """
                DELETE FROM index_group_members WHERE index_id IN
                  (SELECT id FROM indices WHERE exposure_id IN ($idset))""")
            # exposure-keyed children
            for tbl in ("assignment_members", "assignments", "index_groups",
                        "indices", "auto_peaks", "peak_curations")
                DBInterface.execute(db,
                    "DELETE FROM $tbl WHERE exposure_id IN ($idset)")
            end
            # exposure_sources keys two exposure columns; clean both.
            DBInterface.execute(db, """
                DELETE FROM exposure_sources
                 WHERE averaged_exposure_id IN ($idset)
                    OR source_exposure_id IN ($idset)""")
            _table_exists(db, "exposure_tags") && DBInterface.execute(db,
                "DELETE FROM exposure_tags WHERE exposure_id IN ($idset)")
            # finally the non-survivor exposure rows themselves.
            DBInterface.execute(db,
                "DELETE FROM exposures WHERE id IN ($idset)")
        end

        # swap the dedup key: drop old (sample_id, filename), add (experiment_id, filename)
        DBInterface.execute(db, "DROP INDEX IF EXISTS exposures_unique_filename")
        DBInterface.execute(db,
            "CREATE UNIQUE INDEX exposures_unique_filename ON exposures(experiment_id, filename)")
        DBInterface.execute(db,
            "CREATE INDEX IF NOT EXISTS exposures_experiment_idx ON exposures(experiment_id)")
        DBInterface.execute(db,
            "CREATE INDEX IF NOT EXISTS exposures_load_idx ON exposures(load_id)")
        _record_migration!(db, MIGRATION_EXPOSURES_EXPERIMENT_ID)
    end
end

function migrate_drop_scan_signature!(db::SQLite.DB)
    _migrated(db, MIGRATION_DROP_SCAN_SIGNATURE) && return nothing
    SQLite.transaction(db) do
        "scan_signature" in cols_of(db, "experiments") &&
            DBInterface.execute(db, "ALTER TABLE experiments DROP COLUMN scan_signature")
        _record_migration!(db, MIGRATION_DROP_SCAN_SIGNATURE)
    end
    nothing
end

function migrate_experiments_geometry!(db::SQLite.DB)
    _migrated(db, MIGRATION_EXPERIMENTS_GEOMETRY) && return nothing
    SQLite.transaction(db) do
        existing = cols_of(db, "experiments")
        adds = [
            ("beam_center_x", "REAL"), ("beam_center_y", "REAL"),
            ("pixel_size_um", "REAL"), ("q_units", "TEXT"),
            ("energy_kev_source", "TEXT DEFAULT 'default'"),
            ("flight_path_m_source", "TEXT DEFAULT 'default'"),
            ("beam_center_x_source", "TEXT DEFAULT 'default'"),
            ("beam_center_y_source", "TEXT DEFAULT 'default'"),
            ("pixel_size_um_source", "TEXT DEFAULT 'default'"),
            ("q_units_source", "TEXT DEFAULT 'default'"),
            ("last_scanned_at", "TEXT"),
            ("ingest_status", "TEXT DEFAULT 'idle'"),
            ("last_scan_tier", "TEXT DEFAULT 'fast'"),
            ("consecutive_empty_ticks", "INTEGER NOT NULL DEFAULT 0"),
        ]
        for (name, decl) in adds
            name in existing || DBInterface.execute(db,
                "ALTER TABLE experiments ADD COLUMN $name $decl")
        end
        _record_migration!(db, MIGRATION_EXPERIMENTS_GEOMETRY)
    end
    nothing
end

function migrate_samples_name_collapse!(db::SQLite.DB)
    _migrated(db, MIGRATION_SAMPLES_NAME_COLLAPSE) && return nothing
    SQLite.transaction(db) do
        existing = cols_of(db, "samples")
        # 1. grouping columns (additive, idempotent)
        for (name, decl) in [
                ("load_id", "INTEGER REFERENCES loads(id)"),
                ("slot_index", "INTEGER"),
                ("grouping_source", "TEXT DEFAULT 'auto_position'"),
                ("name_source", "TEXT DEFAULT 'auto'")]
            name in existing || DBInterface.execute(db, "ALTER TABLE samples ADD COLUMN $name $decl")
        end
        # 2. collapse the two text columns to one. Guard so the migration is idempotent and
        #    safe on DBs that never had both columns.
        if "display_name" in existing
            # DROP the unique index FIRST: SQLite re-points indexes onto a renamed column,
            # so a surviving samples_unique_name would re-impose label uniqueness post-rename.
            DBInterface.execute(db, "DROP INDEX IF EXISTS samples_unique_name")
            if "name" in existing
                DBInterface.execute(db, "ALTER TABLE samples DROP COLUMN name")
            end
            DBInterface.execute(db, "ALTER TABLE samples RENAME COLUMN display_name TO name")
        end
        # 3. NO new UNIQUE on the label.
        DBInterface.execute(db,
            "CREATE INDEX IF NOT EXISTS samples_experiment_idx ON samples(experiment_id)")
        DBInterface.execute(db,
            "CREATE INDEX IF NOT EXISTS samples_load_slot_idx ON samples(load_id, slot_index)")
        _record_migration!(db, MIGRATION_SAMPLES_NAME_COLLAPSE)
    end
    nothing
end

function migrate_drop_merged_into!(db::SQLite.DB)
    "merged_into_id" in cols_of(db, "samples") || return nothing
    SQLite.transaction(db) do
        DBInterface.execute(db, "ALTER TABLE samples DROP COLUMN merged_into_id")
    end
    nothing
end

function migrate_experiments_config_cols!(db::SQLite.DB)
    # Phase E1: editable config columns on experiments (description + 3 pattern fields).
    # Additive — safe to run against old DBs (ADD COLUMN does not rewrite existing rows).
    _migrated(db, MIGRATION_EXPERIMENTS_CONFIG_COLS) && return nothing
    SQLite.transaction(db) do
        existing = cols_of(db, "experiments")
        for (name, decl) in [
            ("description",         "TEXT"),
            ("image_pattern",       "TEXT"),
            ("metadata_pattern",    "TEXT"),
            ("integration_pattern", "TEXT"),
        ]
            name in existing || DBInterface.execute(db,
                "ALTER TABLE experiments ADD COLUMN $name $decl")
        end
        _record_migration!(db, MIGRATION_EXPERIMENTS_CONFIG_COLS)
    end
    nothing
end

function migrate_exposures_beam_monitors!(db::SQLite.DB)
    _migrated(db, MIGRATION_EXPOSURES_BEAM_MONITORS) && return nothing
    SQLite.transaction(db) do
        existing = cols_of(db, "exposures")
        for (name, decl) in (("transmitted_counts", "REAL"), ("incident_monitor", "REAL"),
                             ("counting_time", "REAL"), ("monitors_read", "INTEGER NOT NULL DEFAULT 0"))
            name in existing || DBInterface.execute(db,
                "ALTER TABLE exposures ADD COLUMN $name $decl")
        end
        _record_migration!(db, MIGRATION_EXPOSURES_BEAM_MONITORS)
    end
    nothing
end

function migrate_exposures_pipe_length!(db::SQLite.DB)
    _migrated(db, MIGRATION_EXPOSURES_PIPE_LENGTH) && return nothing
    SQLite.transaction(db) do
        "pipe_length_m" in cols_of(db, "exposures") ||
            DBInterface.execute(db, "ALTER TABLE exposures ADD COLUMN pipe_length_m REAL")
        DBInterface.execute(db,
            "UPDATE exposures SET monitors_read = 0 WHERE prp_path IS NOT NULL AND pipe_length_m IS NULL")
        _record_migration!(db, MIGRATION_EXPOSURES_PIPE_LENGTH)
    end
    nothing
end

function migrate_exposures_q_step!(db::SQLite.DB)
    _migrated(db, MIGRATION_EXPOSURES_Q_STEP) && return nothing
    SQLite.transaction(db) do
        "q_step" in cols_of(db, "exposures") ||
            DBInterface.execute(db, "ALTER TABLE exposures ADD COLUMN q_step REAL")
        _record_migration!(db, MIGRATION_EXPOSURES_Q_STEP)
    end
    nothing
end

function migrate_exposures_trace_median!(db::SQLite.DB)
    _migrated(db, MIGRATION_EXPOSURES_TRACE_MEDIAN) && return nothing
    SQLite.transaction(db) do
        "trace_median" in cols_of(db, "exposures") ||
            DBInterface.execute(db, "ALTER TABLE exposures ADD COLUMN trace_median REAL")
        _record_migration!(db, MIGRATION_EXPOSURES_TRACE_MEDIAN)
    end
    nothing
end

function migrate_exposures_ingest_state!(db::SQLite.DB)
    _migrated(db, MIGRATION_EXPOSURES_INGEST_STATE) && return nothing
    SQLite.transaction(db) do
        existing = cols_of(db, "exposures")
        for name in ("analysis_version", "analysis_error", "analysis_trace_stat")
            name in existing || DBInterface.execute(db, "ALTER TABLE exposures ADD COLUMN $name TEXT")
        end
        _record_migration!(db, MIGRATION_EXPOSURES_INGEST_STATE)
    end
    nothing
end

function migrate_experiments_ingest_error!(db::SQLite.DB)
    _migrated(db, MIGRATION_EXPERIMENTS_INGEST_ERROR) && return nothing
    SQLite.transaction(db) do
        "ingest_error" in cols_of(db, "experiments") ||
            DBInterface.execute(db, "ALTER TABLE experiments ADD COLUMN ingest_error TEXT")
        _record_migration!(db, MIGRATION_EXPERIMENTS_INGEST_ERROR)
    end
    nothing
end

"""
    backfill_load_sessions!(db)

Assign `session_id` to existing loads that still have NULL. For each experiment
with any NULL-session_id load, reads the ordered load start times, applies
`_assign_sessions` (spec §5: >3 h gap → new macro-session), and UPDATEs each
load row. Idempotent: experiments whose loads all already have a non-NULL
`session_id` are skipped entirely. Sentinel-gated.

Called from `migrate_schema!` after `migrate_loads_table!`.
"""
function backfill_load_sessions!(db::SQLite.DB)
    _migrated(db, MIGRATION_LOAD_SESSIONS_BACKFILL) && return nothing
    # Tolerate the loads table not yet existing (very old DBs where
    # migrate_loads_table! hasn't installed it — the sentinel below still fires
    # so a subsequent open_db doesn't re-run and hit the missing table again).
    has_loads = !isempty(sqlrows(DBInterface.execute(db,
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='loads'")))
    if has_loads
        # Experiments that have at least one load with NULL session_id.
        exp_ids = Int[Int(r.experiment_id) for r in sqlrows(DBInterface.execute(db,
            "SELECT DISTINCT experiment_id FROM loads WHERE session_id IS NULL"))]
        for eid in exp_ids
            # Read ALL of this experiment's loads (no NULL filter): session assignment
            # is holistic over the ordered start-times, and `_assign_sessions` is a pure
            # deterministic function of that full set — so re-deriving every load is safe
            # even if some already held a session_id (a mixed state no writer produces).
            rows = sqlrows(DBInterface.execute(db,
                "SELECT id, start_time FROM loads WHERE experiment_id = ? ORDER BY start_time",
                [eid]))
            isempty(rows) && continue
            # Parse start_time strings to DateTime (NULL → missing).
            parsed = Union{Dates.DateTime, Missing}[
                (r.start_time === nothing || r.start_time === missing) ?
                    missing : Dates.DateTime(String(r.start_time), "yyyy-mm-ddTHH:MM:SS")
                for r in rows]
            session_ids = _assign_sessions(parsed)
            SQLite.transaction(db) do
                for (row, sid) in zip(rows, session_ids)
                    DBInterface.execute(db,
                        "UPDATE loads SET session_id = ? WHERE id = ?",
                        [sid, Int(row.id)])
                end
            end
        end
    end
    SQLite.transaction(db) do
        _record_migration!(db, MIGRATION_LOAD_SESSIONS_BACKFILL)
    end
    nothing
end

"""
    migrate_speculative_peak_durability!(db)

One-shot durability migration for speculative indices (rules in
`packages/HimalayaUI/src/AGENTS.md`):

1. Creates `speculative_peak_intents` — the durable record of the user's
   ratio→q assignments, written at creation and never wiped by analysis.
   Deliberately NOT in the `SCHEMA` const: `create_schema!` runs before
   `migrate_pk_to_autoincrement!` rename-rebuilds `indices`, and SQLite's
   RENAME tracking would rewrite this table's FK to `_migrate_old_indices`
   (the same RENAME hazard `_fix_fk_references_after_autoincrement_migration!` heals).
2. Rescales legacy speculative `basis` values from the un-normalized fit
   convention to the normalized one (factor = first raw phase ratio: √2
   Pn3m/Im3m, √3 Fm3m/Fd3m, √6 Ia3d; 1 elsewhere — no-op).
3. Backfills intents from surviving `index_peaks` rows (best available proxy
   for the creation-time assignment). Rows whose peak was deleted join to
   NULL q and are skipped — orphaned before this table existed, unrecoverable.
4. NULLs `exposures.analysis_inputs_hash` where a speculative index has zero
   `index_peaks` rows, reopening the `indexpeaks_skipped` memoization gate so
   the next analyze runs the re-attach loop (the heal pass). One-shot repair,
   not a standing gate bypass.

Steps 2-4 + the sentinel marker run in one transaction, sentinel-gated.
Must be called at the END of `migrate_schema!`: after
`migrate_pk_to_autoincrement!` + FK-heal (indices table settled), after
`migrate_r2_split_peaks!` (backfill join needs the post-R2 peak tables).
"""
function migrate_speculative_peak_durability!(db::SQLite.DB)
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS speculative_peak_intents (
            index_id       INTEGER NOT NULL REFERENCES indices(id) ON DELETE CASCADE,
            ratio_position INTEGER NOT NULL,
            q              REAL    NOT NULL,
            PRIMARY KEY (index_id, ratio_position)
        )""")

    already = sqlrows(DBInterface.execute(db,
        "SELECT 1 FROM schema_migrations WHERE name = ?",
        [MIGRATION_SPECULATIVE_PEAK_DURABILITY]))
    isempty(already) || return nothing

    SQLite.transaction(db) do
        spec_rows = sqlrows(DBInterface.execute(db,
            "SELECT id, phase, basis FROM indices WHERE kind = 'speculative'"))
        for r in spec_rows
            ismissing(r.basis) && continue
            P = resolve_phase(String(r.phase))
            P === nothing && continue
            factor = Float64(first(Himalaya.phaseratios(P)))
            factor == 1.0 && continue
            DBInterface.execute(db,
                "UPDATE indices SET basis = ? WHERE id = ?",
                [Float64(r.basis) * factor, Int(r.id)])
        end

        # `ip.ratio_position IS NOT NULL`: the column is nullable in index_peaks
        # but NOT NULL in intents — one anomalous legacy row must become a
        # skipped row, not an aborted migration (the sentinel would never be
        # written and every subsequent open_db would throw).
        DBInterface.execute(db, """
            INSERT OR IGNORE INTO speculative_peak_intents (index_id, ratio_position, q)
            SELECT ip.index_id, ip.ratio_position, COALESCE(ap.q, pc.q)
            FROM index_peaks ip
            JOIN indices i ON i.id = ip.index_id AND i.kind = 'speculative'
            LEFT JOIN auto_peaks ap     ON ap.id = ip.peak_id AND ip.peak_kind = 'auto'
            LEFT JOIN peak_curations pc ON pc.id = ip.peak_id AND ip.peak_kind = 'curation'
                                        AND pc.kind = 'add'
            WHERE ip.ratio_position IS NOT NULL
              AND COALESCE(ap.q, pc.q) IS NOT NULL""")

        DBInterface.execute(db, """
            UPDATE exposures SET analysis_inputs_hash = NULL
             WHERE id IN (SELECT DISTINCT i.exposure_id FROM indices i
                           WHERE i.kind = 'speculative'
                             AND NOT EXISTS (SELECT 1 FROM index_peaks ip
                                              WHERE ip.index_id = i.id))""")

        # Plain INSERT (never OR IGNORE): the basis rescale above is NOT
        # idempotent under double application (×√2 twice = silent corruption),
        # and the sentinel gate is read outside this transaction. The PK
        # conflict on `schema_migrations.name` is the backstop that makes a
        # concurrent-open race roll back the losing racer's whole tx.
        DBInterface.execute(db,
            "INSERT INTO schema_migrations (name, applied_at) VALUES (?, ?)",
            [MIGRATION_SPECULATIVE_PEAK_DURABILITY, now_iso()])
    end
    nothing
end

"""
    migrate_hex_sqrt11!(db)

#304: `√11` was removed from `phaseratios(Hexagonal)` — it is not a permitted 2D
hexagonal reflection (`N = h²+hk+k²` has no integer solution for 11). The series
went 14 → 13 entries and every position past 5 renumbered, so persisted
`ratio_position`s for Hexagonal indices now point at the wrong reflection.

Three steps, sentinel-gated, in one transaction:

1. DELETE Hexagonal `index_peaks`. That table is DERIVED — `persist_analysis!`
   wipes it for every kind and rebuilds (auto from `indexpeaks`, speculative from
   `speculative_peak_intents`), so renumbering it would be pointless work. It is
   also not safely renumberable in SQL: auto rows written before #308 store
   `ratio_position` as a serialized `UInt8` BLOB (healed later by
   `migrate_index_peaks_ratio_int!`), so `ratio_position > 6` is a BLOB-vs-INTEGER
   affinity comparison that matches every blob row. Until the next analyze the index reads
   as claiming nothing — an already-supported state, and honest.
2. Renumber `speculative_peak_intents`, the only DURABLE `ratio_position` store
   (analyze reads it and never rewrites it). `p == 6` claimed the √11 that never
   existed → drop; `p > 6` shifts down one so each intent keeps its radicand.
3. NULL `analysis_inputs_hash` to reopen the `indexpeaks_skipped` memoization
   gate, so the next analyze re-indexes and re-scores.

Must run AFTER `migrate_speculative_peak_durability!` (needs the intents table).

**Refuses to run against a core that still has the √11.** `HimalayaUI`'s
`Project.toml` has no `[sources]` (it declares `julia = "1.9"`, and `[sources]`
is 1.11+), so the loaded `Himalaya` is whatever the gitignored manifest resolved
— which for a bare `Pkg.instantiate()` is the registry copy. Renumbering durable
intents into the 13-entry scheme and then resolving them against a 14-entry core
maps an intent for √12 at old position 7 onto position 6, which that core reads
as √11: silent wrong physics, and the burnt sentinel means it never self-corrects.
The `[compat] Himalaya = "0.6"` bound is the primary defence; this guard is the
backstop for a stale manifest that predates it. It returns WITHOUT writing the
sentinel, so a corrected environment still applies the migration later.

Rejected exposures are skipped by `analyze --all` (`cli.jl`), so their Hexagonal
`index_peaks` rows are deleted here with no rebuild path and stay claimless until
someone un-rejects them. Nothing unique is lost — auto rows are derived and
intents survive — but the rows do not come back on their own.
"""
function migrate_hex_sqrt11!(db::SQLite.DB;
                             series::AbstractVector = Himalaya.phaseratios(Himalaya.Hexagonal))
    _migrated(db, MIGRATION_HEX_SQRT11) && return nothing

    # `series` is injectable for one reason: this branch is the last thing
    # standing between a stale manifest and irreversible corruption of durable
    # intents, and it cannot be reached from a test otherwise (the test process
    # loads exactly one Himalaya). Production never passes it.
    if 11 in round.(Int, series .^ 2)
        @warn """
              migrate_hex_sqrt11! deferred: the loaded Himalaya still lists √11 in \
              phaseratios(Hexagonal), so renumbering now would point durable intents \
              at the wrong reflection. Resolve Himalaya >= 0.6 (dev the local core, \
              or update the manifest) and reopen the database.
              """ himalaya_series_length = length(series)
        return nothing
    end

    SQLite.transaction(db) do
        # `IN (...)` rather than `= 'Hexagonal'`: nothing in the schema constrains
        # this column, and `resolve_phase` strips a `Himalaya.` prefix precisely
        # because the qualified spelling is a known hazard. Every current writer
        # uses `string(nameof(P))` (bare), and prod is clean — but a qualified row
        # would be silently skipped while the sentinel burnt, which is exactly the
        # unrecoverable shape this migration must not have.
        DBInterface.execute(db, """
            DELETE FROM index_peaks
             WHERE index_id IN (SELECT id FROM indices
                                 WHERE phase IN ('Hexagonal', 'Himalaya.Hexagonal'))""")

        DBInterface.execute(db, """
            DELETE FROM speculative_peak_intents
             WHERE ratio_position = 6
               AND index_id IN (SELECT id FROM indices
                                 WHERE phase IN ('Hexagonal', 'Himalaya.Hexagonal'))""")

        # Two-pass via a +1000 offset. PRIMARY KEY (index_id, ratio_position) is
        # enforced per row and SQLite gives no ORDER BY on UPDATE, so a single
        # decrement can hit 8→7 while 7 still exists. A PK abort here would roll
        # back the sentinel and throw on every subsequent open_db.
        #
        # `typeof(...) = 'integer'` guards the same wedge from the other side.
        # Every intent is written as a plain Int today, but `index_peaks` proves
        # this column family can receive serialized UInt8 BLOBs (#308) — and under
        # BLOB affinity `> 6` matches every blob row while `+ 1000` collapses them
        # all to one value, so the PK aborts and open_db throws forever after.
        DBInterface.execute(db, """
            UPDATE speculative_peak_intents SET ratio_position = ratio_position + 1000
             WHERE typeof(ratio_position) = 'integer'
               AND ratio_position > 6
               AND index_id IN (SELECT id FROM indices
                                 WHERE phase IN ('Hexagonal', 'Himalaya.Hexagonal'))""")
        DBInterface.execute(db, """
            UPDATE speculative_peak_intents SET ratio_position = ratio_position - 1001
             WHERE typeof(ratio_position) = 'integer'
               AND ratio_position > 1000
               AND index_id IN (SELECT id FROM indices
                                 WHERE phase IN ('Hexagonal', 'Himalaya.Hexagonal'))""")

        # NOT scoped to exposures with a Hexagonal index: `remove_subsets` is
        # phase-agnostic (`issubset` keys on basis + peak set, not phase), so a
        # Hexagonal candidate whose score moved can now survive against a
        # same-basis candidate of another phase on an exposure that has no
        # Hexagonal row today — and eliminated candidates aren't persisted, so
        # that set can't be enumerated here.
        #
        # Zero-index exposures ARE provably safe to skip: elimination requires a
        # surviving eliminator, so zero persisted rows means nothing met
        # `minpeaks` at all, and dropping a ratio can only shrink a candidate's
        # peak count (every remaining ratio keeps its value, so no peak gains a
        # match it lacked). What failed `minpeaks` before still fails after.
        DBInterface.execute(db, """
            UPDATE exposures SET analysis_inputs_hash = NULL
             WHERE analysis_inputs_hash IS NOT NULL
               AND EXISTS (SELECT 1 FROM indices i WHERE i.exposure_id = exposures.id)""")

        # Plain INSERT (never OR IGNORE): the intent renumber is NOT idempotent
        # under double application (a second pass shifts everything down again),
        # and the sentinel gate is read outside this transaction. The PK conflict
        # on `schema_migrations.name` is the backstop that makes a concurrent-open
        # race roll back the losing racer's whole tx.
        DBInterface.execute(db,
            "INSERT INTO schema_migrations (name, applied_at) VALUES (?, ?)",
            [MIGRATION_HEX_SQRT11, now_iso()])
    end
    nothing
end

function _legacy_representative_id(exps)
    kept = [e for e in exps if is_kept_exposure(e)]
    isempty(kept) && return nothing
    picked = [e for e in kept if coalesce(e.selected, 0) == 1]
    maximum(e -> Int(e.id), isempty(picked) ? kept : picked)
end

function migrate_drop_exposure_selected!(db::SQLite.DB)
    _table_exists(db, "exposures") && "selected" in cols_of(db, "exposures") &&
        DBInterface.execute(db, "ALTER TABLE exposures DROP COLUMN selected")
    nothing
end

function migrate_sample_calls!(db::SQLite.DB)
    _migrated(db, MIGRATION_SAMPLE_CALLS) && return nothing
    SQLite.transaction(db) do
        last_event = Dict{Int, Int}()
        for r in sqlrows(DBInterface.execute(db, """
                SELECT entity_id AS eid, MAX(id) AS last FROM user_actions
                 WHERE entity_type = 'exposure'
                   AND action IN ('assignment_add', 'assignment_remove', 'assignment_set_state')
                 GROUP BY entity_id"""))
            last_event[Int(r.eid)] = Int(r.last)
        end
        members = Dict{Int, Vector{Int}}()
        for r in sqlrows(DBInterface.execute(db,
                "SELECT exposure_id AS eid, index_id AS iid FROM assignment_members ORDER BY index_id"))
            push!(get!(members, Int(r.eid), Int[]), Int(r.iid))
        end
        states = Dict{Int, String}()
        for r in sqlrows(DBInterface.execute(db,
                "SELECT exposure_id AS eid, state FROM assignments"))
            states[Int(r.eid)] = String(r.state)
        end
        carries(i) = haskey(members, i) || get(states, i, "indexed") != "indexed"
        by_sample = Dict{Int, Vector{NamedTuple}}()
        sel = "selected" in cols_of(db, "exposures") ? "selected" : "0 AS selected"
        for e in sqlrows(DBInterface.execute(db,
                "SELECT id, sample_id, $sel, status FROM exposures
                  WHERE sample_id IS NOT NULL ORDER BY id"))
            push!(get!(by_sample, Int(e.sample_id), NamedTuple[]), e)
        end
        for sid in sort!(collect(keys(by_sample)))
            exps = by_sample[sid]
            candidates = [Int(e.id) for e in exps if carries(Int(e.id))]
            isempty(candidates) && continue
            edited = [i for i in candidates if haskey(last_event, i)]
            rep = _legacy_representative_id(exps)
            keep = !isempty(edited) ? argmax(i -> last_event[i], edited) :
                   (rep !== nothing && rep in candidates) ? rep : minimum(candidates)
            for i in candidates, ix in get(members, i, Int[])
                i == keep && get(states, i, "indexed") == "indexed" && continue
                _log_migration_event!(db, "assignment_remove", "exposure", i, Dict(:index_id => ix))
            end
            _log_migration_event!(db, "sample_call_set_state", "sample", sid,
                                  Dict(:state => get(states, keep, "indexed")))
        end
        DBInterface.execute(db,
            "INSERT INTO schema_migrations (name, applied_at) VALUES (?, ?)",
            [MIGRATION_SAMPLE_CALLS, now_iso()])
    end
    nothing
end

function _log_migration_event!(db::SQLite.DB, kind::String, entity_type::String,
                               entity_id::Int, payload::Dict)
    payload_json = JSON3.write(payload)
    res = DBInterface.execute(db,
        """INSERT INTO user_actions (user_id, action, entity_type, entity_id, payload)
           VALUES (NULL, ?, ?, ?, ?)""",
        [kind, entity_type, entity_id, payload_json])
    update_view_for_event!(db, kind, entity_id, JSON3.read(payload_json),
                           Int(DBInterface.lastrowid(res)))
    nothing
end

function migrate_index_peaks_ratio_int!(db::SQLite.DB)
    _migrated(db, MIGRATION_RATIO_POSITION_INT) && return nothing
    unknown = SQLite.transaction(db) do
        rows = sqlrows(DBInterface.execute(db, """
            SELECT rowid AS rid,
                   hex(substr(ratio_position, 1, 3)) AS magic,
                   hex(substr(ratio_position, -1)) AS last
              FROM index_peaks
             WHERE typeof(ratio_position) = 'blob'"""))
        unknown = 0
        for r in rows
            if r.magic == "374A4C"
                DBInterface.execute(db,
                    "UPDATE index_peaks SET ratio_position = ? WHERE rowid = ?",
                    [parse(Int, r.last; base = 16), r.rid])
            else
                unknown += 1
            end
        end
        if unknown == 0
            DBInterface.execute(db,
                "INSERT INTO schema_migrations (name, applied_at) VALUES (?, ?)",
                [MIGRATION_RATIO_POSITION_INT, now_iso()])
        end
        unknown
    end
    unknown > 0 && @warn "migrate_index_peaks_ratio_int!: left index_peaks.ratio_position blobs it could not decode" count = unknown
    nothing
end

"""
    migrate_assignments!(db)

Backfill the durable per-exposure assignment from the legacy active group:
create an `assignments` row (state='indexed') and copy the active group's
members into `assignment_members`, for every exposure that has an active group.
Sentinel-gated, idempotent, own transaction. Raw INSERTs (never apply_event!) —
this is a data backfill, not a user action.

`active=1` matches WITHOUT a kind filter, by design: this captures both
user-confirmed (active custom) groups AND the auto group of a never-confirmed
exposure. This is a one-time HISTORICAL backfill — it preserves what a
pre-Plan-A DB already displayed at the moment of upgrade. NOTE: it deliberately
DIVERGES from current analyze behavior. A fresh analyze no longer seeds the
assignment from the auto group (auto-grouping is not a durable concept — see
`_persist_analysis_inner!`), so new exposures read as unindexed until a human
curates them; only legacy upgrades carry the auto guess forward. Consequence to
be aware of: after upgrade, an analyzed-but-never-confirmed legacy exposure
reports a non-null `confirmed_index` (the auto guess) where the legacy
`kind='custom'` snapshot reported `nothing`.

LOG-DERIVABILITY CAVEAT (load-bearing for disaster recovery): this backfill is
the ONLY record of a pre-Plan-A confirmation's membership — those exposures have
`index_confirmed` events (now no-op guards in `update_view_for_event!`) but no
`assignment_add` event, so their assignment is NOT reproducible from the event
log. In normal operation nothing is lost (`rebuild_views_from_log!` does not
truncate, and this migration's sentinel persists). But a true from-empty rebuild
(drop `assignments`/`assignment_members`, re-fold the log) MUST first clear the
`assignments_v1` sentinel and re-run this migration, or every pre-Plan-A
confirmation is silently lost. Post-Plan-A confirmations ride `assignment_add`
and round-trip through the log normally.
"""
function migrate_assignments!(db::SQLite.DB)
    already = sqlrows(DBInterface.execute(db,
        "SELECT 1 FROM schema_migrations WHERE name = ?", [MIGRATION_ASSIGNMENTS]))
    isempty(already) || return nothing

    SQLite.transaction(db) do
        DBInterface.execute(db,
            """INSERT OR IGNORE INTO assignments (exposure_id, state)
               SELECT DISTINCT exposure_id, 'indexed'
               FROM index_groups WHERE active = 1""")
        DBInterface.execute(db,
            """INSERT OR IGNORE INTO assignment_members (exposure_id, index_id)
               SELECT g.exposure_id, m.index_id
               FROM index_groups g
               JOIN index_group_members m ON m.group_id = g.id
               WHERE g.active = 1""")
        DBInterface.execute(db,
            "INSERT INTO schema_migrations (name, applied_at) VALUES (?, ?)",
            [MIGRATION_ASSIGNMENTS, now_iso()])
    end
    nothing
end

"""
    migrate_r4_rebase_entity_type!(db)

Rewrite historical user_actions rows with `entity_type='peak'` or `'index'`
to `entity_type='exposure'`, with `entity_id` resolved to the parent exposure.
Without this, the `idx_events_by_exposure` index is useless for fold-by-exposure
queries on legacy events.

Idempotent: returns early if no legacy rows remain.
"""
function migrate_r4_rebase_entity_type!(db::SQLite.DB)
    legacy = first(sqlrows(DBInterface.execute(db, """
        SELECT COUNT(*) AS n FROM user_actions
        WHERE entity_type IN ('peak', 'index')
    """))).n
    Int(legacy) == 0 && return

    SQLite.transaction(db) do
        # Resolve peak → exposure via auto_peaks/peak_curations split
        # (legacy `peaks` table no longer exists post-R2.2).
        DBInterface.execute(db, """
            UPDATE user_actions SET
                entity_type = 'exposure',
                entity_id = COALESCE(
                    (SELECT exposure_id FROM auto_peaks      WHERE auto_peaks.id     = user_actions.entity_id),
                    (SELECT exposure_id FROM peak_curations  WHERE peak_curations.id = user_actions.entity_id)
                )
            WHERE entity_type = 'peak'
              AND (
                EXISTS (SELECT 1 FROM auto_peaks     WHERE auto_peaks.id     = user_actions.entity_id)
                OR
                EXISTS (SELECT 1 FROM peak_curations WHERE peak_curations.id = user_actions.entity_id)
              )
        """)
        # Index events: resolve via indices.exposure_id.
        DBInterface.execute(db, """
            UPDATE user_actions SET
                entity_type = 'exposure',
                entity_id = (SELECT exposure_id FROM indices WHERE indices.id = user_actions.entity_id)
            WHERE entity_type = 'index'
              AND EXISTS (SELECT 1 FROM indices WHERE indices.id = user_actions.entity_id)
        """)
        # Stragglers (peak/index id no longer resolves) keep their original
        # entity_type — they won't appear in fold-by-exposure queries but
        # they're not lost; queryable by raw kind/id.
    end
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
# Single source of truth for the entity tables that participate in the
# AUTOINCREMENT migration. Used by both `migrate_pk_to_autoincrement!`
# (which renames+rebuilds them) and `_fix_fk_references_after_autoincrement_migration!`
# (which heals corrupted FK references that point at the rename-staging name).
# Keep these in lockstep — drift causes the heal loop to silently no-op for
# any new entity (review issue #19).
const _AUTOINCREMENT_ENTITIES = ("experiments", "samples", "exposures", "peaks", "indices")

# Superset of tables whose `_migrate_old_<name>` rename-staging name may
# appear in CREATE statements after a partial-failure migration. Includes
# `_AUTOINCREMENT_ENTITIES` plus `index_peaks` (R2 widen migration uses the
# same temp-name pattern; review issue #22). Heal loop iterates this list.
const _MIGRATION_TEMP_ENTITIES = (_AUTOINCREMENT_ENTITIES..., "index_peaks")

function migrate_pk_to_autoincrement!(db::SQLite.DB)
    tables = collect(_AUTOINCREMENT_ENTITIES)

    # Sentinel: skip iff every table in `tables` that EXISTS already has
    # AUTOINCREMENT. "peaks" may no longer exist in R2.2+ DBs (removed from
    # SCHEMA); skip it if absent so migrate_pk_to_autoincrement! stays
    # idempotent on fresh DBs.
    needs_migration = false
    for t in tables
        rows = sqlrows(DBInterface.execute(db,
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
        Set(String[String(r.name) for r in sqlrows(DBInterface.execute(db_tmp,
            "SELECT name FROM sqlite_master WHERE type='table'"))])
    end
    tables_to_migrate = filter(t ->
        !isempty(sqlrows(DBInterface.execute(db,
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
                new_cols = Set(String[String(c.name) for c in sqlrows(
                    DBInterface.execute(db, "PRAGMA table_info($t)"))])
                old_cols = String[String(c.name) for c in sqlrows(
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

After `migrate_pk_to_autoincrement!` runs, SQLite's ALTER TABLE RENAME
tracking rewrites FK references in EVERY table whose CREATE statement
named the renamed entity, including tables that were already populated.
The references point at `_migrate_old_<entity>` which was dropped at the
end of the rename transaction — leaving behind dead FK names that surface
later as `SQLiteException("no such table: main._migrate_old_*")` on
unrelated INSERTs (the prepare step walks the FK graph and trips when it
hits a stale reference).

The previous version of this function only rebuilt empty tables it knew
were freshly created (`auto_peaks`, `peak_curations`). That missed
populated tables like the since-dropped `sample_tags` and `exposure_tags`,
`exposure_sources`, `index_groups`, `index_group_members` whose FKs had
the same corruption. A real prod DB at `/tmp/himalaya-dev.db` with rows
in the since-dropped `sample_tags` and `sample_messages` exhibited the bug — caught only by
running the actual server (smoke checklist), not by the test suite which
uses freshly-built DBs that don't go through this exact sequence.

This implementation rewrites the FK references in-place via
`PRAGMA writable_schema = ON` + a textual REPLACE on `sqlite_master.sql`.
Faster and data-preserving compared to a copy-rebuild. Safe because we
only swap the dead `_migrate_old_<x>` table name back to its real name
`<x>` — the FK graph is otherwise unchanged. A `VACUUM` forces SQLite to
re-read the schema cache.
"""
function _fix_fk_references_after_autoincrement_migration!(db::SQLite.DB)
    broken = sqlrows(DBInterface.execute(db,
        "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE '%_migrate_old_%'"))
    isempty(broken) && return

    DBInterface.execute(db, "PRAGMA writable_schema = ON")
    try
        # Wrap the rewrite + drift assertion in an explicit transaction so a
        # partial failure (e.g. SQLITE_CORRUPT, disk full) OR a future entity
        # drift (assertion fires) rolls back atomically rather than leaving
        # the DB half-healed (review issues #20 + #25). The assertion MUST
        # run inside the tx — it's the only drift detector for issue #19.
        # `writable_schema` is connection-local, so toggling it outside the
        # transaction is safe.
        SQLite.transaction(db) do
            # Replace each `_migrate_old_<entity>` with `<entity>`. Driven from
            # `_MIGRATION_TEMP_ENTITIES` (single source of truth — covers
            # both the AUTOINCREMENT rename and R2 widen rename patterns)
            # so adding a new entity to either migration automatically
            # heals here too (review issues #19 + #22).
            #
            # Quoted form ("_migrate_old_x") and bare form (_migrate_old_x)
            # both occur in CREATE statements depending on whether the FK
            # target was rendered with or without quoting.
            #
            # ORDERING INVARIANT (review suggestion #17): the unanchored
            # substring REPLACE on `_migrate_old_<entity>` is safe today
            # because no entity name in `_MIGRATION_TEMP_ENTITIES` is a
            # substring of another (e.g. no `peak` AND `peaks`). If you add
            # an entity here, verify this invariant or anchor the replace.
            for entity in _MIGRATION_TEMP_ENTITIES
                old_q = "\"_migrate_old_$entity\""
                old_b = "_migrate_old_$entity"
                DBInterface.execute(db,
                    "UPDATE sqlite_master SET sql = REPLACE(sql, ?, ?) WHERE type='table'",
                    [old_q, entity])
                DBInterface.execute(db,
                    "UPDATE sqlite_master SET sql = REPLACE(sql, ?, ?) WHERE type='table'",
                    [old_b, entity])
            end
            # Defense-in-depth: assert the heal worked. Catches a future drift
            # where an entity gets added to migrate_pk_to_autoincrement! /
            # migrate_r2_widen_index_peaks_pk! but the constant list above
            # isn't updated (review issue #19). Inside the tx so the throw
            # rolls back the partial REPLACE pass (review issue #25).
            remaining = sqlrows(DBInterface.execute(db,
                "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE '%_migrate_old_%'"))
            isempty(remaining) || error(
                "_fix_fk_references_after_autoincrement_migration!: tables still " *
                "carry `_migrate_old_*` FKs after heal: " *
                join(String[String(r.name) for r in remaining], ", ") *
                " (likely missing entry in _MIGRATION_TEMP_ENTITIES " *
                "or _AUTOINCREMENT_ENTITIES — review suggestion #18)")
        end
    finally
        DBInterface.execute(db, "PRAGMA writable_schema = OFF")
    end
    # Force SQLite to invalidate its in-memory schema cache so subsequent
    # prepares see the corrected FKs. VACUUM is the standard idiom but it
    # acquires an EXCLUSIVE lock and can fail with SQLITE_BUSY under
    # multi-process access. Fall back to bumping schema_version (also
    # invalidates the cache) so the heal still takes effect (review issue #21).
    #
    # The schema_version PRAGMA is writable only when writable_schema=ON on
    # most builds (review issue #26). The `finally` above already toggled it
    # off, so the fallback re-enables it for the bump and toggles off again.
    try
        DBInterface.execute(db, "VACUUM")
    catch err
        @warn "VACUUM after FK heal failed; falling back to schema_version bump" exception=err
        DBInterface.execute(db, "PRAGMA writable_schema = ON")
        try
            DBInterface.execute(db, "PRAGMA schema_version = schema_version + 1")
        finally
            DBInterface.execute(db, "PRAGMA writable_schema = OFF")
        end
    end
end

"""
    migrate_samples_naming!(db) :: Nothing

Convert legacy `samples (label, name)` shape to `(name, display_name)`. Idempotent
on the canonical shape (sentinel: `display_name` exists AND `label` does not).
Atomic: all four ALTER/UPDATE/DROP statements + the duplicate-suffix pass +
the UNIQUE INDEX creation + the idempotent_responses purge run inside one
SQLite.transaction so a partial migration is impossible.

Pre-existing `(experiment_id, name)` duplicates (the missing UNIQUE never
blocked them) are renamed to `<name>-2`, `<name>-3`, … ordered by ascending id
(oldest sample wins the bare name; deterministic across reruns).
"""
function migrate_samples_naming!(db::SQLite.DB)::Nothing
    cols = Set(r.name for r in sqlrows(DBInterface.execute(db,
        "PRAGMA table_info('samples')")))
    # No samples table yet (partial legacy fixture or pre-create_schema! call) — skip.
    # create_schema! will create it with the canonical shape; nothing to rename.
    isempty(cols) && return nothing
    # Post-redesign shape: a single `name` label, no `display_name`/`label`. Nothing to rename.
    # Without this, the name-collapse migration (migrate_samples_name_collapse!) is silently
    # reverted on every open_db. (Phase-A P0-1.)
    if "name" in cols && !("display_name" in cols) && !("label" in cols)
        return nothing
    end
    if "display_name" in cols && !("label" in cols)
        return nothing  # already migrated
    end
    SQLite.transaction(db) do
        if !("display_name" in cols)
            try DBInterface.execute(db, "ALTER TABLE samples ADD COLUMN display_name TEXT")
            catch e; occursin("duplicate column", sprint(showerror, e)) || rethrow(); end
        end
        if "label" in cols
            if "name" in cols
                # Full Plan-7 shape: both label (stable id) and name (friendly) present.
                # label → name (stable id), name → display_name (friendly text).
                DBInterface.execute(db, """UPDATE samples
                                              SET display_name = name,
                                                  name         = COALESCE(NULLIF(label, ''), name)""")
            else
                # Older shape: only label present (no name column yet).
                # Add name column and populate from label; display_name stays NULL.
                try DBInterface.execute(db, "ALTER TABLE samples ADD COLUMN name TEXT")
                catch e; occursin("duplicate column", sprint(showerror, e)) || rethrow(); end
                DBInterface.execute(db, "UPDATE samples SET name = COALESCE(NULLIF(label, ''), CAST(id AS TEXT))")
            end
            try DBInterface.execute(db, "ALTER TABLE samples DROP COLUMN label")
            catch e; occursin("no such column", sprint(showerror, e)) || rethrow(); end
        end
        # Duplicate-label disambiguation retired (ingestion redesign): sample labels may
        # legitimately repeat across loads (e.g. two `HA85 (S01P15)`); identity is (load_id, slot_index),
        # not the label. See migrate_samples_name_collapse! and spec §10.
        DBInterface.execute(db,
            "CREATE UNIQUE INDEX IF NOT EXISTS samples_unique_name ON samples(experiment_id, name)")
        # Old idempotent_responses rows carry pre-rename payload shape; purge to
        # prevent stale-shape replays on retried client_op_id keys post-deploy.
        try DBInterface.execute(db, "DELETE FROM idempotent_responses")
        catch e; occursin("no such table", lowercase(sprint(showerror, e))) || rethrow(); end
    end
    nothing
end

"""
    migrate_exposures_unique_filename!(db)

Add `UNIQUE INDEX exposures_unique_filename ON exposures(sample_id, filename)`
following the dedupe-then-enforce pattern from `migrate_samples_naming!`.
Renames pre-existing duplicates (oldest id keeps the bare filename;
second-and-later get `<filename>-2`, …) before creating the index, so the
`CREATE UNIQUE INDEX` always succeeds against clean data.

Idempotent on re-run. Wrapped in `SQLite.transaction` so a partial run
cannot leave duplicates renamed without uniqueness enforcement.

Direct-invocation pattern from CLAUDE.md (FK-heal regression tests) —
tests can call this without going through `open_db`.
"""
function migrate_exposures_unique_filename!(db::SQLite.DB)
    SQLite.transaction(db) do
        # Track existing (sample_id, filename) pairs so a row literally named
        # "<x>-2" doesn't collide with our rename target.
        existing = Set{Tuple{Int64,String}}(
            (Int64(r.sample_id), String(r.filename))
            for r in sqlrows(DBInterface.execute(db,
                "SELECT sample_id, filename FROM exposures WHERE sample_id IS NOT NULL AND filename IS NOT NULL")))

        dups = sqlrows(DBInterface.execute(db, """
            SELECT sample_id, filename FROM exposures
            WHERE sample_id IS NOT NULL AND filename IS NOT NULL
            GROUP BY sample_id, filename HAVING COUNT(*) > 1"""))

        for d in dups
            ids = sqlrows(DBInterface.execute(db,
                "SELECT id FROM exposures WHERE sample_id = ? AND filename = ? ORDER BY id ASC",
                [d.sample_id, d.filename]))
            for (i, row) in enumerate(ids)
                i == 1 && continue  # oldest keeps the bare name
                suffix_n = i
                new_name = "$(d.filename)-$(suffix_n)"
                while (Int64(d.sample_id), new_name) in existing
                    suffix_n += 1
                    new_name = "$(d.filename)-$(suffix_n)"
                end
                push!(existing, (Int64(d.sample_id), new_name))
                @warn "Renamed duplicate exposure" sample_id=d.sample_id old=d.filename new=new_name id=row.id
                DBInterface.execute(db,
                    "UPDATE exposures SET filename = ? WHERE id = ?",
                    [new_name, row.id])
            end
        end

        DBInterface.execute(db,
            "CREATE UNIQUE INDEX IF NOT EXISTS exposures_unique_filename ON exposures(sample_id, filename)")
    end
    nothing
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
    info = sqlrows(DBInterface.execute(db,
        "SELECT name, pk FROM pragma_table_info('index_peaks')"))
    already_widened = any(r -> String(r.name) == "peak_kind" && Int(r.pk) > 0, info)
    already_widened && return
    # Skip on fresh DBs that already have the new shape via CREATE.
    isempty(info) && return

    DBInterface.execute(db, "PRAGMA foreign_keys = OFF")
    try
        SQLite.transaction(db) do
            DBInterface.execute(db, "ALTER TABLE index_peaks RENAME TO _migrate_old_index_peaks")
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
                FROM _migrate_old_index_peaks
            """)
            DBInterface.execute(db, "DROP TABLE _migrate_old_index_peaks")
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
    peaks_exists = !isempty(sqlrows(DBInterface.execute(db,
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
    peaks_ever_written = !isempty(sqlrows(DBInterface.execute(db,
        "SELECT 1 FROM sqlite_sequence WHERE name = 'peaks'"))) ||
        (first(sqlrows(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM peaks"))).n > 0)
    peaks_ever_written || return

    # Sentinel: if auto_peaks already has rows, we're partway through —
    # the only safe action is to bail and require operator intervention.
    auto_peaks_exists = !isempty(sqlrows(DBInterface.execute(db,
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='auto_peaks'")))
    if auto_peaks_exists
        auto_count = first(sqlrows(DBInterface.execute(db,
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
                    findpeaks_index INTEGER,
                    found_by        TEXT NOT NULL DEFAULT 'hybrid' CHECK (found_by IN ('hybrid', 'model', 'both')),
                    model_id        INTEGER
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
                sqlrows(DBInterface.execute(db, "PRAGMA table_info(peaks)")) .|>
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
            manual_rows = sqlrows(DBInterface.execute(db,
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
            orphans = sqlrows(DBInterface.execute(db, """
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
        name::Union{AbstractString,Nothing} = nothing,
        path::AbstractString,
        data_dir::Union{AbstractString,Nothing} = nothing,
        analysis_dir::Union{AbstractString,Nothing} = nothing,
        description = nothing,
        last_scanned_at = nothing, ingest_status = "idle",
        setup...)
    result = DBInterface.execute(db, """
        INSERT INTO experiments (name, path, description, last_scanned_at, ingest_status)
        VALUES (?,?,?,?,?)
    """, [name, String(path), description, last_scanned_at, ingest_status])
    exp_id = Int(DBInterface.lastrowid(result))
    if data_dir !== nothing
        create_setup!(db; experiment_id = exp_id, data_dir = data_dir,
            analysis_dir = something(analysis_dir, data_dir), setup...)
    end
    return exp_id
end

function create_load!(db::SQLite.DB;
        experiment_id::Integer,
        load_index::Integer,
        session_id::Union{Integer, Nothing}    = nothing,
        start_time::Union{String, Nothing}     = nothing,
        end_time::Union{String, Nothing}       = nothing,
        frame_count::Integer                   = 0,
        note::Union{String, Nothing}           = nothing,
        setup_id                               = default_setup_id(db, experiment_id))
    res = DBInterface.execute(db, """
        INSERT INTO loads (experiment_id, setup_id, load_index, session_id, start_time, end_time, frame_count, note)
        VALUES (?,?,?,?,?,?,?,?)
    """, [experiment_id, setup_id, load_index, session_id, start_time, end_time, frame_count, note])
    Int(DBInterface.lastrowid(res))
end

function create_sample!(db::SQLite.DB; experiment_id::Integer, name::AbstractString,
        notes=nothing, load_id=nothing, slot_index=nothing,
        grouping_source="auto_position", name_source="auto",
        setup_id = default_setup_id(db, experiment_id))
    res = DBInterface.execute(db, """
        INSERT INTO samples (experiment_id, setup_id, name, notes, load_id, slot_index, grouping_source, name_source)
        VALUES (?,?,?,?,?,?,?,?)
    """, [experiment_id, setup_id, name, notes, load_id, slot_index, grouping_source, name_source])
    return Int(DBInterface.lastrowid(res))
end

function create_exposure!(db::SQLite.DB;
        experiment_id::Int,                              # required (Phase A)
        sample_id::Union{Int,Nothing}      = nothing,    # optional (transient pre-group state)
        filename::Union{String,Nothing}    = nothing,
        kind::String                       = "file",
        status::Union{String,Nothing}      = nothing,
        image_path::Union{String,Nothing}  = nothing,
        # --- new PRP fields (Phase A) ---
        prp_path          = nothing,
        timestamp         = nothing,
        exposure_time     = nothing,
        horizontal_position = nothing,
        scan_id           = nothing,
        frame_no          = nothing,
        load_id           = nothing,
        content_fingerprint = nothing,
        transmitted_counts = nothing,
        incident_monitor  = nothing,
        counting_time     = nothing,
        pipe_length_m     = nothing,
        monitors_read::Bool = prp_path !== nothing,
        setup_id = default_setup_id(db, experiment_id))
    result = DBInterface.execute(db, """
        INSERT INTO exposures
            (setup_id, experiment_id, sample_id, filename, kind, status, image_path,
             prp_path, timestamp, exposure_time, horizontal_position, scan_id, frame_no,
             load_id, content_fingerprint, transmitted_counts, incident_monitor, counting_time,
             pipe_length_m, monitors_read)
        VALUES (?,?,?,?,?,?,?, ?,?,?,?,?,?,?,?, ?,?,?,?,?)
    """, [setup_id, experiment_id, sample_id, filename, kind, status, image_path,
          prp_path, timestamp, exposure_time, horizontal_position, scan_id, frame_no,
          load_id, content_fingerprint, transmitted_counts, incident_monitor, counting_time,
          pipe_length_m, Int(monitors_read)])
    return Int(DBInterface.lastrowid(result))
end

function get_experiment(db::SQLite.DB, id::Int)
    rows = sqlrows(DBInterface.execute(db,
        "SELECT * FROM experiments WHERE id = ?", [id]))
    isempty(rows) && error("experiment $id not found")
    first(rows)
end

function get_samples(db::SQLite.DB, experiment_id::Int)
    sqlrows(DBInterface.execute(db,
        "SELECT * FROM samples WHERE experiment_id = ? ORDER BY id", [experiment_id]))
end

function get_exposures(db::SQLite.DB, sample_id::Int)
    sqlrows(DBInterface.execute(db,
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
    DBInterface.execute(db, "PRAGMA busy_timeout = 5000")

    # WAL lets concurrent readers proceed alongside one writer — load-bearing
    # for parallel request handling (#115). The default rollback journal
    # serializes every reader behind any in-flight writer. WAL persists in
    # the DB file header, so this PRAGMA is effectively a one-time migration
    # for existing DBs and a no-op on every subsequent open. Skipped for
    # `:memory:` and shared-cache URI DBs, which can't run WAL.
    #
    # `sqlrows(...)` drains the result iterator so the prepared
    # statement is dropped before `finalize_statements!`. Without the drain,
    # the iterator keeps the PRAGMA's statement attached, and any subsequent
    # DDL (e.g. `DROP TABLE` in tests, or migration loops on legacy DBs)
    # fails with `database table is locked`. The trailing
    # `finalize_statements!` then clears the cache so callers see a
    # quiescent connection.
    if db_path != ":memory:" && !startswith(db_path, "file:")
        sqlrows(DBInterface.execute(db, "PRAGMA journal_mode = WAL"))
        SQLite.finalize_statements!(db)
    end

    # SQLite hardcodes O_CREAT mode 0644 in os_unix.c — process umask only
    # masks bits OUT, so umask 0002 can't promote 0644 to 0664. For
    # multi-user deploys (curators in a shared group writing the same DB),
    # we need group-write on the file. WAL creates `-wal` and `-shm`
    # sidecars on first write; chmod those too so other group members can
    # write through them.
    #
    # Skip the chmod when another user owns the file: in the shared-group
    # deploy, only the owner can chmod (chmod fails — typically EPERM —
    # on a cross-user file), and they presumably already set group-write
    # on first open. On a file WE own, any failure (read-only mount,
    # immutable bit, EROFS, fs quirk) is a real configuration problem
    # and must propagate — the prior broad `e isa IOError || e isa
    # SystemError` catch silently masked all of those, leaving SQLite's
    # hardcoded 0644 in place.
    # Base.Libc.getuid() wraps :jl_getuid (portable across platforms) and
    # returns Culong, the same type as stat(p).uid — no promotion gymnastics
    # at the comparison.
    my_uid = Base.Libc.getuid()
    for p in (db_path, db_path * "-wal", db_path * "-shm")
        isfile(p) || continue
        if Sys.isunix() && stat(p).uid != my_uid
            continue  # cross-user shared-group case: owner sets perms
        end
        chmod(p, 0o664)
    end
    db
end
