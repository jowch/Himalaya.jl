using SQLite, DBInterface, Tables
using Dates: now, UTC, format, @dateformat_str

# Sentinel marker name for the I3.1 comparison→series data migration (#171).
const MIGRATION_COMPARISONS_TO_SERIES = "comparisons_to_series"

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
    name          TEXT,
    display_name  TEXT,
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
    id                   INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id            INTEGER REFERENCES samples(id),
    filename             TEXT,
    kind                 TEXT DEFAULT 'file',
    selected             BOOLEAN DEFAULT FALSE,
    status               TEXT CHECK (status IN ('accepted', 'rejected')),
    image_path           TEXT,
    trace_hash           TEXT,
    analysis_inputs_hash TEXT
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
    kind        TEXT NOT NULL DEFAULT 'auto',
    inputs_hash TEXT
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
    # Each ALTER TABLE adds a column that may already exist on legacy DBs;
    # SQLite has no `ADD COLUMN IF NOT EXISTS`, so we tolerate "duplicate column"
    # specifically. CREATE INDEX/TABLE statements below are `IF NOT EXISTS`-
    # guarded and any unexpected failure must surface — never swallow them
    # under a bare catch (issue #6 from PR review).
    alter_stmts = [
        "ALTER TABLE exposures ADD COLUMN status TEXT CHECK (status IN ('accepted', 'rejected'))",
        "ALTER TABLE exposures ADD COLUMN image_path TEXT",
        "ALTER TABLE experiments ADD COLUMN config TEXT",
        "ALTER TABLE experiments ADD COLUMN experiment_type TEXT",
        "ALTER TABLE experiments ADD COLUMN energy_kev REAL",
        "ALTER TABLE experiments ADD COLUMN flight_path_m REAL",
        "ALTER TABLE users ADD COLUMN first_name TEXT",
        "ALTER TABLE users ADD COLUMN last_name TEXT",
        "ALTER TABLE indices ADD COLUMN kind TEXT NOT NULL DEFAULT 'auto'",
        "ALTER TABLE exposures ADD COLUMN trace_hash TEXT",
        "ALTER TABLE exposures ADD COLUMN analysis_inputs_hash TEXT",
        "ALTER TABLE indices ADD COLUMN inputs_hash TEXT",
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
    # like sample_messages/sample_tags/etc. Idempotent — no-op when clean.
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
    legacy = Tables.rowtable(DBInterface.execute(db,
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

    # Compare page (Plan §Phase 1, Task 1.1): comparisons / comparison_members /
    # comparison_messages. Must run AFTER R4 — none of the compare tables touch
    # user_actions, but ordering keeps every cross-cutting fix-up bounded by the
    # earlier R-numbered migrations. See docs/superpowers/specs/2026-05-02-compare-page-design.md.
    migrate_compare!(db)

    # Compare page Phase 13: per-user pinned comparisons.
    migrate_comparison_pins!(db)

    # Issue #67: drop NOT NULL on the four `comparisons` columns the route
    # used to seed with empty-string placeholders. Must run AFTER
    # migrate_comparison_pins! (and migrate_compare!) — both create FK refs
    # to comparisons that the relax migration heals when it RENAMEs the
    # table during the rebuild.
    migrate_compare_relax_nullability!(db)

    # 2026-05-14 — Compare UX refinement (spec §6.4): persist view choices
    # on the comparison so they round-trip across viewers.
    migrate_compare_view_choices!(db)

    # Series model (#164, master plan §5.1): new `series*` tables + the
    # `schema_migrations` sentinel. Placed after the compare migrations so
    # the compare/series schema stays grouped and #171's
    # `migrate_comparisons_to_series!` has a natural slot after this. New
    # tables only — the `comparison*` tables are never renamed (§2.1).
    migrate_series!(db)

    # I3.1 (#171): copy the comparison corpus into the series* tables. MUST run
    # after migrate_series! (tables + schema_migrations sentinel exist) and after
    # the two compare migrations above (reads their columns). Own transaction;
    # sentinel-gated; raw-INSERT user_actions, never apply_event! (no broadcast).
    migrate_comparisons_to_series!(db)

    # PR #107 left the on-disk experiment.toml AND the in-DB experiments.config
    # blob using the legacy `[manifest].label/name` shape. The deprecation
    # error in `_build_config` (config.jl) hard-fails any route that calls
    # `config_from_db` (trace plot, analyze_exposure!, reanalyze) for those
    # experiments. Migrate the in-DB blob in place; the on-disk file is the
    # operator's responsibility (`himalaya migrate-toml <dir>`) but is no
    # longer load-bearing at runtime since `experiments.config` is the
    # source of truth for `analyze_exposure!`.
    migrate_experiment_config_label_to_name!(db)
end

"""
    migrate_experiment_config_label_to_name!(db)

Rewrite each `experiments.config` blob from the legacy `[manifest].label/name`
shape to the canonical `[manifest].name/display_name` shape (PR #107). Pure
text rewrite via `migrate_manifest_toml_text`; only touches blobs that
actually contain `[manifest].label`. Idempotent — already-migrated blobs
return unchanged and we skip the UPDATE.

If a blob is in the corrupt "both `label` AND `display_name`" state, the
helper raises; we let it propagate so the operator can see which experiment
is broken rather than silently masking it.
"""
function migrate_experiment_config_label_to_name!(db::SQLite.DB)
    rows = try
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, config FROM experiments WHERE config IS NOT NULL"))
    catch err
        # In normal `migrate_schema!` flow the `config` column already
        # exists — `alter_stmts` adds it earlier in the same function. This
        # tolerance is purely for test fixtures that build a `SQLite.DB`
        # directly and call this migration in isolation (no
        # `migrate_schema!` invocation). Production `open_db` callers never
        # hit this branch.
        msg = lowercase(sprint(showerror, err))
        (occursin("no such column", msg) || occursin("no such table", msg)) ||
            rethrow()
        return nothing
    end
    for row in rows
        blob = String(row.config)
        new_text, changed = migrate_manifest_toml_text(blob)
        changed || continue
        DBInterface.execute(db,
            "UPDATE experiments SET config = ? WHERE id = ?",
            [new_text, row.id])
        @info "Healed experiments.config blob (legacy [manifest].label → name/display_name)" experiment_id = Int(row.id)
    end
    nothing
end

"""
    migrate_compare!(db)

Install the Compare-page tables (`comparisons`, `comparison_members`,
`comparison_messages`) and their supporting indexes. Idempotent — every
statement is `IF NOT EXISTS`-guarded so reopening an already-migrated DB
is a no-op.

Why each FK action:
- `comparisons.forked_from_id ON DELETE SET NULL` — forks survive parent
  deletion as independent artifacts; the `forked_at_hash` (immutable) is
  retained for historical reference.
- `comparisons.created_by ON DELETE SET NULL` — user-FK rule.
- `comparison_members.comparison_id ON DELETE CASCADE` — members are part
  of the artifact; deleting the comparison drops its membership.
- `comparison_members.exposure_id ON DELETE SET NULL` — exposure deletion
  leaves a visible orphan placeholder rather than silently mutating the
  figure (preserves chat references).
- `comparison_messages.comparison_id ON DELETE CASCADE` — chat is part
  of the comparison's discussion thread.
- `comparison_messages.author_id ON DELETE SET NULL` — user-FK rule.

`comparisons.id` uses `INTEGER PRIMARY KEY AUTOINCREMENT` because comparisons
are `@`-mention targets (mention-target rule, see CLAUDE.md). Members and
messages use plain `INTEGER PRIMARY KEY` — neither is `@`-mentioned, and
`comparison_messages` matches the existing `sample_messages` shape.
"""
function migrate_compare!(db::SQLite.DB)
    # `title`, `content_hash`, `created_at`, `updated_at` are nullable here
    # (issue #67): the route uses `INSERT INTO comparisons DEFAULT VALUES`
    # to mint the AUTOINCREMENT id and the dispatcher fills these via
    # `COALESCE(col, ?)`. Pre-#67 these were `TEXT NOT NULL` and the route
    # seeded `''` placeholders, which trapped #54 because empty-string is
    # NOT NULL but IS a sentinel that the dispatcher's COALESCE preserved.
    # `migrate_compare_relax_nullability!` rebuilds legacy DBs to this shape.
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS comparisons (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            title           TEXT,
            description     TEXT,
            content_hash    TEXT,
            created_by      INTEGER REFERENCES users(id) ON DELETE SET NULL,
            created_at      TEXT,
            updated_at      TEXT,
            forked_from_id  INTEGER REFERENCES comparisons(id) ON DELETE SET NULL,
            forked_at_hash  TEXT
        )""")
    DBInterface.execute(db, """
        CREATE INDEX IF NOT EXISTS idx_comparisons_forked_from
            ON comparisons(forked_from_id)""")

    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS comparison_members (
            id              INTEGER PRIMARY KEY,
            comparison_id   INTEGER NOT NULL REFERENCES comparisons(id) ON DELETE CASCADE,
            exposure_id     INTEGER REFERENCES exposures(id) ON DELETE SET NULL,
            display_order   INTEGER NOT NULL,
            band_height     REAL    NOT NULL DEFAULT 1.0,
            y_offset        REAL    NOT NULL DEFAULT 0,
            normalization   TEXT    NOT NULL DEFAULT 'none',
            color_override  TEXT,
            label_override  TEXT,
            q_window_min    REAL,
            q_window_max    REAL,
            peak_display    TEXT    CHECK (peak_display IS NULL OR json_valid(peak_display)),
            snapshot        TEXT    NOT NULL CHECK (json_valid(snapshot)),
            created_by      INTEGER REFERENCES users(id) ON DELETE SET NULL,
            created_at      TEXT    NOT NULL
        )""")
    DBInterface.execute(db, """
        CREATE INDEX IF NOT EXISTS idx_comparison_members_by_comparison
            ON comparison_members(comparison_id, display_order)""")

    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS comparison_messages (
            id            INTEGER PRIMARY KEY,
            comparison_id INTEGER NOT NULL REFERENCES comparisons(id) ON DELETE CASCADE,
            author_id     INTEGER REFERENCES users(id) ON DELETE SET NULL,
            body          TEXT NOT NULL,
            created_at    DATETIME DEFAULT CURRENT_TIMESTAMP
        )""")
    DBInterface.execute(db, """
        CREATE INDEX IF NOT EXISTS idx_comparison_messages_comparison
            ON comparison_messages(comparison_id, created_at)""")
    nothing
end

"""
    migrate_comparison_pins!(db)

Install the `comparison_pins` table for per-user pinned comparisons
(Plan §Phase 13, Task 13.2). Idempotent.

Composite PK on `(user_id, comparison_id)` enforces "one pin per (user,
comparison)" — pinning the same comparison twice is a no-op INSERT OR
IGNORE. Both FK columns cascade on delete: if a user is removed, their
pins disappear; if a comparison is deleted, the pin disappears with it.
This matches the expectation that a "pin" is purely metadata about an
otherwise-unaffected user/comparison pair.

Pin/unpin routes wrap in `with_idempotency` (the wrapper provides the
outer SQLite tx that `apply_event!(InTransaction(), …)` requires); when
no `X-Client-Op-Id` header is present the wrapper falls through, and
repeated POSTs simply re-affirm the "pinned" state via
`INSERT OR IGNORE` / idempotent `DELETE`. View-table state stays
correct either way. The `pinned_at` timestamp captures user-perceived
ordering for display in the sidebar (most-recently-pinned first).
"""
function migrate_comparison_pins!(db::SQLite.DB)
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS comparison_pins (
            user_id        INTEGER NOT NULL REFERENCES users(id)       ON DELETE CASCADE,
            comparison_id  INTEGER NOT NULL REFERENCES comparisons(id) ON DELETE CASCADE,
            pinned_at      TEXT    NOT NULL DEFAULT CURRENT_TIMESTAMP,
            PRIMARY KEY (user_id, comparison_id)
        )""")
    DBInterface.execute(db, """
        CREATE INDEX IF NOT EXISTS idx_comparison_pins_by_user
            ON comparison_pins(user_id, pinned_at DESC)""")
    nothing
end

"""
    migrate_compare_relax_nullability!(db)

Issue #67: drop `NOT NULL` on `title`, `content_hash`, `created_at`,
`updated_at` of the `comparisons` table. Pre-#67 these four columns were
`TEXT NOT NULL` and the route's mint-the-id INSERT seeded them with
empty-string placeholders that the dispatcher then folded over. Empty
string is NOT NULL but IS a sentinel — it trapped #54 because plain
`COALESCE(created_at, ?)` preserved `''` forever. PR #66 patched the
symptom with `NULLIF(created_at, '')`; this migration removes the trap
structurally so the route can use a clean `INSERT INTO comparisons
DEFAULT VALUES` and the dispatcher's `COALESCE` Just Works on NULLs.

Idempotent — sentinel parses `sqlite_master.sql` and no-ops once the
relaxed shape is already installed (fresh DBs created by the post-#67
`migrate_compare!` are already relaxed).

Data heal: rows with `created_at = ''` (the #54 bug) have it rewritten
to `CURRENT_TIMESTAMP` so the frontend's row-shape contract (non-empty
ISO timestamp, asserted in test_route_response_shapes.jl) holds for
legacy rows too.
"""
function migrate_compare_relax_nullability!(db::SQLite.DB)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT sql FROM sqlite_master WHERE type='table' AND name='comparisons'"))
    isempty(rows) && return  # comparisons table absent (pre-migrate_compare!)
    sql = String(rows[1].sql)
    # Sentinel: `title TEXT NOT NULL` is the canonical legacy marker; absent
    # iff the relaxed schema is already installed.
    occursin(r"\btitle\s+TEXT\s+NOT\s+NULL\b"i, sql) || return

    # FK enforcement must be disabled OUTSIDE the transaction (SQLite docs)
    # — comparisons is referenced by comparison_members (CASCADE) and
    # comparison_pins (CASCADE), which would otherwise trip during the
    # rename + rebuild sequence.
    DBInterface.execute(db, "PRAGMA foreign_keys = OFF")
    try
        SQLite.transaction(db) do
            DBInterface.execute(db,
                "ALTER TABLE comparisons RENAME TO _migrate_old_comparisons")
            DBInterface.execute(db, """
                CREATE TABLE comparisons (
                    id              INTEGER PRIMARY KEY AUTOINCREMENT,
                    title           TEXT,
                    description     TEXT,
                    content_hash    TEXT,
                    created_by      INTEGER REFERENCES users(id) ON DELETE SET NULL,
                    created_at      TEXT,
                    updated_at      TEXT,
                    forked_from_id  INTEGER REFERENCES comparisons(id) ON DELETE SET NULL,
                    forked_at_hash  TEXT
                )""")
            DBInterface.execute(db, """
                INSERT INTO comparisons
                  (id, title, description, content_hash, created_by,
                   created_at, updated_at, forked_from_id, forked_at_hash)
                SELECT id, title, description, content_hash, created_by,
                       created_at, updated_at, forked_from_id, forked_at_hash
                  FROM _migrate_old_comparisons""")
            DBInterface.execute(db, "DROP TABLE _migrate_old_comparisons")
            # The forked_from_id index travels with the table — recreate it.
            DBInterface.execute(db, """
                CREATE INDEX IF NOT EXISTS idx_comparisons_forked_from
                    ON comparisons(forked_from_id)""")
        end
        # After the transaction commits, heal FK references corrupted by
        # SQLite's ALTER TABLE RENAME tracking — comparison_members and
        # comparison_pins (and any future table that REFERENCES comparisons)
        # had their stored FK target rewritten to `_migrate_old_comparisons`.
        _heal_renamed_table_fk_refs!(db, "comparisons")
        # One-time data heal: convert the #54 stale `''` created_at /
        # updated_at rows to a real timestamp so the frontend's row-shape
        # assertion holds for legacy rows. Format MUST match
        # `comparison_now_iso()` — fresh rows use
        # `yyyy-mm-ddTHH:MM:SS.sssZ`, and ComparisonSidebar.tsx sorts on
        # `updated_at` as a string. SQLite's `CURRENT_TIMESTAMP` returns
        # `YYYY-MM-DD HH:MM:SS` (space, no `Z`); space sorts BEFORE `T`,
        # so healed legacy rows would lexically precede fresh rows of the
        # same instant. Substitute a Julia-formatted string instead.
        # Idempotent — no-op once the rewrite has run.
        now_iso = format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ")
        DBInterface.execute(db, """
            UPDATE comparisons
               SET created_at = ?
             WHERE created_at = ''""", [now_iso])
        DBInterface.execute(db, """
            UPDATE comparisons
               SET updated_at = ?
             WHERE updated_at = ''""", [now_iso])
    finally
        DBInterface.execute(db, "PRAGMA foreign_keys = ON")
    end
end

"""
    migrate_compare_view_choices!(db)

Add `view_grouping_mode`, `view_show_peak_ticks`, `view_show_peak_labels`
columns to `comparisons` so the author's view choices round-trip across
viewers (spec §6.4). All NULL on existing rows; the frontend falls back to
per-tab Zustand defaults when NULL.

Idempotent: each `ALTER TABLE ... ADD COLUMN` is wrapped in a try/catch
that treats "duplicate column name" as success — the same pattern used
by other additive migrations.
"""
function migrate_compare_view_choices!(db::SQLite.DB)
    for stmt in (
        "ALTER TABLE comparisons ADD COLUMN view_grouping_mode TEXT",
        "ALTER TABLE comparisons ADD COLUMN view_show_peak_ticks INTEGER",
        "ALTER TABLE comparisons ADD COLUMN view_show_peak_labels INTEGER",
    )
        try
            DBInterface.execute(db, stmt)
        catch err
            msg = sprint(showerror, err)
            occursin("duplicate column name", lowercase(msg)) || rethrow()
        end
    end
end

"""
    migrate_series!(db)

Install the Series-model tables (`series`, `series_members`, `series_samples`,
`series_messages`, `series_pins`) and the `schema_migrations` sentinel table.
Idempotent — every statement is `IF NOT EXISTS`-guarded, so reopening an
already-migrated DB is a no-op.

The `series` model is a *new* table set, never a rename of `comparison*`:
`rebuild_views_from_log!` re-folds historical events through dispatcher
branches that `INSERT INTO` the named `comparison*` tables, so renaming them
would break event-log replay on every existing DB (master plan §2.1). The
`comparison*` tables therefore stay permanently as replay machinery; series
data is copied forward later by `migrate_comparisons_to_series!` (#171).

Each table is created in its final shape in one `CREATE TABLE` — unlike the
`comparison*` tables, which reached their shape through three migrations.

Why each FK action:
- `series.created_by ON DELETE SET NULL` / `series_members.created_by` /
  `series_messages.author_id` — user-FK rule.
- `series.forked_from_id ON DELETE SET NULL` — forks survive parent deletion
  as independent artifacts (the `comparisons.forked_from_id` precedent).
- `series_members.exposure_id ON DELETE SET NULL` — exposure deletion leaves
  a visible orphan placeholder rather than silently mutating the figure
  (the `comparison_members.exposure_id` precedent).
- All four child tables `… series_id … ON DELETE CASCADE` — children are
  part of the artifact; deleting the series drops them, so the future
  `series_deleted` dispatcher branch stays a one-line `DELETE FROM series`.
- `series_samples.sample_id ON DELETE CASCADE` (not `SET NULL`) — a
  `series_samples` row is a pure pointer with no snapshot, so a NULL
  `sample_id` is unrenderable (master plan §5.1).

`series.id` is `INTEGER PRIMARY KEY AUTOINCREMENT` because series are
`@`-mention targets (mention-target rule, see CLAUDE.md), exactly as
`comparisons.id` is. The child tables use plain `INTEGER PRIMARY KEY` —
none is `@`-mentioned, and `series_samples.id` is replay-volatile (master
plan §11): never key client state on it.

`schema_migrations` is created empty. The `migrate_comparisons_to_series!`
copy (#171) writes its marker row last, inside its own transaction.

Timestamp caveat (inherited from `comparison_messages`, issue #76):
`series_messages.created_at` is `DATETIME DEFAULT CURRENT_TIMESTAMP`, which
yields the space-separated `YYYY-MM-DD HH:MM:SS` form — NOT the ISO
`yyyy-mm-ddTHH:MM:SS.sssZ` form that `series.created_at` carries when a
route writes it via `comparison_now_iso()`. A future series route or
dispatcher must not string-sort `series_messages.created_at` against the
`series`/`series_members` timestamps; sort `series_messages` on its own
column only.
"""
function migrate_series!(db::SQLite.DB)
    # `series`: the comparison columns (nullable, post-#67 shape) + the
    # view-choice columns + the recipe columns (`ordering_variable`,
    # `order_rule`, `state`). `content_hash` is NULL while `state='draft'`.
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS series (
            id                    INTEGER PRIMARY KEY AUTOINCREMENT,
            title                 TEXT,
            description           TEXT,
            content_hash          TEXT,
            created_by            INTEGER REFERENCES users(id)  ON DELETE SET NULL,
            created_at            TEXT,
            updated_at            TEXT,
            forked_from_id        INTEGER REFERENCES series(id) ON DELETE SET NULL,
            forked_at_hash        TEXT,
            view_grouping_mode    TEXT,
            view_show_peak_ticks  INTEGER,
            view_show_peak_labels INTEGER,
            ordering_variable     TEXT,
            order_rule            TEXT NOT NULL DEFAULT 'manual'
                                    CHECK (order_rule IN ('ascending','descending','manual')),
            state                 TEXT NOT NULL DEFAULT 'committed'
                                    CHECK (state IN ('draft','committed'))
        )""")
    DBInterface.execute(db, """
        CREATE INDEX IF NOT EXISTS idx_series_forked_from
            ON series(forked_from_id)""")

    # `series_members` — the plate. The `comparison_members` shape with
    # `comparison_id` renamed to `series_id` (CASCADE). `exposure_id` keeps
    # `ON DELETE SET NULL` (orphan-placeholder rule).
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS series_members (
            id              INTEGER PRIMARY KEY,
            series_id       INTEGER NOT NULL REFERENCES series(id)    ON DELETE CASCADE,
            exposure_id     INTEGER          REFERENCES exposures(id) ON DELETE SET NULL,
            display_order   INTEGER NOT NULL,
            band_height     REAL    NOT NULL DEFAULT 1.0,
            y_offset        REAL    NOT NULL DEFAULT 0,
            normalization   TEXT    NOT NULL DEFAULT 'none',
            color_override  TEXT,
            label_override  TEXT,
            q_window_min    REAL,
            q_window_max    REAL,
            peak_display    TEXT    CHECK (peak_display IS NULL OR json_valid(peak_display)),
            snapshot        TEXT    NOT NULL CHECK (json_valid(snapshot)),
            created_by      INTEGER REFERENCES users(id) ON DELETE SET NULL,
            created_at      TEXT    NOT NULL
        )""")
    DBInterface.execute(db, """
        CREATE INDEX IF NOT EXISTS idx_series_members_by_series
            ON series_members(series_id, display_order)""")

    # `series_samples` — the recipe membership: an explicit ordered sample
    # list. `series_id`/`sample_id` are both NOT NULL (a pointer row with a
    # NULL target is unrenderable; #164 spec text omits NOT NULL — see plan
    # preamble). `UNIQUE(series_id, position)` also serves ordered reads, so
    # no separate index is created.
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS series_samples (
            id          INTEGER PRIMARY KEY,
            series_id   INTEGER NOT NULL REFERENCES series(id)  ON DELETE CASCADE,
            sample_id   INTEGER NOT NULL REFERENCES samples(id) ON DELETE CASCADE,
            position    INTEGER NOT NULL,
            pinned      INTEGER NOT NULL DEFAULT 0 CHECK (pinned   IN (0,1)),
            excluded    INTEGER NOT NULL DEFAULT 0 CHECK (excluded IN (0,1)),
            UNIQUE(series_id, position)
        )""")

    # `series_messages` — the `comparison_messages` shape, `series_id` (CASCADE).
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS series_messages (
            id         INTEGER PRIMARY KEY,
            series_id  INTEGER NOT NULL REFERENCES series(id) ON DELETE CASCADE,
            author_id  INTEGER REFERENCES users(id) ON DELETE SET NULL,
            body       TEXT NOT NULL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP
        )""")
    DBInterface.execute(db, """
        CREATE INDEX IF NOT EXISTS idx_series_messages_by_series
            ON series_messages(series_id, created_at)""")

    # `series_pins` — the `comparison_pins` shape. Composite PK enforces one
    # pin per (user, series); both FKs CASCADE.
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS series_pins (
            user_id    INTEGER NOT NULL REFERENCES users(id)  ON DELETE CASCADE,
            series_id  INTEGER NOT NULL REFERENCES series(id) ON DELETE CASCADE,
            pinned_at  TEXT    NOT NULL DEFAULT CURRENT_TIMESTAMP,
            PRIMARY KEY (user_id, series_id)
        )""")
    DBInterface.execute(db, """
        CREATE INDEX IF NOT EXISTS idx_series_pins_by_user
            ON series_pins(user_id, pinned_at DESC)""")

    # Migration-version sentinel. No such table exists today; #171's copy
    # needs a real marker. Created empty here — #164 writes no rows.
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS schema_migrations (
            name        TEXT PRIMARY KEY,
            applied_at  TEXT
        )""")
    nothing
end

# Build the `series_created` recipe-snapshot payload for one comparison.
# Recipe samples are the DISTINCT, non-NULL sample_ids of the comparison's
# members' exposures, in first-seen display_order, at sequential 0-based
# positions. Orphan members (NULL exposure_id) and exposures with NULL
# sample_id contribute no recipe row. Master plan §6.1 step 2.
function _series_created_payload_from_comparison(db::SQLite.DB, cmp)
    sample_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT e.sample_id AS sample_id
           FROM comparison_members cm
           JOIN exposures e ON e.id = cm.exposure_id
           WHERE cm.comparison_id = ? AND e.sample_id IS NOT NULL
           ORDER BY cm.display_order ASC, cm.id ASC""", [Int(cmp.id)]))
    seen = Set{Int}()
    samples = Vector{Dict{Symbol,Any}}()
    for r in sample_rows
        sid = Int(r.sample_id)
        sid in seen && continue
        push!(seen, sid)
        push!(samples, Dict{Symbol,Any}(
            :sample_id => sid,
            :position  => length(samples),
            :pinned    => false,
            :excluded  => false,
        ))
    end
    Dict{Symbol,Any}(
        :title                 => ismissing(cmp.title) ? nothing : String(cmp.title),
        :description           => ismissing(cmp.description) ? nothing : String(cmp.description),
        # CRITICAL (P1 — id-space mismatch): NULL the fork lineage on migrated series.
        # `comparisons.forked_from_id` references comparisons(id); `series.forked_from_id`
        # references series(id) — DIFFERENT id-spaces. Copying the comparison's id straight
        # through would write a dangling reference (an arbitrary/nonexistent series id). FK
        # enforcement is ON at this point (migrate_compare_relax_nullability!'s `finally`
        # re-enabled `PRAGMA foreign_keys=ON` before this migration runs), so such a write
        # would FK-throw at INSERT and abort the migration; the `PRAGMA foreign_key_check`
        # test confirms the copy leaves no dangling ref. A comparison_id→series_id remap is
        # unsafe (id-order does NOT guarantee parent.id < child.id). Fork lineage across the
        # comparison→series cutover is not meaningful, so drop it: BOTH fields NULL.
        # (Human-approved option a.)
        :forked_from_id        => nothing,
        :forked_at_hash        => nothing,
        :ordering_variable     => nothing,   # comparisons have no recipe ordering
        :order_rule            => nothing,   # dispatcher COALESCEs to 'manual'
        :view_grouping_mode    => ismissing(cmp.view_grouping_mode) ? nothing : String(cmp.view_grouping_mode),
        :view_show_peak_ticks  => ismissing(cmp.view_show_peak_ticks) ? nothing : cmp.view_show_peak_ticks,
        :view_show_peak_labels => ismissing(cmp.view_show_peak_labels) ? nothing : cmp.view_show_peak_labels,
        :samples               => samples,
    )
end

# Build the `series_plate_committed` member-list payload for one comparison.
# Members carry NO ids (the dispatcher mints fresh PKs). Field shape matches
# `_series_member_payload` output. Master plan §6.1 step 2.
function _series_plate_committed_payload_from_comparison(db::SQLite.DB, cmp)
    member_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT exposure_id, display_order, band_height, y_offset, normalization,
                  color_override, label_override, q_window_min, q_window_max,
                  peak_display, snapshot
           FROM comparison_members WHERE comparison_id = ?
           ORDER BY display_order ASC, id ASC""", [Int(cmp.id)]))
    members = Vector{Dict{Symbol,Any}}()
    for m in member_rows
        push!(members, Dict{Symbol,Any}(
            :id             => nothing,
            :exposure_id    => ismissing(m.exposure_id)    ? nothing : Int(m.exposure_id),
            :display_order  => Int(m.display_order),
            :band_height    => Float64(m.band_height),
            :y_offset       => Float64(m.y_offset),
            :normalization  => String(m.normalization),
            :color_override => ismissing(m.color_override) ? nothing : String(m.color_override),
            :label_override => ismissing(m.label_override) ? nothing : String(m.label_override),
            :q_window_min   => ismissing(m.q_window_min)   ? nothing : Float64(m.q_window_min),
            :q_window_max   => ismissing(m.q_window_max)   ? nothing : Float64(m.q_window_max),
            # peak_display/snapshot are stored JSON text; parse so the payload
            # round-trips as structured JSON (matches the route's parsed shape).
            :peak_display   => ismissing(m.peak_display) ? nothing : JSON3.read(String(m.peak_display)),
            :snapshot       => JSON3.read(String(m.snapshot)),
        ))
    end
    Dict{Symbol,Any}(:members => members)
end

# Raw-INSERT one synthesized user_actions row (no broadcast, NULL client ids,
# carrying the comparison's user_id + historical timestamp), then fold its
# payload through the live dispatcher. Returns the new event_id.
# NOTE: the timestamp column is named `timestamp`, NOT `created_at`.
function _synthesize_series_event!(db::SQLite.DB, kind::String, series_id::Integer,
                                   payload::Dict, user_id, ts)
    payload_json = JSON3.write(payload)
    res = DBInterface.execute(db,
        """INSERT INTO user_actions
             (user_id, action, entity_type, entity_id, payload,
              undoes_event_id, client_id, client_op_id, timestamp)
           VALUES (?, ?, 'series', ?, ?, NULL, NULL, NULL, ?)""",
        [user_id, kind, Int(series_id), payload_json, ts])
    event_id = Int(DBInterface.lastrowid(res))
    # Canonicalize exactly as apply_event! does so the dispatcher branches see a
    # JSON3.Object, not a Dict.
    payload_canonical = JSON3.read(payload_json)
    update_view_for_event!(db, kind, Int(series_id), payload_canonical, event_id)
    event_id
end

"""
    migrate_comparisons_to_series!(db)

Event-sourced copy of every `comparisons` row into the `series*` tables
(master plan §6.1, issue #171). Runs at `migrate_schema!` time, after
`migrate_series!` (so the `series*` tables and `schema_migrations` sentinel
exist) and after `migrate_compare_view_choices!`/`migrate_compare_relax_nullability!`
(it reads the columns those add). Wrapped in its own `SQLite.transaction`.

Idempotent via the `schema_migrations` sentinel (gate at start, marker row
written LAST inside the same transaction — the gate flips only on a fully
committed copy). NEVER calls `apply_event!` (its public method broadcasts; a
migration must not fan out N replay-as-reruns). Raw-`INSERT`s `user_actions`
rows with `client_op_id`/`client_id` NULL — so no `idempotent_responses` cache
rows are written and the partial unique index does not apply.
"""
function migrate_comparisons_to_series!(db::SQLite.DB)
    # Gate: skip if already applied. A single read before the transaction is
    # sufficient — migrate_schema! runs single-threaded inside open_db, before
    # serve accepts connections, so there is no concurrent opener to race.
    already = Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM schema_migrations WHERE name = ?",
        [MIGRATION_COMPARISONS_TO_SERIES]))
    isempty(already) || return nothing

    SQLite.transaction(db) do
        # NOTE: `forked_from_id`/`forked_at_hash` are deliberately NOT selected —
        # the migrated series always NULLs its fork lineage (id-space mismatch; see
        # `_series_created_payload_from_comparison`). Selecting them would be a dead read.
        cmps = Tables.rowtable(DBInterface.execute(db,
            """SELECT id, title, description, content_hash, created_by, created_at,
                      updated_at,
                      view_grouping_mode, view_show_peak_ticks, view_show_peak_labels
               FROM comparisons ORDER BY id ASC"""))

        for cmp in cmps
            uid = ismissing(cmp.created_by) ? nothing : Int(cmp.created_by)
            # Historical timestamp carried onto the synthesized event rows'
            # `timestamp` column. Falls back to now() only if the comparison
            # somehow has a NULL created_at (post-#67 it is nullable).
            ts  = ismissing(cmp.created_at) ? comparison_now_iso() : String(cmp.created_at)

            # Mint the series id exactly as the live route does.
            res = DBInterface.execute(db, "INSERT INTO series DEFAULT VALUES")
            new_id = Int(DBInterface.lastrowid(res))

            created_payload = _series_created_payload_from_comparison(db, cmp)
            _synthesize_series_event!(db, "series_created", new_id, created_payload, uid, ts)

            plate_payload = _series_plate_committed_payload_from_comparison(db, cmp)
            _synthesize_series_event!(db, "series_plate_committed", new_id, plate_payload, uid, ts)

            # Messages: carry author_id + body + created_at (don't let the
            # DEFAULT CURRENT_TIMESTAMP overwrite history). series_messages.id
            # auto-assigns.
            DBInterface.execute(db,
                """INSERT INTO series_messages (series_id, author_id, body, created_at)
                   SELECT ?, author_id, body, created_at
                   FROM comparison_messages WHERE comparison_id = ?
                   ORDER BY created_at ASC, id ASC""", [new_id, Int(cmp.id)])

            # Pins: carry user_id + pinned_at. OR IGNORE defends against a dup
            # (user, series) PK (can't actually occur — comparison_pins PK is
            # (user, comparison) and one comparison maps to one series).
            DBInterface.execute(db,
                """INSERT OR IGNORE INTO series_pins (user_id, series_id, pinned_at)
                   SELECT user_id, ?, pinned_at
                   FROM comparison_pins WHERE comparison_id = ?""", [new_id, Int(cmp.id)])
        end

        # Sentinel marker LAST, inside the same transaction.
        DBInterface.execute(db,
            "INSERT INTO schema_migrations (name, applied_at) VALUES (?, ?)",
            [MIGRATION_COMPARISONS_TO_SERIES, comparison_now_iso()])
    end
    nothing
end

# Heal FK references in `sqlite_master.sql` that point at `_migrate_old_<entity>`
# back to `<entity>`. SQLite's ALTER TABLE RENAME tracking rewrites stored FK
# refs in EVERY table whose CREATE statement named the renamed entity — those
# stale refs surface later as `no such table: main._migrate_old_*` on unrelated
# INSERTs (the prepare step walks the FK graph). Mirrors the strategy in
# `_fix_fk_references_after_autoincrement_migration!` but scoped to one entity
# (the `comparisons` rebuild for #67).
#
# When to use this vs. _fix_fk_references_after_autoincrement_migration!:
# - one-off scoped rebuild of a single entity (e.g. comparisons for #67) →
#   call this helper directly with the entity name.
# - entity participates in the autoincrement / R2 widen rebuild → add it to
#   `_MIGRATION_TEMP_ENTITIES` instead so the existing iteration covers it.
function _heal_renamed_table_fk_refs!(db::SQLite.DB, entity::String)
    old_name = "_migrate_old_$entity"
    broken = Tables.rowtable(DBInterface.execute(db,
        "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE ?",
        ["%$old_name%"]))
    isempty(broken) && return
    DBInterface.execute(db, "PRAGMA writable_schema = ON")
    try
        SQLite.transaction(db) do
            # Quoted ("…") and bare forms both occur depending on how the FK
            # target was rendered in the original CREATE.
            DBInterface.execute(db,
                "UPDATE sqlite_master SET sql = REPLACE(sql, ?, ?) WHERE type='table'",
                ["\"$old_name\"", entity])
            DBInterface.execute(db,
                "UPDATE sqlite_master SET sql = REPLACE(sql, ?, ?) WHERE type='table'",
                [old_name, entity])
            remaining = Tables.rowtable(DBInterface.execute(db,
                "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE ?",
                ["%$old_name%"]))
            isempty(remaining) || error(
                "_heal_renamed_table_fk_refs!($entity): tables still carry " *
                "`$old_name` FKs after heal: " *
                join(String[String(r.name) for r in remaining], ", "))
        end
    finally
        DBInterface.execute(db, "PRAGMA writable_schema = OFF")
    end
    # Force SQLite to invalidate its in-memory schema cache. VACUUM is the
    # standard idiom but acquires an EXCLUSIVE lock and can fail with
    # SQLITE_BUSY under multi-process access — fall back to bumping
    # schema_version (also invalidates the cache).
    try
        DBInterface.execute(db, "VACUUM")
    catch err
        @warn "VACUUM after $entity FK heal failed; bumping schema_version" exception=err
        DBInterface.execute(db, "PRAGMA writable_schema = ON")
        try
            DBInterface.execute(db, "PRAGMA schema_version = schema_version + 1")
        finally
            DBInterface.execute(db, "PRAGMA writable_schema = OFF")
        end
    end
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
    legacy = first(Tables.rowtable(DBInterface.execute(db, """
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
populated tables like `sample_messages`, `sample_tags`, `exposure_tags`,
`exposure_sources`, `index_groups`, `index_group_members` whose FKs had
the same corruption. A real prod DB at `/tmp/himalaya-dev.db` with rows
in `sample_tags`/`sample_messages` exhibited the bug — caught only by
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
    broken = Tables.rowtable(DBInterface.execute(db,
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
            remaining = Tables.rowtable(DBInterface.execute(db,
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
    cols = Set(r.name for r in Tables.rowtable(DBInterface.execute(db,
        "PRAGMA table_info('samples')")))
    # No samples table yet (partial legacy fixture or pre-create_schema! call) — skip.
    # create_schema! will create it with the canonical shape; nothing to rename.
    isempty(cols) && return nothing
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
        # Duplicate suffix pass (oldest id keeps bare name).
        # Collision-safe: track existing (experiment_id, name) pairs so a user-named
        # sample literally called "<name>-2" doesn't conflict with our rename target.
        existing = Set{Tuple{Int64,String}}(
            (Int64(r.experiment_id), String(r.name))
            for r in Tables.rowtable(DBInterface.execute(db,
                "SELECT experiment_id, name FROM samples WHERE name IS NOT NULL AND experiment_id IS NOT NULL")))
        dups = Tables.rowtable(DBInterface.execute(db, """
            SELECT experiment_id, name FROM samples
            GROUP BY experiment_id, name HAVING COUNT(*) > 1"""))
        for d in dups
            ids = Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM samples WHERE experiment_id = ? AND name = ? ORDER BY id ASC",
                [d.experiment_id, d.name]))
            for (i, row) in enumerate(ids)
                i == 1 && continue  # oldest keeps the bare name
                # Pick the next suffix that isn't already taken by a user-named sample.
                suffix_n = i
                new_name = "$(d.name)-$(suffix_n)"
                while (Int64(d.experiment_id), new_name) in existing
                    suffix_n += 1
                    new_name = "$(d.name)-$(suffix_n)"
                end
                push!(existing, (Int64(d.experiment_id), new_name))
                @warn "Renamed duplicate sample" experiment_id=d.experiment_id old=d.name new=new_name id=row.id
                DBInterface.execute(db, "UPDATE samples SET name = ? WHERE id = ?",
                    [new_name, row.id])
            end
        end
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
            for r in Tables.rowtable(DBInterface.execute(db,
                "SELECT sample_id, filename FROM exposures WHERE sample_id IS NOT NULL AND filename IS NOT NULL")))

        dups = Tables.rowtable(DBInterface.execute(db, """
            SELECT sample_id, filename FROM exposures
            WHERE sample_id IS NOT NULL AND filename IS NOT NULL
            GROUP BY sample_id, filename HAVING COUNT(*) > 1"""))

        for d in dups
            ids = Tables.rowtable(DBInterface.execute(db,
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
    info = Tables.rowtable(DBInterface.execute(db,
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
        name::Union{String,Nothing}         = nothing,
        display_name::Union{String,Nothing} = nothing,
        notes::Union{String,Nothing}        = nothing)
    result = DBInterface.execute(db,
        "INSERT INTO samples (experiment_id, name, display_name, notes) VALUES (?, ?, ?, ?)",
        [experiment_id, name, display_name, notes])
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

    # WAL lets concurrent readers proceed alongside one writer — load-bearing
    # for parallel request handling (#115). The default rollback journal
    # serializes every reader behind any in-flight writer. WAL persists in
    # the DB file header, so this PRAGMA is effectively a one-time migration
    # for existing DBs and a no-op on every subsequent open. Skipped for
    # `:memory:` and shared-cache URI DBs, which can't run WAL.
    #
    # `Tables.rowtable(...)` drains the result iterator so the prepared
    # statement is dropped before `finalize_statements!`. Without the drain,
    # the iterator keeps the PRAGMA's statement attached, and any subsequent
    # DDL (e.g. `DROP TABLE` in tests, or migration loops on legacy DBs)
    # fails with `database table is locked`. The trailing
    # `finalize_statements!` then clears the cache so callers see a
    # quiescent connection.
    if db_path != ":memory:" && !startswith(db_path, "file:")
        Tables.rowtable(DBInterface.execute(db, "PRAGMA journal_mode = WAL"))
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
