# Series Database Schema (I2.1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the series data model to the SQLite schema — five new `series*` tables plus a `schema_migrations` sentinel — via a new idempotent `migrate_series!` registered in the `migrate_schema!` sequence.

**Architecture:** `migrate_series!` mirrors the existing `migrate_compare!` (`packages/HimalayaUI/src/db.jl:449`): a standalone function of `CREATE TABLE IF NOT EXISTS` statements, registered as one line in `migrate_schema!`. Because `series*` are brand-new tables, each is created in its **final shape in a single `CREATE TABLE`** — unlike the `comparison*` tables, which reached their shape through three migrations (`migrate_compare!` + `migrate_compare_view_choices!` + `migrate_compare_relax_nullability!`). The `comparison*` tables are untouched. No routes, event kinds, or data migration — those are #165–#168 and #171.

**Tech Stack:** Julia, SQLite.jl, stdlib `Test`.

**Source spec:** GitHub issue #164 (table shapes + acceptance criteria), master plan §5.1/§9 (`docs/superpowers/plans/2026-05-17-himalaya-ui-redesign.md`).

**One deliberate read of intent:** issue #164 writes `series_samples` FK columns as `series_id REFERENCES series ...` / `sample_id REFERENCES samples ...` without an explicit `NOT NULL`. This plan adds `NOT NULL` to both: a recipe row with a NULL `series_id` is unreachable, and master plan §5.1 states an orphan (NULL `sample_id`) is "unrenderable" — the reason `sample_id` is `CASCADE` not `SET NULL`. This matches the `comparison_members.comparison_id NOT NULL` precedent. Flag for reviewer confirmation.

---

## File Structure

- **Modify** `packages/HimalayaUI/src/db.jl` — add the `migrate_series!` function (after `migrate_compare_view_choices!`, ~line 663) and register one call in `migrate_schema!` (after the `migrate_compare_view_choices!(db)` call, ~line 368).
- **Modify** `packages/HimalayaUI/test/test_db.jl` — add `migrate_series!` to the `using HimalayaUI: …` import list; add two `@testset`s at end of file.

No new files. `migrate_series!` is a migration function, not part of the base `SCHEMA` constant — same as `migrate_compare!`.

---

## Task 1: `series` parent table + `schema_migrations` sentinel

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (add `migrate_series!`, register in `migrate_schema!`)
- Test: `packages/HimalayaUI/test/test_db.jl`

- [ ] **Step 1: Add `migrate_series!` to the test import list**

In `packages/HimalayaUI/test/test_db.jl`, the file starts with:

```julia
using HimalayaUI: create_schema!, migrate_schema!, create_experiment!, create_sample!,
                  create_exposure!, get_experiment, get_samples, get_exposures,
                  migrate_r2_split_peaks!
```

Add `migrate_series!` to the list:

```julia
using HimalayaUI: create_schema!, migrate_schema!, create_experiment!, create_sample!,
                  create_exposure!, get_experiment, get_samples, get_exposures,
                  migrate_r2_split_peaks!, migrate_series!
```

- [ ] **Step 2: Write the failing test for `series` + `schema_migrations`**

Append to the end of `packages/HimalayaUI/test/test_db.jl`:

```julia
@testset "migrate_series! creates series + schema_migrations" begin
    db = SQLite.DB()  # in-memory
    create_schema!(db)
    migrate_series!(db)

    tables = Set(r.name for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT name FROM sqlite_master WHERE type='table'")))
    @test "series" in tables
    @test "schema_migrations" in tables

    cols = Set(r.name for r in Tables.rowtable(DBInterface.execute(db,
        "PRAGMA table_info(series)")))
    # comparison columns + view-choice columns + recipe columns
    for c in ("id", "title", "description", "content_hash", "created_by",
              "created_at", "updated_at", "forked_from_id", "forked_at_hash",
              "view_grouping_mode", "view_show_peak_ticks", "view_show_peak_labels",
              "ordering_variable", "order_rule", "state")
        @test c in cols
    end

    # schema_migrations ships empty — #171 writes the marker row, not #164.
    @test only(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM schema_migrations"))).n == 0

    # order_rule / state CHECK constraints reject bad values.
    DBInterface.execute(db, "INSERT INTO series DEFAULT VALUES")  # defaults pass CHECK
    @test_throws SQLite.SQLiteException DBInterface.execute(db,
        "INSERT INTO series (order_rule) VALUES ('sideways')")
    @test_throws SQLite.SQLiteException DBInterface.execute(db,
        "INSERT INTO series (state) VALUES ('archived')")

    # Idempotent re-run.
    migrate_series!(db)
    @test "series" in Set(r.name for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT name FROM sqlite_master WHERE type='table'")))
end

@testset "migrate_schema! installs series tables on a legacy DB" begin
    db = SQLite.DB()
    create_schema!(db)
    migrate_schema!(db)  # migrate_series! is registered in the sequence

    tables = Set(r.name for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT name FROM sqlite_master WHERE type='table'")))
    @test "series" in tables
    @test "schema_migrations" in tables

    # The comparison* tables are untouched by this change.
    for c in ("comparisons", "comparison_members", "comparison_messages",
              "comparison_pins")
        @test c in tables
    end

    migrate_schema!(db)  # second run is idempotent
end
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_db.jl`
Expected: FAIL — `UndefVarError: migrate_series! not defined` (the import in Step 1 references a symbol that does not yet exist).

- [ ] **Step 4: Add the `migrate_series!` function**

In `packages/HimalayaUI/src/db.jl`, immediately after the end of `migrate_compare_view_choices!` (the closing `end` at ~line 663), add:

```julia
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

    # Migration-version sentinel. No such table exists today; #171's copy
    # needs a real marker. Created empty here — #164 writes no rows.
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS schema_migrations (
            name        TEXT PRIMARY KEY,
            applied_at  TEXT
        )""")
    nothing
end
```

- [ ] **Step 5: Register `migrate_series!` in `migrate_schema!`**

In `migrate_schema!` (`packages/HimalayaUI/src/db.jl`), find the `migrate_compare_view_choices!` call (~line 368):

```julia
    # 2026-05-14 — Compare UX refinement (spec §6.4): persist view choices
    # on the comparison so they round-trip across viewers.
    migrate_compare_view_choices!(db)
```

Add immediately after it (before the `migrate_experiment_config_label_to_name!` block):

```julia
    # Series model (#164, master plan §5.1): new `series*` tables + the
    # `schema_migrations` sentinel. Placed after the compare migrations so
    # the compare/series schema stays grouped and #171's
    # `migrate_comparisons_to_series!` has a natural slot after this. New
    # tables only — the `comparison*` tables are never renamed (§2.1).
    migrate_series!(db)
```

- [ ] **Step 6: Run the test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_db.jl`
Expected: PASS — both `@testset`s green; no failures elsewhere in `test_db.jl`.

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_db.jl
git commit -m "feat: add series + schema_migrations tables (I2.1, #164)"
```

---

## Task 2: The four `series*` child tables

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (extend `migrate_series!`)
- Test: `packages/HimalayaUI/test/test_db.jl`

- [ ] **Step 1: Write the failing test for the child tables + constraints**

Append to the end of `packages/HimalayaUI/test/test_db.jl`:

```julia
@testset "migrate_series! creates the four series child tables" begin
    db = SQLite.DB()
    create_schema!(db)
    migrate_series!(db)

    tables = Set(r.name for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT name FROM sqlite_master WHERE type='table'")))
    for t in ("series_members", "series_samples", "series_messages", "series_pins")
        @test t in tables
    end

    # series_samples CHECK constraints on pinned / excluded.
    DBInterface.execute(db, "INSERT INTO series DEFAULT VALUES")
    DBInterface.execute(db, """
        INSERT INTO series_samples (series_id, sample_id, position)
        VALUES (1, 1, 0)""")  # defaults: pinned=0, excluded=0 — pass CHECK
    @test_throws SQLite.SQLiteException DBInterface.execute(db, """
        INSERT INTO series_samples (series_id, sample_id, position, pinned)
        VALUES (1, 1, 1, 2)""")
    @test_throws SQLite.SQLiteException DBInterface.execute(db, """
        INSERT INTO series_samples (series_id, sample_id, position, excluded)
        VALUES (1, 1, 2, 9)""")

    # UNIQUE(series_id, position) — position 0 is already taken above.
    @test_throws SQLite.SQLiteException DBInterface.execute(db, """
        INSERT INTO series_samples (series_id, sample_id, position)
        VALUES (1, 1, 0)""")

    # series_members.snapshot must be valid JSON.
    @test_throws SQLite.SQLiteException DBInterface.execute(db, """
        INSERT INTO series_members (series_id, display_order, snapshot, created_at)
        VALUES (1, 0, 'not-json', '2026-05-17T00:00:00.000Z')""")

    migrate_series!(db)  # idempotent re-run
end

@testset "series child tables cascade-delete with the series" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))  # FK enforcement ON
        exp_id  = create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
        samp_id = create_sample!(db, exp_id, "s1")

        DBInterface.execute(db, "INSERT INTO series DEFAULT VALUES")
        sid = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT last_insert_rowid() AS id"))).id

        DBInterface.execute(db, """
            INSERT INTO series_samples (series_id, sample_id, position)
            VALUES (?, ?, 0)""", [sid, samp_id])
        DBInterface.execute(db, """
            INSERT INTO series_members (series_id, display_order, snapshot, created_at)
            VALUES (?, 0, '{}', '2026-05-17T00:00:00.000Z')""", [sid])
        DBInterface.execute(db, """
            INSERT INTO series_messages (series_id, body) VALUES (?, 'hi')""", [sid])
        DBInterface.execute(db, """
            INSERT INTO series_pins (user_id, series_id, pinned_at)
            VALUES (1, ?, '2026-05-17T00:00:00.000Z')""", [sid])

        DBInterface.execute(db, "DELETE FROM series WHERE id = ?", [sid])

        for t in ("series_samples", "series_members", "series_messages", "series_pins")
            n = only(Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS n FROM $t WHERE series_id = ?", [sid]))).n
            @test n == 0
        end
    end
end
```

> Note: `series_pins.user_id` references `users(id)`; with FK enforcement ON,
> `open_db` runs `migrate_schema!`, which seeds no users. If the cascade test
> fails on the `series_pins` insert with a foreign-key error, insert a user
> first: `DBInterface.execute(db, "INSERT INTO users (username) VALUES ('t')")`
> then use that id. Check the `users` table shape via `PRAGMA table_info(users)`
> before assuming the column name.

- [ ] **Step 2: Run the test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_db.jl`
Expected: FAIL — `series_members`/`series_samples`/`series_messages`/`series_pins` do not exist (`no such table`).

- [ ] **Step 3: Extend `migrate_series!` with the four child tables**

In `packages/HimalayaUI/src/db.jl`, inside `migrate_series!`, insert the four child-table statements **after** the `schema_migrations` `CREATE TABLE` and **before** the final `nothing`:

```julia
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
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_db.jl`
Expected: PASS — all four `@testset`s from Tasks 1–2 green, no regressions in `test_db.jl`.

- [ ] **Step 5: Run the full backend suite once to confirm no regression**

Run: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1`
Then: `grep -E "Test Summary|did not pass|fail" /tmp/jl-test.out`
Expected: no `did not pass` / `fail` lines. The suite is slow (5–10 min) — capture once, grep the file (see `packages/HimalayaUI/test/AGENTS.md`).

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_db.jl
git commit -m "feat: add series_members/samples/messages/pins tables (I2.1, #164)"
```

---

## Acceptance Criteria (issue #164)

- [ ] A fresh DB and `migrate_schema!` on an existing DB both produce all five `series*` tables and `schema_migrations` — Task 1 Step 2 (`migrate_schema! installs series tables on a legacy DB`) + Task 2 Step 1.
- [ ] All `CHECK` constraints, `UNIQUE(series_id, position)`, and `ON DELETE CASCADE` clauses are present — Task 1 (`order_rule`/`state` CHECK) + Task 2 (`pinned`/`excluded` CHECK, `UNIQUE`, `json_valid`, cascade-delete).
- [ ] `migrate_series!` is idempotent on re-run — Task 1 + Task 2 each re-run `migrate_series!` and re-assert.
- [ ] The `comparison*` tables are unchanged — Task 1 Step 2 asserts all four `comparison*` tables still present; no statement in `migrate_series!` touches them.
- [ ] A Julia schema test covers the new tables and constraints — four `@testset`s in `test_db.jl`.
