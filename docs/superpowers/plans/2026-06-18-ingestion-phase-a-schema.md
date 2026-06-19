# Ingestion Redesign — Phase A: Schema & DB Layer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Lay the database foundation for the ingestion redesign: the additive migrations (a `loads` table, denormalized `exposures.experiment_id`, typed geometry/scan columns on `experiments`, the `samples` two-label→one-label collapse) plus the `create_*` function and `analyze_exposure!` signature changes they force, so later phases can scan, group, and analyze without touching the data model again.

**Architecture:** Extend the existing sentinel-gated migration system in `packages/HimalayaUI/src/db.jl` (each migration is a named row in `schema_migrations`, applied once, idempotent DDL inside a transaction). Mirror the existing `migrate_exposures_unique_filename!` dedupe-then-enforce discipline. All changes are additive except the `samples` machine-`name` drop, which is a careful drop-index → drop-column → rename sequence. No grouping/scanning logic here (that is Phase B) — this plan only changes shapes and the functions that write them.

**Tech Stack:** Julia, SQLite.jl / DBInterface.jl, stdlib `Test`. Backend package at `packages/HimalayaUI/`.

**Spec:** `docs/superpowers/specs/2026-06-15-ingestion-redesign-design.md` §4, §9.1 (`analyze_exposure!` rewire), §10 (schema). Read it before starting.

**Source of truth for current code:** the build kit anchors in this plan were line-verified 2026-06-18, but line numbers drift — confirm each anchor with a quick `grep`/read before editing.

---

## File Structure

| File | Responsibility | This plan |
|---|---|---|
| `packages/HimalayaUI/src/db.jl` | schema DDL, migration driver, `create_*` writers | MODIFY (all tasks) |
| `packages/HimalayaUI/src/pipeline.jl` | `analyze_exposure!` (experiment resolution) | MODIFY (Task 9) |
| `packages/HimalayaUI/src/routes_experiments.jl` | the one API caller of `analyze_exposure!` | MODIFY (Task 9) |
| `packages/HimalayaUI/test/test_ingestion_schema.jl` | new standalone test file for this phase | CREATE (all tasks) |
| `packages/HimalayaUI/test/runtests.jl` | test registry | MODIFY (Task 1) |

**Out of scope (later plans):** `prp.jl` / `geometry.jl` / `grouping.jl` / `scan_and_group!` (Phase B), the scan/rescan/SSE routes + scheduler (Phase C), the structural-edit event kinds (Phase D), the whole frontend + the `display_name→name` TS sweep (Phase E).

---

## Conventions for every task

- **Run a single test file** during TDD from the repo root:
  `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
  (Standalone: the file ends with its own `@testset`s and runs directly. The full suite via `Pkg.test` stays for CI.)
- Each test builds a **temp DB** (`mktempdir`) so tests never touch a real database.
- **Commit after each task** once its test passes.
- The migration sentinel name pattern is `"MIGRATION_<THING>"` (string consts at the top of `db.jl`).

---

## Task 0: Test harness + helpers

**Files:**
- Create: `packages/HimalayaUI/test/test_ingestion_schema.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Create the standalone test file with a pre-redesign DB fixture helper**

```julia
# packages/HimalayaUI/test/test_ingestion_schema.jl
using Test
using SQLite, DBInterface, Tables
using HimalayaUI

"Build a fresh DB at the CURRENT (post-migration) schema and return its path."
function fresh_db()
    dir = mktempdir()
    path = joinpath(dir, "test.db")
    db = SQLite.DB(path)
    HimalayaUI.create_schema!(db)    # canonical tables
    HimalayaUI.migrate_schema!(db)   # NOTE: takes an SQLite.DB, not a path (db.jl:247)
    SQLite.close(db)
    return path
end

"Insert a minimal experiment honouring its NOT NULL columns (name/path/data_dir/analysis_dir, db.jl:21-23). Returns the id."
function seed_experiment(db; name="e")
    DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES (?, ?, ?, ?)",
        [name, "/exp/$name", "/exp/$name/data", "/exp/$name/analysis"]))
end

"Open a path, run f(db), always close."
function with_db(f, path)
    db = SQLite.DB(path)
    try
        return f(db)
    finally
        SQLite.close(db)
    end
end

"Column names of a table, lowercased."
cols(db, table) = lowercase.(getproperty.(Tables.rowtable(
    DBInterface.execute(db, "PRAGMA table_info($table)")), :name))

"Names of indexes on a table."
indexes(db, table) = String.(getproperty.(Tables.rowtable(
    DBInterface.execute(db, "PRAGMA index_list($table)")), :name))

@testset "ingestion schema (Phase A)" begin
    # tasks append @testset blocks below
end
```

> **Notes (verified against live `db.jl`):** there is **no `init_db!`** — the production entry is `open_db` (`db.jl:1894`), which is a singleton-connection opener unsuitable for per-test temp DBs, so tests build the schema explicitly with `create_schema!(db)` + `migrate_schema!(db)` as above. `migrate_schema!` takes an `SQLite.DB`, **not** a path (`db.jl:247`). The test helper `cols(db, table)` below is **distinct** from the production helper `cols_of` defined in `db.jl` (Task 3) — don't conflate them. Every test that needs an experiment must use `seed_experiment(db)` (never a bare `INSERT INTO experiments (name) …`, which violates the 3 NOT NULL columns).

- [ ] **Step 2: Register the file in runtests.jl**

Find the list of `include(...)` lines in `packages/HimalayaUI/test/runtests.jl` and add, in alphabetical/topical position:

```julia
include("test_ingestion_schema.jl")
```

- [ ] **Step 3: Run it to confirm the harness loads**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: PASS (an empty `@testset` passes; it proves `init_db!`/imports resolve).

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/test/test_ingestion_schema.jl packages/HimalayaUI/test/runtests.jl
git commit -m "test: scaffold ingestion schema test harness"
```

---

## Task 1: Migration sentinels + driver ordering

Add the named sentinels for the three new migrations and wire them into the driver **after** `migrate_samples_naming!` (the collapse depends on it having run).

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (sentinel consts near the other `MIGRATION_*` consts; the `migrate_schema!` driver body)
- Test: `packages/HimalayaUI/test/test_ingestion_schema.jl`

- [ ] **Step 1: Write the failing test**

Append to the `@testset "ingestion schema (Phase A)"`:

```julia
@testset "migrations recorded" begin
    path = fresh_db()
    with_db(path) do db
        applied = Set(String.(getproperty.(Tables.rowtable(
            DBInterface.execute(db, "SELECT name FROM schema_migrations")), :name)))
        @test HimalayaUI.MIGRATION_LOADS_TABLE in applied
        @test HimalayaUI.MIGRATION_EXPOSURES_EXPERIMENT_ID in applied
        @test HimalayaUI.MIGRATION_EXPERIMENTS_GEOMETRY in applied
        @test HimalayaUI.MIGRATION_SAMPLES_NAME_COLLAPSE in applied
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: FAIL (the three sentinels are not in `schema_migrations`).

- [ ] **Step 3: Add the sentinels + driver calls**

In `db.jl`, next to the existing `const MIGRATION_* = "..."` declarations, add the **four** sentinels for this phase. Match the existing lowercase-short value convention (e.g. `"comparisons_to_series"`, `"assignments_v1"`):

```julia
const MIGRATION_LOADS_TABLE             = "loads_table_v1"
const MIGRATION_EXPOSURES_EXPERIMENT_ID = "exposures_experiment_id_v1"
const MIGRATION_EXPERIMENTS_GEOMETRY    = "experiments_geometry_v1"
const MIGRATION_SAMPLES_NAME_COLLAPSE   = "samples_name_collapse_v1"
```

The helpers `_migrated`/`_record_migration!` **do not exist** today (`migrate_assignments!` inlines the gate-read + sentinel INSERT). Define them once in `db.jl` (use `comparison_now_iso()`, which lives in `comparisons.jl:169` but resolves at runtime across the package):

```julia
_migrated(db::SQLite.DB, name::AbstractString) = !isempty(Tables.rowtable(DBInterface.execute(db,
    "SELECT 1 FROM schema_migrations WHERE name = ?", [name])))
_record_migration!(db::SQLite.DB, name::AbstractString) = DBInterface.execute(db,
    "INSERT INTO schema_migrations (name, applied_at) VALUES (?, ?)", [name, comparison_now_iso()])
```

In `migrate_schema!`, **after** the existing `migrate_samples_naming!(db)` call, add the **four** calls in dependency order (loads → exposures denorm → experiments geometry → samples collapse). This is the canonical final sequence:

```julia
    migrate_loads_table!(db)
    migrate_exposures_experiment_id!(db)
    migrate_experiments_geometry!(db)
    migrate_samples_name_collapse!(db)
```

Add **four** stub functions (replaced in Tasks 2/3/4/5) that only write their sentinel so this task's test passes:

```julia
function migrate_loads_table!(db::SQLite.DB)
    _migrated(db, MIGRATION_LOADS_TABLE) && return nothing
    SQLite.transaction(db) do; _record_migration!(db, MIGRATION_LOADS_TABLE); end
end
function migrate_exposures_experiment_id!(db::SQLite.DB)
    _migrated(db, MIGRATION_EXPOSURES_EXPERIMENT_ID) && return nothing
    SQLite.transaction(db) do; _record_migration!(db, MIGRATION_EXPOSURES_EXPERIMENT_ID); end
end
function migrate_experiments_geometry!(db::SQLite.DB)
    _migrated(db, MIGRATION_EXPERIMENTS_GEOMETRY) && return nothing
    SQLite.transaction(db) do; _record_migration!(db, MIGRATION_EXPERIMENTS_GEOMETRY); end
end
function migrate_samples_name_collapse!(db::SQLite.DB)
    _migrated(db, MIGRATION_SAMPLES_NAME_COLLAPSE) && return nothing
    SQLite.transaction(db) do; _record_migration!(db, MIGRATION_SAMPLES_NAME_COLLAPSE); end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_ingestion_schema.jl
git commit -m "feat(db): scaffold ingestion migration sentinels + driver order"
```

---

## Task 2: `loads` table

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (`migrate_loads_table!`)
- Test: `packages/HimalayaUI/test/test_ingestion_schema.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "loads table" begin
    path = fresh_db()
    with_db(path) do db
        @test "loads" in lowercase.(String.(getproperty.(Tables.rowtable(
            DBInterface.execute(db, "SELECT name FROM sqlite_master WHERE type='table'")), :name)))
        c = cols(db, "loads")
        for col in ["id","experiment_id","load_index","session_id","start_time","end_time","frame_count","note"]
            @test col in c
        end
        # FK to experiments enforced
        DBInterface.execute(db, "PRAGMA foreign_keys=ON")
        @test_throws Exception DBInterface.execute(db,
            "INSERT INTO loads (experiment_id, load_index) VALUES (999999, 0)")
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: FAIL ("no such table: loads").

- [ ] **Step 3: Implement the migration**

Replace the `migrate_loads_table!` stub:

```julia
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_ingestion_schema.jl
git commit -m "feat(db): add loads table migration"
```

---

## Task 3: `exposures` — new columns, `experiment_id` backfill, dedup index swap

Mirror the dedupe-then-enforce discipline of `migrate_exposures_unique_filename!` (`db.jl:1516-1554`). Fail-fast on any exposure whose `experiment_id` cannot be derived (spec §10, P1-5).

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (`migrate_exposures_experiment_id!`)
- Test: `packages/HimalayaUI/test/test_ingestion_schema.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "exposures experiment_id + dedup" begin
    path = fresh_db()
    with_db(path) do db
        c = cols(db, "exposures")
        for col in ["experiment_id","prp_path","timestamp","exposure_time",
                    "horizontal_position","scan_id","frame_no","load_id","content_fingerprint"]
            @test col in c
        end
        # new dedup index present, old one gone
        idx = indexes(db, "exposures")
        @test any(i -> occursin("exposures_unique_filename", i), idx)
        # the new index is on (experiment_id, filename): same filename, different sample, same experiment → rejected
        DBInterface.execute(db, "PRAGMA foreign_keys=ON")
        eid = seed_experiment(db)   # honours NOT NULL name/path/data_dir/analysis_dir
        s1 = DBInterface.lastrowid(DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, name) VALUES (?, 'A')", [eid]))
        s2 = DBInterface.lastrowid(DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, name) VALUES (?, 'B')", [eid]))
        DBInterface.execute(db,
            "INSERT INTO exposures (experiment_id, sample_id, filename) VALUES (?, ?, 'f.tif')", [eid, s1])
        @test_throws Exception DBInterface.execute(db,
            "INSERT INTO exposures (experiment_id, sample_id, filename) VALUES (?, ?, 'f.tif')", [eid, s2])
    end
end
```

> If `experiments`/`samples` minimal INSERTs above fail because of NOT NULL columns, add the required columns to the INSERTs after checking `cols(db, "experiments")` / `cols(db, "samples")`.

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: FAIL (no `experiment_id` column).

- [ ] **Step 3: Implement the migration**

Replace the `migrate_exposures_experiment_id!` stub:

```julia
function migrate_exposures_experiment_id!(db::SQLite.DB)
    _migrated(db, MIGRATION_EXPOSURES_EXPERIMENT_ID) && return nothing
    SQLite.transaction(db) do
        existing = cols_of(db, "exposures")   # see helper note below
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
        orphans = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, filename FROM exposures WHERE experiment_id IS NULL"))
        isempty(orphans) || error(
            "migrate_exposures_experiment_id!: $(length(orphans)) exposures have no derivable " *
            "experiment_id (orphaned sample_id). Resolve manually before migrating. ids: " *
            join(getproperty.(orphans, :id), ", "))

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
```

> **Helper note:** add a `cols_of(db, table)` next to the migration helpers if one doesn't exist:
> ```julia
> cols_of(db, table) = String.(getproperty.(Tables.rowtable(
>     DBInterface.execute(db, "PRAGMA table_info($table)")), :name))
> ```
> Guarding each ADD COLUMN with `name in existing` keeps the migration idempotent on partially-migrated DBs.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_ingestion_schema.jl
git commit -m "feat(db): denormalize exposures.experiment_id + swap dedup key"
```

---

## Task 4: `experiments` — typed geometry + scan columns

All additive `ALTER TABLE ADD COLUMN`. Per spec §10: geometry typed cols + a `*_source` per geometry field, plus scan/scheduler bookkeeping. `flight_path_m`/`energy_kev` already exist; add the rest. There is **no** `detector_distance_mm` column (spec §10: one quantity → `flight_path_m`).

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (extend an existing additive-columns migration, or add `migrate_experiments_geometry!` gated by a new sentinel `MIGRATION_EXPERIMENTS_GEOMETRY`)
- Test: `packages/HimalayaUI/test/test_ingestion_schema.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "experiments geometry+scan columns" begin
    path = fresh_db()
    with_db(path) do db
        c = cols(db, "experiments")
        for col in ["beam_center_x","beam_center_y","pixel_size_um","q_units",
                    "energy_kev_source","flight_path_m_source","beam_center_x_source",
                    "beam_center_y_source","pixel_size_um_source","q_units_source",
                    "last_scanned_at","scan_signature","ingest_status",
                    "last_scan_tier","consecutive_empty_ticks"]
            @test col in c
        end
        @test !("detector_distance_mm" in c)   # spec §10: no separate column
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: FAIL.

- [ ] **Step 3: Implement the migration**

The sentinel `MIGRATION_EXPERIMENTS_GEOMETRY`, its driver call (3rd, after `migrate_exposures_experiment_id!`), and a stub already exist from Task 1. **Replace the stub** with the real implementation:

```julia
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
            ("last_scanned_at", "TEXT"), ("scan_signature", "TEXT"),
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
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_ingestion_schema.jl
git commit -m "feat(db): add typed geometry + scan columns to experiments"
```

---

## Task 5: `samples` — name collapse + grouping columns

The careful part (spec §10, P0-2/P0-4). Order: add grouping columns → DROP INDEX `samples_unique_name` → DROP COLUMN `name` (manifest machine id) → RENAME `display_name` → `name`. Must run **after** `migrate_samples_naming!` and retire its duplicate-suffix pass (Task 5b). No new UNIQUE on the label.

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (`migrate_samples_name_collapse!`)
- Test: `packages/HimalayaUI/test/test_ingestion_schema.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "samples name collapse + grouping cols" begin
    path = fresh_db()
    with_db(path) do db
        c = cols(db, "samples")
        @test "name" in c
        @test !("display_name" in c)            # collapsed away
        for col in ["load_id","slot_index","grouping_source","name_source"]
            @test col in c
        end
        # label is NOT unique: two samples in one experiment may share a name
        DBInterface.execute(db, "PRAGMA foreign_keys=ON")
        eid = seed_experiment(db)   # honours NOT NULL name/path/data_dir/analysis_dir
        DBInterface.execute(db, "INSERT INTO samples (experiment_id, name) VALUES (?, 'HA85 (S01P15)')", [eid])
        @test_nowarn DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, name) VALUES (?, 'HA85 (S01P15)')", [eid])
        @test !any(i -> occursin("samples_unique_name", i), indexes(db, "samples"))
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: FAIL (either `display_name` still present, or the second insert throws on the surviving UNIQUE).

- [ ] **Step 3a: Stop `migrate_samples_naming!` from re-reverting the collapse (P0-1)**

`migrate_samples_naming!` (`db.jl:1431`) is **unsentineled** and runs on every `open_db`. Its guard is `if "display_name" in cols && !("label" in cols); return nothing; end` (`db.jl:1437`). After this task drops `display_name` and renames to `name`, that guard is **false**, so the body runs again and **re-adds `display_name`** (`db.jl:1442`) and **recreates `samples_unique_name`** (`db.jl:1492`) every startup. Add a post-redesign early-return immediately after the `isempty(cols) && return nothing` line (`db.jl:1436`):

```julia
    # Post-redesign shape: a single `name` label, no `display_name`/`label`. Nothing to rename.
    # Without this, the name-collapse migration (migrate_samples_name_collapse!) is silently
    # reverted on every open_db. (Phase-A P0-1.)
    if "name" in cols && !("display_name" in cols) && !("label" in cols)
        return nothing
    end
```

- [ ] **Step 3b: Implement the collapse migration**

Replace the `migrate_samples_name_collapse!` stub. (`migrate_schema!` calls `migrate_samples_naming!` before this; this runs after.)

```julia
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
end
```

> **SQLite version note:** `ALTER TABLE ... DROP COLUMN` / `RENAME COLUMN` require SQLite ≥ 3.35 / ≥ 3.25. Confirm the bundled `SQLite.jl` ships a new-enough engine (`DBInterface.execute(db, "SELECT sqlite_version()")` in the REPL). If not, fall back to the 12-step table-rebuild (create new table, copy, drop, rename) — but verify first; modern SQLite.jl is fine.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_ingestion_schema.jl
git commit -m "feat(db): collapse samples display_name->name + add grouping columns"
```

---

## Task 5b: Retire the duplicate-suffix pass in `migrate_samples_naming!`

`migrate_samples_naming!` (`db.jl:1463-1489`) appends a numeric suffix to disambiguate a UNIQUE machine name. Once labels may legitimately repeat (Task 5), that pass corrupts auto-generated labels. Neutralize it for fresh runs while leaving it harmless on already-migrated DBs.

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (`migrate_samples_naming!`, the dup-suffix loop ~`db.jl:1463-1489`)
- Test: `packages/HimalayaUI/test/test_ingestion_schema.jl`

- [ ] **Step 1: Write the test (duplicate labels survive the full migration)**

This is the meaningful check: a legacy DB whose two samples share a *display_name* (the future label) but have distinct unique machine *names* must, after `migrate_schema!`, end with two samples both named the shared label, **unsuffixed**, and with no surviving unique index. It guards both the dup-suffix removal (Step 3) and the P0-1 guard (Task 5 Step 3a).

```julia
@testset "duplicate labels survive naming+collapse" begin
    dir = mktempdir(); path = joinpath(dir, "duplabel.db")
    db = SQLite.DB(path)
    DBInterface.execute(db, "CREATE TABLE schema_migrations (name TEXT PRIMARY KEY, applied_at TEXT)")
    DBInterface.execute(db, "CREATE TABLE experiments (id INTEGER PRIMARY KEY, name TEXT, path TEXT NOT NULL, data_dir TEXT, analysis_dir TEXT)")
    DBInterface.execute(db, "CREATE TABLE samples (id INTEGER PRIMARY KEY, experiment_id INTEGER, name TEXT, display_name TEXT, notes TEXT)")
    DBInterface.execute(db, "CREATE UNIQUE INDEX samples_unique_name ON samples(experiment_id, name)")
    DBInterface.execute(db, "CREATE TABLE exposures (id INTEGER PRIMARY KEY, sample_id INTEGER, filename TEXT)")
    DBInterface.execute(db, "CREATE UNIQUE INDEX exposures_unique_filename ON exposures(sample_id, filename)")
    DBInterface.execute(db, "INSERT INTO experiments (id, name, path) VALUES (1,'e','/exp/e')")
    # distinct machine names (unique), SAME display_name (the future label)
    DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (1,1,'m1','HA85 (S01P15)')")
    DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (2,1,'m2','HA85 (S01P15)')")
    SQLite.close(db)

    dbm = SQLite.DB(path); HimalayaUI.migrate_schema!(dbm); SQLite.close(dbm)

    with_db(path) do db2
        names = String.(getproperty.(Tables.rowtable(DBInterface.execute(db2,
            "SELECT name FROM samples ORDER BY id")), :name))
        @test names == ["HA85 (S01P15)", "HA85 (S01P15)"]   # NOT suffixed to '...-2'
        @test !("display_name" in cols(db2, "samples"))
        @test !any(i -> occursin("samples_unique_name", i), indexes(db2, "samples"))
    end
end
```

- [ ] **Step 2: Run it**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: PASS once Task 5 (collapse + the P0-1 guard) is in place. **This is a safety-net + integration test, not a strict red-green:** the dup-suffix pass operated on the *unique machine `name`* (never the label, which lived in `display_name` pre-collapse), so it was effectively dead w.r.t. labels — Step 3 removes it for clarity, and this test ensures the whole naming→collapse pipeline never suffixes a legitimately-repeated label.

- [ ] **Step 3: Read and neutralize the pass**

Read `db.jl:1431-1498`. The dedup pass is the loop that detects duplicate `(experiment_id, display_name)` and rewrites later ones with a suffix. Replace that loop body with a no-op + comment (keep the function and its sentinel intact so the migration record is unchanged):

```julia
    # Duplicate-label disambiguation retired (ingestion redesign): sample labels may
    # legitimately repeat across loads (e.g. two `HA85 (S01P15)`); identity is (load_id, slot_index),
    # not the label. See migrate_samples_name_collapse! and spec §10.
```

- [ ] **Step 4: Run the full Phase-A test file**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: PASS (all testsets).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_ingestion_schema.jl
git commit -m "refactor(db): retire duplicate-suffix pass in migrate_samples_naming!"
```

---

## Task 6: `create_experiment!` signature

Add kwargs for the new geometry + scan columns so callers (the new `POST /api/experiments` in Phase C) can persist derived geometry. All optional with sensible defaults; existing callers keep working.

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (`create_experiment!`, ~`db.jl:1808-1826`)
- Test: `packages/HimalayaUI/test/test_ingestion_schema.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "create_experiment! geometry kwargs" begin
    path = fresh_db()
    with_db(path) do db
        id = HimalayaUI.create_experiment!(db; name="geo",
            path="/exp/geo", data_dir="/exp/geo/data", analysis_dir="/exp/geo/analysis",
            energy_kev=9.0, energy_kev_source="prp",
            flight_path_m=1.8095, flight_path_m_source="setup",
            beam_center_x=421.409, beam_center_y=836.946, beam_center_x_source="setup",
            pixel_size_um=172.0, q_units="A^-1")
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT flight_path_m, flight_path_m_source, beam_center_x FROM experiments WHERE id=?", [id])))
        @test row.flight_path_m ≈ 1.8095
        @test row.flight_path_m_source == "setup"
        @test row.beam_center_x ≈ 421.409
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: FAIL (unknown kwargs / columns not written).

- [ ] **Step 3: Extend `create_experiment!`**

The **live** signature (`db.jl:1808-1826`, verified) is a **strict superset** target: keep every existing kwarg/column (`name`, the NOT NULL `path`/`data_dir`/`analysis_dir`, `manifest_path`, `config`, `experiment_type`, `energy_kev`, `flight_path_m`) and ADD the new geometry/scan ones (all defaulting so current callers are unaffected):

```julia
function create_experiment!(db::SQLite.DB;
        name::Union{String,Nothing} = nothing,
        path::String, data_dir::String, analysis_dir::String,   # NOT NULL, keep required
        manifest_path::Union{String,Nothing} = nothing,
        config::Union{String,Nothing} = nothing,
        experiment_type::Union{String,Nothing} = nothing,
        energy_kev::Union{Float64,Nothing} = nothing,
        flight_path_m::Union{Float64,Nothing} = nothing,
        # --- new (Phase A) ---
        energy_kev_source = "default", flight_path_m_source = "default",
        beam_center_x = nothing, beam_center_x_source = "default",
        beam_center_y = nothing, beam_center_y_source = "default",
        pixel_size_um = nothing, pixel_size_um_source = "default",
        q_units = nothing, q_units_source = "default",
        last_scanned_at = nothing, scan_signature = nothing, ingest_status = "idle")
    result = DBInterface.execute(db, """
        INSERT INTO experiments
            (name, path, data_dir, analysis_dir, manifest_path, config, experiment_type,
             energy_kev, energy_kev_source, flight_path_m, flight_path_m_source,
             beam_center_x, beam_center_x_source, beam_center_y, beam_center_y_source,
             pixel_size_um, pixel_size_um_source, q_units, q_units_source,
             last_scanned_at, scan_signature, ingest_status)
        VALUES (?,?,?,?,?,?,?, ?,?,?,?, ?,?,?,?, ?,?,?,?, ?,?,?)
    """, [name, path, data_dir, analysis_dir, manifest_path, config, experiment_type,
          energy_kev, energy_kev_source, flight_path_m, flight_path_m_source,
          beam_center_x, beam_center_x_source, beam_center_y, beam_center_y_source,
          pixel_size_um, pixel_size_um_source, q_units, q_units_source,
          last_scanned_at, scan_signature, ingest_status])
    return Int(DBInterface.lastrowid(result))
end
```

> The test in Step 1 must pass the required `path`/`data_dir`/`analysis_dir` (use `seed_experiment`-style values). `manifest_path` stays accepted (its full removal is a later-phase cleanup) so existing `test_db.jl` callers (`:42,:145,:198,:253,:305`) keep working.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_ingestion_schema.jl
git commit -m "feat(db): create_experiment! accepts typed geometry + scan fields"
```

---

## Task 7: `create_sample!` signature (drop `display_name`, add grouping fields)

Must align with the post-collapse schema (Task 5): no `display_name` column. Spec §10/P0-4: the INSERT and the column drop are conceptually one commit; here we make the writer match the migrated schema.

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (`create_sample!`, ~`db.jl:1828-1837`)
- Test: `packages/HimalayaUI/test/test_ingestion_schema.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "create_sample! grouping fields, no display_name" begin
    path = fresh_db()
    with_db(path) do db
        eid = seed_experiment(db)   # honours NOT NULL name/path/data_dir/analysis_dir
        lid = DBInterface.lastrowid(DBInterface.execute(db,
            "INSERT INTO loads (experiment_id, load_index) VALUES (?, 1)", [eid]))
        sid = HimalayaUI.create_sample!(db; experiment_id=eid, name="HA85 (S01P15)",
            load_id=lid, slot_index=15, grouping_source="auto_position", name_source="auto")
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT name, load_id, slot_index, grouping_source, name_source FROM samples WHERE id=?", [sid])))
        @test row.name == "HA85 (S01P15)"
        @test row.load_id == lid
        @test row.slot_index == 15
        @test row.name_source == "auto"
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: FAIL (unknown kwargs, or INSERT names `display_name`).

- [ ] **Step 3: Rewrite `create_sample!`**

Read the current `create_sample!` (`db.jl:1828-1837` — note its INSERT currently names `(experiment_id, name, display_name, notes)`). Replace with:

```julia
function create_sample!(db::SQLite.DB; experiment_id::Integer, name::AbstractString,
        notes=nothing, load_id=nothing, slot_index=nothing,
        grouping_source="auto_position", name_source="auto")
    res = DBInterface.execute(db, """
        INSERT INTO samples (experiment_id, name, notes, load_id, slot_index, grouping_source, name_source)
        VALUES (?,?,?,?,?,?,?)
    """, [experiment_id, name, notes, load_id, slot_index, grouping_source, name_source])
    return Int(DBInterface.lastrowid(res))
end
```

> **Caller sweep (same commit):** grep `create_sample!(` across `packages/HimalayaUI/`. Update every caller that passes `display_name=` (e.g. `cli.jl:199`, `cli.jl` reingest path) to pass `name=` instead. Also update the unconditional clobber `UPDATE samples SET display_name=…` at `cli.jl:202-204` to `SET name=…, notes=…` **and** gate it `WHERE name_source != 'human'` (spec §4 never-clobber). Update any test fixtures passing `display_name:`. The CLI's grouping fields can be omitted (defaults apply) until Phase B rewrites that path.

- [ ] **Step 4: Run the full Phase-A file + a quick compile check of cli.jl**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Then: `julia --project=packages/HimalayaUI -e 'using HimalayaUI'`  (Expected: loads with no `display_name` errors.)
Expected: PASS + clean load.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/src/cli.jl packages/HimalayaUI/test/test_ingestion_schema.jl
git commit -m "feat(db): create_sample! single label + grouping fields; sweep callers"
```

---

## Task 8: `create_exposure!` signature (mandatory `experiment_id`, PRP fields, optional `sample_id`)

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (`create_exposure!`, ~`db.jl:1839-1851`)
- Test: `packages/HimalayaUI/test/test_ingestion_schema.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "create_exposure! experiment_id + prp fields, nullable sample_id" begin
    path = fresh_db()
    with_db(path) do db
        eid = seed_experiment(db)   # honours NOT NULL name/path/data_dir/analysis_dir
        # sample_id omitted (transient pre-group state must be allowed)
        xid = HimalayaUI.create_exposure!(db; experiment_id=eid, filename="HA_85_001.tif",
            prp_path="/d/HA_85_001.prp", timestamp="2026-04-12T10:02:00",
            exposure_time=2.0, horizontal_position=12.4, scan_id=2404, frame_no=1)
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT experiment_id, sample_id, horizontal_position, frame_no FROM exposures WHERE id=?", [xid])))
        @test row.experiment_id == eid
        @test row.sample_id === missing
        @test row.horizontal_position ≈ 12.4
        @test row.frame_no == 1
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: FAIL.

- [ ] **Step 3: Rewrite `create_exposure!`**

The **live** signature (`db.jl:1839-1851`, verified) keeps `filename`/`kind`/`selected`/`status`/`image_path` — those must survive (both `cli.jl` callers at `:98-101` and `:216-219` pass `image_path=`). The changes: make `sample_id` **optional** (transient pre-group state) and add `experiment_id` (required) + the PRP fields. Strict superset:

```julia
function create_exposure!(db::SQLite.DB;
        experiment_id::Int,                              # new, required
        sample_id::Union{Int,Nothing}      = nothing,    # was required Int; now optional
        filename::Union{String,Nothing}    = nothing,
        kind::String                       = "file",
        selected::Bool                     = false,
        status::Union{String,Nothing}      = nothing,
        image_path::Union{String,Nothing}  = nothing,
        # --- new (Phase A) ---
        prp_path = nothing, timestamp = nothing, exposure_time = nothing,
        horizontal_position = nothing, scan_id = nothing, frame_no = nothing,
        load_id = nothing, content_fingerprint = nothing)
    result = DBInterface.execute(db, """
        INSERT INTO exposures
            (experiment_id, sample_id, filename, kind, selected, status, image_path,
             prp_path, timestamp, exposure_time, horizontal_position, scan_id, frame_no,
             load_id, content_fingerprint)
        VALUES (?,?,?,?,?,?,?, ?,?,?,?,?,?,?,?)
    """, [experiment_id, sample_id, filename, kind, Int(selected), status, image_path,
          prp_path, timestamp, exposure_time, horizontal_position, scan_id, frame_no,
          load_id, content_fingerprint])
    return Int(DBInterface.lastrowid(result))
end
```

> The two `cli.jl` callers currently pass `sample_id=` (positionally-required before); they keep working since `sample_id` is still accepted. They must additionally pass `experiment_id=` once Phase B wires the scan, but for Phase A they compile (add `experiment_id=` at both sites in the Task's commit if they call `create_exposure!` directly during init/reingest — grep to confirm; the init path at `cli.jl:98-101` does).

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_ingestion_schema.jl
git commit -m "feat(db): create_exposure! takes experiment_id + PRP fields, nullable sample_id"
```

---

## Task 9: `analyze_exposure!` resolves experiment via `exposures.experiment_id`

Spec §9.1: today `analyze_exposure!` (`pipeline.jl:913-931`) resolves the experiment via `exposures → samples → experiments` (`pipeline.jl:920-922`), so an exposure with `sample_id = NULL` (the transient scan-before-group state) falls out and errors. Read `exposures.experiment_id` directly and resolve `analysis_dir` from it, dropping the `analysis_dir::String` positional param. Update both call sites.

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl` (`analyze_exposure!`, ~913-931)
- Modify: `packages/HimalayaUI/src/cli.jl` (`_analyze_experiment!` call site)
- Modify: `packages/HimalayaUI/src/routes_experiments.jl:102` (the `/analyze` call site)
- Test: `packages/HimalayaUI/test/test_ingestion_schema.jl`

- [ ] **Step 1: Write the failing test**

This test asserts the *resolution path* without running a full analysis (which needs real `.dat` files). It checks that `analyze_exposure!` can resolve the experiment for an exposure with `sample_id = NULL`, by stubbing the analysis at the resolution boundary. Concretely, assert the new internal resolver helper:

```julia
@testset "experiment resolution by exposures.experiment_id" begin
    path = fresh_db()
    with_db(path) do db
        eid = HimalayaUI.create_experiment!(db; name="e",
            analysis_dir="/tmp/analysis-xyz")
        xid = HimalayaUI.create_exposure!(db; experiment_id=eid, filename="f.tif")  # sample_id NULL
        # new helper resolves analysis_dir from the exposure's experiment_id directly
        @test HimalayaUI._resolve_analysis_dir(db, xid) == "/tmp/analysis-xyz"
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: FAIL (`_resolve_analysis_dir` undefined).

- [ ] **Step 3: Add the resolver + rewire `analyze_exposure!`**

In `pipeline.jl`, add the helper and change the JOIN. Read the current function first; the JOIN at ~920-922 looks like `FROM exposures e JOIN samples s ON s.id = e.sample_id JOIN experiments x ON x.id = s.experiment_id`. Replace the experiment-resolution with a direct read:

```julia
function _resolve_analysis_dir(db::SQLite.DB, exposure_id::Integer)
    row = first(Tables.rowtable(DBInterface.execute(db, """
        SELECT x.analysis_dir AS analysis_dir
          FROM exposures e JOIN experiments x ON x.id = e.experiment_id
         WHERE e.id = ?
    """, [Int(exposure_id)])))
    return row.analysis_dir
end
```

Then change `analyze_exposure!` to resolve `analysis_dir` internally when not supplied, dropping the positional. Keep a **backward-compatible** method so other callers compile:

```julia
# new primary method: resolves analysis_dir from the exposure's experiment
function analyze_exposure!(db::SQLite.DB, exposure_id::Integer; kwargs...)
    analysis_dir = _resolve_analysis_dir(db, exposure_id)
    return analyze_exposure!(db, exposure_id, analysis_dir; kwargs...)
end
```

Then fix the existing `analyze_exposure!(db, exposure_id, analysis_dir::String; ...)` body so it no longer requires `sample_id`. The current lookup (`pipeline.jl:919-924`) JOINs `exposures → samples` and returns **zero rows** for a pre-group exposure (`sample_id IS NULL`) → "exposure not found", defeating the whole point. Replace that lookup:

```julia
# BEFORE (pipeline.jl:919-924) — requires sample_id via the samples JOIN:
#   SELECT e.filename, s.experiment_id ... FROM exposures e
#     JOIN samples s ON s.id = e.sample_id WHERE e.id = ?
# AFTER — sample_id-independent:
row = first(Tables.rowtable(DBInterface.execute(db,
    "SELECT filename, experiment_id FROM exposures WHERE id = ?", [Int(exposure_id)])))
```
Use `row.filename` / `row.experiment_id` downstream as before. Confirm no other line in the body reads a `samples`-derived column.

- [ ] **Step 4: Update the two call sites to the new 2-arg form**

- `routes_experiments.jl:102`: `analyze_exposure!(db, Int(ex.id), String(analysis_dir))` → `analyze_exposure!(db, Int(ex.id))`. Then check whether `analysis_dir` is still used elsewhere in that handler; if not, drop its now-dead computation.
- `cli.jl:374` (`_analyze_experiment!`): the same `analyze_exposure!(db, …, analysis_dir)` call → drop the 3rd arg. Grep `analyze_exposure!(` to confirm these are the only two call sites.

- [ ] **Step 5: Run test + load check**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Then: `julia --project=packages/HimalayaUI -e 'using HimalayaUI'`
Expected: PASS + clean load.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/pipeline.jl packages/HimalayaUI/src/routes_experiments.jl packages/HimalayaUI/src/cli.jl packages/HimalayaUI/test/test_ingestion_schema.jl
git commit -m "feat(pipeline): resolve experiment via exposures.experiment_id (sample_id-independent)"
```

---

## Task 10: Full migration round-trip on a legacy fixture

Prove the migrations run cleanly **in order** against a database created at the *pre-redesign* schema (samples with `name`+`display_name`, exposures with `UNIQUE(sample_id, filename)`), not just a fresh one.

**Files:**
- Test: `packages/HimalayaUI/test/test_ingestion_schema.jl`

- [ ] **Step 1: Write the test**

```julia
@testset "legacy DB upgrades cleanly" begin
    dir = mktempdir(); path = joinpath(dir, "legacy.db")
    # Hand-build a minimal pre-redesign schema, then run the full migrator.
    db = SQLite.DB(path)
    DBInterface.execute(db, "CREATE TABLE schema_migrations (name TEXT PRIMARY KEY, applied_at TEXT)")
    DBInterface.execute(db, """CREATE TABLE experiments (id INTEGER PRIMARY KEY AUTOINCREMENT,
        name TEXT, path TEXT NOT NULL, data_dir TEXT, analysis_dir TEXT, description TEXT,
        energy_kev REAL, flight_path_m REAL, manifest_path TEXT)""")
    DBInterface.execute(db, """CREATE TABLE samples (id INTEGER PRIMARY KEY AUTOINCREMENT,
        experiment_id INTEGER, name TEXT, display_name TEXT, notes TEXT)""")
    DBInterface.execute(db, "CREATE UNIQUE INDEX samples_unique_name ON samples(experiment_id, name)")
    DBInterface.execute(db, """CREATE TABLE exposures (id INTEGER PRIMARY KEY AUTOINCREMENT,
        sample_id INTEGER, filename TEXT)""")
    DBInterface.execute(db, "CREATE UNIQUE INDEX exposures_unique_filename ON exposures(sample_id, filename)")
    # seed one experiment + sample + exposure so backfill has work
    DBInterface.execute(db, "INSERT INTO experiments (id, name, path) VALUES (1, 'legacy', '/exp/legacy')")
    DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (1, 1, 'JC001', 'JC C04')")
    DBInterface.execute(db, "INSERT INTO exposures (sample_id, filename) VALUES (1, 'f.tif')")
    SQLite.close(db)

    dbm = SQLite.DB(path); HimalayaUI.migrate_schema!(dbm); SQLite.close(dbm)   # takes a DB, not a path

    with_db(path) do db2
        @test "loads" in lowercase.(String.(getproperty.(Tables.rowtable(DBInterface.execute(db2,
            "SELECT name FROM sqlite_master WHERE type='table'")), :name)))
        @test "name" in cols(db2, "samples") && !("display_name" in cols(db2, "samples"))
        @test "experiment_id" in cols(db2, "exposures")
        # the backfill set exposures.experiment_id from the sample JOIN
        row = first(Tables.rowtable(DBInterface.execute(db2,
            "SELECT experiment_id FROM exposures WHERE filename='f.tif'")))
        @test row.experiment_id == 1
        # the surviving label is the old display_name value
        srow = first(Tables.rowtable(DBInterface.execute(db2, "SELECT name FROM samples WHERE id=1")))
        @test srow.name == "JC C04"
    end
end
```

> The legacy fixture builds a *minimal* pre-redesign `experiments`/`samples`/`exposures` with the period's UNIQUE indexes, seeds one row each, then runs the real `migrate_schema!(db)` to prove the upgrade path. `samples` here has both `name` (machine id `'JC001'`) and `display_name` (`'JC C04'`) so the collapse keeps `JC C04`.

- [ ] **Step 2: Run it**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_schema.jl`
Expected: PASS. If the collapse fails on `DROP COLUMN`, see the SQLite-version note in Task 5.

- [ ] **Step 3: Run the full backend suite once to catch regressions**

Run: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1; tail -40 /tmp/jl-test.out`
Expected: all green. Investigate any failure caused by the `create_*` signature or `display_name` changes (likely test fixtures that need the `name=` sweep from Task 7).

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/test/test_ingestion_schema.jl
git commit -m "test(db): full migration round-trip on a legacy fixture"
```

---

## Self-Review

**Spec coverage (§10 + §9.1):**
- loads table → Task 2 ✓
- exposures denorm + new cols + dedup swap + NULL fail-fast → Task 3 ✓
- experiments typed geometry + `*_source` + scan/scheduler cols, no `detector_distance_mm` → Task 4 ✓
- samples drop-index→drop-name→rename, grouping cols, no new UNIQUE, dup-suffix retired → Tasks 5, 5b ✓
- create_experiment!/create_sample!/create_exposure! signatures → Tasks 6, 7, 8 ✓ (P0-4 atomicity: Task 7 ships the writer + caller sweep in one commit alongside the Task 5 migration)
- analyze_exposure! sample_id-independent resolution + both call sites → Task 9 ✓
- legacy upgrade path → Task 10 ✓
- never-clobber gate on `name_source` → folded into Task 7's caller sweep ✓

**Deferred to later plans (not gaps):** `prp.jl`/`geometry.jl`/`grouping.jl`/`scan_and_group!` (Phase B); routes + SSE + scheduler (Phase C); event kinds (Phase D); frontend + the TS `display_name→name` sweep (Phase E). `manifest_path` is left accepted-but-ignored here; its full removal is an explicit later cleanup.

**Placeholder scan:** every code step shows complete code; every test step shows the assertion; every run step shows the command + expected result. The two "read the current function first" notes are accuracy guards (line numbers drift), not placeholders — the new code is given in full.

**Type/name consistency:** `create_sample!(; experiment_id, name, …)`, `create_exposure!(; experiment_id, filename, …)`, `create_experiment!(; name, …)`, `_resolve_analysis_dir(db, exposure_id)`, `analyze_exposure!(db, exposure_id)` are used identically in their definitions and tests. Sentinel names (`MIGRATION_LOADS_TABLE`, `MIGRATION_EXPOSURES_EXPERIMENT_ID`, `MIGRATION_SAMPLES_NAME_COLLAPSE`, `MIGRATION_EXPERIMENTS_GEOMETRY`) match between Tasks 1/4 and the test in Task 1.

**Review status (2026-06-18):** this plan was validated by a 2-reviewer workflow (himalaya-reviewer + plan-quality) which returned `not-ready` on the first draft and caught 5 P0s; all are now folded in: P0-1 the self-reverting `migrate_samples_naming!` guard (Task 5 Step 3a), P0-2/P0-3 the `create_experiment!`/`create_exposure!` superset rewrites that keep `path`/`data_dir`/`analysis_dir` and `kind`/`selected`/`status`/`image_path` (Tasks 6, 8 — verified against live `db.jl:1808-1851`), P0-4 the nonexistent `init_db!` + wrong `migrate_schema!` arity (Tasks 0, 10 — entry is `create_schema!`+`migrate_schema!(db)`), P0-5 the untested geometry sentinel (folded into Task 1). P1s folded: the `seed_experiment` helper for NOT NULL columns, the real `_migrated`/`_record_migration!`/`comparison_now_iso` facts (Task 1), the sample_id-independent `analyze_exposure!` body + exact call sites (Task 9), and the duplicate-label integration test (Task 5b). Remaining open: confirm the SQLite engine version supports `DROP COLUMN`/`RENAME COLUMN` (Task 5 note), and that the `cli.jl` `create_exposure!` callers pass `experiment_id` (Task 8 note).
