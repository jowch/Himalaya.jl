# Sample Naming Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor `samples` from `(label, name)` into `(name [stable identifier], display_name [editable])`, swap which manifest column flows to which DB field, add full pre-write manifest validation, hard-cutover the `experiment.toml` schema with a `migrate-toml` helper, and propagate the rename through the queue framework, contract tests, and live specs. Resolves issue #88; supersedes #83. Prerequisite for #89 (slug permalinks).

**Architecture:** Single PR, seven independently-green commits. Backend bottom-up (schema → validator → parser → CLI helper → routes/pipeline) then atomic frontend rename (type + helper + queue mutator + UI flow + every test fixture in one commit). The frontend rename is atomic because the frozen-shape contract test ties Sample type, payload, and cache row together — splitting them breaks intermediate-state contracts.

**Tech Stack:** Julia 1.x + SQLite.jl + Oxygen.jl 1.10 + Tables.jl (backend); React 18 + TypeScript strict + TanStack Query + Zustand + Vitest + Playwright (frontend); stdlib `Test` (Julia tests).

**Spec:** [`docs/superpowers/specs/2026-05-08-sample-naming-refactor-design.md`](../specs/2026-05-08-sample-naming-refactor-design.md). The spec is the authority on design intent; this plan is the build sequence. Each phase references spec §sections — open the spec when implementing a step.

**Worktree:** `.claude/worktrees/sample-naming-refactor` (branch `sample-naming-refactor`).

**Prod-DB safety snapshot:** `/opt/Himalaya.jl/data/himalaya.db.pre-issue88.20260509T012719Z` (139 samples / 678 exposures / 421 auto_peaks; integrity OK). Use **only** for the live-smoke step, never as a unit-test fixture.

**Conventions for this branch:**
- TDD: failing test → minimal implementation → green → commit. Each task ends with one commit.
- Julia tests: stdlib `Test` (`@testset`, `@test`, `@test_throws`). Internal helpers via `HimalayaUI.<name>`.
- Vitest tests: import from `../../src/...`, use `vi.fn()` for mocks, no Tailwind class assertions, `screen.getByText("X").closest("li")` not `document.querySelector`.
- The full HimalayaUI Julia suite takes 5–10 min — capture once with `> /tmp/jl-test.out 2>&1`, then `grep -E "Test Summary|did not pass|fail"`. Don't re-run with different greps.
- Frontend dev loop expects `tsc --noEmit + vite build` clean before PR — run `npm run build` after Phase 6.

---

## File structure

**Created:**
```
packages/HimalayaUI/src/
  validate.jl                 — ManifestViolation + validate_manifest (pure function)

packages/HimalayaUI/test/
  test_validate.jl            — pure-function tests for the validator
  test_migrate_toml.jl        — TOML rewrite happy path + edge cases

packages/HimalayaUI/frontend/src/lib/sample/
  displayName.ts              — single helper used by every sample-display call site

packages/HimalayaUI/frontend/test/lib/sample/
  displayName.test.ts         — helper unit tests
```

**Modified (backend):**
```
packages/HimalayaUI/src/
  db.jl                       — schema (samples DDL), migrate_samples_naming!, idempotent_responses purge
  config.jl                   — ExperimentConfig.col_label → col_display_name; reject [manifest].label
  manifest.jl                 — ManifestSample.label → display_name; column-meaning swap
  cli.jl                      — _cli_init_inner! extraction; create_sample! kwargs; --sample filter; cli_migrate_toml dispatch
  routes_samples.jl           — PATCH allowlist swap
  routes_experiments.jl       — PATCH allowlist drop :name; reingest 400/500 split
  routes_export.jl            — JSON :label key + CSV header rename
  comparisons.jl              — picker_samples SELECT column rename
  HimalayaUI.jl               — include "validate.jl"; export migrate-related symbols
configs/
  simple.toml                 — [manifest].label → name; [manifest].name → display_name
```

**Modified (frontend):**
```
packages/HimalayaUI/frontend/src/
  api.ts                      — Sample.label → display_name; updateSample patch type; updateExperiment cleanup
  lib/queue/mutators/trivial.ts          — updateSampleMutator full rewrite
  lib/queue/persistence.ts    — SCHEMA_VERSION 1 → 2
  lib/comparison/labels.ts    — JSDoc + inline expression
  components/SampleMetadataCard.tsx       — input now edits display_name; breadcrumb shows stable name
  components/NavModal.tsx     — primary/secondary lines + search-both-fields
  components/ComparisonPickerBody.tsx     — search-both-fields
  components/SamplePickerRow.tsx
  components/ExposureListRow.tsx
  components/PlotCard.tsx
  components/TitleButton.tsx
  components/MentionPicker.tsx
  components/MentionChip.tsx
  queries.ts
```

**Modified (tests — comprehensive sweep):**
```
packages/HimalayaUI/test/
  test_db.jl                              — migration regression tests
  test_config.jl                          — old-key rejection
  test_manifest.jl                        — column-meaning swap
  test_routes_samples.jl                  — PATCH allowlist
  test_routes_experiments.jl              — 400/500 split
  test_routes_export.jl                   — column rename
  test_pipeline.jl                        — _cli_init_inner! transaction; validate at the door
  test_idempotency_replay_invariant.jl    — every fixture's create_sample! call + line 207 PATCH body

packages/HimalayaUI/frontend/test/
  queue/cache-shape.test.ts               — SAMPLE_KEYS + 3 fixture sites
  queue/sseEventPayload.contract.test.ts  — fixture + assertion
  queue/rollbackSymmetry.test.ts          — shared fixture + updateSample case
  queue/authHeaders.test.ts               — request mock
  queue/mutatorOnSseWins.test.ts          — entire SSE-wins block
  queue/treats404AsSuccess.test.tsx
  queries-samples.test.tsx
  SamplePickerRow.test.tsx
  ComparisonPicker.test.tsx
  ComparisonPickerPanel.test.tsx
  ComparisonPickerBody.test.tsx
  ComparePageEdit.test.tsx
  ComparePage.test.tsx                    — legacy comments
  NavModal.test.tsx
  TitleButton.test.tsx
  smoke.test.tsx
  ExposureListRow.test.tsx                — flip "falls back to" test name

packages/HimalayaUI/frontend/e2e/
  inspect.spec.ts
  smoke.spec.ts
  compare.spec.ts
  bundle-b-paper-cuts.spec.ts
  figure-export.spec.ts
  multiplayer-replay-rerun.spec.ts
  live/sample-rename-preserves-fields.spec.ts   — semantically inverted; full rewrite

docs/
  experiment-config.md
  CLAUDE.md                               — gotchas for samples table
  superpowers/specs/2026-04-22-himalaya-web-app-design.md   — verify samples shape not stale
AGENTS.md (root + nested)
CHANGELOG.md (or release notes)
```

---

# Phase 1 — Backend: schema migration

Sentinel-guarded, atomic. Runs after `migrate_pk_to_autoincrement!` + `_fix_fk_references_after_autoincrement_migration!` so the FK heal sees a stable shape. Final step purges `idempotent_responses` to align the wire-format cache.

Spec reference: §3.

---

### Task 1.1: Update `create_schema!` so fresh DBs use canonical column order

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl:26-32`

- [ ] **Step 1: Edit the `samples` CREATE TABLE block**

```julia
CREATE TABLE IF NOT EXISTS samples (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id INTEGER REFERENCES experiments(id),
    name          TEXT,
    display_name  TEXT,
    notes         TEXT
);
```

- [ ] **Step 2: Add the UNIQUE INDEX immediately below the CREATE TABLE**

```julia
CREATE UNIQUE INDEX IF NOT EXISTS samples_unique_name ON samples(experiment_id, name);
```

- [ ] **Step 3: Verify backend module still loads**

Run: `julia --project=packages/HimalayaUI -e 'using HimalayaUI; println("HimalayaUI ok")'`
Expected: `HimalayaUI ok` (no schema errors)

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/src/db.jl
git commit -m "refactor(samples): canonical column order in fresh schema (issue #88)"
```

---

### Task 1.2: Write the failing migration test (synthetic legacy DB → migrated shape)

**Files:**
- Modify: `packages/HimalayaUI/test/test_db.jl` (append a new `@testset`)

- [ ] **Step 1: Write the failing test**

Add to `test_db.jl`:

```julia
@testset "migrate_samples_naming! — legacy (label, name) → (name, display_name)" begin
    mktempdir() do tmp
        dbpath = joinpath(tmp, "h.db")
        db     = SQLite.DB(dbpath)
        # Synthetic LEGACY shape (pre-rename), bypassing open_db's migrations.
        DBInterface.execute(db, "PRAGMA foreign_keys = ON")
        DBInterface.execute(db, """CREATE TABLE experiments (
            id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, path TEXT,
            data_dir TEXT, analysis_dir TEXT, manifest_path TEXT, config TEXT,
            experiment_type TEXT, energy_kev REAL, flight_path_m REAL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP)""")
        DBInterface.execute(db, """CREATE TABLE samples (
            id INTEGER PRIMARY KEY AUTOINCREMENT, experiment_id INTEGER REFERENCES experiments(id),
            label TEXT, name TEXT, notes TEXT)""")
        DBInterface.execute(db,
            "INSERT INTO experiments (id, name, path) VALUES (?, ?, ?)", [1, "exp", "/tmp"])
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, label, name, notes) VALUES (?, ?, ?, ?)",
            [1, "JC001", "DOPC + chol", "first run"])
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, label, name, notes) VALUES (?, ?, ?, ?)",
            [1, "", "fallback only", nothing])
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, label, name, notes) VALUES (?, ?, ?, ?)",
            [1, "JC002", "second run", nothing])

        HimalayaUI.migrate_samples_naming!(db)

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name, display_name, notes FROM samples ORDER BY id"))
        @test rows[1].name == "JC001";          @test rows[1].display_name == "DOPC + chol"
        @test rows[2].name == "fallback only";  @test rows[2].display_name == "fallback only"
        @test rows[3].name == "JC002";          @test rows[3].display_name == "second run"

        # Old `label` column is gone.
        cols = [r.name for r in Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info('samples')"))]
        @test "label" ∉ cols
        @test "display_name" ∈ cols

        # Unique index exists.
        idxs = [r.name for r in Tables.rowtable(DBInterface.execute(db,
            "PRAGMA index_list('samples')"))]
        @test "samples_unique_name" ∈ idxs
    end
end
```

- [ ] **Step 2: Run the test to confirm it fails**

```bash
julia --project=packages/HimalayaUI -e '
using HimalayaUI, Test, SQLite, DBInterface, Tables; include("packages/HimalayaUI/test/test_db.jl")' 2>&1 | tail -20
```

Expected: `UndefVarError: migrate_samples_naming! not defined` (or `MethodError`).

---

### Task 1.3: Implement `migrate_samples_naming!`

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (add helper next to `migrate_pk_to_autoincrement!`)
- Modify: `packages/HimalayaUI/src/db.jl` inside `migrate_schema!` to call the new helper

- [ ] **Step 1: Add the migration helper**

Append to `db.jl` (placement: after `_fix_fk_references_after_autoincrement_migration!`):

```julia
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
    if "display_name" in cols && !("label" in cols)
        return nothing  # already migrated
    end
    SQLite.transaction(db) do
        if !("display_name" in cols)
            try DBInterface.execute(db, "ALTER TABLE samples ADD COLUMN display_name TEXT")
            catch e; occursin("duplicate column", sprint(showerror, e)) || rethrow(); end
        end
        if "label" in cols
            DBInterface.execute(db, """UPDATE samples
                                          SET display_name = name,
                                              name         = COALESCE(NULLIF(label, ''), name)""")
            try DBInterface.execute(db, "ALTER TABLE samples DROP COLUMN label")
            catch e; occursin("no such column", sprint(showerror, e)) || rethrow(); end
        end
        # Duplicate suffix pass (oldest id keeps bare name).
        dups = Tables.rowtable(DBInterface.execute(db, """
            SELECT experiment_id, name FROM samples
            GROUP BY experiment_id, name HAVING COUNT(*) > 1"""))
        for d in dups
            ids = Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM samples WHERE experiment_id = ? AND name = ? ORDER BY id ASC",
                [d.experiment_id, d.name]))
            for (i, row) in enumerate(ids)
                i == 1 && continue  # oldest keeps the bare name
                new_name = "$(d.name)-$(i)"
                @warn "Renamed duplicate sample" experiment_id=d.experiment_id old=d.name new=new_name id=row.id
                DBInterface.execute(db, "UPDATE samples SET name = ? WHERE id = ?",
                    [new_name, row.id])
            end
        end
        DBInterface.execute(db,
            "CREATE UNIQUE INDEX IF NOT EXISTS samples_unique_name ON samples(experiment_id, name)")
        # Old idempotent_responses rows carry pre-rename payload shape; purge to
        # prevent stale-shape replays on retried client_op_id keys post-deploy.
        DBInterface.execute(db, "DELETE FROM idempotent_responses")
    end
    nothing
end
```

- [ ] **Step 2: Wire the helper into `migrate_schema!`**

Find the `migrate_schema!` function and insert the call **after** `migrate_pk_to_autoincrement!` and `_fix_fk_references_after_autoincrement_migration!`:

```julia
# Inside migrate_schema!, after the existing FK heal call:
migrate_samples_naming!(db)
```

- [ ] **Step 3: Run the test from Task 1.2**

```bash
julia --project=packages/HimalayaUI -e '
using HimalayaUI, Test, SQLite, DBInterface, Tables; include("packages/HimalayaUI/test/test_db.jl")' \
> /tmp/jl-mig.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-mig.out
```

Expected: PASS, all assertions green.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_db.jl
git commit -m "feat(db): migrate_samples_naming! with atomic transaction (issue #88)"
```

---

### Task 1.4: Add the duplicate-suffix regression test

**Files:**
- Modify: `packages/HimalayaUI/test/test_db.jl` (append)

- [ ] **Step 1: Write the failing test**

```julia
@testset "migrate_samples_naming! — duplicate names suffixed by ascending id" begin
    mktempdir() do tmp
        db = SQLite.DB(joinpath(tmp, "h.db"))
        DBInterface.execute(db, "PRAGMA foreign_keys = ON")
        DBInterface.execute(db, """CREATE TABLE experiments (
            id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, path TEXT,
            data_dir TEXT, analysis_dir TEXT, manifest_path TEXT, config TEXT,
            experiment_type TEXT, energy_kev REAL, flight_path_m REAL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP)""")
        DBInterface.execute(db, """CREATE TABLE samples (
            id INTEGER PRIMARY KEY AUTOINCREMENT, experiment_id INTEGER REFERENCES experiments(id),
            label TEXT, name TEXT, notes TEXT)""")
        DBInterface.execute(db, "INSERT INTO experiments (id, name) VALUES (1, 'exp')")
        # Three rows that collide post-COALESCE on (1, "JC001").
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, label, name) VALUES (?, ?, ?, ?)",
            [10, 1, "JC001", "v1"])
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, label, name) VALUES (?, ?, ?, ?)",
            [11, 1, "JC001", "v2"])
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, label, name) VALUES (?, ?, ?, ?)",
            [12, 1, "JC001", "v3"])

        HimalayaUI.migrate_samples_naming!(db)

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name FROM samples ORDER BY id"))
        @test rows[1].name == "JC001"
        @test rows[2].name == "JC001-2"
        @test rows[3].name == "JC001-3"
    end
end
```

- [ ] **Step 2: Run** — expect PASS (already implemented in Task 1.3).

```bash
julia --project=packages/HimalayaUI -e '
using HimalayaUI, Test, SQLite, DBInterface, Tables; include("packages/HimalayaUI/test/test_db.jl")' \
> /tmp/jl-mig.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-mig.out
```

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/test/test_db.jl
git commit -m "test(db): duplicate-name suffix regression for migrate_samples_naming! (issue #88)"
```

---

### Task 1.5: Add the idempotency regression test

**Files:**
- Modify: `packages/HimalayaUI/test/test_db.jl` (append)

- [ ] **Step 1: Write the failing test**

```julia
@testset "migrate_samples_naming! — idempotent on second run" begin
    mktempdir() do tmp
        db = SQLite.DB(joinpath(tmp, "h.db"))
        # Build canonical post-migration shape directly.
        DBInterface.execute(db, """CREATE TABLE samples (
            id INTEGER PRIMARY KEY AUTOINCREMENT, experiment_id INTEGER,
            name TEXT, display_name TEXT, notes TEXT)""")
        DBInterface.execute(db,
            "CREATE UNIQUE INDEX samples_unique_name ON samples(experiment_id, name)")
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (1, 1, 'JC001', 'DOPC')")
        # Migration should be a no-op (sentinel triggers).
        HimalayaUI.migrate_samples_naming!(db)
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name, display_name FROM samples"))
        @test length(rows) == 1
        @test rows[1].name == "JC001"
        @test rows[1].display_name == "DOPC"
    end
end
```

- [ ] **Step 2: Run** — expect PASS.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/test/test_db.jl
git commit -m "test(db): idempotency regression for migrate_samples_naming! (issue #88)"
```

---

### Task 1.6: Verify migration integrity against the prod-backup snapshot

This is a **smoke test, not a committed unit test** — verifies the migration handles the real 139-row prod data shape.

**Files:** none modified.

- [ ] **Step 1: Copy the backup to a scratch path so we never touch the snapshot**

```bash
cp /opt/Himalaya.jl/data/himalaya.db.pre-issue88.20260509T012719Z /tmp/himalaya-smoke.db
```

- [ ] **Step 2: Run open_db (which triggers migrate_schema!) against the copy**

```bash
julia --project=packages/HimalayaUI -e '
using HimalayaUI, SQLite, DBInterface, Tables
db = HimalayaUI.open_db("/tmp/himalaya-smoke.db")
ic  = first(Tables.rowtable(DBInterface.execute(db, "PRAGMA integrity_check"))).integrity_check
n   = first(Tables.rowtable(DBInterface.execute(db, "SELECT COUNT(*) AS c FROM samples"))).c
cols = [r.name for r in Tables.rowtable(DBInterface.execute(db, "PRAGMA table_info(samples)"))]
println("integrity=$ic samples=$n cols=$cols")
println("has_label=$("label" in cols) has_display_name=$("display_name" in cols)")
' 2>&1 | tail -15
```

Expected: `integrity=ok samples=139 cols=[id, experiment_id, name, notes, display_name]` (or similar with `display_name` appended), `has_label=false has_display_name=true`. Any `@warn` about renamed duplicates is OK.

---

# Phase 2 — Backend: validation module

Pure function, no DB. Lives in its own file. Spec reference: §4.

---

### Task 2.1: Write failing tests for `validate_manifest`

**Files:**
- Create: `packages/HimalayaUI/test/test_validate.jl`

- [ ] **Step 1: Create the test file**

```julia
using Test
using HimalayaUI: ManifestSample, ManifestViolation, validate_manifest

@testset "validate_manifest — empty name" begin
    samples = [ManifestSample("", "DOPC", "", "", ["JC001"])]
    vs = validate_manifest(samples)
    @test length(vs) == 1
    @test vs[1].kind == :empty_name
    @test vs[1].sample_index == 1
end

@testset "validate_manifest — bad name characters" begin
    samples = [ManifestSample("JC 001", "DOPC", "", "", ["JC001"])]  # space
    vs = validate_manifest(samples)
    @test length(vs) == 1
    @test vs[1].kind == :bad_name_chars
    @test occursin("JC 001", vs[1].detail)
end

@testset "validate_manifest — duplicate name" begin
    samples = [
        ManifestSample("JC001", "first",  "", "", ["JC001"]),
        ManifestSample("JC001", "second", "", "", ["JC002"]),
    ]
    vs = validate_manifest(samples)
    dup = filter(v -> v.kind == :duplicate_name, vs)
    @test length(dup) == 1
end

@testset "validate_manifest — duplicate filenames within sample" begin
    samples = [ManifestSample("JC001", "x", "", "", ["JC001", "JC001"])]
    vs = validate_manifest(samples)
    @test any(v -> v.kind == :duplicate_filename_in_sample, vs)
end

@testset "validate_manifest — overlapping filenames across samples" begin
    samples = [
        ManifestSample("JC001", "x", "", "", ["JC001", "JC002"]),
        ManifestSample("JC002", "y", "", "", ["JC002", "JC003"]),
    ]
    vs = validate_manifest(samples)
    @test any(v -> v.kind == :overlapping_filenames, vs)
end

@testset "validate_manifest — clean manifest produces no violations" begin
    samples = [
        ManifestSample("JC001", "first",  "", "", ["JC001"]),
        ManifestSample("JC002", "second", "", "", ["JC002"]),
    ]
    @test isempty(validate_manifest(samples))
end

@testset "validate_manifest — collects all violations, no fail-fast" begin
    samples = [
        ManifestSample("",       "x", "", "", ["JC001"]),         # empty
        ManifestSample("JC 002", "y", "", "", ["JC001", "JC001"]) # bad chars + dup filenames
    ]
    vs = validate_manifest(samples)
    @test length(vs) >= 3   # at least: empty_name, bad_name_chars, duplicate_filename_in_sample
end
```

- [ ] **Step 2: Wire it into the runtests.jl harness**

Append to `packages/HimalayaUI/test/runtests.jl`:

```julia
include("test_validate.jl")
```

- [ ] **Step 3: Run, expect failures (module not yet exported)**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI"; test_args=["test_validate"])' \
> /tmp/jl-validate.out 2>&1
tail -30 /tmp/jl-validate.out
```

Expected: `UndefVarError: validate_manifest not defined` (and ManifestViolation).

---

### Task 2.2: Implement `validate.jl`

**Files:**
- Create: `packages/HimalayaUI/src/validate.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl` (include the new file; export the public symbols)

- [ ] **Step 1: Write the module**

`packages/HimalayaUI/src/validate.jl`:

```julia
# validate.jl — pure-function manifest validation. No DB, no IO.

const _VALID_NAME_REGEX = r"^[A-Za-z0-9._-]+$"

"""
    ManifestViolation(kind, sample_index, sample_name, detail)

One validation problem found in a parsed manifest. `kind` is one of
`:empty_name`, `:bad_name_chars`, `:duplicate_name`,
`:duplicate_filename_in_sample`, `:overlapping_filenames`.
"""
struct ManifestViolation
    kind         ::Symbol
    sample_index ::Int        # 1-based row in manifest
    sample_name  ::String     # may be "" for :empty_name
    detail       ::String     # human-readable specifics
end

"""
    validate_manifest(samples) -> Vector{ManifestViolation}

Apply all five rules, collect every violation, return them in stable order
(rule 1, rule 2, rule 3, rule 4, rule 5). No fail-fast — operator sees every
fix needed in one pass.
"""
function validate_manifest(samples::Vector{ManifestSample})::Vector{ManifestViolation}
    out = ManifestViolation[]

    # Rule 1: name non-empty.
    for (i, s) in enumerate(samples)
        if isempty(s.name)
            push!(out, ManifestViolation(:empty_name, i, "", "row $i has an empty name"))
        end
    end

    # Rule 2: name matches convention.
    for (i, s) in enumerate(samples)
        if !isempty(s.name) && !occursin(_VALID_NAME_REGEX, s.name)
            push!(out, ManifestViolation(:bad_name_chars, i, s.name,
                "row $i name '$(s.name)' has characters outside [A-Za-z0-9._-]"))
        end
    end

    # Rule 3: name unique within manifest (by first-occurrence index).
    seen = Dict{String,Int}()
    for (i, s) in enumerate(samples)
        isempty(s.name) && continue
        if haskey(seen, s.name)
            push!(out, ManifestViolation(:duplicate_name, i, s.name,
                "row $i duplicates name '$(s.name)' first seen at row $(seen[s.name])"))
        else
            seen[s.name] = i
        end
    end

    # Rule 4: filenames within a sample have no internal duplicates.
    for (i, s) in enumerate(samples)
        seen_f = Dict{String,Int}()
        for f in s.filenames
            if haskey(seen_f, f)
                push!(out, ManifestViolation(:duplicate_filename_in_sample, i, s.name,
                    "row $i sample '$(s.name)' duplicates filename '$f'"))
            else
                seen_f[f] = i
            end
        end
    end

    # Rule 5: filenames don't overlap across samples.
    owner = Dict{String,Tuple{Int,String}}()
    for (i, s) in enumerate(samples)
        for f in s.filenames
            if haskey(owner, f) && owner[f][1] != i
                prev_i, prev_name = owner[f]
                push!(out, ManifestViolation(:overlapping_filenames, i, s.name,
                    "row $i sample '$(s.name)' filename '$f' also claimed by row $prev_i sample '$prev_name'"))
            else
                owner[f] = (i, s.name)
            end
        end
    end

    out
end
```

- [ ] **Step 2: Wire into the module entry**

Edit `packages/HimalayaUI/src/HimalayaUI.jl`. Find the existing `include("manifest.jl")` line; add immediately after:

```julia
include("validate.jl")
```

Add `ManifestViolation` and `validate_manifest` to the export list (or to `using` re-exports — match the file's existing pattern).

- [ ] **Step 3: Run the tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI"; test_args=["test_validate"])' \
> /tmp/jl-validate.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-validate.out
```

Expected: all `validate_manifest` testsets pass.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/src/validate.jl packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/test/test_validate.jl packages/HimalayaUI/test/runtests.jl
git commit -m "feat(validate): pure-function manifest validation (issue #88)"
```

---

# Phase 3 — Backend: parser cutover

Swaps `[manifest].label`/`name` → `name`/`display_name` in the TOML schema and the parser. Old key triggers a hard error pointing at `migrate-toml`. Spec reference: §5 (parser changes).

---

### Task 3.1: Update `simple.toml` template

**Files:**
- Modify: `packages/HimalayaUI/configs/simple.toml`

- [ ] **Step 1: Edit the `[manifest]` block**

```toml
[manifest]
delimiter      = "\t"
skip_rows      = 1
header_row     = 0
sample_id      = 1
name           = 2
display_name   = 3
filenames      = 9
notes_sample   = 10
notes_exposure = 11
```

- [ ] **Step 2: Commit**

```bash
git add packages/HimalayaUI/configs/simple.toml
git commit -m "refactor(config): swap [manifest].label/name to name/display_name in simple.toml (issue #88)"
```

---

### Task 3.2: Write failing tests for the parser cutover

**Files:**
- Modify: `packages/HimalayaUI/test/test_config.jl` (append)
- Modify: `packages/HimalayaUI/test/test_manifest.jl` (append)

- [ ] **Step 1: test_config.jl — old key rejected**

```julia
@testset "load_config rejects deprecated [manifest].label key" begin
    mktempdir() do tmp
        path = joinpath(tmp, "experiment.toml")
        write(path, """
[experiment]
name = "demo"
[manifest]
sample_id = 1
label     = 2
name      = 3
filenames = 9
[layout]
data_dir = "data"
analysis_dir = "analysis/automatic_analysis"
[files]
integration = "{name}.dat"
image       = "{name}.tiff"
""")
        @test_throws ErrorException HimalayaUI.load_config(path)
        try
            HimalayaUI.load_config(path)
        catch e
            @test occursin("deprecated key '[manifest].label'", sprint(showerror, e))
            @test occursin("migrate-toml", sprint(showerror, e))
        end
    end
end

@testset "load_config accepts new [manifest].name + display_name" begin
    mktempdir() do tmp
        path = joinpath(tmp, "experiment.toml")
        write(path, """
[experiment]
name = "demo"
[manifest]
sample_id    = 1
name         = 2
display_name = 3
filenames    = 9
[layout]
data_dir = "data"
analysis_dir = "analysis/automatic_analysis"
[files]
integration = "{name}.dat"
image       = "{name}.tiff"
""")
        cfg = HimalayaUI.load_config(path)
        @test cfg.col_name == 2
        @test cfg.col_display_name == 3
    end
end
```

- [ ] **Step 2: test_manifest.jl — column-meaning swap**

```julia
@testset "parse_manifest — column 2 is identifier (name), column 3 is display_name" begin
    cfg = HimalayaUI.load_builtin_config("simple")
    csv = """sample_id\tname\tdisplay_name\tcol4\tcol5\tcol6\tcol7\tcol8\tfilenames\tnotes_sample\tnotes_exposure
1\tJC001\tDOPC + chol\t\t\t\t\t\tJC001\t\t
2\tJC002\tPOPC\t\t\t\t\t\tJC002\t\t
"""
        # skip_rows=1 with header_row=0 means line 1 is skipped, parsing starts at line 2.
    samples = HimalayaUI.parse_manifest(cfg, IOBuffer(csv))
    @test length(samples) == 2
    @test samples[1].name == "JC001"
    @test samples[1].display_name == "DOPC + chol"
    @test samples[2].name == "JC002"
    @test samples[2].display_name == "POPC"
end
```

- [ ] **Step 3: Run, expect failures**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI"; test_args=["test_config","test_manifest"])' \
> /tmp/jl-parser.out 2>&1
tail -30 /tmp/jl-parser.out
```

Expected: `cfg.col_display_name` undefined; `samples[1].display_name` undefined.

---

### Task 3.3: Implement parser cutover

**Files:**
- Modify: `packages/HimalayaUI/src/config.jl`
- Modify: `packages/HimalayaUI/src/manifest.jl`

- [ ] **Step 1: `config.jl` — rename struct field + reject old key**

In the `ExperimentConfig` struct, rename `col_label` to `col_display_name`.

In `_build_config`, replace the `_coerce_col(get(mf, "label", 2), …)` line with a deprecation guard, and read the new keys:

```julia
if haskey(mf, "label")
    error("Deprecated key '[manifest].label' in experiment.toml. " *
          "The manifest column meanings were swapped: column 2 is now `name` " *
          "(stable identifier), column 3 is now `display_name` (user-friendly label). " *
          "Run `himalaya migrate-toml <experiment-dir>` to upgrade automatically.")
end
# ... build config ...
ExperimentConfig(
    # ... existing fields ...
    _coerce_col(get(mf, "sample_id",      1),  "manifest.sample_id"),
    _coerce_col(get(mf, "name",           2),  "manifest.name"),
    _coerce_col(get(mf, "display_name",   3),  "manifest.display_name"),
    _coerce_col(get(mf, "filenames",      9),  "manifest.filenames"),
    # ... rest ...
)
```

In `config_to_toml`, emit `name` and `display_name` keys instead of `label`/`name`:

```julia
"manifest" => Dict(
    # ... existing keys ...
    "sample_id"      => col_value(cfg.col_sample_id),
    "name"           => col_value(cfg.col_name),
    "display_name"   => col_value(cfg.col_display_name),
    "filenames"      => col_value(cfg.col_filenames),
    # ... rest ...
),
```

- [ ] **Step 2: `manifest.jl` — rename ManifestSample field + swap parser reads**

```julia
struct ManifestSample
    name           ::String   # was: label
    display_name   ::String   # was: name
    notes_sample   ::String
    notes_exposure ::String
    filenames      ::Vector{String}
end
```

In `parse_manifest`, rename `sym_label` → `sym_name_col` and `sym_name` → `sym_display_name_col`:

```julia
sym_id              = resolve_col(cfg.col_sample_id)
sym_name_col        = resolve_col(cfg.col_name)
sym_display_name_col = resolve_col(cfg.col_display_name)
sym_filenames       = resolve_col(cfg.col_filenames)
sym_notes_s         = resolve_col(cfg.col_notes_sample)
sym_notes_e         = resolve_col(cfg.col_notes_exposure)

# ... in the row loop:
push!(samples, ManifestSample(
    safe_get(row, sym_name_col),         # was: safe_get(row, sym_label)  (column 2 → name)
    safe_get(row, sym_display_name_col), # was: safe_get(row, sym_name)   (column 3 → display_name)
    safe_get(row, sym_notes_s),
    safe_get(row, sym_notes_e),
    expand_filename_field(filenames_str),
))
```

- [ ] **Step 3: Run the tests from Task 3.2**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI"; test_args=["test_config","test_manifest"])' \
> /tmp/jl-parser.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-parser.out
```

Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/src/config.jl packages/HimalayaUI/src/manifest.jl \
        packages/HimalayaUI/test/test_config.jl packages/HimalayaUI/test/test_manifest.jl
git commit -m "refactor(parser): swap [manifest] column meanings; reject old label key (issue #88)"
```

---

# Phase 4 — Backend: `migrate-toml` CLI

Section-aware regex-anchored substitution so the "axis label units" comment in `[beamline]` is not misfired. Spec reference: §5 (migrate-toml helper).

---

### Task 4.1: Write failing tests for `cli_migrate_toml`

**Files:**
- Create: `packages/HimalayaUI/test/test_migrate_toml.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl` (include the new file)

- [ ] **Step 1: Create the test file**

```julia
using Test
using HimalayaUI: cli_migrate_toml

const _OLD_TOML = """
[experiment]
name        = "demo"

[beamline]
# axis label units (ASCII)
q_units = "A-1"

[manifest]
delimiter = "\\t"
skip_rows = 1
sample_id = 1
label     = 2  # column for sample identifier
name      = 3
filenames = 9

[files]
integration = "{name}.dat"
image       = "{name}.tiff"
"""

const _NEW_TOML = """
[experiment]
name        = "demo"

[beamline]
# axis label units (ASCII)
q_units = "A-1"

[manifest]
delimiter = "\\t"
skip_rows = 1
sample_id = 1
name      = 2  # column for sample identifier
display_name = 3
filenames = 9

[files]
integration = "{name}.dat"
image       = "{name}.tiff"
"""

@testset "migrate-toml — happy path" begin
    mktempdir() do dir
        path = joinpath(dir, "experiment.toml")
        write(path, _OLD_TOML)
        cli_migrate_toml([dir])
        out = read(path, String)
        @test occursin(r"^name\s*=\s*2\s+#"m, out)        # comment preserved
        @test occursin(r"^display_name\s*=\s*3"m, out)
        @test !occursin(r"^\s*label\s*="m, out)            # old key gone
        # The beamline comment "axis label units" must not be misfired:
        @test occursin("# axis label units", out)
    end
end

@testset "migrate-toml — idempotent on already-migrated file" begin
    mktempdir() do dir
        path = joinpath(dir, "experiment.toml")
        write(path, _NEW_TOML)
        @test_logs (:info, r"already migrated") cli_migrate_toml([dir])
        @test read(path, String) == _NEW_TOML  # unchanged
    end
end

@testset "migrate-toml — errors on file with both old and new keys" begin
    mktempdir() do dir
        path = joinpath(dir, "experiment.toml")
        write(path, """
[manifest]
sample_id    = 1
label        = 2
display_name = 3
filenames    = 9
""")
        @test_throws ErrorException cli_migrate_toml([dir])
    end
end

@testset "migrate-toml — errors on missing experiment.toml" begin
    mktempdir() do dir
        @test_throws ErrorException cli_migrate_toml([dir])
    end
end
```

Append `include("test_migrate_toml.jl")` to `runtests.jl`.

- [ ] **Step 2: Run, expect failure (function undefined)**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI"; test_args=["test_migrate_toml"])' \
> /tmp/jl-migt.out 2>&1
tail -20 /tmp/jl-migt.out
```

Expected: `UndefVarError: cli_migrate_toml not defined`.

---

### Task 4.2: Implement `cli_migrate_toml`

**Files:**
- Modify: `packages/HimalayaUI/src/cli.jl` (add the function and dispatch in `main`)

- [ ] **Step 1: Add the function**

Append to `cli.jl`:

```julia
"""
    cli_migrate_toml(args)

Rewrite `<dir>/experiment.toml` from the legacy `[manifest].label/name` shape
to the canonical `[manifest].name/display_name` shape. Section-aware: only
substitutes inside the `[manifest]` block so the "axis label units" comment
in `[beamline]` is not misfired. Idempotent. Atomic file write.
"""
function cli_migrate_toml(args::Vector{<:AbstractString})
    isempty(args) && error("Usage: himalaya migrate-toml <experiment-dir>")
    dir  = args[1]
    path = joinpath(dir, "experiment.toml")
    isfile(path) || error("experiment.toml not found in $dir")

    lines = readlines(path; keep=true)
    section = ""
    has_label = false
    has_display_name = false
    has_old_name = false  # `name` inside [manifest] when display_name is absent → legacy-shape `name`

    # First pass: classify what's in [manifest].
    for line in lines
        m = match(r"^\s*\[([A-Za-z0-9_]+)\]\s*$", line)
        if m !== nothing
            section = m.captures[1]; continue
        end
        if section == "manifest"
            occursin(r"^\s*label\s*="m,        line) && (has_label = true)
            occursin(r"^\s*display_name\s*="m, line) && (has_display_name = true)
            occursin(r"^\s*name\s*="m,         line) && (has_old_name = true)
        end
    end

    if has_display_name && !has_label
        @info "experiment.toml at $path already migrated"
        return nothing
    end
    if has_display_name && has_label
        error("experiment.toml at $path has both `label` and `display_name` in [manifest]; " *
              "manual edit needed")
    end
    if !has_label
        error("experiment.toml at $path has no `[manifest].label` to migrate")
    end

    # Second pass: rewrite. Only inside [manifest], and only on lines that look like assignments.
    section = ""
    out_lines = String[]
    for line in lines
        m = match(r"^\s*\[([A-Za-z0-9_]+)\]\s*$", line)
        if m !== nothing
            section = m.captures[1]
            push!(out_lines, line); continue
        end
        if section == "manifest"
            # Per-line state machine: a line is EITHER a label= rewrite OR a name= rewrite,
            # never both — so `name` never catches the line we just rewrote from `label`.
            if (m2 = match(r"^(\s*)label(\s*=\s*\S+)(.*)$", line)) !== nothing
                push!(out_lines, m2.captures[1] * "name" * m2.captures[2] * m2.captures[3] * "\n")
            elseif (m3 = match(r"^(\s*)name(\s*=\s*\S+)(.*)$", line)) !== nothing
                push!(out_lines, m3.captures[1] * "display_name" * m3.captures[2] * m3.captures[3] * "\n")
            else
                push!(out_lines, line)
            end
        else
            push!(out_lines, line)
        end
    end

    # Atomic write.
    tmp = path * ".tmp"
    open(tmp, "w") do io
        for l in out_lines; print(io, l); end
    end
    mv(tmp, path; force=true)
    @info "Migrated $path"
    nothing
end
```

- [ ] **Step 2: Wire dispatch in `main`**

Find the `main(ARGS)` function in `cli.jl`. Add a branch for `"migrate-toml"` next to the other subcommands:

```julia
elseif cmd == "migrate-toml"
    cli_migrate_toml(rest)
```

- [ ] **Step 3: Run the tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI"; test_args=["test_migrate_toml"])' \
> /tmp/jl-migt.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-migt.out
```

Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/src/cli.jl packages/HimalayaUI/test/test_migrate_toml.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "feat(cli): himalaya migrate-toml helper for issue #88"
```

---

# Phase 5 — Backend: route + pipeline cutover

The big backend commit. `create_sample!` kwargs swap, all 8 `cli.jl` sites, route allowlists, JSON+CSV export rename, picker_samples SELECT, `_cli_init_inner!` extraction with transaction wrap, validate at the door. Spec reference: §6, §4 (validation surface).

This phase is broken into smaller tasks because the surface area is broad.

---

### Task 5.1: Update `create_sample!` kwargs and tighten downstream callers (compile-only)

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl:1119-1128`
- Modify: `packages/HimalayaUI/src/cli.jl` (lines 64, 164, 170 — the three current `label = …` keyword args)

- [ ] **Step 1: Edit `create_sample!`**

```julia
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
```

- [ ] **Step 2: Edit `cli.jl:64` (init insert)**

```julia
s_id = create_sample!(db;
    experiment_id = exp_id,
    name          = ms.name,
    display_name  = ms.display_name,
    notes         = ms.notes_sample)
```

- [ ] **Step 3: Edit `cli.jl:164` (reingest insert)**

Same pattern.

- [ ] **Step 4: Edit `cli.jl:170` (reingest UPDATE)**

```julia
DBInterface.execute(db,
    "UPDATE samples SET display_name = ?, notes = ? WHERE id = ?",
    [ms.display_name, ms.notes_sample, existing[1].id])
```

- [ ] **Step 5: Module compiles**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI; println("ok")' 2>&1 | tail -3
```

Expected: `ok`.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/src/cli.jl
git commit -m "refactor(samples): create_sample! kwargs swap label→display_name (issue #88)"
```

---

### Task 5.2: Sweep remaining `sample.label` references in cli.jl (filter, log, help)

**Files:**
- Modify: `packages/HimalayaUI/src/cli.jl:311, 335, 338, 357, 374, 383`

- [ ] **Step 1: Replace `sm.label` with `sm.name` in the `--sample` filter**

`cli.jl:311`:

```julia
sample_filter !== nothing && filter!(sm -> sm.name == sample_filter, samples)
```

- [ ] **Step 2: Replace `sample.label` in println log lines (cli.jl:335,338)**

```julia
println("  Skipping $(sample.name) / $(exp_row.filename) (rejected)")
print("  Analyzing $(sample.name) / $(exp_row.filename) ... ")
```

- [ ] **Step 3: Update `--sample/-s` flag help text (cli.jl:357 + 374)**

Change `"sample label"` → `"sample name (stable identifier)"`.

- [ ] **Step 4: Replace `sm.label` in cli.jl:383**

```julia
findfirst(sm -> sm.name == p[:sample], samples)
```

- [ ] **Step 5: Module compiles, no `sm.label` remains in cli.jl**

```bash
grep -n "sm\.label\|sample\.label" packages/HimalayaUI/src/cli.jl || echo "clean"
```

Expected: `clean`.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/cli.jl
git commit -m "refactor(cli): rename sm.label → sm.name in filter/log/help (issue #88)"
```

---

### Task 5.3: Sweep `samples.label` from non-CLI backend SQL + log strings

**Files:**
- Modify: `packages/HimalayaUI/src/comparisons.jl:507` (picker_samples explicit SELECT)
- Modify: `packages/HimalayaUI/src/routes_experiments.jl:95` (error log string)

- [ ] **Step 1: `comparisons.jl:507`**

Find the `SELECT id, experiment_id, name, label, notes FROM samples …` query. Change to:

```julia
SELECT id, experiment_id, name, display_name, notes FROM samples …
```

(Verify which row fields the picker projects to JSON; if it spreads the row into a dict, the projected key changes from `label` to `display_name` automatically. If not, also update the projection map.)

- [ ] **Step 2: `routes_experiments.jl:95`**

Replace `"$(sm.label)/$(ex.filename): …"` with `"$(sm.name)/$(ex.filename): …"`.

- [ ] **Step 3: Verify**

```bash
grep -rn "samples\.label\|sm\.label\|sample\.label" packages/HimalayaUI/src/ \
  | grep -v "label_color\|label_override\|axis label\|colour label" || echo "clean"
```

Expected: `clean` (or only the genuinely-unrelated matches like `colour label`, `label_override` for comparison-member overrides, etc).

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/src/comparisons.jl packages/HimalayaUI/src/routes_experiments.jl
git commit -m "refactor(backend): sweep sample.label references in SQL + log strings (issue #88)"
```

---

### Task 5.4: Update `routes_export.jl` JSON + CSV column rename

**Files:**
- Modify: `packages/HimalayaUI/src/routes_export.jl` (sites at lines 56, 67, 73, 81 per spec)
- Modify: `packages/HimalayaUI/test/test_routes_export.jl`

- [ ] **Step 1: Test first (failing). Update test_routes_export.jl**

Find the test that builds samples with `label="D1", name="UX1"` and asserts CSV header `"sample_label,sample_name,…"`. Rewrite:

```julia
# Setup:
HimalayaUI.create_sample!(db; experiment_id=exp_id,
    name="D1", display_name="UX1")
# ...
# Assertions:
@test s.name == "D1"
@test s.display_name == "UX1"
# CSV header:
@test startswith(csv_body, "sample_name,sample_display_name,exposure_filename,phases")
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI"; test_args=["test_routes_export"])' \
> /tmp/jl-export.out 2>&1
tail -20 /tmp/jl-export.out
```

Expected: failures on the column header + sample.label assertions.

- [ ] **Step 3: Implement — update routes_export.jl**

In the JSON dict construction, replace `:label => sm.label` with `:display_name => sm.display_name`.

In the CSV header definition, replace `"sample_label"` with `"sample_display_name"` (and update the per-row `sm[:label]` access to `sm[:display_name]`).

- [ ] **Step 4: Re-run test, expect pass**

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_export.jl packages/HimalayaUI/test/test_routes_export.jl
git commit -m "refactor(export): JSON + CSV column rename label→display_name (issue #88)"
```

---

### Task 5.5: PATCH `/api/samples/:id` allowlist swap + display_name trim

**Files:**
- Modify: `packages/HimalayaUI/src/routes_samples.jl:24`
- Modify: `packages/HimalayaUI/test/test_routes_samples.jl`

- [ ] **Step 1: Failing test — old `name` key returns 400; new `display_name` accepted + trimmed**

In `test_routes_samples.jl`, add:

```julia
@testset "PATCH /api/samples/:id rejects name (now immutable), accepts display_name" begin
    # ... existing fixture setup creating a sample ...
    # Old shape returns 400.
    r = HTTP.request("PATCH", "$base/api/samples/$sid",
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:name => "renamed")); status_exception=false)
    @test r.status == 400
    # New shape accepted, leading/trailing whitespace trimmed.
    r2 = HTTP.request("PATCH", "$base/api/samples/$sid",
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:display_name => "  spaced  ")))
    @test r2.status == 200
    body = JSON3.read(String(r2.body))
    @test body[:display_name] == "spaced"
end
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement — swap allowlist, add trim**

In `routes_samples.jl:24` PATCH handler:

```julia
fields, vals = Symbol[], Any[]
for k in (:display_name, :notes)
    if haskey(body, k)
        v = body[k]
        # Trim leading/trailing whitespace on display_name only.
        v isa AbstractString && k === :display_name && (v = strip(String(v)))
        push!(fields, k); push!(vals, v)
    end
end
```

- [ ] **Step 4: Re-run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_samples.jl packages/HimalayaUI/test/test_routes_samples.jl
git commit -m "refactor(routes): PATCH /api/samples allowlist display_name+notes (issue #88)"
```

---

### Task 5.6: PATCH `/api/experiments/:id` drops `:name`; reingest 400/500 split

**Files:**
- Modify: `packages/HimalayaUI/src/routes_experiments.jl:42-75` (PATCH allowlist)
- Modify: `packages/HimalayaUI/src/routes_experiments.jl:108-131` (reingest catch)
- Create: `ManifestValidationError <: Exception` (probably in `validate.jl` or a new tiny `errors.jl`)
- Modify: `packages/HimalayaUI/src/cli.jl::_reingest_inner!` to throw `ManifestValidationError` after running `validate_manifest`
- Modify: `packages/HimalayaUI/test/test_routes_experiments.jl`

- [ ] **Step 1: Failing test for the PATCH drop**

```julia
@testset "PATCH /api/experiments/:id no longer accepts name" begin
    # ... fixture ...
    r = HTTP.request("PATCH", "$base/api/experiments/$eid",
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:name => "newname")); status_exception=false)
    @test r.status == 400
end

@testset "POST /api/experiments/:id/reingest returns 400 on validation error" begin
    # ... build an experiment whose manifest has duplicates ...
    r = HTTP.request("POST", "$base/api/experiments/$eid/reingest", []; status_exception=false)
    @test r.status == 400
    body = JSON3.read(String(r.body))
    @test body[:error] == "manifest_invalid"
    @test !isempty(body[:violations])
end
```

- [ ] **Step 2: Run, expect failures**

- [ ] **Step 3: Add `ManifestValidationError` to `validate.jl`**

```julia
"""
    ManifestValidationError(violations)

Raised when `validate_manifest` returns non-empty violations and the caller
rejects the manifest. The HTTP layer maps this to 400 + structured JSON;
the CLI layer prints all violations and exits non-zero.
"""
struct ManifestValidationError <: Exception
    violations::Vector{ManifestViolation}
end
function Base.showerror(io::IO, e::ManifestValidationError)
    print(io, "Manifest invalid ($(length(e.violations)) violations):\n")
    for v in e.violations
        print(io, "  [$(v.kind)] $(v.detail)\n")
    end
end
```

Add to exports.

- [ ] **Step 4: Wire into `_reingest_inner!` (cli.jl:124)**

After `samples = parse_manifest(cfg, manifest_path)`, immediately:

```julia
violations = validate_manifest(samples)
isempty(violations) || throw(ManifestValidationError(violations))
```

- [ ] **Step 5: Update `routes_experiments.jl` PATCH (drop `:name`)**

```julia
fields = Symbol[]
vals   = Any[]
# (No mutable fields right now — route is defensive surface for the future.)
if isempty(fields)
    return HTTP.Response(400,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "no updatable fields provided")))
end
```

(Either remove the `for k in (:name,)` block entirely, or keep an empty allowlist.)

- [ ] **Step 6: Update `routes_experiments.jl` reingest (400/500 split)**

```julia
@post "/api/experiments/{id}/reingest" function(req::HTTP.Request, id::Int)
    db   = current_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT path FROM experiments WHERE id = ?", [id]))
    isempty(rows) && return HTTP.Response(404, ...)
    exp_path = String(rows[1].path)
    try
        res = reingest!(db, id, exp_path)
        log_action!(db, req; action = "reingest",
            entity_type = "experiment", entity_id = id)
        return HTTP.Response(200, ..., JSON3.write(Dict(...)))
    catch e
        if e isa ManifestValidationError
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "manifest_invalid",
                                 :violations => [Dict(:kind=>v.kind, :sample_index=>v.sample_index,
                                                      :sample_name=>v.sample_name, :detail=>v.detail)
                                                 for v in e.violations])))
        end
        return HTTP.Response(500,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => sprint(showerror, e))))
    end
end
```

- [ ] **Step 7: Re-run, expect PASS**

- [ ] **Step 8: Commit**

```bash
git add packages/HimalayaUI/src/routes_experiments.jl packages/HimalayaUI/src/cli.jl \
        packages/HimalayaUI/src/validate.jl packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/test/test_routes_experiments.jl
git commit -m "refactor(routes): drop name from experiments PATCH; 400/500 split on reingest (issue #88)"
```

---

### Task 5.7: Extract `_cli_init_inner!` and wrap in transaction

**Files:**
- Modify: `packages/HimalayaUI/src/cli.jl` (extract the manifest-driven section of `cli_init_with_db!`)
- Modify: `packages/HimalayaUI/test/test_pipeline.jl`

- [ ] **Step 1: Failing test — partial init rolls back fully on validation failure**

```julia
@testset "cli_init_with_db! is atomic on validation failure" begin
    mktempdir() do tmp
        # Create an experiment dir whose manifest has a duplicate name.
        # ... build TOML + CSV ...
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        @test_throws HimalayaUI.ManifestValidationError HimalayaUI.cli_init_with_db!(db, exp_dir)
        # No experiment row leaked.
        rows = Tables.rowtable(DBInterface.execute(db, "SELECT * FROM experiments"))
        @test isempty(rows)
    end
end
```

- [ ] **Step 2: Run, expect failure (the experiment row leaks because there's no transaction wrap yet)**

- [ ] **Step 3: Extract `_cli_init_inner!`**

Refactor `cli_init_with_db!` so the manifest-driven section becomes:

```julia
function cli_init_with_db!(db::SQLite.DB, exp_dir::AbstractString; analyze::Bool=true)
    exp_dir = abspath(exp_dir)
    SQLite.transaction(db) do
        _cli_init_inner!(db, exp_dir)
    end
    if analyze
        # auto-analyze stays OUTSIDE the transaction — a crash mid-analyze must not
        # roll back the experiment registration.
        ...
    end
end

function _cli_init_inner!(db::SQLite.DB, exp_dir::String)
    # Existing experiment row creation + manifest parse + validate + sample/exposure inserts.
    # Includes: load_config, create_experiment!, parse_manifest, validate_manifest (throws),
    #           per-sample create_sample! + per-exposure create_exposure!.
end
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Run the existing read-only-experiment-dir snapshot regression**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI"; test_args=["test_pipeline"])' \
> /tmp/jl-pipeline.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-pipeline.out
```

Expected: PASS. The snapshot test ("cli_init_with_db! does not write to experiment directory") must stay green.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/cli.jl packages/HimalayaUI/test/test_pipeline.jl
git commit -m "refactor(cli): extract _cli_init_inner! with transaction wrap (closes #83 footgun #4)"
```

---

### Task 5.8: Update test_idempotency_replay_invariant.jl fixtures + PATCH body

**Files:**
- Modify: `packages/HimalayaUI/test/test_idempotency_replay_invariant.jl` (lines 78, 118, 166, 199, 207, 240)

- [ ] **Step 1: Replace every `create_sample!(db; experiment_id=exp_id, label="D1")` call**

Sweep:

```bash
grep -n "create_sample!" packages/HimalayaUI/test/test_idempotency_replay_invariant.jl
```

For each match, change `label="D1"` to `name="D1"`. (Drop the `label` kwarg; keep / add `name`.)

- [ ] **Step 2: Update the PATCH body at line 207**

```julia
body_json = JSON3.write(Dict(:display_name => "renamed"))
```

- [ ] **Step 3: Run the test file**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI"; test_args=["test_idempotency_replay_invariant"])' \
> /tmp/jl-idemp.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-idemp.out
```

Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/test/test_idempotency_replay_invariant.jl
git commit -m "test(idempotency): swap fixtures + PATCH body to display_name (issue #88)"
```

---

### Task 5.9: Run the full backend Julia suite once

**Files:** none modified. This is a checkpoint.

- [ ] **Step 1: Capture full suite output**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' \
> /tmp/jl-test.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-test.out
tail -50 /tmp/jl-test.out
```

Expected: all testsets pass. If any fail, fix before moving to Phase 6.

---

# Phase 6 — Frontend: type + helper + atomic rename

ONE commit. Sample type rename, helper creation, mutator rewrite, SampleMetadataCard rewrite, persistence schema bump, search both fields, every test fixture in lockstep. Spec reference: §7.

This phase has many files but ONE commit because the frozen-shape contract test (`cache-shape.test.ts`) ties Sample type, payload shape, and cache row together.

---

### Task 6.1: Failing tests for the displayName helper

**Files:**
- Create: `packages/HimalayaUI/frontend/test/lib/sample/displayName.test.ts`

- [ ] **Step 1: Write the test file**

```ts
import { describe, it, expect } from "vitest";
import { sampleDisplayName } from "../../../src/lib/sample/displayName";

describe("sampleDisplayName", () => {
  it("prefers display_name when set", () => {
    expect(sampleDisplayName({ id: 1, name: "JC001", display_name: "DOPC" })).toBe("DOPC");
  });
  it("falls back to name when display_name is null", () => {
    expect(sampleDisplayName({ id: 1, name: "JC001", display_name: null })).toBe("JC001");
  });
  it("falls back to name when display_name is empty string (uses ||, not ??)", () => {
    expect(sampleDisplayName({ id: 1, name: "JC001", display_name: "" })).toBe("JC001");
  });
  it("falls back to `Sample #id` when both are null/empty", () => {
    expect(sampleDisplayName({ id: 7, name: null, display_name: null })).toBe("Sample #7");
    expect(sampleDisplayName({ id: 7, name: "", display_name: "" })).toBe("Sample #7");
  });
});
```

- [ ] **Step 2: Run, expect failure (file doesn't exist)**

```bash
(cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/lib/sample/displayName.test.ts) 2>&1 | tail -10
```

Expected: import error.

---

### Task 6.2: Implement the helper

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/sample/displayName.ts`

- [ ] **Step 1: Write the helper**

```ts
// displayName.ts — single source of truth for "what string do we render for a sample".
import type { Sample } from "../../api";

// Uses `||` not `??` so an empty-string display_name falls through rather
// than rendering as a blank tile or a leading separator (e.g. " · run.dat"
// in comparison labels). Matches the existing logic in lib/comparison/labels.ts.
export const sampleDisplayName = (
  s: Pick<Sample, "id" | "name" | "display_name">,
): string => s.display_name || s.name || `Sample #${s.id}`;
```

(Note: this depends on the Sample type rename in Task 6.3 to compile cleanly. The vitest run will fail until that lands — that's intentional for the atomic commit.)

---

### Task 6.3: Sample type rename in api.ts

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/api.ts:26-32` (Sample interface)
- Modify: `packages/HimalayaUI/frontend/src/api.ts:103` (updateSample fetcher)
- Modify: `packages/HimalayaUI/frontend/src/api.ts:96` (updateExperiment cleanup)

- [ ] **Step 1: Sample interface**

```ts
export interface Sample {
  id: number;
  experiment_id: number;
  name: string | null;
  display_name: string | null;
  notes: string | null;
  tags: SampleTag[];
}
```

- [ ] **Step 2: updateSample fetcher**

```ts
export const updateSample = (
  id: number,
  patch: { display_name?: string; notes?: string },
  opts?: AuthOpts,
) => request<Sample>("PATCH", `/api/samples/${id}`, patch, opts);
```

- [ ] **Step 3: updateExperiment cleanup**

```ts
// updateExperiment is currently typed against fields the route doesn't accept
// post-issue #88; the route is defensive surface only. Keep the fetcher but
// type it as accepting nothing useful, so callers can't construct a payload
// that would 400.
export const updateExperiment = (
  id: number,
  patch: Record<string, never>,
  opts?: AuthOpts,
) => request<Experiment>("PATCH", `/api/experiments/${id}`, patch, opts);
```

(If `updateExperiment` is unused after this, delete it. Verify with `grep`.)

---

### Task 6.4: Rewrite `mutators/trivial.ts::updateSampleMutator`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts:40, 76-103`

- [ ] **Step 1: UpdateSampleInput type (line 40)**

```ts
export interface UpdateSampleInput {
  display_name?: string;
  notes?: string;
}
```

- [ ] **Step 2: SSE-wins merge (lines 76-85). Drop the `label` line; add `display_name`**

```ts
// Existing comment at lines 76-79 about defensive SSE-wins merge stays.
const patch: Partial<Sample> = {};
if (response.name !== undefined)         patch.name         = response.name;
if (response.display_name !== undefined) patch.display_name = response.display_name;
if (response.notes !== undefined)        patch.notes        = response.notes;
// (label key removed entirely — column no longer exists)
```

- [ ] **Step 3: patchOf helper (lines 98-103)**

```ts
function patchOf(p: UpdateSampleInput): Partial<Sample> {
  const out: Partial<Sample> = {};
  if (p.display_name !== undefined) out.display_name = p.display_name;
  if (p.notes        !== undefined) out.notes        = p.notes;
  return out;
}
```

---

### Task 6.5: Bump SCHEMA_VERSION

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/persistence.ts:6`

- [ ] **Step 1: Edit**

```ts
export const SCHEMA_VERSION = 2;
```

(Add a one-line comment: `// 2: sample.label dropped, payload.name → payload.display_name (issue #88).`)

---

### Task 6.6: SampleMetadataCard rewrite (the rename UI)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/SampleMetadataCard.tsx`

- [ ] **Step 1: Rewrite the input binding**

The component currently binds an `<input>` to `sample.name`. After the refactor, the input edits `sample.display_name` and posts `{ display_name }` on blur:

- Replace `useState(sample.name ?? "")` initialization with `useState(sample.display_name ?? "")`.
- Replace `useEffect(() => setDraft(sample.name ?? ""), [sample.name])` with `…[sample.display_name])`.
- Replace `updateSample(sample.id, { name: draft }, …)` with `…{ display_name: draft }, …`.
- Keep the `data-testid="sample-name-input"` selector unchanged. Add a one-line comment: `// data-testid is historical (issue #88); the field now edits display_name.`

- [ ] **Step 2: Update the breadcrumb (lines 65-76)**

The breadcrumb currently renders `sample.label` as the secondary crumb. After the swap, render `sample.name` (stable identifier) so the breadcrumb shows the identifier and the editable input below shows the friendly text:

```tsx
{(experimentName || sample.name) && (
  <div className="…">
    {experimentName && sample.name && <span> · </span>}
    {sample.name && <span className="font-medium text-fg-dim">{sample.name}</span>}
  </div>
)}
```

---

### Task 6.7: Display-priority sweep — every component using `s.label`

**Files:** see the §7 list in the spec. Do them all in this commit.

- [ ] **Step 1: Replace inline `s.name ?? s.label`/`s.label ?? s.name` chains with `sampleDisplayName(s)`**

For each of the following files, import `sampleDisplayName` from `../../lib/sample/displayName` (adjust relative path) and substitute:

- `queries.ts:208`
- `lib/comparison/labels.ts:38` — keep this one inline (the file has its own `||` documentation); just change the expression to `s.display_name || s.name`. Update the JSDoc at lines 7, 9 too.
- `components/SamplePickerRow.tsx:18`
- `components/NavModal.tsx:127, 204`
- `components/ExposureListRow.tsx:57`
- `components/PlotCard.tsx:83`
- `components/TitleButton.tsx:28`
- `components/MentionPicker.tsx:46, 116`
- `components/MentionChip.tsx:49`

- [ ] **Step 2: NavModal secondary line (line 128)**

```ts
secondary: s.name && s.display_name && s.name !== s.display_name ? s.name : ""
```

- [ ] **Step 3: Search-both-fields in NavModal (line 109) and ComparisonPickerBody (line 119)**

```ts
// NavModal:
[s.name ?? "", s.display_name ?? ""].some(h => h.toLowerCase().includes(needle))
// ComparisonPickerBody (in the haystack array build):
r.sample.name ?? "", r.sample.display_name ?? "", r.sample.notes ?? ""
```

- [ ] **Step 4: ComparisonPickerBody.tsx inline mock fixtures (lines 16, 22, 28)**

Strip `label: null` from each Sample literal. Add `display_name: null` if appropriate.

- [ ] **Step 5: SampleMetadataCard.tsx lines 65, 70, 73, 74**

These are the breadcrumb sites — already covered in Task 6.6.

---

### Task 6.8: Test fixture sweep — Vitest unit tests

**Files:** see the §9-style list in this plan's "File structure" section. Comprehensive sweep.

- [ ] **Step 1: cache-shape.test.ts**

Line 53 — SAMPLE_KEYS Set:

```ts
const SAMPLE_KEYS = new Set([
  "id", "experiment_id", "name", "display_name", "notes", "tags",
]);
```

Lines 329, 335, 355 — strip `label: "D1"` from fixtures, add `display_name: "..."`. Line 342 — `name: "new"` → `display_name: "new"`.

- [ ] **Step 2: sseEventPayload.contract.test.ts**

Line 336 — fixture: drop `label`, add `display_name`. Line 342 — payload `{ name: "new" }` → `{ display_name: "new" }`. Assertion `after.name` → `after.display_name`.

- [ ] **Step 3: rollbackSymmetry.test.ts**

Line 56 (shared SAMPLE fixture), lines 158-163 + 213-227 (updateSample case): same pattern.

- [ ] **Step 4: authHeaders.test.ts**

Lines 157-162 — request mock body `name` → `display_name`.

- [ ] **Step 5: mutatorOnSseWins.test.ts**

Lines 393-428 — entire SSE-wins block: fixtures + payloads + assertions.

- [ ] **Step 6: treats404AsSuccess.test.tsx**

Line 163 — Sample fixture.

- [ ] **Step 7: queries-samples.test.tsx**

Lines 38, 40, 55, 74 — Sample fixtures.

- [ ] **Step 8: Other component tests**

`SamplePickerRow.test.tsx:7`, `ComparisonPicker.test.tsx:45,50`, `ComparisonPickerPanel.test.tsx:9,16`, `ComparisonPickerBody.test.tsx:12,17,115,123`, `ComparePageEdit.test.tsx:638,644`, `NavModal.test.tsx:15,16,20`, `TitleButton.test.tsx:15`, `smoke.test.tsx:56` — drop `label`, add `display_name`.

- [ ] **Step 9: ExposureListRow.test.tsx**

Lines 36, 56, 60. Rename the test name `"falls back to sample.label when sample.name is missing"` to `"falls back to sample.name when sample.display_name is missing"`. Update the assertion to match the new helper behavior.

- [ ] **Step 10: ComparePage.test.tsx**

Lines 153, 296 — legacy comments referencing `sample.label · exposure.filename`. Update the prose.

---

### Task 6.9: Test fixture sweep — Playwright mocked specs

**Files:** every spec in `packages/HimalayaUI/frontend/e2e/` that hand-rolls Sample objects in `page.route` mocks.

- [ ] **Step 1: Per-file edit**

For each of `inspect.spec.ts:8`, `smoke.spec.ts:8,9`, `compare.spec.ts:28,29`, `bundle-b-paper-cuts.spec.ts:25,26`, `figure-export.spec.ts:31`, `multiplayer-replay-rerun.spec.ts:105`:

- Drop `label: …`.
- Add `display_name: …` matching the original `label` value (so the rendered text is unchanged).

---

### Task 6.10: Live spec rewrite — sample-rename-preserves-fields

**Files:**
- Modify: `packages/HimalayaUI/frontend/e2e/live/sample-rename-preserves-fields.spec.ts`

- [ ] **Step 1: Sample type fixture (line 28)**

Drop `label`, add `display_name`.

- [ ] **Step 2: PATCH body (line 92)**

```ts
const r = await ctx.request.patch(`${API}/api/samples/${fx.sampleId}`, {
  data: { display_name: fx.originalName, notes: fx.originalNotes },
  headers: { "X-Username": "playwright" },
});
```

- [ ] **Step 3: Polled assertion (lines 116-124)**

Change `(await getSample()).name` to `(await getSample()).display_name`.

- [ ] **Step 4: afterAll cleanup**

Same swap.

(Do NOT rename the `data-testid="sample-name-input"` selector; that's covered by the historical-comment in Task 6.6.)

---

### Task 6.11: Run frontend full suite, then commit Phase 6

**Files:** none modified. Checkpoint.

- [ ] **Step 1: Vitest**

```bash
(cd packages/HimalayaUI/frontend && npm test) > /tmp/vitest.out 2>&1
grep -E "FAIL|fail|Tests:" /tmp/vitest.out
```

Expected: all green.

- [ ] **Step 2: tsc + vite build**

```bash
(cd packages/HimalayaUI/frontend && npm run build) 2>&1 | tail -10
```

Expected: clean build.

- [ ] **Step 3: Playwright mocked**

```bash
(cd packages/HimalayaUI/frontend && npm run e2e) 2>&1 | tail -15
```

Expected: all green. (If the dev server can't bind to 5173, run `lsof -ti:5173 | xargs kill -9` first.)

- [ ] **Step 4: Single atomic commit for the entire frontend rename**

```bash
git add packages/HimalayaUI/frontend
git commit -m "$(cat <<'EOF'
refactor(frontend): atomic sample.label → display_name rename (issue #88)

- Sample type: drop label, add display_name
- New helper lib/sample/displayName.ts (uses || not ??)
- mutators/trivial.ts::updateSampleMutator full rewrite (UpdateSampleInput,
  patchOf, SSE-wins merge — was missed in original spec)
- SampleMetadataCard.tsx: input edits display_name; breadcrumb shows stable name
- lib/queue/persistence.ts SCHEMA_VERSION bumped 1 → 2 so cross-deploy
  persisted ops drop cleanly
- NavModal + ComparisonPickerBody search both name and display_name
- All 17 component sites + 18 Vitest fixtures + 6 mocked Playwright specs +
  1 live Playwright spec switched in lockstep with the backend wire format

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

# Phase 7 — Docs

Spec reference: §8 step 7.

---

### Task 7.1: Update `docs/experiment-config.md`

**Files:**
- Modify: `docs/experiment-config.md`

- [ ] **Step 1: Update the `[manifest]` schema reference**

Find the section documenting `[manifest]` keys. Replace `label` and `name` with `name` and `display_name`. Add a paragraph:

> The `[manifest].name` column is the **stable scientific identifier** (e.g. `JC001`). It must match `[A-Za-z0-9._-]+`, be non-empty, and be unique within the experiment. It cannot be edited via the UI; rename happens through the manifest CSV + reingest.
>
> The `[manifest].display_name` column is the **friendly user-facing label** (e.g. `DOPC + cholesterol`). It is initialised from the manifest at first ingest and editable via the UI thereafter; reingest never clobbers it.
>
> Migrating an existing `experiment.toml` from the legacy `label/name` shape to the new `name/display_name` shape: run `himalaya migrate-toml <experiment-dir>`. Section-aware regex-anchored substitution; idempotent.

- [ ] **Step 2: Commit**

```bash
git add docs/experiment-config.md
git commit -m "docs(experiment-config): name/display_name schema + migrate-toml (issue #88)"
```

---

### Task 7.2: Update CLAUDE.md gotchas

**Files:**
- Modify: `CLAUDE.md` (root)
- Modify: `packages/HimalayaUI/src/AGENTS.md` (if it documents the samples shape)
- Modify: `packages/HimalayaUI/frontend/src/AGENTS.md` (if it documents Sample type)

- [ ] **Step 1: Add a "samples naming" gotcha if not present**

Add (under HimalayaUI gotchas):

> **`samples.name` is the stable identifier; `samples.display_name` is editable.** Set at ingest from manifest column 2 (was `label`); never UI-mutable. `display_name` is from manifest column 3 (was `name`); user-editable via PATCH; reingest never clobbers. UNIQUE INDEX on `(experiment_id, name)`. Use the `sampleDisplayName(s)` helper in `lib/sample/displayName.ts` everywhere — uses `||` not `??` so empty strings fall through.

- [ ] **Step 2: Verify nested AGENTS.md files don't duplicate stale schema info**

```bash
grep -rn "samples\.label\|samples_label\|samples (label" docs/ AGENTS.md packages/HimalayaUI/ \
  --include="*.md" 2>/dev/null
```

Update any stale references.

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md AGENTS.md packages/HimalayaUI/src/AGENTS.md packages/HimalayaUI/frontend/src/AGENTS.md
git commit -m "docs: update CLAUDE.md/AGENTS.md for samples name/display_name (issue #88)"
```

---

### Task 7.3: CHANGELOG entry

**Files:**
- Modify: `CHANGELOG.md` (or `RELEASE_NOTES.md` — whatever this repo uses; if neither, add to top of `docs/superpowers/specs/2026-05-08-sample-naming-refactor-design.md` "Wire-format changes" appendix)

- [ ] **Step 1: Add a Breaking-Change entry**

```markdown
## Unreleased

### Breaking changes

- **`samples.label` removed; replaced by `name` (stable id, was column 3) and
  `display_name` (editable label, was column 2). Manifest column meanings swapped.**
  Existing experiments need `himalaya migrate-toml <experiment-dir>` to upgrade
  their `experiment.toml`. Existing DBs auto-migrate on first `open_db` after
  deploy. Issue #88.
- **`/api/export` CSV header changed: `sample_label` → `sample_display_name`.**
  Downstream pipelines parsing this CSV need to update their column names.
- **`PATCH /api/samples/:id` no longer accepts `name`; use `display_name`.**
  `samples.name` is now the stable identifier and is set only at ingest time.
- **`PATCH /api/experiments/:id` no longer accepts any field.** The route is
  retained as defensive surface for future fields. Experiment renames must go
  through `experiment.toml` + reingest.
- **First boot after migration purges `idempotent_responses`.** In-flight
  `client_op_id` retries from before the deploy will get fresh executions
  rather than cached pre-rename response bodies.
```

- [ ] **Step 2: Commit**

```bash
git add CHANGELOG.md
git commit -m "docs(changelog): wire-format changes for issue #88"
```

---

### Task 7.4: Verify the web-app design spec hasn't gone stale

**Files:**
- Modify: `docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md` (only if it documents the `samples` shape)

- [ ] **Step 1: Grep for `samples.label`**

```bash
grep -n "samples\.label\|samples (label\|sample.label" \
  docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md
```

- [ ] **Step 2: If matches, update them; otherwise skip**

Update prose describing the samples table shape. Note that the original spec is the historical record — small surgical updates only, with a footnote about issue #88.

- [ ] **Step 3: Commit (if changed)**

```bash
git add docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md
git commit -m "docs(spec): update samples shape reference for issue #88"
```

---

# Phase 8 — Pre-PR verification

This phase is mandatory; the PR cannot be opened until everything below is green.

---

### Task 8.1: Live smoke against the prod-backup snapshot

- [ ] **Step 1: Confirm `/tmp/himalaya-smoke.db` from Task 1.6 is still in place**

```bash
ls -la /tmp/himalaya-smoke.db
```

If absent, re-copy: `cp /opt/Himalaya.jl/data/himalaya.db.pre-issue88.20260509T012719Z /tmp/himalaya-smoke.db`.

- [ ] **Step 2: Open the migrated DB, list samples, verify shape**

```bash
HIMALAYA_DB_PATH=/tmp/himalaya-smoke.db julia --project=packages/HimalayaUI -e '
using HimalayaUI, SQLite, DBInterface, Tables
db = HimalayaUI.open_db()
n  = first(Tables.rowtable(DBInterface.execute(db, "SELECT COUNT(*) AS c FROM samples"))).c
println("samples: $n")
println("first 5: ",
    Tables.rowtable(DBInterface.execute(db,
        "SELECT id, name, display_name FROM samples LIMIT 5")))
' 2>&1 | tail -10
```

Expected: 139 samples; `name` populated, `display_name` populated where `samples.name` was originally non-empty.

- [ ] **Step 3: Pick one experiment and copy its TOML to a scratch path; run migrate-toml on the copy**

```bash
exp_dir=$(HIMALAYA_DB_PATH=/tmp/himalaya-smoke.db julia --project=packages/HimalayaUI -e '
using HimalayaUI, SQLite, DBInterface, Tables
db = HimalayaUI.open_db()
println(first(Tables.rowtable(DBInterface.execute(db, "SELECT path FROM experiments LIMIT 1"))).path)
' 2>&1 | tail -1)
mkdir -p /tmp/smoke-exp && cp "$exp_dir/experiment.toml" /tmp/smoke-exp/
julia --project=packages/HimalayaUI -e '
using HimalayaUI; HimalayaUI.cli_migrate_toml(["/tmp/smoke-exp"])
' 2>&1 | tail -5
diff "$exp_dir/experiment.toml" /tmp/smoke-exp/experiment.toml || true
```

Expected: the diff shows `[manifest].label` line replaced with `name`; old `name` replaced with `display_name`. The `[beamline]` "axis label units" comment is unchanged.

---

### Task 8.2: Full backend test suite

- [ ] **Step 1: Run**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' \
> /tmp/jl-final.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-final.out
tail -50 /tmp/jl-final.out
```

Expected: all testsets pass.

---

### Task 8.3: Full frontend suite

- [ ] **Step 1: Vitest + tsc + vite build + Playwright mocked**

```bash
(cd packages/HimalayaUI/frontend && npm test && npm run build && npm run e2e) \
  > /tmp/fe-final.out 2>&1
grep -E "FAIL|fail|Tests:" /tmp/fe-final.out | tail -10
```

Expected: clean.

---

### Task 8.4: Acceptance criteria checklist

Walk through §11 of the spec against the worktree. Open the spec, tick each box once verified.

- [ ] All 11 acceptance-criteria boxes ticked.

---

### Task 8.5: Open the PR

- [ ] **Step 1: Push the branch**

```bash
git push -u origin sample-naming-refactor
```

- [ ] **Step 2: Open the PR with the issue cross-reference**

Use `gh pr create` with title `refactor(samples): stable name + editable display_name (closes #88, supersedes #83)` and the body summarising the seven phases. Include the prod-backup snapshot path in the description so reviewers can reproduce the smoke test.

---

## Self-review checklist (run before declaring this plan ready)

Walk through the spec section-by-section against this plan:

- [ ] **§1 Problem** — covered by the entire plan.
- [ ] **§2 New data model** — Tasks 1.1 (DDL), 1.3 (migration), 6.3 (Sample type), 5.1 (create_sample!).
- [ ] **§3 Schema migration** — Tasks 1.1–1.6.
- [ ] **§4 Validation module** — Tasks 2.1, 2.2, 5.6 (wired into route + cli).
- [ ] **§5 experiment.toml cutover + migrate-toml** — Tasks 3.1–3.3 (parser), 4.1–4.2 (helper).
- [ ] **§6 Backend route + DB-helper changes** — Tasks 5.1 (create_sample!), 5.2 (cli sites), 5.3 (comparisons + log), 5.4 (export), 5.5 (PATCH samples), 5.6 (PATCH experiments + reingest), 5.7 (transaction).
- [ ] **§7 Frontend changes** — Tasks 6.1–6.10. (Helper, mutator, metadata card, persistence bump, all sites, all fixtures, live spec.)
- [ ] **§8 Build sequence** — mirrored across Phases 1–7.
- [ ] **§9 Testing strategy** — Six-layer audit covered: L1 (Task 5.5), L2 (Task 6.8 step 2), L3 (no change required, verified in spec review), L4 (Task 6.8 step 1), L5+L6 (Task 6.4). Backend full suite: Task 5.9, 8.2. Frontend suite: Task 6.11, 8.3. Live smoke: Task 1.6, 8.1.
- [ ] **§10 Out of scope** — explicitly NOT addressed: experiment.name unique constraint, exposures.filename unique constraint, slug permalinks (#89), comparisons slug column, mention-chip click-through, experiment-name editable/stable split. Plan correctly omits all of these.
- [ ] **§11 Acceptance criteria** — Task 8.4 walks through every box.
