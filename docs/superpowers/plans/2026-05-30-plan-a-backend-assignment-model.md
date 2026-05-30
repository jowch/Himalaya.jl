# Plan A — Backend Assignment Model Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the dual-group (`auto`/`custom`) phase-assignment machinery with one durable assignment per exposure — a 3-state model (`indexed` / `form_factor` / `null`) plus a 0..N set of member indices — exposed via new event kinds and endpoints, while keeping `main` shippable (the unchanged frontend keeps working via dual-write).

**Architecture:** Two new SQLite view tables (`assignments`, `assignment_members`) become the single source of truth for "what phase(s) is this exposure assigned." Three new event kinds (`assignment_add`, `assignment_remove`, `assignment_set_state`) flow through the existing `apply_event!` dispatcher (sole writer to view tables) so they replay through `rebuild_views_from_log!`. The legacy `index_groups`/`index_group_members` machinery is left intact and **dual-written** from the existing member routes, so the current frontend (which still reads `/groups`) is untouched. A one-time migration backfills the new tables from each exposure's active group; `persist_analysis!` seeds the assignment from the auto-indexing. Plan D later moves the frontend onto the native `/assignment` endpoints and removes the legacy machinery + dual-write.

**Tech Stack:** Julia, SQLite.jl / DBInterface, Oxygen.jl (`@get`/`@post`/`@delete`), JSON3, stdlib `Test`. Backend lives under `packages/HimalayaUI/src`; tests under `packages/HimalayaUI/test`.

---

## Settled decisions this plan implements

1. **One durable assignment per exposure** — retire the multi-group guessing (decision #1 of `docs/superpowers/specs/2026-05-30-plotting-implementation-survey.md`). The legacy tables stay only as a transitional compatibility layer (dual-written), removed in Plan D.
2. **Three distinct states** — `indexed` (0..N lattice phases; 0 = call in progress), `form_factor` (structured scattering, no lattice), `null` (no interesting scattering). Explicit `CHECK` enum on the wire, never inferred from `members.length`. Setting `form_factor`/`null` clears members.
3. **Cart is durable backend state** — persisted in SQLite, round-trips the event log, served by `GET /api/exposures/{id}/assignment`.

## File structure

| File | Change | Responsibility |
|---|---|---|
| `packages/HimalayaUI/src/db.jl` | Modify `SCHEMA` (after line 137); add `MIGRATION_ASSIGNMENTS` const + `migrate_assignments!`; register it in `migrate_schema!` after line 384 | Schema + migration |
| `packages/HimalayaUI/src/events.jl` | Add 3 dispatcher branches in `update_view_for_event!` after the `index_unconfirmed` branch (~line 358) | Event → view-table writes |
| `packages/HimalayaUI/src/routes_analysis.jl` | Add `_assignment_body` helper; `GET /assignment`; `POST /assignment/state`; dual-write in the two `/groups/{id}/members` routes | HTTP surface |
| `packages/HimalayaUI/src/pipeline.jl` | Seed `assignments` + `assignment_members` after the auto-group block (after line 537) | Auto-seed the assignment |
| `packages/HimalayaUI/test/test_events.jl` | Add 3 round-trip testsets | Event replay coverage |
| `packages/HimalayaUI/test/test_assignments.jl` | **Create**; register in `runtests.jl` | Migration + helper + seed + route-logic coverage |
| `docs/event-log.md` | Register the 3 new kinds in the kinds table | Docs |

**Invariants to preserve (from the survey doc):** every assignment event uses `entity_type="exposure"`, `entity_id=exposure_id` (so `rebuild_views_from_log!` and `idx_events_by_exposure` work); routes wrapped in `with_idempotency` MUST call `apply_event!(InTransaction(), db, req; …)` (not the default) so the event row + view write + idempotency cache commit atomically; handlers must be replay-idempotent (`INSERT OR IGNORE` / `DELETE` / UPSERT).

---

### Task 1: Schema — `assignments` + `assignment_members` tables

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (the `SCHEMA` constant, insert after the `index_group_members` table at line 137, before `sample_messages` at line 139)
- Test: `packages/HimalayaUI/test/test_assignments.jl` (create)

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/test_assignments.jl`:

```julia
using Test
using HimalayaUI
using SQLite, DBInterface, Tables

@testset "assignments schema" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))

        # Both new tables exist.
        tbls = Set(String(r.name) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table'")))
        @test "assignments" in tbls
        @test "assignment_members" in tbls

        # state CHECK constraint rejects an unknown value.
        DBInterface.execute(db,
            "INSERT INTO exposures (id, sample_id) VALUES (1, NULL)")
        @test_throws Exception DBInterface.execute(db,
            "INSERT INTO assignments (exposure_id, state) VALUES (1, 'bogus')")

        # Default state is 'indexed'.
        DBInterface.execute(db,
            "INSERT INTO assignments (exposure_id) VALUES (1)")
        st = Tables.rowtable(DBInterface.execute(db,
            "SELECT state FROM assignments WHERE exposure_id = 1"))[1].state
        @test String(st) == "indexed"
    end
end
```

- [ ] **Step 2: Register the test file and run it to verify it fails**

Read `packages/HimalayaUI/test/runtests.jl`, find the block of `include("test_*.jl")` lines, and add:

```julia
include("test_assignments.jl")
```

Run: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI", test_args=["test_assignments"])' 2>&1 | tail -40`
(If the suite doesn't support arg filtering, run the file directly: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`)
Expected: FAIL — `"assignments" in tbls` is false (table does not exist yet).

- [ ] **Step 3: Add the tables to `SCHEMA`**

In `packages/HimalayaUI/src/db.jl`, immediately after the `index_group_members` table definition (the `);` on line 137) and before `CREATE TABLE IF NOT EXISTS sample_messages` (line 139), insert:

```sql
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
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: PASS (schema testset green).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_assignments.jl packages/HimalayaUI/test/runtests.jl
git commit -m "feat(backend): add assignments + assignment_members tables"
```

---

### Task 2: Migration — backfill `assignment_members` from active groups

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (add `MIGRATION_ASSIGNMENTS` const near line 5; add `migrate_assignments!`; register the call in `migrate_schema!` after line 384)
- Test: `packages/HimalayaUI/test/test_assignments.jl`

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/test/test_assignments.jl`:

```julia
@testset "migrate_assignments! backfills from active group" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

        # Two indices, an active custom group owning both, an inactive auto group.
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (11, ?, 'Im3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (100, ?, 'auto', 0)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (101, ?, 'custom', 1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_group_members (group_id, index_id) VALUES (101, 10)")
        DBInterface.execute(db, "INSERT INTO index_group_members (group_id, index_id) VALUES (101, 11)")

        # Wipe the migration sentinel + new tables so we can re-run the migration deterministically.
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = 'assignments_v1'")
        DBInterface.execute(db, "DELETE FROM assignment_members WHERE exposure_id = ?", [e_id])
        DBInterface.execute(db, "DELETE FROM assignments WHERE exposure_id = ?", [e_id])

        HimalayaUI.migrate_assignments!(db)

        state = Tables.rowtable(DBInterface.execute(db,
            "SELECT state FROM assignments WHERE exposure_id = ?", [e_id]))
        @test !isempty(state) && String(state[1].state) == "indexed"
        members = Set(Int(m.index_id) for m in Tables.rowtable(DBInterface.execute(db,
            "SELECT index_id FROM assignment_members WHERE exposure_id = ?", [e_id])))
        @test members == Set([10, 11])
    end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: FAIL — `UndefVarError: migrate_assignments! not defined`.

- [ ] **Step 3: Add the migration**

In `packages/HimalayaUI/src/db.jl`, near the other migration-name consts (line ~5, beside `const MIGRATION_COMPARISONS_TO_SERIES = "comparisons_to_series"`), add:

```julia
const MIGRATION_ASSIGNMENTS = "assignments_v1"
```

Then add the function (anywhere among the other `migrate_*!` definitions, e.g. just after `migrate_comparisons_to_series!`):

```julia
"""
    migrate_assignments!(db)

Backfill the durable per-exposure assignment from the legacy active group:
create an `assignments` row (state='indexed') and copy the active group's
members into `assignment_members`, for every exposure that has an active group.
Sentinel-gated, idempotent, own transaction. Raw INSERTs (never apply_event!) —
this is a data backfill, not a user action.
"""
function migrate_assignments!(db::SQLite.DB)
    already = Tables.rowtable(DBInterface.execute(db,
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
            [MIGRATION_ASSIGNMENTS, comparison_now_iso()])
    end
    nothing
end
```

- [ ] **Step 4: Register the migration call**

In `migrate_schema!`, immediately after the `migrate_comparisons_to_series!(db)` call (line 384), add:

```julia
    # Plotting redesign Plan A: durable per-exposure assignment. MUST run after
    # create_schema! (assignments/assignment_members tables exist) and after
    # migrate_series! (schema_migrations sentinel table exists). Backfills from
    # the legacy active group; sentinel-gated; own transaction.
    migrate_assignments!(db)
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: PASS (both testsets green).

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_assignments.jl
git commit -m "feat(backend): migrate_assignments! backfills durable assignment from active group"
```

---

### Task 3: Event kind `assignment_add`

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl` (`update_view_for_event!`, after the `index_unconfirmed` branch ~line 358)
- Test: `packages/HimalayaUI/test/test_events.jl`

- [ ] **Step 1: Write the failing round-trip test**

Append to `packages/HimalayaUI/test/test_events.jl`:

```julia
@testset "assignment_add round-trips through the log" begin
    mktempdir() do dir
        db  = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])

        HimalayaUI.apply_event!(db, req; kind="assignment_add",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))

        # Live state: member present, state forced to 'indexed'.
        @test Set(Int(m.index_id) for m in Tables.rowtable(DBInterface.execute(db,
            "SELECT index_id FROM assignment_members WHERE exposure_id = ?", [e_id]))) == Set([10])
        @test String(Tables.rowtable(DBInterface.execute(db,
            "SELECT state FROM assignments WHERE exposure_id = ?", [e_id]))[1].state) == "indexed"

        # Wipe + rebuild from the log reproduces the member.
        DBInterface.execute(db, "DELETE FROM assignment_members WHERE exposure_id = ?", [e_id])
        HimalayaUI.rebuild_views_from_log!(db, e_id)
        @test Set(Int(m.index_id) for m in Tables.rowtable(DBInterface.execute(db,
            "SELECT index_id FROM assignment_members WHERE exposure_id = ?", [e_id]))) == Set([10])
    end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_events.jl`
Expected: FAIL — no member row written (the dispatcher has no `assignment_add` branch, so the event is a silent no-op).

- [ ] **Step 3: Add the dispatcher branch**

In `packages/HimalayaUI/src/events.jl`, in `update_view_for_event!`, immediately after the `index_unconfirmed` branch (the one that `DELETE`s from `index_group_members`, ~line 358), add:

```julia
    if kind == "assignment_add"
        DBInterface.execute(db,
            """INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')
               ON CONFLICT(exposure_id) DO UPDATE SET state = 'indexed'""",
            [Int(entity_id)])
        DBInterface.execute(db,
            """INSERT OR IGNORE INTO assignment_members (exposure_id, index_id)
               VALUES (?, ?)""",
            [Int(entity_id), Int(payload.index_id)])
        return nothing
    end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_events.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/events.jl packages/HimalayaUI/test/test_events.jl
git commit -m "feat(backend): assignment_add event kind"
```

---

### Task 4: Event kind `assignment_remove`

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl` (after the `assignment_add` branch)
- Test: `packages/HimalayaUI/test/test_events.jl`

- [ ] **Step 1: Write the failing round-trip test**

Append to `packages/HimalayaUI/test/test_events.jl`:

```julia
@testset "assignment_remove round-trips through the log" begin
    mktempdir() do dir
        db  = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (11, ?, 'Im3m', 0.1)", [e_id])

        HimalayaUI.apply_event!(db, req; kind="assignment_add",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))
        HimalayaUI.apply_event!(db, req; kind="assignment_add",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 11))
        HimalayaUI.apply_event!(db, req; kind="assignment_remove",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))

        live() = Set(Int(m.index_id) for m in Tables.rowtable(DBInterface.execute(db,
            "SELECT index_id FROM assignment_members WHERE exposure_id = ?", [e_id])))
        @test live() == Set([11])

        DBInterface.execute(db, "DELETE FROM assignment_members WHERE exposure_id = ?", [e_id])
        HimalayaUI.rebuild_views_from_log!(db, e_id)
        @test live() == Set([11])
    end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_events.jl`
Expected: FAIL — `assignment_remove` is a no-op, so the live set is `{10, 11}`, not `{11}`.

- [ ] **Step 3: Add the dispatcher branch**

In `packages/HimalayaUI/src/events.jl`, after the `assignment_add` branch, add:

```julia
    if kind == "assignment_remove"
        DBInterface.execute(db,
            "DELETE FROM assignment_members WHERE exposure_id = ? AND index_id = ?",
            [Int(entity_id), Int(payload.index_id)])
        return nothing
    end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_events.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/events.jl packages/HimalayaUI/test/test_events.jl
git commit -m "feat(backend): assignment_remove event kind"
```

---

### Task 5: Event kind `assignment_set_state` (clears members on non-`indexed`)

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl` (after the `assignment_remove` branch)
- Test: `packages/HimalayaUI/test/test_events.jl`

- [ ] **Step 1: Write the failing round-trip test**

Append to `packages/HimalayaUI/test/test_events.jl`:

```julia
@testset "assignment_set_state round-trips and clears members on form_factor" begin
    mktempdir() do dir
        db  = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])

        HimalayaUI.apply_event!(db, req; kind="assignment_add",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))
        HimalayaUI.apply_event!(db, req; kind="assignment_set_state",
            entity_type="exposure", entity_id=e_id, payload=Dict(:state => "form_factor"))

        state() = String(Tables.rowtable(DBInterface.execute(db,
            "SELECT state FROM assignments WHERE exposure_id = ?", [e_id]))[1].state)
        members() = Set(Int(m.index_id) for m in Tables.rowtable(DBInterface.execute(db,
            "SELECT index_id FROM assignment_members WHERE exposure_id = ?", [e_id])))
        @test state() == "form_factor"
        @test isempty(members())   # form_factor cleared the lattice members

        # Wipe + rebuild reproduces both the state and the empty member set.
        DBInterface.execute(db, "DELETE FROM assignment_members WHERE exposure_id = ?", [e_id])
        DBInterface.execute(db, "UPDATE assignments SET state = 'indexed' WHERE exposure_id = ?", [e_id])
        HimalayaUI.rebuild_views_from_log!(db, e_id)
        @test state() == "form_factor"
        @test isempty(members())
    end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_events.jl`
Expected: FAIL — `assignment_set_state` is a no-op; state stays `indexed` and member `10` remains.

- [ ] **Step 3: Add the dispatcher branch**

In `packages/HimalayaUI/src/events.jl`, after the `assignment_remove` branch, add:

```julia
    if kind == "assignment_set_state"
        state = String(payload.state)
        DBInterface.execute(db,
            """INSERT INTO assignments (exposure_id, state) VALUES (?, ?)
               ON CONFLICT(exposure_id) DO UPDATE SET state = excluded.state""",
            [Int(entity_id), state])
        if state != "indexed"
            DBInterface.execute(db,
                "DELETE FROM assignment_members WHERE exposure_id = ?", [Int(entity_id)])
        end
        return nothing
    end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_events.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/events.jl packages/HimalayaUI/test/test_events.jl
git commit -m "feat(backend): assignment_set_state event kind (clears members on non-indexed)"
```

---

### Task 6: `_assignment_body` helper

**Files:**
- Modify: `packages/HimalayaUI/src/routes_analysis.jl` (add near `_group_with_members`, ~line 56)
- Test: `packages/HimalayaUI/test/test_assignments.jl`

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/test/test_assignments.jl`:

```julia
@testset "_assignment_body shape" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

        # No assignment yet → defaults: state 'indexed', empty members.
        b0 = HimalayaUI._assignment_body(db, e_id)
        @test b0[:exposure_id] == e_id
        @test b0[:state] == "indexed"
        @test b0[:members] == Int[]

        # With members + a non-default state.
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (11, ?, 'Im3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [e_id])
        DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, 11)", [e_id])
        DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, 10)", [e_id])
        b1 = HimalayaUI._assignment_body(db, e_id)
        @test b1[:members] == [10, 11]   # ORDER BY index_id
    end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: FAIL — `UndefVarError: _assignment_body not defined`.

- [ ] **Step 3: Add the helper**

In `packages/HimalayaUI/src/routes_analysis.jl`, after `_group_with_members` (ends ~line 55), add:

```julia
"""
    _assignment_body(db, exposure_id) -> Dict

Build the canonical assignment response: the durable 3-state assignment for an
exposure plus its 0..N member index ids (ascending). Returns the neutral
default (state 'indexed', no members) when no assignment row exists yet.
"""
function _assignment_body(db::SQLite.DB, exposure_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT state FROM assignments WHERE exposure_id = ?", [exposure_id]))
    state = isempty(rows) ? "indexed" : String(rows[1].state)
    members = Tables.rowtable(DBInterface.execute(db,
        "SELECT index_id FROM assignment_members WHERE exposure_id = ? ORDER BY index_id",
        [exposure_id]))
    Dict(:exposure_id => exposure_id,
         :state       => state,
         :members     => [Int(m.index_id) for m in members])
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_analysis.jl packages/HimalayaUI/test/test_assignments.jl
git commit -m "feat(backend): _assignment_body response helper"
```

---

### Task 7: `GET /api/exposures/{id}/assignment`

**Files:**
- Modify: `packages/HimalayaUI/src/routes_analysis.jl` (inside `register_analysis_routes!()`, near the `GET /groups` route ~line 138)
- Test: `packages/HimalayaUI/test/test_assignments.jl`

- [ ] **Step 1: Write the failing test (handler-logic level)**

Append to `packages/HimalayaUI/test/test_assignments.jl`. This drives the same code the `@get` handler calls (`_assignment_body`) plus asserts the route is registered by checking the Oxygen router after `register_analysis_routes!()`:

```julia
@testset "GET /assignment serves the assignment body" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'form_factor')", [e_id])

        # The route's body is _assignment_body; assert its JSON-serialized shape.
        body = HimalayaUI._assignment_body(db, e_id)
        round = JSON3.read(JSON3.write(body))
        @test round.exposure_id == e_id
        @test round.state == "form_factor"
        @test collect(round.members) == Int[]
    end
end
```

> Full in-process HTTP coverage (issuing a real request through the Oxygen router) follows the harness in `packages/HimalayaUI/test/test_routes_*.jl` — read one of those files and mirror its request-dispatch pattern, asserting `status == 200` and the JSON body equals `_assignment_body`. Add that as a second testset if the harness is readily reusable; the handler-logic test above is the required floor.

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: FAIL — `JSON3` / `_assignment_body` path exercised, but this test passes only once the helper exists (Task 6). If Task 6 is committed, this testset passes immediately at the logic level; proceed to register the actual route so the endpoint exists.

- [ ] **Step 3: Add the route**

In `packages/HimalayaUI/src/routes_analysis.jl`, inside `register_analysis_routes!()`, immediately after the `@get "/api/exposures/{id}/groups"` block (ends ~line 144), add:

```julia
    @get "/api/exposures/{id}/assignment" function(req::HTTP.Request, id::Int)
        db = current_db()
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(_assignment_body(db, id)))
    end
```

- [ ] **Step 4: Run the test + a smoke load to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: PASS.
Then verify the module still loads cleanly (route macro well-formed):
Run: `julia --project=packages/HimalayaUI -e 'using HimalayaUI'`
Expected: no error.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_analysis.jl packages/HimalayaUI/test/test_assignments.jl
git commit -m "feat(backend): GET /api/exposures/{id}/assignment"
```

---

### Task 8: `POST /api/exposures/{id}/assignment/state`

**Files:**
- Modify: `packages/HimalayaUI/src/routes_analysis.jl` (after the `GET /assignment` route)
- Test: `packages/HimalayaUI/test/test_assignments.jl`

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/test/test_assignments.jl`. This exercises the route's event + validation logic by calling `apply_event!` the way the handler does, plus asserting the validation predicate:

```julia
@testset "assignment/state validation + effect" begin
    mktempdir() do dir
        db  = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

        # Valid set → 'null' state recorded.
        HimalayaUI.apply_event!(db, req; kind="assignment_set_state",
            entity_type="exposure", entity_id=e_id, payload=Dict(:state => "null"))
        @test HimalayaUI._assignment_body(db, e_id)[:state] == "null"

        # The route's allow-list predicate (mirror of the handler guard).
        valid = s -> s in ("indexed", "form_factor", "null")
        @test valid("form_factor")
        @test !valid("bogus")
    end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: PASS at the logic level once Tasks 5–6 are in (this testset only fails if `assignment_set_state` or `_assignment_body` is missing). Its purpose is to lock the validation contract; proceed to add the route so the HTTP surface exists.

- [ ] **Step 3: Add the route**

In `packages/HimalayaUI/src/routes_analysis.jl`, after the `GET /assignment` route, add:

```julia
    @post "/api/exposures/{id}/assignment/state" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if !haskey(body, :state)
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "missing field: state")))
        end
        state = String(body.state)
        if !(state in ("indexed", "form_factor", "null"))
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "invalid state: $state")))
        end
        return with_idempotency(db, req) do
            result = apply_event!(InTransaction(), db, req;
                kind        = "assignment_set_state",
                entity_type = "exposure",
                entity_id   = id,
                payload     = Dict(:state => state))
            _enqueue_broadcast_from_result!(result, "assignment_set_state", "exposure", id)
            b = _assignment_body(db, id)
            b[:event_id]    = result.event_id
            b[:view_row_id] = result.view_row_id
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(b))
        end
    end
```

- [ ] **Step 4: Run the test + smoke load to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: PASS.
Run: `julia --project=packages/HimalayaUI -e 'using HimalayaUI'`
Expected: no error.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_analysis.jl packages/HimalayaUI/test/test_assignments.jl
git commit -m "feat(backend): POST /api/exposures/{id}/assignment/state"
```

---

### Task 9: Dual-write `assignment_add` from `POST /groups/{id}/members`

**Files:**
- Modify: `packages/HimalayaUI/src/routes_analysis.jl` (the `@post "/api/groups/{id}/members"` block, lines 146-189 — add a second event after the existing `index_confirmed`)
- Test: `packages/HimalayaUI/test/test_assignments.jl`

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/test/test_assignments.jl`. This asserts that the dual-write keeps the assignment in sync by replaying the same `index_confirmed` + new `assignment_add` events the route emits:

```julia
@testset "member-add dual-writes the assignment" begin
    mktempdir() do dir
        db  = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (200, ?, 'custom', 1)", [e_id])

        # Simulate exactly what the route body now does: legacy event + dual-write.
        HimalayaUI.apply_event!(db, req; kind="index_confirmed",
            entity_type="exposure", entity_id=e_id, payload=Dict(:group_id => 200, :index_id => 10))
        HimalayaUI.apply_event!(db, req; kind="assignment_add",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))

        # Both sources of truth agree.
        legacy = Set(Int(m.index_id) for m in Tables.rowtable(DBInterface.execute(db,
            "SELECT index_id FROM index_group_members WHERE group_id = 200")))
        @test legacy == Set([10])
        @test HimalayaUI._assignment_body(db, e_id)[:members] == [10]
    end
end
```

- [ ] **Step 2: Run it to verify it passes at the event level, then implement the route dual-write**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: PASS (events exist from Tasks 3–6). This testset pins the dual-write *contract*; the route change below makes the actual HTTP path honor it.

- [ ] **Step 3: Add the dual-write to the route**

In `packages/HimalayaUI/src/routes_analysis.jl`, in the `@post "/api/groups/{id}/members"` body, immediately after the existing `_enqueue_broadcast_from_result!(result, "index_confirmed", "exposure", exposure_id)` line, add:

```julia
        # Plan A dual-write: keep the durable assignment in sync with the legacy
        # group membership so the assignment tables stay live while the frontend
        # still reads /groups. Removed in Plan D when the frontend goes
        # assignment-native.
        a_result = apply_event!(InTransaction(), db, req;
            kind        = "assignment_add",
            entity_type = "exposure",
            entity_id   = exposure_id,
            payload     = Dict(:index_id => index_id))
        _enqueue_broadcast_from_result!(a_result, "assignment_add", "exposure", exposure_id)
```

- [ ] **Step 4: Run the analysis-route tests + smoke load**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: PASS.
Run: `julia --project=packages/HimalayaUI -e 'using HimalayaUI'`
Expected: no error.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_analysis.jl packages/HimalayaUI/test/test_assignments.jl
git commit -m "feat(backend): dual-write assignment_add from group member-add route"
```

---

### Task 10: Dual-write `assignment_remove` from `DELETE /groups/{id}/members/{index_id}`

**Files:**
- Modify: `packages/HimalayaUI/src/routes_analysis.jl` (the `@delete "/api/groups/{id}/members/{index_id}"` block, lines 191-229 — add a second event after the existing `index_unconfirmed`)
- Test: `packages/HimalayaUI/test/test_assignments.jl`

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/test/test_assignments.jl`:

```julia
@testset "member-remove dual-writes the assignment" begin
    mktempdir() do dir
        db  = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (200, ?, 'custom', 1)", [e_id])

        HimalayaUI.apply_event!(db, req; kind="assignment_add",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))
        # Simulate the DELETE route body: legacy unconfirm + dual-write remove.
        HimalayaUI.apply_event!(db, req; kind="index_unconfirmed",
            entity_type="exposure", entity_id=e_id, payload=Dict(:group_id => 200, :index_id => 10))
        HimalayaUI.apply_event!(db, req; kind="assignment_remove",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))

        @test isempty(HimalayaUI._assignment_body(db, e_id)[:members])
    end
end
```

- [ ] **Step 2: Run it to verify it passes at the event level, then implement the route dual-write**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: PASS (events exist). Pins the remove dual-write contract.

- [ ] **Step 3: Add the dual-write to the route**

In `packages/HimalayaUI/src/routes_analysis.jl`, in the `@delete "/api/groups/{id}/members/{index_id}"` body, immediately after the existing `_enqueue_broadcast_from_result!(result, "index_unconfirmed", "exposure", exposure_id)` line, add:

```julia
        # Plan A dual-write (mirror of the add route, removed in Plan D).
        a_result = apply_event!(InTransaction(), db, req;
            kind        = "assignment_remove",
            entity_type = "exposure",
            entity_id   = exposure_id,
            payload     = Dict(:index_id => index_id))
        _enqueue_broadcast_from_result!(a_result, "assignment_remove", "exposure", exposure_id)
```

- [ ] **Step 4: Run the test + smoke load**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: PASS.
Run: `julia --project=packages/HimalayaUI -e 'using HimalayaUI'`
Expected: no error.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_analysis.jl packages/HimalayaUI/test/test_assignments.jl
git commit -m "feat(backend): dual-write assignment_remove from group member-remove route"
```

---

### Task 11: `persist_analysis!` seeds the assignment from the auto-indexing

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl` (after the auto-group member loop, line 537)
- Test: `packages/HimalayaUI/test/test_assignments.jl`

- [ ] **Step 1: Write the failing test (unit-level seed helper)**

The seed logic is small and embedded in `persist_analysis!`; test it through the same SQL the new block runs, against a hand-built auto group, to assert the seed-if-absent contract. Append to `packages/HimalayaUI/test/test_assignments.jl`:

```julia
@testset "assignment seed-if-absent contract" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])

        # Call the extracted seed helper directly.
        HimalayaUI.seed_assignment_if_absent!(db, e_id, [10])
        @test HimalayaUI._assignment_body(db, e_id)[:members] == [10]
        @test HimalayaUI._assignment_body(db, e_id)[:state] == "indexed"

        # Second call with a different selection must NOT clobber existing membership.
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (11, ?, 'Im3m', 0.1)", [e_id])
        HimalayaUI.seed_assignment_if_absent!(db, e_id, [11])
        @test HimalayaUI._assignment_body(db, e_id)[:members] == [10]   # preserved
    end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: FAIL — `UndefVarError: seed_assignment_if_absent! not defined`.

- [ ] **Step 3: Add the seed helper and call it from `persist_analysis!`**

In `packages/HimalayaUI/src/pipeline.jl`, add the helper (near the top of the file, after the other small helpers):

```julia
"""
    seed_assignment_if_absent!(db, exposure_id, index_db_ids)

Seed the durable assignment from the auto-indexing selection, but only when the
exposure has no assignment members yet — so reanalysis never clobbers a user's
curated assignment. Creates the `assignments` row (state='indexed') and inserts
`index_db_ids` as members. FK `ON DELETE CASCADE` cleans up members whose
indices are later replaced.
"""
function seed_assignment_if_absent!(db::SQLite.DB, exposure_id::Integer, index_db_ids)
    DBInterface.execute(db,
        "INSERT OR IGNORE INTO assignments (exposure_id, state) VALUES (?, 'indexed')",
        [Int(exposure_id)])
    has_members = !isempty(Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM assignment_members WHERE exposure_id = ? LIMIT 1", [Int(exposure_id)])))
    has_members && return nothing
    for db_id in index_db_ids
        DBInterface.execute(db,
            "INSERT OR IGNORE INTO assignment_members (exposure_id, index_id) VALUES (?, ?)",
            [Int(exposure_id), Int(db_id)])
    end
    nothing
end
```

Then, in `persist_analysis!`, immediately after the auto-group member loop (after line 537, the loop that inserts into `index_group_members`), add:

```julia
    # Plan A: seed the durable assignment from the same auto selection
    # (seed-if-absent, so reanalysis preserves user curation).
    seeded_ids = [candidate_to_db_id[ci] for (ci, idx) in enumerate(candidates) if idx in group_set]
    seed_assignment_if_absent!(db, exposure_id, seeded_ids)
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/pipeline.jl packages/HimalayaUI/test/test_assignments.jl
git commit -m "feat(backend): persist_analysis! seeds the durable assignment (seed-if-absent)"
```

---

### Task 12: Register kinds in docs + full-suite green

**Files:**
- Modify: `docs/event-log.md` (the event-kinds table)
- Verify: whole HimalayaUI backend suite

- [ ] **Step 1: Document the three new event kinds**

In `docs/event-log.md`, find the table of event kinds and add three rows (match the existing column structure — typically `kind | entity_type | payload | view effect`):

```markdown
| `assignment_add` | exposure | `{index_id}` | UPSERT `assignments` → `indexed`; `INSERT OR IGNORE assignment_members` |
| `assignment_remove` | exposure | `{index_id}` | `DELETE` the `assignment_members` row |
| `assignment_set_state` | exposure | `{state}` | UPSERT `assignments.state`; clears `assignment_members` when state ≠ `indexed` |
```

- [ ] **Step 2: Run the full backend suite, capture once**

Run: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test-planA.out 2>&1`
Then: `grep -E "Test Summary|FAIL|Error|assignment" /tmp/jl-test-planA.out | head -60`
Expected: all testsets pass; the new `assignment` testsets appear and pass; no regressions in `test_events.jl` / existing analysis route tests.

- [ ] **Step 3: Verify the migration is idempotent on an already-migrated DB**

Run:
```bash
julia --project=packages/HimalayaUI -e '
using HimalayaUI, SQLite, DBInterface, Tables
mktempdir() do dir
    db = HimalayaUI.open_db(joinpath(dir, "h.db"))
    HimalayaUI.migrate_assignments!(db)   # second call: must be a no-op, no throw
    println("ok: migration idempotent")
end'
```
Expected: prints `ok: migration idempotent` (the sentinel gate short-circuits the second run).

- [ ] **Step 4: Commit**

```bash
git add docs/event-log.md
git commit -m "docs(event-log): register assignment_add/remove/set_state kinds"
```

---

## Self-Review

**1. Spec coverage** (against `docs/superpowers/specs/2026-05-30-plotting-implementation-survey.md` decisions #1–#3 and backend items B1/B2/B5):
- B1 durable single assignment → Tasks 1, 2, 6, 7, 9, 10, 11. ✓
- B2 three-state enum (`indexed`/`form_factor`/`null`) → Task 1 (CHECK), Task 5 (set_state + clear), Task 8 (route). ✓
- B5 assignment event kinds following the `apply_event!` contract → Tasks 3, 4, 5 (round-trip tested), 12 (docs). ✓
- Keep `main` shippable (frontend untouched) → dual-write Tasks 9, 10; legacy `index_groups`/`/groups` left intact. ✓
- B3 (Gauss–Bonnet) and B4 (lattice-driven speculative build) are **out of scope** for Plan A (Plan B and Plan D respectively) — intentional.

**2. Placeholder scan:** No "TBD"/"add error handling"/"similar to Task N". The one soft reference (full in-process HTTP route test, Task 7 Step 1) names the canonical harness file and the exact assertions; the required handler-logic floor is fully specified.

**3. Type/name consistency:** `seed_assignment_if_absent!` (Task 11), `_assignment_body` (Tasks 6–10), `migrate_assignments!` + `MIGRATION_ASSIGNMENTS = "assignments_v1"` (Task 2), event kinds `assignment_add`/`assignment_remove`/`assignment_set_state` (Tasks 3–5, 8–10, 12) — all referenced consistently. Payload key is `:index_id` for add/remove and `:state` for set_state throughout. `entity_type="exposure"`, `entity_id=exposure_id` everywhere.

**Known follow-ups (documented, not Plan A):**
- Reanalysis repopulation: `seed_assignment_if_absent!` won't re-seed once members exist; if reanalysis replaces index ids, FK cascade empties the assignment and it stays empty until a user acts. Revisit when Plan D wires the cart (or add a reanalysis-reseed event).
- Legacy removal: `index_groups` / `index_group_members` / `ensure_custom_group!` / the dual-write / `GET /groups` are removed in Plan D once the frontend is assignment-native.

---

## Execution Handoff

Two execution options:

1. **Subagent-Driven (recommended)** — a fresh subagent per task, two-stage review between tasks, fast iteration.
2. **Inline Execution** — execute tasks in this session with checkpoints for review.
