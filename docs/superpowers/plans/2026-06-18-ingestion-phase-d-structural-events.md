# Ingestion Redesign — Phase D: Structural-Edit Events Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire the four grouping-edit operations — rename a sample, move one exposure to a different sample, merge two samples, split a sample at an exposure boundary — into the event-log / idempotency / SSE rails. After this phase, every structural edit from the grouping-review surface is durable, broadcast over SSE, and individually undoable (single-entity events) or session-locally undoable (multi-entity orchestrations).

**Architecture:**

- `exposure_moved` and `sample_renamed` are genuine single-entity event kinds. `exposure_moved` gets a dispatcher branch in `update_view_for_event!` (`events.jl:305`) that writes the view table directly; `sample_renamed` is a dispatcher no-op (the route writes `samples.name` directly, mirroring `update_sample` at `events.jl:405`). Both inherit `with_idempotency` + SSE broadcast and round-trip through `rebuild_views_from_log!` (`events.jl:1148`, which already accepts an `entity_type` kwarg — confirmed below).
- **Every structural event payload carries `experiment_id`** (spec §9.3, resolved 2026-06-19) so the frontend SSE arms can invalidate the right experiment's `loads(id)`/`samples(id)` caches without reading `remote.entity_id` (which is the exposure/sample id, not the experiment). Canonical payloads: `exposure_moved` = `{ sample_id, from_sample_id, experiment_id }` (dispatcher writes `exposures.sample_id` from `sample_id`, the destination; `from_sample_id` lets the frontend invalidate the source sample too); `sample_renamed` = `{ name, experiment_id }`; `sample_created` = `{ experiment_id }`; `sample_split` = `{ new_sample_id, exposure_ids, experiment_id }`; `grouping_flag_dismissed` = `{ flag_kind, merge_with_sample_id?, experiment_id }`. **There is NO `sample_merged` kind** — a merge fans out as per-exposure `exposure_moved` frames (+ a `sample_renamed`/`update_sample` if the survivor is renamed) and retires the loser via `merged_into_id`.
- `merge` and `split` are route-level orchestrations. Each route body is wrapped in one `with_idempotency` block; inside it issues per-exposure `apply_event!(InTransaction(), …)` calls (each individually idempotent via the `(client_op_id, action, entity_id)` UNIQUE on `user_actions`) plus explicit sample create/retire writes, and calls `_enqueue_broadcast_from_result!` per event INSIDE the block (the call only appends to a task-local queue; `with_idempotency` flushes it after the tx commits — see Conventions). `split` emits a `sample_created` event for the new sample + a `sample_split` event on the source; `merge` does NOT emit `sample_created`. Undo is session-local (`useUndoStack` on the frontend); there is no server-side multi-row undo API.
- `grouping_flag_dismissed` is a durable, single-entity, **undoable** event keyed to the flagged sample. Its dispatcher is a no-op — the merge/split flag is never stored on `samples`; it is re-derived on every roll-up read by `derive_sample_flags` (§9.1), and a non-undone `grouping_flag_dismissed` event *suppresses* the matching flag at read time. An `undoes_event_id`-stamped undo re-shows it.
- All routes live in a new file `routes_grouping.jl`, registered at server startup. They clone the idempotent-mutation pattern from `routes_samples.jl:150–186` (the `PATCH /api/samples/{id}` handler), which calls `_enqueue_broadcast_from_result!` inside its `with_idempotency` block.
- The grouping view (Load ▸ Sample ▸ Exposures) is re-derived from `(load_id, slot_index)`, never folded from the log — consistent with the "rebuild-not-log-derivable" precedent from the plotting redesign. The per-sample `flag` on the roll-up is produced by Phase B's `derive_sample_flags` (§9.1) over the rows the roll-up reads, with dismissed flags suppressed.

**Tech Stack:** Julia, SQLite.jl / DBInterface.jl / HTTP.jl / JSON3.jl / Oxygen.jl, stdlib `Test`. Backend package at `packages/HimalayaUI/`.

**Spec:** `docs/superpowers/specs/2026-06-15-ingestion-redesign-design.md` §9.3 (event kinds, merge orchestration, undo, series/tag dedup, FK re-point tables) and §9.2 (REST endpoints — move/rename/merge/split routes). Read §9.3 in full before starting.

**Source of truth for current code:** the build-kit anchors in this plan were line-verified 2026-06-18, but line numbers drift — confirm each anchor with a quick `grep`/read before editing.

---

## File Structure

| File | Responsibility | This plan |
|---|---|---|
| `packages/HimalayaUI/src/events.jl` | `update_view_for_event!` dispatcher (Task 1 `exposure_moved` branch; Task 8 optional no-op markers; Task 8b `grouping_flag_dismissed` no-op). `rebuild_views_from_log!` already accepts `entity_type` — Task 2 needs no edit here. | MODIFY (Tasks 1, 8, 8b — NOT Task 2) |
| `packages/HimalayaUI/src/db.jl` | `retire_sample!` writer; `MIGRATION_SAMPLES_MERGED_INTO` migration; `get_loads_rollup` (calls Phase B's `derive_sample_flags` + suppresses dismissed flags) | MODIFY (Tasks 3, 4) |
| `packages/HimalayaUI/src/routes_grouping.jl` | rename / move / merge / split / dismiss-flag / loads REST routes | CREATE (Tasks 5–9) |
| `packages/HimalayaUI/src/HimalayaUI.jl` | `include("routes_grouping.jl")` before `include("server.jl")` | MODIFY (Task 5) |
| `packages/HimalayaUI/src/server.jl` | `register_grouping_routes!()` call inside `register_routes!()` | MODIFY (Task 5) |
| `packages/HimalayaUI/test/test_ingestion_structural_events.jl` | new standalone test file for this phase | CREATE (all tasks) |
| `packages/HimalayaUI/test/runtests.jl` | test registry | MODIFY (Task 0) |

**Out of scope (other plans):**

- Phase A — schema migrations (`loads` table, `exposures.experiment_id`, `samples` name collapse).
- Phase B — `prp.jl`, `geometry.jl`, `grouping.jl`, `scan_and_group!`.
- Phase C — scan/rescan routes, `broadcast_progress!`, auto-rescan scheduler.
- Phase E — frontend `/experiments/:id` page, grouping-review surface, `useUndoStack` extraction, `display_name → name` TS sweep.

This plan depends on Phases A, B, and C having run (the `loads` table and `exposures.experiment_id` column must exist; `scan_and_group!` must populate them). If running in isolation on a test DB, the test harness seeds the required rows directly.

> **Cross-plan seam — `derive_sample_flags` (resolved; aligned with Phase B Task 12).** Spec §9.1 / §8.8 pin **`grouping.jl` exposes `derive_sample_flags(load_rows) → Dict{Int, GroupingFlag}`** (keyed by `sample_id`) — a pure read-time function returning the per-sample merge/split suggestion that `get_loads_rollup` (Task 4) attaches to each sample and the frontend consumes verbatim. **Phase B Task 12 (`2026-06-18-ingestion-phase-b-ingest-core.md`) implements it**, and its input contract is confirmed to match this plan's Task 4 call site: it consumes the **nested roll-up rows** (a `Vector` of loads, each with `samples`, each sample with `sample_id`/`name`/`slot_index` + its `exposures` carrying `filename`/`horizontal_position`/`timestamp` — it reads the roll-up field name `horizontal_position`, not the pre-persistence `_mm` form), NOT the raw `ExposureMeta` stream from `scan_directory`. Task 4 passes exactly that shape. The `isdefined(HimalayaUI, :derive_sample_flags)` guard at the call site is **defensive graceful-degradation for incremental build order** (so Phase D's roll-up still loads and tests cleanly if it is built before Phase B Task 12 lands in the same tree, emitting `flag = null` until the function is present) — not a sign the function is missing. **One minor field-name reconciliation:** Phase B Task 12's *docstring* still lists the exposure leaf as `(exposure_id, …, frame_no, status)`, the pre-rewrite shape; this plan's Task 4 emits the §8.8 leaf `(id, filename, horizontal_position, timestamp)`. `derive_sample_flags` only reads `sample_id`/`slot_index`/`horizontal_position`/`timestamp`, all present and identically named, so the logic is unaffected — Phase B's docstring should be synced to the §8.8 leaf names for clarity.

---

## Conventions for every task

- **Run a single test file** from the repo root:
  `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
  (Standalone: the file uses `HimalayaUI.open_db` to build a temp DB, runs all `@testset` blocks, and exits. The full suite `Pkg.test("HimalayaUI")` is for CI.)
- Each test builds its DB in a **`mktempdir()` block** so tests never touch a real database.
- **Commit after each task** once its test passes. Exact `git add` lists are given per task.
- The `InTransaction` variant of `apply_event!` does NOT broadcast — routes wrapped in `with_idempotency` must call `_enqueue_broadcast_from_result!` explicitly. **Call it INSIDE the `with_idempotency do … end` block**, not after it. `_enqueue_broadcast_from_result!` only appends the frame to a task-local queue (`POST_COMMIT_BROADCAST_KEY`); `with_idempotency` itself flushes that queue via `_flush_post_commit_broadcasts!` *after* its `SQLite.transaction` commits (`idempotency.jl:189`), and clears it on rollback or cache-replay (`idempotency.jl:103,147,169`). By the time the `do` block returns, the queue has already been flushed/cleared, so a call placed after the block would enqueue a frame that never fires. This matches the canonical live route `routes_samples.jl:179` (the `PATCH /api/samples/{id}` handler calls `_enqueue_broadcast_from_result!` inside its `with_idempotency` block). On an idempotent cache-replay the body does not execute, so nothing is enqueued and the cached HTTP response is returned without re-broadcasting — correct multiplayer semantics (the original event already broadcast on first execution).

---

## Task 0: Test harness scaffold

**Files:**
- Create: `packages/HimalayaUI/test/test_ingestion_structural_events.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Create the standalone test file**

```julia
# packages/HimalayaUI/test/test_ingestion_structural_events.jl
using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# Standalone-run support: under runtests.jl the shared HTTP harness
# (with_test_server / start_test_server! / bind_db!) is already loaded
# (test_http.jl is included first); on a standalone run pull it in.
# Mirrors test_assignment_reattach.jl:28–30.
if !isdefined(@__MODULE__, :with_test_server)
    include("test_http.jl")
end

# ── helpers ─────────────────────────────────────────────────────────────────

"Open a fresh DB with the full current schema."
function fresh_db()
    dir = mktempdir()
    HimalayaUI.open_db(joinpath(dir, "h.db"))
end

"Seed one experiment + one load + two samples + two exposures (one per sample).
Caller already holds `db`; returns (exp_id, load_id, s1_id, s2_id, e1_id, e2_id)."
function seed_two_samples(db)
    exp_id = HimalayaUI.create_experiment!(db; path="/d", data_dir="/d", analysis_dir="/a")
    # load_id: Phase A adds the loads table; if running in isolation seed it manually.
    load_id = DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO loads (experiment_id, load_index) VALUES (?, 1)", [exp_id]))
    s1_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, load_id=load_id,
        slot_index=1, name="HA85 (S01P01)")
    s2_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, load_id=load_id,
        slot_index=2, name="HA85 (S01P02)")
    e1_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id,
        sample_id=s1_id, filename="f01.tif")
    e2_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id,
        sample_id=s2_id, filename="f02.tif")
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id)
end

"POST-like request with X-Username, for direct apply_event! unit tests
(Tasks 1–2). Route tests drive real HTTP via with_test_server instead."
user_req(name="alice") = HTTP.Request("POST", "/x", ["X-Username" => name], UInt8[])

@testset "structural events (Phase D)" begin
    # task testsets append below
end
```

> **Note:** `create_sample!` and `create_exposure!` signatures used here reflect the Phase A rewrites (they accept `load_id`, `slot_index`, `experiment_id` kwargs). If Phase A is not yet merged, update the calls to match the current signatures — the Phase A plan documents the exact new signatures.
>
> **Why `with_test_server`, not a route-dispatch shim:** routes are registered via Oxygen `@patch`/`@post`/`@get` macros into Oxygen's global router; there is no `(req, db)` dispatch function to call directly. The only way to exercise a route is over real in-process HTTP. `with_test_server(db)` (defined in `test_http.jl`) calls `bind_db!(db)` + `start_test_server!(db, port)` so the handler's `current_db()` resolves to the test `db`, then yields `(port, base)` for `HTTP.get/patch/post("$base/api/…")`. This is the canonical pattern in `test_routes_samples.jl`. Because `start_test_server!` calls `register_routes!()` (server.jl:200) — which after Task 9 includes `register_grouping_routes!()` — the grouping routes are registered automatically; **route tests must NOT call `register_grouping_routes!()` themselves**, and the seed/db they pass to `with_test_server` is the one the handlers read.

- [ ] **Step 2: Register in runtests.jl**

Find the `include(...)` list in `packages/HimalayaUI/test/runtests.jl` (line 33 area) and add in alphabetical position:

```julia
include("test_ingestion_structural_events.jl")
```

- [ ] **Step 3: Run to verify harness loads**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: PASS (empty `@testset` always passes; proves imports resolve).

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/test/test_ingestion_structural_events.jl packages/HimalayaUI/test/runtests.jl
git commit -m "test: scaffold structural-events test harness (Phase D)"
```

---

## Task 1: `exposure_moved` dispatcher branch + round-trip

`exposure_moved` is keyed to the exposure (`entity_type = "exposure"`, `entity_id = exposure_id`). Its canonical payload is `{ sample_id, from_sample_id, experiment_id }` (spec §9.3): `sample_id` is the destination (the new owner), `from_sample_id` is the source, `experiment_id` lets the frontend invalidate the right experiment's caches. **The dispatcher reads ONLY `payload.sample_id`** and writes `exposures.sample_id` directly — it is unchanged by the added `from_sample_id`/`experiment_id` fields (extra payload keys are inert in `update_view_for_event!`). It is replay-idempotent (an `UPDATE` is naturally idempotent).

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl` (`update_view_for_event!`, ~line 305; `rebuild_views_from_log!`, ~line 1148)
- Test: `packages/HimalayaUI/test/test_ingestion_structural_events.jl`

- [ ] **Step 1: Write the failing test**

Append inside the `@testset "structural events (Phase D)"` block:

```julia
@testset "exposure_moved — dispatcher writes exposures.sample_id" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)
    req = user_req()

    # Move e1 (currently in s1) to s2. Use the InTransaction variant — the
    # broadcasting default-method variant calls broadcast_event! → current_db()
    # (events.jl:1106) on a non-nothing user_id, which errors with no server
    # bound. InTransaction directly exercises the dispatcher and does NOT
    # broadcast (events.jl:101). It assumes an open tx, so wrap it.
    result = SQLite.transaction(db) do
        HimalayaUI.apply_event!(HimalayaUI.InTransaction(), db, req;
            kind = "exposure_moved",
            entity_type = "exposure",
            entity_id = e1_id,
            payload = Dict(:sample_id => s2_id,
                           :from_sample_id => s1_id,
                           :experiment_id => exp_id))
    end
    @test result.event_id > 0

    row = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT sample_id FROM exposures WHERE id = ?", [e1_id])))
    @test Int(row.sample_id) == s2_id

    # The dispatcher ignores the extra payload keys — it reads only sample_id.
    # Confirm the durable payload round-trips all three fields for the SSE arm.
    pl = JSON3.read(String(first(Tables.rowtable(DBInterface.execute(db,
        "SELECT payload FROM user_actions WHERE id = ?", [result.event_id]))).payload))
    @test Int(pl.sample_id) == s2_id
    @test Int(pl.from_sample_id) == s1_id
    @test Int(pl.experiment_id) == exp_id
end

@testset "exposure_moved — rebuild_views_from_log! round-trip" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)
    req = user_req()

    SQLite.transaction(db) do
        HimalayaUI.apply_event!(HimalayaUI.InTransaction(), db, req;
            kind = "exposure_moved", entity_type = "exposure", entity_id = e1_id,
            payload = Dict(:sample_id => s2_id,
                           :from_sample_id => s1_id, :experiment_id => exp_id))
    end

    # Reset view state (undo the dispatcher's write) and replay from log.
    DBInterface.execute(db, "UPDATE exposures SET sample_id = ? WHERE id = ?", [s1_id, e1_id])
    @test Int(first(Tables.rowtable(DBInterface.execute(db,
        "SELECT sample_id FROM exposures WHERE id = ?", [e1_id]))).sample_id) == s1_id

    HimalayaUI.rebuild_views_from_log!(db, e1_id; entity_type="exposure")
    row = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT sample_id FROM exposures WHERE id = ?", [e1_id])))
    @test Int(row.sample_id) == s2_id
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: FAIL (`exposure_moved` falls through `update_view_for_event!` without writing anything).

- [ ] **Step 3: Add the dispatcher branch**

Read `events.jl:305–424` to confirm current branch layout. Add the `exposure_moved` branch **after** `assignment_set_state` and **before** the `kind == "update_sample"` no-op line:

```julia
    if kind == "exposure_moved"
        # payload.sample_id: the destination sample. UPDATE is naturally idempotent —
        # replaying sets the same value, so rebuild_views_from_log! is safe.
        DBInterface.execute(db,
            "UPDATE exposures SET sample_id = ? WHERE id = ?",
            [Int(payload.sample_id), Int(entity_id)])
        return nothing
    end
```

Also add an intent-marker no-op for `sample_renamed` now (it has no view write — the PATCH route writes `samples.name` directly, matching the `update_sample` precedent at `events.jl:405`):

```julia
    kind == "sample_renamed" && return nothing
```

This line goes **before** the existing `kind == "update_sample" && return nothing` line.

> **On the no-op branch:** it is NOT required for correctness. The dispatcher's default is already a silent no-op (`events.jl:516–517` returns `nothing` for any unrecognized kind), so an unknown kind does not throw and `rebuild_views_from_log!` tolerates it as-is. This line is an explicit intent marker — it documents that `sample_renamed` is deliberately a no-op (not an oversight), matching the existing retired-kind guards (`index_confirmed`/`index_unconfirmed` at `events.jl:363`, the M2.1 trivial-route block at `events.jl:405–411`). Adding it is a style choice for readability; the `exposure_moved` branch above is the only functionally load-bearing change in this task.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/events.jl packages/HimalayaUI/test/test_ingestion_structural_events.jl
git commit -m "feat(events): exposure_moved + sample_renamed dispatcher branches"
```

---

## Task 2: `sample_renamed` — event recorded, SSE-ready, log is the audit trail

`sample_renamed` is keyed to the sample (`entity_type = "sample"`, `entity_id = sample_id`). The dispatcher is a no-op (the route writes `samples.name` directly, identical to `update_sample` today at `routes_samples.jl:168–170`). The event's purpose is a durable audit row + SSE broadcast. `rebuild_views_from_log!` must tolerate `entity_type = "sample"`.

> **VERIFIED — no events.jl change needed in this task.** `rebuild_views_from_log!` (`events.jl:1148–1149`) **already** declares `entity_type::String = "exposure"` as a keyword argument and threads it into the WHERE clause (`events.jl:1153`: `WHERE entity_type = ? AND entity_id = ?` bound to `[entity_type, Int(entity_id)]`). Passing `entity_type="sample"` already routes to the sample partition of the log. The hedge in earlier drafts ("the function may need an `entity_type` kwarg added") is stale — confirm by reading `events.jl:1148–1164`, then skip the events.jl edit. This task's only new code is the test asserting the existing kwarg path tolerates `entity_type="sample"`.

**Files:**
- Test: `packages/HimalayaUI/test/test_ingestion_structural_events.jl` (only)
- (No `events.jl` change — the `entity_type` kwarg already exists.)

- [ ] **Step 1: Write the failing test**

```julia
@testset "sample_renamed — event durable + rebuild tolerates entity_type='sample'" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)
    req = user_req()

    # InTransaction variant — no broadcast (see Task 1 rationale).
    result = SQLite.transaction(db) do
        HimalayaUI.apply_event!(HimalayaUI.InTransaction(), db, req;
            kind = "sample_renamed",
            entity_type = "sample",
            entity_id = s1_id,
            payload = Dict(:name => "JC C04 (S01P01)", :experiment_id => exp_id))
    end
    @test result.event_id > 0

    # The event is durable.
    row = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT action, entity_type, entity_id FROM user_actions WHERE id = ?",
        [result.event_id])))
    @test String(row.action) == "sample_renamed"
    @test String(row.entity_type) == "sample"
    @test Int(row.entity_id) == s1_id

    # rebuild_views_from_log! with entity_type="sample" must not throw.
    @test_nowarn HimalayaUI.rebuild_views_from_log!(db, s1_id; entity_type="sample")
end
```

- [ ] **Step 2: Confirm the existing kwarg path (read events.jl, no edit)**

Read `events.jl:1148–1164` and confirm: the signature is `rebuild_views_from_log!(db::SQLite.DB, entity_id::Integer; entity_type::String = "exposure")` and the WHERE binds `[entity_type, Int(entity_id)]`. Because the kwarg already exists and is threaded, `entity_type="sample"` already works — **no events.jl edit is required**. This is a characterization test, not a TDD red→green cycle: it pins behavior that already holds so a future refactor of `rebuild_views_from_log!` cannot silently break the sample-partition path.

- [ ] **Step 3: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: PASS immediately (the kwarg path already exists). If it FAILS with a `MethodError` on the `entity_type` kwarg, the live function has drifted from `events.jl:1148` — re-read and only then widen the signature.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/test/test_ingestion_structural_events.jl
git commit -m "test(events): pin sample_renamed durable audit + rebuild entity_type='sample' path"
```

---

## Task 3: `retire_sample!` writer + `merged_into_id` column migration

Merge requires a way to mark the loser as retired without hard-deleting it (spec §9.3: nullable `merged_into_id`, not a hard delete; `series_samples` `ON DELETE CASCADE` would silently drop the series membership). Add the column as a migration and add a `retire_sample!` helper.

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (new sentinel + migration + `retire_sample!`)
- Test: `packages/HimalayaUI/test/test_ingestion_structural_events.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "samples.merged_into_id column exists after migration" begin
    db = fresh_db()
    col_names = String.(getproperty.(Tables.rowtable(
        DBInterface.execute(db, "PRAGMA table_info(samples)")), :name))
    @test "merged_into_id" in col_names
end

@testset "retire_sample! sets merged_into_id" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

    HimalayaUI.retire_sample!(db, s2_id; merged_into_id=s1_id)

    row = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT merged_into_id FROM samples WHERE id = ?", [s2_id])))
    @test Int(row.merged_into_id) == s1_id

    # The loser row still exists (no hard delete).
    @test !isempty(Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM samples WHERE id = ?", [s2_id])))
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: FAIL (`merged_into_id` column absent, `retire_sample!` undefined).

- [ ] **Step 3: Add the sentinel, migration, and writer**

In `db.jl`, next to the other `const MIGRATION_*` declarations, add:

```julia
const MIGRATION_SAMPLES_MERGED_INTO = "MIGRATION_SAMPLES_MERGED_INTO"
```

In `migrate_schema!`, add a call after the existing migration chain:

```julia
    migrate_samples_merged_into!(db)
```

Add the migration function:

```julia
function migrate_samples_merged_into!(db::SQLite.DB)
    _migrated(db, MIGRATION_SAMPLES_MERGED_INTO) && return nothing
    SQLite.transaction(db) do
        existing = cols_of(db, "samples")
        "merged_into_id" in existing || DBInterface.execute(db,
            "ALTER TABLE samples ADD COLUMN merged_into_id INTEGER REFERENCES samples(id)")
        _record_migration!(db, MIGRATION_SAMPLES_MERGED_INTO)
    end
end
```

Add the writer below `create_sample!` (`db.jl:1828`):

```julia
"""
    retire_sample!(db, loser_id; merged_into_id)

Mark a sample as retired by pointing its `merged_into_id` at the survivor.
Does NOT hard-delete the row — FKs from `series_samples` etc. are re-pointed
before this is called (see merge route). Sets `merged_into_id` only; no SSE.
"""
function retire_sample!(db::SQLite.DB, loser_id::Integer;
                         merged_into_id::Integer)
    DBInterface.execute(db,
        "UPDATE samples SET merged_into_id = ? WHERE id = ?",
        [Int(merged_into_id), Int(loser_id)])
    nothing
end
```

> **Helper availability:** `_migrated`, `_record_migration!`, and `cols_of` are defined in Phase A's plan. If Phase A is not yet merged, either cherry-pick or define them here (see Phase A plan, Task 1 and Task 3 helper notes).

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_ingestion_structural_events.jl
git commit -m "feat(db): samples.merged_into_id migration + retire_sample! writer"
```

---

## Task 4: Grouping-view query

The grouping-review surface needs a `GET /api/experiments/{id}/loads` endpoint that returns the Load ▸ Sample ▸ Exposures roll-up (spec §9.2 and §9.3: "keyed `queryKeys.loads(id)`, distinct from the flat corpus"). The view is re-derived from `(load_id, slot_index)` — NOT folded from the log — so no event log is involved. Task 4 only adds the query helper used by Task 5's route.

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (add `get_loads_rollup`)
- Test: `packages/HimalayaUI/test/test_ingestion_structural_events.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "get_loads_rollup returns nested Load→Sample→Exposure tree (§8.8 shape)" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

    rollup = HimalayaUI.get_loads_rollup(db, exp_id)
    @test length(rollup) == 1
    load = first(rollup)
    # Load-level field names match the §8.8 contract exactly.
    @test load.load_id == load_id
    @test haskey(load, :load_index)
    @test haskey(load, :session_id)
    @test haskey(load, :start_time)
    @test haskey(load, :end_time)
    @test haskey(load, :frame_count)
    @test haskey(load, :note)
    @test length(load.samples) == 2

    sample_ids = Set(s.sample_id for s in load.samples)
    @test s1_id in sample_ids && s2_id in sample_ids
    for s in load.samples
        # Sample-level field names match the §8.8 LoadSample contract exactly.
        @test haskey(s, :name)
        @test haskey(s, :slot_index)
        @test haskey(s, :grouping_source)
        @test haskey(s, :name_source)
        @test haskey(s, :merged_into_id)
        @test haskey(s, :flag)         # GroupingFlag or nothing
        @test length(s.exposures) == 1
        ex = first(s.exposures)
        # Exposure-leaf field names match the §8.8 LoadExposure contract exactly.
        @test haskey(ex, :id)
        @test haskey(ex, :filename)
        @test haskey(ex, :horizontal_position)
        @test haskey(ex, :timestamp)
    end
end

@testset "get_loads_rollup attaches a present merge flag, and suppresses a dismissed one" begin
    # Phase B Task 12's derive_sample_flags is the authority for `flag`. This
    # test is written to pass regardless of incremental build order:
    #  - if derive_sample_flags is not yet loaded in this tree, every flag is
    #    `nothing` (the guard below skips the present-flag assert).
    #  - once Phase B Task 12 is present, the rollup attaches its output, and a
    #    non-undone grouping_flag_dismissed event on s1 suppresses s1's flag.
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

    have_flags = isdefined(HimalayaUI, :derive_sample_flags)

    rollup = HimalayaUI.get_loads_rollup(db, exp_id)
    s1 = first(filter(s -> s.sample_id == s1_id, first(rollup).samples))
    # Before any dismiss: flag is whatever derive_sample_flags returns
    # (nothing when the function is absent).
    have_flags || @test s1.flag === nothing

    # Record a non-undone grouping_flag_dismissed event on s1 directly.
    req = user_req()
    SQLite.transaction(db) do
        HimalayaUI.apply_event!(HimalayaUI.InTransaction(), db, req;
            kind = "grouping_flag_dismissed", entity_type = "sample", entity_id = s1_id,
            payload = Dict(:flag_kind => "merge", :merge_with_sample_id => s2_id,
                           :experiment_id => exp_id))
    end

    rollup2 = HimalayaUI.get_loads_rollup(db, exp_id)
    s1b = first(filter(s -> s.sample_id == s1_id, first(rollup2).samples))
    # A dismissed flag is suppressed regardless of what derive_sample_flags returns.
    @test s1b.flag === nothing
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: FAIL (`get_loads_rollup` undefined).

- [ ] **Step 3: Implement `get_loads_rollup`**

Add below `get_exposures` in `db.jl` (~line 1865):

```julia
"""
    get_loads_rollup(db, experiment_id) -> Vector{NamedTuple}

Return the Load ▸ Sample ▸ Exposures roll-up for the grouping-review surface,
in the exact nested shape spec §8.8 pins (field names are load-bearing — the
frontend `Load`/`LoadSample`/`LoadExposure` interfaces mirror them verbatim).

Each Load NamedTuple:
    (load_id, load_index, session_id, start_time, end_time, frame_count, note, samples)
each LoadSample NamedTuple:
    (sample_id, name, slot_index, grouping_source, name_source,
     merged_into_id, flag, exposures)
each LoadExposure NamedTuple:
    (id, filename, horizontal_position, timestamp)

`flag` is the per-sample merge/split suggestion: it is NOT a stored column —
it is produced by Phase B Task 12's `derive_sample_flags` (§9.1) over the nested
rows this function just read, THEN any flag whose sample has a non-undone
`grouping_flag_dismissed` event is suppressed to `nothing`. `derive_sample_flags`
is resolved at call time via `isdefined` — defensive graceful-degradation for
incremental build order (if this roll-up is built before Phase B Task 12 lands
in the same tree, every `flag` is `nothing` and the suppression set is still
computed — there is just nothing to suppress). See the §9.1 seam note at the
top of this plan.

Only non-retired samples (merged_into_id IS NULL) are included.
Re-derived from (load_id, slot_index) — NOT replayed from the event log.
"""
function get_loads_rollup(db::SQLite.DB, experiment_id::Integer)
    loads = Tables.rowtable(DBInterface.execute(db,
        """SELECT id AS load_id, load_index, session_id,
                  start_time, end_time, frame_count, note
           FROM loads
           WHERE experiment_id = ?
           ORDER BY load_index""", [Int(experiment_id)]))

    # Samples whose flag is dismissed by a NON-UNDONE grouping_flag_dismissed
    # event. A *dismiss* is a grouping_flag_dismissed row with NO undoes_event_id
    # (an undo is also a grouping_flag_dismissed row but CARRIES undoes_event_id,
    # so `ua.undoes_event_id IS NULL` excludes undo rows from being counted as
    # dismisses-of-their-own — Task 8b). A dismiss is "undone" iff a later event
    # carries `undoes_event_id = <dismiss id>`; the NOT EXISTS drops those.
    # (entity_type='sample', so this never collides with exposure-keyed events.)
    dismissed_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT DISTINCT ua.entity_id AS sample_id
           FROM user_actions ua
           JOIN samples s ON s.id = ua.entity_id
           WHERE ua.action = 'grouping_flag_dismissed'
             AND ua.entity_type = 'sample'
             AND ua.undoes_event_id IS NULL
             AND s.experiment_id = ?
             AND NOT EXISTS (
               SELECT 1 FROM user_actions u2
               WHERE u2.undoes_event_id = ua.id
             )""", [Int(experiment_id)]))
    dismissed = Set(Int(r.sample_id) for r in dismissed_rows)

    nested = map(loads) do ld
        samples = Tables.rowtable(DBInterface.execute(db,
            """SELECT id AS sample_id, name, slot_index,
                      grouping_source, name_source, merged_into_id
               FROM samples
               WHERE load_id = ? AND (merged_into_id IS NULL)
               ORDER BY slot_index""", [Int(ld.load_id)]))

        samples_with_exposures = map(samples) do sm
            exposures = Tables.rowtable(DBInterface.execute(db,
                """SELECT id, filename, horizontal_position, timestamp
                   FROM exposures
                   WHERE sample_id = ?
                   ORDER BY frame_no, id""", [Int(sm.sample_id)]))
            # `flag` is filled in below from derive_sample_flags; seed nothing.
            merge(sm, (flag = nothing, exposures = exposures))
        end

        merge(ld, (samples = samples_with_exposures,))
    end

    # Attach the derived per-sample flag (Phase B §9.1), resolved at call time.
    # derive_sample_flags is pure over the nested rows above and returns
    # Dict{Int sample_id => GroupingFlag (a NamedTuple) | nothing}. We pass it
    # the load rows it needs (each sample's slot_index/name + its exposures'
    # horizontal_position/timestamp/filename) — exactly the nested shape it
    # documents. Suppress any flag whose sample is in the dismissed set.
    if isdefined(@__MODULE__, :derive_sample_flags)
        flags = derive_sample_flags(nested)        # Dict{Int, Any}
        nested = map(nested) do ld
            new_samples = map(ld.samples) do sm
                f = get(flags, Int(sm.sample_id), nothing)
                Int(sm.sample_id) in dismissed && (f = nothing)
                merge(sm, (flag = f,))
            end
            merge(ld, (samples = new_samples,))
        end
    end

    return nested
end
```

> **`derive_sample_flags` argument shape (aligned with Phase B Task 12).** This call passes the **already-nested `Vector` of Load NamedTuples** (`nested`) — each Load carries `.samples`, each sample carries `.sample_id`/`.name`/`.slot_index` and `.exposures` with `.filename`/`.horizontal_position`/`.timestamp`. **Phase B Task 12 implements `derive_sample_flags(load_rows) → Dict{Int, GroupingFlag}` consuming exactly this nested shape** (confirmed against Phase B's input-contract section — it reads the roll-up field name `horizontal_position`), NOT the raw `ExposureMeta` stream from `scan_directory`. The seam is aligned; no reconciliation needed beyond the minor docstring-leaf-name sync noted at the top of this plan (Phase B Task 12 only reads `sample_id`/`slot_index`/`horizontal_position`/`timestamp`, all identically named in both).

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_ingestion_structural_events.jl
git commit -m "feat(db): get_loads_rollup for grouping-review surface"
```

---

## Task 5: Rename route (`PATCH /api/samples/{id}/name`)

Clones the `PATCH /api/samples/{id}` pattern from `routes_samples.jl:150–186`. Writes `samples.name` directly, records a `sample_renamed` event for audit + SSE, enforces never-clobber only on the field being written here (always writes, since this is an explicit rename action — spec §4 never-clobber is about scan refreshes, not user edits).

**Files:**
- Create: `packages/HimalayaUI/src/routes_grouping.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_structural_events.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "PATCH /api/samples/{id}/name — renames and records event" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

    # Drive the real registered route over in-process HTTP. with_test_server
    # binds `db` (so the handler's current_db() resolves to it) and starts a
    # server whose start_test_server! → register_routes! → register_grouping_routes!
    # registers the grouping routes (wired in Task 5 Step 4 / Task 9). Do NOT
    # call register_grouping_routes!() here — registration is global and
    # idempotent re-registration into Oxygen's router is not what we want.
    with_test_server(db) do port, base
        resp = HTTP.patch("$base/api/samples/$(s1_id)/name";
            body = JSON3.write(Dict(:name => "JC C04 (S01P01)")),
            headers = ["Content-Type"   => "application/json",
                       "X-Username"      => "alice",
                       "X-Client-Op-Id"  => "test-op-rename-1"],
            status_exception = false)
        @test resp.status == 200

        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT name, name_source FROM samples WHERE id = ?", [s1_id])))
        @test String(row.name) == "JC C04 (S01P01)"
        @test String(row.name_source) == "user"

        # A sample_renamed event was written.
        evts = Tables.rowtable(DBInterface.execute(db,
            "SELECT action FROM user_actions WHERE entity_type='sample' AND entity_id=?", [s1_id]))
        @test any(r -> String(r.action) == "sample_renamed", evts)
    end
end
```

> **Route-test pattern (matches `test_routes_samples.jl`):** there is NO `handle_grouping_route(req, db)` function and no `(req, db)` dispatch — Oxygen registers handlers via `@patch`/`@post`/`@get` macros into its global router, so routes are only reachable over HTTP. `with_test_server(db) do port, base … end` (from `test_http.jl`, available standalone via the include guard in Task 0) is the canonical driver: it `bind_db!(db)`s and `start_test_server!`s, then yields `(port, base)` for `HTTP.patch/post/get("$base/api/…"; status_exception=false)`. Pass `status_exception=false` so non-2xx responses return a `Response` instead of throwing, letting `@test resp.status == …` assert error codes. **The grouping routes must be wired into `register_routes!` before these route tests can pass — see Step 4 below, which does the wiring as part of this first route task so Tasks 5–9 all have the routes registered.**

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: FAIL (`register_grouping_routes!` / the route does not exist).

- [ ] **Step 3: Create `routes_grouping.jl` with the rename route**

```julia
# packages/HimalayaUI/src/routes_grouping.jl
using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_grouping_routes!()

    # ── Rename a sample ──────────────────────────────────────────────────────
    # PATCH /api/samples/{id}/name
    # Body: { "name": "new label" }
    # Sets samples.name = body.name and name_source = 'user'. Records a
    # sample_renamed event for audit + SSE broadcast.
    # Clones the with_idempotency + apply_event! + _enqueue_broadcast pattern
    # from routes_samples.jl:167–185.
    @patch "/api/samples/{id}/name" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if !haskey(body, :name) || !(body.name isa AbstractString)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "name is required and must be a string")))
        end
        new_name = strip(String(body.name))
        isempty(new_name) && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "name must not be blank")))

        # Resolve the experiment for the SSE payload (frontend invalidates by it).
        exp_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT experiment_id FROM samples WHERE id = ?", [id]))
        isempty(exp_rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "sample not found")))
        exp_id = Int(exp_rows[1].experiment_id)

        return with_idempotency(db, req) do
            DBInterface.execute(db,
                "UPDATE samples SET name = ?, name_source = 'user' WHERE id = ?",
                [new_name, id])
            result = apply_event!(InTransaction(), db, req;
                kind        = "sample_renamed",
                entity_type = "sample",
                entity_id   = id,
                payload     = Dict(:name => new_name, :experiment_id => exp_id))
            _enqueue_broadcast_from_result!(result, "sample_renamed", "sample", id)

            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT * FROM samples WHERE id = ?", [id]))
            isempty(rows) && return HTTP.Response(404,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "sample not found")))
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(row_to_json(rows[1])))
        end
    end

end  # register_grouping_routes!
```

> **Confirm current_db / json / row_to_json availability:** these are module-level helpers defined in `server.jl` (`current_db()`, `events.jl` (`apply_event!`/`_enqueue_broadcast_from_result!`), `json.jl` (`json(req)`, `row_to_json(row)`). They are visible to `routes_grouping.jl` because the whole package is one module (`HimalayaUI`) with each file `include`d into it — confirm by noting `routes_samples.jl` uses the same bare names. No `using`/import beyond the header line is needed.

- [ ] **Step 4: Wire the file into the module and register the routes**

Route tests reach the route only through `register_routes!` (called by `start_test_server!`), so the wiring must happen NOW — as part of this first route task — not deferred to Task 9. Two edits:

(a) In `packages/HimalayaUI/src/HimalayaUI.jl`, add the include **after the last `include("routes_*.jl")` line (currently `include("routes_resolve.jl")` at line 30) and BEFORE `include("server.jl")` (line 31)** — `server.jl`'s `register_routes!` references `register_grouping_routes!`, so the function must be defined first:

```julia
include("routes_resolve.jl")
include("routes_grouping.jl")   # ← add this line
include("server.jl")
```

(b) In `packages/HimalayaUI/src/server.jl`, add the registration call inside `register_routes!()` (the single registration site, `server.jl:49`; the existing `register_*_routes!()` calls are at lines 120–131). Append after `register_resolve_routes!()`:

```julia
    register_resolve_routes!()
    register_grouping_routes!()   # ← add this line
end
```

Both production `serve` and the test harness (`start_test_server!` → `register_routes!`, `server.jl:200`) now register the grouping routes.

- [ ] **Step 5: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/routes_grouping.jl packages/HimalayaUI/src/HimalayaUI.jl packages/HimalayaUI/src/server.jl packages/HimalayaUI/test/test_ingestion_structural_events.jl
git commit -m "feat(routes): PATCH /api/samples/{id}/name — rename + event + SSE; wire register_grouping_routes!"
```

---

## Task 6: Move route (`POST /api/exposures/{id}/move`)

Moves one exposure to a different sample within the same experiment. Emits an `exposure_moved` event (which the Task 1 dispatcher translates into `UPDATE exposures SET sample_id = ?`). Validates same-experiment constraint (spec §9.2: "within-experiment only"). Does not trigger reanalysis (spec §5: "regroup is bookkeeping only").

**Files:**
- Modify: `packages/HimalayaUI/src/routes_grouping.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_structural_events.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "POST /api/exposures/{id}/move — moves exposure + records event" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

    with_test_server(db) do port, base
        resp = HTTP.post("$base/api/exposures/$(e1_id)/move";
            body = JSON3.write(Dict(:sample_id => s2_id)),
            headers = ["Content-Type"  => "application/json",
                       "X-Username"     => "alice",
                       "X-Client-Op-Id" => "test-op-move-1"],
            status_exception = false)
        @test resp.status == 200

        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id FROM exposures WHERE id = ?", [e1_id])))
        @test Int(row.sample_id) == s2_id

        evts = Tables.rowtable(DBInterface.execute(db,
            "SELECT action FROM user_actions WHERE entity_type='exposure' AND entity_id=?", [e1_id]))
        @test any(r -> String(r.action) == "exposure_moved", evts)
    end
end

@testset "POST /api/exposures/{id}/move — rejects cross-experiment move" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)
    # Create a sample in a different experiment.
    exp2_id = HimalayaUI.create_experiment!(db; path="/d2", data_dir="/d2", analysis_dir="/a2")
    s_other = HimalayaUI.create_sample!(db; experiment_id=exp2_id, name="Other")

    with_test_server(db) do port, base
        resp = HTTP.post("$base/api/exposures/$(e1_id)/move";
            body = JSON3.write(Dict(:sample_id => s_other)),
            headers = ["Content-Type"  => "application/json",
                       "X-Username"     => "alice",
                       "X-Client-Op-Id" => "test-op-move-2"],
            status_exception = false)
        @test resp.status == 422   # cross-experiment forbidden
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: FAIL (route returns 404 — the move path is not yet defined inside `register_grouping_routes!`, though the function is now wired from Task 5).

- [ ] **Step 3: Add the move route to `routes_grouping.jl`**

Inside `register_grouping_routes!()`, after the rename route:

```julia
    # ── Move one exposure to a different sample ───────────────────────────────
    # POST /api/exposures/{id}/move
    # Body: { "sample_id": Int }
    # Within-experiment only (spec §9.2 / §5).
    @post "/api/exposures/{id}/move" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if !haskey(body, :sample_id)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "sample_id is required")))
        end
        dest_sample_id = Int(body.sample_id)

        # Validate same-experiment constraint; also capture the source sample id
        # and the experiment id for the SSE payload.
        src_rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT e.sample_id AS from_sample_id,
                      e.experiment_id AS src_exp,
                      s.experiment_id AS dst_exp
               FROM exposures e
               JOIN samples s ON s.id = ?
               WHERE e.id = ?""", [dest_sample_id, id]))
        if isempty(src_rows)
            return HTTP.Response(404,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "exposure or destination sample not found")))
        end
        src_row = first(src_rows)
        if !ismissing(src_row.src_exp) && !ismissing(src_row.dst_exp) &&
                Int(src_row.src_exp) != Int(src_row.dst_exp)
            return HTTP.Response(422,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error =>
                    "cross-experiment moves are not allowed; move would change the resolved analysis_dir")))
        end
        from_sample_id = ismissing(src_row.from_sample_id) ? nothing : Int(src_row.from_sample_id)
        exp_id = ismissing(src_row.dst_exp) ? nothing : Int(src_row.dst_exp)

        return with_idempotency(db, req) do
            result = apply_event!(InTransaction(), db, req;
                kind        = "exposure_moved",
                entity_type = "exposure",
                entity_id   = id,
                payload     = Dict(:sample_id => dest_sample_id,
                                   :from_sample_id => from_sample_id,
                                   :experiment_id => exp_id))
            _enqueue_broadcast_from_result!(result, "exposure_moved", "exposure", id)

            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT * FROM exposures WHERE id = ?", [id]))
            isempty(rows) && return HTTP.Response(404,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "exposure not found")))
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(row_to_json(rows[1])))
        end
    end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_grouping.jl packages/HimalayaUI/test/test_ingestion_structural_events.jl
git commit -m "feat(routes): POST /api/exposures/{id}/move — move + cross-experiment guard"
```

---

## Task 7: Merge route (`POST /api/samples/{id}/merge`)

The most complex task. The whole route body is one `with_idempotency` block. Inside it: (a) pre-merge dedup writes (drop loser's `series_samples` memberships in series the survivor already belongs to; drop loser's `sample_tags` where survivor holds the same key); (b) re-point all child FK tables to the survivor; (c) retire the loser; (d) if a rename payload is provided, record a `sample_renamed` event; (e) enqueue per-event broadcasts. SSE cache handling on the frontend is invalidate-only (refetch `loads(id)` — spec §9.3).

**Files:**
- Modify: `packages/HimalayaUI/src/routes_grouping.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_structural_events.jl`

- [ ] **Step 1: Write the failing tests**

```julia
@testset "POST /api/samples/{id}/merge — basic re-point + retire" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

    with_test_server(db) do port, base
        # Merge s2 (loser) into s1 (survivor).
        resp = HTTP.post("$base/api/samples/$(s2_id)/merge";
            body = JSON3.write(Dict(:survivor_id => s1_id)),
            headers = ["Content-Type"  => "application/json",
                       "X-Username"     => "alice",
                       "X-Client-Op-Id" => "test-op-merge-1"],
            status_exception = false)
        @test resp.status == 200

        # e2 (was in s2) now belongs to s1.
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id FROM exposures WHERE id = ?", [e2_id])))
        @test Int(row.sample_id) == s1_id

        # s2 is retired, not deleted.
        s2_row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT merged_into_id FROM samples WHERE id = ?", [s2_id])))
        @test Int(s2_row.merged_into_id) == s1_id

        # exposure_moved events were recorded for the moved exposures.
        evts = Tables.rowtable(DBInterface.execute(db,
            "SELECT action FROM user_actions WHERE entity_type='exposure' AND entity_id=?", [e2_id]))
        @test any(r -> String(r.action) == "exposure_moved", evts)
    end
end

@testset "POST /api/samples/{id}/merge — series_samples dedup: drops loser membership" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

    # Create a series containing both s1 and s2.
    series_id = DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO series DEFAULT VALUES"))
    DBInterface.execute(db,
        "INSERT INTO series_samples (series_id, sample_id, position) VALUES (?, ?, 1)",
        [series_id, s1_id])
    DBInterface.execute(db,
        "INSERT INTO series_samples (series_id, sample_id, position) VALUES (?, ?, 2)",
        [series_id, s2_id])

    with_test_server(db) do port, base
        resp = HTTP.post("$base/api/samples/$(s2_id)/merge";
            body = JSON3.write(Dict(:survivor_id => s1_id)),
            headers = ["Content-Type"  => "application/json",
                       "X-Username"     => "alice",
                       "X-Client-Op-Id" => "test-op-merge-dedup-ss"],
            status_exception = false)
        @test resp.status == 200

        # s2's membership in the series where s1 already appears must be dropped.
        surviving = Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id FROM series_samples WHERE series_id = ?", [series_id]))
        remaining_ids = Set(Int(r.sample_id) for r in surviving)
        @test s2_id ∉ remaining_ids   # loser's duplicate membership dropped
        @test s1_id in remaining_ids  # survivor's membership intact
    end
end

@testset "POST /api/samples/{id}/merge — sample_tags dedup: survivor wins" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

    # Both samples have a tag with key "condition"; values differ.
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value) VALUES (?, 'condition', 'acid')", [s1_id])
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value) VALUES (?, 'condition', 'base')", [s2_id])
    # s2 has a unique tag that should survive the merge.
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value) VALUES (?, 'ph', '8.0')", [s2_id])

    with_test_server(db) do port, base
        resp = HTTP.post("$base/api/samples/$(s2_id)/merge";
            body = JSON3.write(Dict(:survivor_id => s1_id)),
            headers = ["Content-Type"  => "application/json",
                       "X-Username"     => "alice",
                       "X-Client-Op-Id" => "test-op-merge-tags"],
            status_exception = false)
        @test resp.status == 200

        tags = Tables.rowtable(DBInterface.execute(db,
            "SELECT key, value FROM sample_tags WHERE sample_id = ?", [s1_id]))
        tag_map = Dict(String(t.key) => String(t.value) for t in tags)
        @test tag_map["condition"] == "acid"  # survivor wins on collision
        @test tag_map["ph"] == "8.0"          # loser's unique tag transferred
    end
end
```

> **`sample_tags` columns:** the dedup tests insert `(sample_id, key, value)` without a `source` column; the live `sample_tags` table has a `source` column with a default of `'manual'` (see `routes_samples.jl:223`). Omitting it relies on that default — confirm the column is nullable-or-defaulted before writing these inserts; if not, add `, source` / `'manual'`.

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: FAIL (route returns 404 — `/api/samples/{id}/merge` path not yet defined inside the already-wired `register_grouping_routes!`).

- [ ] **Step 3: Add the merge route to `routes_grouping.jl`**

Inside `register_grouping_routes!()`, after the move route:

```julia
    # ── Merge samples: loser → survivor ──────────────────────────────────────
    # POST /api/samples/{id}/merge
    # Body: { "survivor_id": Int, "name"?: String }
    # id = loser. All child rows re-pointed, then loser retired.
    # Whole body is one with_idempotency block (spec §9.3).
    # Undo is session-local (frontend useUndoStack). No server-side undo API.
    # SSE: applyRemoteToCache must refetch loads(id) on merge events (invalidate-only).
    @post "/api/samples/{id}/merge" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if !haskey(body, :survivor_id)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "survivor_id is required")))
        end
        loser_id    = id
        survivor_id = Int(body.survivor_id)
        loser_id == survivor_id && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "loser and survivor must be different")))

        # Validate both samples exist and belong to the same experiment.
        both = Tables.rowtable(DBInterface.execute(db,
            """SELECT id, experiment_id FROM samples WHERE id IN (?, ?)""",
            [loser_id, survivor_id]))
        length(both) != 2 && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "one or both samples not found")))
        if Int(both[1].experiment_id) != Int(both[2].experiment_id)
            return HTTP.Response(422,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "cannot merge samples from different experiments")))
        end
        exp_id = Int(both[1].experiment_id)

        return with_idempotency(db, req) do
            broadcast_results = []

            # (i) series_samples dedup — drop loser's membership in any series
            # where the survivor already appears (spec §9.3). Done BEFORE re-point
            # so the surviving UNIQUE(series_id, position) isn't invalidated first.
            DBInterface.execute(db,
                """DELETE FROM series_samples
                   WHERE sample_id = ?
                     AND series_id IN (
                       SELECT series_id FROM series_samples WHERE sample_id = ?
                     )""",
                [loser_id, survivor_id])

            # (ii) sample_tags dedup — drop loser's tag where survivor holds same key;
            # survivor wins on collision (spec §9.3).
            DBInterface.execute(db,
                """DELETE FROM sample_tags
                   WHERE sample_id = ?
                     AND key IN (
                       SELECT key FROM sample_tags WHERE sample_id = ?
                     )""",
                [loser_id, survivor_id])

            # (iii) re-point all child FK tables to survivor.
            # exposures: emit one exposure_moved event per exposure moved.
            loser_exposures = Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM exposures WHERE sample_id = ?", [loser_id]))
            for ex_row in loser_exposures
                ex_id = Int(ex_row.id)
                # Use a derived op-id so each per-exposure sub-event has its own
                # idempotency key within the enclosing with_idempotency block.
                sub_req = HTTP.Request("POST", req.target,
                    vcat(req.headers,
                         ["X-Client-Op-Id" => "merge-$(loser_id)-$(survivor_id)-move-$(ex_id)"]),
                    UInt8[])
                result = apply_event!(InTransaction(), db, sub_req;
                    kind        = "exposure_moved",
                    entity_type = "exposure",
                    entity_id   = ex_id,
                    payload     = Dict(:sample_id => survivor_id,
                                       :from_sample_id => loser_id,
                                       :experiment_id => exp_id))
                push!(broadcast_results, (result, "exposure_moved", "exposure", ex_id))
            end

            # series_samples re-point (survivors of the dedup step).
            DBInterface.execute(db,
                "UPDATE series_samples SET sample_id = ? WHERE sample_id = ?",
                [survivor_id, loser_id])

            # sample_tags re-point (survivors of the dedup step).
            DBInterface.execute(db,
                "UPDATE sample_tags SET sample_id = ? WHERE sample_id = ?",
                [survivor_id, loser_id])

            # sample_messages re-point (no uniqueness constraint — plain UPDATE).
            DBInterface.execute(db,
                "UPDATE sample_messages SET sample_id = ? WHERE sample_id = ?",
                [survivor_id, loser_id])

            # (iv) Retire the loser (sets merged_into_id; does not hard-delete).
            retire_sample!(db, loser_id; merged_into_id=survivor_id)

            # (v) Optional rename of survivor.
            if haskey(body, :name) && body.name isa AbstractString
                new_name = strip(String(body.name))
                if !isempty(new_name)
                    DBInterface.execute(db,
                        "UPDATE samples SET name = ?, name_source = 'user' WHERE id = ?",
                        [new_name, survivor_id])
                    rename_result = apply_event!(InTransaction(), db, req;
                        kind        = "sample_renamed",
                        entity_type = "sample",
                        entity_id   = survivor_id,
                        payload     = Dict(:name => new_name, :experiment_id => exp_id))
                    push!(broadcast_results,
                        (rename_result, "sample_renamed", "sample", survivor_id))
                end
            end

            # Enqueue all broadcasts after the transaction commits.
            for (result, kind, etype, eid) in broadcast_results
                _enqueue_broadcast_from_result!(result, kind, etype, eid)
            end

            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:loser_id => loser_id, :survivor_id => survivor_id)))
        end
    end
```

> **`X-Client-Op-Id` for sub-events:** the sub-request headers add a per-exposure op-id derived from the merge pair + exposure id. This ensures each `apply_event!(InTransaction(), …)` call has its own `(client_op_id, action, entity_id)` UNIQUE slot so idempotent replay works correctly for each individual exposure move. The enclosing `with_idempotency` block caches the entire HTTP response.

> **Verify `_enqueue_broadcast_from_result!` visibility:** it is currently defined in `events.jl` (`events.jl:241`). Confirm it is exported or accessible from `routes_grouping.jl` by the same module-level inclusion used by `routes_samples.jl`. If it is not exported, either export it or prefix it with `HimalayaUI.` in the test.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_grouping.jl packages/HimalayaUI/test/test_ingestion_structural_events.jl
git commit -m "feat(routes): POST /api/samples/{id}/merge — re-point + dedup + retire"
```

---

## Task 8: Split route (`POST /api/samples/{id}/split`)

Split creates a new sample from a subset of the source sample's exposures (the user specifies the exposures to move). Unlike merge, there is no FK re-point of non-exposure tables (the original sample stays; the new one is fresh). Follows the `comparison_created` / `series_created` precedent: the route mints the new sample id and passes it as `entity_id` in a `sample_created` event (payload `{ experiment_id }`). **`split` is the sole emitter of `sample_created`** (the contract carries a dedicated frontend `sample_created` arm that does invalidate-only refresh of `samples`/`corpusSamples`/`loads`); **`merge` does NOT emit `sample_created`** — it only fans out per-exposure `exposure_moved` (+ optional `sample_renamed`). The split route emits, in order: `sample_created` (new sample), one `exposure_moved` per moved exposure (source = the original `id`), then `sample_split` (entity = source). Undo is session-local.

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl` (add `sample_created` dispatcher branch)
- Modify: `packages/HimalayaUI/src/routes_grouping.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_structural_events.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "POST /api/samples/{id}/split — creates new sample + moves exposures" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)
    # Add a second exposure to s1 so we have two to split.
    e3_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id,
        sample_id=s1_id, filename="f03.tif")

    with_test_server(db) do port, base
        # Split: move e3 out of s1 into a new sample.
        resp = HTTP.post("$base/api/samples/$(s1_id)/split";
            body = JSON3.write(Dict(:exposure_ids => [e3_id], :name => "HA85 (S01P01b)")),
            headers = ["Content-Type"  => "application/json",
                       "X-Username"     => "alice",
                       "X-Client-Op-Id" => "test-op-split-1"],
            status_exception = false)
        @test resp.status == 201

        resp_body = JSON3.read(String(resp.body))
        new_sample_id = Int(resp_body.new_sample_id)

        # e3 now belongs to the new sample.
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id FROM exposures WHERE id = ?", [e3_id])))
        @test Int(row.sample_id) == new_sample_id

        # e1 still belongs to s1.
        row2 = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id FROM exposures WHERE id = ?", [e1_id])))
        @test Int(row2.sample_id) == s1_id

        # New sample has the given name.
        new_row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT name, experiment_id FROM samples WHERE id = ?", [new_sample_id])))
        @test String(new_row.name) == "HA85 (S01P01b)"
        @test Int(new_row.experiment_id) == exp_id

        # exposure_moved event recorded for each moved exposure.
        evts = Tables.rowtable(DBInterface.execute(db,
            "SELECT action FROM user_actions WHERE entity_type='exposure' AND entity_id=?", [e3_id]))
        @test any(r -> String(r.action) == "exposure_moved", evts)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: FAIL (route returns 404 — `/api/samples/{id}/split` path not yet defined inside the already-wired `register_grouping_routes!`).

- [ ] **Step 3: Add intent-marker no-op dispatcher branches in `events.jl` (optional)**

Read `events.jl:400–424` (the M2.1 trivial-route no-op block) and the default no-op at `events.jl:516–517`. After the `kind == "analyze_run"` no-op line, add:

```julia
    # Structural grouping edits — dispatcher no-ops for kinds whose view
    # writes happen directly in the route (per the update_sample precedent).
    kind == "sample_created" && return nothing   # split mints the sample row in-route
    kind == "sample_split"   && return nothing   # split emits per-exposure exposure_moved
    # grouping_flag_dismissed is durable + undoable but view-less: the merge/split
    # flag is never stored — get_loads_rollup suppresses it at READ time by joining
    # non-undone dismiss events (Task 4). So the dispatcher is a deliberate no-op,
    # and undo is just another event carrying undoes_event_id (Task 8b).
    kind == "grouping_flag_dismissed" && return nothing
```

> **There is NO `sample_merged` kind.** Spec §9.3 is explicit: a merge fans out as per-exposure `exposure_moved` (+ optional `sample_renamed`) and retires the loser via `merged_into_id`; the frontend has no `sample_merged` arm. So do **not** add a `sample_merged` no-op — it would be a dead branch for a kind that is never emitted.
>
> **These branches are NOT required for correctness.** The dispatcher's default is already a silent no-op (`events.jl:516–517` returns `nothing` for any unrecognized kind), so `rebuild_views_from_log!` does NOT throw on these kinds even without the branches. They are intent markers — they document that `sample_created`/`sample_split` are deliberately view-less audit-only events (the route writes the sample row and the per-exposure `exposure_moved` events carry the only replayable view state), matching the existing retired-kind and trivial-route guards (`events.jl:363,405–411`). Only `exposure_moved` (Task 1), `sample_created`, and `sample_split` (split route), plus `sample_renamed` (Task 1 marker) and `grouping_flag_dismissed` (Task 8b marker) correspond to events this phase actually emits.

- [ ] **Step 4: Add the split route to `routes_grouping.jl`**

Inside `register_grouping_routes!()`, after the merge route:

```julia
    # ── Split a sample: move a subset of exposures to a new sample ───────────
    # POST /api/samples/{id}/split
    # Body: { "exposure_ids": [Int], "name"?: String }
    # Creates a new sample in the same experiment+load (inherits slot_index = NULL),
    # moves the specified exposures to it, records an audit event per move.
    # Returns 201 with { new_sample_id: Int }.
    @post "/api/samples/{id}/split" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if !haskey(body, :exposure_ids) || !(body.exposure_ids isa AbstractVector)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "exposure_ids (array) is required")))
        end
        exposure_ids = Int.(body.exposure_ids)
        isempty(exposure_ids) && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure_ids must not be empty")))

        # Validate source sample exists.
        src_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT experiment_id, load_id FROM samples WHERE id = ?", [id]))
        isempty(src_rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "sample not found")))
        src = first(src_rows)
        exp_id  = Int(src.experiment_id)
        load_id = ismissing(src.load_id) ? nothing : Int(src.load_id)

        # Validate all exposures belong to the source sample (not a foreign sample).
        # The IN-list is string-interpolated rather than parameterized. This is
        # injection-SAFE: `exposure_ids = Int.(body.exposure_ids)` above coerces
        # every element to Int (a non-Int element throws before we get here), so
        # the interpolated text is a comma-joined list of integer literals — no
        # user string ever reaches the SQL. The `?` placeholder still binds
        # `sample_id` via the parameter array.
        owned = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM exposures WHERE sample_id = ? AND id IN ($(join(exposure_ids, ',')))",
            [id]))
        length(owned) != length(exposure_ids) && return HTTP.Response(422,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "one or more exposure_ids do not belong to this sample")))

        new_name = (haskey(body, :name) && body.name isa AbstractString) ?
            strip(String(body.name)) : nothing

        return with_idempotency(db, req) do
            broadcast_results = []

            # Mint the new sample row (route owns the id — comparison_created precedent).
            new_sample_id = create_sample!(db;
                experiment_id    = exp_id,
                load_id          = load_id,
                name             = something(new_name, "Split from sample $(id)"),
                grouping_source  = "manual",
                name_source      = isnothing(new_name) ? "auto" : "user")

            # Record a sample_created audit event for the new sample (no view
            # write — no-op dispatcher). Canonical payload is `{ experiment_id }`
            # only (spec §9.3): the frontend's sample_created arm does an
            # invalidate-only refresh keyed on experiment_id and does not read
            # any other field. (Split is the emitter; merge does NOT emit it.)
            created_result = apply_event!(InTransaction(), db, req;
                kind        = "sample_created",
                entity_type = "sample",
                entity_id   = new_sample_id,
                payload     = Dict(:experiment_id => exp_id))
            push!(broadcast_results,
                (created_result, "sample_created", "sample", new_sample_id))

            # Move each exposure (source = the original sample `id`).
            for ex_id in exposure_ids
                sub_req = HTTP.Request("POST", req.target,
                    vcat(req.headers,
                         ["X-Client-Op-Id" => "split-$(id)-$(new_sample_id)-move-$(ex_id)"]),
                    UInt8[])
                result = apply_event!(InTransaction(), db, sub_req;
                    kind        = "exposure_moved",
                    entity_type = "exposure",
                    entity_id   = ex_id,
                    payload     = Dict(:sample_id => new_sample_id,
                                       :from_sample_id => id,
                                       :experiment_id => exp_id))
                push!(broadcast_results, (result, "exposure_moved", "exposure", ex_id))
            end

            # Record a sample_split audit event on the source.
            split_result = apply_event!(InTransaction(), db, req;
                kind        = "sample_split",
                entity_type = "sample",
                entity_id   = id,
                payload     = Dict(:new_sample_id => new_sample_id,
                                   :exposure_ids  => exposure_ids,
                                   :experiment_id => exp_id))
            push!(broadcast_results, (split_result, "sample_split", "sample", id))

            for (result, kind, etype, eid) in broadcast_results
                _enqueue_broadcast_from_result!(result, kind, etype, eid)
            end

            HTTP.Response(201,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:new_sample_id => new_sample_id)))
        end
    end
```

> **`something` usage:** `something(x, default)` is available in Julia Base. If `new_name` is `nothing`, it passes `nothing` to the first arg of `something` and falls back to the default string. Verify this compiles by reading `Base.something`.

- [ ] **Step 5: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/events.jl packages/HimalayaUI/src/routes_grouping.jl packages/HimalayaUI/test/test_ingestion_structural_events.jl
git commit -m "feat(routes): POST /api/samples/{id}/split + sample_created/split dispatcher no-ops"
```

---

## Task 8b: `grouping_flag_dismissed` durable event — dismiss + undo routes

"Keep separate" on a merge/split suggestion writes a durable, undoable `grouping_flag_dismissed` event keyed to the flagged sample (spec §9.3). There is no view write — `get_loads_rollup` (Task 4) suppresses the flag at read time by joining non-undone dismiss events. The dispatcher no-op was added in Task 8 Step 3. This task adds the route + its inverse (undo).

Two routes:
- `POST /api/samples/{id}/dismiss-flag` — body `{ flag_kind: "merge" | "split", merge_with_sample_id? }`. Records `grouping_flag_dismissed` (entity = sample `id`), payload `{ flag_kind, merge_with_sample_id?, experiment_id }`. After this, the roll-up suppresses that sample's flag.
- `POST /api/samples/{id}/dismiss-flag/undo` — records a `grouping_flag_dismissed` event whose `undoes_event_id` points at the most-recent non-undone dismiss for this sample, which un-suppresses the flag (the Task 4 `NOT EXISTS` join now sees a surviving undo). This mirrors the `peak_unexcluded` "find the prior event to record as undoes" precedent at `routes_peaks.jl:296–313`.

> **Why a second event, not a DELETE.** The event log is append-only; undo is modeled by a new event stamped with `undoes_event_id` (the `peak_excluded` → `peak_unexcluded` precedent), never by deleting the original. Task 4's suppression query is `NOT EXISTS (… u2.undoes_event_id = ua.id)`, so an undo event cancels the dismiss it points at. A redo would be another dismiss; the latest non-undone dismiss wins.

**Files:**
- Modify: `packages/HimalayaUI/src/routes_grouping.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_structural_events.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "POST /api/samples/{id}/dismiss-flag — suppresses flag; undo re-shows" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

    with_test_server(db) do port, base
        # Dismiss a merge flag on s1.
        resp = HTTP.post("$base/api/samples/$(s1_id)/dismiss-flag";
            body = JSON3.write(Dict(:flag_kind => "merge", :merge_with_sample_id => s2_id)),
            headers = ["Content-Type"  => "application/json",
                       "X-Username"     => "alice",
                       "X-Client-Op-Id" => "test-op-dismiss-1"],
            status_exception = false)
        @test resp.status == 200

        # A non-undone dismiss event exists on s1.
        evts = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, action FROM user_actions WHERE entity_type='sample' AND entity_id=?", [s1_id]))
        @test any(r -> String(r.action) == "grouping_flag_dismissed", evts)

        # The rollup suppresses s1's flag.
        s1 = first(filter(s -> s.sample_id == s1_id,
                          first(HimalayaUI.get_loads_rollup(db, exp_id)).samples))
        @test s1.flag === nothing

        # Undo the dismiss.
        resp2 = HTTP.post("$base/api/samples/$(s1_id)/dismiss-flag/undo";
            body = JSON3.write(Dict()),
            headers = ["Content-Type"  => "application/json",
                       "X-Username"     => "alice",
                       "X-Client-Op-Id" => "test-op-dismiss-undo-1"],
            status_exception = false)
        @test resp2.status == 200

        # The dismiss is now undone — the suppression set no longer contains s1.
        # (s1.flag is whatever derive_sample_flags returns; absent Phase B it is
        # nothing, but it is no longer suppressed-by-dismiss — assert via the
        # suppression query directly so the test is independent of Phase B.)
        suppressed = Tables.rowtable(DBInterface.execute(db,
            """SELECT 1 FROM user_actions ua
               WHERE ua.action='grouping_flag_dismissed' AND ua.entity_type='sample'
                 AND ua.entity_id=? AND ua.undoes_event_id IS NULL
                 AND NOT EXISTS (SELECT 1 FROM user_actions u2 WHERE u2.undoes_event_id=ua.id)""",
            [s1_id]))
        @test isempty(suppressed)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: FAIL (route returns 404 — `/api/samples/{id}/dismiss-flag` not yet defined inside the already-wired `register_grouping_routes!`).

- [ ] **Step 3: Add the dismiss + undo routes to `routes_grouping.jl`**

Inside `register_grouping_routes!()`, after the split route:

```julia
    # ── Dismiss a grouping flag ("Keep separate") ────────────────────────────
    # POST /api/samples/{id}/dismiss-flag
    # Body: { "flag_kind": "merge"|"split", "merge_with_sample_id"?: Int }
    # Durable, undoable. No view write — get_loads_rollup suppresses the flag
    # at read time by joining non-undone dismiss events (Task 4).
    @post "/api/samples/{id}/dismiss-flag" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        flag_kind = haskey(body, :flag_kind) && body.flag_kind isa AbstractString ?
                    String(body.flag_kind) : nothing
        if flag_kind === nothing || !(flag_kind in ("merge", "split"))
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "flag_kind must be \"merge\" or \"split\"")))
        end

        exp_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT experiment_id FROM samples WHERE id = ?", [id]))
        isempty(exp_rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "sample not found")))
        exp_id = Int(exp_rows[1].experiment_id)

        payload = Dict{Symbol, Any}(:flag_kind => flag_kind, :experiment_id => exp_id)
        if haskey(body, :merge_with_sample_id) && body.merge_with_sample_id !== nothing
            payload[:merge_with_sample_id] = Int(body.merge_with_sample_id)
        end

        return with_idempotency(db, req) do
            result = apply_event!(InTransaction(), db, req;
                kind        = "grouping_flag_dismissed",
                entity_type = "sample",
                entity_id   = id,
                payload     = payload)
            _enqueue_broadcast_from_result!(result, "grouping_flag_dismissed", "sample", id)
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:sample_id => id, :flag_kind => flag_kind)))
        end
    end

    # ── Undo a grouping-flag dismiss ─────────────────────────────────────────
    # POST /api/samples/{id}/dismiss-flag/undo
    # Records a grouping_flag_dismissed event stamped with undoes_event_id
    # pointing at the latest non-undone dismiss for this sample, un-suppressing
    # the flag. Mirrors the peak_excluded→peak_unexcluded undo precedent
    # (routes_peaks.jl:296–313).
    @post "/api/samples/{id}/dismiss-flag/undo" function(req::HTTP.Request, id::Int)
        db = current_db()

        exp_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT experiment_id FROM samples WHERE id = ?", [id]))
        isempty(exp_rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "sample not found")))
        exp_id = Int(exp_rows[1].experiment_id)

        # Find the latest non-undone dismiss to reverse.
        prior = Tables.rowtable(DBInterface.execute(db,
            """SELECT id FROM user_actions ua
               WHERE ua.action = 'grouping_flag_dismissed'
                 AND ua.entity_type = 'sample' AND ua.entity_id = ?
                 AND NOT EXISTS (SELECT 1 FROM user_actions u2 WHERE u2.undoes_event_id = ua.id)
               ORDER BY ua.id DESC LIMIT 1""", [id]))
        isempty(prior) && return HTTP.Response(409,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "no active grouping-flag dismiss to undo")))
        undoes = Int(prior[1].id)

        return with_idempotency(db, req) do
            result = apply_event!(InTransaction(), db, req;
                kind            = "grouping_flag_dismissed",
                entity_type     = "sample",
                entity_id       = id,
                payload         = Dict(:experiment_id => exp_id, :undo => true),
                undoes_event_id = undoes)
            _enqueue_broadcast_from_result!(result, "grouping_flag_dismissed", "sample", id)
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:sample_id => id, :undone_event_id => undoes)))
        end
    end
```

> **The undo event is itself a `grouping_flag_dismissed` row** (kind reused, distinguished only by `undoes_event_id`). Task 4's suppression `NOT EXISTS` join keys off `u2.undoes_event_id = ua.id`, so the undo row both (a) cancels the dismiss it points at and (b) is itself never treated as an active dismiss because the suppression `SELECT` only joins `samples` on `ua.entity_id` and counts *any* dismiss with no surviving undo — and the undo row, being a dismiss-kind with `undoes_event_id` set, WOULD otherwise count as an active dismiss. **To avoid the undo row re-suppressing the flag, Task 4's suppression query must exclude rows that carry `undoes_event_id`** (an undo is not a dismiss-of-its-own). Confirm Task 4's query has `AND ua.undoes_event_id IS NULL` — see the note appended to Task 4 below.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_grouping.jl packages/HimalayaUI/test/test_ingestion_structural_events.jl
git commit -m "feat(routes): grouping_flag_dismissed dismiss + undo routes (durable, undoable)"
```

---

## Task 9: Loads roll-up GET endpoint

Add the loads roll-up GET endpoint to the already-wired `register_grouping_routes!`.

> **Wiring already done in Task 5.** `include("routes_grouping.jl")` (HimalayaUI.jl, before `include("server.jl")`) and the `register_grouping_routes!()` call inside `register_routes!()` (server.jl) were both added in Task 5 Step 4 so the route tests in Tasks 5–8 could reach the routes over HTTP. Task 9 therefore makes NO new wiring change — it only adds one more `@get` route inside the existing function and tests it.

**Files:**
- Modify: `packages/HimalayaUI/src/routes_grouping.jl` (add loads endpoint inside `register_grouping_routes!`)
- Test: `packages/HimalayaUI/test/test_ingestion_structural_events.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "GET /api/experiments/{id}/loads — returns load tree" begin
    db = fresh_db()
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

    with_test_server(db) do port, base
        resp = HTTP.get("$base/api/experiments/$(exp_id)/loads";
            status_exception = false)
        @test resp.status == 200

        data = JSON3.read(String(resp.body))
        @test length(data) == 1
        @test Int(data[1].load_id) == load_id
        @test haskey(data[1], :session_id)
        @test haskey(data[1], :note)
        @test length(data[1].samples) == 2
        sm = data[1].samples[1]
        @test haskey(sm, :flag)              # GroupingFlag JSON object or null
        ex = sm.exposures[1]
        @test haskey(ex, :id)                # exposure leaf uses `id`, not `exposure_id`
        @test haskey(ex, :horizontal_position)
        @test haskey(ex, :timestamp)
    end
end
```

> **No `register_grouping_routes!()` call in the test:** the route is registered by `start_test_server!` → `register_routes!` (wired in Task 5). This GET path does not exist yet, so before Step 3 the SPA catch-all (`server.jl:63` `@get "/**"`) returns a 404 for unregistered `api/` paths — the test fails on `resp.status == 200`.

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Expected: FAIL (the `/api/experiments/{id}/loads` `@get` route is not yet defined; the `api/`-guarded SPA catch-all returns 404).

- [ ] **Step 3: Add the loads endpoint to `routes_grouping.jl`**

Add inside `register_grouping_routes!()`, before the closing `end`:

```julia
    # ── Loads roll-up for grouping-review surface ─────────────────────────────
    # GET /api/experiments/{id}/loads
    # Returns the Load ▸ Sample ▸ Exposures tree re-derived from (load_id, slot_index).
    # Distinct from GET /api/experiments/{id}/samples (the flat corpus query).
    @get "/api/experiments/{id}/loads" function(req::HTTP.Request, id::Int)
        db = current_db()
        rollup = get_loads_rollup(db, id)
        # Serialize: NamedTuples → JSON-friendly Dict tree.
        # Field names here are the §8.8 contract verbatim (Load / LoadSample /
        # LoadExposure). `flag` serializes the GroupingFlag NamedTuple (or null);
        # JSON3 emits a NamedTuple as a JSON object, so a
        # (kind="merge", merge_with_sample_id=…, merge_with_label=…) flag lands
        # as {"kind":"merge",…} and `nothing` lands as JSON null.
        out = map(rollup) do ld
            Dict(
                :load_id     => ld.load_id,
                :load_index  => ld.load_index,
                :session_id  => ld.session_id,
                :start_time  => ld.start_time,
                :end_time    => ld.end_time,
                :frame_count => ld.frame_count,
                :note        => ld.note,
                :samples     => map(ld.samples) do sm
                    Dict(
                        :sample_id       => sm.sample_id,
                        :name            => sm.name,
                        :slot_index      => sm.slot_index,
                        :grouping_source => sm.grouping_source,
                        :name_source     => sm.name_source,
                        :merged_into_id  => sm.merged_into_id,
                        :flag            => sm.flag,
                        :exposures       => map(sm.exposures) do ex
                            Dict(
                                :id                  => ex.id,
                                :filename            => ex.filename,
                                :horizontal_position => ex.horizontal_position,
                                :timestamp           => ex.timestamp)
                        end)
                end)
        end
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
    end
```

- [ ] **Step 4: Run test + load check**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_structural_events.jl`
Then: `julia --project=packages/HimalayaUI -e 'using HimalayaUI'`
Expected: PASS + clean load (the include + `register_grouping_routes!()` wiring from Task 5 is already in place; this confirms the new `@get` route loads cleanly).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_grouping.jl packages/HimalayaUI/test/test_ingestion_structural_events.jl
git commit -m "feat(routes): GET /api/experiments/{id}/loads grouping roll-up"
```

---

## Task 10: Full backend suite + regression check

Prove this phase adds no regressions against the full `HimalayaUI` test suite.

**Files:**
- Test only.

- [ ] **Step 1: Run the full backend suite**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test-phaseD.out 2>&1
tail -40 /tmp/jl-test-phaseD.out
```

Expected: all green. Common failure modes to investigate:

- `test_events.jl` round-trip property tests that fold every event kind through `rebuild_views_from_log!` — the new kinds (`exposure_moved`, `sample_renamed`, `sample_created`, `sample_split`, `grouping_flag_dismissed`) must round-trip. They do: `exposure_moved` is an idempotent `UPDATE`, the rest are no-ops (the flag suppression is a read-time derivation in `get_loads_rollup`, never a view write, so there is nothing to re-fold), and the dispatcher's default is itself a silent no-op (`events.jl:516–517`), so an unknown kind never throws. There is no `@warn "unknown event kind"` path to chase.
- `test_routes_samples.jl` tests that reference the `PATCH /api/samples/{id}` route — the rename route is a NEW path (`/api/samples/{id}/name`), not a replacement; the old PATCH must still work. Confirm no collision.
- Any test that hard-codes the `update_view_for_event!` dispatch table — the new branches add lines but do not remove any.

- [ ] **Step 2: Confirm the new routes don't shadow existing ones**

```bash
grep -rn "api/samples/{id}" packages/HimalayaUI/src/ | grep -v "test"
```

Expected: the new `/api/samples/{id}/name` and `/api/samples/{id}/merge` and `/api/samples/{id}/split` paths are distinct from the existing `/api/samples/{id}` PATCH and `/api/samples/{id}/tags` POST. No shadowing.

- [ ] **Step 3: Commit if clean**

If the suite is green, commit a final "Phase D complete" marker:

```bash
git add packages/HimalayaUI/test/test_ingestion_structural_events.jl
git commit -m "test(phase-d): full suite green — structural-edit events complete"
```

---

## Self-Review

### Spec coverage (§9.3 + §9.2 requirements → task)

| Spec requirement | Task |
|---|---|
| `exposure_moved` — single-entity event, payload `{ sample_id, from_sample_id, experiment_id }`, dispatcher writes `exposures.sample_id` from `sample_id`, replay-idempotent | Task 1 |
| `sample_renamed` — single-entity event, payload `{ name, experiment_id }`, durable audit row + SSE, no view write | Task 2 |
| All structural payloads carry `experiment_id` | Tasks 1, 2, 5, 6, 7, 8, 8b |
| `rebuild_views_from_log!` tolerates `entity_type="sample"` | Task 2 |
| `merged_into_id` nullable column on samples (not hard-delete) | Task 3 |
| `retire_sample!` writer | Task 3 |
| Loads roll-up nested `Load[]`/`LoadSample`/`LoadExposure` shape (§8.8 field names) re-derived from `(load_id, slot_index)`, not log-folded | Task 4 |
| Per-sample `flag` via Phase B Task 12 `derive_sample_flags` + dismiss suppression | Task 4 (seam aligned with Phase B Task 12) |
| `PATCH /api/samples/{id}/name` rename route (idempotent-mutation pattern) | Task 5 |
| `POST /api/exposures/{id}/move` move route + cross-experiment guard | Task 6 |
| `POST /api/samples/{id}/merge` merge route — whole body one `with_idempotency` block | Task 7 |
| Merge: `series_samples` dedup (drop loser membership where survivor present) | Task 7 |
| Merge: `sample_tags` dedup (survivor wins on key collision) | Task 7 |
| Merge: re-point `exposures`, `series_samples`, `sample_tags`, `sample_messages` | Task 7 |
| Merge does NOT emit `sample_created`; no `sample_merged` kind | Task 7, Task 8 |
| `POST /api/samples/{id}/split` split route — mints sample id, emits `sample_created` + per-exposure `exposure_moved` + `sample_split` | Task 8 |
| `sample_created` / `sample_split` / `grouping_flag_dismissed` intent-marker no-ops in dispatcher | Task 8 |
| `grouping_flag_dismissed` durable, undoable event — dismiss + undo routes | Task 8b |
| `GET /api/experiments/{id}/loads` endpoint distinct from flat `/api/samples` | Task 9 |
| Routes registered at server startup (`register_routes!`, wired in Task 5) | Task 5 |
| No regressions against existing suite | Task 10 |

### Contract conformance pass (2026-06-18, second round)

Applied the canonical §8.8/§9.2/§9.3 contract verbatim after a fresh live-source verification. Changes:

- **`experiment_id` added to every structural payload** (`exposure_moved`, `sample_renamed`, `sample_created`, `sample_split`, `grouping_flag_dismissed`), plus `from_sample_id` on `exposure_moved`. **The dispatcher is unchanged** — `update_view_for_event!`'s `exposure_moved` branch reads only `payload.sample_id` (confirmed in `events.jl:305+`); extra payload keys are inert. Route bodies now resolve `experiment_id`/`from_sample_id` from the DB before the event write; the dispatcher tests assert the full payload round-trips.
- **`get_loads_rollup` rewritten to the exact nested `Load[]` shape** (§8.8): load = `(load_id, load_index, session_id, start_time, end_time, frame_count, note, samples)`; sample = `(sample_id, name, slot_index, grouping_source, name_source, merged_into_id, flag, exposures)`; exposure = `(id, filename, horizontal_position, timestamp)`. The exposure leaf is now `id` (not `exposure_id`) and drops `frame_no`/`status` per the contract. The Task 9 serializer and both tests were updated to match.
- **Per-sample `flag`** is produced by calling Phase B Task 12's `derive_sample_flags` over the nested rows the roll-up read, then suppressing any flag whose sample has a non-undone `grouping_flag_dismissed` event (a `NOT EXISTS` join on `user_actions`, with `ua.undoes_event_id IS NULL` so undo rows aren't themselves counted as dismisses). **Phase B Task 12 implements `derive_sample_flags`; arg shape confirmed to match this call site (nested roll-up rows).** The `isdefined` guard is defensive graceful-degradation for incremental build order (yields `flag=null` only if the roll-up is built before Phase B Task 12 lands in the same tree).
- **New Task 8b — `grouping_flag_dismissed`** durable, undoable event: dispatcher no-op (suppression is read-time, Task 4), `POST /api/samples/{id}/dismiss-flag` + `POST /api/samples/{id}/dismiss-flag/undo` (the undo records a dismiss-kind row stamped with `undoes_event_id` pointing at the latest non-undone dismiss — the `peak_excluded`→`peak_unexcluded` precedent at `routes_peaks.jl:296–313`). TDD test asserts dismiss→suppressed, undo→re-shown.
- **Split emits `sample_created`; merge does not; there is no `sample_merged` kind.** The split route emits `sample_created` (payload `{ experiment_id }` only) + per-exposure `exposure_moved` + `sample_split`; the Task 8 dispatcher markers drop the spurious `sample_merged` no-op and add `grouping_flag_dismissed`.

### Review applied (2026-06-18)

This plan was revised against live source after a structural-events review. Fixes folded in:

- **[P0] Route tests rewritten to the live harness.** Every route testset (Tasks 5–9) now drives real in-process HTTP via `with_test_server(db) do port, base … end` with `HTTP.patch/post/get("$base/api/…"; status_exception=false)`, matching `test_routes_samples.jl`. The fabricated `handle_grouping_route(req, db)` shim — which does not exist and was never implemented — was removed everywhere. The in-test `register_grouping_routes!()` calls were dropped (registration happens via `start_test_server!` → `register_routes!`).
- **[P0] Test DB now reaches the handlers.** `with_test_server(db)` calls `bind_db!`/`start_test_server!`, so the handler's `current_db()` resolves to the test `db`. The earlier framing (independent `fresh_db()` handed to a shim) would have read/written a different connection than the handlers.
- **[P0] Wiring corrected and reordered.** The `include("routes_grouping.jl")` goes in `HimalayaUI.jl` (before `include("server.jl")`), not in `server.jl`; the registration call is `register_grouping_routes!()` inside `register_routes!()` (server.jl:49). Both edits moved into Task 5 Step 4 so the route tests in Tasks 5–9 can reach the routes; Task 9 no longer re-wires.
- **[P1] Broadcast-ordering contradiction resolved (opposite of the finding's first suggestion).** Live source shows `_enqueue_broadcast_from_result!` only appends to a task-local queue that `with_idempotency` flushes *after* commit (and clears on rollback/cache-replay). The canonical route `routes_samples.jl:179` calls it INSIDE the `with_idempotency` block. So the routes' existing placement is correct; the **Conventions bullet** (which wrongly said "after the block") was the defect and is now corrected with the cache-hit-no-rebroadcast behavior explained. Moving the call after the block — the finding's stated preference — would enqueue a frame that never fires, so it was deliberately NOT applied. Flagged for human awareness.
- **[P1] Dispatcher unit tests switched to the non-broadcasting variant.** Tasks 1–2 now call `apply_event!(HimalayaUI.InTransaction(), db, req; …)` wrapped in `SQLite.transaction(db)`. The broadcasting default variant calls `broadcast_event!` → `current_db()` (events.jl:1106) on a non-nothing `user_id`, which errors with no server bound.
- **[P2] Confirmed findings 6 & 7 against `events.jl`** (details in the per-task notes and the summary returned to the reviewer).
- **[P2] Task 0 harness cleaned up:** `seed_two_samples` no longer returns `db` (callers already hold it); the unused `sys_req()` helper removed; the `with_test_server` include-guard added (mirrors `test_assignment_reattach.jl`).
- **[P2] Split route IN-clause** annotated as injection-safe (Int-coerced) and given its missing `[id]` bind parameter.

**§9.3 SSE note (deferred to Phase E):** the spec requires that `applyRemoteToCache.ts` on the frontend **refetches** `loads(id)`/`samples` on `sample_created` (split) events (invalidate-only — the new id isn't payload-derivable), reconciles a foreign merge from the `exposure_moved` burst (there is **no** `sample_merged` arm), and does surgical splices for `exposure_moved` + `sample_renamed` (using `payload.experiment_id`/`from_sample_id`, never `remote.entity_id`). `grouping_flag_dismissed` invalidates `loads(id)` so the suppressed/re-shown flag refreshes. The backend SSE broadcast is wired here (via `_enqueue_broadcast_from_result!`, with `experiment_id` in every payload); the frontend cache handlers are Phase E work.

### Placeholder scan

Every code step in this plan provides complete code. No step says "add error handling" or "similar to Task N" without giving the full implementation. The two "read the live file first" notes in Tasks 5 and 9 are accuracy guards for line numbers that drift — the new code is given in full in each case.

### Type and name consistency

- `retire_sample!(db, loser_id; merged_into_id)` — defined Task 3, used Task 7. Signature matches call site.
- `get_loads_rollup(db, experiment_id)` — defined Task 4, route-called Task 9. Return type (nested NamedTuple tree) matches the Task 9 JSON serializer field-for-field (§8.8 names).
- `derive_sample_flags(nested)` — **called** Task 4, **defined in Phase B Task 12** (`grouping.jl`, §9.1); arg shape (nested roll-up rows) confirmed to match this call site. Resolved at runtime via `isdefined` as defensive graceful-degradation for incremental build order.
- Event kind strings (`"exposure_moved"`, `"sample_renamed"`, `"sample_created"`, `"sample_split"`, `"grouping_flag_dismissed"`) are identical between dispatcher branches (Task 1 `exposure_moved`/`sample_renamed`; Task 8 `sample_created`/`sample_split`/`grouping_flag_dismissed` markers), route call sites (Tasks 5–8b), and test assertions. **There is no `"sample_merged"` kind** (merge fans out as `exposure_moved`).
- `register_grouping_routes!()` — defined in `routes_grouping.jl` (Task 5), called inside `register_routes!()` in `server.jl` (wired Task 5 Step 4). Tests do NOT call it directly — they reach the routes over HTTP via `with_test_server`, which triggers `register_routes!` through `start_test_server!`.
- `_enqueue_broadcast_from_result!` — verified at `events.jl:240`. Not exported; the module's `include`-based compilation makes it visible to `routes_grouping.jl` by the bare name (same as `routes_samples.jl:179`). The tests never reference it directly.

### Known accuracy risks

- **`derive_sample_flags` cross-plan seam (resolved).** Spec §9.1/§8.8 pin it as a Phase B export (`grouping.jl`); **Phase B Task 12 (`2026-06-18-ingestion-phase-b-ingest-core.md`) implements it** — `derive_sample_flags(load_rows) → Dict{Int, GroupingFlag}` keyed by `sample_id`, returning the §8.8 `merge`/`split`/`null` shape, consuming the nested roll-up rows. Arg shape confirmed to match Task 4's call site. Task 4's `isdefined(HimalayaUI, :derive_sample_flags)` guard is defensive graceful-degradation for incremental build order (yields `flag=null` only if Phase D's roll-up is built before Phase B Task 12 lands in the same tree), not a sign the function is absent. **One minor sync left:** Phase B Task 12's docstring still names the exposure leaf `(exposure_id, …, frame_no, status)` (pre-rewrite); this plan emits the §8.8 leaf `(id, filename, horizontal_position, timestamp)`. The function reads only `sample_id`/`slot_index`/`horizontal_position`/`timestamp` (identically named in both), so the logic is unaffected — Phase B's docstring should be synced to the §8.8 leaf names for clarity.
- **Phase A columns the roll-up reads** (`loads.session_id`/`note`/`start_time`/`end_time`/`frame_count`, `samples.load_id`/`slot_index`/`grouping_source`/`name_source`/`merged_into_id`, `exposures.experiment_id`/`horizontal_position`/`timestamp`/`frame_no`) do **not** exist in the live schema (verified 2026-06-18: live `samples` is `(id, experiment_id, name, display_name, notes)`; live `exposures` has no `experiment_id`/`horizontal_position`/`timestamp`/`frame_no`; there is no `loads` table). All are Phase A deliverables. Task 4's query and the move/merge routes' `e.experiment_id` reference assume Phase A merged. The standalone test harness seeds them directly (Task 0).
- **`create_sample!` / `create_exposure!` signatures** in the test harness (Task 0) assume the Phase A rewrites (accepting `load_id`, `slot_index`, `experiment_id`). If Phase A is not merged, update the seed calls to the current signatures from `db.jl:1828–1851` (verified: current `create_sample!` takes `experiment_id, name, display_name, notes`; current `create_exposure!` takes `sample_id, filename, kind, selected, status, image_path`).
- **Route test pattern** (Tasks 5–9): the tests drive real in-process HTTP via `with_test_server(db) do port, base … end` (the canonical `test_routes_samples.jl` pattern). There is NO `handle_grouping_route(req, db)` shim — Oxygen routes are macro-registered into a global router with no `(req, db)` dispatch entry point. The harness picks up the routes because `start_test_server!` → `register_routes!` → `register_grouping_routes!` (wired in Task 5).
- **`_enqueue_broadcast_from_result!` visibility:** verified at `events.jl:240`. It is NOT exported, but `routes_grouping.jl` is `include`d into the `HimalayaUI` module so it calls the bare name — same as `routes_samples.jl:179`. The route tests call it only indirectly (through HTTP), so no `HimalayaUI.` prefix is needed in tests.
- **`series` table default row:** the merge test inserts `INSERT INTO series DEFAULT VALUES`. Confirm the `series` table's `NOT NULL` columns (read `db.jl:700–760`) — if any column lacks a default, add the required values to the INSERT.
