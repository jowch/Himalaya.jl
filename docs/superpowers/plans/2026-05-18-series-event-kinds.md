# Series Event-Kind Cluster (#166 / #167 / #168) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire the six series event kinds (`series_created`, `series_recipe_updated`, `series_deleted`, `series_plate_committed`, `series_pinned`, `series_unpinned`) through every layer — the Julia `update_view_for_event!` dispatcher, the frontend mutation queue, and paired contract tests — so the already-built `/api/series*` routes stop producing degenerate placeholder rows.

**Architecture:** Three GitHub issues, landed as one branch in order #166 → #167 → #168 so the three shared switch-statement files (`events.jl`, `applyRemoteToCache.ts`, `mutatorRegistry.ts`) never merge-conflict. The routes (`routes_series.jl`) and business logic (`series.jl`) already exist and already call `apply_event!`; this plan fills the missing dispatcher branches and the missing frontend queue wiring. The view-producing series dispatcher branches are **pure-replace** (parent upsert + child `DELETE`-by-`series_id` + `INSERT`-all from a full-snapshot payload), explicitly *not* the member-id-discriminating shape of `_update_view_for_comparison_submitted!` — this makes every fold idempotent and order-independent so it replays correctly from an empty view.

**Tech Stack:** Julia (SQLite.jl, Oxygen.jl, stdlib `Test`), TypeScript/React (TanStack Query, Vitest).

---

## Decision Record — the `api.ts` / query-key scaffolding inference

**This section documents an inference made during brainstorming. Read it before Task 4.**

The issue map (`docs/superpowers/plans/2026-05-17-himalaya-ui-redesign-issue-map.md`) cards I2.3/I2.4/I2.5 list their frontend files as only `applyRemoteToCache.ts`, `mutatorRegistry.ts`, and new `lib/queue/mutators/series*.ts`. The §3 contention table for `api.ts` / `queries.ts` names only I1.2 / I3.3 / I3.5a — **not** this event-kind cluster.

That file list is incomplete. The `Mutator` interface (`lib/queue/types.ts:159`) has a mandatory `request` field that must call an `api.ts` fetcher, and a verified grep confirms **zero** series frontend code exists today (no `Series` type, no fetcher, no query key, no mutator). The master plan §5.2 mandates the per-kind mutator + `applyRemoteToCache` layers as Phase 2 deliverables, but neither §5.2 nor the issue map legislates the `api.ts` seam those layers structurally depend on.

**Inference / resolution:** the cluster adds the **minimal queue-side scaffolding** its mutators and cache handlers strictly require — a `Series` TS type (plus `SeriesMember`, `SeriesSample`, and the input/body types), the three mutating fetchers (`saveSeries`, `commitSeriesPlate`, `deleteSeries`), and a query-key shape (`queryKeys.series(id)`, `queryKeys.seriesList`, `queryKeys.seriesPins`). It deliberately does **not** add read hooks (`useSeriesList`, `useSeries`) or a `SeriesSummary` listing type — those belong to I3.3 (folio UI), which will build on top of this scaffolding.

**Consequence for I3.3:** this is unlisted contention. When I3.3 is scheduled it must *rebase onto* the cluster's `api.ts` / `queries.ts` additions rather than create them. A corresponding memory has been recorded (`series-event-cluster-frontend-scaffolding`).

---

## Post-review revisions

This plan was reviewed by four agents (himalaya-reviewer, queue-reviewer, frontend-reviewer, a spec-fidelity reviewer). Changes folded in from that pass:

- **`commitSeriesPlate.synthesizeFromSse` returns `post_state` raw** — no `{...base, ...post_state}` merge. The merge would spread `QueueResponseMeta` (`event_id` / `client_op_id` / `analysis_inputs_hash`) into the object `onSuccess` writes to the detail cache via `setQueryData`, polluting the cached `Series` row and diverging from the clean `post_state` splice in `applyRemoteToCache` (Task 8). This is the `GroupMutationResponse` anti-pattern (`api.ts:322`).
- **`series_created` live-path `UPDATE` sets `content_hash = NULL` explicitly** (Task 1) — hardens the draft invariant against a defensive re-fold onto an already-committed row.
- **SSE-broadcast layer is now asserted** — Task 1 adds a `_capture_series_sse` in-process-subscriber helper and asserts a `series_created` frame; Task 7 asserts a `series_plate_committed` frame carries `post_state`. The six-layer rule names SSE broadcast as a distinct layer; it was previously assumed-wired, not tested.
- **`commitSeriesPlate.request` uses the shared `CommitSeriesPlateBody` type** rather than re-declaring it inline (Task 8); test files use `as any` (the established convention) not `as never`; the Task 8 cache test uses a complete typed `Series` fixture (Task 8).
- **Verified-and-kept:** the himalaya-reviewer flagged the Task 1 replay-path `INSERT` as having a column/placeholder count mismatch — this was re-counted and is **correct** (15 columns, 13 `?` + `NULL` + `'draft'` = 15 value slots, 13 binds). No change. `applyPostStateOnly`'s `Array.isArray(post_state.indices)` guard remains safe after `SseEvent.post_state` is widened to include `Series` — a `Series` has no `indices` field, so the guard bails correctly.

---

## Background facts (verified against the worktree)

- **Routes already wired.** `routes_series.jl` POST/PATCH/DELETE/commit/pin/unpin already call `apply_event!(InTransaction(), …)` with the correct `kind` / `entity_type` and (for commit) `post_state = out`. **This plan does not modify `routes_series.jl` or `series.jl`.**
- **Dispatcher gap.** `update_view_for_event!` (`events.jl:305`) has no `series_*` branches; its default is a silent no-op (`events.jl:426-427`). So the routes currently write a `user_actions` row but no view rows, and `fetch_series_with_plate` returns a degenerate placeholder.
- **Payload is canonicalized.** `apply_event!` round-trips the payload through JSON3 before dispatch (`events.jl:159-161`), so every dispatcher branch sees a `JSON3.Object` — `obj.field` and `obj[:field]` both work; a missing key returns `nothing`.
- **`content_hash` rule.** `content_hash` is NULL while `state='draft'` and is computed *only* on `series_plate_committed`, from the plate (`series_members`) — `compute_series_content_hash` (`series.jl:275`) already excludes the recipe. `series_created` and `series_recipe_updated` therefore never touch `content_hash`.
- **`series` table columns** (`db.jl:728`): `id, title, description, content_hash, created_by, created_at, updated_at, forked_from_id, forked_at_hash, view_grouping_mode, view_show_peak_ticks, view_show_peak_labels, ordering_variable, order_rule, state`. `order_rule` defaults `'manual'`, `state` defaults `'committed'`.
- **`series_samples` columns** (`db.jl:782`): `id, series_id, sample_id, position, pinned, excluded` — `UNIQUE(series_id, position)`, `pinned`/`excluded` are `CHECK (… IN (0,1))`.
- **`series_pins`** (`db.jl:808`): `user_id, series_id, pinned_at`, `PRIMARY KEY (user_id, series_id)`.
- **`series_samples.id` is replay-volatile** (master plan §11) — `DELETE`+`INSERT` re-mints ids on every fold, so client state must key recipe rows on `(series_id, position)`, never on `series_samples.id`.

## Running tests

Backend (focused dev loop — the full suite is 5–10 min, see `test/AGENTS.md`):

```bash
julia --project=packages/HimalayaUI -e '
using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_routes_series.jl")
' 2>&1 | grep -E "Test Summary|Pass|Fail|Error|did not pass"
```

(`with_test_server` is defined in `test_http.jl`; including it first puts that helper in scope. Including `test_http.jl` also runs its own testsets — expected and acceptable for the dev loop.)

Frontend (from `packages/HimalayaUI/frontend/`):

```bash
npx vitest run test/queue/<file>.test.ts     # one file
npm run build                                # tsc --noEmit + vite — final gate
```

**Final acceptance gate (run once at the end, Task 10):** the full Julia suite (`julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'`) and `npm run build` + `npx vitest run` must be green.

---

## File Structure

**Modified — Julia:**
- `packages/HimalayaUI/src/events.jl` — six new dispatcher branches in `update_view_for_event!`; four new helper functions (`_update_view_for_series_created!`, `_update_view_for_series_recipe_updated!`, `_update_view_for_series_plate_committed!`, `_insert_series_sample!`, `_insert_series_member!`).
- `packages/HimalayaUI/test/test_routes_series.jl` — new behavioural testsets per kind (route round-trip + `rebuild_views_from_log!` fold-from-empty).

**Modified — frontend:**
- `frontend/src/api.ts` — `Series` / `SeriesMember` / `SeriesSample` types, input/body types, `saveSeries` / `commitSeriesPlate` / `deleteSeries` fetchers.
- `frontend/src/queries.ts` — `queryKeys.series` / `seriesList` / `seriesPins`.
- `frontend/src/lib/queue/types.ts` — `OpKind` gains `series_save` / `series_commit` / `series_delete`; `SseEvent.post_state` widened to include `Series`.
- `frontend/src/lib/queue/mutatorRegistry.ts` — `resolveMutator` + `resolveMutatorForEvent` cases.
- `frontend/src/lib/queue/applyRemoteToCache.ts` — six new `case` arms.

**Created — frontend:**
- `frontend/src/lib/queue/mutators/saveSeries.ts` — covers `series_created` + `series_recipe_updated` (one mutator, `payload.id` discriminates).
- `frontend/src/lib/queue/mutators/commitSeriesPlate.ts` — covers `series_plate_committed`.
- `frontend/src/lib/queue/mutators/deleteSeries.ts` — covers `series_deleted`.
- `frontend/test/queue/applyRemoteToCache.series.test.ts` — cache-handler contract tests.
- `frontend/test/queue/seriesMutators.test.ts` — mutator `onSuccess` / `synthesizeFromSse` tests.

Series dispatcher helpers live in `events.jl` next to `_update_view_for_comparison_created!` / `_insert_comparison_member!` — following the existing pattern, not `series.jl` (which holds route-facing business logic).

---

## Task 1: `series_created` dispatcher branch (#166)

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl` — add branch in `update_view_for_event!` (after line 422); add `_update_view_for_series_created!` and `_insert_series_sample!` (after `_insert_comparison_member!`, ~line 681).
- Test: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Add the SSE-capture helper**

The "SSE `broadcast_event!`" layer is one of the six layers and must be asserted, not assumed. Add this helper to `test_routes_series.jl` at top level, just after the `_series_test_db` function (before `@testset "Series routes"`). It is the in-process subscriber pattern from `test/AGENTS.md` (a self-contained local copy of `test_idempotency_replay_invariant.jl::_capture_sse_during`, so `test_routes_series.jl` has no cross-file include-order dependency):

```julia
# Capture the SSE frames of one kind broadcast while `f` runs — the
# in-process subscriber pattern (test/AGENTS.md). `do`-block form: `f` is the
# FIRST argument, so `_capture_series_sse("k") do … end` ⇒ (f, "k").
function _capture_series_sse(f::Function, kind::String)
    pending = Channel{String}(64)
    sub = (pending = pending,)
    lock(HimalayaUI.SSE_LOCK) do
        push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
    end
    try
        f()
        sleep(0.3)   # let the post-commit broadcast queue flush
    finally
        lock(HimalayaUI.SSE_LOCK) do
            filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
        end
        close(pending)
    end
    [fr for fr in pending
        if !startswith(fr, ":") && occursin("\"kind\":\"$kind\"", fr)]
end
```

- [ ] **Step 2: Write the failing test**

Add this testset inside the `@testset "Series routes" begin … end` block in `test_routes_series.jl` (e.g. after the `GET /api/series/{id}` testset):

```julia
    @testset "POST /api/series — series_created writes the recipe + draft state" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                body = JSON3.write(Dict(
                    :title   => "Heat ramp",
                    :samples => [Dict(:sample_id => 100, :position => 0, :pinned => true)],
                ))
                resp = HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"], body)
                @test resp.status == 201
                got = JSON3.read(resp.body, Dict{Symbol, Any})
                sid = got[:id]
                # The dispatcher upserts state='draft' and leaves content_hash NULL.
                @test got[:state] == "draft"
                @test got[:content_hash] == ""           # fetch_series_with_plate maps NULL → ""
                @test got[:title] == "Heat ramp"
                # The recipe snapshot landed in series_samples.
                @test length(got[:samples]) == 1
                @test got[:samples][1]["sample_id"] == 100
                @test got[:samples][1]["pinned"] == true
                @test isempty(got[:members])             # series_created carries zero members

                # rebuild_views_from_log! round-trip: empty the view rows, re-fold.
                DBInterface.execute(db, "DELETE FROM series_samples WHERE series_id = ?", [sid])
                DBInterface.execute(db, "DELETE FROM series WHERE id = ?", [sid])
                HimalayaUI.rebuild_views_from_log!(db, sid; entity_type = "series")
                refold = HimalayaUI.fetch_series_with_plate(db, sid)
                @test refold !== nothing
                @test refold[:state] == "draft"
                @test length(refold[:samples]) == 1
                @test refold[:samples][1][:sample_id] == 100

                # SSE layer: a second create, observed through the in-process
                # subscriber, must broadcast exactly one series_created frame
                # carrying entity_type='series'.
                frames = _capture_series_sse("series_created") do
                    HTTP.post("$base/api/series",
                        ["X-Username" => "alice", "Content-Type" => "application/json"],
                        JSON3.write(Dict(:title => "Frame check",
                            :samples => [Dict(:sample_id => 100, :position => 0)])))
                end
                @test length(frames) == 1
                @test occursin("\"entity_type\":\"series\"", frames[1])
            end
            close(db)
        end
    end
```

- [ ] **Step 3: Run the test to verify it fails**

Run the focused backend command (see "Running tests"). Expected: FAIL — `got[:state]` is `"committed"` (schema default on the degenerate placeholder) and `got[:samples]` is empty, because no `series_created` dispatcher branch exists.

- [ ] **Step 4: Add the dispatcher branch**

In `events.jl`, inside `update_view_for_event!`, immediately after the `comparison_unpinned` branch (after line 422) and before the `noop_test` line, insert:

```julia
    # Series event kinds (#166 / I2.3). The view-producing branches are
    # pure-replace: parent upsert + child DELETE-by-series_id + INSERT-all from
    # a full-snapshot payload — explicitly NOT comparison_submitted-shaped, so
    # every fold is idempotent and order-independent (folds from an empty view).
    if kind == "series_created"
        return _update_view_for_series_created!(db, entity_id, payload, event_id)
    end
    if kind == "series_recipe_updated"
        return _update_view_for_series_recipe_updated!(db, entity_id, payload, event_id)
    end
    if kind == "series_deleted"
        DBInterface.execute(db, "DELETE FROM series WHERE id = ?", [Int(entity_id)])
        # The schema's four ON DELETE CASCADE clauses drop series_members,
        # series_samples, series_messages and series_pins automatically.
        return nothing
    end
```

- [ ] **Step 5: Add the helper functions**

In `events.jl`, after `_insert_comparison_member!` (ends ~line 681), add:

```julia
# Helper: INSERT one series_samples row from a recipe payload entry. The
# payload entry is a JSON3.Object from `_series_sample_payload` —
# `{sample_id, position, pinned, excluded}`. `pinned`/`excluded` are coerced
# to 0/1 for the CHECK (… IN (0,1)) columns.
function _insert_series_sample!(db, series_id, s)
    pinned   = (haskey(s, :pinned)   && s.pinned   == true) ? 1 : 0
    excluded = (haskey(s, :excluded) && s.excluded == true) ? 1 : 0
    DBInterface.execute(db,
        """INSERT INTO series_samples
             (series_id, sample_id, position, pinned, excluded)
           VALUES (?, ?, ?, ?, ?)""",
        [Int(series_id), Int(s.sample_id), Int(s.position), pinned, excluded])
    nothing
end

"""
    _update_view_for_series_created!(db, entity_id, payload, event_id)

`series_created` dispatcher (#166). Upserts the `series` row at `entity_id`
(the route mints it with `INSERT INTO series DEFAULT VALUES` to capture the
AUTOINCREMENT id; a plain INSERT would collide on the live path — so this
SELECTs and UPDATEs an existing row, else INSERTs with an explicit id for the
replay-from-empty path). Sets `state='draft'`; `content_hash` stays NULL (a
draft has no committed plate — master plan §5.1). Then pure-replaces
`series_samples` from the full payload snapshot. Touches no `series_members`.
"""
function _update_view_for_series_created!(db, entity_id, payload, event_id)
    sid     = Int(entity_id)
    user_id = user_id_for_event(db, event_id)
    now_str = comparison_now_iso()

    title          = haskey(payload, :title) && payload.title !== nothing ?
                     String(payload.title) : nothing
    description    = haskey(payload, :description) && payload.description !== nothing ?
                     String(payload.description) : nothing
    forked_from_id = haskey(payload, :forked_from_id) && payload.forked_from_id !== nothing ?
                     Int(payload.forked_from_id) : nothing
    forked_at_hash = haskey(payload, :forked_at_hash) && payload.forked_at_hash !== nothing ?
                     String(payload.forked_at_hash) : nothing
    ordering_var   = haskey(payload, :ordering_variable) && payload.ordering_variable !== nothing ?
                     String(payload.ordering_variable) : nothing
    order_rule     = haskey(payload, :order_rule) && payload.order_rule !== nothing ?
                     String(payload.order_rule) : "manual"
    vgm  = haskey(payload, :view_grouping_mode) && payload.view_grouping_mode !== nothing ?
           String(payload.view_grouping_mode) : nothing
    vspt = haskey(payload, :view_show_peak_ticks)  ? payload.view_show_peak_ticks  : nothing
    vspl = haskey(payload, :view_show_peak_labels) ? payload.view_show_peak_labels : nothing

    existing = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM series WHERE id = ?", [sid]))
    if isempty(existing)
        # Replay path: the original row was deleted. INSERT with an explicit id;
        # SQLite's AUTOINCREMENT counter advances past it automatically.
        DBInterface.execute(db,
            """INSERT INTO series
               (id, title, description, content_hash, created_by, created_at,
                updated_at, forked_from_id, forked_at_hash,
                view_grouping_mode, view_show_peak_ticks, view_show_peak_labels,
                ordering_variable, order_rule, state)
               VALUES (?, ?, ?, NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'draft')""",
            [sid, title, description, user_id, now_str, now_str,
             forked_from_id, forked_at_hash, vgm, vspt, vspl,
             ordering_var, order_rule])
    else
        # Live path: the route's `INSERT … DEFAULT VALUES` placeholder exists.
        # COALESCE stamps created_at on first fold; a prior fold's value
        # survives. `content_hash = NULL` is set explicitly (not just relied on
        # from the placeholder default) so the draft invariant — a draft has a
        # NULL content_hash (master plan §5.1) — holds even on a defensive
        # re-fold onto a row left committed by a prior `series_plate_committed`.
        DBInterface.execute(db,
            """UPDATE series
               SET title = ?, description = ?, content_hash = NULL,
                   created_by = ?,
                   created_at = COALESCE(created_at, ?), updated_at = ?,
                   forked_from_id = ?, forked_at_hash = ?,
                   view_grouping_mode = ?, view_show_peak_ticks = ?,
                   view_show_peak_labels = ?,
                   ordering_variable = ?, order_rule = ?, state = 'draft'
               WHERE id = ?""",
            [title, description, user_id, now_str, now_str,
             forked_from_id, forked_at_hash, vgm, vspt, vspl,
             ordering_var, order_rule, sid])
    end

    # Pure-replace the recipe rows from the full-snapshot payload.
    DBInterface.execute(db, "DELETE FROM series_samples WHERE series_id = ?", [sid])
    samples = haskey(payload, :samples) ? payload.samples : []
    for s in samples
        _insert_series_sample!(db, sid, s)
    end
    return sid
end
```

(`_update_view_for_series_recipe_updated!` is added in Task 2 — Task 1's branch references it but Task 1's test only exercises `series_created`; the function must nonetheless be defined for the file to load. Add a minimal stub now and replace it in Task 2, OR implement Tasks 1 and 2 before the first run. **Recommended:** add the stub below now.)

```julia
# Stub — full implementation lands in Task 2 (#166, series_recipe_updated).
function _update_view_for_series_recipe_updated!(db, entity_id, payload, event_id)
    error("series_recipe_updated dispatcher not yet implemented")
end
```

- [ ] **Step 6: Run the test to verify it passes**

Run the focused backend command. Expected: PASS — the new `POST /api/series` testset (route round-trip, fold-from-empty, and SSE-frame assertions) is green.

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/src/events.jl packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(cat <<'EOF'
Add the series_created dispatcher branch + SSE-capture helper (#166)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: `series_recipe_updated` dispatcher branch (#166)

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl` — replace the `_update_view_for_series_recipe_updated!` stub.
- Test: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Add this testset inside `@testset "Series routes"`:

```julia
    @testset "PATCH /api/series/{id} — series_recipe_updated pure-replaces the recipe" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            DBInterface.execute(db,
                "INSERT INTO samples (id, experiment_id, name) VALUES (101, 10, 'sB')")
            with_test_server(db) do port, base
                # Create a draft with one recipe row.
                createBody = JSON3.write(Dict(
                    :title   => "Ramp",
                    :samples => [Dict(:sample_id => 100, :position => 0)],
                ))
                created = JSON3.read(HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    createBody).body, Dict{Symbol, Any})
                sid = created[:id]

                # PATCH with a completely different recipe — pure-replace, not a diff.
                patchBody = JSON3.write(Dict(
                    :ordering_variable => "temperature",
                    :order_rule        => "ascending",
                    :samples           => [
                        Dict(:sample_id => 101, :position => 0),
                        Dict(:sample_id => 100, :position => 1),
                    ],
                ))
                resp = HTTP.patch("$base/api/series/$sid",
                    ["X-Username" => "alice", "Content-Type" => "application/json"], patchBody)
                @test resp.status == 200
                got = JSON3.read(resp.body, Dict{Symbol, Any})
                @test got[:ordering_variable] == "temperature"
                @test got[:order_rule] == "ascending"
                @test length(got[:samples]) == 2
                @test got[:samples][1]["sample_id"] == 101   # ordered by position
                @test got[:samples][2]["sample_id"] == 100
                @test got[:state] == "draft"                  # recipe edit never commits
                @test got[:content_hash] == ""                # recipe edit never hashes

                # rebuild_views_from_log! round-trip: fold series_created THEN
                # series_recipe_updated from empty; final state == post-PATCH.
                DBInterface.execute(db, "DELETE FROM series_samples WHERE series_id = ?", [sid])
                DBInterface.execute(db, "DELETE FROM series WHERE id = ?", [sid])
                HimalayaUI.rebuild_views_from_log!(db, sid; entity_type = "series")
                refold = HimalayaUI.fetch_series_with_plate(db, sid)
                @test length(refold[:samples]) == 2
                @test refold[:samples][1][:sample_id] == 101
                @test refold[:ordering_variable] == "temperature"
            end
            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the focused backend command. Expected: FAIL — the route hits the `_update_view_for_series_recipe_updated!` stub, which `error()`s; the PATCH returns 500.

- [ ] **Step 3: Replace the stub with the real implementation**

In `events.jl`, replace the Task-1 stub:

```julia
"""
    _update_view_for_series_recipe_updated!(db, entity_id, payload, event_id)

`series_recipe_updated` dispatcher (#166). Pure-replace: updates the recipe
scalars (`ordering_variable`, `order_rule`) on the `series` row, then
`DELETE`s every `series_samples` row for the series and re-`INSERT`s the full
payload snapshot. Never touches `content_hash` (recipe is excluded from the
hash — master plan §5.1), `state`, or `series_members`. `ordering_variable`
is a bare write (a PATCH is a full recipe replace, so an omitted field nulls
it); `order_rule` is `COALESCE`d because the column is `NOT NULL` and
`CHECK`-constrained.
"""
function _update_view_for_series_recipe_updated!(db, entity_id, payload, event_id)
    sid     = Int(entity_id)
    now_str = comparison_now_iso()
    ordering_var = haskey(payload, :ordering_variable) && payload.ordering_variable !== nothing ?
                   String(payload.ordering_variable) : nothing
    order_rule   = haskey(payload, :order_rule) && payload.order_rule !== nothing ?
                   String(payload.order_rule) : nothing

    DBInterface.execute(db,
        """UPDATE series
           SET ordering_variable = ?,
               order_rule = COALESCE(?, order_rule),
               updated_at = ?
           WHERE id = ?""",
        [ordering_var, order_rule, now_str, sid])

    DBInterface.execute(db, "DELETE FROM series_samples WHERE series_id = ?", [sid])
    samples = haskey(payload, :samples) ? payload.samples : []
    for s in samples
        _insert_series_sample!(db, sid, s)
    end
    return nothing
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run the focused backend command. Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/events.jl packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(cat <<'EOF'
Add the series_recipe_updated dispatcher branch (#166)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: `series_deleted` dispatcher branch (#166)

The `series_deleted` branch itself was already added in Task 1 (the inline `DELETE FROM series`). This task adds its behavioural test.

**Files:**
- Test: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Add this testset inside `@testset "Series routes"`:

```julia
    @testset "DELETE /api/series/{id} — series_deleted cascades + folds to absent" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                createBody = JSON3.write(Dict(
                    :title   => "Doomed",
                    :samples => [Dict(:sample_id => 100, :position => 0)],
                ))
                sid = JSON3.read(HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    createBody).body, Dict{Symbol, Any})[:id]

                resp = HTTP.delete("$base/api/series/$sid", ["X-Username" => "alice"])
                @test resp.status == 200
                # Four-table cascade: the series row and its recipe rows are gone.
                @test HimalayaUI.fetch_series_with_plate(db, sid) === nothing
                remaining = Tables.rowtable(DBInterface.execute(db,
                    "SELECT COUNT(*) AS n FROM series_samples WHERE series_id = ?", [sid]))
                @test remaining[1].n == 0

                # rebuild_views_from_log! round-trip: fold series_created THEN
                # series_deleted from empty → series stays absent.
                HimalayaUI.rebuild_views_from_log!(db, sid; entity_type = "series")
                @test HimalayaUI.fetch_series_with_plate(db, sid) === nothing
            end
            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it passes**

Run the focused backend command. Expected: PASS — the `series_deleted` branch already exists from Task 1. (If it fails, the Task-1 branch was not added correctly.)

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(cat <<'EOF'
Add the series_deleted dispatcher round-trip test (#166)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Frontend series scaffolding — types, fetchers, query keys (#166)

See the **Decision Record** at the top of this plan. This task adds the minimal `api.ts` / `queries.ts` scaffolding the queue layers structurally require.

**Files:**
- Modify: `frontend/src/api.ts` (append after the comparison-pins block, ~line 666)
- Modify: `frontend/src/queries.ts` (`queryKeys` object, after `comparisonMessages`, ~line 88)
- Modify: `frontend/src/lib/queue/types.ts` (`OpKind` union; `SseEvent.post_state`)

- [ ] **Step 1: Add the series types + fetchers to `api.ts`**

Append to `frontend/src/api.ts` (the `MemberSnapshot` type at line 400 is reused):

```typescript
// ─── Series (#166 / #167 / #168 — series event-kind cluster) ────────────────
//
// Minimal queue-side scaffolding only — the read hooks (useSeriesList /
// useSeries) and the listing summary type belong to I3.3 (folio UI). See the
// Decision Record in docs/superpowers/plans/2026-05-18-series-event-kinds.md.
// Shapes mirror `fetch_series_with_plate` in series.jl.

/** The recipe membership — one `series_samples` row. */
export interface SeriesSample {
  id: number;
  series_id: number;
  sample_id: number;
  position: number;
  pinned: boolean;
  excluded: boolean;
}

/** The plate — one `series_members` row. Mirrors `ComparisonMember`. */
export interface SeriesMember {
  id: number;
  series_id: number;
  exposure_id: number | null;
  display_order: number;
  band_height: number;
  y_offset: number;
  normalization: string;
  color_override: string | null;
  label_override: string | null;
  q_window_min: number | null;
  q_window_max: number | null;
  peak_display: unknown;
  snapshot: MemberSnapshot | null;
  is_stale: boolean;
  created_by: number | null;
  created_at: string | null;
}

/** Full nested response from `GET/POST/PATCH /api/series*`. */
export interface Series {
  id: number;
  title: string;
  description: string | null;
  content_hash: string;
  created_by: number | null;
  created_at: string | null;
  updated_at: string | null;
  forked_from_id: number | null;
  forked_at_hash: string | null;
  forked_from_title: string | null;
  view_grouping_mode: string | null;
  view_show_peak_ticks: boolean | null;
  view_show_peak_labels: boolean | null;
  ordering_variable: string | null;
  order_rule: string;
  state: string;
  members: SeriesMember[];
  samples: SeriesSample[];
}

/** Per-recipe-row input for `POST /api/series` and `PATCH /api/series/:id`. */
export interface SeriesSampleInput {
  sample_id: number;
  position?: number;
  pinned?: boolean;
  excluded?: boolean;
}

/** Per-plate-member input for `POST /api/series/:id/commit`. Members carry no
 *  id — the dispatcher mints them. `snapshot` is server-filled when omitted. */
export interface SeriesMemberInput {
  exposure_id: number | null;
  display_order: number;
  band_height?: number;
  y_offset?: number;
  normalization?: string;
  color_override?: string | null;
  label_override?: string | null;
  q_window_min?: number | null;
  q_window_max?: number | null;
  peak_display?: unknown;
  snapshot?: MemberSnapshot;
}

/** Body for `POST /api/series` (create) and `PATCH /api/series/:id` (recipe). */
export interface SaveSeriesBody {
  title: string;
  description?: string | null;
  samples: SeriesSampleInput[];
  ordering_variable?: string | null;
  order_rule?: "ascending" | "descending" | "manual";
  forked_from_id?: number | null;
  forked_at_hash?: string | null;
  view_grouping_mode?: string | null;
  view_show_peak_ticks?: boolean | null;
  view_show_peak_labels?: boolean | null;
}

/** Body for `POST /api/series/:id/commit`. */
export interface CommitSeriesPlateBody {
  members: SeriesMemberInput[];
  expected_content_hash?: string;
}

/**
 * Save a series — create (no id ⇒ `POST /api/series`) or recipe-edit
 * (id present ⇒ `PATCH /api/series/:id`). Branches on id so the queue
 * mutator's `request` stays a single-payload call (mirrors `saveComparison`).
 */
export async function saveSeries(
  body: SaveSeriesBody,
  seriesId: number | undefined,
  opts?: AuthOpts,
): Promise<Series> {
  return seriesId === undefined
    ? request<Series>("POST", "/api/series", body, opts)
    : request<Series>("PATCH", `/api/series/${seriesId}`, body, opts);
}

/** Commit the plate (the old "submit"). A 409 surfaces as a generic `ApiError`;
 *  the typed conflict modal is I3.5b's concern, not the cluster's. */
export const commitSeriesPlate = (
  id: number, body: CommitSeriesPlateBody, opts?: AuthOpts,
) => request<Series>("POST", `/api/series/${id}/commit`, body, opts);

export const deleteSeries = (id: number, opts?: AuthOpts) =>
  request<{ id: number; deleted: boolean; event_id: number }>(
    "DELETE", `/api/series/${id}`, undefined, opts);
```

- [ ] **Step 2: Add the series query keys to `queries.ts`**

In `frontend/src/queries.ts`, inside the `queryKeys` object, after the `comparisonMessages` entry (~line 88), add:

```typescript
  // Series (#166/#167/#168). Detail root `["series", id]`, listing root
  // `["series-list"]` — distinct roots so a listing invalidation never
  // clobbers a detail entry (mirrors the comparison/comparisons split). Read
  // hooks (useSeriesList / useSeries) are added by I3.3.
  series:     (id: number | undefined) => ["series", id ?? "none"] as const,
  seriesList: ["series-list"] as const,
  seriesPins: ["series-pins"] as const,
```

- [ ] **Step 3: Extend `OpKind` and `SseEvent.post_state` in `types.ts`**

In `frontend/src/lib/queue/types.ts`, change the import on line 2 and the `OpKind` union (line 52) and `SseEvent.post_state` (line 119):

Import line — add `Series`:

```typescript
import type { Comparison, Series } from "../../api";
```

`OpKind` — replace the final line (`  | "comparison_save" | "comparison_delete";`) with:

```typescript
  | "comparison_save" | "comparison_delete"
  // Series (#166/#167/#168). `series_save` covers create (POST /api/series)
  // AND recipe edit (PATCH /api/series/:id) — `payload.id` discriminates,
  // mirroring `comparison_save`. `series_commit` → `series_plate_committed`.
  // `series_delete` → `series_deleted`.
  | "series_save" | "series_commit" | "series_delete";
```

`SseEvent.post_state` — replace line 119:

```typescript
  post_state?: CurationPostState | Comparison | Series;
```

- [ ] **Step 4: Verify the types compile**

Run from `packages/HimalayaUI/frontend/`:

```bash
npx tsc --noEmit
```

Expected: PASS (no errors). This task adds only declarations — no behaviour yet.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/api.ts packages/HimalayaUI/frontend/src/queries.ts packages/HimalayaUI/frontend/src/lib/queue/types.ts
git commit -m "$(cat <<'EOF'
Add minimal frontend series scaffolding for the queue layer (#166)

api.ts Series types + fetchers, queryKeys, OpKind union. See the Decision
Record in docs/superpowers/plans/2026-05-18-series-event-kinds.md — the
issue map omits this api.ts seam; a Mutator structurally requires it.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: `saveSeries` + `deleteSeries` mutators + registry wiring (#166)

**Files:**
- Create: `frontend/src/lib/queue/mutators/saveSeries.ts`
- Create: `frontend/src/lib/queue/mutators/deleteSeries.ts`
- Modify: `frontend/src/lib/queue/mutatorRegistry.ts`
- Test: `frontend/test/queue/seriesMutators.test.ts` (create)

- [ ] **Step 1: Write the failing test**

Create `frontend/test/queue/seriesMutators.test.ts`:

```typescript
import { describe, it, expect } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { saveSeriesMutator } from "../../src/lib/queue/mutators/saveSeries";
import { deleteSeriesMutator } from "../../src/lib/queue/mutators/deleteSeries";
import { queryKeys } from "../../src/queries";
import type { Series } from "../../src/api";

function fullSeries(id: number): Series {
  return {
    id, title: "S", description: null, content_hash: "sha256:abc",
    created_by: 1, created_at: null, updated_at: null,
    forked_from_id: null, forked_at_hash: null, forked_from_title: null,
    view_grouping_mode: null, view_show_peak_ticks: null,
    view_show_peak_labels: null, ordering_variable: null, order_rule: "manual",
    state: "draft", members: [], samples: [],
  };
}

describe("saveSeriesMutator", () => {
  it("onSuccess writes a full Series response into the detail cache", () => {
    const qc = new QueryClient();
    const resp = fullSeries(7);
    // `as any` for the flat-payload arg — the established test convention
    // (see deleteComparison.test.tsx); a precise FlatPayload literal is noise.
    saveSeriesMutator.onSuccess(
      { kind: "series_save", payload: {}, clientOpId: "op1",
        title: "S", samples: [], username: "alice", clientId: "c1" } as any,
      resp, qc);
    expect(qc.getQueryData<Series>(queryKeys.series(7))).toEqual(resp);
  });

  it("synthesizeFromSse yields a partial shape carrying the entity id", () => {
    const synth = saveSeriesMutator.synthesizeFromSse?.(
      { id: 99, kind: "series_created", entity_type: "series", entity_id: 7,
        payload: { title: "S", samples: [] } },
      { event_id: 99, client_op_id: "op1", analysis_inputs_hash: undefined });
    expect(synth?.id).toBe(7);
  });
});

describe("deleteSeriesMutator", () => {
  it("onSuccess removes the detail cache and filters the listing", () => {
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.series(7), fullSeries(7));
    qc.setQueryData(queryKeys.seriesList, [{ id: 7 }, { id: 8 }]);
    deleteSeriesMutator.onSuccess(
      { kind: "series_delete", payload: {}, clientOpId: "op1",
        id: 7, username: "alice", clientId: "c1" } as any,
      { id: 7, deleted: true, event_id: 99 }, qc);
    expect(qc.getQueryData(queryKeys.series(7))).toBeUndefined();
    expect(qc.getQueryData(queryKeys.seriesList)).toEqual([{ id: 8 }]);
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run from `packages/HimalayaUI/frontend/`:

```bash
npx vitest run test/queue/seriesMutators.test.ts
```

Expected: FAIL — `Cannot find module './mutators/saveSeries'`.

- [ ] **Step 3: Create `saveSeries.ts`**

Create `frontend/src/lib/queue/mutators/saveSeries.ts`:

```typescript
/**
 * saveSeries mutator (#166). One mutator handles BOTH create and recipe edit;
 * `payload.id` discriminates — absent ⇒ create (POST /api/series), present ⇒
 * recipe edit (PATCH /api/series/:id). Mirrors `saveComparison`.
 *
 * No optimistic write — the builder UI (I3.5b) shows the local draft; this
 * mutator's job is the server round-trip + cache reconciliation. `onSuccess`
 * splices a full Series response into the detail cache and invalidates the
 * listing; on the SSE-wins path the response is a partial shape, detected and
 * routed to invalidate-fallback.
 */
import * as api from "../../../api";
import type { AuthOpts, Series, SaveSeriesBody, SeriesSampleInput } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";

export interface SaveSeriesInput {
  /** Absent ⇒ create; present ⇒ recipe edit of an existing series. */
  id?: number;
  title: string;
  description?: string | null;
  samples: SeriesSampleInput[];
  ordering_variable?: string | null;
  order_rule?: "ascending" | "descending" | "manual";
  forked_from_id?: number | null;
  forked_at_hash?: string | null;
  view_grouping_mode?: string | null;
  view_show_peak_ticks?: boolean | null;
  view_show_peak_labels?: boolean | null;
}

interface SaveSeriesScope {
  username: string | undefined;
  clientId: string;
}

function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

function buildBody(p: SaveSeriesInput): SaveSeriesBody {
  const out: SaveSeriesBody = { title: p.title, samples: p.samples };
  if (p.description !== undefined) out.description = p.description;
  if (p.ordering_variable !== undefined) out.ordering_variable = p.ordering_variable;
  if (p.order_rule !== undefined) out.order_rule = p.order_rule;
  if (p.forked_from_id !== undefined) out.forked_from_id = p.forked_from_id;
  if (p.forked_at_hash !== undefined) out.forked_at_hash = p.forked_at_hash;
  if (p.view_grouping_mode !== undefined) out.view_grouping_mode = p.view_grouping_mode;
  if (p.view_show_peak_ticks !== undefined) out.view_show_peak_ticks = p.view_show_peak_ticks;
  if (p.view_show_peak_labels !== undefined) out.view_show_peak_labels = p.view_show_peak_labels;
  return out;
}

export const saveSeriesMutator: Mutator<SaveSeriesInput, SaveSeriesScope, Series> = {
  kind: "series_save",
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: (p) => api.saveSeries(buildBody(p), p.id, buildAuthOpts(p)),
  onSuccess: (_p, response, qc) => {
    // SSE-wins path: `synthesizeFromSse` yields a partial shape (no `samples`
    // array, no `state`). Probe for the full Series shape; fall back to
    // invalidate so the next read fetches the canonical projection.
    const looksFull = Array.isArray((response as { samples?: unknown }).samples)
      && typeof (response as { state?: unknown }).state === "string";
    if (looksFull) {
      qc.setQueryData(queryKeys.series(response.id), response);
    } else if (typeof response?.id === "number") {
      qc.invalidateQueries({ queryKey: queryKeys.series(response.id) });
    }
    qc.invalidateQueries({ queryKey: queryKeys.seriesList });
  },
  synthesizeFromSse: (remote, base) => {
    const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
    // Partial Series shape — `onSuccess`'s looksFull detector trips the
    // invalidate fallback (no `state` in the SSE payload). `id` is required so
    // the invalidate branch targets the right cache key.
    return { ...base, ...payload, id: remote.entity_id } as unknown as Series;
  },
};
```

- [ ] **Step 4: Create `deleteSeries.ts`**

Create `frontend/src/lib/queue/mutators/deleteSeries.ts`:

```typescript
/**
 * deleteSeries mutator (#166). Mirrors `deleteComparison`: no optimistic
 * write; `onSuccess` removes the detail cache (NOT invalidate — refetching a
 * deleted resource 404s) and filters the id out of the listing cache.
 */
import * as api from "../../../api";
import type { AuthOpts } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";

export interface DeleteSeriesInput {
  id: number;
}

interface DeleteSeriesScope {
  username: string | undefined;
  clientId: string;
}

function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

type DeleteResponse = { id: number; deleted: boolean; event_id: number };

export const deleteSeriesMutator: Mutator<DeleteSeriesInput, DeleteSeriesScope, DeleteResponse> = {
  kind: "series_delete",
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: (p) => api.deleteSeries(p.id, buildAuthOpts(p)),
  onSuccess: (p, _response, qc) => {
    qc.removeQueries({ queryKey: queryKeys.series(p.id) });
    // Filter the id out of the cached listing. Typed structurally (`{id}`) —
    // the cluster does not own the listing summary type (I3.3 does).
    qc.setQueriesData<{ id: number }[]>(
      { queryKey: queryKeys.seriesList },
      (old) => (old ? old.filter((s) => s.id !== p.id) : old),
    );
  },
  // 404 = "already gone" → desired end state (idempotent delete under retry).
  treats404AsSuccess: true,
  synthesizeFromSse: (remote, base) => ({
    ...base,
    id: remote.entity_id,
    deleted: true,
  } as DeleteResponse),
};
```

- [ ] **Step 5: Wire the registry**

In `frontend/src/lib/queue/mutatorRegistry.ts`:

Add imports after line 29 (`import { deleteComparisonMutator } …`):

```typescript
import { saveSeriesMutator } from "./mutators/saveSeries";
import { deleteSeriesMutator } from "./mutators/deleteSeries";
```

In `resolveMutator`, after the `case "comparison_delete":` arm (line 110), add:

```typescript
    case "series_save":
      return saveSeriesMutator;
    case "series_delete":
      return deleteSeriesMutator;
```

In `resolveMutatorForEvent`, after the `case "comparison_deleted":` arm (line 145), add:

```typescript
    case "series_created":
    case "series_recipe_updated":
      return saveSeriesMutator;
    case "series_deleted":
      return deleteSeriesMutator;
```

(`series_pinned` / `series_unpinned` are deliberately *not* cased — they have no mutator and fall through to the `default: return undefined` arm, exactly like `comparison_pinned`. `series_commit` / `series_plate_committed` are added in Task 8.)

- [ ] **Step 6: Run the test to verify it passes**

```bash
npx vitest run test/queue/seriesMutators.test.ts
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/saveSeries.ts \
        packages/HimalayaUI/frontend/src/lib/queue/mutators/deleteSeries.ts \
        packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts \
        packages/HimalayaUI/frontend/test/queue/seriesMutators.test.ts
git commit -m "$(cat <<'EOF'
Add saveSeries + deleteSeries mutators and registry wiring (#166)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: `applyRemoteToCache` handlers for the #166 kinds

**Files:**
- Modify: `frontend/src/lib/queue/applyRemoteToCache.ts`
- Test: `frontend/test/queue/applyRemoteToCache.series.test.ts` (create)

- [ ] **Step 1: Write the failing test**

Create `frontend/test/queue/applyRemoteToCache.series.test.ts`:

```typescript
/**
 * applyRemoteToCache — series foreign-event handlers (#166/#167/#168).
 *
 * series_created / series_recipe_updated carry no post_state (master plan
 * §5.2) and the SSE payload's recipe rows are id-less + replay-volatile (§11),
 * so the handler is invalidate-only. series_deleted removes the detail cache
 * and filters the listing.
 */
import { describe, it, expect } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import { queryKeys } from "../../src/queries";
import type { SseEvent } from "../../src/lib/queue/types";

describe("applyRemoteToCache: series #166 kinds", () => {
  it("series_created invalidates the detail + listing caches", () => {
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.series(5), { id: 5 });
    qc.setQueryData(queryKeys.seriesList, [{ id: 5 }]);
    const remote: SseEvent = {
      id: 1, kind: "series_created", entity_type: "series", entity_id: 5,
      payload: { title: "S", samples: [] },
    };
    applyRemoteToCache(remote, qc);
    expect(qc.getQueryState(queryKeys.series(5))?.isInvalidated).toBe(true);
    expect(qc.getQueryState(queryKeys.seriesList)?.isInvalidated).toBe(true);
  });

  it("series_recipe_updated invalidates the detail cache", () => {
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.series(5), { id: 5 });
    const remote: SseEvent = {
      id: 2, kind: "series_recipe_updated", entity_type: "series", entity_id: 5,
      payload: { samples: [] },
    };
    applyRemoteToCache(remote, qc);
    expect(qc.getQueryState(queryKeys.series(5))?.isInvalidated).toBe(true);
  });

  it("series_deleted removes the detail cache and filters the listing", () => {
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.series(5), { id: 5 });
    qc.setQueryData(queryKeys.seriesList, [{ id: 5 }, { id: 6 }]);
    const remote: SseEvent = {
      id: 3, kind: "series_deleted", entity_type: "series", entity_id: 5,
      payload: { id: 5 },
    };
    applyRemoteToCache(remote, qc);
    expect(qc.getQueryData(queryKeys.series(5))).toBeUndefined();
    expect(qc.getQueryData(queryKeys.seriesList)).toEqual([{ id: 6 }]);
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
npx vitest run test/queue/applyRemoteToCache.series.test.ts
```

Expected: FAIL — the unknown-kind `default` arm invalidates `peaks`/`indices`/`groups`, so `series(5)` is not invalidated and `seriesList` is untouched.

- [ ] **Step 3: Add the `case` arms**

In `frontend/src/lib/queue/applyRemoteToCache.ts`, add these arms inside the `switch (remote.kind)`, immediately after the `comparison_pinned` / `comparison_unpinned` arm (after line 278). These three arms use only structural `{ id: number }` typing — no new import is needed (the `Series` import is added in Task 8, which is the first arm that uses it):

```typescript
    case "series_created":
    case "series_recipe_updated": {
      // Neither kind carries post_state (master plan §5.2); the SSE payload's
      // series_samples entries are id-less and series_samples.id is
      // replay-volatile (§11), so there is no safe surgical splice.
      // Invalidate-only — the next read refetches the canonical projection.
      qc.invalidateQueries({ queryKey: queryKeys.series(id) });
      qc.invalidateQueries({ queryKey: queryKeys.seriesList });
      break;
    }
    case "series_deleted": {
      // Remove the detail cache — refetching a deleted resource 404s and
      // leaves stale `isError` state. Filter the id out of the listing.
      qc.removeQueries({ queryKey: queryKeys.series(id) });
      qc.setQueriesData<{ id: number }[]>(
        { queryKey: queryKeys.seriesList },
        (old) => (old ? old.filter((s) => s.id !== id) : old),
      );
      break;
    }
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
npx vitest run test/queue/applyRemoteToCache.series.test.ts
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts \
        packages/HimalayaUI/frontend/test/queue/applyRemoteToCache.series.test.ts
git commit -m "$(cat <<'EOF'
Add applyRemoteToCache handlers for the #166 series kinds

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: `series_plate_committed` dispatcher branch (#167)

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl` — branch + `_update_view_for_series_plate_committed!` + `_insert_series_member!`.
- Test: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Add this testset inside `@testset "Series routes"`:

```julia
    @testset "POST /api/series/{id}/commit — series_plate_committed commits the plate" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                sid = JSON3.read(HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:title => "Ramp",
                                     :samples => [Dict(:sample_id => 100, :position => 0)]))
                    ).body, Dict{Symbol, Any})[:id]

                # Pass an explicit snapshot so the route's _series_member_payload
                # uses it directly (no dependence on compute_member_snapshot
                # succeeding against the unanalyzed exposure fixture).
                commitBody = JSON3.write(Dict(:members => [
                    Dict(:exposure_id => 1000, :display_order => 0,
                         :snapshot => Dict(:effective_peaks => [],
                                           :confirmed_index => nothing,
                                           :analysis_inputs_hash => nothing)),
                ]))
                resp = HTTP.post("$base/api/series/$sid/commit",
                    ["X-Username" => "alice", "Content-Type" => "application/json"], commitBody)
                @test resp.status == 200
                got = JSON3.read(resp.body, Dict{Symbol, Any})
                @test got[:state] == "committed"
                @test startswith(got[:content_hash], "sha256:")   # hashed from the plate
                @test length(got[:members]) == 1
                @test got[:members][1]["exposure_id"] == 1000
                @test length(got[:samples]) == 1                  # recipe untouched

                # rebuild_views_from_log! round-trip from empty.
                DBInterface.execute(db, "DELETE FROM series_members WHERE series_id = ?", [sid])
                DBInterface.execute(db, "DELETE FROM series_samples WHERE series_id = ?", [sid])
                DBInterface.execute(db, "DELETE FROM series WHERE id = ?", [sid])
                HimalayaUI.rebuild_views_from_log!(db, sid; entity_type = "series")
                refold = HimalayaUI.fetch_series_with_plate(db, sid)
                @test refold[:state] == "committed"
                @test length(refold[:members]) == 1
                @test refold[:content_hash] == got[:content_hash]   # hash is fold-stable

                # SSE layer: series_plate_committed is the one series event
                # carrying a post_state envelope. Create + commit a fresh
                # series through the in-process subscriber and assert the
                # frame carries entity_type='series' and a post_state field.
                sid2 = JSON3.read(HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:title => "Frame check",
                        :samples => [Dict(:sample_id => 100, :position => 0)]))
                    ).body, Dict{Symbol, Any})[:id]
                frames = _capture_series_sse("series_plate_committed") do
                    HTTP.post("$base/api/series/$sid2/commit",
                        ["X-Username" => "alice", "Content-Type" => "application/json"],
                        JSON3.write(Dict(:members => [
                            Dict(:exposure_id => 1000, :display_order => 0,
                                 :snapshot => Dict(:effective_peaks => [],
                                                   :confirmed_index => nothing,
                                                   :analysis_inputs_hash => nothing))])))
                end
                @test length(frames) == 1
                @test occursin("\"entity_type\":\"series\"", frames[1])
                @test occursin("\"post_state\"", frames[1])
            end
            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the focused backend command. Expected: FAIL — no `series_plate_committed` branch, so the commit route's dispatcher no-ops; `got[:state]` stays `"draft"` and `members` is empty.

- [ ] **Step 3: Add the dispatcher branch**

In `events.jl`, inside `update_view_for_event!`, add after the `series_deleted` branch (added in Task 1):

```julia
    if kind == "series_plate_committed"
        return _update_view_for_series_plate_committed!(db, entity_id, payload, event_id)
    end
```

- [ ] **Step 4: Add the helper functions**

In `events.jl`, after `_update_view_for_series_recipe_updated!` (added in Task 2), add:

```julia
# Helper: INSERT one series_members row. Mirrors `_insert_comparison_member!`
# (the series_members and comparison_members column shapes are identical) —
# reuses `_member_field` / `_member_json`. Members in a series_plate_committed
# payload carry NO ids; the dispatcher mints fresh PKs.
function _insert_series_member!(db, series_id, m, user_id, now_str)
    snap = _member_json(m, :snapshot)
    snap === nothing &&
        error("series member missing required `snapshot` field")
    DBInterface.execute(db,
        """INSERT INTO series_members
             (series_id, exposure_id, display_order, band_height, y_offset,
              normalization, color_override, label_override,
              q_window_min, q_window_max, peak_display, snapshot,
              created_by, created_at)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        [Int(series_id),
         _member_field(m, :exposure_id),
         Int(_member_field(m, :display_order)),
         Float64(_member_field(m, :band_height; default=1.0)),
         Float64(_member_field(m, :y_offset; default=0.0)),
         String(_member_field(m, :normalization; default="none")),
         _member_field(m, :color_override),
         _member_field(m, :label_override),
         _member_field(m, :q_window_min),
         _member_field(m, :q_window_max),
         _member_json(m, :peak_display),
         snap,
         user_id, now_str])
    nothing
end

"""
    _update_view_for_series_plate_committed!(db, entity_id, payload, event_id)

`series_plate_committed` dispatcher (#167). Pure-replace the plate:
`DELETE` every `series_members` row, then `INSERT` the full payload member
list (members carry no ids — `_insert_series_member!` mints them). Sets
`state='committed'` and computes `content_hash` from the plate via
`compute_series_content_hash` (which excludes the recipe — master plan §5.1).
Never touches `series_samples`.
"""
function _update_view_for_series_plate_committed!(db, entity_id, payload, event_id)
    sid     = Int(entity_id)
    user_id = user_id_for_event(db, event_id)
    now_str = comparison_now_iso()

    DBInterface.execute(db, "DELETE FROM series_members WHERE series_id = ?", [sid])
    members = haskey(payload, :members) ? payload.members : []
    for m in members
        _insert_series_member!(db, sid, m, user_id, now_str)
    end

    new_hash = compute_series_content_hash(db, sid)
    DBInterface.execute(db,
        "UPDATE series SET content_hash = ?, state = 'committed', updated_at = ? WHERE id = ?",
        [new_hash, now_str, sid])
    return nothing
end
```

- [ ] **Step 5: Run the test to verify it passes**

Run the focused backend command. Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/events.jl packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(cat <<'EOF'
Add the series_plate_committed dispatcher branch (#167)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: `commitSeriesPlate` mutator + `applyRemoteToCache` post_state handler (#167)

**Files:**
- Create: `frontend/src/lib/queue/mutators/commitSeriesPlate.ts`
- Modify: `frontend/src/lib/queue/mutatorRegistry.ts`
- Modify: `frontend/src/lib/queue/applyRemoteToCache.ts`
- Test: `frontend/test/queue/seriesMutators.test.ts`, `frontend/test/queue/applyRemoteToCache.series.test.ts`

- [ ] **Step 1: Write the failing tests**

Append to `frontend/test/queue/seriesMutators.test.ts`:

```typescript
import { commitSeriesPlateMutator } from "../../src/lib/queue/mutators/commitSeriesPlate";

describe("commitSeriesPlateMutator", () => {
  it("synthesizeFromSse builds a full Series from the post_state envelope", () => {
    const post = fullSeries(7);
    post.state = "committed";
    const synth = commitSeriesPlateMutator.synthesizeFromSse?.(
      { id: 99, kind: "series_plate_committed", entity_type: "series",
        entity_id: 7, payload: { members: [] }, post_state: post },
      { event_id: 99, client_op_id: "op1", analysis_inputs_hash: undefined });
    expect(synth?.state).toBe("committed");
    expect(synth?.id).toBe(7);
  });
});
```

First add a `Series` type import at the top of `frontend/test/queue/applyRemoteToCache.series.test.ts` (alongside the existing `SseEvent` import):

```typescript
import type { Series } from "../../src/api";
```

Then append, inside the `describe` block:

```typescript
  it("series_plate_committed splices post_state into the detail cache", () => {
    const qc = new QueryClient();
    // A structurally-complete Series — post_state IS the full
    // fetch_series_with_plate projection, so the cache write is a real
    // round-trip, not a partial-object cast.
    const post: Series = {
      id: 5, title: "S", description: null, content_hash: "sha256:x",
      created_by: 1, created_at: null, updated_at: null,
      forked_from_id: null, forked_at_hash: null, forked_from_title: null,
      view_grouping_mode: null, view_show_peak_ticks: null,
      view_show_peak_labels: null, ordering_variable: null,
      order_rule: "manual", state: "committed", members: [], samples: [],
    };
    const remote: SseEvent = {
      id: 4, kind: "series_plate_committed", entity_type: "series", entity_id: 5,
      payload: { members: [] },
      post_state: post,
    };
    applyRemoteToCache(remote, qc);
    expect(qc.getQueryData(queryKeys.series(5))).toEqual(post);
  });
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
npx vitest run test/queue/seriesMutators.test.ts test/queue/applyRemoteToCache.series.test.ts
```

Expected: FAIL — `Cannot find module './mutators/commitSeriesPlate'`; `series_plate_committed` hits the `default` arm so `series(5)` is not written.

- [ ] **Step 3: Create `commitSeriesPlate.ts`**

Create `frontend/src/lib/queue/mutators/commitSeriesPlate.ts`:

```typescript
/**
 * commitSeriesPlate mutator (#167). Commits the plate via
 * `POST /api/series/:id/commit` — the old "submit". No optimistic write
 * (spinner). `series_plate_committed` is the one series event carrying a
 * `post_state` envelope (the full `fetch_series_with_plate` projection), so
 * `synthesizeFromSse` can return a complete Series on the SSE-wins path.
 *
 * A 409 (content_hash conflict) surfaces as a generic `ApiError`; the typed
 * conflict modal is I3.5b's concern.
 */
import * as api from "../../../api";
import type {
  AuthOpts, Series, SeriesMemberInput, CommitSeriesPlateBody,
} from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";

export interface CommitSeriesPlateInput {
  id: number;
  members: SeriesMemberInput[];
  expected_content_hash?: string;
}

interface CommitSeriesPlateScope {
  username: string | undefined;
  clientId: string;
}

function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

export const commitSeriesPlateMutator: Mutator<CommitSeriesPlateInput, CommitSeriesPlateScope, Series> = {
  kind: "series_commit",
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: (p) => {
    // Annotate with the shared CommitSeriesPlateBody type — do not re-declare
    // it inline (drift risk).
    const body: CommitSeriesPlateBody = { members: p.members };
    if (p.expected_content_hash !== undefined) {
      body.expected_content_hash = p.expected_content_hash;
    }
    return api.commitSeriesPlate(p.id, body, buildAuthOpts(p));
  },
  onSuccess: (_p, response, qc) => {
    const looksFull = Array.isArray((response as { members?: unknown }).members)
      && typeof (response as { state?: unknown }).state === "string";
    if (looksFull) {
      qc.setQueryData(queryKeys.series(response.id), response);
    } else if (typeof response?.id === "number") {
      qc.invalidateQueries({ queryKey: queryKeys.series(response.id) });
    }
    qc.invalidateQueries({ queryKey: queryKeys.seriesList });
  },
  synthesizeFromSse: (remote, _base) => {
    // post_state IS the full fetch_series_with_plate projection (master plan
    // §5.2). Return it RAW — do NOT spread `_base` (QueueResponseMeta) in:
    // `onSuccess`'s looksFull branch writes this object straight to the detail
    // cache via setQueryData, and merging event_id / client_op_id /
    // analysis_inputs_hash would pollute the cached Series row and diverge
    // from the clean post_state splice in the applyRemoteToCache arm. `_base`
    // is intentionally unused (underscore-prefixed for noUnusedParameters).
    if (remote.post_state != null) {
      return remote.post_state as Series;
    }
    return undefined;
  },
};
```

- [ ] **Step 4: Wire the registry**

In `frontend/src/lib/queue/mutatorRegistry.ts`, add the import after the `deleteSeriesMutator` import:

```typescript
import { commitSeriesPlateMutator } from "./mutators/commitSeriesPlate";
```

In `resolveMutator`, after the `case "series_delete":` arm (added in Task 5):

```typescript
    case "series_commit":
      return commitSeriesPlateMutator;
```

In `resolveMutatorForEvent`, after the `case "series_deleted":` arm (added in Task 5):

```typescript
    case "series_plate_committed":
      return commitSeriesPlateMutator;
```

- [ ] **Step 5: Add the `applyRemoteToCache` arm**

In `frontend/src/lib/queue/applyRemoteToCache.ts`, first add `Series` to the type-import block (lines 3–6 — `Peak, GroupEntry, … Comparison,`):

```typescript
import type {
  Peak, GroupEntry, Exposure, Sample, SampleMessage,
  ComparisonMessage, ComparisonSummary, Comparison, Series,
} from "../../api";
```

Then add this arm after the `series_deleted` arm (added in Task 6):

```typescript
    case "series_plate_committed": {
      // The one series event carrying a post_state envelope (master plan
      // §5.2): post_state IS the fetch_series_with_plate projection. Splice it
      // straight into the detail cache; invalidate the listing (denormalised
      // member_count / has_stale_members / last_event_at fields).
      if (remote.post_state != null) {
        qc.setQueryData(queryKeys.series(id), remote.post_state as Series);
      } else {
        qc.invalidateQueries({ queryKey: queryKeys.series(id) });
      }
      qc.invalidateQueries({ queryKey: queryKeys.seriesList });
      break;
    }
```

- [ ] **Step 6: Run the tests to verify they pass**

```bash
npx vitest run test/queue/seriesMutators.test.ts test/queue/applyRemoteToCache.series.test.ts
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/commitSeriesPlate.ts \
        packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts \
        packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts \
        packages/HimalayaUI/frontend/test/queue/seriesMutators.test.ts \
        packages/HimalayaUI/frontend/test/queue/applyRemoteToCache.series.test.ts
git commit -m "$(cat <<'EOF'
Add commitSeriesPlate mutator + series_plate_committed cache handler (#167)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: `series_pinned` / `series_unpinned` dispatcher branches (#168)

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl` — two inline branches in `update_view_for_event!`.
- Test: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Add this testset inside `@testset "Series routes"`:

```julia
    @testset "POST/DELETE /api/series/{id}/pin — series pins" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            DBInterface.execute(db, "INSERT INTO series (id, title, state) VALUES (5, 'S5', 'draft')")
            with_test_server(db) do port, base
                # Pin → series_pins row written; GET reflects it.
                respPin = HTTP.post("$base/api/series/5/pin", ["X-Username" => "alice"])
                @test respPin.status == 200
                pins = JSON3.read(HTTP.get("$base/api/users/me/series-pins",
                    ["X-Username" => "alice"]).body)
                @test pins == [5]

                # Unpin → series_pins row removed.
                respUnpin = HTTP.delete("$base/api/series/5/pin", ["X-Username" => "alice"])
                @test respUnpin.status == 200
                pins2 = JSON3.read(HTTP.get("$base/api/users/me/series-pins",
                    ["X-Username" => "alice"]).body)
                @test pins2 == []

                # rebuild_views_from_log! round-trip: pin events live on
                # entity_type='user'. Re-fold for the user → series_pins
                # rebuilt. (Pin then unpin → empty; pin again to assert a row.)
                HTTP.post("$base/api/series/5/pin", ["X-Username" => "alice"])
                uid = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id FROM users WHERE username = 'alice'"))[1].id
                DBInterface.execute(db, "DELETE FROM series_pins WHERE user_id = ?", [uid])
                HimalayaUI.rebuild_views_from_log!(db, uid; entity_type = "user")
                refold = Tables.rowtable(DBInterface.execute(db,
                    "SELECT series_id FROM series_pins WHERE user_id = ?", [uid]))
                @test length(refold) == 1
                @test refold[1].series_id == 5
            end
            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the focused backend command. Expected: FAIL — no `series_pinned` / `series_unpinned` branches, so the pin routes write a `user_actions` row but no `series_pins` row; `pins` comes back `[]`.

- [ ] **Step 3: Add the dispatcher branches**

In `events.jl`, inside `update_view_for_event!`, add after the `series_plate_committed` branch (added in Task 7):

```julia
    # Per-user series pins (#168 / I2.5). Stored with entity_type='user',
    # entity_id=user_id (the comparison_pinned precedent) — the affected series
    # rides in the payload as `series_id`. Five-layer: no post_state, no
    # mutator. The dispatcher derives user_id by joining user_actions on
    # event_id, exactly as the comparison pin branches do.
    if kind == "series_pinned"
        DBInterface.execute(db,
            """INSERT OR REPLACE INTO series_pins (user_id, series_id, pinned_at)
               VALUES ((SELECT user_id FROM user_actions WHERE id = ?),
                       ?, CURRENT_TIMESTAMP)""",
            [event_id, Int(payload.series_id)])
        return nothing
    end
    if kind == "series_unpinned"
        DBInterface.execute(db,
            """DELETE FROM series_pins
               WHERE user_id = (SELECT user_id FROM user_actions WHERE id = ?)
                 AND series_id = ?""",
            [event_id, Int(payload.series_id)])
        return nothing
    end
```

- [ ] **Step 4: Run the test to verify it passes**

Run the focused backend command. Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/events.jl packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(cat <<'EOF'
Add the series_pinned / series_unpinned dispatcher branches (#168)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: `applyRemoteToCache` pin handlers + full-suite acceptance (#168)

**Files:**
- Modify: `frontend/src/lib/queue/applyRemoteToCache.ts`
- Test: `frontend/test/queue/applyRemoteToCache.series.test.ts`

- [ ] **Step 1: Write the failing test**

Append to `frontend/test/queue/applyRemoteToCache.series.test.ts`, inside the `describe` block:

```typescript
  it("series_pinned / series_unpinned invalidate the pins cache", () => {
    for (const kind of ["series_pinned", "series_unpinned"] as const) {
      const qc = new QueryClient();
      qc.setQueryData(queryKeys.seriesPins, [5]);
      const remote: SseEvent = {
        id: 5, kind, entity_type: "user", entity_id: 1,
        payload: { series_id: 5 },
      };
      applyRemoteToCache(remote, qc);
      expect(qc.getQueryState(queryKeys.seriesPins)?.isInvalidated).toBe(true);
    }
  });
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
npx vitest run test/queue/applyRemoteToCache.series.test.ts
```

Expected: FAIL — `series_pinned` / `series_unpinned` hit the `default` arm, which does not touch `seriesPins`.

- [ ] **Step 3: Add the `case` arm**

In `frontend/src/lib/queue/applyRemoteToCache.ts`, add after the `series_plate_committed` arm (added in Task 8):

```typescript
    case "series_pinned":
    case "series_unpinned": {
      // Pin/unpin fan out cross-tab. The seriesPins cache is global per-tab
      // (the current user's pin set); the SSE self-echo filter discards the
      // originating tab's own frame. Invalidate so the next read gets the
      // canonical list — mirrors comparison_pinned / comparison_unpinned.
      qc.invalidateQueries({ queryKey: queryKeys.seriesPins });
      break;
    }
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
npx vitest run test/queue/applyRemoteToCache.series.test.ts
```

Expected: PASS.

- [ ] **Step 5: Run the full frontend gate**

From `packages/HimalayaUI/frontend/`:

```bash
npm run build && npx vitest run
```

Expected: `tsc --noEmit` + vite build green; the full Vitest suite green.

- [ ] **Step 6: Run the full backend suite**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-test.out
tail -50 /tmp/jl-test.out
```

Expected: green. (Two known flaky tests — the fast-skip P99 latency test and a `GET /api/health` port race — may fail intermittently regardless of this change; re-run once to confirm they are not regressions.)

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts \
        packages/HimalayaUI/frontend/test/queue/applyRemoteToCache.series.test.ts
git commit -m "$(cat <<'EOF'
Add applyRemoteToCache pin handlers; close the series event-kind cluster (#168)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Acceptance criteria (issues #166 / #167 / #168)

- [ ] **#166** — `series_created`, `series_recipe_updated`, `series_deleted` each have a dispatcher branch, a frontend mutator/cache handler, registry cases, and paired contract + `rebuild_views_from_log!` round-trip tests.
- [ ] The view-producing branches are **pure-replace** (parent upsert + child `DELETE`/`INSERT` from a full-snapshot payload) — verified by Task 1/2's fold-from-empty round-trips.
- [ ] `series_created` folds correctly from an empty view (Task 1, Step 1 round-trip).
- [ ] All four non-pin `series_*` events carry `entity_type='series'` on the wire (already set by `routes_series.jl`; unchanged here).
- [ ] **#167** — `series_plate_committed` commits the plate, sets `state='committed'`, computes `content_hash` from the plate; its `applyRemoteToCache` handler splices the `post_state` envelope; `commitSeriesPlateMutator` defines `synthesizeFromSse`.
- [ ] **#168** — `series_pinned` / `series_unpinned` write `series_pins` under `entity_type='user'`; five-layer (no mutator, no post_state); pins replay correctly under `entity_type='user'` (Task 9 round-trip).
- [ ] `resolveMutatorForEvent` has a case for every non-pin series kind; pins fall through to `undefined`.
- [ ] Recipe rows are keyed on `(series_id, position)`, never `series_samples.id` (master plan §11) — the frontend handlers are invalidate-only / post_state-splice, never id-splice.
- [ ] The SSE `broadcast_event!` layer is asserted, not assumed — `series_created` and `series_plate_committed` frames are captured via the in-process subscriber (Tasks 1, 7).
- [ ] Full Julia suite + `npm run build` + Vitest green (Task 10).

## Out of scope

- The native-series round-trip test #170 (I2.7) — a *separate* dedicated test; the per-kind round-trips here satisfy the six-layer "round-trip" layer, not #170.
- The comparison→series data migration #171 (I3.1).
- Series read hooks (`useSeriesList`, `useSeries`) and the `SeriesSummary` listing type — I3.3.
- Builder mutation wiring / negative-id optimistic placeholders / the typed 409 conflict modal — I3.5b.
- The historical `comparison_*` dispatcher branches — frozen, never aliased.
