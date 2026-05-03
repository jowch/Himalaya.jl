# Server-Reconciliation Mutation Queue Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement [docs/superpowers/specs/2026-05-02-mutation-queue-design.md](../specs/2026-05-02-mutation-queue-design.md): a server-reconciliation mutation queue that eliminates the autoReanalyze double round-trip, replaces invalidate-and-refetch with cache-merged optimistic state, and resolves concurrent conflicts via replay-as-rerun.

**Architecture:** Four milestones, each independently mergeable. M0 lands schema + idempotency primitives + `analyze_exposure!` fast-skip + `analyze_run` broadcast suppression (mostly backend). M1 lands the queue framework with deferred-promise mutator resolution and replay coordinator (frontend infrastructure, no consumers). M2 migrates all 16 mutations to the queue, with the **peak-op slice as one atomic PR** (cannot be subdivided). M3 mechanical cleanup.

**Tech stack:** Julia 1.9+, SQLite.jl, Oxygen.jl, HTTP.jl (SSE), JSON3; React 18, TanStack Query v5, Vitest, Playwright.

**Read first:**
- [docs/superpowers/specs/2026-05-02-mutation-queue-design.md](../specs/2026-05-02-mutation-queue-design.md) — design rationale, 14 architectural decisions, fallback triggers
- [docs/event-log.md](../../event-log.md) — `apply_event!` dispatcher contract, hash invariants, SSE semantics; the queue's server-side primitives compose with these unchanged
- [docs/superpowers/specs/2026-05-02-sse-client-id-design.md](../specs/2026-05-02-sse-client-id-design.md) — PR #32, per-tab `client_id`. Direct prerequisite for the queue's SSE confirmation primitive.
- [docs/superpowers/plans/2026-05-01-multiplayer-instrumentation.md](2026-05-01-multiplayer-instrumentation.md) §"Migration Architecture" — idempotency contract, ALTER TABLE try/catch precedent, sentinel patterns

**Prerequisites:**
- PR #32 must be merged before M0 begins. The plan extends `apply_event!`'s `client_id` plumbing pattern to a parallel `client_op_id` parameter; the file diff is meaningfully simpler when `client_id` is in main rather than a sibling worktree.

---

## File Map

**M0: Backend (mostly) + minimal frontend**

Modified backend:
- `packages/HimalayaUI/src/db.jl` — append `client_op_id TEXT` to `user_actions` `CREATE TABLE` DDL; new `idempotent_responses` table DDL; two ALTER TABLEs in `migrate_schema!`'s `stmts` array; new `idempotent_responses` table in fresh-DB schema
- `packages/HimalayaUI/src/actions.jl` — new `get_client_op_id(req)` helper (mirrors `get_client_id` from PR #32)
- `packages/HimalayaUI/src/events.jl` — `apply_event!` extracts + persists `client_op_id`; `broadcast_event!` adds `client_op_id` and `ts` to JSON; new `with_idempotency` wrapper + `OP_LOCKS` registry; broadcast suppression for both-skipped `analyze_run` events
- `packages/HimalayaUI/src/pipeline.jl` — `analyze_exposure!` gains `trace_known_unchanged` keyword arg; refactored to do DB-only reads first and skip `hash_trace_file`/`load_dat` on hot path; new helpers `any_add_curations`, `hash_peak_set_from_db`, `read_trace_hash`, `read_inputs_hash`, `count_auto_peaks`, `count_indices`
- `packages/HimalayaUI/src/events.jl` — also gains `gc_idempotent_responses!` function (sweeps `idempotent_responses` rows older than TTL + corresponding `OP_LOCKS` entries)
- `packages/HimalayaUI/src/server.jl` — `start_gc_timer!` helper called from `serve(db)` to invoke the sweep every ~30 minutes

Modified frontend:
- `packages/HimalayaUI/frontend/src/api.ts` — `AuthOpts.clientOpId`, `X-Client-Op-Id` header on mutations
- `packages/HimalayaUI/frontend/src/queries.ts` — `authOpts(username, clientId, clientOpId?)`; mutators that need a fresh op id mint one per call (NOT module-level — unlike `client_id` which is per-tab, `client_op_id` is per-mutation)

New frontend:
- `packages/HimalayaUI/frontend/src/lib/clientOpId.ts` — `newClientOpId()` mints `crypto.randomUUID()` per call

New tests:
- `packages/HimalayaUI/test/test_idempotency.jl` (new file; `include` in `runtests.jl`) — `with_idempotency` cache hit/miss, concurrent retry, failure-not-cached, multi-event route cache hit
- `packages/HimalayaUI/test/test_fast_skip.jl` (new file; `include` in `runtests.jl`) — `analyze_exposure!` fast-skip path, broadcast suppression, slow path correctness with add-curations
- `packages/HimalayaUI/frontend/test/clientOpId.test.ts` — mint generates unique UUIDs

Modified tests:
- `packages/HimalayaUI/test/test_db.jl` — fresh-DB `idempotent_responses` table presence, legacy-schema `client_op_id` migration, partial-index presence
- `packages/HimalayaUI/test/test_actions.jl` — `get_client_op_id` extraction (file created in PR #32)
- `packages/HimalayaUI/test/test_events.jl` — `apply_event!` writes/reads `client_op_id` column
- `packages/HimalayaUI/test/test_sse.jl` — `broadcast_event!` frame includes `client_op_id` and `ts`; both-skipped `analyze_run` does NOT broadcast
- `packages/HimalayaUI/test/test_pipeline.jl` — existing `analyze_exposure!` tests stay green; verify the refactor doesn't regress slow-path callers (reingest, scheduled scan, manual reanalyze)
- `packages/HimalayaUI/frontend/test/api.test.ts` — `X-Client-Op-Id` header on mutations + GET-omit case

Modified docs:
- `docs/event-log.md` — new section on `client_op_id` semantics, `idempotent_responses` table, broadcast suppression rule
- `CLAUDE.md` §"HimalayaUI gotchas" — note the per-mutation UUID pattern (vs PR #32's per-tab UUID)

**M1: Frontend infrastructure (no consumers)**

New frontend:
- `packages/HimalayaUI/frontend/src/lib/queue/types.ts` — `OpKind`, `Mutator<T,R>`, `RollbackContext`, `PendingDeferred<T>`, `SseEvent`
- `packages/HimalayaUI/frontend/src/lib/queue/deferred.ts` — `pendingDeferreds` Map, `makeDeferred`, registry helpers
- `packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts` — `handleRemoteEvent`, per-kind `applyRemoteToCache`
- `packages/HimalayaUI/frontend/src/lib/queue/persistence.ts` — sessionStorage mirror, rehydrate flow
- `packages/HimalayaUI/frontend/src/lib/queue/hooks.ts` — `useExposureHasPendingPeakOps`, `useQueueOpStatus`
- `packages/HimalayaUI/frontend/src/lib/queue/useQueueMutation.ts` — framework wrapper: deferred-promise resolution, optimistic-effect plumbing, AbortSignal handling, sessionStorage persistence wiring, failure-class routing
- `packages/HimalayaUI/frontend/src/lib/queue/errors.ts` — `isValidationError`, `isInfrastructureError`, per-kind error message lookup
- `packages/HimalayaUI/frontend/src/lib/queue/index.ts` — barrel export

Modified frontend:
- `packages/HimalayaUI/frontend/src/lib/sseSubscriber.ts` — wire the queue's `handleRemoteEvent` ahead of the existing fallback path. The `client_id`-based skip filter remains as the fallback for non-queue echoes (system events, pre-feature events).

New tests:
- `packages/HimalayaUI/frontend/test/queue/replayCoordinator.test.ts` — pop-on-client-op-id, rollback-apply-replay sequence, MutationCache iteration order
- `packages/HimalayaUI/frontend/test/queue/persistence.test.ts` — rehydrate cycle, schema-version drop with toast
- `packages/HimalayaUI/frontend/test/queue/hooks.test.ts` — `useExposureHasPendingPeakOps` reflects pending mutations
- `packages/HimalayaUI/frontend/test/queue/deferred.test.ts` — registry mint/resolve/leak-on-abort
- `packages/HimalayaUI/frontend/test/queue/useQueueMutation.test.ts` — wrapper-level integration: per-call clientOpId minting, optimistic context plumbing, error-class routing, AbortSignal threading
- `packages/HimalayaUI/frontend/test/queue/helpers.ts` — test helpers (`makeFakeMutation`, `remoteForeignEvent`)

**M2: Mutator migration**

The peak-op slice is one atomic PR (Task M2.2 below). Other slices (M2.1, M2.3, M2.4, M2.5) land independently in any order.

Modified backend (peak-op slice):
- `packages/HimalayaUI/src/routes_peaks.jl` — three curation routes wrap body in `with_idempotency`, use the `defer_broadcast=true` kwarg from M0.6 inside the transaction, call `analyze_exposure!(db, id, dir; trace_known_unchanged=true)` synchronously, then emit one enriched broadcast with `post_state` populated
- (No new changes to `events.jl` — M0.6 already shipped `defer_broadcast` and `post_state`; M2 strictly consumes them)

Modified frontend:
- `packages/HimalayaUI/frontend/src/queries.ts` — replace each `useMutation` with a queue-shaped mutator hook; delete `autoReanalyze` chain; `useReanalyzeExposure` becomes null-mutator queue op
- `packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx` — gate on `useExposureHasPendingPeakOps`; reduce debounce to ~150ms
- `packages/HimalayaUI/frontend/src/components/SpeculativeBuilder.tsx` — `useSpeculativeSnap` gates on `useExposureHasPendingPeakOps`; modal renders last response with "updating to latest…" subtext during gate
- `packages/HimalayaUI/frontend/src/components/{ChatCard,SamplePropertyEditor,...}.tsx` — switch to queue-shaped hooks for tags, messages, sample updates, status (per migration order)
- `packages/HimalayaUI/frontend/src/components/ui/Toast.tsx` (new or extension of existing toast) — Validation-class error toast
- `packages/HimalayaUI/frontend/src/components/InfrastructureBanner.tsx` (new) — persistent banner for retry-exceeded ops

New tests:
- `packages/HimalayaUI/frontend/test/queue/peakAddMutator.test.ts` and parallel files per mutator — round-trip optimistic effect → request → confirmation → cache settled
- `packages/HimalayaUI/frontend/e2e/multiplayer-replay-rerun.spec.ts` — two-context Playwright test: User A and User B both add peaks; both succeed; both exclude same peak; replay-rerun resolves invisibly
- `packages/HimalayaUI/test/test_routes_peaks.jl` — synchronous reanalyze inside curation handlers; response shape includes `event_id`, `view_row_id`, `analysis_inputs_hash`, and the inserted curation row (under key `peak`)

**M3: Cleanup**

- `packages/HimalayaUI/frontend/src/lib/sseSubscriber.ts` — drop now-redundant `client_id` skip filter (or retain as defense-in-depth — decide during implementation)
- `packages/HimalayaUI/frontend/src/queries.ts` — audit remaining `qc.invalidateQueries` calls; replace with `setQueryData` from response where applicable
- `docs/event-log.md` — document the queue's relationship to the dispatcher contract
- `CLAUDE.md` §"HimalayaUI gotchas" — queue invariants, `useExposureHasPendingPeakOps` pattern, mutator framework conventions

---

## Migration Architecture

This plan introduces two schema changes (one column, one new table) plus an `analyze_exposure!` refactor. All schema work follows the idempotency contract from [Plan 7's migration architecture](2026-05-01-multiplayer-instrumentation.md):

- `migrate_schema!` must be safe on fresh / pre-migration / partially-migrated / already-migrated DBs.
- `ALTER TABLE ADD COLUMN` is wrapped in `try ... catch end` (SQLite raises on duplicate column; precedent at [db.jl:142-160](../../../packages/HimalayaUI/src/db.jl)).
- `CREATE TABLE IF NOT EXISTS` and `CREATE INDEX IF NOT EXISTS` for the response cache table — same idempotency without a try/catch wrapper.
- No backfill needed: existing `user_actions` rows get NULL `client_op_id` (correct semantic value for pre-feature events). The `idempotent_responses` table starts empty.
- Forward-only; rollback is via DB backup before deploy.

The `analyze_exposure!` refactor is behaviorally identical for callers that don't pass `trace_known_unchanged=true` — existing callers (`reingest`, scheduled scans, manual `useReanalyzeExposure`) preserve current behavior. The fast-skip path activates only when callers opt in by passing `trace_known_unchanged=true`. Today (post-M0) only Task M2.2's peak-op handlers do; other callers stay on the slow path.

---

## Task M0.1: `client_op_id` schema + idempotent_responses table

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl`
- Modify: `packages/HimalayaUI/test/test_db.jl`

- [ ] **Step 1: Write failing tests**

Add to `packages/HimalayaUI/test/test_db.jl`:

```julia
@testset "user_actions.client_op_id column exists on fresh DB" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        cols = Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(user_actions)"))
        @test any(c -> c.name == "client_op_id", cols)
        SQLite.close(db)
    end
end

@testset "idempotent_responses table exists on fresh DB" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        tables = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table' AND name='idempotent_responses'"))
        @test length(tables) == 1
        cols = Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(idempotent_responses)"))
        @test Set(c.name for c in cols) == Set(["client_op_id", "status_code", "body", "created_at"])
        # client_op_id is the PK
        @test any(c -> c.name == "client_op_id" && c.pk == 1, cols)
        SQLite.close(db)
    end
end

@testset "open_db adds client_op_id to legacy user_actions schema" begin
    mktempdir() do tmp
        path = joinpath(tmp, "test.db")
        # Build a legacy user_actions table without client_op_id (post-PR-32 schema).
        # Mirror current create_schema! DDL minus the new column.
        db = SQLite.DB(path)
        DBInterface.execute(db, """
            CREATE TABLE user_actions (
                id              INTEGER PRIMARY KEY,
                user_id         INTEGER REFERENCES users(id) ON DELETE SET NULL,
                timestamp       DATETIME DEFAULT CURRENT_TIMESTAMP,
                action          TEXT,
                entity_type     TEXT,
                entity_id       INTEGER,
                note            TEXT,
                payload         TEXT,
                undoes_event_id INTEGER REFERENCES user_actions(id),
                client_id       TEXT
            )
        """)
        SQLite.close(db)
        db = HimalayaUI.open_db(path)
        cols = Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(user_actions)"))
        @test any(c -> c.name == "client_op_id", cols)
        SQLite.close(db)
    end
end

@testset "open_db creates idempotent_responses on legacy DB" begin
    mktempdir() do tmp
        path = joinpath(tmp, "test.db")
        db = SQLite.DB(path)
        # No idempotent_responses table.
        DBInterface.execute(db, "CREATE TABLE foo (x INTEGER)")
        SQLite.close(db)
        db = HimalayaUI.open_db(path)
        tables = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table' AND name='idempotent_responses'"))
        @test length(tables) == 1
        SQLite.close(db)
    end
end

@testset "client_op_id partial index present" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        idx = Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='index' AND name='idx_events_by_client_op_id'"))
        @test length(idx) == 1
        @test occursin("WHERE client_op_id IS NOT NULL", String(idx[1].sql))
        SQLite.close(db)
    end
end
```

Run: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'` — expect failures.

- [ ] **Step 2: Add column, table, and index**

In `packages/HimalayaUI/src/db.jl`:

1. Add `client_op_id TEXT` to the `user_actions` `CREATE TABLE` DDL inside `create_schema!`.

2. Add `CREATE TABLE IF NOT EXISTS idempotent_responses (...)` to `create_schema!`:

```julia
DBInterface.execute(db, """
    CREATE TABLE IF NOT EXISTS idempotent_responses (
        client_op_id  TEXT PRIMARY KEY,
        status_code   INTEGER NOT NULL,
        body          TEXT NOT NULL,
        created_at    DATETIME DEFAULT CURRENT_TIMESTAMP
    )
""")
DBInterface.execute(db, """
    CREATE INDEX IF NOT EXISTS idx_idempotent_responses_created
      ON idempotent_responses(created_at)
""")
```

3. Add to `migrate_schema!`'s `stmts` array (existing try/catch loop handles duplicate column):

```julia
"ALTER TABLE user_actions ADD COLUMN client_op_id TEXT",
"CREATE INDEX IF NOT EXISTS idx_events_by_client_op_id ON user_actions(client_op_id) WHERE client_op_id IS NOT NULL",
```

4. Add a separate idempotent statement for the response cache (creates the table if missing on legacy DBs):

```julia
"""CREATE TABLE IF NOT EXISTS idempotent_responses (
    client_op_id  TEXT PRIMARY KEY,
    status_code   INTEGER NOT NULL,
    body          TEXT NOT NULL,
    created_at    DATETIME DEFAULT CURRENT_TIMESTAMP
)""",
"CREATE INDEX IF NOT EXISTS idx_idempotent_responses_created ON idempotent_responses(created_at)",
```

- [ ] **Step 3: Verify**

Run the test above; all four testsets pass. Run full suite — no regressions because new columns are nullable and new table is unread by existing code.

---

## Task M0.2: `client_op_id` plumbing through `apply_event!` and SSE frame

**Files:**
- Modify: `packages/HimalayaUI/src/actions.jl` — new `get_client_op_id(req)` helper
- Modify: `packages/HimalayaUI/src/events.jl` — `apply_event!` reads + persists; `broadcast_event!` adds field
- Modify: `packages/HimalayaUI/test/test_actions.jl`
- Modify: `packages/HimalayaUI/test/test_events.jl`
- Modify: `packages/HimalayaUI/test/test_sse.jl`

- [ ] **Step 1: Write failing tests**

Add to `packages/HimalayaUI/test/test_actions.jl`:

```julia
@testset "get_client_op_id extracts X-Client-Op-Id header" begin
    req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "uuid-789"], UInt8[])
    @test HimalayaUI.get_client_op_id(req) == "uuid-789"
end

@testset "get_client_op_id returns nothing when absent" begin
    req = HTTP.Request("POST", "/", Pair{String,String}[], UInt8[])
    @test HimalayaUI.get_client_op_id(req) === nothing
end
```

Add to `packages/HimalayaUI/test/test_events.jl`:

```julia
@testset "apply_event! persists client_op_id from request header" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        # ... setup exposure ...
        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "uuid-abc"], UInt8[])
        result = HimalayaUI.apply_event!(db, req;
            kind="noop_test", entity_type="exposure", entity_id=exp_id, payload=Dict(:q => 1.0))
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT client_op_id FROM user_actions WHERE id = ?", [result.event_id]))
        @test String(rows[1].client_op_id) == "uuid-abc"
    end
end

@testset "apply_event! tolerates missing X-Client-Op-Id" begin
    # ... same setup, no header → client_op_id is missing/NULL ...
    @test ismissing(rows[1].client_op_id)
end

@testset "rebuild_views_from_log! tolerates NULL client_op_id" begin
    # rebuild against a fixture exposure with mixed NULL/non-NULL client_op_id
    # → rebuilt view state matches live state
end
```

Add to `packages/HimalayaUI/test/test_sse.jl` (extend existing 4 broadcast_event! testsets):

```julia
# For each existing testset that calls broadcast_event!, capture the frame
# and assert client_op_id and ts are present in JSON.
@test haskey(frame_json, "client_op_id")
@test haskey(frame_json, "ts")
# When client_op_id is non-null:
@test frame_json["client_op_id"] == "uuid-789"
```

Run tests — expect failures.

- [ ] **Step 2: Implement**

In `actions.jl`, add `get_client_op_id` mirroring `get_client_id`:

```julia
function get_client_op_id(req::HTTP.Request)::Union{String, Nothing}
    for (name, value) in HTTP.headers(req)
        if lowercase(name) == "x-client-op-id"
            return String(value)
        end
    end
    return nothing
end
```

In `events.jl`, extend `apply_event!`:

```julia
function apply_event!(db::SQLite.DB, req;
                      kind::String,
                      entity_type::String,
                      entity_id::Integer,
                      payload = nothing,
                      undoes_event_id::Union{Int,Nothing} = nothing)
    username     = get_username(req)
    client_id    = get_client_id(req)
    client_op_id = get_client_op_id(req)  # NEW
    user_id      = username === nothing ? nothing : get_or_create_user!(db, username)
    payload_json = payload === nothing ? nothing : JSON3.write(payload)

    # ... transaction ...
    DBInterface.execute(db, """
        INSERT INTO user_actions
            (user_id, action, entity_type, entity_id, payload, undoes_event_id, client_id, client_op_id)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
        [user_id, kind, entity_type, Int(entity_id), payload_json, undoes_event_id, client_id, client_op_id])
    # ...
end
```

In `broadcast_event!`, add `client_op_id` and `ts` to the JSON dict:

```julia
function broadcast_event!(event_id, kind, entity_type, entity_id, user_id, client_id, client_op_id, payload_json)
    actor = user_id === nothing ? nothing : lookup_username(current_db(), user_id)
    msg = JSON3.write(Dict(
        :id           => Int(event_id),
        :kind         => kind,
        :entity_type  => entity_type,
        :entity_id    => Int(entity_id),
        :actor        => actor,
        :client_id    => client_id,
        :client_op_id => client_op_id,                          # NEW
        :ts           => string(now(UTC)),                       # NEW
        :payload      => payload_json === nothing ? nothing : JSON3.read(payload_json),
    ))
    # ... rest unchanged ...
end
```

Update the `apply_event!` callsite to pass `client_op_id` to `broadcast_event!`.

- [ ] **Step 3: Verify**

Tests pass. Existing tests stay green. `rebuild_views_from_log!` tolerates NULL `client_op_id` — verified by extending the existing rebuild round-trip test with a fixture that has mixed NULL/non-NULL `client_op_id` values.

---

## Task M0.3: `with_idempotency` wrapper + `OP_LOCKS` registry

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl`
- New: `packages/HimalayaUI/test/test_idempotency.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl` — `include("test_idempotency.jl")`

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/test/test_idempotency.jl`:

```julia
using Test, HTTP, SQLite, DBInterface, Tables
using HimalayaUI: with_idempotency, open_db

@testset "with_idempotency: passthrough when X-Client-Op-Id absent" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        req = HTTP.Request("POST", "/", Pair{String,String}[], UInt8[])
        called = Ref(0)
        result = with_idempotency(db, req) do
            called[] += 1
            HTTP.Response(200; body = "{\"x\":1}")
        end
        @test called[] == 1
        @test result.status == 200
        # Cache table should be empty.
        rows = Tables.rowtable(DBInterface.execute(db, "SELECT * FROM idempotent_responses"))
        @test isempty(rows)
    end
end

@testset "with_idempotency: cache hit returns prior response" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "uuid-1"], UInt8[])
        called = Ref(0)
        # First call.
        r1 = with_idempotency(db, req) do
            called[] += 1
            HTTP.Response(201; body = "{\"id\":42}")
        end
        # Second call with same op_id — body should NOT execute.
        r2 = with_idempotency(db, req) do
            called[] += 1
            HTTP.Response(500; body = "{\"x\":\"should not appear\"}")
        end
        @test called[] == 1
        @test r1.status == 201
        @test r2.status == 201
        @test String(r2.body) == "{\"id\":42}"
    end
end

@testset "with_idempotency: failures NOT cached" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "uuid-2"], UInt8[])
        called = Ref(0)
        # First call returns 4xx.
        r1 = with_idempotency(db, req) do
            called[] += 1
            HTTP.Response(400; body = "{\"error\":\"bad\"}")
        end
        # Second call should re-execute body (failures aren't cached).
        r2 = with_idempotency(db, req) do
            called[] += 1
            HTTP.Response(200; body = "{\"id\":99}")
        end
        @test called[] == 2
        @test r2.status == 200
        @test String(r2.body) == "{\"id\":99}"
        # Cache now contains the successful retry.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM idempotent_responses WHERE client_op_id = 'uuid-2'"))
        @test length(rows) == 1
    end
end

@testset "with_idempotency: concurrent retries serialize via OP_LOCKS" begin
    # Two parallel tasks with same op_id, both miss the cache initially.
    # Per-op-id ReentrantLock serializes; one writes, other reads cache.
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "uuid-3"], UInt8[])
        called = Ref(0)
        body_lock = ReentrantLock()
        function delay_body()
            lock(body_lock) do
                called[] += 1
            end
            sleep(0.05)
            HTTP.Response(201; body = "{\"called\":$(called[])}")
        end
        # Spawn two tasks racing.
        t1 = @async with_idempotency(db, req) do; delay_body(); end
        t2 = @async with_idempotency(db, req) do; delay_body(); end
        r1, r2 = fetch(t1), fetch(t2)
        # Body executed exactly once.
        @test called[] == 1
        # Both responses have the same body.
        @test String(r1.body) == String(r2.body)
    end
end

@testset "with_idempotency: multi-event route cache hit returns identical body" begin
    # Simulate a route that emits N events (speculative POST shape).
    # First call commits 3 events with same client_op_id; second call returns cached.
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "uuid-4"], UInt8[])
        called = Ref(0)
        function multi_event_body()
            called[] += 1
            # Emit 3 apply_event! calls with the same client_op_id (extracted from req).
            for i in 1:3
                HimalayaUI.apply_event!(db, req;
                    kind = "speculative_created",
                    entity_type = "exposure",
                    entity_id = 1,  # assume seeded
                    payload = Dict(:n => i))
            end
            HTTP.Response(201; body = "{\"created\":3}")
        end
        r1 = with_idempotency(db, req) do; multi_event_body(); end
        r2 = with_idempotency(db, req) do; multi_event_body(); end
        @test called[] == 1
        # 3 user_actions rows from first call only.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM user_actions WHERE client_op_id = 'uuid-4'"))
        @test length(rows) == 3
        # Cached body matches.
        @test String(r1.body) == String(r2.body)
    end
end
```

Add `include("test_idempotency.jl")` to `runtests.jl`. Run — expect failures.

- [ ] **Step 2: Implement**

In `events.jl`, add the `OP_LOCKS` registry and `with_idempotency`:

```julia
const OP_LOCKS    = Dict{String, ReentrantLock}()
const OP_LOCKS_MU = ReentrantLock()

function _op_lock(op_id::String)::ReentrantLock
    lock(OP_LOCKS_MU) do
        get!(OP_LOCKS, op_id, ReentrantLock())
    end
end

function _lookup_cached_response(db::SQLite.DB, op_id::String)::Union{HTTP.Response, Nothing}
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT status_code, body FROM idempotent_responses WHERE client_op_id = ?",
        [op_id]))
    isempty(rows) && return nothing
    row = rows[1]
    return HTTP.Response(Int(row.status_code),
                         ["Content-Type" => "application/json"];
                         body = String(row.body))
end

function with_idempotency(f, db::SQLite.DB, req::HTTP.Request)
    op_id = get_client_op_id(req)
    op_id === nothing && return f()

    # Fast path: lock-free cache check.
    cached = _lookup_cached_response(db, op_id)
    cached === nothing || return cached

    # Acquire per-op-id lock; re-check cache inside the lock.
    return lock(_op_lock(op_id)) do
        cached2 = _lookup_cached_response(db, op_id)
        cached2 === nothing || return cached2

        response = f()
        if response.status < 400
            DBInterface.execute(db,
                "INSERT INTO idempotent_responses (client_op_id, status_code, body) VALUES (?, ?, ?)",
                [op_id, Int(response.status), String(response.body)])
        end
        return response
    end
end
```

- [ ] **Step 3: Verify**

All five testsets pass. The concurrent-retries test specifically validates the per-op-id lock serializes body execution.

> **Note on `OP_LOCKS` lifecycle:** entries accumulate per unique `client_op_id`. They're tiny (a `ReentrantLock` is bytes), but a long-running process eventually accrues many. Add a TTL sweep alongside the `idempotent_responses` GC sweep in Task M0.7. The two cleanup loops can share a single periodic timer.

---

## Task M0.4: `analyze_run` broadcast suppression

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl`
- Modify: `packages/HimalayaUI/test/test_sse.jl`

- [ ] **Step 1: Write failing tests**

Add to `packages/HimalayaUI/test/test_sse.jl`:

```julia
@testset "analyze_run with both skip flags true does NOT broadcast" begin
    # Capture broadcasts via the test SSE harness.
    captured = String[]
    HimalayaUI._test_subscribe!(captured)

    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        # ... seed exposure ...
        req = HimalayaUI._system_request()
        HimalayaUI.apply_event!(db, req;
            kind = "analyze_run",
            entity_type = "exposure",
            entity_id = exp_id,
            payload = Dict(
                :findpeaks_skipped => true,
                :indexpeaks_skipped => true,
                :duration_ms => 0))
        # No frame should have been broadcast.
        @test isempty(captured)
    end
end

@testset "analyze_run with one skip flag true DOES broadcast" begin
    # Only the both-skipped case is suppressed.
    captured = String[]
    HimalayaUI._test_subscribe!(captured)

    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        # ... seed exposure ...
        req = HimalayaUI._system_request()
        HimalayaUI.apply_event!(db, req;
            kind = "analyze_run",
            entity_type = "exposure",
            entity_id = exp_id,
            payload = Dict(
                :findpeaks_skipped => true,
                :indexpeaks_skipped => false))
        @test length(captured) == 1
    end
end

@testset "non-analyze_run events broadcast normally regardless of payload" begin
    # peak_added with arbitrary payload skip flags must broadcast.
    # ...
end
```

- [ ] **Step 2: Implement**

In `events.jl`, modify the broadcast call site in `apply_event!`:

```julia
if isdefined(@__MODULE__, :broadcast_event!)
    suppress = kind == "analyze_run" &&
               payload !== nothing &&
               get(payload, :findpeaks_skipped, false) === true &&
               get(payload, :indexpeaks_skipped, false) === true
    if !suppress
        try
            broadcast_event!(event_id, kind, entity_type, Int(entity_id),
                             user_id, client_id, client_op_id, payload_json)
        catch err
            @warn "broadcast_event! failed (event still durable in user_actions)" exception=err
        end
    end
end
```

> **Subtle point:** `payload` is the original Julia Dict passed to `apply_event!`, NOT the round-tripped JSON3.Object. Both forms support `get(p, :key, default)`, but the `=== true` comparison guards against the false-y JSON3 case where a missing key returns `nothing`.

- [ ] **Step 3: Verify**

All three testsets pass.

---

## Task M0.5: `analyze_exposure!` fast-skip refactor

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl`
- Modify: `packages/HimalayaUI/src/hash.jl` — new `hash_peak_set_from_db`
- New: `packages/HimalayaUI/test/test_fast_skip.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl`

This is the load-bearing M0 piece. Failure here means the synchronous-reanalyze pattern in M2 relocates latency rather than eliminating it.

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/test/test_fast_skip.jl`:

```julia
using Test, SQLite, HimalayaUI

@testset "analyze_exposure! fast-skip: zero file I/O on no-change with trace_known_unchanged=true" begin
    # Setup: exposure with auto_peaks already populated, no add curations,
    # stored hashes match what the trace + curations would produce.
    mktempdir() do tmp
        db, exp_id, analysis_dir = setup_clean_analyzed_exposure(tmp)
        # Track file reads.
        n_reads_before = count_dat_file_opens()
        HimalayaUI.analyze_exposure!(db, exp_id, analysis_dir; trace_known_unchanged=true)
        n_reads_after = count_dat_file_opens()
        @test n_reads_after == n_reads_before  # Zero file I/O.
    end
end

@testset "analyze_exposure! fast-skip: latency target <100µs P99" begin
    # Run 100 invocations on a stable exposure.
    mktempdir() do tmp
        db, exp_id, analysis_dir = setup_clean_analyzed_exposure(tmp)
        ts = Float64[]
        for _ in 1:100
            t = @elapsed HimalayaUI.analyze_exposure!(db, exp_id, analysis_dir; trace_known_unchanged=true)
            push!(ts, t)
        end
        sort!(ts)
        p99 = ts[99]
        @test p99 < 100e-6  # 100 microseconds
    end
end

@testset "analyze_exposure! fast-skip: load_dat skipped when no add-curations" begin
    # Add an exclude curation, run analyze with trace_known_unchanged=true,
    # verify no load_dat call.
end

@testset "analyze_exposure! fast-skip: load_dat REQUIRED when add-curation present" begin
    # Add a 'add'-curation, run analyze, verify load_dat IS called
    # (because we need sharpness from the trace for the new q).
end

@testset "analyze_exposure! fast-skip: trace_known_unchanged=false preserves slow path" begin
    # Default behavior unchanged: hash_trace_file fires.
    mktempdir() do tmp
        db, exp_id, analysis_dir = setup_clean_analyzed_exposure(tmp)
        n_reads_before = count_dat_file_opens()
        HimalayaUI.analyze_exposure!(db, exp_id, analysis_dir)  # no kw arg
        n_reads_after = count_dat_file_opens()
        @test n_reads_after > n_reads_before  # Slow path read the file.
    end
end

@testset "hash_peak_set_from_db produces same hash as hash_peak_set on equivalent inputs" begin
    # Three configurations exercised:
    # (a) auto peaks present, no curations
    # (b) auto peaks present, one exclude curation
    # (c) auto peaks present, multiple exclude curations
    # Each: seed the DB, then compute both ways, assert equality.
    mktempdir() do tmp
        db, exp_id, analysis_dir = setup_clean_analyzed_exposure(tmp)
        # Seed scenario (a)
        h_db = HimalayaUI.hash_peak_set_from_db(db, exp_id)
        q, I, σ = HimalayaUI.load_dat(resolve_trace_path(db, exp_id, analysis_dir))
        eff = HimalayaUI.effective_peaks(db, exp_id, q, I)
        h_ref = HimalayaUI.hash_peak_set(eff)
        @test h_db == h_ref

        # (b) Add an exclude curation; recompute.
        peak_q = first_auto_peak_q(db, exp_id)
        DBInterface.execute(db,
            "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', ?)",
            [exp_id, peak_q])
        h_db2 = HimalayaUI.hash_peak_set_from_db(db, exp_id)
        eff2 = HimalayaUI.effective_peaks(db, exp_id, q, I)
        h_ref2 = HimalayaUI.hash_peak_set(eff2)
        @test h_db2 == h_ref2
        @test h_db2 != h_db  # exclusion changed the hash

        # (c) Add a second exclude; recompute.
        peak_q2 = second_auto_peak_q(db, exp_id)
        DBInterface.execute(db,
            "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', ?)",
            [exp_id, peak_q2])
        h_db3 = HimalayaUI.hash_peak_set_from_db(db, exp_id)
        eff3 = HimalayaUI.effective_peaks(db, exp_id, q, I)
        h_ref3 = HimalayaUI.hash_peak_set(eff3)
        @test h_db3 == h_ref3
    end
end
```

Helper `count_dat_file_opens()` can be implemented by hooking `Base.open` for `.dat` files within the test scope, or by counting via a wrapper used in `load_dat`. (The exact mechanism can be a test-only helper; what matters is the assertion.)

- [ ] **Step 2: Implement helpers**

In `pipeline.jl` (or new helpers section):

```julia
function read_trace_hash(db::SQLite.DB, exposure_id::Int)::Union{String, Nothing}
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT trace_hash FROM exposures WHERE id = ?", [exposure_id]))
    isempty(rows) && return nothing
    ismissing(rows[1].trace_hash) ? nothing : String(rows[1].trace_hash)
end

read_inputs_hash(db, exposure_id) = # ... mirror shape for analysis_inputs_hash ...

function count_auto_peaks(db::SQLite.DB, exposure_id::Int)::Int
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM auto_peaks WHERE exposure_id = ?", [exposure_id]))
    Int(rows[1].c)
end

count_indices(db, exposure_id) = # ... mirror for indices ...

function any_add_curations(db::SQLite.DB, exposure_id::Int)::Bool
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM peak_curations WHERE exposure_id = ? AND kind = 'add' LIMIT 1",
        [exposure_id]))
    !isempty(rows)
end
```

In `hash.jl`:

```julia
"""
    hash_peak_set_from_db(db, exposure_id) -> String

Compute the analysis_inputs_hash from auto_peaks + exclude-curations alone,
without loading the trace. Equivalent to `hash_peak_set(effective_peaks(...))`
when no `add` curations exist for the exposure. Caller MUST verify
`!any_add_curations(db, exposure_id)` before relying on this result.
"""
function hash_peak_set_from_db(db::SQLite.DB, exposure_id::Int)::String
    # SELECT auto_peaks (q, sharpness) JOIN peak_curations(kind='exclude') for this exposure
    # → filter out excluded; sort by q; compute hash over (q::Float64, sharpness::Float64) tuples.
    # Mirror the canonical encoding used by hash_peak_set in pipeline.jl.
end
```

- [ ] **Step 3: Implement fast-skip in `analyze_exposure!`**

Restructure `analyze_exposure!` per the spec sketch:

```julia
function analyze_exposure!(db::SQLite.DB, exposure_id::Int, analysis_dir::String;
                            trace_known_unchanged::Bool=false)
    dat_path           = resolve_trace_path(db, exposure_id, analysis_dir)
    stored_trace_hash  = read_trace_hash(db, exposure_id)
    stored_inputs_hash = read_inputs_hash(db, exposure_id)
    autopeaks_count    = count_auto_peaks(db, exposure_id)
    indices_count      = count_indices(db, exposure_id)

    if trace_known_unchanged
        new_trace_hash = stored_trace_hash
        findpeaks_skipped = autopeaks_count > 0
    else
        new_trace_hash = hash_trace_file(dat_path)
        findpeaks_skipped = (stored_trace_hash == new_trace_hash) && (autopeaks_count > 0)
    end

    if findpeaks_skipped && !any_add_curations(db, exposure_id)
        new_inputs_hash = hash_peak_set_from_db(db, exposure_id)
        indexpeaks_skipped = (stored_inputs_hash == new_inputs_hash) && (indices_count > 0)
        if indexpeaks_skipped
            # Both skipped, no disk I/O.
            apply_event!(db, _system_request();
                kind        = "analyze_run",
                entity_type = "exposure",
                entity_id   = exposure_id,
                payload     = Dict(
                    :findpeaks_skipped  => true,
                    :indexpeaks_skipped => true,
                    :trace_hash         => new_trace_hash,
                    :inputs_hash        => new_inputs_hash,
                    :duration_ms        => 0,
                    :post_state_size_bytes => 0,
                ))
            return
        end
    end

    # Slow path: existing code.
    q, I, σ = load_dat(dat_path)
    # ... existing implementation ...
end
```

The slow path should also emit an `analyze_run` event with `post_state_size_bytes` populated (compute from JSON-serialized indices array length); see Task M0.6.

- [ ] **Step 4: Verify**

All test cases pass. Existing pipeline tests stay green. The latency test (P99 < 100µs) is the load-bearing assertion — if it fails, fall back per the spec's fallback trigger.

---

## Task M0.6: `post_state` payload, broadcast deferral, size observability

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl` — emit `post_state_size_bytes` on analyze_run
- Modify: `packages/HimalayaUI/src/events.jl` — `apply_event!` gains `defer_broadcast::Bool=false` and `post_state::Union{Dict,Nothing}=nothing`; `broadcast_event!` accepts and emits `post_state`
- Modify: `packages/HimalayaUI/test/test_sse.jl`
- Modify: `packages/HimalayaUI/test/test_events.jl`

This task is load-bearing for M2's replay-without-refetch correctness AND addresses a contract change to `apply_event!` that several routes depend on. Two coupled additions:

1. **`post_state` field on the SSE frame**, optionally populated by callers that know the post-analyze state.
2. **`defer_broadcast::Bool=false` keyword** on `apply_event!`. When true, the broadcast is *not* fired automatically after commit; the caller is responsible for calling `broadcast_event!` later (typically after running additional logic like `analyze_exposure!`). When false (default), behavior is identical to the pre-M0.6 contract — every existing call site is preserved without modification.

Why both in one task: the `post_state` data isn't known until *after* `analyze_exposure!` runs, but `analyze_exposure!` runs after `apply_event!`'s transaction commits. Without `defer_broadcast`, the route handler would have to either (a) issue a second broadcast (two frames per logical event — confusing for clients) or (b) allow the inner broadcast to fire without `post_state` and accept refetch-fallback. `defer_broadcast` lets the route serialize as: open transaction → apply_event! (no broadcast) → analyze_exposure! → commit → emit one broadcast with post_state. One frame per event, post_state inline.

- [ ] **Step 1: Write failing tests**

```julia
# test_events.jl
@testset "apply_event! with defer_broadcast=true does NOT fire broadcast" begin
    captured = String[]
    HimalayaUI._test_subscribe!(captured)
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        # ... seed exposure ...
        result = HimalayaUI.apply_event!(db, _test_request();
            kind="noop_test", entity_type="exposure", entity_id=exp_id,
            payload=Dict(:q => 1.0),
            defer_broadcast=true)
        @test isempty(captured)
        @test result.event_id > 0  # event still durable in user_actions
    end
end

@testset "apply_event! defer_broadcast=false (default) preserves existing behavior" begin
    # Every existing call site uses the default; broadcasts fire as before.
    captured = String[]
    HimalayaUI._test_subscribe!(captured)
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        # ... seed exposure ...
        HimalayaUI.apply_event!(db, _test_request();
            kind="noop_test", entity_type="exposure", entity_id=exp_id,
            payload=Dict(:q => 1.0))  # no defer_broadcast kwarg
        @test length(captured) == 1
    end
end

@testset "broadcast_event! emits post_state when provided" begin
    captured = String[]
    HimalayaUI._test_subscribe!(captured)
    HimalayaUI.broadcast_event!(1, "peak_added", "exposure", 42,
        nothing, "tab-id", "op-id", JSON3.write(Dict(:q => 1.0));
        post_state = Dict(:analysis_inputs_hash => "abc123", :indices => []))
    @test length(captured) == 1
    frame_json = JSON3.read(captured[1])
    @test haskey(frame_json, :post_state)
    @test frame_json.post_state.analysis_inputs_hash == "abc123"
end

@testset "broadcast_event! omits post_state when not provided" begin
    captured = String[]
    HimalayaUI._test_subscribe!(captured)
    HimalayaUI.broadcast_event!(1, "peak_added", "exposure", 42,
        nothing, "tab-id", "op-id", JSON3.write(Dict(:q => 1.0)))
    frame_json = JSON3.read(captured[1])
    @test !haskey(frame_json, :post_state)
end

# test_pipeline.jl
@testset "analyze_run slow path records post_state_size_bytes" begin
    mktempdir() do tmp
        db, exp_id, analysis_dir = setup_clean_analyzed_exposure(tmp)
        # Force slow path: change the trace.
        # ... mutate trace file ...
        HimalayaUI.analyze_exposure!(db, exp_id, analysis_dir)
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT payload FROM user_actions WHERE entity_id = ? AND action = 'analyze_run' ORDER BY id DESC LIMIT 1",
            [exp_id]))
        payload = JSON3.read(String(rows[1].payload))
        @test get(payload, :post_state_size_bytes, 0) > 0
    end
end
```

- [ ] **Step 2: Extend `broadcast_event!` and `apply_event!`**

```julia
function broadcast_event!(event_id, kind, entity_type, entity_id,
                          user_id, client_id, client_op_id, payload_json;
                          post_state::Union{Dict, Nothing} = nothing)
    actor = user_id === nothing ? nothing : lookup_username(current_db(), user_id)
    fields = Dict(
        :id           => Int(event_id),
        :kind         => kind,
        :entity_type  => entity_type,
        :entity_id    => Int(entity_id),
        :actor        => actor,
        :client_id    => client_id,
        :client_op_id => client_op_id,
        :ts           => string(now(UTC)),
        :payload      => payload_json === nothing ? nothing : JSON3.read(payload_json),
    )
    post_state === nothing || (fields[:post_state] = post_state)
    msg = JSON3.write(fields)
    frame = "event: curation\ndata: $msg\n\n"
    lock(SSE_LOCK) do
        # ... existing pruning logic unchanged ...
    end
    nothing
end
```

```julia
function apply_event!(db::SQLite.DB, req;
                      kind::String,
                      entity_type::String,
                      entity_id::Integer,
                      payload = nothing,
                      undoes_event_id::Union{Int,Nothing} = nothing,
                      defer_broadcast::Bool = false,
                      post_state::Union{Dict, Nothing} = nothing)
    # ... existing extraction + transaction logic ...

    if !defer_broadcast && isdefined(@__MODULE__, :broadcast_event!)
        # Apply the both-skipped suppression rule from M0.4.
        suppress = kind == "analyze_run" &&
                   payload !== nothing &&
                   get(payload, :findpeaks_skipped, false) === true &&
                   get(payload, :indexpeaks_skipped, false) === true
        if !suppress
            try
                broadcast_event!(event_id, kind, entity_type, Int(entity_id),
                                 user_id, client_id, client_op_id, payload_json;
                                 post_state = post_state)
            catch err
                @warn "broadcast_event! failed (event still durable in user_actions)" exception=err
            end
        end
    end

    return (event_id = event_id, view_row_id = view_row_id_ref[])
end
```

- [ ] **Step 3: Update `analyze_exposure!` slow path to record `post_state_size_bytes`**

```julia
# In the slow-path branch of analyze_exposure!:
indices_array = serialize_indices_for_post_state(db, exposure_id)
size_bytes = length(JSON3.write(indices_array))
apply_event!(db, _system_request();
    kind = "analyze_run",
    entity_type = "exposure",
    entity_id = exposure_id,
    payload = Dict(
        :findpeaks_skipped     => false,
        :indexpeaks_skipped    => false,
        :duration_ms           => duration_ms,
        :post_state_size_bytes => size_bytes,
    ))
```

The slow path does NOT pass `post_state` as a kwarg — `analyze_exposure!` is called from many places (`reingest`, scheduled scans, etc.); it's the route handler's job to attach the full `post_state` to the curation event when desired (see M2.2). `analyze_exposure!` only records the *size* for observability.

`serialize_indices_for_post_state(db, exposure_id)` is a small helper that JSON-shapes the current indices for an exposure — same shape `GET /api/exposures/:id/indices` returns, reusing the existing serialization in `json.jl` where possible.

- [ ] **Step 4: Verify**

All five testsets pass. Critically, the "default preserves existing behavior" test asserts that no existing callsite of `apply_event!` regresses — backward compatibility for the unchanged routes.

> **Why this lives in M0, not M2.2:** the `defer_broadcast` semantics is a contract change to a function that *every* route depends on. Putting it in M2.2 would mean the contract change ships in the same PR as three frontend mutator migrations + an autoReanalyze deletion + a banner update + new components. Pulling it into M0 keeps the contract change small, well-tested, and reviewable independently. M2.2 then just *uses* the new kwarg.

---

## Task M0.7: GC sweep for `idempotent_responses` and `OP_LOCKS`

**Files:**
- Modify: `packages/HimalayaUI/src/server.jl` or new `packages/HimalayaUI/src/gc.jl`
- Modify: `packages/HimalayaUI/test/test_idempotency.jl`

- [ ] **Step 1: Write failing tests**

```julia
@testset "idempotent_responses GC sweeps rows older than TTL" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        DBInterface.execute(db, """
            INSERT INTO idempotent_responses (client_op_id, status_code, body, created_at)
            VALUES ('old-1', 200, '{}', datetime('now', '-2 hours'))
        """)
        DBInterface.execute(db, """
            INSERT INTO idempotent_responses (client_op_id, status_code, body, created_at)
            VALUES ('new-1', 200, '{}', datetime('now', '-5 minutes'))
        """)
        HimalayaUI.gc_idempotent_responses!(db; ttl_seconds = 3600)
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT client_op_id FROM idempotent_responses ORDER BY client_op_id"))
        @test [String(r.client_op_id) for r in rows] == ["new-1"]
    end
end

@testset "OP_LOCKS GC removes entries whose responses have been swept" begin
    # Assert OP_LOCKS shrinks after gc.
end
```

- [ ] **Step 2: Implement**

```julia
function gc_idempotent_responses!(db::SQLite.DB; ttl_seconds::Int = 3600)
    DBInterface.execute(db,
        "DELETE FROM idempotent_responses WHERE created_at < datetime('now', ?)",
        ["-$(ttl_seconds) seconds"])
    # Also sweep OP_LOCKS: any op_id whose response was just deleted (or was
    # never recorded) is unlikely to retry. Cheap to just clear entries
    # whose corresponding response is absent.
    lock(OP_LOCKS_MU) do
        live = Set{String}()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT client_op_id FROM idempotent_responses"))
        for r in rows
            push!(live, String(r.client_op_id))
        end
        for k in keys(OP_LOCKS)
            k in live || delete!(OP_LOCKS, k)
        end
    end
end
```

Schedule the sweep on a Timer (every 30 minutes is generous):

```julia
# server.jl
const GC_TIMER = Ref{Union{Timer, Nothing}}(nothing)

function start_gc_timer!(db)
    GC_TIMER[] = Timer(0.0; interval = 30 * 60) do _
        try
            gc_idempotent_responses!(db; ttl_seconds = 3600)
        catch err
            @warn "idempotent_responses GC sweep failed" exception = err
        end
    end
end
```

Call `start_gc_timer!` from `serve(db)`.

- [ ] **Step 3: Verify**

Both testsets pass.

---

## Task M0.8: Frontend `clientOpId` minting + header threading

**Files:**
- New: `packages/HimalayaUI/frontend/src/lib/clientOpId.ts`
- Modify: `packages/HimalayaUI/frontend/src/api.ts` — `AuthOpts.clientOpId`, header threading
- Modify: `packages/HimalayaUI/frontend/src/queries.ts` — `authOpts` extension; mutators mint per-call
- New: `packages/HimalayaUI/frontend/test/clientOpId.test.ts`
- Modify: `packages/HimalayaUI/frontend/test/api.test.ts`

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/clientOpId.test.ts`:

```typescript
import { describe, it, expect } from "vitest";
import { newClientOpId } from "../src/lib/clientOpId";

describe("newClientOpId", () => {
  it("mints a UUID", () => {
    const id = newClientOpId();
    expect(id).toMatch(/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/);
  });
  it("returns unique IDs across calls", () => {
    const ids = new Set<string>();
    for (let i = 0; i < 100; i++) ids.add(newClientOpId());
    expect(ids.size).toBe(100);
  });
});
```

Add to `packages/HimalayaUI/frontend/test/api.test.ts` (parallel to existing `X-Username` tests):

```typescript
it("threads X-Client-Op-Id on POST mutations", async () => {
  // ...mock fetch, call api.addPeak with clientOpId, assert header...
});
it("omits X-Client-Op-Id on GET requests", async () => {
  // ...
});
it("threads X-Client-Op-Id without X-Username (system flows)", async () => {
  // ...
});
```

- [ ] **Step 2: Implement**

Create `lib/clientOpId.ts`:

```typescript
export function newClientOpId(): string {
  return crypto.randomUUID();
}
```

Extend `AuthOpts` in `api.ts`:

```typescript
export interface AuthOpts {
  username?: string;
  clientId?: string;
  clientOpId?: string;
}
```

Update `request()` to thread `X-Client-Op-Id` on mutations:

```typescript
if (method !== "GET" && opts?.clientOpId !== undefined) {
  headers.set("X-Client-Op-Id", opts.clientOpId);
}
```

In `queries.ts`, extend `authOpts`:

```typescript
function authOpts(
  username: string | undefined,
  clientId: string,
  clientOpId?: string,
): api.AuthOpts {
  const out: api.AuthOpts = { clientId };
  if (username !== undefined) out.username = username;
  if (clientOpId !== undefined) out.clientOpId = clientOpId;
  return out;
}
```

For mutations that do NOT yet go through the queue (M0 still has all 16 using legacy `useMutation`), mint a fresh `clientOpId` per call:

```typescript
return useMutation({
  mutationFn: async (q: number) => {
    const opts = authOpts(username, CLIENT_ID, newClientOpId());
    const peak = await api.addPeak(exposureId, q, opts);
    await autoReanalyze(exposureId, username, CLIENT_ID);  // unchanged in M0
    return peak;
  },
  onSuccess: () => invalidateExposure(qc, exposureId),
});
```

> **Why mint per-call rather than per-mutation-instance:** TanStack Query's `useMutation` reuses the same `mutationFn` reference across multiple `mutate()` invocations. A `clientOpId` minted at hook-creation time would be reused for every click, defeating retry idempotency. Mint inside `mutationFn` so each call gets its own.

- [ ] **Step 3: Verify**

All tests pass. Run `npm test` for the frontend suite — no regressions.

---

## M0 verification

After all M0 tasks land, manual verification:

- [ ] Backend test suite: `Pkg.test("HimalayaUI")` — all green, including new test_idempotency.jl, test_fast_skip.jl, extended test_db.jl/test_actions.jl/test_events.jl/test_sse.jl.
- [ ] Frontend test suite: `npm test` — all green, including new clientOpId.test.ts and extended api.test.ts.
- [ ] Frontend build: `npm run build` — TS strict + Vite clean.
- [ ] Live wire smoke: spin up the test server, fire a curation request with `X-Client-Op-Id`, retry the same request, observe the second response is byte-identical to the first and the body wasn't re-executed (verifiable via test logs or DB inspection).
- [ ] Latency check: run the M0.5 P99 latency test on real fixture data; confirm <100µs.

If the latency check fails, **STOP** and apply the corresponding fallback trigger from the spec — revert M2's synchronous-reanalyze pattern before proceeding.

> **M0 is *user-invisible* but not *wire-invisible*.** Two observable changes ship in M0 that are worth flagging:
>
> 1. **`analyze_run` broadcast suppression.** Previously every analysis call (including idempotent re-runs from the autoReanalyze chain) broadcast a frame to all subscribers. After M0.4, `analyze_run` events with both `findpeaks_skipped` and `indexpeaks_skipped` true are suppressed. Pre-M0 clients that depended on the spurious frame for some side effect would notice; HimalayaUI's frontend ignores them, so this is desirable.
>
> 2. **SSE frame shape.** M0.2 adds `client_op_id` and `ts`; M0.6 adds optional `post_state`. Clients that strictly typed-validate frame shape (rather than ignoring unknown fields) would error. The current `sseSubscriber.ts` uses indexed property access; new fields are tolerated.
>
> Neither is a regression for the intended client (HimalayaUI's own frontend), but the "behavioral no-op" claim earlier in the plan is precisely that M0 changes *no observable user behavior*. Wire changes are explicit; flag in release notes.

---

## Task M1.1: Queue framework types and registry

**Files:**
- New: `packages/HimalayaUI/frontend/src/lib/queue/types.ts`
- New: `packages/HimalayaUI/frontend/src/lib/queue/deferred.ts`
- New: `packages/HimalayaUI/frontend/src/lib/queue/index.ts`
- New: `packages/HimalayaUI/frontend/test/queue/deferred.test.ts`

- [ ] **Step 1: Write failing tests**

```typescript
// queue/deferred.test.ts
describe("pendingDeferreds", () => {
  it("registers and resolves", async () => {
    const d = makeDeferred<{ id: number }>("op-1");
    setTimeout(() => d.resolve({ id: 42 }), 0);
    const result = await d.promise;
    expect(result).toEqual({ id: 42 });
  });

  it("registry lookup finds by client_op_id", () => {
    const d = makeDeferred<unknown>("op-2");
    expect(pendingDeferreds.get("op-2")).toBe(d);
  });

  it("delete removes from registry", async () => {
    const d = makeDeferred<unknown>("op-3");
    d.resolve(null);
    await d.promise.finally(() => pendingDeferreds.delete("op-3"));
    expect(pendingDeferreds.get("op-3")).toBeUndefined();
  });

  it("AbortSignal rejection clears the deferred", async () => {
    const ctrl = new AbortController();
    const d = makeDeferred<unknown>("op-4");
    ctrl.signal.addEventListener("abort", () => d.reject(new DOMException("aborted", "AbortError")));
    ctrl.abort();
    await expect(d.promise).rejects.toThrow("aborted");
  });
});
```

- [ ] **Step 2: Implement types and deferred registry**

```typescript
// queue/types.ts

// Optimistic-id invariant:
// Mutators that need an entity id before the server has assigned one
// (peak_added inserting a peak_curations row; speculative_created inserting
// an indices row) use NEGATIVE placeholder ids in optimistic cache writes.
// Real DB ids are always positive (INTEGER PRIMARY KEY AUTOINCREMENT in
// SQLite never returns ≤ 0 for fresh inserts). Consumers that read these
// caches (PeakRow, MentionChip, etc.) must tolerate negative ids: don't
// `id > 0` filter, don't strict-parse, don't dereference into URLs without
// a sign check. The placeholder is replaced with the real id when the
// mutator's onSuccess runs against the server response.

export type OpKind =
  | "peak_added" | "peak_excluded" | "peak_unexcluded" | "peak_removed"
  | "index_confirmed" | "index_unconfirmed"
  | "speculative_created" | "speculative_deleted"
  | "set_exposure_status" | "select_exposure"
  | "add_tag" | "remove_tag"
  | "post_message" | "update_sample"
  | "reanalyze_exposure"
  | "delete_index";

export interface OpPayload<T = unknown> {
  kind: OpKind;
  payload: T;
  clientOpId: string;
  exposureId?: number;
}

export interface RollbackContext {
  restore: () => void;
}

export interface PendingDeferred<T> {
  promise: Promise<T>;
  resolve: (value: T) => void;
  reject: (reason: unknown) => void;
}

export interface SseEvent {
  id: number;
  kind: string;
  entity_type: string;
  entity_id: number;
  actor?: string | null;
  client_id?: string | null;
  client_op_id?: string | null;
  ts?: string;
  payload?: unknown;
  post_state?: { analysis_inputs_hash: string; indices: unknown[] };
}

export interface Mutator<TPayload, TResponse> {
  kind: OpKind;
  onMutate: (payload: TPayload, qc: QueryClient) => RollbackContext;
  request: (payload: TPayload, signal: AbortSignal) => Promise<TResponse>;
  onSuccess: (payload: TPayload, response: TResponse, qc: QueryClient) => void;
  affectsExposurePeaks?: (payload: TPayload, exposureId: number) => boolean;
}
```

```typescript
// queue/deferred.ts
const pendingDeferreds = new Map<string, PendingDeferred<unknown>>();

export function makeDeferred<T>(clientOpId: string): PendingDeferred<T> {
  let resolve!: (v: T) => void;
  let reject!: (e: unknown) => void;
  const promise = new Promise<T>((res, rej) => { resolve = res; reject = rej; });
  const d: PendingDeferred<T> = { promise, resolve, reject };
  pendingDeferreds.set(clientOpId, d as PendingDeferred<unknown>);
  return d;
}

export function clearDeferred(clientOpId: string): void {
  pendingDeferreds.delete(clientOpId);
}

export function getDeferred(clientOpId: string): PendingDeferred<unknown> | undefined {
  return pendingDeferreds.get(clientOpId);
}
```

`queue/index.ts` re-exports everything for ergonomic imports.

- [ ] **Step 3: Verify**

All four testcases pass.

---

## Task M1.2: Replay coordinator

**Files:**
- New: `packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts`
- New: `packages/HimalayaUI/frontend/test/queue/replayCoordinator.test.ts`

- [ ] **Step 1: Write failing tests**

```typescript
// queue/replayCoordinator.test.ts
import { describe, it, expect, beforeEach } from "vitest";
import { QueryClient, MutationCache } from "@tanstack/react-query";
import { handleRemoteEvent } from "../src/lib/queue/replayCoordinator";
import { makeDeferred, getDeferred } from "../src/lib/queue/deferred";

describe("handleRemoteEvent", () => {
  let qc: QueryClient;
  let mc: MutationCache;

  beforeEach(() => {
    qc = new QueryClient();
    mc = qc.getMutationCache();
  });

  it("SSE with matching client_op_id resolves the pending deferred", async () => {
    const d = makeDeferred<{ ok: true }>("op-1");
    handleRemoteEvent({
      id: 1, kind: "peak_added", entity_type: "exposure", entity_id: 42,
      client_op_id: "op-1", payload: { q: 1.0 },
    }, qc, mc);
    const result = await d.promise;
    expect(result).toMatchObject({ event_id: 1 });
    expect(getDeferred("op-1")).toBeUndefined();  // deferred has been resolved
  });

  it("SSE without matching client_op_id triggers rollback-apply-replay", () => {
    // Spy onMutate / context.restore on a fake pending mutation.
    const restore = vi.fn();
    const onMutate = vi.fn();
    const fakeMutation = makeFakeMutation({ status: "pending", context: { restore }, onMutate });
    mc.add(fakeMutation);
    handleRemoteEvent({
      id: 2, kind: "peak_excluded", entity_type: "exposure", entity_id: 42,
      client_op_id: "other-tab-op", payload: { q: 2.0 },
    }, qc, mc);
    expect(restore).toHaveBeenCalledTimes(1);
    expect(onMutate).toHaveBeenCalledTimes(1);
  });

  it("Rollback iterates pending mutations in reverse order", () => {
    const order: number[] = [];
    for (let i = 0; i < 3; i++) {
      mc.add(makeFakeMutation({
        status: "pending",
        context: { restore: () => order.push(i) },
        onMutate: () => {},
      }));
    }
    handleRemoteEvent(remoteForeignEvent(), qc, mc);
    expect(order).toEqual([2, 1, 0]);
  });

  it("Re-run iterates in original (insertion) order", () => {
    const order: number[] = [];
    for (let i = 0; i < 3; i++) {
      mc.add(makeFakeMutation({
        status: "pending",
        context: { restore: () => {} },
        onMutate: () => order.push(i),
      }));
    }
    handleRemoteEvent(remoteForeignEvent(), qc, mc);
    expect(order).toEqual([0, 1, 2]);
  });

  it("MutationCache.getAll() insertion order is preserved (load-bearing TanStack invariant)", () => {
    // This test exists specifically to fail loudly if a future TanStack
    // version changes MutationCache's underlying storage from Set to
    // something with different iteration semantics. See Open Q #7.
    const ids: string[] = ["a", "b", "c", "d", "e"];
    for (const id of ids) {
      mc.add(makeFakeMutation({ status: "pending", context: {}, mutationKey: [id] }));
    }
    const observed = mc.getAll().map(m => (m.options.mutationKey as string[])[0]);
    expect(observed).toEqual(ids);
  });

  it("AbortSignal during replay does not double-rollback", () => {
    const restore = vi.fn();
    const m = makeFakeMutation({ status: "pending", context: { restore }, onMutate: () => {} });
    mc.add(m);
    handleRemoteEvent(remoteForeignEvent(), qc, mc);
    handleRemoteEvent(remoteForeignEvent(), qc, mc);
    // Two replay events but rollback context's restore should be invoked
    // once per event (each event independently rolls back, applies, replays).
    // The "double-rollback" guard means: within ONE event, restore fires at
    // most once. Verified by: restore called exactly twice across two events,
    // not four.
    expect(restore).toHaveBeenCalledTimes(2);
  });
});
```

`makeFakeMutation` is a small test helper that fabricates a `Mutation`-shaped object with the fields the coordinator reads; it lives in `test/queue/helpers.ts`. `remoteForeignEvent()` returns a canned event with a `client_op_id` not matching any live deferred.

- [ ] **Step 2: Implement**

```typescript
// queue/replayCoordinator.ts
import type { QueryClient, MutationCache } from "@tanstack/react-query";
import type { SseEvent } from "./types";
import { getDeferred, clearDeferred } from "./deferred";
import { queryKeys } from "../../queries";

export function handleRemoteEvent(
  remote: SseEvent,
  qc: QueryClient,
  mc: MutationCache,
): void {
  // SSE-confirms-my-pending-op path.
  if (remote.client_op_id) {
    const deferred = getDeferred(remote.client_op_id);
    if (deferred) {
      deferred.resolve(synthesizeResponseFromSse(remote));
      clearDeferred(remote.client_op_id);
      return;
    }
  }

  // Remote op from another tab/user.
  const pending = mc.getAll().filter(m => m.state.status === "pending");

  // Rollback in reverse queue order.
  for (const m of [...pending].reverse()) {
    const ctx = m.state.context as { restore?: () => void } | undefined;
    ctx?.restore?.();
  }

  applyRemoteToCache(remote, qc);

  // Re-run optimistic effects in queue order.
  for (const m of pending) {
    m.options.onMutate?.(m.state.variables);
  }
}

function synthesizeResponseFromSse(remote: SseEvent): unknown {
  return {
    event_id:             remote.id,
    client_op_id:         remote.client_op_id,
    analysis_inputs_hash: remote.post_state?.analysis_inputs_hash,
    // Per-route fields lifted from payload (q, exclude target, etc.)
    ...((remote.payload as object) ?? {}),
  };
}

function applyRemoteToCache(remote: SseEvent, qc: QueryClient): void {
  const id = remote.entity_id;
  const payload = remote.payload as Record<string, unknown> | undefined;

  // Apply post_state (atomic indices + exposure.analysis_inputs_hash) when
  // present — this is the replay-without-refetch path.
  const applyPostState = () => {
    if (!remote.post_state) return;
    qc.setQueryData(queryKeys.indices(id), remote.post_state.indices);
    qc.setQueryData(queryKeys.exposure(id), (old: Exposure | undefined) =>
      old ? { ...old, analysis_inputs_hash: remote.post_state!.analysis_inputs_hash } : old);
  };

  // Note on event-kind coverage:
  // Today only events emitted via `apply_event!` reach SSE — that's
  // peak_added/excluded/unexcluded/removed, index_confirmed/unconfirmed,
  // speculative_created/deleted, and analyze_run. The cases below for
  // post_message, set_exposure_status, select_exposure, add_tag, remove_tag,
  // and update_sample are forward-scaffolded: those routes still use
  // `log_action!` (no broadcast) today and won't trigger these branches
  // until the trivial-mutator slice (M2.1) migrates them to `apply_event!`
  // for cross-tab parity. Until then the branches are unreachable; the
  // `default:` invalidate fallback would also handle them correctly if
  // they fired.

  switch (remote.kind) {
    case "peak_added": {
      // Append the optimistic curation row; server-side row id is on the event
      // (the inner broadcast was deferred; this enriched broadcast carries it
      // via remote.id and the matching curation row in payload).
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) => [
        ...old,
        { id: -remote.id, q: payload?.q as number, kind: "add", excluded: false },
      ]);
      applyPostState();
      break;
    }
    case "peak_excluded": {
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) =>
        old.map(p => Math.abs(p.q - (payload?.q as number)) < 1e-6
          ? { ...p, excluded: true } : p));
      applyPostState();
      break;
    }
    case "peak_unexcluded": {
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) =>
        old.map(p => Math.abs(p.q - (payload?.q as number)) < 1e-6
          ? { ...p, excluded: false } : p));
      applyPostState();
      break;
    }
    case "peak_removed": {
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) =>
        old.filter(p => Math.abs(p.q - (payload?.q as number)) >= 1e-6));
      applyPostState();
      break;
    }
    case "index_confirmed": {
      const groupId = payload?.group_id as number;
      const indexId = payload?.index_id as number;
      qc.setQueryData<Group[]>(queryKeys.groups(id), (old = []) =>
        old.map(g => g.id === groupId
          ? { ...g, members: [...g.members, indexId] } : g));
      break;
    }
    case "index_unconfirmed": {
      const groupId = payload?.group_id as number;
      const indexId = payload?.index_id as number;
      qc.setQueryData<Group[]>(queryKeys.groups(id), (old = []) =>
        old.map(g => g.id === groupId
          ? { ...g, members: g.members.filter(m => m !== indexId) } : g));
      break;
    }
    case "speculative_created":
    case "speculative_deleted": {
      // Speculatives can change index ordering; safer to invalidate indices
      // and groups rather than reconstruct piecewise. Refetch is bounded.
      qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
      qc.invalidateQueries({ queryKey: queryKeys.groups(id) });
      break;
    }
    case "set_exposure_status": {
      qc.setQueryData(queryKeys.exposure(id), (old: Exposure | undefined) =>
        old ? { ...old, status: payload?.status as string } : old);
      break;
    }
    case "post_message": {
      // Append to message list cache.
      const sampleId = payload?.sample_id as number;
      qc.setQueryData<Message[]>(queryKeys.messages(sampleId), (old = []) => [
        ...old,
        payload as Message,
      ]);
      break;
    }
    case "add_tag":
    case "remove_tag": {
      // Tags live on samples or exposures; the entity_type tells us which.
      // Invalidate the parent collection — tag updates are infrequent.
      const parentKey = remote.entity_type === "sample"
        ? queryKeys.samples(payload?.experiment_id as number)
        : ["sample", payload?.sample_id, "exposures"];
      qc.invalidateQueries({ queryKey: parentKey });
      break;
    }
    case "update_sample": {
      qc.setQueryData(queryKeys.sample(id), (old: Sample | undefined) =>
        old ? { ...old, ...(payload ?? {}) } : old);
      break;
    }
    case "select_exposure": {
      // Sample-scoped LWW. Update all exposures' selected flag.
      const sampleId = payload?.sample_id as number;
      qc.invalidateQueries({ queryKey: queryKeys.exposures(sampleId) });
      break;
    }
    case "analyze_run": {
      // No view-table writes; only post_state matters here.
      applyPostState();
      break;
    }
    default: {
      // Unknown kind: invalidate-and-refetch as fallback so the cache
      // converges even if we can't reason about the event shape.
      qc.invalidateQueries({ queryKey: queryKeys.peaks(id) });
      qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
      qc.invalidateQueries({ queryKey: queryKeys.groups(id) });
    }
  }
}
```

> **Why some events use `setQueryData` and others `invalidateQueries`:** the per-kind logic mirrors the spec's "replay-without-refetch where possible, fall back to refetch where the event payload is insufficient." Curation events for which the route handler issued an enriched broadcast (with `post_state`) get full cache-merge. Speculative events affect index ordering globally and refetching is simpler. Tag events are infrequent enough that refetch is cheaper than per-tag merge logic. Update this map when adding a new event kind.

- [ ] **Step 3: Verify**

All six testcases pass. The MutationCache iteration-order test (case 5) is the load-bearing assertion against TanStack version drift — see Open Question #7 in the spec.

---

## Task M1.3: Persistence layer

**Files:**
- New: `packages/HimalayaUI/frontend/src/lib/queue/persistence.ts`
- New: `packages/HimalayaUI/frontend/test/queue/persistence.test.ts`

- [ ] **Step 1: Write failing tests**

```typescript
// Six cases:
// 1. Pending mutations are mirrored to sessionStorage on enqueue
// 2. Confirmed mutations are removed from sessionStorage
// 3. Rehydrate runs onMutate against fresh cache on tab reload
// 4. Rehydrate re-fires request with original clientOpId
// 5. Schema-version mismatch drops the persisted op with a toast
// 6. AbortSignal-cancelled mutations clear from sessionStorage
```

- [ ] **Step 2: Implement**

```typescript
// queue/persistence.ts
import type { MutationCache, QueryClient } from "@tanstack/react-query";
import type { OpKind } from "./types";

const STORAGE_KEY = "himalaya-ui:queue";
const SCHEMA_VERSION = 1;

interface PersistedOp {
  schemaVersion: number;
  kind: OpKind;
  payload: unknown;
  clientOpId: string;
}

export function attachPersistence(mc: MutationCache): () => void {
  const unsubscribe = mc.subscribe(event => {
    if (event.type === "added" || event.type === "updated") {
      mirrorToSessionStorage(mc);
    }
  });
  return unsubscribe;
}

function mirrorToSessionStorage(mc: MutationCache): void {
  const pending = mc.getAll().filter(m => m.state.status === "pending");
  const ops: PersistedOp[] = pending.map(m => ({
    schemaVersion: SCHEMA_VERSION,
    kind: (m.state.variables as { kind: OpKind }).kind,
    payload: m.state.variables,
    clientOpId: (m.state.variables as { clientOpId: string }).clientOpId,
  }));
  sessionStorage.setItem(STORAGE_KEY, JSON.stringify(ops));
}

export async function rehydrate(qc: QueryClient, mutators: Map<OpKind, Mutator<unknown, unknown>>): Promise<void> {
  const raw = sessionStorage.getItem(STORAGE_KEY);
  if (!raw) return;
  let ops: PersistedOp[] = [];
  try {
    ops = JSON.parse(raw);
  } catch {
    sessionStorage.removeItem(STORAGE_KEY);
    return;
  }
  let droppedCount = 0;
  for (const op of ops) {
    if (op.schemaVersion !== SCHEMA_VERSION) {
      droppedCount++;
      continue;
    }
    const mutator = mutators.get(op.kind);
    if (!mutator) {
      droppedCount++;
      continue;
    }
    // Re-fire the mutation; server-side idempotency returns cached response.
    qc.getMutationCache().build(qc, {
      mutationKey: [op.kind],
      mutationFn: (vars) => mutator.request(vars, new AbortController().signal),
      onMutate: (vars) => mutator.onMutate(vars, qc),
      onSuccess: (data, vars) => mutator.onSuccess(vars, data, qc),
    }).execute(op.payload);
  }
  if (droppedCount > 0) {
    showToast(`${droppedCount} edits from the previous session were lost`, "warning");
  }
  sessionStorage.removeItem(STORAGE_KEY);
}
```

- [ ] **Step 3: Verify**

All six testcases pass.

---

## Task M1.4: `useExposureHasPendingPeakOps` and `useQueueOpStatus` hooks

**Files:**
- New: `packages/HimalayaUI/frontend/src/lib/queue/hooks.ts`
- New: `packages/HimalayaUI/frontend/test/queue/hooks.test.ts`

- [ ] **Step 1: Write failing tests**

```typescript
// queue/hooks.test.ts
import { describe, it, expect } from "vitest";
import { renderHook, act } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useExposureHasPendingPeakOps, useQueueOpStatus } from "../src/lib/queue/hooks";

function withQueryClient(qc: QueryClient) {
  return ({ children }: { children: React.ReactNode }) =>
    <QueryClientProvider client={qc}>{children}</QueryClientProvider>;
}

describe("useExposureHasPendingPeakOps", () => {
  it("returns false when no peak ops are pending for the exposure", () => {
    const qc = new QueryClient();
    const { result } = renderHook(() => useExposureHasPendingPeakOps(42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe(false);
  });

  it("returns true while a peak_added op is pending for the same exposure", () => {
    const qc = new QueryClient();
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_added", exposureId: 42, q: 1.0, clientOpId: "op-1" },
    }));
    const { result } = renderHook(() => useExposureHasPendingPeakOps(42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe(true);
  });

  it("returns false for a different exposure id", () => {
    const qc = new QueryClient();
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_added", exposureId: 99, q: 1.0, clientOpId: "op-1" },
    }));
    const { result } = renderHook(() => useExposureHasPendingPeakOps(42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe(false);
  });

  it("returns false for a non-peak op kind on the same exposure", () => {
    const qc = new QueryClient();
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "add_tag", exposureId: 42, key: "k", value: "v", clientOpId: "op-1" },
    }));
    const { result } = renderHook(() => useExposureHasPendingPeakOps(42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe(false);
  });

  it("returns false again after the op confirms", async () => {
    const qc = new QueryClient();
    const m = makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_added", exposureId: 42, q: 1.0, clientOpId: "op-1" },
    });
    qc.getMutationCache().add(m);
    const { result, rerender } = renderHook(() => useExposureHasPendingPeakOps(42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe(true);
    act(() => { m.state.status = "success"; qc.getMutationCache().notify(m); });
    rerender();
    expect(result.current).toBe(false);
  });
});

describe("useQueueOpStatus", () => {
  it("returns 'pending' for matching kind", () => {
    const qc = new QueryClient();
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "reanalyze_exposure", exposureId: 42, clientOpId: "op-1" },
    }));
    const { result } = renderHook(() => useQueueOpStatus("reanalyze_exposure", 42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe("pending");
  });
  it("returns 'idle' for non-matching kind", () => { /* analogous */ });
  it("ignores exposureId scope when not provided", () => { /* analogous */ });
});
```

- [ ] **Step 2: Implement**

```typescript
// queue/hooks.ts
import { useMutationState } from "@tanstack/react-query";
import type { OpKind, OpPayload } from "./types";

const PEAK_AFFECTING_KINDS: Set<OpKind> = new Set([
  "peak_added", "peak_excluded", "peak_unexcluded", "peak_removed",
  "reanalyze_exposure",
]);

export function useExposureHasPendingPeakOps(exposureId: number | undefined): boolean {
  const pending = useMutationState({
    filters: {
      predicate: (m) => {
        if (m.state.status !== "pending") return false;
        const op = m.state.variables as OpPayload | undefined;
        if (!op) return false;
        if (!PEAK_AFFECTING_KINDS.has(op.kind)) return false;
        return op.exposureId === exposureId;
      },
    },
  });
  return pending.length > 0;
}

export function useQueueOpStatus(kind: OpKind, exposureId?: number): "idle" | "pending" {
  const pending = useMutationState({
    filters: {
      predicate: (m) => {
        if (m.state.status !== "pending") return false;
        const op = m.state.variables as OpPayload | undefined;
        if (!op || op.kind !== kind) return false;
        if (exposureId !== undefined && op.exposureId !== exposureId) return false;
        return true;
      },
    },
  });
  return pending.length > 0 ? "pending" : "idle";
}
```

- [ ] **Step 3: Verify**

All three testcases pass.

---

## Task M1.5: `useQueueMutation` framework wrapper

This is the integration point where the deferred-promise primitive (M1.1), the replay coordinator (M1.2), the persistence layer (M1.3), and the hooks (M1.4) come together as a single hook that consumers (M2) call once per mutator.

**Files:**
- New: `packages/HimalayaUI/frontend/src/lib/queue/useQueueMutation.ts`
- New: `packages/HimalayaUI/frontend/test/queue/useQueueMutation.test.ts`

- [ ] **Step 1: Specify the public surface**

```typescript
// queue/useQueueMutation.ts (interface only)

export interface UseQueueMutationResult<TInput> {
  mutate: (input: TInput) => void;
  isPending: boolean;
  error: unknown;
  reset: () => void;
}

/**
 * Wires a Mutator into TanStack Query's useMutation, with:
 * - per-call client_op_id minted via crypto.randomUUID()
 * - optimistic effect via mutator.onMutate(payload, qc), context restored on
 *   error or replay-rollback
 * - HTTP request via mutator.request(payload, signal); deferred-promise
 *   resolution lets either HTTP or SSE confirm
 * - sessionStorage persistence via attachPersistence (M1.3)
 * - failure-class routing (Validation toast, Infrastructure banner)
 *
 * The hook accepts a static "scope" parameter (e.g., { exposureId,
 * sampleId, experimentId, username, clientId }) that the mutator's
 * onMutate / request need but that aren't part of the per-call input.
 * mutate(input) merges the input into a payload of shape { ...scope,
 * ...input, kind, clientOpId }.
 */
export function useQueueMutation<TInput, TResponse>(
  mutator: Mutator<OpPayload<TInput>, TResponse>,
  scope: Record<string, unknown>,
): UseQueueMutationResult<TInput>;
```

- [ ] **Step 2: Write failing tests**

```typescript
// queue/useQueueMutation.test.ts
describe("useQueueMutation", () => {
  it("mints a fresh clientOpId per mutate() call", () => {
    // Two mutate() calls should produce two distinct deferreds in the registry.
  });
  it("calls mutator.onMutate with merged scope+input payload", () => { /* ... */ });
  it("HTTP success → mutator.onSuccess called with response", async () => { /* ... */ });
  it("HTTP error 4xx → context.restore called (rollback) AND Validation toast emitted", async () => { /* ... */ });
  it("HTTP timeout → retry with exponential backoff up to 30s", async () => { /* ... */ });
  it("SSE event with matching client_op_id confirms before HTTP returns", async () => { /* ... */ });
  it("AbortSignal threads through to mutator.request", async () => { /* ... */ });
  it("isPending reflects MutationCache state", () => { /* ... */ });
  it("rehydrate-on-mount replays persisted ops with original clientOpId", async () => { /* ... */ });
});
```

Each test asserts a specific contract — the canonical "happy path" is the HTTP-succeeds case; everything else exercises one failure mode at a time.

- [ ] **Step 3: Implement**

```typescript
// queue/useQueueMutation.ts (full)
import { useMutation, useQueryClient } from "@tanstack/react-query";
import { newClientOpId } from "../clientOpId";
import { makeDeferred, clearDeferred } from "./deferred";
import { isValidationError, isInfrastructureError } from "./errors";
import { showToast } from "../toast";
import type { Mutator, OpPayload, RollbackContext } from "./types";

export function useQueueMutation<TInput, TResponse>(
  mutator: Mutator<OpPayload<TInput>, TResponse>,
  scope: Record<string, unknown>,
): UseQueueMutationResult<TInput> {
  const qc = useQueryClient();

  const mutation = useMutation<TResponse, unknown, OpPayload<TInput>, RollbackContext>({
    mutationKey: [mutator.kind],
    mutationFn: async (payload, { signal }) => {
      const deferred = makeDeferred<TResponse>(payload.clientOpId);
      // Wire HTTP into the deferred.
      mutator.request(payload, signal)
        .then(response => deferred.resolve(response))
        .catch(err => deferred.reject(err));
      // AbortSignal cleanup: reject deferred so registry doesn't leak.
      signal.addEventListener("abort", () =>
        deferred.reject(new DOMException("aborted", "AbortError")));
      try {
        return await deferred.promise;
      } finally {
        clearDeferred(payload.clientOpId);
      }
    },
    onMutate: (payload) => mutator.onMutate(payload, qc),
    onSuccess: (response, payload) => mutator.onSuccess(payload, response, qc),
    onError: (err, payload, context) => {
      // Roll back optimistic effect.
      context?.restore?.();
      // Route to failure class.
      if (isValidationError(err)) {
        showToast(buildValidationMessage(mutator.kind, err), "error");
      } else if (isInfrastructureError(err)) {
        // Infrastructure banner is mounted at App.tsx; it reads from
        // useMutationState. Nothing per-mutation to do here beyond the
        // built-in TanStack retry behavior.
      }
      // Conflict-class errors don't surface here — they manifest as
      // mutator.onMutate producing a no-op or different effect during
      // replay-rerun.
    },
    retry: (failureCount, err) => {
      if (isValidationError(err)) return false;
      return failureCount < 5;  // exponential backoff up to ~30s
    },
    retryDelay: (attempt) => Math.min(1000 * 2 ** attempt, 30000),
  });

  // Type-erased mutate that mints clientOpId at call time.
  const mutate = (input: TInput) => {
    const payload: OpPayload<TInput> = {
      kind: mutator.kind,
      clientOpId: newClientOpId(),
      ...scope,
      ...input,
    } as OpPayload<TInput>;
    mutation.mutate(payload);
  };

  return { mutate, isPending: mutation.isPending, error: mutation.error, reset: mutation.reset };
}
```

`isValidationError` / `isInfrastructureError` are small helpers in `lib/queue/errors.ts` that classify based on HTTP status (4xx → validation; 5xx / network → infrastructure). `buildValidationMessage` looks up per-kind copy from a small table; add new kinds to the table when introducing new mutators.

- [ ] **Step 4: Verify**

All nine testcases pass. Manual: M2.1's `useAddSampleTag` reference implementation (already templated in M2.1) compiles against this hook signature; if it doesn't, M2.1 needs updating before tasks downstream.

> **Why this is its own task and not folded into M1.1-M1.4:** the hook is the integration surface every consumer (16 mutations) talks to. Specifying it as a footnote in M2.1 (as the prior plan draft did) hides the contract from M1's test surface and risks the framework drifting between what M1 builds and what M2 expects. Putting it in M1.5 means M1 ships a complete, testable framework that M2 strictly consumes.

---

## M1 verification

- [ ] All M1 unit tests green (M1.1–M1.5).
- [ ] No consumer migrations in M1; existing component tests stay green.
- [ ] `npm run build` clean.
- [ ] Module is intentionally dead code at end of M1; M2 begins consumption.

---

## Task M2.1: Trivial mutator slice

The trivial slice (eight mutations) migrates first to validate the framework on simple cases. Each migration is independently mergeable.

Mutations in this slice:
- `useUpdateSample`, `useAddSampleTag`, `useRemoveSampleTag`
- `useAddExposureTag`, `useRemoveExposureTag`
- `usePostSampleMessage`
- `useSetExposureStatus`
- `useSelectExposure`

> **Backend coupling for cross-tab sync.** The corresponding backend routes (`routes_samples.jl`, `routes_exposures.jl` for tags + status + select, `routes_messages.jl` for chat) currently use `log_action!` (no SSE broadcast). The frontend mutator's optimistic UI works without backend changes — clicks feel instant via local cache writes. **But cross-tab sync requires migrating these routes to `apply_event!`** so SSE frames fire and other tabs see updates without manual refetch. Each per-mutation PR in this slice should include the corresponding backend migration:
>
> - `useUpdateSample` → migrate `PATCH /api/samples/:id` body to `apply_event!(kind="update_sample", entity_type="sample", ...)`. New dispatcher branch in `update_view_for_event!` is a no-op (no view writes; the action is already a direct UPDATE on `samples`). The migration is for SSE delivery, not view materialization.
> - `useAddSampleTag` / `useRemoveSampleTag` / `useAddExposureTag` / `useRemoveExposureTag` → similar; new dispatcher branches return `nothing` (tag tables are written directly by the route, not via the dispatcher; we want the SSE side effect).
> - `usePostSampleMessage` → migrate `POST /api/samples/:id/messages` to `apply_event!(kind="add_message", ...)` for cross-tab chat updates.
> - `useSetExposureStatus` → migrate to `apply_event!(kind="set_exposure_status", ...)`.
> - `useSelectExposure` → migrate to `apply_event!(kind="select_exposure", ...)`. Stays LWW per spec; the broadcast is purely for cross-tab "Bob is now looking at exposure 7" awareness.
>
> Without these backend migrations, the trivial slice's optimistic UI works in the editing tab but other tabs stay stale until refetch. The `applyRemoteToCache` switch in M1.2 is forward-scaffolded for these kinds; it activates when each route migrates.

**Files (per mutation):**
- Modify: `packages/HimalayaUI/frontend/src/queries.ts` — replace `useMutation` with queue-shaped mutator
- New: per-mutation mutator file under `packages/HimalayaUI/frontend/src/lib/queue/mutators/<name>.ts`
- New: per-mutation test file under `packages/HimalayaUI/frontend/test/queue/mutators/<name>.test.ts`
- Modify: corresponding `routes_*.jl` route — swap `log_action!` for `apply_event!` (per the backend coupling note above)
- Modify: `packages/HimalayaUI/src/events.jl` — `update_view_for_event!` dispatcher gets new branches that return `nothing` (no view write) for the new kinds; the kinds are still recorded in `user_actions` and broadcast via SSE
- Modify: `packages/HimalayaUI/test/test_events.jl` — `rebuild_views_from_log!` round-trip test extended for each new kind (no-op dispatcher branch must not break the property)
- Modify: any consumer components that used the old hook return shape (most don't change)

- [ ] **Step 1: Write failing tests for `useAddSampleTag` (canonical first migration)**

```typescript
describe("useAddSampleTag mutator", () => {
  it("optimistically appends tag to cache", () => { /* ... */ });
  it("HTTP success replaces optimistic with server response", () => { /* ... */ });
  it("HTTP error rolls back optimistic and surfaces validation toast", () => { /* ... */ });
  it("survives remote SSE event during in-flight (replay-rerun)", () => { /* ... */ });
  it("rehydrates from sessionStorage on reload", () => { /* ... */ });
});
```

- [ ] **Step 2: Implement `useAddSampleTag` mutator**

```typescript
// lib/queue/mutators/addSampleTag.ts
export const addSampleTagMutator: Mutator<AddTagPayload, TagResponse> = {
  kind: "add_tag",
  onMutate: (payload, qc) => {
    const prev = qc.getQueryData<Sample[]>(queryKeys.samples(payload.experimentId));
    qc.setQueryData(queryKeys.samples(payload.experimentId), (old: Sample[] = []) => {
      // Append optimistic tag
      return old.map(s => s.id === payload.sampleId
        ? { ...s, tags: [...(s.tags ?? []), { id: -1, key: payload.key, value: payload.value }] }
        : s);
    });
    return { restore: () => qc.setQueryData(queryKeys.samples(payload.experimentId), prev) };
  },
  request: (payload, signal) => api.addSampleTag(payload.sampleId, payload.key, payload.value,
    authOpts(payload.username, payload.clientId, payload.clientOpId)),
  onSuccess: (payload, response, qc) => {
    // Replace optimistic id with server-assigned id.
    qc.setQueryData(queryKeys.samples(payload.experimentId), (old: Sample[] = []) => {
      return old.map(s => s.id === payload.sampleId
        ? { ...s, tags: (s.tags ?? []).map(t => t.id === -1 ? response : t) }
        : s);
    });
  },
};

export function useAddSampleTag(experimentId: number, sampleId: number) {
  return useQueueMutation(addSampleTagMutator, { experimentId, sampleId });
}
```

> **`useQueueMutation` is the framework-provided wrapper** that wires the mutator into `useMutation` with deferred-promise resolution, sessionStorage persistence, replay-rerun, and failure-class routing. Specified and tested in Task M1.5 (`packages/HimalayaUI/frontend/src/lib/queue/useQueueMutation.ts`). M2.1 onward is strictly consumption.

- [ ] **Step 3: Migrate the remaining seven trivial mutations following the same pattern**

Each migration is its own commit. After each, run the relevant component tests to ensure consumer integration is clean.

- [ ] **Step 4: Verify**

All trivial-mutator tests green. Existing component tests for tags/messages/sample-fields stay green (the optimistic UI is invisible — same observable outcome, faster perception).

---

## Task M2.2: Peak-op atomic slice

This is one PR. It cannot be subdivided. The PR includes backend handler changes, three frontend mutator migrations, autoReanalyze deletion, banner gating update, SSE payload extension, and the failure-class indicator components.

**Files:**
- Modify: `packages/HimalayaUI/src/routes_peaks.jl` — three handlers (POST /:id/peaks, DELETE /peaks/:id, PATCH /peaks/:id) extended with `with_idempotency` + synchronous `analyze_exposure!` + extended response shape
- Modify: `packages/HimalayaUI/src/events.jl` — `broadcast_event!` accepts `post_state` parameter; route handlers pass it for peak events
- New: `packages/HimalayaUI/frontend/src/lib/queue/mutators/peakAdd.ts`, `peakRemove.ts`, `peakSetExcluded.ts`
- Modify: `packages/HimalayaUI/frontend/src/queries.ts` — replace three `useMutation` hooks; delete `autoReanalyze`
- Modify: `packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx` — gate on `useExposureHasPendingPeakOps`; reduce debounce to 150ms
- New: `packages/HimalayaUI/frontend/src/components/ui/Toast.tsx` (or extend existing toast), `InfrastructureBanner.tsx`
- Modify: `packages/HimalayaUI/frontend/src/App.tsx` — mount InfrastructureBanner, ToastContainer
- Modify: `packages/HimalayaUI/test/test_routes_peaks.jl`
- New: per-mutator test files; e2e two-context replay-rerun test

- [ ] **Step 1: Backend route handler changes**

Uses the `defer_broadcast::Bool=true` semantics shipped in M0.6. The route handler defers the broadcast inside the transaction, runs `analyze_exposure!`, then emits one frame with `post_state` populated.

```julia
@post "/api/exposures/{id}/peaks" function(req::HTTP.Request, id::Int)
    db = current_db()
    return with_idempotency(db, req) do
        body = json(req)
        q = Float64(body["q"])
        result = SQLite.transaction(db) do
            event_result = apply_event!(db, req;
                kind            = "peak_added",
                entity_type     = "exposure",
                entity_id       = id,
                payload         = Dict(:q => q),
                defer_broadcast = true)  # M0.6: skip inner broadcast
            analyze_exposure!(db, id, current_analysis_dir(db);
                              trace_known_unchanged = true)
            event_result
        end

        # Build post-analysis payload (after commit; safe because the
        # post-state is read from durable rows, not in-flight state).
        new_hash      = read_inputs_hash(db, id)
        indices_array = serialize_indices_for_post_state(db, id)

        # One enriched broadcast per logical event.
        broadcast_event!(result.event_id, "peak_added", "exposure", id,
                         get_user_id(req, db), get_client_id(req), get_client_op_id(req),
                         JSON3.write(Dict(:q => q));
                         post_state = Dict(
                             :analysis_inputs_hash => new_hash,
                             :indices              => indices_array))

        return HTTP.Response(201;
            headers = ["Content-Type" => "application/json"],
            body = JSON3.write(Dict(
                :peak                  => fetch_curation_row(db, result.view_row_id),
                :event_id              => result.event_id,
                :view_row_id           => result.view_row_id,
                :analysis_inputs_hash  => new_hash,
            )))
    end
end
```

The same shape applies to `DELETE /api/peaks/:id` (kind `peak_removed`) and `PATCH /api/peaks/:id` (kind `peak_excluded` or `peak_unexcluded`); each handler uses `defer_broadcast=true` and emits its own enriched broadcast after `analyze_exposure!` returns.

- [ ] **Step 2: Update Julia tests**

```julia
# test_routes_peaks.jl
@testset "POST /api/exposures/:id/peaks runs synchronous reanalyze" begin
    # Make request, observe response includes analysis_inputs_hash that matches
    # the post-state. Assert auto_peaks/indices are updated synchronously
    # before the response.
end

@testset "POST /api/exposures/:id/peaks SSE frame includes post_state" begin
    # Spy on broadcasts; assert exactly one frame per request, with post_state.
end
```

- [ ] **Step 3: Implement frontend mutators**

Three mutators following the M2.1 pattern. The `onMutate` for each writes optimistic peaks; `onSuccess` writes peaks + indices + analysis_inputs_hash atomically from the response.

- [ ] **Step 4: Delete `autoReanalyze` chain**

Remove [queries.ts:103-114](../../../packages/HimalayaUI/frontend/src/queries.ts) and the `await autoReanalyze(...)` calls in the three peak mutations.

- [ ] **Step 5: Update `StaleIndicesBanner`**

```typescript
// StaleIndicesBanner.tsx
const DEFAULT_STALE_DEBOUNCE_MS = 150;  // was 2000

export function StaleIndicesBanner({ exposureId, debounceMs = DEFAULT_STALE_DEBOUNCE_MS }) {
  const indicesQ = useIndices(exposureId);
  const exposureQ = useExposure(exposureId);
  const reanalyze = useReanalyzeExposure(exposureId ?? 0);
  const hasPendingPeakOps = useExposureHasPendingPeakOps(exposureId);  // NEW

  // Hide entirely while peak ops are pending; cross-entity refetch races
  // are masked.
  if (hasPendingPeakOps) return null;

  // ... existing hash-mismatch logic, with the shorter debounce ...
}
```

- [ ] **Step 6: Implement failure-class indicators**

`Toast.tsx`: Validation-class red toast with mutator-provided message; auto-dismiss after 5s.

`InfrastructureBanner.tsx`: persistent banner showing "Trying to save your changes…" past 500ms, "Couldn't save your changes — refresh and retry" past 30s. Driven by `useMutationState` filtered to long-pending mutations.

- [ ] **Step 7: Write end-to-end test**

The two-context test exercises real cross-context SSE delivery, so it needs a real backend (not `page.route` mocks). The existing `inspect.spec.ts` and `smoke.spec.ts` use mocks; this spec is the first that needs an actual `himalaya serve` instance.

**Test fixture setup:**

A `playwright.config.ts` extension that adds the backend to its `webServer` array (Playwright supports multiple entries; both start before tests, both stop on teardown). Running them as separate processes outside Playwright introduces flakiness on CI (race between Vite startup and backend readiness, no automatic teardown if a test crashes); the `webServer`-array approach gives Playwright lifecycle control.

```typescript
// playwright.config.ts
export default defineConfig({
  webServer: [
    { command: "npm run dev -- --host 127.0.0.1", url: "http://127.0.0.1:5173", reuseExistingServer: !process.env.CI },
    { command: "node ./e2e/start-backend.mjs", url: `http://127.0.0.1:${process.env.E2E_BACKEND_PORT ?? 8081}/api/experiments`, reuseExistingServer: !process.env.CI },
  ],
  // ...
});
```

`./e2e/start-backend.mjs`:
1. Reserves a port via `get-port` (npm package), defaulting to 8081 if free. Writes the chosen port to `process.env.E2E_BACKEND_PORT` for tests.
2. Spawns `julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(["serve", tmpdir, "--port", ENV["E2E_BACKEND_PORT"]])'` with a fixture experiment directory.
3. Forwards the child's stderr to its own stderr so failures surface in Playwright's logs.
4. Cleans up the child + tmpdir on `SIGTERM` / `SIGINT`.

The fixture experiment dir contains a small synthetic `experiment.toml` + a few seeded .dat files with predictable peaks. `e2e/fixtures/experiment-replay-rerun/` is the right place under the existing repo layout.

**Vite proxy:** `vite.config.ts` proxies `/api/*` to `:8080` for dev. For e2e the proxy target needs to be the dynamic `E2E_BACKEND_PORT`. Either:
- Override the proxy target via env var read in `vite.config.ts`, OR
- Build the frontend statically (`npm run build`) and serve via the Julia backend's `dynamicfiles("dist")` path — backend serves both API and frontend on one port, no proxy needed.

The second option is simpler for CI and matches production deployment. Use it.

**Port-conflict handling:** the `get-port` package is not destructive (no `lsof | xargs kill` per CLAUDE.md's anti-pattern note). On a developer machine that happens to have :8081 occupied, the helper picks an alternate.

`packages/HimalayaUI/frontend/e2e/multiplayer-replay-rerun.spec.ts`:

```typescript
import { test, expect, type Browser, type Page } from "@playwright/test";

async function openExposure(browser: Browser, username: string, exposureId: number): Promise<Page> {
  const ctx = await browser.newContext();
  const page = await ctx.newPage();
  await page.goto(`/?exposure=${exposureId}`);
  // OnboardingFlow auto-advances if username param is set.
  await page.fill('[data-testid="username-input"]', username);
  await page.click('[data-testid="onboarding-continue"]');
  await page.waitForSelector('[data-testid="trace-viewer"]');
  return page;
}

test("two users add peaks at non-overlapping q values; both succeed and both visible", async ({ browser }) => {
  const pageA = await openExposure(browser, "alice", 1);
  const pageB = await openExposure(browser, "bob", 1);

  // Initial peak count snapshot.
  const initialPeaksA = await pageA.locator('[data-peak-id]').count();
  const initialPeaksB = await pageB.locator('[data-peak-id]').count();
  expect(initialPeaksB).toBe(initialPeaksA);

  // Alice double-clicks at q ≈ 0.234 to add a peak.
  await pageA.locator('[data-testid="trace-viewer"]').dblclick({ position: { x: 100, y: 200 } });
  // Bob double-clicks at q ≈ 0.567.
  await pageB.locator('[data-testid="trace-viewer"]').dblclick({ position: { x: 250, y: 200 } });

  // Wait for both SSE frames to deliver.
  await expect.poll(() => pageA.locator('[data-peak-id]').count())
    .toBe(initialPeaksA + 2);
  await expect.poll(() => pageB.locator('[data-peak-id]').count())
    .toBe(initialPeaksB + 2);

  // Verify the indices panel reflects both adds in both contexts.
  // Index lists should match between contexts (modulo render order).
  const indicesA = await pageA.locator('[data-index-id]').evaluateAll(els => els.map(e => e.getAttribute('data-index-id')).sort());
  const indicesB = await pageB.locator('[data-index-id]').evaluateAll(els => els.map(e => e.getAttribute('data-index-id')).sort());
  expect(indicesA).toEqual(indicesB);
});

test("two users exclude the same peak; replay-rerun produces a no-op for the second", async ({ browser }) => {
  const pageA = await openExposure(browser, "alice", 1);
  const pageB = await openExposure(browser, "bob", 1);

  // Both pages target the same auto peak.
  const targetPeakSelector = '[data-peak-id]:not([data-excluded="true"]):first-child';

  // Alice excludes; Bob excludes the same peak ~simultaneously.
  // (await both clicks; race semantics handled by replay-as-rerun)
  await Promise.all([
    pageA.locator(targetPeakSelector).click({ modifiers: ["Alt"] }),  // Alt+click = exclude per existing UX
    pageB.locator(targetPeakSelector).click({ modifiers: ["Alt"] }),
  ]);

  // Both pages converge on the peak being excluded.
  await expect.poll(() => pageA.locator(`${targetPeakSelector.replace(":first-child","")}[data-excluded="true"]`).count())
    .toBeGreaterThan(0);
  await expect.poll(() => pageB.locator(`${targetPeakSelector.replace(":first-child","")}[data-excluded="true"]`).count())
    .toBeGreaterThan(0);

  // No Validation toast appears on either page (both ops "succeeded" — the
  // second is a no-op via replay-rerun, not a conflict).
  expect(await pageA.locator('[role="alert"][data-class="validation"]').count()).toBe(0);
  expect(await pageB.locator('[role="alert"][data-class="validation"]').count()).toBe(0);

  // The event log shows both operations recorded (one peak_excluded; the
  // second is dispatcher-deduplicated). Verifiable via DB inspection if the
  // test fixture exposes it.
});
```

The selectors (`[data-testid=...]`, `[data-peak-id]`, `[data-index-id]`, `[data-excluded]`) follow the project's convention of stable `data-*` attributes (per CLAUDE.md "E2E selectors"). If existing components don't surface the right hooks, add them in M2.2 step 3 (frontend mutator implementation).

- [ ] **Step 8: Verify**

Full test suite green: backend, frontend unit, e2e. Manual smoke: open two browser tabs, exclude a peak in one, observe the other tab's banner doesn't flicker and the indices update without a refetch round-trip.

---

## Task M2.3: Index/group ops slice

Three mutations: `useAddIndexToGroup`, `useRemoveIndexFromGroup`, `useDeleteIndex`. Each migration is independently mergeable; they share structure but have different optimistic shapes.

**Files (per mutation):**
- New: `packages/HimalayaUI/frontend/src/lib/queue/mutators/{addIndexToGroup,removeIndexFromGroup,deleteIndex}.ts`
- New: `packages/HimalayaUI/frontend/test/queue/mutators/{...}.test.ts`
- Modify: `packages/HimalayaUI/frontend/src/queries.ts` — replace each `useMutation` with the queue-shaped hook
- Modify: consumers — `IndicesCard.tsx`, `PhasePanel.tsx` (most use the existing hook return shape unchanged)

- [ ] **Step 1: Write failing tests for `useAddIndexToGroup`**

```typescript
describe("useAddIndexToGroup mutator", () => {
  it("optimistically appends index_id to group's members", () => {
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.groups(42), [
      { id: 1, kind: "custom", members: [10, 20] },
    ]);
    const ctx = addIndexToGroupMutator.onMutate(
      { kind: "index_confirmed", clientOpId: "op-1", exposureId: 42, groupId: 1, indexId: 30 } as OpPayload<{ groupId: number; indexId: number }>,
      qc,
    );
    const groups = qc.getQueryData<Group[]>(queryKeys.groups(42));
    expect(groups?.[0].members).toContain(30);
    ctx.restore();
    expect(qc.getQueryData<Group[]>(queryKeys.groups(42))?.[0].members).toEqual([10, 20]);
  });
  it("HTTP success replaces optimistic state with server response", async () => { /* analogous to M2.1 */ });
  it("HTTP 4xx rolls back optimistic and surfaces validation toast", async () => { /* ... */ });
  it("survives remote SSE event during in-flight (replay-rerun)", async () => {
    // Spawn a pending mutation; fire a remote `index_unconfirmed` for a different
    // index in the same group; assert the pending op's optimistic effect re-applies
    // against the post-remote base.
  });
  it("rehydrates from sessionStorage on reload", async () => { /* ... */ });
});
```

Tests for `useRemoveIndexFromGroup` and `useDeleteIndex` follow the same five-test template, swapping the optimistic effect:
- `useRemoveIndexFromGroup`: optimistic effect *removes* an index from the group's members; rollback re-inserts.
- `useDeleteIndex`: optimistic effect removes the index from the indices cache AND from any group's members; rollback restores both. (Server-side: only speculative indices can be deleted; auto indices return 422 → Validation toast.)

- [ ] **Step 2: Implement mutators**

```typescript
// queue/mutators/addIndexToGroup.ts
export const addIndexToGroupMutator: Mutator<
  OpPayload<{ groupId: number; indexId: number }>,
  GroupMemberResponse
> = {
  kind: "index_confirmed",
  onMutate: (payload, qc) => {
    const prev = qc.getQueryData<Group[]>(queryKeys.groups(payload.exposureId));
    qc.setQueryData<Group[]>(queryKeys.groups(payload.exposureId), (old = []) =>
      old.map(g => g.id === payload.groupId
        ? { ...g, members: [...g.members, payload.indexId] }
        : g));
    return { restore: () => qc.setQueryData(queryKeys.groups(payload.exposureId), prev) };
  },
  request: (payload, signal) => api.addIndexToGroup(payload.groupId, payload.indexId,
    authOpts(payload.username, payload.clientId, payload.clientOpId), { signal }),
  onSuccess: (payload, response, qc) => {
    // Server returns the canonical group membership; replace cache entry.
    qc.setQueryData<Group[]>(queryKeys.groups(payload.exposureId), (old = []) =>
      old.map(g => g.id === payload.groupId ? response.group : g));
  },
  affectsExposurePeaks: () => false,
};

export function useAddIndexToGroup(exposureId: number, groupId: number) {
  const username = useAppState(s => s.username);
  return useQueueMutation(addIndexToGroupMutator, {
    exposureId, groupId, username, clientId: CLIENT_ID,
  });
}
```

`removeIndexFromGroup` and `deleteIndex` follow the same pattern with their own optimistic shapes.

> **Why `affectsExposurePeaks: false`:** index/group ops don't change the effective peak set, so they don't gate the speculative-snap query nor the stale banner. Only mutators that affect peaks set this to `true`.

- [ ] **Step 3: Verify**

All 15 testcases (5 per mutation × 3 mutations) pass. Existing `IndicesCard.tsx` and `PhasePanel.tsx` component tests stay green.

---

## Task M2.4: Speculative ops slice

Two mutations + one query gate. The speculative POST is the canonical multi-event idempotency test for `with_idempotency`'s response-body cache.

**Files:**
- New: `packages/HimalayaUI/frontend/src/lib/queue/mutators/createSpeculative.ts`
- New: `packages/HimalayaUI/frontend/test/queue/mutators/createSpeculative.test.ts`
- Modify: `packages/HimalayaUI/src/routes_analysis.jl` — wrap `POST /api/exposures/:id/speculative` body in `with_idempotency` (the route emits N `apply_event!` calls; `with_idempotency` records the response for verbatim retry replay)
- Modify: `packages/HimalayaUI/test/test_routes_analysis.jl` — verify multi-event route returns identical body on retry
- Modify: `packages/HimalayaUI/frontend/src/queries.ts` — `useCreateSpeculative` uses `useQueueMutation`; `useSpeculativeSnap` gates on `useExposureHasPendingPeakOps`
- Modify: `packages/HimalayaUI/frontend/src/components/SpeculativeBuilder.tsx` — render last response with "updating to latest…" subtext during gate

- [ ] **Step 1: Write failing backend tests**

```julia
# test_routes_analysis.jl
@testset "POST /api/exposures/:id/speculative is idempotent under retry" begin
    mktempdir() do tmp
        # ... seed exposure with auto peaks ...
        body = Dict("phase" => "Pn3m", "anchor_peak_id" => peak_id,
                    "anchor_ratio" => 1.0, "additional" => [])
        op_id = "uuid-spec-1"
        # First request creates N speculative indices.
        r1 = HTTP.request("POST", "http://localhost:$port/api/exposures/$exp_id/speculative",
            ["X-Client-Op-Id" => op_id, "Content-Type" => "application/json"],
            JSON3.write(body))
        @test r1.status == 201
        body1 = String(r1.body)
        n_after_first = count_speculative_indices(db, exp_id)

        # Second request with same op_id returns identical body, no new events.
        r2 = HTTP.request("POST", "http://localhost:$port/api/exposures/$exp_id/speculative",
            ["X-Client-Op-Id" => op_id, "Content-Type" => "application/json"],
            JSON3.write(body))
        @test r2.status == 201
        @test String(r2.body) == body1
        @test count_speculative_indices(db, exp_id) == n_after_first
    end
end
```

- [ ] **Step 2: Wrap `routes_analysis.jl`'s speculative POST in `with_idempotency`**

```julia
@post "/api/exposures/{id}/speculative" function(req::HTTP.Request, id::Int)
    db = current_db()
    return with_idempotency(db, req) do
        body = json(req)
        # ... existing validation + iteration over (phase, anchor) combinations ...
        SQLite.transaction(db) do
            for combo in combinations
                apply_event!(db, req;
                    kind            = "speculative_created",
                    entity_type     = "exposure",
                    entity_id       = id,
                    payload         = Dict(...),
                    defer_broadcast = false)  # speculatives don't need post_state
                # ...
            end
        end
        return HTTP.Response(201; body = JSON3.write(response_body))
    end
end
```

The `with_idempotency` wrapper handles both single-event and multi-event routes uniformly — the cached response body is byte-identical regardless of how many events the body emits.

- [ ] **Step 3: Write failing frontend tests**

```typescript
describe("useSpeculativeSnap gates on pending peak ops", () => {
  it("query disabled while a peak_added op is pending for the same exposure", () => {
    const qc = new QueryClient();
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      mutationKey: ["peak_added"],
      variables: { kind: "peak_added", exposureId: 42, q: 1.0 },
    }));
    const { result } = renderHook(() => useSpeculativeSnap(42, "Pn3m", 5, 1.0));
    expect(result.current.fetchStatus).toBe("idle");  // gated
  });
  it("query enabled when no pending peak ops", () => { /* ... */ });
  it("re-enables after pending peak op settles", async () => { /* ... */ });
});

describe("useCreateSpeculative mutator", () => {
  it("optimistically inserts placeholder speculative into indices cache", () => { /* ... */ });
  it("HTTP success replaces placeholder with server-returned indices", async () => { /* ... */ });
  it("idempotent retry returns cached body without re-creating", async () => { /* ... */ });
});

describe("SpeculativeBuilder during gate", () => {
  it("renders last successful snap response", () => { /* ... */ });
  it("shows 'updating to latest…' subtext while gated", () => { /* ... */ });
});
```

- [ ] **Step 4: Implement frontend changes**

```typescript
// queue/mutators/createSpeculative.ts
export const createSpeculativeMutator: Mutator<
  OpPayload<{ phase: string; anchorPeakId: number; anchorRatio: number; additional: SpeculativeAdditional[] }>,
  SpeculativeCreateResponse
> = {
  kind: "speculative_created",
  onMutate: (payload, qc) => {
    // Optimistic: insert a placeholder speculative into indices cache so
    // PhasePanel reflects "you just created a Pn3m" instantly.
    const prev = qc.getQueryData<Index[]>(queryKeys.indices(payload.exposureId));
    qc.setQueryData<Index[]>(queryKeys.indices(payload.exposureId), (old = []) => [
      ...old,
      { id: -Date.now(), phase: payload.phase, basis: payload.anchorRatio, /* placeholder */ } as Index,
    ]);
    return { restore: () => qc.setQueryData(queryKeys.indices(payload.exposureId), prev) };
  },
  request: (payload, signal) => api.createSpeculative(payload.exposureId, {
    phase: payload.phase, anchor_peak_id: payload.anchorPeakId,
    anchor_ratio: payload.anchorRatio, additional: payload.additional,
  }, authOpts(payload.username, payload.clientId, payload.clientOpId), { signal }),
  onSuccess: (payload, response, qc) => {
    // Replace placeholder with server's actual indices array.
    qc.setQueryData(queryKeys.indices(payload.exposureId), response.indices);
    qc.setQueryData(queryKeys.groups(payload.exposureId), response.groups);
  },
  affectsExposurePeaks: () => false,
};
```

```typescript
// queries.ts
export function useSpeculativeSnap(exposureId, phase, anchorPeakId, anchorRatio) {
  const blocked = useExposureHasPendingPeakOps(exposureId);
  return useQuery({
    queryKey: ["exposure", exposureId ?? "none", "speculative-snap", phase ?? "", anchorPeakId ?? -1, anchorRatio] as const,
    queryFn: () => api.getSpeculativeSnap(exposureId as number, phase as string, anchorPeakId as number, anchorRatio),
    enabled: exposureId !== undefined && phase !== undefined && anchorPeakId !== undefined && !blocked,
  });
}
```

```typescript
// SpeculativeBuilder.tsx — pseudocode for the gate UI
const snap = useSpeculativeSnap(exposureId, phase, anchorPeakId, anchorRatio);
const blocked = useExposureHasPendingPeakOps(exposureId);
const lastGood = useRef<SnapResponse | null>(null);
if (snap.data) lastGood.current = snap.data;

return (
  <Card>
    {/* ... selector UI ... */}
    {snap.data && <SnapPreview data={snap.data} />}
    {blocked && lastGood.current && (
      <div>
        <SnapPreview data={lastGood.current} />
        <span className="text-xs text-fg-muted">updating to latest…</span>
      </div>
    )}
  </Card>
);
```

- [ ] **Step 5: Verify**

Backend test passes: speculative POST returns identical body on retry; `count_speculative_indices` is unchanged.
Frontend tests pass: `useSpeculativeSnap` correctly gates and re-enables; `SpeculativeBuilder` shows the "updating to latest…" subtext during the gate window.

---

## Task M2.5: Null-mutator ops slice (`useReanalyzeExposure`)

`useReanalyzeExposure` migrates to queue framework with a null `onMutate`. Button loading state derives from `useQueueOpStatus("reanalyze_exposure", exposureId)`. FIFO ordering with pending peak ops is asserted via the queue's MutationCache iteration order (already validated in M1.2).

**Files:**
- New: `packages/HimalayaUI/frontend/src/lib/queue/mutators/reanalyzeExposure.ts`
- New: `packages/HimalayaUI/frontend/test/queue/mutators/reanalyzeExposure.test.ts`
- Modify: `packages/HimalayaUI/frontend/src/queries.ts` — `useReanalyzeExposure` swaps to `useQueueMutation`
- Modify: `packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx` — button uses `useQueueOpStatus`

- [ ] **Step 1: Write failing tests**

```typescript
describe("useReanalyzeExposure", () => {
  it("button isPending true while op is in queue", () => {
    const { result } = renderHook(() => useReanalyzeExposure(42));
    act(() => result.current.mutate({}));
    expect(result.current.isPending).toBe(true);
  });

  it("button isPending false after op confirms", async () => { /* ... */ });

  it("FIFO ordering: pending peak op enqueued first runs first", async () => {
    // Deterministic ordering via mocked HTTP that resolves in order of
    // mutate() invocation.
    const httpResolveOrder: string[] = [];
    mockApi.addPeak = vi.fn(async (...args) => {
      httpResolveOrder.push("addPeak");
      await new Promise(r => setTimeout(r, 10));
      return { /* ... */ };
    });
    mockApi.reanalyzeExposure = vi.fn(async (...args) => {
      httpResolveOrder.push("reanalyze");
      return { /* ... */ };
    });
    const { result: addPeakHook } = renderHook(() => useAddPeak(42));
    const { result: reanalyzeHook } = renderHook(() => useReanalyzeExposure(42));
    act(() => {
      addPeakHook.current.mutate({ q: 1.0 });
      reanalyzeHook.current.mutate({});
    });
    await waitFor(() => expect(httpResolveOrder).toEqual(["addPeak", "reanalyze"]));
  });

  it("Re-analyze with M0 fast-skip is fast when peak ops already produced latest hash", async () => {
    // After a peak_added op confirms, immediately fire reanalyze; assert the
    // server-side fast-skip path produces a near-no-op response (analyze_run
    // payload has both skip flags true).
  });
});
```

- [ ] **Step 2: Implement mutator**

```typescript
// queue/mutators/reanalyzeExposure.ts
export const reanalyzeExposureMutator: Mutator<
  OpPayload<{}>,
  ReanalyzeResponse
> = {
  kind: "reanalyze_exposure",
  onMutate: () => ({ restore: () => {} }),  // null optimistic effect
  request: (payload, signal) => api.reanalyzeExposure(payload.exposureId,
    authOpts(payload.username, payload.clientId, payload.clientOpId), { signal }),
  onSuccess: (payload, response, qc) => {
    qc.setQueryData(queryKeys.exposure(payload.exposureId), response.exposure);
    qc.setQueryData(queryKeys.indices(payload.exposureId), response.indices);
  },
  affectsExposurePeaks: () => true,  // Reanalyze touches peak set indirectly
};

export function useReanalyzeExposure(exposureId: number) {
  const username = useAppState(s => s.username);
  return useQueueMutation(reanalyzeExposureMutator, {
    exposureId, username, clientId: CLIENT_ID,
  });
}
```

- [ ] **Step 3: Update `StaleIndicesBanner`'s button**

```typescript
const reanalyze = useReanalyzeExposure(exposureId ?? 0);
const status = useQueueOpStatus("reanalyze_exposure", exposureId);

<Button
  variant="primary"
  disabled={status === "pending"}
  onClick={() => reanalyze.mutate({})}
>
  {status === "pending" ? "Re-analyzing…" : "Re-analyze"}
</Button>
```

- [ ] **Step 4: Verify**

All four testcases pass. The FIFO ordering test specifically validates that the queue serializes by mutationCache insertion order, which means the M1.2 invariant test transitively gates this behavior.

> **Note on "FIFO" semantics:** TanStack Query mutations run concurrently by default — the FIFO ordering here refers to *enqueue order in the MutationCache*, which the replay coordinator iterates in. The actual HTTP requests fire concurrently; FIFO ordering matters only for replay-rerun (the order onMutate is re-invoked). For `useReanalyzeExposure` specifically, the M0 fast-skip + server-side `apply_event!` ordering guarantees that even concurrent requests serialize at the SQLite write level.

---

## Task M3.1: Code health cleanup

- [ ] Audit `qc.invalidateQueries` calls in `queries.ts`'s `onSuccess` handlers; replace with `setQueryData` from response where the response carries the new state. Invalidation remains correct for SSE-triggered cache updates and for queries that the response doesn't reach (e.g., cross-entity).
- [ ] In `sseSubscriber.ts`, decide: drop the `client_id`-based skip filter (since queue confirmation now handles same-tab echoes via `client_op_id`) or retain as defense-in-depth. Default: retain with a comment explaining it's a fallback for non-queue echoes.
- [ ] Remove any TODO comments from M0–M2 that have been addressed.

## Task M3.2: Documentation

- [ ] Update `docs/event-log.md` with a new section describing the queue's relationship to the dispatcher contract: events still flow through `apply_event!`; routes now wrap in `with_idempotency`; `analyze_run` events are suppressed when both skip flags fire; SSE frames include `client_op_id`, `ts`, `post_state`.
- [ ] Update `CLAUDE.md` §"HimalayaUI gotchas":
  - Per-mutation `client_op_id` (vs PR #32's per-tab `client_id`) — mint inside `mutationFn`, not at hook creation.
  - `useExposureHasPendingPeakOps` pattern for gating server-computed reads and UI affordances.
  - Queue invariants: `MutationCache` insertion order, deferred-promise registry, sessionStorage rehydrate flow.
  - "Synchronous reanalyze in curation routes" pattern: pass `trace_known_unchanged=true`, return `analysis_inputs_hash` in response.
- [ ] Update `CLAUDE.md` "Current state" section to reflect the queue having shipped.

---

## Verification checklist (full plan)

After M0–M3 ship:

- [ ] All four milestone PRs merged. Each was independently mergeable; on-disk DB heals on `open_db`.
- [ ] Backend test suite: `Pkg.test("HimalayaUI")` — green. New tests: idempotency, fast-skip, broadcast suppression, post_state, GC sweep.
- [ ] Frontend test suite: `npm test` — green. New tests: queue framework (deferred, replay coordinator, persistence, hooks), per-mutator round-trip tests for all 16 mutations.
- [ ] Frontend build: `npm run build` — TS strict + Vite clean.
- [ ] Two-context Playwright multiplayer replay-rerun test green.
- [ ] Manual smoke checklist:
  - [ ] Click "exclude peak" — peak disappears within one frame, no banner flicker, indices update without observable round-trip.
  - [ ] Two browser tabs (different `client_id`s, same `username`): exclude peak in tab A → tab B reflects within ~1s without manual refresh.
  - [ ] Trigger network failure (offline mode in DevTools) during a curation — Validation toast or Infrastructure banner appears appropriately.
  - [ ] Reload tab mid-curation — sessionStorage queue rehydrates, op replays via HTTP, cache settles correctly.
  - [ ] Speculative builder modal: open while a peak op is pending — modal shows "updating to latest…" briefly, then snap suggestions populate.
- [ ] Latency observability: `analyze_run` event payloads in `user_actions` show fast-skip path P99 < 100µs in steady state.
- [ ] `post_state` size telemetry: query `user_actions` for the last 100 `analyze_run` events; compute distribution of `payload->>'post_state_size_bytes'`. Confirm P50 < 3KB and P99 < 8KB. Feed back into Open Question #4 disposition: if outlier sizes are infrequent (P99 < 8KB), `post_state` enrichment is sustainable; if P99 ≥ 8KB, switch to compact "you should refetch the following keys" payload per the OQ #4 fallback. Document the decision in M3.2's docs update.
- [ ] Cached-response staleness check: a Julia regression test asserts `idempotent_responses` returns its body verbatim even after a referenced entity (peak, index) is deleted by another user — verifies the spec's documented eventual-consistency model.
- [ ] Cascading-rejection check: a Vitest scenario seeds N queued ops referencing one entity, then deletes that entity via a remote SSE event; assert each pending op surfaces its own Validation toast (not a single cascade summary).

If any of the manual smoke checks fail, refer to the spec's [Fallback triggers](../specs/2026-05-02-mutation-queue-design.md#fallback-triggers).

---

## Out of scope

- Cross-session undo / Cmd+Z (deferred follow-up; primitives in place).
- `Last-Event-ID` server-side SSE replay (made unnecessary by HTTP-authoritative confirmation).
- Offline editing / localStorage queue persistence.
- Polished UI affordances (per-element pending pulse, aggregate "saving…" widget, ambient-change highlight). Failure-class indicators (Validation toast, Infrastructure banner) ship in M2.2; the polish indicators are a separate follow-up PR.
- Time-travel UI / state reconstruction at arbitrary historical event id.
- ML-driven conflict resolution.
- Authentication tightening beyond `X-Username` (orthogonal).
- `If-Match` headers + 409 retry on delta-shaped routes (R5b — redirected; replay-as-rerun handles the conflict-resolution problem with different tradeoffs).
