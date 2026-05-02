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
- `packages/HimalayaUI/src/server.jl` — periodic GC sweep of `idempotent_responses` rows older than TTL (1 hour)

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
- `packages/HimalayaUI/test/test_pipeline.jl` — fast-skip path is microseconds, slow path triggers when `add` curation present, `trace_known_unchanged=true` skips `hash_trace_file`
- `packages/HimalayaUI/frontend/test/api.test.ts` — `X-Client-Op-Id` header on mutations + GET-omit case

Modified docs:
- `docs/event-log.md` — new section on `client_op_id` semantics, `idempotent_responses` table, broadcast suppression rule
- `CLAUDE.md` §"HimalayaUI gotchas" — note the per-mutation UUID pattern (vs PR #32's per-tab UUID)

**M1: Frontend infrastructure (no consumers)**

New frontend:
- `packages/HimalayaUI/frontend/src/lib/queue/types.ts` — `OpKind`, `Mutator<T,R>`, `RollbackContext`, `PendingDeferred<T>`, `SseEvent`
- `packages/HimalayaUI/frontend/src/lib/queue/deferred.ts` — `pendingDeferreds` Map, `makeDeferred`, registry helpers
- `packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts` — `handleRemoteEvent`, `applyRemoteToCache`
- `packages/HimalayaUI/frontend/src/lib/queue/persistence.ts` — sessionStorage mirror, rehydrate flow
- `packages/HimalayaUI/frontend/src/lib/queue/hooks.ts` — `useExposureHasPendingPeakOps`, `useQueueOpStatus`
- `packages/HimalayaUI/frontend/src/lib/queue/index.ts` — barrel export

Modified frontend:
- `packages/HimalayaUI/frontend/src/lib/sseSubscriber.ts` — wire the queue's `handleRemoteEvent` ahead of the existing fallback path. The `client_id`-based skip filter remains as the fallback for non-queue echoes (system events, pre-feature events).

New tests:
- `packages/HimalayaUI/frontend/test/queue/replayCoordinator.test.ts` — pop-on-client-op-id, rollback-apply-replay sequence, MutationCache iteration order
- `packages/HimalayaUI/frontend/test/queue/persistence.test.ts` — rehydrate cycle, schema-version drop with toast
- `packages/HimalayaUI/frontend/test/queue/hooks.test.ts` — `useExposureHasPendingPeakOps` reflects pending mutations
- `packages/HimalayaUI/frontend/test/queue/deferred.test.ts` — registry mint/resolve/leak-on-abort

**M2: Mutator migration**

The peak-op slice is one atomic PR (Migration step 2 below). Other slices land independently in any order.

Modified backend (peak-op slice):
- `packages/HimalayaUI/src/routes_peaks.jl` — three curation routes wrap body in `with_idempotency`, call `analyze_exposure!(db, id, dir; trace_known_unchanged=true)` synchronously inside the transaction, return extended response shape
- `packages/HimalayaUI/src/events.jl` — `broadcast_event!` adds `post_state` field to SSE frame for events that carry it (curation events emitted from peak-edit routes after synchronous reanalyze)

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
- `packages/HimalayaUI/test/test_routes_peaks.jl` — synchronous reanalyze inside curation handlers; response shape includes `analysis_inputs_hash` and `peak_or_curation_row`

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

The `analyze_exposure!` refactor is behaviorally identical for callers that don't pass `trace_known_unchanged=true` — existing callers (`reingest`, scheduled scans, manual `useReanalyzeExposure`) preserve current behavior. The fast-skip path activates only when curation routes opt in (M2 peak-op slice).

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
    # Set up exposure with auto_peaks + exclude curations.
    # Compute hash via DB-only path; compute hash via trace-loading path.
    # Assert they match.
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

## Task M0.6: `post_state` payload + size observability

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl` — emit `post_state_size_bytes` on analyze_run
- Modify: `packages/HimalayaUI/src/events.jl` — `broadcast_event!` includes `post_state` for events that opt in
- Modify: `packages/HimalayaUI/src/routes_peaks.jl` — peak-edit routes (when wired in M2) build `post_state` from the new analysis state
- Modify: `packages/HimalayaUI/test/test_sse.jl`

`post_state` enrichment is the load-bearing M0 prerequisite for M2's replay-without-refetch correctness. Without it, replay coordinator falls back to invalidate-and-refetch and the cross-entity refetch race re-introduces banner flicker.

- [ ] **Step 1: Write failing tests**

```julia
@testset "post_state included in SSE frame for curation events with reanalyze" begin
    # A peak_added event from the M2 wiring (handler builds post_state) → frame has it.
    # Pre-M2 events (no post_state passed) → field absent.
end

@testset "post_state_size_bytes recorded on slow-path analyze_run" begin
    # After analyze_exposure! slow path, the analyze_run event payload
    # has post_state_size_bytes > 0.
end
```

- [ ] **Step 2: Implement**

Extend `broadcast_event!` to accept an optional `post_state` parameter:

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
    # ... rest unchanged ...
end
```

`apply_event!` doesn't always know what `post_state` should be — the caller does. Add an optional `post_state` keyword arg to `apply_event!`:

```julia
function apply_event!(db, req;
                      kind, entity_type, entity_id,
                      payload = nothing,
                      undoes_event_id = nothing,
                      post_state = nothing)
    # ...
    broadcast_event!(...; post_state = post_state)
    # ...
end
```

In `analyze_exposure!`'s slow path, compute `post_state_size_bytes` and embed in the `analyze_run` payload:

```julia
indices_array = serialize_indices_for_post_state(db, exposure_id)
size_bytes = length(JSON3.write(indices_array))
apply_event!(db, _system_request();
    kind = "analyze_run",
    # ...
    payload = Dict(
        :findpeaks_skipped     => false,
        :indexpeaks_skipped    => false,
        :duration_ms           => duration_ms,
        :post_state_size_bytes => size_bytes,
    ))
```

The full `post_state` field (with the actual indices array) is built and passed by the *route handler* in M2, not by `analyze_exposure!`.

- [ ] **Step 3: Verify**

Tests pass.

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
// queue/replayCoordinator.test.ts — six cases:
// 1. SSE with matching client_op_id resolves the pending deferred
// 2. SSE without matching client_op_id triggers rollback-apply-replay
// 3. Rollback iterates pending mutations in reverse order
// 4. Re-run iterates in original (insertion) order
// 5. MutationCache.getAll() insertion order is preserved (Set iteration)
// 6. AbortSignal during replay does not double-rollback
```

- [ ] **Step 2: Implement**

```typescript
// queue/replayCoordinator.ts
import type { QueryClient } from "@tanstack/react-query";
import type { MutationCache } from "@tanstack/react-query";
import type { SseEvent } from "./types";
import { getDeferred } from "./deferred";

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

  // Apply remote to cache (cache-merge from event payload + post_state).
  applyRemoteToCache(remote, qc);

  // Re-run optimistic effects in queue order.
  for (const m of pending) {
    m.options.onMutate?.(m.state.variables);
  }
}

function synthesizeResponseFromSse(remote: SseEvent): unknown {
  // Build a response object matching what the HTTP would have returned.
  // Per route, this means {event_id, view_row_id, analysis_inputs_hash, post_state}.
  return {
    event_id: remote.id,
    client_op_id: remote.client_op_id,
    analysis_inputs_hash: remote.post_state?.analysis_inputs_hash,
    // Per-route fields filled in via remote.payload.
    ...((remote.payload as object) ?? {}),
  };
}

function applyRemoteToCache(remote: SseEvent, qc: QueryClient): void {
  // Per kind, write to the appropriate cache entries via setQueryData.
  // - peak_added: append to peaks(exposureId) cache; if post_state, replace
  //   indices and exposure.analysis_inputs_hash atomically.
  // - peak_excluded: similar, removing the peak from cache.
  // - index_confirmed: update group_members cache.
  // - ...etc per kind.
  // Falls back to qc.invalidateQueries for events without sufficient payload.
  switch (remote.kind) {
    // ...
  }
}
```

- [ ] **Step 3: Verify**

All six testcases pass. The MutationCache iteration-order test is the load-bearing assertion against TanStack version drift — see Open Question #7 in the spec.

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
// Three cases:
// 1. Returns false when no peak ops are pending for the exposure
// 2. Returns true while a peak_added/peak_excluded/etc op is pending
// 3. Returns false again after the op confirms
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

## M1 verification

- [ ] All M1 unit tests green.
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

**Files (per mutation):**
- Modify: `packages/HimalayaUI/frontend/src/queries.ts` — replace `useMutation` with queue-shaped mutator
- New: per-mutation mutator file under `packages/HimalayaUI/frontend/src/lib/queue/mutators/<name>.ts`
- New: per-mutation test file under `packages/HimalayaUI/frontend/test/queue/mutators/<name>.test.ts`
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

> **`useQueueMutation` is the framework-provided wrapper** that wires the mutator into `useMutation` with deferred-promise resolution, sessionStorage persistence, and replay-rerun. Defined in `lib/queue/index.ts` as part of M1; tests for it are part of M1.2.

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

Update each of the three peak-edit routes per the spec sketch:

```julia
@post "/api/exposures/{id}/peaks" function(req::HTTP.Request, id::Int)
    db = current_db()
    return with_idempotency(db, req) do
        body = json(req)
        q = Float64(body["q"])
        result = SQLite.transaction(db) do
            event_result = apply_event!(db, req;
                kind        = "peak_added",
                entity_type = "exposure",
                entity_id   = id,
                payload     = Dict(:q => q))
            analyze_exposure!(db, id, current_analysis_dir(db);
                              trace_known_unchanged = true)
            event_result
        end

        # Build extended response with post-analysis state.
        new_hash = read_inputs_hash(db, id)
        indices_array = serialize_indices_for_post_state(db, id)

        # Re-broadcast with post_state.
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

> **Subtle point on broadcast ordering:** the `apply_event!` call inside the transaction already broadcasts a frame *without* `post_state` (the `apply_event!` call site doesn't know about post-analyze state). The handler then issues a *second* broadcast with `post_state` after the analysis. Two options here:
>
> 1. Suppress the inner broadcast for events whose route will issue an enriched broadcast. Add a flag.
> 2. Always include `post_state` in the inner broadcast. Move post-analysis to before the broadcast — but `apply_event!` and `analyze_exposure!` need to be in the same transaction, and the broadcast fires after commit.
>
> The cleanest option is (1): add an `outer_broadcast::Bool = false` keyword to `apply_event!`. When true (the M2 wiring), the inner broadcast is suppressed; the route handler issues the enriched one explicitly. Default false preserves all other call sites. Document the rule.

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

`packages/HimalayaUI/frontend/e2e/multiplayer-replay-rerun.spec.ts`:

```typescript
test("two contexts add peaks to same exposure; both succeed", async ({ browser }) => {
  const ctxA = await browser.newContext();
  const ctxB = await browser.newContext();
  // ... two pages ...
  // Both add peaks at non-overlapping q values.
  // Wait for both SSE frames to land.
  // Assert both peaks appear in both contexts; indices reflect both.
});

test("two contexts exclude same peak; replay-rerun resolves invisibly", async ({ browser }) => {
  // ...
});
```

- [ ] **Step 8: Verify**

Full test suite green: backend, frontend unit, e2e. Manual smoke: open two browser tabs, exclude a peak in one, observe the other tab's banner doesn't flicker and the indices update without a refetch round-trip.

---

## Task M2.3: Index/group ops slice

Three mutations: `useAddIndexToGroup`, `useRemoveIndexFromGroup`, `useDeleteIndex`. Pattern follows M2.1; each migration is independent.

- [ ] **Step 1**: Write failing tests per mutator.
- [ ] **Step 2**: Implement mutators following the M2.1 pattern.
- [ ] **Step 3**: Verify; existing component tests stay green.

---

## Task M2.4: Speculative ops slice

Two mutations: `useCreateSpeculative`, plus the speculative-snap query gating.

- [ ] **Step 1: Write failing tests**

```typescript
describe("useSpeculativeSnap gates on pending peak ops", () => {
  it("query disabled while peak op pending", () => { /* ... */ });
  it("modal renders last response with 'updating to latest' subtext during gate", () => { /* ... */ });
  it("re-enables after pending peak ops settle", () => { /* ... */ });
});
```

- [ ] **Step 2: Migrate `useCreateSpeculative` to queue-shaped mutator**

The speculative create-and-delete multi-event case is the canonical test for `with_idempotency`'s response-body-cache path. Confirm that the route handler (already wrapping in `with_idempotency` per Task M0.3) returns the cached body byte-identically on retry.

- [ ] **Step 3: Update `useSpeculativeSnap`**

```typescript
export function useSpeculativeSnap(exposureId, phase, anchorPeakId, anchorRatio) {
  const blocked = useExposureHasPendingPeakOps(exposureId);
  return useQuery({
    queryKey: ...,
    queryFn: ...,
    enabled: enabled && !blocked,
  });
}
```

- [ ] **Step 4: Update `SpeculativeBuilder.tsx`**

While the snap query is disabled (gated), render the last-known response with subtle "updating to latest…" subtext rather than a skeleton or blank state.

- [ ] **Step 5: Verify**

---

## Task M2.5: Null-mutator ops slice

`useReanalyzeExposure`. Migrates to queue framework with a null `onMutate`. Button loading state derives from `useQueueOpStatus("reanalyze_exposure", exposureId)`.

- [ ] **Step 1: Write tests**

```typescript
it("Re-analyze enters queue behind pending peak ops", () => { /* ... */ });
it("button shows pending state while op is in queue", () => { /* ... */ });
it("FIFO ordering: pending peak ops complete before reanalyze fires", () => { /* ... */ });
```

- [ ] **Step 2: Implement**

```typescript
export const reanalyzeExposureMutator: Mutator<ReanalyzePayload, ReanalyzeResponse> = {
  kind: "reanalyze_exposure",
  onMutate: () => ({ restore: () => {} }),  // null optimistic effect
  request: (payload, signal) => api.reanalyzeExposure(payload.exposureId,
    authOpts(payload.username, payload.clientId, payload.clientOpId)),
  onSuccess: (payload, response, qc) => {
    qc.setQueryData(queryKeys.exposure(payload.exposureId), response.exposure);
    qc.setQueryData(queryKeys.indices(payload.exposureId), response.indices);
  },
};
```

- [ ] **Step 3: Update `StaleIndicesBanner`'s button**

```typescript
const status = useQueueOpStatus("reanalyze_exposure", exposureId);
<Button disabled={status === "pending"} onClick={() => reanalyze.mutate({...})}>
  {status === "pending" ? "Re-analyzing…" : "Re-analyze"}
</Button>
```

- [ ] **Step 4: Verify**

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
- [ ] `post_state` size telemetry: M2 reads `post_state_size_bytes` from `analyze_run` events; confirms typical sizes are within budget (1-3KB) and outlier sizes (≥8KB) are infrequent enough that the SSE frame growth is sustainable.

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
