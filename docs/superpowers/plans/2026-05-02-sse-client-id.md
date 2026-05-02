# Per-Tab SSE Client Identity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement [docs/superpowers/specs/2026-05-02-sse-client-id-design.md](../specs/2026-05-02-sse-client-id-design.md): separate audit identity (`actor`/`username`) from SSE routing identity (`client_id`) so multi-tab same-user workflows (Compare + Inspect, multi-monitor split) see live updates across tabs without sacrificing self-echo suppression within a tab.

**Architecture:** Additive change at every layer. New nullable column `user_actions.client_id`, new `X-Client-Id` header (mutations only), new per-tab UUID minted into `sessionStorage` on first load. SSE filter swaps from `actor === username` to `client_id === clientId`. `actor` remains on the SSE frame for future presence/attribution UI.

**Tech stack:** Julia 1.9+, SQLite.jl, Oxygen.jl, HTTP.jl (SSE), JSON3; React 18, TanStack Query v5, Vitest, Playwright.

**Read first:**
- [docs/superpowers/specs/2026-05-02-sse-client-id-design.md](../specs/2026-05-02-sse-client-id-design.md) — design rationale and open questions
- [docs/event-log.md](../../event-log.md) §"Client side" — current SSE filter behavior
- [docs/superpowers/plans/2026-05-01-multiplayer-instrumentation.md](2026-05-01-multiplayer-instrumentation.md) §"Migration Architecture" — idempotency contract, ALTER TABLE try/catch precedent, sentinel patterns

---

## File Map

**Modified backend:**
- `packages/HimalayaUI/src/db.jl` — append `client_id TEXT` to `user_actions` `CREATE TABLE` DDL in `create_schema!`; append one `ALTER TABLE` to the `stmts` array in `migrate_schema!`
- `packages/HimalayaUI/src/actions.jl` — new `get_client_id(req)` helper
- `packages/HimalayaUI/src/events.jl` — `apply_event!` extracts + persists `client_id`; `broadcast_event!` adds parameter + JSON field

**Modified frontend:**
- `packages/HimalayaUI/frontend/src/api.ts` — `AuthOpts.clientId`, `X-Client-Id` header on mutations
- `packages/HimalayaUI/frontend/src/queries.ts` — `authOpts(username, clientId)` extension; module-level `CLIENT_ID` const; 17 `authOpts(username)` call sites updated to `authOpts(username, CLIENT_ID)`
- `packages/HimalayaUI/frontend/src/lib/sseSubscriber.ts` — filter swap (`actor` → `client_id`), `CurationEvent.client_id`, ctx shape change
- `packages/HimalayaUI/frontend/src/App.tsx` — capture `clientId` as a stable `const` (no ref); remove `usernameRef` declaration, sync useEffect, and the `useRef` import (all unused after the swap)

**New frontend:**
- `packages/HimalayaUI/frontend/src/lib/clientId.ts` — mint + persist UUID in sessionStorage

**New tests:**
- `packages/HimalayaUI/test/test_actions.jl` (does NOT currently exist; must be created and `include`d in `runtests.jl`) — `get_client_id` cases
- `packages/HimalayaUI/frontend/test/clientId.test.ts`

**Modified tests:**
- `packages/HimalayaUI/test/test_events.jl` — `apply_event!` writes/reads column; rebuild round-trip with NULL `client_id`
- `packages/HimalayaUI/test/test_sse.jl` — extend ALL existing `broadcast_event!` testsets to pass + assert `client_id` (4 call sites: lines 25, 61, 96, 140 — verify with grep before editing)
- `packages/HimalayaUI/test/test_db.jl` — fresh-DB column presence + legacy-schema migration
- `packages/HimalayaUI/test/runtests.jl` — `include("test_actions.jl")`
- `packages/HimalayaUI/frontend/test/api.test.ts` — `X-Client-Id` parallel assertions on every existing `X-Username` mutation test (6 mutation cases at lines 7, 53, 68, 118, 144, 161 + the GET-omit case at line 16 + new clientId-without-username case)
- `packages/HimalayaUI/frontend/test/sse.test.tsx` (NOTE: existing file is `sse.test.tsx`, not `sseSubscriber.test.ts`) — update all 6 existing cases to use new `{ clientId, qc }` ctx shape; add new `client_id` filter cases

**Modified docs:**
- `docs/event-log.md` §"Client side" — replace `actor`-based filter description, drop multi-tab caveat
- `CLAUDE.md` §"HimalayaUI frontend gotchas" — replace multi-tab staleness note with new filter semantics + sessionStorage lifetime

---

## Migration Architecture

This plan introduces one schema change. It must follow the idempotency contract from [Plan 7's migration architecture](2026-05-01-multiplayer-instrumentation.md), summarized:

- `migrate_schema!` must be safe on fresh / pre-migration / partially-migrated / already-migrated DBs.
- `ALTER TABLE ADD COLUMN` is wrapped in `try ... catch end` (SQLite raises on duplicate column; precedent at [db.jl:142-160](../../../packages/HimalayaUI/src/db.jl)).
- No backfill needed — existing rows get NULL `client_id`, which is the correct semantic value (system-emitted or pre-feature events broadcast to all tabs).
- Forward-only; rollback is via DB backup before deploy.

---

## Task 1: Schema change

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl`
- Modify: `packages/HimalayaUI/test/test_db.jl`

- [ ] **Step 1: Write failing tests**

Add to `packages/HimalayaUI/test/test_db.jl`:

```julia
@testset "user_actions.client_id column exists on fresh DB" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        cols = Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(user_actions)"))
        @test any(c -> c.name == "client_id", cols)
        SQLite.close(db)
    end
end

@testset "open_db adds client_id to legacy user_actions schema" begin
    mktempdir() do tmp
        path = joinpath(tmp, "test.db")
        # Build a legacy user_actions table without client_id. Mirror the
        # CURRENT create_schema! DDL ([db.jl:146-156]) verbatim minus the
        # new column — column names, defaults, and FKs must match exactly,
        # otherwise the test exercises a schema that never shipped.
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
                undoes_event_id INTEGER REFERENCES user_actions(id)
            )
        """)
        SQLite.close(db)

        # Re-open via open_db: migrate_schema!'s ALTER TABLE adds the column
        db = HimalayaUI.open_db(path)
        cols = Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(user_actions)"))
        @test any(c -> c.name == "client_id", cols)
        SQLite.close(db)
    end
end
```

> **Pre-migration fixture pattern:** the project does not use `ALTER TABLE … DROP COLUMN` (requires SQLite ≥ 3.35; not portable across SQLite.jl builds). Construct legacy schemas via raw `CREATE TABLE` directly. Mirror the current `user_actions` DDL minus the new column.

Run: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'` — expect failure.

- [ ] **Step 2: Add column to `create_schema!` and `migrate_schema!`**

In `packages/HimalayaUI/src/db.jl`:

1. Add `client_id TEXT` to the `user_actions` `CREATE TABLE` DDL inside `create_schema!` so fresh DBs get the column without relying on the migration loop.

2. Append one entry to the `stmts` array in `migrate_schema!` ([db.jl:200-215](../../../packages/HimalayaUI/src/db.jl)):

```julia
"ALTER TABLE user_actions ADD COLUMN client_id TEXT",
```

The surrounding try/catch loop handles the duplicate-column case for already-migrated DBs (precedent: every other ALTER in this list).

> **No separate `migrate_rN_*` helper.** Per-column migration helpers are reserved for table-rebuild migrations (R2's peaks split, the PK widening). Single ADD COLUMN goes in the `stmts` array. See [docs/superpowers/specs/2026-05-02-sse-client-id-design.md](../specs/2026-05-02-sse-client-id-design.md) §"Schema".

- [ ] **Step 3: Verify**

Run the test above; it should pass. Also run the full HimalayaUI test suite — no other test should regress because `client_id` is nullable and unread by existing code.

---

## Task 2: Backend extraction + persistence

**Files:**
- Modify: `packages/HimalayaUI/src/actions.jl`
- Modify: `packages/HimalayaUI/src/events.jl`
- New: `packages/HimalayaUI/test/test_actions.jl`
- Modify: `packages/HimalayaUI/test/test_events.jl`
- Modify: `packages/HimalayaUI/test/test_sse.jl` (broadcast frame field)
- Modify: `packages/HimalayaUI/test/runtests.jl` (`include("test_actions.jl")`)

- [ ] **Step 1: Write failing tests**

`test_actions.jl` does NOT currently exist. Create it and add `include("test_actions.jl")` to `runtests.jl`. New file:

```julia
@testset "get_client_id" begin
    @test HimalayaUI.get_client_id(
        HTTP.Request("POST", "/x", ["X-Client-Id" => "abc-123"], UInt8[])
    ) == "abc-123"
    @test HimalayaUI.get_client_id(
        HTTP.Request("POST", "/x", Pair{String,String}[], UInt8[])
    ) === nothing
    @test HimalayaUI.get_client_id(
        HTTP.Request("POST", "/x", ["X-Client-Id" => ""], UInt8[])
    ) === nothing
end
```

In `packages/HimalayaUI/test/test_events.jl` — use the standard fixture chain (mirrors [test_events.jl:97-99](../../../packages/HimalayaUI/test/test_events.jl)). FK enforcement is on, so `sample_id=1` won't work without a real experiment+sample row.

> **Important:** call `HimalayaUI.bind_db!(db)` after `open_db`. `apply_event!` triggers `broadcast_event!` which calls `current_db()` for the `lookup_username` step. Without `bind_db!`, `current_db()` errors; the error is caught and logged as `@warn` (events.jl:60-64), so tests still pass but emit noise. Mirror the test_sse.jl convention.

```julia
@testset "apply_event! persists client_id when X-Client-Id header present" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        eid    = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="e1")
        req = HTTP.Request("POST", "/x",
            ["X-Username" => "alice", "X-Client-Id" => "tab-xyz"], UInt8[])
        HimalayaUI.apply_event!(db, req;
            kind="noop_test", entity_type="exposure", entity_id=eid, payload=Dict(:k=>1))
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT client_id FROM user_actions ORDER BY id DESC LIMIT 1"))
        @test rows[1].client_id == "tab-xyz"
    end
end

@testset "apply_event! writes NULL client_id when header absent" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        eid    = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="e1")
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        HimalayaUI.apply_event!(db, req;
            kind="noop_test", entity_type="exposure", entity_id=eid, payload=Dict(:k=>1))
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT client_id FROM user_actions ORDER BY id DESC LIMIT 1"))
        @test ismissing(rows[1].client_id)
    end
end

@testset "rebuild_views_from_log! tolerates rows with NULL client_id" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        eid    = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="e1")
        DBInterface.execute(db, """
            INSERT INTO user_actions (user_id, action, entity_type, entity_id, payload, client_id)
            VALUES (NULL, 'noop_test', 'exposure', ?, '{}', NULL)
        """, [eid])
        @test_nowarn HimalayaUI.rebuild_views_from_log!(db, Int(eid))
    end
end
```

In `packages/HimalayaUI/test/test_sse.jl` — there are **4 call sites** to `broadcast_event!` (lines 25, 61, 96, 140; grep to confirm). All four must be updated to pass `client_id` as the new 7th positional argument, otherwise the suite fails to compile after the signature change.

For the first testset (`"broadcast_event! puts a curation frame onto subscriber channels"`), change the call and add an assertion:

```julia
HimalayaUI.broadcast_event!(
    1, "test_broadcast", "exposure", 42,
    user_id, JSON3.write(Dict(:foo => "bar")), "tab-xyz")

@test occursin("\"client_id\":\"tab-xyz\"", frame)
```

For the other three testsets (`"with nil user_id sets actor to null"` at line 49+, `"prunes closed subscriber"` at line 96+, `"prunes slow subscriber"` at line 140+), pass `nothing` as the 7th argument. For the nil-user_id testset, also add `@test occursin("\"client_id\":null", frame)` since that case exercises the NULL-client_id wire format.

The existing test uses `HimalayaUI.bind_db!(db)` (NOT `set_current_db!` — does not exist) and accesses `HimalayaUI.SSE_SUBSCRIBERS[]` / `HimalayaUI.SSE_LOCK` directly. Mirror that style in any additions.

Run tests; expect failures.

- [ ] **Step 2: Add `get_client_id` to `actions.jl`**

```julia
"""
    get_client_id(req) -> Union{String, Nothing}

Return the `X-Client-Id` header value if present and non-empty, else nothing.
This is the per-tab SSE routing identity, distinct from `X-Username` (audit
identity). See docs/superpowers/specs/2026-05-02-sse-client-id-design.md.
"""
function get_client_id(req::HTTP.Request)
    v = HTTP.header(req, "X-Client-Id", "")
    isempty(v) ? nothing : String(v)
end
```

- [ ] **Step 3: Plumb through `apply_event!`**

In `packages/HimalayaUI/src/events.jl::apply_event!`:

1. At the **top of the function body**, alongside the existing `username = get_username(req)` line ([events.jl:25](../../../packages/HimalayaUI/src/events.jl:25)), add:

```julia
client_id = get_client_id(req)
```

   It must be extracted **before** `SQLite.transaction` opens, so the value is in scope for both the INSERT (inside the transaction) and the post-commit `broadcast_event!` call.
2. Add `client_id` to the INSERT. Current ([events.jl:32-35](../../../packages/HimalayaUI/src/events.jl:32)) is 6 columns: `(user_id, action, entity_type, entity_id, payload, undoes_event_id)`. New shape (binding order MUST match column order):

```julia
res = DBInterface.execute(db,
    """INSERT INTO user_actions
       (user_id, action, entity_type, entity_id, payload, undoes_event_id, client_id)
       VALUES (?, ?, ?, ?, ?, ?, ?)""",
    [user_id, kind, entity_type, Int(entity_id), payload_json, undoes_event_id, client_id])
```

3. Pass `client_id` to `broadcast_event!` as a new 7th positional argument (after `payload_json`).

4. **`_system_request()` requires no change.** It produces requests with no `X-Client-Id` header, so `get_client_id` returns `nothing` and pipeline-emitted events broadcast with NULL `client_id` — correct (system events refresh all tabs). See spec §"Behavior notes".

5. **No SSE handler changes needed.** The `GET /api/events` handler at [server.jl:47-77](../../../packages/HimalayaUI/src/server.jl) emits only heartbeats and broadcast frames — no initial/replay/hello frame to backfill with `client_id`. Verified before writing this plan.

- [ ] **Step 4: Extend `broadcast_event!`**

Two locations to update — the function definition and its single call site in `apply_event!`:

```julia
# events.jl: function definition (currently 6 args, becomes 7)
function broadcast_event!(event_id::Integer, kind::String, entity_type::String,
                          entity_id::Integer, user_id::Union{Integer, Nothing},
                          payload_json::Union{String, Nothing},
                          client_id::Union{String, Nothing})
    actor = user_id === nothing ? nothing : lookup_username(current_db(), user_id)
    msg = JSON3.write(Dict(
        :id          => Int(event_id),
        :kind        => kind,
        :entity_type => entity_type,
        :entity_id   => Int(entity_id),
        :actor       => actor,
        :client_id   => client_id,
        :payload     => payload_json === nothing ? nothing : JSON3.read(payload_json),
    ))
    # ... rest unchanged ...
end
```

Update the call site at [events.jl:61](../../../packages/HimalayaUI/src/events.jl:61): add `client_id` as the 7th positional argument.

Update the docstring at [events.jl:184](../../../packages/HimalayaUI/src/events.jl:184) to document the new parameter.

- [ ] **Step 5: Verify**

Run HimalayaUI test suite. The `rebuild_views_from_log!` round-trip property test should still pass — `client_id` doesn't affect view dispatch. If it fails, the dispatcher is reading something it shouldn't.

---

## Task 3: Frontend mint + plumb

**Files:**
- New: `packages/HimalayaUI/frontend/src/lib/clientId.ts`
- New: `packages/HimalayaUI/frontend/test/clientId.test.ts`
- Modify: `packages/HimalayaUI/frontend/src/api.ts`
- Modify: `packages/HimalayaUI/frontend/src/queries.ts`
- Modify: `packages/HimalayaUI/frontend/test/api.test.ts`

- [ ] **Step 1: Write failing tests for `clientId.ts`**

`packages/HimalayaUI/frontend/test/clientId.test.ts`:

```typescript
import { describe, it, expect, beforeEach, vi } from "vitest";

describe("getClientId", () => {
  beforeEach(() => sessionStorage.clear());

  it("mints a UUID on first call and persists to sessionStorage", async () => {
    const { getClientId } = await import("../src/lib/clientId");
    const id = getClientId();
    expect(id).toMatch(/^[0-9a-f-]{36}$/);
    expect(sessionStorage.getItem("himalaya.client_id")).toBe(id);
  });

  it("returns the same id on subsequent calls within a session", async () => {
    const { getClientId } = await import("../src/lib/clientId");
    expect(getClientId()).toBe(getClientId());
  });

  it("reuses an existing id in sessionStorage", async () => {
    sessionStorage.setItem("himalaya.client_id", "preset-uuid");
    vi.resetModules();
    const { getClientId } = await import("../src/lib/clientId");
    expect(getClientId()).toBe("preset-uuid");
  });
});
```

> **Note:** `crypto.randomUUID` is available in JSDOM 22+. If the test environment lacks it, add a `vi.stubGlobal("crypto", { randomUUID: () => "fake-uuid-..." })` in `test/setup.ts` — but check first; do not add unnecessary stubs.

- [ ] **Step 2: Implement `clientId.ts`**

```typescript
const KEY = "himalaya.client_id";

export function getClientId(): string {
  let id = sessionStorage.getItem(KEY);
  if (!id) {
    id = crypto.randomUUID();
    sessionStorage.setItem(KEY, id);
  }
  return id;
}
```

- [ ] **Step 3: Write failing tests for `X-Client-Id` header**

In `packages/HimalayaUI/frontend/test/api.test.ts`, add parallel assertions next to each existing `X-Username` mutation assertion (6 mutation `it()` blocks at lines 7, 53, 68, 118, 144, 161, plus the GET-omit case at line 16). **Do NOT add `clientId` to GET-mode tests** (`getTrace`, `listIndices`, `listGroups`, `listSampleMessages`) — they assert header *exclusion*, and `X-Client-Id` follows the same exclusion rule as `X-Username`. Example:

```typescript
it("sends X-Client-Id header on mutating calls when opts.clientId provided", async () => {
  // ... existing fetch mock ...
  await api.addPeak(1, 0.5, { username: "alice", clientId: "tab-xyz" });
  expect((init.headers as Record<string, string>)["X-Client-Id"]).toBe("tab-xyz");
});

it("omits X-Client-Id on GET requests", async () => {
  // mirror the existing X-Username GET test
});

it("sends X-Client-Id without X-Username when only clientId is provided", async () => {
  // Anonymous tab (pre-onboarding) still gets a client_id for SSE filtering.
  await api.addPeak(1, 0.5, { clientId: "tab-xyz" });
  const headers = init.headers as Record<string, string>;
  expect(headers["X-Client-Id"]).toBe("tab-xyz");
  expect(headers["X-Username"]).toBeUndefined();
});
```

Add `clientId: "tab-xyz"` to the `AuthOpts` argument in every existing mutation test that asserts on `X-Username`. The new assertion goes in the same `it` block.

- [ ] **Step 4: Extend `AuthOpts` and `request()` in `api.ts`**

```typescript
export interface AuthOpts {
  username?: string;
  clientId?: string;
}

// Inside request():
if (opts?.username && method !== "GET") headers["X-Username"] = opts.username;
if (opts?.clientId && method !== "GET") headers["X-Client-Id"] = opts.clientId;
```

- [ ] **Step 5: Extend `authOpts()` in `queries.ts`**

```typescript
function authOpts(
  username: string | undefined,
  clientId: string | undefined,
): api.AuthOpts {
  const out: api.AuthOpts = {};
  if (username !== undefined) out.username = username;
  if (clientId !== undefined) out.clientId = clientId;
  return out;
}
```

> Note: do NOT shorthand to `{ username, clientId }` — `exactOptionalPropertyTypes: true` rejects `undefined` values for optional fields. The conditional-build pattern is mandatory; see CLAUDE.md gotcha.

Update all 17 call sites. The pattern shifts from `authOpts(username)` to `authOpts(username, CLIENT_ID)`.

Define `CLIENT_ID` once at module load in `queries.ts`:

```typescript
import { getClientId } from "./lib/clientId";
const CLIENT_ID = getClientId();
```

Do **not** call `getClientId()` inside each hook — it's stable for the tab session. The function is sessionStorage-memoized, but the JS call still has cost on every render. Module-level constant is the right shape.

- [ ] **Step 6: Verify**

Run `npm test`. All `api.test.ts` cases pass. All existing `queries.ts` consumers continue to compile under TS strict.

---

## Task 4: SSE filter swap

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/sseSubscriber.ts`
- Modify: `packages/HimalayaUI/frontend/src/App.tsx`
- Modify: `packages/HimalayaUI/frontend/test/sse.test.tsx` (existing file with 6 cases — extend, do NOT create a parallel `sseSubscriber.test.ts`)

- [ ] **Step 1: Update existing tests + add new filter tests**

`packages/HimalayaUI/frontend/test/sse.test.tsx` already exists with 6 cases that all pass `{ username, qc }` as the ctx. Step 2's ctx-shape change (`{ clientId, qc }`) breaks every one of them. Two cases need more than a mechanical ctx update:

1. **Line 21 `"skips self-echo events whose actor matches local username"` — DELETE this test.** Its semantics invert under the new filter: an event with `actor: "alice"` from another tab of the same user should now invalidate, not be skipped. Replace with the new "processes events from another tab even when same actor" case below.
2. **Line 48 `"propagates anonymous events (actor: null) from other tabs"` — KEEP unchanged (only update ctx shape).** The intent (events with no actor pass through) survives the filter swap. The NEW system-events test below (`"processes system events (NULL client_id) regardless of tab"`) covers the NULL-client_id branch specifically; keeping the original case wording avoids a confusing near-duplicate test name.

For the other 4 cases (lines 6, 31, 38, 58), pass a placeholder like `clientId: "test-tab"` — they exercise unrelated branches (entity_id missing, malformed JSON, non-exposure entity_type, base invalidation) and don't care about filter identity.

Then add the new filter-specific cases below to the same file:

```typescript
import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { handleCurationEvent } from "../src/lib/sseSubscriber";
import { queryKeys } from "../src/queries";

function makeCtx(clientId: string | undefined) {
  const qc = new QueryClient();
  vi.spyOn(qc, "invalidateQueries");
  return { qc, ctx: { qc, clientId } };
}

describe("handleCurationEvent", () => {
  it("skips events authored by this tab (client_id match)", () => {
    const { qc, ctx } = makeCtx("tab-A");
    handleCurationEvent(JSON.stringify({
      kind: "peak_added", entity_type: "exposure", entity_id: 5,
      actor: "alice", client_id: "tab-A",
    }), ctx);
    expect(qc.invalidateQueries).not.toHaveBeenCalled();
  });

  it("processes events from another tab even when same actor", () => {
    const { qc, ctx } = makeCtx("tab-A");
    handleCurationEvent(JSON.stringify({
      kind: "peak_added", entity_type: "exposure", entity_id: 5,
      actor: "alice", client_id: "tab-B",
    }), ctx);
    expect(qc.invalidateQueries).toHaveBeenCalledWith({ queryKey: queryKeys.peaks(5) });
  });

  it("processes system events (NULL client_id) regardless of tab", () => {
    const { qc, ctx } = makeCtx("tab-A");
    handleCurationEvent(JSON.stringify({
      kind: "analyze_run", entity_type: "exposure", entity_id: 5,
      actor: null, client_id: null,
    }), ctx);
    expect(qc.invalidateQueries).toHaveBeenCalled();
  });

  it("ignores malformed JSON", () => {
    const { qc, ctx } = makeCtx("tab-A");
    handleCurationEvent("{not json", ctx);
    expect(qc.invalidateQueries).not.toHaveBeenCalled();
  });
});
```

> The existing file is `sse.test.tsx` (NOT `sseSubscriber.test.ts`). Confirmed by `ls packages/HimalayaUI/frontend/test/sse*`. Both `handleCurationEvent` (the function under test) and the new test cases live alongside the existing 6.

- [ ] **Step 2: Update `CurationEvent` interface and filter**

In `packages/HimalayaUI/frontend/src/lib/sseSubscriber.ts`:

```typescript
interface CurationEvent {
  kind?: string;
  entity_type?: string;
  entity_id?: number;
  actor?: string | null;
  client_id?: string | null;
}

export function handleCurationEvent(
  data: string,
  ctx: { clientId: string | undefined; qc: QueryClient },
): void {
  let event: CurationEvent | null = null;
  try {
    event = JSON.parse(data) as CurationEvent;
  } catch {
    return;
  }
  if (!event || typeof event.entity_id !== "number") return;
  // Self-echo filter: skip events authored by this tab. Other tabs of the
  // same user (or other users) pass through. System events (client_id=null)
  // also pass through — they originated outside any tab.
  if (event.client_id && event.client_id === ctx.clientId) return;
  if (event.entity_type !== "exposure") return;
  const id = event.entity_id;
  ctx.qc.invalidateQueries({ queryKey: queryKeys.peaks(id) });
  ctx.qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
  ctx.qc.invalidateQueries({ queryKey: queryKeys.groups(id) });
  ctx.qc.invalidateQueries({ queryKey: queryKeys.exposure(id) });
}
```

- [ ] **Step 3: Update `App.tsx` to thread `clientId` into ctx**

`clientId` is **stable for the tab session** (sessionStorage), so a `useRef` is unnecessary. Capture once at component mount:

```typescript
import { getClientId } from "./lib/clientId";
// inside the App component:
const clientId = getClientId();  // stable, no ref needed
// pass clientId directly into the SSE handler ctx instead of usernameRef.
```

The existing `usernameRef` ([App.tsx:21](../../../packages/HimalayaUI/frontend/src/App.tsx:21)) exists because `username` changes during onboarding (`undefined → "alice"`). After this change, `username` is no longer read by the SSE handler. `usernameRef` is the **only** `useRef` in App.tsx (verified). Remove:

1. The `usernameRef` declaration and its sync useEffect.
2. The `useRef` import — App.tsx line 2 currently reads `import { useEffect, useRef } from "react";`. Change to `import { useEffect } from "react";` (keep `useEffect`, drop `useRef`).
3. Any reads of `usernameRef.current` (today: only the EventSource handler at lines 26-35, which is being rewritten to use `clientId` directly).

- [ ] **Step 4: Verify**

Run `npm test` — all unit tests pass. Then manual verification:

1. Start `himalaya serve` and `npm run dev`. Open two tabs as `alice`. Edit a peak in tab A; tab B's PhasePanel and trace plot should refresh without manual reload.
2. With one tab, edit a peak; verify TanStack Query devtools show no SSE-driven duplicate invalidation beyond the local mutation's `onSuccess` invalidations. (Devtools may show multiple events per cache key — what matters is that the SSE handler does NOT contribute additional invalidations for the originating tab's own edits.)
3. (If Compare exists by then.) Open Compare with exposure Y in tab A and Inspect for Y in tab B. Edit a peak in tab B; verify Compare's Y panel updates.

---

## Task 5: Documentation

**Files:**
- Modify: `docs/event-log.md`
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update `docs/event-log.md` §"Client side"**

Replace the `actor`-based filter description (lines ~174-197) with the `client_id`-based one. Drop the multi-tab caveat about `alice-laptop` vs `alice-desktop`. Note that `actor` remains on the SSE frame for future presence/attribution UI.

- [ ] **Step 2: Update `CLAUDE.md` §"HimalayaUI frontend gotchas"**

The current `StaleIndicesBanner` section doesn't mention this directly, but a new gotcha is warranted near the SSE / Zustand notes:

```markdown
**Per-tab SSE identity.** SSE self-echo filtering uses a per-tab `client_id`
minted into `sessionStorage` on first load (see `lib/clientId.ts`). Audit
identity (`actor` / `X-Username`) is unchanged. Two tabs of the same user
are treated as distinct subscribers — edits in one tab refresh the other.
The `client_id` lives for the tab session: survives reload, scoped to one
tab. See docs/event-log.md §"Client side".
```

---

## Verification checklist

- [ ] `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'` passes (Julia: `test_db`, `test_actions`, `test_events`).
- [ ] `npm test` (from `packages/HimalayaUI/frontend/`) passes (`api.test.ts`, `clientId.test.ts`, `sse.test.tsx`).
- [ ] `npm run build` passes (TypeScript strict + Vite build).
- [ ] Manual two-tab verification per Task 4 Step 4.
- [ ] `git diff` shows: 1 schema migration, 1 backend helper + plumbing, 1 frontend module, 17 `authOpts` call-site updates, 1 SSE filter swap, 2 doc updates. No surprise files.

## Out of scope

Per the spec:
- `If-Match` / 409 conflict resolution (Plan 7 R5b — gated on R4 instrumentation).
- Long-lived `client_id` in `localStorage`. Tab-session scoping is sufficient.
- Presence UI ("Bob is editing"). The `actor` field is still on the wire if a future spec wants it.
- E2E multi-tab Playwright scenarios. Defer until Compare's E2E suite needs them.
