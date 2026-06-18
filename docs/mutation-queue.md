# Mutation queue architecture

This is the load-bearing reference for HimalayaUI's mutation pipeline: how a
user gesture in the browser becomes a durable row in SQLite, how concurrent
edits from other tabs are reconciled into the local cache without a refetch
round-trip, and what invariants every queue mutator must respect.

The queue is implemented across:

- **Frontend:** `packages/HimalayaUI/frontend/src/lib/queue/`
- **Backend:** `packages/HimalayaUI/src/idempotency.jl`, `events.jl`

It composes with — does not replace — the dispatcher contract in
[event-log.md](event-log.md). Read that first if you haven't.

---

## 1. The two coordinating tokens

Two UUIDs flow through every mutation. They live at different scopes and
solve different problems:

| Token | Scope | Purpose | Source |
|---|---|---|---|
| `client_id` | Per-tab (sessionStorage) | SSE self-echo filtering: "is this frame from me?" | `lib/clientId.ts` |
| `client_op_id` | Per-mutation (fresh per `mutate()` call) | Idempotency key + deferred-promise registry key | `lib/clientOpId.ts` |

The same user across two tabs produces *two distinct* `client_id`s — by
design. Edits in one tab refresh the other.

The same mutation retried produces the *same* `client_op_id` (the framework
keeps the same payload across HTTP retries). A different `mutate()` call
mints a fresh one.

> **Invariant: mint `client_op_id` inside `mutationFn`, not at hook
> construction.** `useQueueMutation` mints it via `newClientOpId()` per
> `mutate()` call (`useQueueMutation.ts:147`). Capturing one in a closure at
> hook mount time would make every retry of every mutation through that hook
> share one id — the server's idempotency cache would treat the second user
> action as a duplicate of the first.

See also: [event-log.md §"Client side"](event-log.md) for `client_id`
lifecycle (per-tab, survives reload).

---

## 2. Anatomy of a queue mutator

Every queue-able operation is a `Mutator<TInput, TScope, TResponse>`
(`lib/queue/types.ts`). Its fields:

```
kind                  // OpKind discriminator (e.g. "peak_added")
onMutate(payload, qc) // optimistic cache write; returns RollbackContext
request(payload, sig) // HTTP call; honours AbortSignal
onSuccess(payload, response, qc)  // canonical cache write
synthesizeFromSse?(remote, base)  // SSE-wins synthetic response (see §3)
affectsExposurePeaks?(payload, exposureId)  // hook-level scoping
treats404AsSuccess?   // treat HTTP 404 as a no-op success (idempotent removes)
```

`useQueueMutation(mutator, scope)` wires it into TanStack Query's
`useMutation`. `scope` is closure-injected at hook construction (e.g.
`exposureId`, `username`, `clientId`); per-call input is merged into the
flat payload before any callback runs. Mutator code never has to cast.

**`OpKind` ≠ event kind in two cases:** `delete_index` emits
`speculative_deleted`; `reanalyze_exposure` emits `analyze_run`. The
`OpKind` describes the *user gesture*; the event name describes the
*durable mutation*. Membership predicates like
`useExposureHasPendingPeakOps` operate on `OpKind`.

---

## 3. The happy-path lifecycle

A successful mutation flows through eight steps:

```
1. User clicks → useQueueMutation.mutate(input) is called
2. Hook mints fresh client_op_id, builds the FlatPayload
3. TanStack's mutationFn:
     a. Registers a deferred-promise in pendingDeferreds keyed on client_op_id
     b. Fires mutator.request(payload, signal) — async HTTP
     c. await deferred.promise   ← does NOT await HTTP directly
4. mutator.onMutate(payload, qc) writes the optimistic cache effect
5. Server: route body wrapped in with_idempotency(db, req) do … end
     a. Outer SQLite.transaction
     b. apply_event!(InTransaction(), …) writes user_actions row + view
     c. Response cached in idempotent_responses keyed by client_op_id
     d. Tx commits → _flush_post_commit_broadcasts!() fires SSE frame
6. SSE frame arrives at all subscribed tabs (including the originating one)
7. Originating tab: handleRemoteEvent matches client_op_id →
     a. applyPostStateOnly(remote, qc) — refresh indices + analysis_inputs_hash
     b. deferred.resolve(synthesizeResponseFromSse(remote))
     c. deferred.controller.abort() — kill the redundant in-flight HTTP
8. mutationFn awakens → mutator.onSuccess writes canonical cache state
   (server-assigned positive id replaces the negative optimistic placeholder)
```

The deferred-promise pattern is what makes this work. `mutationFn` doesn't
`return await mutator.request(…)`; it returns the deferred's promise.
Whoever resolves it first — HTTP response or SSE frame — wins the race. The
loser is aborted (HTTP) or becomes a no-op (SSE arriving second hits an
empty registry).

### Why SSE-first matters

If we awaited HTTP directly and triggered `onSuccess` from the response, the
indices cache wouldn't update until the next `post_state` SSE frame landed
in the *foreign-event* path — but that path skips when `client_id` matches.
The stale-indices indicator would then stick at the pre-mutation hash until a
polling refetch (or page navigation) replaced it. Resolving the deferred
from SSE — and applying `post_state` *before* resolving — closes that
window deterministically.

### Why we abort HTTP after SSE wins

The HTTP response, if it ever arrives, would be redundant (the SSE-
synthesized response already carries the same `event_id` /
`analysis_inputs_hash` / kind-specific fields). Aborting saves a network
round-trip's worth of bytes, and it eliminates a subtle ordering hazard:
without the abort, an HTTP-resolution that races the SSE-resolution can
land its `onSuccess` write *after* a foreign event has already been
applied, double-applying the optimistic effect.

> **Invariant: resolve, then abort.** `replayCoordinator.ts:53-57`. Promises
> only settle once — calling `abort()` first triggers the registered
> `onAbort` listener in `mutationFn`, which rejects the deferred, beating
> our `resolve()`.

---

## 4. Optimistic placeholder ids are negative

Mutators that need an entity id before the server has assigned one (e.g.
`peak_added` inserting a `peak_curations` row, `speculative_created`
inserting an `indices` row) write a **negative integer** id into the
optimistic cache row. SQLite's `INTEGER PRIMARY KEY AUTOINCREMENT` never
returns ≤ 0 for fresh inserts, so positive vs negative is an unambiguous
"server-confirmed" vs "still pending" flag.

`onSuccess` overwrites the placeholder with the real id from the server
response (or, on the SSE-first path, from the synthesized response built
out of the SSE payload).

**Consumer rules.** Any cache reader that filters or compares peak ids
must tolerate negatives:

- Don't `id > 0`-filter — strips legitimate optimistic rows.
- Don't strict-parse against a backend URL pattern without a sign check.
- Don't dereference into URLs (`/api/peaks/${id}`) without first confirming
  the id is positive — a negative id will 404.

This is documented as a wall comment at the top of `types.ts`. Breaking it
surfaces as flicker or duplicate rows during the optimistic-confirm
transition; the regression cost is high relative to the gain from a
cleaner-looking filter.

---

## 5. Foreign events: replay-as-rerun

When an SSE frame arrives whose `client_op_id` does *not* match any
registered deferred, it's either from another user, another tab, or a
system-initiated event. The cache has to incorporate the foreign effect
*plus* whatever optimistic mutations are still pending in this tab — and
do it in an order that preserves causal correctness.

The algorithm (`replayCoordinator.ts:76-117`):

```
1. pending := MutationCache.getAll().filter(status == "pending")
2. for m in pending.reverse():     // last-in first-out
       m.context.restore()         // roll back optimistic effect
3. applyRemoteToCache(remote, qc)  // fold in the foreign effect
4. for m in pending:                // insertion order
       fresh := m.options.onMutate(m.state.variables)
       m.state.context = fresh     // swap in fresh restore closure
```

Two ordering invariants are load-bearing:

> **Invariant: rollback in REVERSE; re-apply in INSERTION order.** Pending
> mutations form a dependency chain — later ops can read state earlier ops
> wrote. Reverse rollback unwinds dependents before their dependencies;
> forward re-apply rebuilds them on top of the new foreign baseline.

> **Invariant: `MutationCache.getAll()` preserves insertion order.** If a
> future TanStack version changes this — e.g. switches to a hash-keyed map
> — the foreign-event branch silently produces wrong cache state. There's a
> regression test (`replayCoordinator.test.ts`, M1.2) that pins this
> against version drift; if you bump TanStack, run it.

### Why we swap fresh restore closures back in

TanStack v5 captures `state.context` once from the original `onMutate`
return. Re-running `onMutate` after rollback produces a *new* rollback
closure that reverts to the post-foreign-event snapshot, but TanStack
doesn't update `state.context` on its own. If a later HTTP-settle's
`onError` fires and uses the original (pre-rollback) context, it
overwrites the cache with stale state. The surgical assignment at
`replayCoordinator.ts:112` keeps the per-mutation rollback consistent.

### Why each iteration is independently try/caught

If one mutator's re-run `onMutate` throws, the *other* pending mutations
must still re-run — otherwise their `state.context` retains the
pre-rollback closure and a later error rolls back to a stale snapshot.
The per-iteration try/catch (`replayCoordinator.ts:107-116`) is what makes
the fan-out resilient to a single bad mutator.

### Self-echo guard

There's a third path between own-op and foreign-event: an SSE frame whose
`client_id` matches this tab but whose `client_op_id` doesn't match any
registered deferred. This happens when HTTP-first won the race — the
mutator's `onSuccess` already wrote the canonical row and cleared the
deferred, but the SSE frame still arrives. Falling into the foreign-event
path here would double-apply the per-kind body (e.g. duplicate peak row
for `peak_added`). Falling into nothing would skip the `post_state` indices
update.

The compromise: `applyPostStateOnly(remote, qc)` and return early. The
indices cache and `analysis_inputs_hash` are propagated; the per-kind body
is skipped because it already ran in `onSuccess`.

---

## 6. The backend half: `with_idempotency`

`with_idempotency(db, req) do … end` (`src/idempotency.jl`) wraps a route
handler body with three guarantees:

1. **Request-level idempotency.** With `X-Client-Op-Id`, the first
   successful response (status < 400) is INSERTed into
   `idempotent_responses` keyed by `client_op_id` *atomically with the
   events the body emitted*. Subsequent calls with the same op-id return
   the cached body without re-executing.
2. **Concurrent retries serialize.** A per-op-id `ReentrantLock`
   (`OP_LOCKS`, guarded by `OP_LOCKS_MU`) ensures two racing tasks with
   the same op-id execute the body exactly once. The cache is checked
   once outside the lock (fast path) and again inside the lock + tx
   (double-check) so the fast path stays cheap.
3. **Failures aren't cached.** Status ≥ 400 responses don't write
   `idempotent_responses`; the next retry re-executes the body. This
   matches Stripe's idempotency-key semantics.

Without `X-Client-Op-Id` the body still runs inside an outer
`SQLite.transaction` — bodies that call `apply_event!(InTransaction(), …)`
need that outer tx whether or not idempotency is engaged.

### `apply_event!` has two methods

- `apply_event!(db, req; …)` — the public method. Opens its own savepoint
  + broadcasts immediately on commit. Use this from non-idempotent routes.
- `apply_event!(InTransaction(), db, req; …)` — the in-tx variant.
  Participates in the outer `SQLite.transaction` and **never broadcasts on
  its own**. Use this from inside `with_idempotency`.

> **Invariant: events emitted inside `with_idempotency` use the
> `InTransaction` variant.** Calling the public method from within would
> open a nested savepoint (correctness-correct, but inefficient) and would
> bypass the post-commit broadcast queue (correctness-broken — subscribers
> would see a frame referring to state that may roll back).

### Post-commit broadcast queue

Events emitted inside the outer tx must defer their SSE frame until *after*
the tx commits, otherwise subscribers could observe state that rolls back.
Each request handler gets its own queue via `task_local_storage()` (key
`POST_COMMIT_BROADCAST_KEY`):

- On commit + status < 400: `_flush_post_commit_broadcasts!()` fires every
  queued frame.
- On rollback OR status ≥ 400 OR cache replay (body didn't execute):
  `_clear_post_commit_broadcasts!()` discards the queue without firing.

Broadcast failures are logged but do not abort the flush — frames are
best-effort, durable rows in `user_actions` are not.

### Single-process scope

`OP_LOCKS` is in-process only. The deployment model is one-experiment-per-
process (see [CLAUDE.md](../CLAUDE.md) "Running the app"); a multi-process or
sharded deployment would need a different primitive. Two paths:

1. Change the cache `INSERT` to `INSERT OR IGNORE` and on `changes()=0`
   re-read the cached row.
2. Replace `OP_LOCKS` with a row-level lock (e.g. `SELECT … FOR UPDATE`
   on a Postgres backend).

The `idempotent_responses(client_op_id)` PK is defense-in-depth for either
path.

### GC

`gc_idempotent_responses!(db; ttl_seconds=3600)` sweeps the cache + the
matching `OP_LOCKS` entries under `OP_LOCKS_MU`. Mid-body ops have no
response row yet (the row is INSERTed *inside* the body's tx), so they're
invisible to the sweep — the sweep can only collect locks for ops that
have COMPLETED past their TTL. This avoids the race where a naive sweep
would delete a lock a fresh retry is about to acquire for a never-completed
op.

`server.jl::start_gc_timer!` wires this into a background timer in
`serve(db)`.

---

## 7. Synchronous reanalyze in curation routes

Peak `add` / `remove` / `exclude` routes pass `trace_known_unchanged=true`
to `analyze_exposure!` so it skips `findpeaks` (the .dat file hasn't
changed since the trace hash was recorded) but still recomputes
`analysis_inputs_hash` and the index set. The route returns the new hash
in the response body, and the queue mutator's `onSuccess` writes it into
the exposure cache directly.

The result: no refetch round-trip, no stale-indices flicker. New
curation routes should follow this pattern; otherwise the stale-indices
indicator spuriously fires until a polling refetch lands.

The `analyze_run` no-op fast path also suppresses both the SSE frame and
the durable `user_actions` row when both `findpeaks_skipped` and
`indexpeaks_skipped` are true (the hashes already prove no-op-ness; a
durable count of "nothing happened" offers no observability value).

---

## 8. The `useExposureHasPendingPeakOps` gate

While a peak-affecting mutation is queued for an exposure, any UI
affordance that depends on the *effective peak set* must be suspended —
otherwise the user sees flicker as the optimistic state, the server
response, and the SSE echo land out of order.

Canonical consumers:

- the stale-indices indicator on the Focus/sample surface — don't show "stale" while a re-analyze is mid-flight.
- `useSpeculativeSnap` — don't snap to a peak that may roll back.

Reuse the hook for any new card or query that reads `peaks(exposureId)`
derivatively. The membership predicate uses `affectsExposurePeaks` on the
mutator definition, so a new peak-touching mutator participates
automatically once it sets that flag.

---

## 9. Persistence across reload

Pending mutations are mirrored into `sessionStorage` via
`lib/queue/persistence.ts`. On reload, the queue rehydrates the pending
set and resumes (HTTP retries pick up where they left off; SSE frames
that landed during the reload window have their `client_op_id` matched
against the rehydrated registry).

`schemaVersion` mismatch between the persisted shape and the current
build drops the queue with a toast — the alternative (silent type drift)
would corrupt the cache.

---

## 10. Failure classes + retry policy

`useQueueMutation` distinguishes two failure classes via
`lib/queue/errors.ts`:

| Class | Examples | Retry? | UX |
|---|---|---|---|
| Validation | 4xx, schema rejection, duplicate ↔ existing row | No | Toast with kind-specific message |
| Infrastructure | 5xx, network error, timeout | Yes (up to 5x, exponential backoff capped at 30s) | Banner mounted at `print/App.tsx` (PrintApp) via `useMutationState` |

A validation error rolls back the optimistic effect via `context.restore`
and the user sees a toast. The mutation's `state` ends in `error`.

An infrastructure error keeps the optimistic effect *visible* and retries
in the background — the banner indicates "trying"; the user can keep
working. After 5 failures, the mutation lands in `error` and the
optimistic effect rolls back.

`InfrastructureBanner` (`print/shell/InfrastructureBanner.tsx`) is a global
status strip that reads from `useMutationState` directly. It hides until
any mutation has been pending > 500ms (avoids flicker on fast successes)
and upgrades to a "stuck" state with a Refresh button after > 30s.

---

## 11. Multi-layer reconciliation contract

Six layers participate in any kind's optimistic-and-confirmation flow:

```
1. Route emit (Julia)            packages/HimalayaUI/src/routes_*.jl
2. SSE frame payload             packages/HimalayaUI/src/events.jl::broadcast_event!
3. applyRemoteToCache merge      lib/queue/applyRemoteToCache.ts
4. Cache row shape               (the resulting Peak / Index / Exposure object)
5. Mutator onMutate              lib/queue/mutators/*.ts (optimistic write)
6. Mutator onSuccess             lib/queue/mutators/*.ts (canonical write)
```

A bug class (e.g. auto/curation peak-id collision, missing `intensity`
field, source-mismatch dedupe) can manifest from any of them. When fixing
one, the regression test must exercise *every* layer where the class can
occur. See [contract-testing.md](contract-testing.md) for the testing
discipline and the canonical paired test files.

---

## File map

**Frontend** (`packages/HimalayaUI/frontend/src/`):

| File | Role |
|---|---|
| `lib/queue/types.ts` | `OpKind`, `Mutator`, `RollbackContext`, `PendingDeferred`, `SseEvent`, `FlatPayload` |
| `lib/queue/useQueueMutation.ts` | Framework wrapper around `useMutation` — deferred resolution, AbortSignal handling, retry policy |
| `lib/queue/deferred.ts` | `pendingDeferreds` Map + lifecycle helpers |
| `lib/queue/replayCoordinator.ts` | `handleRemoteEvent`: own-op confirmation + foreign-event replay-as-rerun |
| `lib/queue/applyRemoteToCache.ts` | Per-kind cache merge for foreign events; `applyPostStateOnly` for own-op + self-echo |
| `lib/queue/persistence.ts` | sessionStorage mirror + rehydrate |
| `lib/queue/hooks.ts` | `useExposureHasPendingPeakOps`, `useQueueOpStatus` |
| `lib/queue/errors.ts` | `isValidationError`, `isInfrastructureError`, `buildValidationMessage` |
| `lib/queue/mutatorRegistry.ts` | Discovery / lookup of mutators by kind |
| `lib/queue/optimisticId.ts` | Negative-id minting helper |
| `lib/queue/peakQTol.ts` | Tolerance for matching peak rows by q value |
| `lib/queue/mutators/` | One file per mutator kind (peak / index / speculative / custom-index / series / assignment / trivial / reanalyze) |
| `lib/clientId.ts` | Per-tab UUID (sessionStorage) |
| `lib/clientOpId.ts` | Per-mutation UUID (`crypto.randomUUID` per call) |
| `print/shell/InfrastructureBanner.tsx` | Global "infrastructure error" status strip |

**Backend** (`packages/HimalayaUI/src/`):

| File | Role |
|---|---|
| `idempotency.jl` | `with_idempotency`, `OP_LOCKS`, `gc_idempotent_responses!` |
| `events.jl` | `apply_event!` (both methods), `broadcast_event!`, post-commit broadcast queue |
| `actions.jl` | `get_client_op_id`, `get_client_id` request-header extractors |
| `routes_*.jl` | Per-resource handlers; mutating ones wrap in `with_idempotency` |
| `pipeline.jl` | `analyze_exposure!` with `trace_known_unchanged` fast-skip |

---

## Tests

Frontend (`packages/HimalayaUI/frontend/test/queue/`):

- `useQueueMutation.test.tsx` — wrapper integration; per-call op-id mint, AbortSignal threading, error-class routing
- `replayCoordinator.test.ts` — own-op confirmation, foreign-event rollback/apply/rerun, MutationCache iteration order
- `cache-shape.test.ts` — `onSuccess` cache row shape per kind
- `sseEventPayload.contract.test.ts` — `applyRemoteToCache` wire-format merge per kind
- `rollbackSymmetry.test.ts` — `onMutate` snapshot ↔ rollback inverse
- `authHeaders.test.ts` — `X-Username` / `X-Client-Id` / `X-Client-Op-Id` header propagation
- `persistence.test.ts` — rehydrate cycle; schema-version drop
- `deferred.test.ts` — registry mint / resolve / leak-on-abort
- `mutatorOnSseWins.test.ts` — SSE-first race resolution per kind
- `hooks.test.tsx` — `useExposureHasPendingPeakOps` membership semantics
- `errors.test.ts` — validation vs infrastructure classification

Backend (`packages/HimalayaUI/test/`):

- `test_idempotency.jl` — `with_idempotency` cache hit / miss / concurrent retry / failure-not-cached / multi-event route
- `test_idempotency_replay_invariant.jl` — SSE-fanout invariant under retry; `_capture_sse_during` test harness
- `test_route_response_shapes.jl` — route emit shapes (the layer-1 half of the contract)
- `test_fast_skip.jl` — `analyze_exposure!` fast-skip + broadcast suppression
- `test_routes_*.jl` — per-resource happy/sad paths

---

## Further reading

- [event-log.md](event-log.md) — the dispatcher contract, hash invariants, and SSE multiplayer semantics. The queue composes with these unchanged.
- [contract-testing.md](contract-testing.md) — the six-layer testing rule.
- [superpowers/specs/2026-05-02-mutation-queue-design.md](superpowers/specs/2026-05-02-mutation-queue-design.md) — original design spec; 14 architectural decisions, fallback triggers.
- [superpowers/plans/2026-05-02-mutation-queue.md](superpowers/plans/2026-05-02-mutation-queue.md) — implementation plan; useful for archaeology, not as live reference.
