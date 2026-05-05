# Multi-layer reconciliation contract testing

This is the testing discipline for HimalayaUI's optimistic-and-confirmation
flow — the path a mutation takes from user gesture to durable SQL row to
SSE-broadcast cache merge across all connected tabs.

The contract is wide: a single mutation kind (e.g. `peak_added`) touches
six independent layers of code. A bug class can manifest at any of them.
The rule that follows from this is simple but easy to forget when fixing
one layer in isolation.

> **Rule: when fixing a bug at one layer, the regression test must
> exercise EVERY layer where the bug class can occur.**

---

## The six layers

```
1. Route emit (Julia)            packages/HimalayaUI/src/routes_*.jl
                                   Issues HTTP response + emits user_actions row
                                   + queues SSE broadcast.
2. SSE frame payload             packages/HimalayaUI/src/events.jl::broadcast_event!
                                   Serializes the event to the wire format —
                                   { kind, entity_id, client_id, client_op_id,
                                     payload, post_state? }.
3. applyRemoteToCache merge      lib/queue/applyRemoteToCache.ts
                                   Foreign-event path: turn an SSE frame back
                                   into a cache-row mutation per kind.
4. Cache row shape               (the resulting Peak / Index / Exposure object)
                                   What downstream consumers (PeakRow,
                                   PhasePanel, MentionChip) read.
5. Mutator onMutate              lib/queue/mutators/*.ts (optimistic write)
                                   Builds the optimistic placeholder row +
                                   returns the rollback closure.
6. Mutator onSuccess             lib/queue/mutators/*.ts (canonical write)
                                   Replaces the placeholder with the real row
                                   from the HTTP response or SSE-synthesized
                                   response.
```

Three review passes during the queue migration (PRs #17, #18, #24) all
surfaced fixes that landed at one layer but not its symmetric layer:

- A mutator's `onMutate` was fixed to scope a row by `source` (auto vs
  manual), but `applyRemoteToCache` for the same kind still matched
  unscoped — so a foreign tab's edit could flip the wrong row.
- A mutator's `onSuccess` was fixed to preserve `intensity` /
  `prominence` / `sharpness` on a `peak_excluded` event, but the SSE
  synthesis didn't carry those fields — so an SSE-first race blanked
  them.
- A mutator's optimistic id was fixed to be negative, but `applyRemoteToCache`
  still used the event id (positive) for foreign inserts — so cache
  rows from foreign tabs collided with auto-detected peaks.

Encoding the rule prevents the next instance.

---

## What "exercise every layer" means in practice

For each kind, three paired tests cover the six layers:

| Test file | Layers exercised |
|---|---|
| `cache-shape.test.ts` | Mutator `onSuccess` row shape (5, 6) |
| `sseEventPayload.contract.test.ts` | `applyRemoteToCache` wire-format merge (2, 3, 4) |
| `rollbackSymmetry.test.ts` | `onMutate` snapshot ↔ rollback inverse (5) |

Plus an authentication header test:

| Test file | Layers exercised |
|---|---|
| `authHeaders.test.ts` | `X-Username` / `X-Client-Id` / `X-Client-Op-Id` propagation across mutators |

And a backend pair:

| Test file | Layers exercised |
|---|---|
| `test_route_response_shapes.jl` | Route emit (1) — what the JSON response actually looks like |
| `test_idempotency_replay_invariant.jl` | SSE-fanout invariant under retry — the layer-1↔layer-2 relationship survives idempotent replay |

All located under `packages/HimalayaUI/frontend/test/queue/` (frontend) and
`packages/HimalayaUI/test/` (backend).

---

## Adding a new kind

When you add a new `OpKind` (a new mutator):

1. **Layer 1 (route emit)** — add a row to `test_route_response_shapes.jl`
   asserting the route's JSON response shape.
2. **Layer 2-4 (SSE merge → cache row)** — add a row to
   `sseEventPayload.contract.test.ts`. The test feeds a synthetic SSE
   frame through `applyRemoteToCache` and asserts the resulting cache
   shape.
3. **Layer 5-6 (optimistic + canonical mutator)** — add rows to
   `cache-shape.test.ts` (asserts `onSuccess` writes the canonical row)
   and `rollbackSymmetry.test.ts` (asserts `onMutate` followed by
   `restore` is a no-op on the cache).

Skipping any of these means a future symmetry bug in your kind will go
unnoticed until it hits production.

---

## Fixing an existing bug

When you fix a bug at one layer:

1. Identify which layer the fix is on — usually 3, 5, or 6.
2. For each *other* layer where the same bug class could occur, ask: does
   the existing regression test cover it? If not, add the row.
3. Land the fix and the test rows in the same PR.

The pattern in PR descriptions: "Fixed at layer X; added regression rows
to layers Y, Z because the same class can manifest there."

---

## Why this is fiddly

The six layers are *independently authored*. The route emit is hand-rolled
JSON; the SSE merge is a hand-rolled per-kind switch; the mutators are
hand-rolled per kind; the cache row is a TypeScript type. Nothing
generates one from another. A schema would make this trivial — and is on
the future-feature list — but until then, the only correctness mechanism
is the test discipline.

The good news: each test row is short (10-20 lines). Adding a complete
six-layer test set for a new kind is ~5 minutes of work. The discipline
is cheap. What's expensive is debugging a production multiplayer flicker
back to a 4-layer fix that should have been a 6-layer fix.

---

## Test layout reference

```
packages/HimalayaUI/frontend/test/queue/
  authHeaders.test.ts                    # headers per kind
  cache-shape.test.ts                    # mutator onSuccess (layer 5-6)
  rollbackSymmetry.test.ts               # mutator onMutate ↔ restore (layer 5)
  sseEventPayload.contract.test.ts       # applyRemoteToCache (layer 2-4)
  mutatorOnSseWins.test.ts               # SSE-first race per kind
  hooks.test.tsx                         # useExposureHasPendingPeakOps etc.
  useQueueMutation.test.tsx              # framework wrapper
  replayCoordinator.test.ts              # foreign-event replay
  deferred.test.ts                       # registry mint / clear
  persistence.test.ts                    # sessionStorage rehydrate
  errors.test.ts                         # validation vs infrastructure
  helpers.ts                             # shared test fixtures

packages/HimalayaUI/test/
  test_route_response_shapes.jl          # route emit (layer 1)
  test_idempotency_replay_invariant.jl   # SSE invariant under retry
  test_routes_*.jl                       # per-resource happy/sad paths
```
