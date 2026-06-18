---
name: new-mutator
description: Scaffold a new queue-shaped mutation in HimalayaUI — Mutator definition, queries.ts hook via useQueueMutation, queries-*.test.tsx round-trip test, optional applyRemoteToCache branch, and (if needed) backend route wrapping in with_idempotency. Reads CLAUDE.md, docs/event-log.md §3a, and existing mutators as live reference. Usage: /new-mutator <kind> [--scope <scope>] [--no-cache-merge]
disable-model-invocation: true
---

# new-mutator

Scaffolds the queue-shaped pattern for a new mutation in HimalayaUI. The mutation queue framework (Plan 8) ships HTTP↔SSE race resolution, optimistic effects, replay-as-rerun, and idempotent retries — but every new mutator has to honor 5+ cross-cutting invariants to compose correctly. This skill walks each step against current source instead of relying on a static template.

## Args

- `<kind>` — the new OpKind in snake_case (e.g. `peak_relabeled`, `index_starred`). Must be added to the `OpKind` union in `lib/queue/types.ts`.
- `--scope <scope>` — what closure-injected fields the mutator needs. Common scopes: `exposure` (`{exposureId, username, clientId}`), `sample` (`{experimentId, sampleId, ...}`), `exposure-tag` (`{sampleId, exposureId, ...}`).
- `--no-cache-merge` — skip the `applyRemoteToCache` branch. Use when the mutation's effect is fully captured by `mutator.onSuccess` writing the response (no cross-tab cache merge needed) OR the default invalidate fallback is acceptable. Default: include a per-kind branch.
- `--readonly-route` — backend route is GET-only, no `with_idempotency` needed. Skips the route-wrapping step.

## Procedure

### 1. Read the current reference implementations

Before writing any code, read these to understand the *current* shape:

```
CLAUDE.md                                                   ← §"HimalayaUI gotchas" + Plan 8 entries
docs/event-log.md                                           ← §3a "Mutation queue + idempotency"
packages/HimalayaUI/frontend/src/lib/queue/types.ts         ← OpKind union, Mutator interface
packages/HimalayaUI/frontend/src/lib/queue/useQueueMutation.ts  ← framework wiring
packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts  ← per-kind SSE merge
packages/HimalayaUI/frontend/src/lib/queue/mutators/peakAdd.ts   ← canonical optimistic-with-placeholder
packages/HimalayaUI/frontend/src/lib/queue/mutators/peakSetExcluded.ts  ← canonical kind-split (excluded vs unexcluded)
packages/HimalayaUI/frontend/src/lib/queue/mutators/indexGroup.ts       ← canonical scope-shared multi-mutator file
packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts          ← canonical 1:1 cache merge
packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts           ← rehydrate dispatch
packages/HimalayaUI/frontend/test/queries-peaks.test.tsx                ← canonical round-trip test
packages/HimalayaUI/src/routes_peaks.jl                                 ← with_idempotency + InTransaction + post-commit broadcast
packages/HimalayaUI/src/idempotency.jl                                  ← with_idempotency contract
```

If the new mutator has a placeholder-id (peak/message/tag insertions), `peakAdd.ts` is the template. If it's a state flip (excluded/select/status), `peakSetExcluded.ts` or `setExposureStatusMutator` in `trivial.ts`. If it bundles related ops on the same scope, `indexGroup.ts`.

### 2. Decide the kind name and OpKind addition

The kind is the SAME string as the backend `apply_event!` `kind=` argument. Past-tense verbs: `peak_added`, `index_confirmed`, `speculative_created`. Match existing conventions.

Add the kind to the `OpKind` union in `lib/queue/types.ts`. TypeScript will tell you everywhere it needs to be handled (`applyRemoteToCache` switch, `errors.ts` VALIDATION_COPY, `mutatorRegistry.ts`).

### 3. Design the payload + scope

**Scope** = closure-injected fields the hook captures from `useAppState`. Typical: `{exposureId, sampleId, username, clientId}`.

**Input** = per-`mutate()` arguments. E.g., `{q: number}` for peak_added, `{key, value}` for add_tag.

**Payload sent to server** = `{kind, clientOpId, ...scope, ...input}`. The backend route deconstructs whatever it needs.

The `Mutator<TInput, TScope, TResponse>` interface (3 generics, post-loose-end-#6) gives inference at the call site without `as unknown as` casts. Don't reintroduce the `Flat = OpPayload<TInput> & TScope & TInput` cast pattern — it was deleted in commit 67a18a3.

### 4. Implement the mutator

Create `packages/HimalayaUI/frontend/src/lib/queue/mutators/<name>.ts` (or add to an existing file if scope-shared, like `indexGroup.ts`).

Required fields:
- `kind: OpKind` — must match backend `apply_event!`
- `onMutate(payload, qc): RollbackContext` — write optimistic effect, snapshot prior cache, return `{restore}` closure
- `request(payload, signal): Promise<TResponse>` — call `api.<route>(...)` with `authOpts(payload.username, payload.clientId, payload.clientOpId)`
- `onSuccess(payload, response, qc)` — replace placeholders with server-assigned ids, write `analysis_inputs_hash` if peak-affecting
- `affectsExposurePeaks?: () => boolean` — `true` for peak/reanalyze ops; gates `useExposureHasPendingPeakOps`. Default false; only set to true if `useSpeculativeSnap` should suspend during this op.

**Cross-cutting invariants** (regression-tested; don't violate):
- **Optimistic placeholder ids are NEGATIVE.** Use `nextOptimisticId()` from `lib/optimisticId.ts` — shared monotonic counter across mutators, no collisions.
- **Q-tolerance comparisons** use `peakQTol(q)` from `lib/peakQTol.ts` (mirrors Julia's `MAX(1e-6, |q|*0.001)`). Never use bare `< 1e-6`.
- **Snapshot prior cache value before any setQueryData.** The restore closure must be self-contained.
- **Don't blanket-invalidate in onSuccess.** SSE post_state via applyRemoteToCache handles cross-tab. Only invalidate if SSE coverage is genuinely missing for this kind.

### 5. Wire the hook in `queries.ts`

```typescript
export function useFooMutation(scope1: number, scope2: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(fooMutator, {
    scope1, scope2, username, clientId: CLIENT_ID,
  });
}
```

TypeScript infers all 3 generics from the mutator. Drop any explicit `<TInput, TScope, TResponse>` annotations at the call site.

### 6. Add the `applyRemoteToCache` branch (unless `--no-cache-merge`)

In `lib/queue/applyRemoteToCache.ts`, add a `case "<kind>":` branch. The branch fires when an SSE frame arrives with this kind from another tab/user (`client_op_id` doesn't match a local pending op). Pattern: cache-merge using payload fields; if the response is too rich for cache merge, invalidate the affected query keys and let the SSE refetch settle them.

### 7. Update `mutatorRegistry.ts`

Add the kind → mutator mapping. For dual-scope kinds (sample-vs-exposure-scoped variants of the same kind, like `add_tag`), discriminate via `payload.experimentId !== undefined`.

### 8. (Backend) Wrap the route in `with_idempotency` (unless `--readonly-route`)

In `packages/HimalayaUI/src/routes_<resource>.jl`:

```julia
@post "/api/..." function(req::HTTP.Request, ...)
    db = current_db()
    return with_idempotency(db, req) do
        # body — INSERTs/UPDATEs that should commit atomically with the event
        result = apply_event!(InTransaction(), db, req;
            kind        = "<your-kind>",
            entity_type = "exposure",  # or sample/etc
            entity_id   = ...,
            payload     = Dict(:field => value))
        # response body — what gets cached for retries
        HTTP.Response(201, ["Content-Type" => "application/json"],
            JSON3.write(Dict(...)))
    end
end
```

**Critical**: use `apply_event!(InTransaction(), db, req; ...)` inside `with_idempotency`'s closure — NOT the default `apply_event!(db, req; ...)`. The default opens a nested savepoint and broadcasts immediately, bypassing the post-commit queue. The `InTransaction()` variant participates in the outer tx without broadcasting. If you need an enriched broadcast (with `post_state`), enqueue it via `_enqueue_post_commit_broadcast!` so it fires after commit.

If the new event affects the analysis hash, also call `analyze_exposure!(db, exp_id, dir; trace_known_unchanged=true, defer_broadcast=true)` synchronously inside the closure — closes the stale-banner window without a refetch round-trip.

### 9. Write the test

Create `packages/HimalayaUI/frontend/test/queries-<topic>.test.tsx` (or extend an existing topic-grouped file). Mirror `queries-peaks.test.tsx`. Required scenarios:

1. **Optimistic effect lands** — never-resolving fetch mock, assert cache mutates synchronously after `mutate(...)`.
2. **HTTP success replaces optimistic with server response** — mock with response body, wait for settle, assert real id replaces placeholder.
3. **HTTP error rolls back** — 4xx mock, assert cache returns to seeded state.
4. (For peak-affecting) **`analysis_inputs_hash` written to exposure cache** on success.

Use `withClient` + `mockOnce` from existing test files; `pendingDeferreds.clear()` in `beforeEach`.

### 10. Verify

```bash
cd packages/HimalayaUI/frontend
npm test 2>&1 | tail -10
npm run build 2>&1 | tail -3
```

Backend (if route changes):

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' 2>&1 | tail -3
```

All tests green; build clean (TS strict + `exactOptionalPropertyTypes`).

## Anti-patterns

- ❌ Using `qc.invalidateQueries` in `onSuccess` for a peak-affecting kind. SSE post_state covers this.
- ❌ Adding `as unknown as` casts back into mutator files. The 3-generic `Mutator` interface eliminates them.
- ❌ Calling `apply_event!(db, req; ...)` (default method) inside `with_idempotency`. Use `apply_event!(InTransaction(), db, req; ...)`.
- ❌ Broadcasting inside the with_idempotency closure. Defer via `_enqueue_post_commit_broadcast!` or use `defer_broadcast=true`.
- ❌ Using `< 1e-6` for q comparisons. Use `peakQTol(q)`.
- ❌ Module-local optimistic-id counter. Use `nextOptimisticId()` from the shared util.
- ❌ Setting `affectsExposurePeaks: () => true` for tag/message/status ops. Only peak ops + reanalyze affect the peak set.
