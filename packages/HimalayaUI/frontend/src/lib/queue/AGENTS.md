# packages/HimalayaUI/frontend/src/lib/queue — Mutation Queue

## OVERVIEW
Optimistic-and-confirmation mutation pipeline. User gesture → durable SQLite row → SSE broadcast → cache merge across tabs.

## STRUCTURE
```
lib/queue/
├── types.ts              # Mutator<TInput, TScope, TResponse> interface + OpKind
├── useQueueMutation.ts   # TanStack Query wiring + idempotency
├── applyRemoteToCache.ts # SSE → cache merge per kind
├── replayCoordinator.ts  # Replay-as-rerun logic
├── mutatorRegistry.ts    # resolveMutator / resolveMutatorForEvent
├── persistence.ts        # mirror MutationCache pending set to sessionStorage
├── queueMeta.ts          # per-mutation metadata
├── replacePlaceholder.ts # swap optimistic placeholder ids for real ids
├── assignmentDropped.ts  # dropped-assignment handling
├── sampleCall.ts         # writes a sample's call to every frame's assignment cache
├── reconnectRefresh.ts   # refetch after an SSE reconnect
├── optimisticId.ts       # Negative placeholder ids
├── peakQTol.ts           # q-value tolerance helper
├── deferred.ts           # Deferred promise utility
├── errors.ts             # Queue error taxonomy
├── testHelpers.ts        # test-only helpers
└── mutators/             # Per-kind implementations
    ├── peakAdd.ts
    ├── peakRemove.ts
    ├── peakSetExcluded.ts
    ├── reanalyzeExposure.ts
    ├── createSpeculative.ts
    ├── customIndex.ts
    ├── indexGroup.ts
    ├── assignment.ts
    ├── bridges.ts
    ├── labels.ts
    └── trivial.ts

# Note: clientOpId.ts (idempotency-key minting) and clientId.ts (per-tab SSE
# self-echo token) live one level up in lib/, not in lib/queue/.
```

## CONVENTIONS
- **Invariant**: mint `client_op_id` inside `mutationFn`, not at hook construction
- **Invariant**: rollback in REVERSE order; re-apply in INSERTION order
- **Invariant**: `MutationCache.getAll()` preserves insertion order (load-bearing for replay)
- **Invariant**: optimistic placeholder ids are negative (`peak.id < 0`)
- `OpKind` ≠ event kind in the Julia backend (mapping is 1-to-1 but names may differ)
- **The queue serializes HTTP.** Every queue mutation, live or replayed after a reload, shares the TanStack scope `"queue"` (`queueMutationOptions`), so one request is in flight at a time and edits land in the order made (#351). Replayed ops queue in stored order, ahead of edits made after the reload.

## ANTI-PATTERNS
- Never capture `client_op_id` in a closure at hook mount time
- Fixing one layer of a queue bug requires regression tests at all six layers (see `docs/contract-testing.md`)

- **Peak list settles by refetch, deferred.** `refreshPeaks.ts` invalidates an exposure's peaks only when no peak mutation for it is pending; an earlier refetch would erase later optimistic patches. Own `onSuccess` and `applyRemoteToCache` peak cases call it. Test: `test/queue/peakRefetchRace.test.tsx`.
- **Never `setQueriesData` with a pass-through updater.** `setQueryData` marks a query fresh and clears its invalidation even when the updater returns the same list, so patching every list "in case" undoes an earlier `invalidateQueries` on the lists it didn't change. That hid a new average after a drop until a reload. Patch only the lists that hold the row: `patchListsWith` in `selectionLists.ts`. Test: `test/applyRemoteToCache.ingest.test.ts`.
