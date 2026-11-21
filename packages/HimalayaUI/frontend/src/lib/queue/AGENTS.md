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
    ├── grouping.ts
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
- **The queue does NOT serialize HTTP.** TanStack `useMutation` fires requests concurrently; there is no `concurrency: 1`. "FIFO ordering" means only the replay-rerun *call order* (the `onMutate` re-invocation sequence after a foreign SSE event), never HTTP dispatch order. Concurrent in-flight requests for the same entity are possible, so correctness rests on idempotency + last-write-wins, not on request ordering.

## ANTI-PATTERNS
- Never capture `client_op_id` in a closure at hook mount time
- Fixing one layer of a queue bug requires regression tests at all six layers (see `docs/contract-testing.md`)
