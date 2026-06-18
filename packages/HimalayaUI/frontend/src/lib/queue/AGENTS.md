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
├── hooks.ts              # useExposureHasPendingPeakOps + per-op gating hooks
├── queueMeta.ts          # per-mutation metadata
├── replacePlaceholder.ts # swap optimistic placeholder ids for real ids
├── assignmentDropped.ts  # dropped-assignment handling
├── optimisticId.ts       # Negative placeholder ids
├── peakQTol.ts           # q-value tolerance helper
├── deferred.ts           # Deferred promise utility
├── errors.ts             # Queue error taxonomy
├── testHelpers.ts        # test-only helpers
├── index.ts              # barrel
└── mutators/             # Per-kind implementations
    ├── peakAdd.ts
    ├── peakRemove.ts
    ├── peakSetExcluded.ts
    ├── reanalyzeExposure.ts
    ├── createSpeculative.ts
    ├── customIndex.ts
    ├── indexGroup.ts
    ├── assignment.ts
    ├── scopeSeries.ts
    ├── createSeries.ts
    ├── saveSeries.ts
    ├── deleteSeries.ts
    ├── commitSeriesPlate.ts
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

## ANTI-PATTERNS
- Never capture `client_op_id` in a closure at hook mount time
- Never read `peaks(id)` without gating on `useExposureHasPendingPeakOps` during in-flight ops
- Fixing one layer of a queue bug requires regression tests at all six layers (see `docs/contract-testing.md`)
