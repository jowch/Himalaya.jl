# packages/HimalayaUI/frontend/src/lib/queue — Mutation Queue

## OVERVIEW
Optimistic-and-confirmation mutation pipeline. User gesture → durable SQLite row → SSE broadcast → cache merge across tabs.

## STRUCTURE
```
lib/queue/
├── types.ts              # Mutator<TInput, TScope, TResponse> + OpKind
├── Mutator.ts            # Base mutator factory
├── useQueueMutation.ts   # TanStack Query wiring + idempotency
├── applyRemoteToCache.ts # SSE → cache merge per kind
├── MutationCache.ts      # In-memory pending mutation store
├── replayCoordinator.ts  # Replay-as-rerun logic
├── clientOpId.ts         # Idempotency key minting
├── clientId.ts           # Per-tab SSE self-echo token
├── optimisticId.ts       # Negative placeholder ids
├── peakQTol.ts           # q-value tolerance helper
├── deferred.ts           # Deferred promise utility
├── errors.ts             # Queue error taxonomy
└── mutators/             # Per-kind implementations
    ├── peakAdd.ts
    ├── peakExclude.ts
    ├── peakUnexclude.ts
    ├── peakRemove.ts
    ├── indexConfirm.ts
    └── indexUnconfirm.ts
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
