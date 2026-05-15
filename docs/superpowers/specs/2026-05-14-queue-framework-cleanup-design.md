# Queue framework cleanups (issue #129)

**Date:** 2026-05-14
**Issue:** [#129](https://github.com/jowch/Himalaya.jl/issues/129)
**Scope:** Frontend `packages/HimalayaUI/frontend/src/lib/queue/`

## Motivation

The mutation queue grew to 11 mutators. Three framework-shape leaks are now visible across multiple concrete examples:

1. **Per-kind synth lives in `replayCoordinator.ts`, not in the mutator.** The shape `synthesizeResponseFromSse` produces for each event kind has to satisfy the cache contract owned by the *corresponding mutator*, but those contracts live in different files. Every mutator with a kind-aware synth branch (saveComparison, createSpeculative, indexGroup) duplicates a "looksFull" detector in `onSuccess` to compensate for synth-vs-HTTP shape divergence.
2. **Five copies of the placeholder-replace loop.** `peakAdd.ts:58`, `trivial.ts:143/219/289/347` all open-code the "filter `id<0` placeholders matching field X, dedupe against server id, append" pattern with subtly different match keys.
3. **Framework metadata noise.** Four mutator `onSuccess` callbacks destructure `event_id` / `view_row_id` / `analysis_inputs_hash` off the response. The framework owns these; the mutator shouldn't have to spell them out (and silently drift if a new field is added).

## Non-goals

- No change to `applyRemoteToCache.ts` (foreign-event path; separate concern).
- No change to the rollback ordering / replay-as-rerun invariants.
- No change to the existing per-mutator HTTP contract (request shape, status codes).
- No rename of `OpKind` or any wire-side event kind.
- The `looksFull` invalidate detectors in `saveComparison` / `createSpeculative` / `indexGroup` stay. They handle synth-vs-HTTP shape divergence — a related but distinct concern from where the synth lives. Removing them would require a separate refactor that gives the synth full Comparison/IndexEntry shapes (with `members[]`, `phase`, lineage fields). That's deferred — this refactor only moves the synth *location*, not its *fidelity*.

## Design

### 1. `synthesizeFromSse` on the `Mutator` interface

Add an optional method to `Mutator<TInput, TScope, TResponse>` in `types.ts`:

```ts
synthesizeFromSse?: (remote: SseEvent, base: QueueResponseMeta) => TResponse | undefined;
```

Where `QueueResponseMeta` is a new exported type representing the framework-owned fields the coordinator already builds:

```ts
export interface QueueResponseMeta {
  event_id: number;
  client_op_id: string | null | undefined;
  analysis_inputs_hash: string | undefined;
}
```

`replayCoordinator.synthesizeResponseFromSse` becomes (~10 lines):

```ts
function synthesizeResponseFromSse(remote: SseEvent): unknown {
  const base: QueueResponseMeta = {
    event_id: remote.id,
    client_op_id: remote.client_op_id,
    analysis_inputs_hash: remote.post_state?.analysis_inputs_hash,
  };
  const mutator = resolveMutatorForEvent(remote.kind, remote.entity_type);
  const synth = mutator?.synthesizeFromSse?.(remote, base);
  return synth ?? { ...base, ...((remote.payload as object) ?? {}) };
}
```

#### Event-kind → mutator lookup

The existing `mutatorRegistry.ts` exports `resolveMutator(op)` keyed by `OpKind` + payload-shape discrimination. SSE frames carry the *event kind* string (not `OpKind`) — most are 1:1 but a few diverge (`comparison_save` → `comparison_created` | `comparison_submitted`, `delete_index` → `speculative_deleted`, `reanalyze_exposure` → `analyze_run`).

Add a sibling `resolveMutatorForEvent(eventKind, entityType): Mutator | undefined` next to `resolveMutator`. The dispatcher is an explicit switch (mirroring `resolveMutator`'s shape), not a derived registry — the `add_tag`/`remove_tag`/`post_message` entity-type splits couldn't be expressed by a per-mutator `eventKinds: string[]` field anyway. A consistency cross-check test ensures every mutator round-trips: for each canonical (eventKind, entityType) the mutator emits, `resolveMutatorForEvent` returns the same mutator.

#### Per-mutator migration

| Event kind(s) | Mutator | Action |
|---|---|---|
| `peak_added` | `peakAddMutator` | Move synth body verbatim to `synthesizeFromSse` |
| `add_tag` | `addSampleTagMutator` + `addExposureTagMutator` | Both register `add_tag`; coordinator must route by entity. → Use entity_type discriminator on the synth method, or split synth to two mutator methods and look up by `(kind, entity_type)` pair. Recommendation: **look up by `(kind, entity_type)` tuple** for these two; both other mutators are entity-uniqued already. |
| `comparison_created`, `comparison_submitted` | `saveComparisonMutator` | Move synth body verbatim; both event kinds route to the same mutator via the `resolveMutatorForEvent` switch (two `case` arms fall through) |
| `comparison_deleted` | `deleteComparisonMutator` | Move synth body verbatim |
| `peak_excluded`, `peak_unexcluded` | `peakSetExcludedMutator` | Move synth body verbatim; both event kinds route via the switch (the `excluded: boolean` flips per kind inside the synth) |

After migration, `replayCoordinator.ts` shrinks by ~70 lines and the kind-aware shape contracts live next to their consumers.

### 2. `replacePlaceholder<T>(list, serverItem, matches)`

New file `lib/queue/replacePlaceholder.ts`:

```ts
export function replacePlaceholder<T extends { id: number }>(
  list: T[],
  serverItem: T,
  matches: (item: T) => boolean,
): T[] {
  let replaced = false;
  const seen = new Set<number>();
  const out: T[] = [];
  for (const item of list) {
    if (item.id === serverItem.id) continue;        // dedupe vs concurrent SSE insertion
    if (!replaced && item.id < 0 && matches(item)) {
      out.push(serverItem);
      replaced = true;
      seen.add(serverItem.id);
      continue;
    }
    if (item.id >= 0) {
      if (seen.has(item.id)) continue;
      seen.add(item.id);
    }
    out.push(item);
  }
  if (!replaced) out.push(serverItem);
  return out;
}
```

#### Site migrations

| File:line | Old loop | New call |
|---|---|---|
| `peakAdd.ts:58` | peak list (flat) | `replacePlaceholder(peaks, response, pk => Math.abs(pk.q - response.q) < peakQTol(response.q) && pk.exposure_id === response.exposure_id)` |
| `trivial.ts:143` | sample.tags (nested) | thin wrapper: spread sample, replace `tags` |
| `trivial.ts:219` | exposure.tags (nested) | same pattern as sample.tags |
| `trivial.ts:289` | sample messages (flat) | `replacePlaceholder(msgs, response, m => m.body === response.body && m.sample_id === response.sample_id)` |
| `trivial.ts:347` | comparison messages (flat) | `replacePlaceholder(msgs, response, m => m.body === response.body && m.comparison_id === response.comparison_id)` |

Tag sites (143, 219) need a parent-record spread because the list lives inside `sample.tags` / `exposure.tags`. The helper still applies to the inner list — just call it inside the parent's setQueryData callback.

### 3. `stripQueueMetadata<T>(response)` *(bonus)*

New file `lib/queue/queueMeta.ts`:

```ts
export interface QueueResponseMeta { ... }  // re-export from types

export function stripQueueMetadata<T extends Partial<QueueResponseMeta>>(
  response: T,
): { meta: QueueResponseMeta; payload: Omit<T, keyof QueueResponseMeta> } {
  const { event_id, view_row_id, analysis_inputs_hash, client_op_id, ...payload } = response as any;
  return {
    meta: { event_id, view_row_id, analysis_inputs_hash, client_op_id },
    payload: payload as Omit<T, keyof QueueResponseMeta>,
  };
}
```

(Note: `view_row_id` is mutator-side metadata too — surveyed in 4 mutators. Add to the meta type.)

#### Site migrations

| File | Old destructure | New form |
|---|---|---|
| `peakAdd.ts:55` | `const { event_id: _e, view_row_id: _v, analysis_inputs_hash, ...rest } = response` | `const { meta, payload } = stripQueueMetadata(response); writeExposureHash(qc, exposureId, meta.analysis_inputs_hash); ...use payload` |
| `peakSetExcluded.ts:59` | same shape | same |
| `indexGroup.ts:72` | same shape | same |
| `indexGroup.ts:120` | same shape | same |

## Test strategy

### Six-layer contract tests stay green

The existing acceptance contract tests are the regression net:

- `rollbackSymmetry.test.ts` — every mutator's `restore` perfectly reverts optimistic writes (unchanged behavior)
- `cache-shape.test.ts` — `onSuccess` cache writes match `queryKeys.*` value types (no leaked queue metadata fields)
- `sseEventPayload.contract.test.ts` — synth response shape matches mutator expectations (moves to per-mutator assertions)
- `authHeaders.test.ts` — `buildAuthOpts` plumbing per mutator (unchanged)

### New tests

- `replacePlaceholder.test.ts` — unit tests covering: placeholder match, no-match append, concurrent SSE dedupe, multiple placeholders (only first replaced), empty list, list with only positive ids.
- `stripQueueMetadata.test.ts` — round-trip preserves payload, meta extraction, missing fields = undefined.
- `sseEventPayload.contract.test.ts` — add assertion: for every event kind currently in the old replayCoordinator switch, the corresponding `mutator.synthesizeFromSse` exists and returns the same shape.
- `mutatorRegistry.test.ts` — assert event-kind→mutator lookup covers every kind emitted by the backend (cross-reference with `events.jl`).

### Manual smoke (one of)

Run `peak_add` end-to-end via Playwright live spec (`e2e/live/peak-add-no-stale-banner.spec.ts`) — exercises the SSE-wins path that this refactor touches most heavily.

## Commit sequence

1. Add `replacePlaceholder` helper + tests (no migration)
2. Migrate 5 sites to `replacePlaceholder` (one commit per logical file: peakAdd, trivial)
3. Add `synthesizeFromSse` to Mutator interface, `resolveMutatorForEvent(kind, entity_type)` dispatcher, refactor replayCoordinator
4. Migrate 5 event kinds to own synth (one commit per mutator: peakAdd, addTag, saveComparison, deleteComparison, peakSetExcluded)
5. Add `stripQueueMetadata` helper + tests
6. Migrate 4 onSuccess consumers

Each step keeps all existing tests green. Step 3 is the riskiest — it changes the dispatch surface — so its commit includes both the new code path and the existing switch fallback until step 4 lands (incremental migration, not a flag-day).

## Risks

- **add_tag has two mutators sharing one event kind.** The dispatcher needs entity-type discrimination. Mitigation: lookup key is `(event_kind, entity_type)` tuple; both tag mutators declare different entity types.
- **`peak_added` synth response includes `exposure_id: remote.entity_id`** — a mutator owns this mapping today implicitly via the coordinator. Moving to `synthesizeFromSse` keeps the mapping but makes it explicit at the mutator site. No behavioral change.
- **Forward-scaffolded kinds (`post_message`, `set_exposure_status`, `add_tag`, `remove_tag`, `update_sample`, `select_exposure`) are unreachable today.** The replayCoordinator block-comment documents this. They stay unreachable — the refactor doesn't activate them; it just lets future activation be a one-line `synthesizeFromSse` addition per mutator.

## Acceptance criteria (from #129)

- [x] `replayCoordinator.synthesizeResponseFromSse` no longer dispatches on a per-kind switch (reduced to fallback lookup + default).
- [x] Mutators with non-trivial synth own their synth as a mutator method.
- [x] Five placeholder-replace loops collapse to one helper call each.
- [x] All queue contract tests green; rollbackSymmetry / cache-shape / sseEventPayload / authHeaders unchanged.
