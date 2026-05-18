# Corpus sample query layer (`useCorpusSamples`) — design

**Issue:** #156 (I1.2) · **Date:** 2026-05-17 · **Milestone:** HimalayaUI workflow-driven redesign

## Summary

Add the frontend corpus sample query layer: a `useCorpusSamples()` TanStack
Query hook backed by the corpus-wide `GET /api/samples` route, with a distinct
corpus query-key namespace (`["corpus", "samples"]`) kept separate from the
existing per-experiment `useSamples(experimentId)` (`queryKeys.samples`).

This is a pure data-layer addition. It backs the Phase 1 corpus surfaces —
the contact sheet (#160) and sample loupe (#161) — and is the query whose
distinct key the corpus `add_tag` work (#159) patches and invalidates.

## Rationale

The redesigned sample table is corpus-wide. `useSamples(experimentId)` exists
today but is keyed per experiment (`queryKeys.samples` → `["experiment", id,
"samples"]`), so its cache entries are experiment-scoped. A corpus surface
needs its own query and its own key namespace so cache invalidation,
optimistic patches, and SSE-driven updates target the right entry. Giving the
corpus query a distinct key now keeps the corpus `add_tag` work (#159) from
having to disambiguate cache keys after the fact.

## Verified context

The corpus `GET /api/samples` route already exists (merged as #154 / I0.1,
`routes_samples.jl:27`). Its projection returns each sample with its `tags`
**plus an extra `q_units` field** sourced from the owning experiment's config
— a field the per-experiment route and the per-experiment `Sample` type do
not carry. The redesign master plan is explicit that this field is
deliberate:

> "The `GET /api/samples` projection **includes `q_units` per sample** (from
> the owning experiment's config) — so the corpus surfaces and the future
> cross-experiment normalization have it without a later route-body change."
> — `docs/superpowers/plans/2026-05-17-himalaya-ui-redesign.md` §3 (line 95)

Phase 3 (§295) names "`q_units` normalization + caption" as a line item that
consumes this field. The corpus hook must therefore preserve `q_units` in its
return type — discarding it would force a later type change.

The route also accepts an optional `?experiment_id=` filter. Per issue #156
scope (key shape fixed at `["corpus", "samples"]`) and YAGNI — no Phase 1
consumer needs narrowing — the hook does **not** expose this filter. It can
be added when a consumer needs it.

## Design

Four edits, all append-only, in files the issue map (§3) flags as
low-contention for I1.2.

### 1. `api.ts` — type + fetcher

A `CorpusSample` interface extends the existing `Sample` with the corpus-only
`q_units` field. The per-experiment `Sample` interface is untouched, so
`useSamples` and its tests are unaffected.

```ts
// Corpus samples carry q_units (resolved from the owning experiment's
// config) — the per-experiment Sample does not. Phase 3 normalization
// reads this field.
export interface CorpusSample extends Sample {
  q_units: string;
}

export const listCorpusSamples = (): Promise<CorpusSample[]> =>
  request<CorpusSample[]>("GET", "/api/samples");
```

### 2. `queries.ts` — query key

A new `["corpus", ...]` namespace in the `queryKeys` object, distinct from
the per-experiment `queryKeys.samples`. The route takes no filter, so the key
is a plain `as const` tuple, not a function.

```ts
// Corpus-wide sample listing (redesign Phase 1). Distinct ["corpus", ...]
// namespace so corpus add_tag (#159) patches/invalidations never collide
// with the per-experiment ["experiment", id, "samples"] entries.
corpusSamples: ["corpus", "samples"] as const,
```

### 3. `queries.ts` — hook

```ts
/**
 * Corpus-wide sample list — every sample across all experiments, each with
 * its tags and a q_units from the owning experiment's config. Backs the
 * Phase 1 contact sheet (#160) and sample loupe (#161). Distinct from the
 * per-experiment useSamples(experimentId).
 */
export function useCorpusSamples() {
  return useQuery({
    queryKey: queryKeys.corpusSamples,
    queryFn: () => api.listCorpusSamples(),
  });
}
```

No `enabled` gate: the hook takes no parameter, so unlike its
`useSamples`/`useExposures` siblings there is nothing to gate the fetch on.

### 4. `test/queries-corpus-samples.test.tsx` — new Vitest file

Follows the `queries-samples.test.tsx` harness (`makeClient`, a `mockOnce`
helper, `renderHook` with a `QueryClientProvider` wrapper):

- **Hook fetches the route:** mock a `GET /api/samples` response, assert
  `result.current.data` resolves to the corpus list and that each row's
  `q_units` field survives into the hook's data.
- **Key shape:** assert `queryKeys.corpusSamples` deep-equals
  `["corpus", "samples"]`, and that it is not equal to `queryKeys.samples(id)`
  for any `id` — proving the corpus and per-experiment namespaces cannot
  prefix-collide.

## Out of scope

- The contact-sheet UI that consumes the hook (#160).
- The sample loupe (#161).
- Corpus `add_tag` / `remove_tag` mutations and their tri-scope cache
  patching (#159).
- Any change to the existing `useSamples(experimentId)` hook or its tests.

## Acceptance criteria

- [ ] `useCorpusSamples()` fetches `GET /api/samples` and returns the corpus
      sample list, with `q_units` preserved on each row.
- [ ] The corpus query uses the key `["corpus", "samples"]`, distinct from
      the per-experiment `queryKeys.samples`.
- [ ] The existing `useSamples(experimentId)` hook and its tests are
      unchanged.
- [ ] Vitest covers the hook and the key shape.
- [ ] `npm run build` (tsc + vite) passes.

## References

- Master plan §3 / §4.2: `docs/superpowers/plans/2026-05-17-himalaya-ui-redesign.md`
- Issue map I1.2: `docs/superpowers/plans/2026-05-17-himalaya-ui-redesign-issue-map.md`
- Corpus route: `packages/HimalayaUI/src/routes_samples.jl:27`
- Existing per-experiment hook: `packages/HimalayaUI/frontend/src/queries.ts:119`
