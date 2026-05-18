# Corpus sample-tag mutations through the queue (#159, I1.3)

**Status:** design approved · **Date:** 2026-05-18
**Issue:** #159 — Wire corpus sample-tag mutations through the queue (I1.3)
**Milestone:** HimalayaUI workflow-driven redesign (Phase 1)

## Problem

The redesign's corpus contact sheet (#160) is the first reader of a
corpus-wide `samples` query — `useCorpusSamples()` / `queryKeys.corpusSamples`
(`["corpus","samples"]`), merged in #156. Tagging a sample from that sheet
must round-trip through the mutation queue: an optimistic patch, an HTTP
write, and a foreign-tab SSE refresh.

Two layers of the queue currently assume tags are *binary*-scoped —
sample-vs-exposure:

1. **`resolveMutator`** (`mutatorRegistry.ts`) discriminates a sample tag
   from an exposure tag on `experimentId !== undefined`. A corpus sample-tag
   op carries `{sampleId}` only — no `experimentId` — so it would misroute to
   the *exposure* mutator.
2. **`applyRemoteToCache`** (`applyRemoteToCache.ts`) maps an `add_tag` /
   `remove_tag` SSE frame to one of two cache keys — `samples(experiment_id)`
   or `exposures(sample_id)`. It never invalidates `corpusSamples`, so a
   foreign tab's tag would never refresh the contact sheet.

Both must learn the corpus scope. This is a full six-layer mutation-queue
change (route → SSE → `applyRemoteToCache` → cache row → `onMutate` →
`onSuccess`), not a one-line addition.

## Scope decision: separate corpus mutators

The corpus `add_tag` / `remove_tag` get **new, dedicated mutators**
(`addCorpusSampleTagMutator`, `removeCorpusSampleTagMutator`) that patch only
`queryKeys.corpusSamples`. The existing experiment-scoped `addSampleTagMutator`
/ `removeSampleTagMutator` are **left untouched**.

Rationale: a unified mutator patching every cached `samples` projection is
"more correct" only in a world where the corpus list and a per-experiment
list render simultaneously — which the Phase 1 cutover (#163, retires the
Inspect page) makes unlikely. Separate mutators are the lowest-risk choice:
no behaviour change to the existing path, existing contract tests stay green,
and `resolveMutator` becomes a clean true tri-way split. The accepted cost is
that tagging from the corpus sheet does not optimistically patch a
per-experiment `samples` list cache; an SSE frame (which always carries
`experiment_id`) still invalidates it for any tab that has it cached.

## Verified finding: `resolveMutatorForEvent` stays 2-arm

The issue text and master plan §11 say to convert `resolveMutator` **and**
`resolveMutatorForEvent` "to tri-scope consistently". `resolveMutatorForEvent`
**cannot** become a literal tri-scope split, and should not pretend to. This
was verified, not assumed:

- There is exactly one sample-tag route — `POST /api/samples/{id}/tags`
  (single) and `POST /api/samples/tags/batch` (batch). Both emit
  `apply_event!(... kind="add_tag", entity_type="sample" ...)` with payload
  `{key, value, tag_id, experiment_id}` (`routes_samples.jl:147`, `:270`).
- The backend has no notion of corpus-vs-experiment scope. "Corpus-ness" is a
  frontend-only construct: which hook ran and thus which scope was
  flat-spread into the queue op. The op reaches the server as a plain
  `POST /api/samples/{id}/tags` regardless.
- Therefore the SSE frame is byte-identical for a corpus-originated and an
  experiment-originated sample tag — same `kind`, same `entity_type="sample"`,
  same payload (the route always resolves `experiment_id` from the sample
  row). No field distinguishes them.
- `resolveMutatorForEvent(kind, entityType)` reads only `remote.kind` and
  `remote.entity_type`. Widening its signature to take the full payload would
  change nothing — the distinction is *absent* from the wire, not merely hard
  to read.
- Its only job (per `synthesizeResponseFromSse` in `replayCoordinator.ts`) is
  to pick a `synthesizeFromSse` for the own-op SSE-confirmation response
  shape. `addSampleTagMutator.synthesizeFromSse` yields a `SampleTag`; the
  corpus mutator's `onSuccess` consumes a `SampleTag`. The actual cache patch
  is done by the *pending mutation's own `onSuccess`* — the corpus mutator,
  correctly. `remove_tag` has no `synthesizeFromSse` on either side and falls
  through to the generic `{...base,...payload}`.

Conclusion: `resolveMutatorForEvent`'s `add_tag` / `remove_tag` cases stay
`entity_type === "sample" ? <sample mutator> : <exposure mutator>`, with a
comment recording this verification so a reviewer does not flag it as an
incomplete conversion. The corpus mutators **do not define
`synthesizeFromSse`** — it would be unreachable dead code.

## Design — the six layers

### Layer 1 — Route emit (`routes_samples.jl`)

`POST /api/samples/{id}/tags` accepts an optional `source` body field:

- Validate string-if-present (mirror the batch route #169:
  `if haskey(body, :source) && !(body.source isa AbstractString)` → 400).
- `source = haskey(body, :source) ? String(body.source) : "manual"`.
- `INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?,?,?,?)`
  with the resolved value; the 201 response emits `:source => source`.

The `add_tag` **event payload is unchanged** — `{key, value, tag_id,
experiment_id}`. Corpus tags are all `manual`; the route param exists for the
future series-scoping path, which will add `source` to the payload when it
needs a non-`manual` source on the wire. Leaving the payload alone means no
`sseEventPayload.contract.test.ts` churn.

### Layer 2 — `resolveMutator` (op path — genuine tri-scope)

`mutatorRegistry.ts`, `add_tag` and `remove_tag` cases. The discriminator
peels `exposureId` off **first** — an exposure op carries `{sampleId,
exposureId}` and, like a corpus op, has no `experimentId`, so the two collide
unless `exposureId` is tested before falling through to corpus:

```
exposureId   !== undefined → exposure mutator
experimentId !== undefined → experiment-sample mutator
else                       → corpus-sample mutator   ← new
```

The local payload cast in `resolveMutator` already declares `exposureId?:
number` — no type change needed, the discriminator just has to *use* it. The
ordering is the load-bearing correctness point — a comment must call it out.

`resolveMutatorForEvent` is **not** changed for `add_tag` / `remove_tag`
beyond a clarifying comment (see "Verified finding" above).

### Layer 3 — Corpus mutators + hooks

`lib/queue/mutators/trivial.ts` — two new mutators:

- `addCorpusSampleTagMutator` (`kind: "add_tag"`): scope `BaseScope &
  { sampleId: number }`, input `{ key, value }`, response `SampleTag`.
  - `onMutate`: read `queryKeys.corpusSamples` (`CorpusSample[]`), append an
    optimistic tag `{ id: nextOptimisticId(), key, value, source: "manual" }`
    to the matching sample's `tags`, return a `restore` closure.
  - `request`: `api.addSampleTag(p.sampleId, p.key, p.value, authOpts)`.
  - `onSuccess`: `replacePlaceholder` the optimistic tag with the server row
    (strip `sample_id`, keep `id/key/value/source`), matched on
    `key === p.key && value === p.value`.
  - **No `synthesizeFromSse`** — the own-op SSE-wins response is shaped by
    `addSampleTagMutator.synthesizeFromSse` via `resolveMutatorForEvent`.
- `removeCorpusSampleTagMutator` (`kind: "remove_tag"`): scope `BaseScope &
  { sampleId: number }`, input `{ tagId }`, response `void`.
  - `onMutate`: filter the tag out of the `corpusSamples` cache,
    return `restore`.
  - `request`: `api.removeSampleTag(p.sampleId, p.tagId, authOpts)`.
  - `onSuccess`: `() => {}`; `treats404AsSuccess: true`.

`queries.ts` — two new hooks, scope carries **no** `experimentId`:

```ts
useAddCorpusSampleTag(sampleId)    // { sampleId, username, clientId }
useRemoveCorpusSampleTag(sampleId) // { sampleId, username, clientId }
```

`useRemoveCorpusSampleTag` wraps `mutate` as `(tagId) => inner.mutate({ tagId })`,
mirroring `useRemoveSampleTag`.

`mutatorRegistry.ts` imports and wires the two new mutators into
`resolveMutator`. No new `api.ts` fetcher — the corpus mutators reuse
`api.addSampleTag` / `api.removeSampleTag` (same route).

### Layer 4 — `applyRemoteToCache` (foreign-event refresh)

The `add_tag` / `remove_tag` branch becomes:

```
entity_type === "sample"  → invalidate samples(payload.experiment_id)
                            AND corpusSamples
entity_type === "exposure"→ invalidate exposures(payload.sample_id)  (unchanged)
```

The route always emits `experiment_id`, so the experiment key is still
invalidated (harmless if not cached); the `corpusSamples` invalidation is the
new behaviour that refreshes the contact sheet from a foreign tab.

### Layers 5–6 — cache row shape & mutator writes

`CorpusSample extends Sample`, so `tags: SampleTag[]` is already the cache row
shape; the corpus mutators write the same `SampleTag` shape the experiment
mutators do. No new type.

## Testing — six-layer contract coverage

Per `docs/contract-testing.md`'s paired-file convention:

| Layer | Test file | Assertion |
|---|---|---|
| 1 Route emit | `test/test_routes_samples.jl` | `POST .../tags` with explicit `source`, with default `source`, non-string `source` → 400 |
| 1 Response shape | `test/test_route_response_shapes.jl` | 201 body carries the resolved `source` |
| 2 SSE payload | `frontend/test/queue/sseEventPayload.contract.test.ts` | `add_tag` sample frame still has `{key,value,tag_id,experiment_id}` (regression: payload unchanged) |
| 3 Foreign merge | `frontend/test/queue/sseEventPayload.contract.test.ts` (or sibling) | foreign `add_tag` / `remove_tag` `entity_type="sample"` frame invalidates **both** `corpusSamples` and `samples(experiment_id)` |
| — Resolver | `frontend/test/queue/mutatorRegistry.test.ts` | tri-way `resolveMutator` split: `{sampleId}`→corpus, `{experimentId,sampleId}`→experiment, `{sampleId,exposureId}`→exposure, for `add_tag` and `remove_tag` |
| 5 `onMutate` | `frontend/test/queue/rollbackSymmetry.test.ts` | corpus `onMutate` snapshot ↔ `restore` inverse on the `corpusSamples` cache |
| 6 `onSuccess` | `frontend/test/queue/cache-shape.test.ts` | corpus `onSuccess` replaces the optimistic placeholder with the canonical `SampleTag` |
| Auth | `frontend/test/queue/authHeaders.test.ts` | `X-Username` / `X-Client-Id` / `X-Client-Op-Id` propagate through the two new mutators |

## Acceptance criteria (from #159)

- [ ] Adding a tag to a corpus sample routes to a sample-tag mutator, not the
      exposure mutator.
- [ ] A corpus tag add/remove round-trips through the queue and patches
      `corpusSamples` optimistically.
- [ ] A corpus tag SSE frame from a foreign tab invalidates `corpusSamples`
      and refreshes the contact sheet.
- [ ] `resolveMutator` is a true tri-scope split; `resolveMutatorForEvent`
      carries the verified-2-arm comment.
- [ ] `applyRemoteToCache`'s `add_tag` / `remove_tag` branch invalidates the
      corpus key for sample-scoped frames.
- [ ] Tag `source` is parameterized at the single-tag route.
- [ ] Six-layer contract tests pass.

## Dependencies & sequencing

- **#156 (corpus query layer)** — merged. `queryKeys.corpusSamples`,
  `useCorpusSamples`, `CorpusSample` exist. Unblocked.
- **Blocks #163** (Phase 1 cutover).
- **Shared-file contention** with the series event-kind cluster #166 / #167 /
  #168 on `mutatorRegistry.ts` and `applyRemoteToCache.ts` (issue-map §3).
  Those issues *append* new `case` arms; this issue *modifies* the existing
  `add_tag` / `remove_tag` branch. Land-order must be coordinated and this
  issue rebased **last** onto whatever that cluster produced — the rebase is
  not purely mechanical, because the modified branch may have moved.
- `routes_samples.jl` — the single-tag route body edit (~line 136) is
  disjoint from the appended batch route (#169, merged); no live contention.

## Out of scope

- The contact-sheet UI itself (#160).
- The corpus query-key definition (#156 — merged).
- Series event kinds (#166–#168) — they touch the same switch files but are
  separate issues.
- Adding `source` to the `add_tag` event payload — deferred to the
  series-scoping path that actually needs a non-`manual` source on the wire.
- A unified mutator patching every `samples` projection — see "Scope
  decision" above.
