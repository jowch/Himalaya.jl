# Corpus-wide `GET /api/samples` route — design

**Issue:** #154 (I0.1) — Wave A foundational, no upstream dependencies.
**Date:** 2026-05-17

## Goal

Add a corpus-wide `GET /api/samples` endpoint returning every sample across all
experiments in one projection, as a sibling to the existing experiment-scoped
`GET /api/experiments/{id}/samples`. The redesigned UI is corpus-scoped —
experiment becomes a filter facet, not a container — so its surfaces need a
single listing route spanning the whole corpus.

This route is the sole data source for the Phase 1 sample table (consumed via
the `useCorpusSamples` hook, issue #156) and, because it carries `q_units` per
sample, the data source for future cross-experiment trace normalization in the
series stage.

## Context: why the experiment-scoped route is *not* the model

`GET /api/experiments/{id}/samples` and its `useSamples` hook are consumed only
by surfaces the redesign retires (`InspectPage`, the Index three-card
workspace, dual-nav scaffolding). The master plan keeps the experiment-scoped
*route* alive ("untouched and green") purely as append-only discipline — not
because anything new calls it. After Phase 5 it has zero frontend consumers.

Consequence: "match the sibling's shape" is **not** a design goal — there is no
second live consumer to share a type with. `useCorpusSamples` is the only
sample-list hook going forward. The projection shape is therefore decided
solely by what the Phase 1 sample table needs.

## Scope

- New route `GET /api/samples` added to `register_samples_routes!()` in
  `packages/HimalayaUI/src/routes_samples.jl`.
- Optional `?experiment_id=` query parameter to filter to one experiment, so
  the corpus route subsumes the experiment-scoped use case.
- The projection includes `tags` and `q_units` per sample.
- A Julia route test covering the full-corpus response, the `?experiment_id=`
  filter, and the malformed-parameter path.
- Append-only change: no schema change, no event kinds, no migration.

## Out of scope

- The existing `GET /api/experiments/{id}/samples` route — left untouched, its
  tests stay green.
- Frontend consumption — the `useCorpusSamples` query hook is issue #156.
- Cross-experiment `q_units` normalization logic — this route only *exposes*
  `q_units`; normalization is a Phase 3 concern.
- Pagination — the corpus is small enough to return whole; revisit only if it
  grows.

## Design

### Route behavior

`GET /api/samples` is added to the existing `register_samples_routes!()`
function. No `server.jl` edit is needed: `register_samples_routes!()` is
already called at `server.jl:122`. (Issue #154's "registered in server.jl"
predates that being a single shared registration function.)

`GET /api/samples/{id}` is a distinct path — no route-matching conflict.

### Response shape

JSON array, one object per sample, ordered by `id`. Each object is
`row_to_json(sample)` plus two added fields:

- `tags` — array of `{id, key, value, source}`, the same nested-tag shape the
  sibling route produces.
- `q_units` — string sourced from the owning experiment's config; ASCII default
  `"A-1"` when the config is missing or the TOML is malformed.

**Why `tags` is included** (Option A, single-call): the Phase 1 sample table's
deliverable spec gives each row a Tags column ("light free-form chips"), and
corpus `add_tag`/`remove_tag` operate on that surface. The table's only
sample-list fetch is this route. Bundling `tags` lets `GET /api/samples` alone
fully feed the table — one request, no client-side join against the future
corpus `/api/sample-tags` route (#159).

### `?experiment_id=` filter

- Valid integer → `WHERE experiment_id = ?`. A nonexistent id yields `[]`
  naturally from SQL.
- Malformed value (non-numeric) → `400 {"error": "..."}`, rather than silently
  ignoring the parameter.
- Absent → whole corpus.

### Edge cases

- Empty corpus, or filtered-to-empty → `[]`.

### Query strategy — never N+1

Three queries total, regardless of corpus size:

1. `SELECT * FROM samples [WHERE experiment_id = ?] ORDER BY id`.
2. `SELECT id, config FROM experiments` → folded into a `Dict{Int,String}` of
   `experiment_id → q_units`.
3. `SELECT id, sample_id, key, value, source FROM sample_tags
   WHERE sample_id IN (?, …) ORDER BY id` → grouped into a
   `Dict{Int,Vector}` keyed by `sample_id`. **Skipped entirely** when query 1
   returns no samples — an empty `IN ()` is invalid SQL.

A single `map` over the samples then joins both dicts in memory:
`d[:tags] = get(tags_by_sample, id, [])`,
`d[:q_units] = get(qunits_by_exp, experiment_id, "A-1")`.

This deliberately avoids the sibling route's per-sample tag query (an N+1).
Since the corpus route is not imitating the sibling, it is free to do this
right.

### Shared helper — `q_units` extraction

The TOML-parse-with-`"A-1"`-fallback logic currently lives inline inside
`_experiment_row_to_json` (`routes_experiments.jl:8`). Extract it to a small
helper:

```julia
_q_units_from_config(cfg_text)::String
```

`_experiment_row_to_json` and the new route both call it. This is a targeted
dedup of code the new route would otherwise copy verbatim — not unrelated
refactoring. The helper lives in `routes_experiments.jl`, next to its existing
caller.

## Testing

Extend `packages/HimalayaUI/test/test_routes_samples.jl` with a second
`@testset` that builds **two** experiments with distinct `q_units` in their
configs, then asserts:

- `GET /api/samples` returns all samples across both experiments, each with the
  correct `q_units` and `tags`.
- `GET /api/samples?experiment_id=<id>` returns only that experiment's samples.
- A sample with tags shows them in the projection.
- A malformed `?experiment_id=` value returns `400`.
- The existing `GET /api/experiments/{id}/samples` testset is untouched and
  stays green.

## Acceptance criteria

- [ ] `GET /api/samples` returns every sample in the database.
- [ ] Each sample in the projection carries a `q_units` field sourced from its
      experiment's config, and a `tags` array.
- [ ] `GET /api/samples?experiment_id=<id>` returns only that experiment's
      samples; a malformed value returns `400`.
- [ ] The existing `GET /api/experiments/{id}/samples` route is unchanged and
      its tests stay green.
- [ ] A Julia route test exercises the corpus, filtered, and malformed-param
      responses.
- [ ] The running frontend is unaffected (the route is additive and
      unconsumed).

## References

- Master plan §3 (Phase 0): `docs/superpowers/plans/2026-05-17-himalaya-ui-redesign.md`
- Issue map I0.1 and the §3 contention table:
  `docs/superpowers/plans/2026-05-17-himalaya-ui-redesign-issue-map.md`
- Architecture decision 1 (experiment as facet): `docs/redesign-notes.md`
- Sibling route: `packages/HimalayaUI/src/routes_samples.jl:4`
- `q_units` inline logic to extract: `packages/HimalayaUI/src/routes_experiments.jl:8`
