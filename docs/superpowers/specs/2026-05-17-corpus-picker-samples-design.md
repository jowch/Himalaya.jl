# Corpus-wide `GET /api/picker-samples` route (I0.3)

**Issue:** #158
**Date:** 2026-05-17
**Status:** Design approved

## Goal

Add a corpus-wide `GET /api/picker-samples` endpoint returning the sample-picker
projection across all experiments, as a sibling to the experiment-scoped
`GET /api/experiments/{eid}/picker-samples`.

The series scoping step (#174) lets a user pick samples for a series from the
whole corpus. The picker projection — a lightweight per-sample shape for
selection UIs — exists only per experiment today.

## Scope

- New route `GET /api/picker-samples` in `routes_picker.jl`.
- New corpus method of `picker_samples` in `comparisons.jl`, via a dispatch
  refactor that shares the existing experiment-agnostic body.
- Julia tests in `test_picker_samples_route.jl`.

Append-only: no schema change, no events, no migration.

## Out of scope

- The experiment-scoped `GET /api/experiments/{eid}/picker-samples` route — its
  observable behavior is unchanged.
- The scoping picker UI (#174).
- Any frontend change. `picker-samples` has no frontend consumer today
  (api.ts / queries.ts do not reference it); the running frontend is unaffected.
- Chunked-`IN` / temp-table pagination of the corpus query. The exposures and
  `sample_tags` bulk queries bind one placeholder per sample, so a corpus
  exceeding `SQLITE_MAX_VARIABLE_NUMBER` (999 on pre-3.32 SQLite, 32766 after)
  would throw. This is not a concern at realistic SAXS corpus sizes; chunking is
  deferred until a deployment approaches that bound.

## Design

### Helper refactor — `comparisons.jl`

`picker_samples(db, experiment_id)` builds the projection in three parts:

1. `SELECT id, experiment_id, name, display_name, notes FROM samples
   WHERE experiment_id = ? ORDER BY id` — the **only** experiment-scoped query.
2. Bulk fetch of exposures and `sample_tags` keyed by `sample_id IN (...)`.
3. Grouping, indexing-exposure resolution, and per-sample dict assembly.

Parts 2–3 are driven purely by `sample_ids` and are experiment-agnostic. Extract
them into a private helper:

```julia
_picker_samples_projection(db, samples) -> Vector{Dict{Symbol, Any}}
```

`samples` is the already-fetched row table from part 1. Then expose two thin
public methods that dispatch on arity:

- `picker_samples(db, experiment_id)` — runs the experiment-scoped `samples`
  query (`WHERE experiment_id = ? ORDER BY id`), delegates to the helper.
  **Behavior unchanged** — all existing `picker_samples helper` testsets stay
  green.
- `picker_samples(db)` — runs the corpus `samples` query (no `WHERE`,
  `ORDER BY experiment_id, id`), delegates to the same helper. Returns
  `Dict{Symbol, Any}[]` on an empty corpus.

The `samples` SELECT column list is identical for both methods, so
`experiment_id` is already present in every sample dict — corpus consumers group
by it client-side. `ORDER BY experiment_id, id` yields stable, experiment-grouped
output.

The exposure-resolution logic is correct corpus-wide without change: exposures
are grouped by `sample_id`, and each exposure belongs to exactly one sample, so
there is no cross-experiment leakage. The existing "multi-experiment isolation"
testset confirms this property.

### Route — `routes_picker.jl`

Add `@get "/api/picker-samples"` inside `register_picker_routes!()`, a sibling to
the existing experiment-scoped block. No path params, no query params. Read-only
— no `with_idempotency`, no SSE, no event-log row, consistent with the other
picker routes.

```julia
@get "/api/picker-samples" function(req::HTTP.Request)
    db = current_db()
    HTTP.Response(200, ["Content-Type" => "application/json"],
                  JSON3.write(picker_samples(db)))
end
```

No registration change is needed — `register_picker_routes!()` is already called
from `server.jl`.

## Response shape

Identical per-sample shape to the experiment-scoped route. For each sample:

- `:sample` — sample row dict (`id, experiment_id, name, display_name, notes`)
  with a `:tags` vector.
- `:indexing_exposure_id` — `Int` (highest-id `selected` exposure, else highest
  id) or `null` for zero-exposure samples; key always present.
- `:all_exposures` — vector of `{:id, :filename, :selected::Bool}`, `id` ASC.

The corpus payload is the concatenation of every experiment's samples, ordered
by `(experiment_id, id)`.

## Testing

Append to `test_picker_samples_route.jl`:

- **`picker_samples(db)` helper testset:** a multi-experiment corpus returns all
  samples ordered by `(experiment_id, id)`; an empty corpus returns `[]`;
  cross-experiment exposure resolution is correct (an exposure with a larger
  global id in another experiment does not leak into a sample's
  `indexing_exposure_id`).
- **`GET /api/picker-samples` route testset:** mirrors the existing
  experiment-route test — status 200, payload shape, `selected` rendered as a
  JSON bool, `indexing_exposure_id` present-but-`null` for a zero-exposure
  sample.

The existing experiment-scoped `picker_samples helper` and route testsets must
remain green unchanged — they are the regression guard for the refactor.

## Dependencies

- No upstream dependencies — Wave A.
- Blocks the series scoping step (#174).
- Shares `routes_picker.jl` with #157 (`GET /api/sample-tags`); append-only
  contention. Land #157 first to avoid the append conflict.

## References

- Master plan §3 (Phase 0): `docs/superpowers/plans/2026-05-17-himalaya-ui-redesign.md`
- Issue map I0.3: `docs/superpowers/plans/2026-05-17-himalaya-ui-redesign-issue-map.md`
- Sibling route: `packages/HimalayaUI/src/routes_picker.jl:71`
- Helper: `packages/HimalayaUI/src/comparisons.jl:506`
