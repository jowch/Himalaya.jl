# Index Scoring Redesign

**Date:** 2026-04-23  
**Status:** Approved — pending implementation

## Context

The current `score(index)` function in `src/index.jl` is an ad hoc formula:

```julia
(numpeaks(index) + totalprom(index)) * (1 - 0.25 * num_gaps) * rsquared
```

This has several known defects:

- **Unit mixing** — `numpeaks` (integer, typically 2–8) is added to `totalprom` (summed prominence, potentially 100× larger). The count term is effectively irrelevant.
- **Gap penalty can go negative** — `1 - 0.25 * num_gaps` becomes negative at 5 gaps.
- **No position sensitivity** — all gap positions are penalized equally, but physically peaks appear in order: if peak 3 is observed, peaks 1 and 2 should also be present. Missing the fundamental is much more suspicious than missing the 8th harmonic.
- **R² is fixed at 1 for single-peak indices** — artificially rewards minimal evidence.
- **No sharpness signal** — peaks from the same phase have consistent shape (wide-with-wide, sharp-with-sharp). Sharpness homogeneity is a meaningful quality signal that is already computed by `findpeaks` but discarded at the `indexpeaks` call site.

The score is load-bearing: `auto_group` in `pipeline.jl` uses it for greedy peak-claiming, so ranking errors propagate directly to auto-analysis results. False negatives (correct phase ranked low) are the primary concern (≈80%), with false positives secondary (≈20%).

## Design

### Score formula

```
score(index) = coverage(index) × consistency(index)
```

Both factors are in `[0, 1]`. The product is in `[0, 1]`.

**Coverage** — position-weighted fraction of expected peaks found:

```
coverage = Σ(1/rank_i  for each found peak i) / Σ(1/rank_j  for all expected peaks j)
```

`rank` is the ordinal position in the phase's ratio series (1 = fundamental). Harmonic weighting encodes the physical ordering constraint: missing the fundamental (weight 1.0) costs far more than missing the 8th peak (weight 0.125).

**Sharpness consistency** — coefficient of variation of peak sharpness across the index, inverted:

```
consistency = clamp(1 − std(sharpnesses) / mean(sharpnesses), 0, 1)
```

For a single-peak index, `consistency = 1.0` (undefined CV → assume perfect). Peaks from the same phase have similar shape; a mixed-sharpness index is penalized.

**R² as a hard gate (UI only):** R² is already stored in the `indices` table. The frontend filters out indices below a configurable threshold (default: 0.98). R² is not part of the score formula — it is almost always high for any assignment Himalaya produces, so it does not discriminate between candidates.

### Struct changes

`Index` currently stores `prom::Real` (summed prominence). Prominence is a peak-quality gate used during peak-finding, not a scoring signal for assignment quality. It is replaced by per-peak sharpness:

```julia
struct Index{P<:Phase}
    basis::Real
    peaks::SparseVector{<:Real, <:Integer}
    sharpness::SparseVector{<:Real, <:Integer}
end
```

`totalprom` and the `prom` field are removed. The `sharpness` sparse vector has the same sparsity structure as `peaks` — nonzero at ratio positions where a peak was assigned.

### `indexpeaks` call site

The pipeline currently discards sharpness:

```julia
# before
candidates = Himalaya.indexpeaks(peaks_result.q, peaks_result.prominence)

# after
candidates = Himalaya.indexpeaks(peaks_result.q, peaks_result.sharpness)
```

The internal `indexpeaks` dispatch chain threads the second positional argument (previously `proms`, now `sharpness`) into the `Index` constructor at the same positions. The argument name in the internal methods changes from `proms` to `sharpness`; the mechanics are identical.

### `score` implementation

Requires `Statistics.std` and `Statistics.mean` — `Statistics` is already a dependency (used by `R²`), but `std` must be confirmed imported in `src/index.jl`.

```julia
function score(index::Index{P}) where P
    n = length(phaseratios(P))
    found_idx, _ = findnz(index.peaks)
    denom = sum(1/r for r in 1:n)
    coverage = sum(1/r for r in found_idx) / denom

    _, sharps = findnz(index.sharpness)
    consistency = length(sharps) > 1 ? clamp(1 - std(sharps)/mean(sharps), 0, 1) : 1.0

    coverage * consistency
end
```

### `remove_subsets`

No changes needed — it already calls `score` to break ties; it will automatically use the new formula.

### `auto_group`

No changes needed — it already sorts by `score` descending.

### Frontend — R² gate

The index candidates panel (PhasePanel or equivalent) adds a client-side filter: indices with `r_squared < threshold` are hidden or visually dimmed. The threshold is configurable (default 0.95). The `r_squared` value is already returned by the `/api/exposures/:id/indices` endpoint and stored on each index row.

## Files to change

| File | Change |
|------|--------|
| `src/index.jl` | Remove `prom` field; add `sharpness` field to `Index`; remove `totalprom`; rewrite `score`; update `indexpeaks` internals |
| `packages/HimalayaUI/src/pipeline.jl` | Pass `peaks_result.sharpness` instead of `peaks_result.prominence` to `indexpeaks` |
| `test/index.jl` | Update tests for new struct shape and `score` formula |
| Frontend `PhasePanel` (or equivalent) | Add R² gate filter control |

## Out of scope

- Changing the peak-finding algorithm or sharpness definition
- Changing how `fit` or R² is computed
- Changing the `auto_group` greedy algorithm
- Per-user R² threshold persistence (use a sensible default)
