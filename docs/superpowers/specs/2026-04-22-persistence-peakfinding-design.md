# Peak-finding rewrite v2 — persistence + sharpness — design

**Date:** 2026-04-22
**Author:** Jonathan Chen (with Claude)
**Status:** awaiting-approval
**Supersedes:** [2026-04-22-peakfinding-rewrite-design.md](2026-04-22-peakfinding-rewrite-design.md)

## Problem

The first rewrite (CWT + Lorentzian fit + σ-based SNR) reached ~60% recall
on hand-labeled real SAXS traces with zero false positives. Iteration on
fit initialisation, threshold defaults, and edge-case handling did not
materially improve recall. Two independent literature passes (one local,
one via Perplexity) converged on the same diagnosis: **the current
pipeline conflates two independent peak-quality criteria — prominence and
sharpness — into a single gating stage, and embeds a parametric peak
shape (Lorentzian) into Stage 2.**

The user has further clarified that:
- Peak quality is judged primarily by (A) topological prominence and (B)
  curvature/sharpness; downstream phase-ratio coherence is the validation,
  not σ-based SNR.
- Per-trace tuning is unacceptable; thresholds must be self-calibrating.
- Hand labels should drive *testing*, not algorithm parameters.
- Peak-shape assumptions (Lorentzian, Gaussian) are physics baggage to
  shed — this is fundamentally a signal-processing problem.

## Goal

Replace `findpeaks` with a shape-agnostic detector that evaluates
prominence and sharpness independently, applies adaptive (kneedle)
thresholds derived from the data's own distribution, and combines them
via an AND-gate. Drop all parametric peak-shape fitting.

## Approach: persistence + sharpness, both kneedle-thresholded

### Pipeline

```
                                  ┌──────────────────────┐
   (q, I, σ)  ──► preprocess ────►│  persistence(y)      │──► (idx, prominence)
                  y = I or I/σ    │  (1D persistent       │
                                  │   homology = prom)   │
                                  └──────────────────────┘
                                              │
                                              ▼
                                  ┌──────────────────────┐
                                  │  sharpness(y)        │──► sharpness vector
                                  │  :savgol or :cwt     │     (one per sample)
                                  └──────────────────────┘
                                              │
                                              ▼
                                  ┌──────────────────────┐
                                  │  knee(values)        │──► auto threshold
                                  │  (kneedle algorithm) │     per axis
                                  └──────────────────────┘
                                              │
                                              ▼
                                  AND-gate: passes both → real peak
```

1. **Preprocess.** Optionally divide intensity by σ (`y = I ./ σ`) so
   prominence is in noise-relative units. Default `true` when σ provided.
2. **Persistence.** Compute prominence for every local max in `y`.
   Mathematically identical to 1D persistent homology. Implementation:
   thin wrapper around `Peaks.peakproms`.
3. **Sharpness.** Compute a sharpness measure for every sample of `y`;
   look up at candidate positions. Two interchangeable methods (kwarg):
   - `:savgol` — magnitude of smoothed 2nd derivative (fast, single-scale).
   - `:cwt` — max Ricker response across a small scale bank (multi-scale).
4. **Adaptive thresholds.** Run the kneedle algorithm on each sorted
   value list to find the elbow separating signal from noise.
5. **AND-gate.** A candidate is accepted iff its prominence ≥ knee(proms)
   AND its sharpness ≥ knee(sharps).

### What is NOT in this design

- Parametric peak shape fitting (Lorentzian, Gaussian, Voigt).
- σ-based SNR scoring (`A / σ_A` from a fit covariance).
- Sub-pixel peak position refinement (deferred — return grid positions).
- Label-trained penalties (PeakSeg-style).
- Bayesian / RJ-MCMC inference.
- Background subtraction (SNIP, asymmetric LS).

These are all live algorithm-design avenues but explicitly excluded from
v2. v2 is the smallest, most readable thing that operationalises the
user's stated A+B criteria.

## API

```julia
findpeaks(q, I, σ; normalize_by_σ = true,
                    sharpness_method = :savgol,  # or :cwt
                    prom_floor   = nothing,    # nothing → kneedle adaptive
                    sharp_floor  = nothing) -> NamedTuple

findpeaks(q, I;     normalize_by_σ = false, kwargs...)  # σ optional
```

| Kwarg | Meaning | Default |
|---|---|---|
| `normalize_by_σ` | If `true`, work on `I ./ σ` so prominence is in noise units | `true` (when σ given) |
| `sharpness_method` | `:savgol` (single-scale 2nd derivative) or `:cwt` (multi-scale Ricker max) | `:savgol` |
| `prom_floor` | Manual prominence threshold; `nothing` ⇒ kneedle | `nothing` |
| `sharp_floor` | Manual sharpness threshold; `nothing` ⇒ kneedle | `nothing` |

**Returns:**

```julia
(; indices, q, prominence, sharpness)
```

Four equal-length vectors. `indices` are positions in the input arrays;
`q` is `q[indices]`; `prominence` and `sharpness` are the per-peak
quality measures.

**Downstream compatibility.** `indexpeaks(peak_qs, peak_proms)` continues
to work — call as `indexpeaks(pk.q, pk.prominence)`. The shape of the
return type is preserved for the consumer.

## Module structure

| File | Responsibility | LOC |
|---|---|---|
| `src/peakfinding.jl` | Top-level `findpeaks`; composes the rest | ~30 |
| `src/persistence.jl` | `persistence(y)` — wrapper over `Peaks.peakproms` | ~5 |
| `src/sharpness.jl` | `sharpness(y; method, ...)`; `sharpness_savgol`, `sharpness_cwt`, `ricker`, `savitzky_golay`, `convolve_reflect` | ~80 |
| `src/threshold.jl` | `knee(v)` — kneedle on sorted descending vector | ~30 |
| `src/Himalaya.jl` | Module entry, exports | (small edit) |

The previous v1 implementation (~430 LOC across `src/peakfinding.jl`)
is fully replaced.

## Algorithm details

### Persistence (`src/persistence.jl`)

```julia
using Peaks: findmaxima, peakproms

"""
    persistence(y) -> (; indices, prominence)

Compute the topological persistence — equivalently, the prominence — of
every local maximum in the 1D signal `y`.

Intuition: imagine flooding the signal from above. Each local max is
born when the water level passes it and dies when its hill merges with
a taller neighbour (or hits the array boundary). The vertical distance
between birth and death is its prominence.
"""
function persistence(y)
    maxima = findmaxima(y)
    with_proms = peakproms(maxima)
    (indices = with_proms.indices, prominence = with_proms.proms)
end
```

### Sharpness (`src/sharpness.jl`)

The user's existing `savitzky_golay(m, n, y; order = 0)` (general-purpose
SG implementation from the v1 source — preserved from git history) is
restored as the underlying engine. The new `sharpness_savgol` is a thin
wrapper.

```julia
const DEFAULT_SCALES = [2.0, 3.0, 5.0, 8.0, 13.0]   # 5 scales, geometric

"""
    sharpness(y; method = :savgol, scales = DEFAULT_SCALES, m = 5)

Compute a per-sample sharpness measure for the 1D signal `y`. High value
= sharply curved local feature.
"""
function sharpness(y; method = :savgol, scales = DEFAULT_SCALES, m = 5)
    method === :savgol     && return sharpness_savgol(y, m)
    method === :cwt        && return sharpness_cwt(y, scales)
    throw(ArgumentError("Unknown sharpness method: $method (use :savgol or :cwt)"))
end

"""
    sharpness_savgol(y, m) -> Vector{Float64}

-d²y/dx² via Savitzky-Golay smoothing (window 2m+1, polynomial order 4).
A sharp local maximum has strongly negative d²y; sign flipped so high =
sharp.
"""
sharpness_savgol(y, m) = -savitzky_golay(m, 4, y; order = 2)

"""
    sharpness_cwt(y, scales) -> Vector{Float64}

Max Ricker-wavelet response across `scales`, normalised by 1/a so
coefficients are scale-comparable.
"""
function sharpness_cwt(y, scales)
    n = length(y)
    out = fill(-Inf, n)
    for a in scales
        m = ceil(Int, 5a)
        kernel = [ricker(t, a) for t in -m:m]
        response = convolve_reflect(y, kernel) ./ a
        out .= max.(out, response)
    end
    out
end

ricker(t, a) = (2 / (sqrt(3a) * π^0.25)) * (1 - (t/a)^2) * exp(-(t^2) / (2a^2))

"""
    convolve_reflect(y, kernel) -> Vector{Float64}

Direct 1D convolution of `y` with `kernel` (length 2m+1), reflecting at
edges (no wraparound, no zero-padding artefacts).
"""
function convolve_reflect(y, kernel)
    n = length(y)
    m = (length(kernel) - 1) ÷ 2
    out = zeros(n)
    for i in 1:n, k in 1:length(kernel)
        idx = i + (k - m - 1)
        idx = idx < 1 ? 2 - idx : idx > n ? 2n - idx : idx
        out[i] += kernel[k] * y[idx]
    end
    out
end

# savitzky_golay(m, n, y; order) — preserved verbatim from the v0.4.5
# `src/peakfinding.jl` (commit fddd611). General-purpose SG: any window,
# any polynomial order, any derivative order (single value or vector).
# Requires `using LinearAlgebra: I` (the identity matrix constructor).
function savitzky_golay(m, n, y; order = 0)
    num_y = length(y)
    z = -m:m
    J = zeros(2m + 1, n + 1)

    for i = 0:n
        @inbounds J[:, i+1] .= z .^ i
    end

    # The convolution term matrix
    C = J' \ I(n .+ 1)[:, order .+ 1] # = inv(J' * J) * J' = pinv(J)
    Y = zeros(num_y, length(order))

    for i in 1:num_y
        if i <= m
            window_indices = abs.(z .+ i) .+ 1
        elseif i > num_y - m
            window_indices = -abs.(z .+ i .- num_y) .+ num_y
        else
            window_indices = z .+ i
        end

        for j in eachindex(order)
            @inbounds Y[i, j] = C[:, j]' * y[window_indices]
        end
    end

    length(order) == 1 ? Y[:, 1] : Y
end
```

### Threshold (`src/threshold.jl`)

```julia
"""
    knee(v) -> Float64

Find the elbow of a sorted descending sequence `v` using the kneedle
algorithm: the value at the point of maximum perpendicular distance from
the chord connecting the first and last entries (after normalising both
axes to [0, 1] so the distance is scale-invariant).

Used as an adaptive threshold: values ≥ knee(v) are signal, values
< knee(v) are noise. No magic numbers, no labels, no per-trace tuning —
the threshold emerges from the shape of the data's own distribution.

Edge cases:
- Empty input → 0.0 (nothing passes).
- Single value → that value.
- Constant sequence → v[1] (no elbow exists; nothing passes).
- No clear elbow (convex-up curve) → v[1] (conservative; nothing passes).

Precondition: `v` must be sorted descending.
"""
function knee(v)
    n = length(v)
    n == 0 && return 0.0
    n == 1 && return float(v[1])
    n == 2 && return float(v[1])

    y_min   = v[end]
    y_range = v[1] - y_min
    y_range == 0 && return float(v[1])

    best_score = 0.0
    best_idx   = 1
    for i in 2:(n-1)
        x = (i - 1) / (n - 1)
        y = (v[i] - y_min) / y_range
        score = 1 - x - y                 # below the chord
        if score > best_score
            best_score = score
            best_idx   = i
        end
    end

    float(v[best_idx])
end
```

### Top-level findpeaks (`src/peakfinding.jl`)

```julia
include("persistence.jl")
include("sharpness.jl")
include("threshold.jl")

"""
    findpeaks(q, I, σ; normalize_by_σ = true, sharpness = :savgol,
                       prom_floor = nothing, sharp_floor = nothing)
    findpeaks(q, I;    normalize_by_σ = false, kwargs...)

Detect peaks in a 1D signal using topological prominence and curvature,
both adaptively thresholded. Shape-agnostic; no per-trace tuning.

Returns a NamedTuple `(; indices, q, prominence, sharpness)` of equal-
length vectors, sorted by ascending q.
"""
function findpeaks(q, I, σ; normalize_by_σ    = true,
                              sharpness_method = :savgol,
                              prom_floor       = nothing,
                              sharp_floor      = nothing)
    @assert length(q) == length(I) == length(σ) "q, I, σ must be equal length"

    y = normalize_by_σ ? I ./ σ : I

    cands           = persistence(y)
    sharps_full     = sharpness(y; method = sharpness_method)
    sharps_at_peaks = sharps_full[cands.indices]

    pf = something(prom_floor,  knee(sort(cands.prominence; rev = true)))
    sf = something(sharp_floor, knee(sort(sharps_at_peaks;  rev = true)))

    keep = (cands.prominence .>= pf) .& (sharps_at_peaks .>= sf)
    idx_keep = cands.indices[keep]
    perm = sortperm(idx_keep)

    (indices    = idx_keep[perm],
     q          = q[idx_keep[perm]],
     prominence = cands.prominence[keep][perm],
     sharpness  = sharps_at_peaks[keep][perm])
end

findpeaks(q, I; normalize_by_σ = false, kwargs...) =
    findpeaks(q, I, ones(length(I)); normalize_by_σ = normalize_by_σ, kwargs...)
```

## Files touched

| File | Change |
|---|---|
| `src/peakfinding.jl` | Full replacement (drops ricker / cwt / find_ridges / fit_peak / old findpeaks) |
| `src/persistence.jl` | Create (5 LOC + docstring) |
| `src/sharpness.jl` | Create (~80 LOC including restored `savitzky_golay`) |
| `src/threshold.jl` | Create (~30 LOC) |
| `src/Himalaya.jl` | Verify exports (`findpeaks` already exported) |
| `Project.toml` | Drop `LsqFit` (no fitting); keep `Peaks` |
| `README.md` | Update example to new return shape (`pk.prominence` instead of `pk.snr`) |
| `test/runtests.jl` | Include the four new test files |
| `test/persistence.jl` | Create — unit tests |
| `test/sharpness.jl` | Create — unit tests for both methods |
| `test/threshold.jl` | Create — unit tests for kneedle |
| `test/peakfinding.jl` | Create — synthetic stress tests (replaces v1 tests) |
| `test/peakfinding_real.jl` | Re-baseline `RECALL_FLOOR` and `SPURIOUS_CEILING` after running |

## Testing

### Unit tests (per module)

**persistence.jl**
- Single Gaussian peak on flat baseline → 1 entry, prominence = peak − baseline.
- Two peaks (one taller) → tall: prominence to global min; short: prominence to saddle.
- Monotonic signal → empty.
- Plateau peak → handled (one entry at midpoint).

**threshold.jl (knee)**
- L-shaped `[100, 90, 80, 5, 4, 3, 2]` → returns ~80.
- Linear descent `[10, 9, 8, 7, 6, 5]` → returns middle.
- Constant `[5, 5, 5, 5]` → returns 5.
- Empty → 0.0.
- Single value → that value.

**sharpness.jl**
- Both methods on a Gaussian peak: peak position has highest sharpness; flat regions ~0.
- `:savgol`: m=5 + order-4 SG coefficients match known reference values.
- `:cwt`: response at peak position highest at scale matching the Gaussian σ.
- Invalid method → `ArgumentError`.

### Synthetic integration tests (peakfinding.jl)

| Scenario | Construction | Assertion |
|---|---|---|
| Single isolated Lorentzian | A=100, γ=0.003, flat baseline, σ=5 | exactly 1 peak, q within Δq |
| 3 peaks at Pn3m ratios | power-law background + 3 Lorentzians | exactly 3 peaks at expected positions |
| Smooth form-factor wiggle | wide Lorentzian background, no narrow peaks | 0 peaks |
| Single-pixel spikes | flat + 3 single-pixel outliers | `:cwt` rejects all; `:savgol` may detect (documented limitation) |
| Two close-but-resolved peaks | two Lorentzians ≥ 3 HWHM apart | exactly 2 peaks |
| σ-normalize on/off | same trace both ways | results consistent in shape |

### Regression tests (peakfinding_real.jl)

Same three traces (`example_tot.dat`, `cubic_tot.dat`, `form-factor_tot.dat`),
same `LABELED_PEAKS`. Same `RECALL_FLOOR` / `SPURIOUS_CEILING` pattern,
**re-baselined** after the new algorithm runs. Floors are documentation
of current behaviour; tightening them is a deliberate choice in the same
commit that improves the algorithm.

### Integration test (test/index.jl)

Existing end-to-end test (load `example_tot.dat`, run `findpeaks`, run
`indexpeaks`, assert top-scoring index has a sensible phase) updated for
the new return shape.

## Out of scope

- Performance benchmarks (no hot-loop concern at n≈1000).
- Sub-pixel peak position refinement.
- Label-trained penalties (PeakSeg).
- Bayesian RJ-MCMC.
- Background subtraction.
- Reworking `score` in `src/index.jl` (separate follow-up).

## Open issues

- Single-pixel spike behaviour with `:savgol` is a known limitation. If
  it manifests in real traces, the path forward is either pre-filtering
  with a hampel/median filter or making `:cwt` the default.
- Re-baselined recall floors may be lower than v1's. If recall regresses
  significantly, the design holds but iteration on (a) sharpness method,
  (b) `:cwt` scale bank, or (c) `normalize_by_σ` default may be needed.
