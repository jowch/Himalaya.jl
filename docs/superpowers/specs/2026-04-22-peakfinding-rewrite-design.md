# Peak-finding rewrite — design

**Date:** 2026-04-22
**Author:** Jonathan Chen (with Claude)
**Status:** awaiting-approval

## Problem

The current `findpeaks` in [src/peakfinding.jl](../../../src/peakfinding.jl) requires
manual per-scatter tuning of a prominence threshold `θ`. Even with tuning, it
returns false-positive peaks from broad form-factor features and fine-scale
wiggles. Two root causes:

1. The "noise" being filtered against is whatever Savitzky-Golay 2nd-derivative
   prominence happens to capture. Prominence has no physical meaning, so no
   single `θ` works across files.
2. The third column of `_tot.dat` (per-q intensity uncertainty σ(q),
   propagated by the azimuthal-integration software) is unused. We have a
   per-point noise scale and we throw it away.

## Goal

Replace `findpeaks` with an algorithm whose threshold is a **physical
SNR in sigmas** that does not need per-trace tuning, and which discriminates
narrow Bragg peaks from broad form-factor features by *construction*, not by
post-hoc filtering.

## Approach: CWT ridge detection + σ-based SNR validation

Two-stage pipeline. Stage 1 nominates scale-correct candidates using a
Continuous Wavelet Transform with Ricker (Mexican-hat) wavelets. Stage 2
fits each candidate locally and computes an honest SNR using the σ(q)
column.

### Pipeline

```
            log₁₀(I)              CWT(Ricker)         ridges      σ-validation
trace ───► transform ───► coefficients across ───► linking ───► fit Gaussian ───► peaks
   │                            ~12 scales         across       per candidate;
  q,σ                                              scales       reject if SNR<nσ
```

1. **Log transform.** Work on `y = log₁₀(I)`. SAXS spans ~3 orders of
   magnitude; without log compression high-q peaks are invisible to the wavelet.
2. **CWT** with Ricker wavelets at ~12 geometrically-spaced widths
   (default `[2, 2.6, 3.5, 4.5, 6, 8, 10, 13, 17, 22, 28, 35]` points).
3. **Ridge linking** across scales. A real peak appears as a chain of CWT
   maxima spanning multiple consecutive scales.
4. **σ-based validation.** For each surviving candidate, fit a Gaussian +
   linear baseline by weighted least squares (weights = 1/σ²), compute
   SNR = A / σ_A from the covariance matrix, reject if SNR < `nσ`.

## API

```julia
findpeaks(q, I, σ; nσ = 5, scales = nothing, min_ridge_length = 3) -> NamedTuple
```

| Argument | Meaning |
|---|---|
| `q`, `I`, `σ` | Three columns of `_tot.dat`. Equal length. |
| `nσ` | SNR threshold in sigmas. Same value across all scatters. Default 5. |
| `scales` | CWT widths in points. `nothing` ⇒ default 12-element geometric series above. |
| `min_ridge_length` | Candidate must persist across this many adjacent scales. Default 3. |

**Returns:** `(; indices, q, snr, width)` — four equal-length vectors. `indices`
are positions in the input arrays; `q` is the sub-pixel position from the
Gaussian fit; `snr = amplitude / σ_amplitude`; `width` is fitted FWHM in q-units.

**Breaks the old signature** `findpeaks(y)`. Acceptable because it doesn't
work well anyway and the package is at 0.4.x.

## Algorithm details

### CWT

Ricker wavelet: `ψ_a(t) = (2 / √(3a)·π^(1/4)) · (1 - t²/a²) · exp(-t²/(2a²))`.
Kernel half-window `5 · a` points; truncate beyond ±5σ. Direct convolution
(n=922 × 12 scales ≈ 1 ms). Reflect-pad at edges.

### Ridge linking

For each scale, find all positive local maxima of CWT coefficients (no height
threshold — Stage 2 does the filtering).

Greedy linking from largest scale to smallest:

- Start a ridge at each maximum at the largest scale.
- For the next-smaller scale, extend a ridge if there's a maximum within
  `±max(1, scale/4)` points of the ridge head.
- Allow up to 1 missing scale ("gap") before terminating a ridge.
- New maxima at smaller scales not extending any ridge start their own.

A ridge survives if it spans ≥ `min_ridge_length` scales (default 3).
Position is taken at the scale `a_opt` where the CWT coefficient is maximal.
Estimated peak width: `2.355 · a_opt · Δq`.

### σ propagation

For reference (not currently used in the algorithm, but worth documenting):
`σ_y = σ / (I · ln 10)` is the propagated 1σ on `log₁₀ I`. Future iterations
may weight the CWT by σ_y; this version uses σ only in Stage 2.

### Stage 2 fit

For each candidate at index `i_c` with width estimate `w_c = a_opt`:

- Window: `±2.5 · w_c` points around `i_c`, clipped to bounds.
- Model: `I(q) = A · exp(-(q - q₀)² / (2σ_w²)) + (m·q + b)`. 5 free parameters.
- Weighted LSQ via `LsqFit.jl`, weights `wᵢ = 1/σᵢ²`.
- Initial guesses: `A = I[i_c] - median(window I)`, `q₀ = q[i_c]`,
  `σ_w = w_c · Δq / 2.355`, baseline from window endpoints.
- `σ_A = stderror(fit)[1]` from the covariance matrix.
- Compute `SNR = A / σ_A`.
- Reject if `SNR < nσ` or fit fails to converge.

## Files touched

| File | Change |
|---|---|
| `src/peakfinding.jl` | Full rewrite. |
| `src/Himalaya.jl` | Drop `using Distributions` (was unused). |
| `Project.toml` | `+ LsqFit`, `- Distributions`. |
| `src/io.jl` | Move to `examples/io.jl` with header comment listing extra deps (`Images`, `TiffImages`, `DelimitedFiles`). |
| `README.md` | Update example block to new `findpeaks` signature. |
| `test/peakfinding.jl` | New file — synthetic ground-truth tests. |
| `test/peakfinding_real.jl` | New file — regression tests on labeled traces. |
| `test/data/` | New directory with three trace files copied from `scratch/`. |
| `test/runtests.jl` | Include the two new test files. |

`score` in [src/index.jl:264](../../../src/index.jl) is *not* touched in this
spec. Replacing its magic constants with an SNR-aware formula is a follow-up.

## Testing

### Synthetic ground-truth (`test/peakfinding.jl`)

Generate q ∈ [0.005, 0.4] on 922 points (matching example_tot.dat geometry).
Construct `I(q) = background + Σ peaks + noise`, then assert `findpeaks`
recovers the planted peaks.

| Scenario | Construction | Assertion |
|---|---|---|
| Single isolated peak, clean | 1 Gaussian, A/σ = 20, flat baseline | exactly 1 peak; q within Δq/2; SNR within 20% of true |
| Multiple peaks at known phase | 3 peaks at Pn3m ratios on power-law background, A/σ = 10 | exactly 3 peaks; q within Δq/2 |
| Below-threshold rejected | 1 peak A/σ = 2 | 0 peaks at default `nσ = 5` |
| Above-threshold detected | 1 peak A/σ = 6 | exactly 1 peak |
| Form-factor wiggle ignored | smooth Lorentzian background of width 100 pts | 0 peaks |
| Single-pixel spikes ignored | flat background + 3 single-pixel outliers | 0 peaks |
| Close-but-resolved | 2 Gaussians 8 pts apart, A/σ = 10 each | exactly 2 peaks |

### Regression on real traces (`test/peakfinding_real.jl`)

Three traces hand-labeled via `scratch/label_peaks.html`:

```julia
const LABELED_PEAKS = Dict(
    "example_tot.dat" => [
        0.053439, 0.065558, 0.075081, 0.091962, 0.107111,
        0.113171, 0.124858, 0.130052, 0.139142, 0.155590,
    ],
    "cubic_tot.dat" => [
        0.046115, 0.058678, 0.065177, 0.071675, 0.079906,
        0.082939, 0.092037, 0.101568, 0.112832, 0.117164,
        0.121496, 0.124529, 0.130594, 0.137959, 0.143591,
        0.145324, 0.152689, 0.159620, 0.171318, 0.177816,
        0.189513, 0.194712, 0.200777,
    ],
    "form-factor_tot.dat" => Float64[],
)
```

For each trace, run `findpeaks(q, I, σ)` with default `nσ = 5` and assert:

- Number of returned peaks matches `length(LABELED_PEAKS[name])`.
- Each labeled q-value has a returned peak within tolerance
  `2 · median(diff(q))` (≈ 2 q-grid samples — labels were placed via
  free-form clicking on data points, not snapped to local maxima, so the
  algorithm's sub-pixel fit may legitimately land 1-2 samples away).
- For `form-factor_tot.dat`: zero returned peaks (the diagnostic test —
  proves we've solved the original problem).

The cubic trace is dense (23 peaks); if matching exactly proves brittle, we
may relax to "≥ 90% of labels recovered, ≤ 10% false positives" — decided
during implementation based on first-run results.

### Integration test (`test/index.jl`, augment existing)

End-to-end: load `example_tot.dat`, run `findpeaks`, run `indexpeaks`,
assert the top-scoring `Index` matches a known-correct phase + basis. Confirms
the API change didn't break the downstream pipeline.

## Out of scope

- Performance / timing tests.
- Reworking `score` in `src/index.jl`.
- Reworking `indexpeaks` to consume SNR more meaningfully.
- Loading `_tot.dat` files from disk (caller's responsibility).
- Multi-trace batch processing.

## Open issues

None at design time. If the cubic regression test proves brittle in
implementation, decide between (a) loosening tolerance, (b) re-labeling with
finer-grain criteria, (c) adjusting default scale range.
