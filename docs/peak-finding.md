# Peak-finding strategy

Himalaya finds Bragg peaks in 1D SAXS integration traces using two
independently-scored criteria combined by an AND-gate. Neither criterion
requires per-trace tuning — both thresholds are set adaptively from the
data's own distribution.

## The two-criterion design

A genuine Bragg peak must satisfy two conditions simultaneously:

**Criterion A — topological prominence.** The peak must stand clearly
above its local floor. This is formalized as *topological persistence*:
imagine flooding the signal from above; a peak is born when the water
level reaches its apex and dies when its hill merges with a taller
neighbour. The prominence is the vertical distance between birth and
death — a physically meaningful, scale-free measure of how much a peak
dominates its neighbourhood. [1, 2] Implemented via
[Peaks.jl](https://github.com/halleysfifthinc/Peaks.jl) `peakproms`.

**Criterion B — sharpness.** The peak must be sharply curved, not a
broad hump. Sharpness is the magnitude of the (smoothed) second
derivative at the peak position. This criterion discriminates narrow
Bragg peaks from slowly-varying form-factor features: a form-factor
envelope may have high prominence, but its curvature is negligible at
the scale of a Bragg peak. Implemented via Savitzky-Golay 2nd derivative
(default) or multi-scale Ricker CWT.

The AND-gate means neither criterion can be satisfied by accident. A
single-pixel noise spike has high sharpness but low prominence (it barely
rises above its neighbours). A broad form-factor feature has prominence
but no sharpness. Only a narrow, tall Bragg peak passes both.

## Adaptive thresholding

For each criterion, the threshold is set via the *kneedle algorithm*
applied to the sorted distribution of scores across all candidates in
that trace. The knee of the curve separates the "steep" high-scoring
candidates (signal) from the "flat" low-scoring tail (noise), without
any fixed cutoff. This means the same algorithm works across different
instruments, sample geometries, and count rates.

Manual overrides (`prom_floor`, `sharp_floor`) are available for
diagnostic or scripted workflows, but should not be needed in normal use.

## Gating no-peak traces

Kneedle assumes the candidate-prominence distribution has a real elbow.
On a trace with no Bragg peaks at all (e.g., a pure form-factor envelope
plus shot noise), kneedle picks a fictitious knee inside what is
effectively a noise decay curve, and a handful of high-q candidates can
squeak through.

Two adaptive backstops sit on top of the kneedle AND-gate:

**High-q post-trim (`q_trim_high`, default `0.05`).** Peaks at the very
top of the q-range are discarded after kneedle. The post-filter
ordering matters — truncating before kneedle would remove the noise
candidates that anchor the lower part of the sorted-prominence curve,
shifting the knee upward and rejecting real peaks. We have observed
this regression on real fixtures.

**Relative-prominence floor (`prom_ratio_floor`, default `30.0`).**
After kneedle, every kept peak must satisfy
`prominence ≥ prom_ratio_floor × median(candidate_prominence)`. The
median is computed over *all* candidates produced by `persistence`, so
on a noise-only trace it tracks the noise floor; on a trace with real
peaks it is dominated by noise candidates and the gate is a no-op
(real peaks sit far above the kneedle threshold, which is itself far
above the ratio floor). The floor is skipped when the candidate count
is below `RATIO_FLOOR_MIN_CANDIDATES` (currently 20) — at low candidate
counts the median is dominated by the peaks themselves and would
suppress real signal.

## Why not σ-based SNR?

The `_tot.dat` files include a per-point intensity uncertainty σ(q)
propagated by the azimuthal integration software. Computing a physical
SNR (`A / σ_A`) from a per-candidate Gaussian fit is an intuitive
approach, but parametric shape fitting is brittle when peaks sit on a
steeply falling background or overlap slightly.

More importantly, phase-ratio coherence (via `indexpeaks`) is a much
stronger filter than per-peak SNR: a single 5σ excess could be noise,
but four peaks at Pn3m's √2 : √3 : √4 : √6 ratios cannot be. The
prominence + sharpness AND-gate gets candidates to `indexpeaks`;
`indexpeaks` validates them by phase geometry. Shape fitting at the
candidate-selection stage would add complexity without improving the
overall recall of correct phase assignments.

## Key defaults and when to change them

**`normalize_by_σ = false` (default).** For Poisson-dominated SAXS,
σ ≈ √I, so dividing I by σ produces √I — a square-root compression that
shrinks the prominence of real peaks relative to noise. The default is
`false` for SAXS data. Set to `true` for data with approximately
homoscedastic noise (e.g., detector read-noise dominated).

**`sharpness_method = :savgol` (default).** The Savitzky-Golay 2nd
derivative is fast and sufficient when peaks are similar in width across
a trace. Switch to `:cwt` (multi-scale Ricker) when peaks span a wide
range of widths within one trace, or when single-pixel spike rejection
is important.

## References

[1] Huber, S. (2021). [Persistent Topology for Peak Detection](https://www.sthu.org/blog/13-perstopology-peakdetection/index.html). *sthu.org*. — Accessible explanation of the flooding intuition applied specifically to 1D peak detection.

[2] Edelsbrunner, H., Letscher, D., & Zomorodian, A. (2002). Topological Persistence and Simplification. *Discrete & Computational Geometry*, 28, 511–533. https://doi.org/10.1007/s00454-002-2885-2 — Foundational paper introducing persistent homology.

## Intentionally out of scope

- **Parametric shape fitting.** Brittle on real SAXS backgrounds; phase-ratio coherence is a stronger downstream filter.
- **Sub-pixel peak positions.** Peaks are returned at grid positions. Parabolic interpolation is a straightforward follow-up if needed.
- **Background subtraction** (SNIP, asymmetric least squares). A potential follow-up for steep-background cases.
- **Multi-trace batch statistics.** Thresholds are per-trace; no pooling across scans.
- **q-banded kneedle** (independent thresholds per q-band). Considered when designing the no-peak-trace gating; rejected because every band of a noise-only trace is itself noise, and per-band kneedle would still find spurious knees.
- **σ(q)-based noise floor** using the third column of `_tot.dat`. Considered for the same purpose; rejected because parametric noise modelling at the candidate-selection stage re-opens the design dimension this document already closed (see "Why not σ-based SNR" above). The relative-prominence floor uses the data's own prominence distribution and avoids this.
- **Low-q trim.** The empirical failure mode is high-q only. A symmetric trim risks clipping Lamellar 1st-order peaks; we wait for a fixture exhibiting low-q noise contamination before adding one.
