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

**Relative-prominence floor (`prom_ratio_floor`, default `15.0`).**
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

## The indexing tolerance is asymmetric

Peak *finding* hands candidates to `indexpeaks`, which decides whether their
q-ratios match a phase. Before any per-ratio residual is computed, `indexpeaks`
prunes the candidate ratios with a single mask (`src/index.jl:142-143`):

```julia
# TODO: work out addition tol on max side
mask = 1 - tol .<= observed_ratios .<= maximum(ratios)
```

The two ends of this window are not treated the same. The lower bound is
**softened** by `tol` (`1 - tol`), so a peak landing just below the fundamental
ratio of 1 still survives the mask. The upper bound is a **hard cutoff** at the
phase's highest ratio: a peak whose ratio exceeds `maximum(ratios)` by even a
hair is rejected outright, with no tolerance band. The in-code `TODO` is an
explicit acknowledgement that the high side still needs its own tolerance.

Two consequences follow, and both are easy to trip over:

- **`tol` is not symmetric at the high end.** A small positive measurement error
  on the top reflection — the very peaks most vulnerable to it, since they sit
  in the noisiest, steeply-falling tail of the trace — can push a real peak past
  the cutoff and out of consideration, even though the same-sized error on a
  low-order peak would be tolerated.
- **Extending a phase's ratio table changes what the *old* top ratio can
  match.** Because the cutoff is `maximum(ratios)`, adding a higher-order
  reflection raises the ceiling. Peaks that previously sat just above the former
  top ratio (and were masked out) now fall inside the window and become
  available — so adding reflections is not a purely additive change to a phase
  definition; it can shift which peaks match the existing top ratio.

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

## Hybrid detector (opt-in)

`findpeaks(q, I, σ; method = :hybrid)`. The default `:persistence` path is
unchanged and remains the default.

**Why.** The persistence detector thresholds prominence in raw intensity.
On a steep background that threshold is set by the large low-q features, so
weak peaks and shoulders on the steep flank fall below it even when they
are clearly significant against the local noise. A detector that works in
noise-relative units after background removal does not have this problem,
but alone it adds more spurious candidates than the persistence detector.

**What it does** (code in `src/matchedfilter.jl`):

1. Run the persistence detector with `q_trim_high = 0.07`.
2. Compute a matched-filter response per kernel width.
   - **Masked data.** A maximal run of points that are exactly `0.0` or
     nonfinite is masked when the median of the nearest 7 finite, non-run
     points on every available side exceeds 3 times their rms `σ` (a beamstop,
     a detector gap, at either end or in the interior). Exactly zero runs in a
     near-zero region (a true zero baseline, a negative tail) fail the test and
     stay data. Masked and nonfinite points are *missing*.
   - **Span.** The trace is cut to the span between the first and last present
     points, so a masked beamstop is treated like the start of the trace.
     Inside the span, missing points are filled by linear interpolation for the
     background only.
   - **Background.** SNIP baseline (half-window 60 points) using the LLS
     operator on data scaled by `s`, the median positive `σ`:
     `v = log(log(sqrt(max(I/s, 0) + 1) + 1) + 1)`, clip, invert, and multiply
     by `s`. The operator compresses the dynamic range without a hole at zero,
     and `s` scales with `I`, so the result does not depend on intensity units.
   - **Filter.** The residual is correlated with unit-norm, zero-mean Ricker
     kernels of widths 2, 3, 4.5, 7, 11, 16 points and divided by
     `sqrt(Σ k² σ²)`. A width is defined at a point when its kernel fits inside
     the span and covers no missing point, so filled data never enters a
     response.
   - **Two scales.** The same response is put on two scales.
     - `z` (used to *add* peaks): each width is normalised by a centre and
       scale estimated from its own values with signal removed. Median and
       1.4826·MAD are computed, points more than 3.5 scales from the centre
       (and a few neighbours) are dropped, and this repeats until the set stops
       changing. A width with fewer than 30 points (or a fifth of its points)
       left is skipped. Removing the signal keeps strong or dense peaks from
       inflating the scale, while the scale still absorbs noise that `σ` does
       not describe (correlated noise, background misfit at wide widths).
     - `zraw` (used to *keep* persistence peaks): each width is divided by one
       noise factor for the whole trace, the clipped scale of the normalised
       second differences `(I[i-1] - 2I[i] + I[i+1]) / sqrt(σ[i-1]² + 4σ[i]² + σ[i+1]²)`.
       Second differences hardly respond to a smooth background or to a peak
       more than a few points wide. So `zraw` stays correct on a trace that is
       mostly signal, where no estimate from the responses themselves can be.
   - **Short traces.** Traces with fewer than 150 points (the widest kernel
     needs 129 points to fit) return the persistence result (with
     `q_trim_high = 0.07`) unchanged, with `z` all `NaN`.
3. Keep a persistence peak unless all three hold:
   - `zraw` is defined at it;
   - the max of `zraw` over the widths and the points within ±3 is below `kb`
     (default 5);
   - no width there is blocked by missing data: its window covers a missing
     point or would cross a span edge next to a beamstop.

   This removes persistence peaks that are not significant against the noise.
   A peak too close to the trace ends for any kernel to be centred on it is
   kept, as the default path would. So is a peak next to missing data (a masked
   run, a beamstop at either end, or a single nonfinite point, which blocks up to
   64 points on each side for the widest kernel), where
   the wide kernels that would see a broad peak are not defined.
4. Add local maxima of `z` (the max over the widths defined at a point) with
   `z ≥ k` (default 20). Each must be:
   - more than 9 points from every other added peak (strongest first) and
     from every kept persistence peak;
   - at or below the same q cut as the persistence path (`q_trim_high`, 0.07);
   - a local maximum at its own best width, with both neighbours defined;
   - not in a valley at a finer width: at every narrower width, `z` within ±1
     point must be above −3.

   The own-width check matters where a wide width becomes defined (at the trace
   ends or next to a masked run): there the max over widths steps up without a
   peak. The valley check rejects the midpoint of two close peaks, which a wide
   kernel sees as one blob while the narrower kernels sit in the valley between
   them.

   Known limitation: a strongly asymmetric peak (one side 8 or more times
   wider than the other) can get a second, added peak about 10 points out on
   its long side, because the long side is concave at the wide widths and not
   a valley at the narrow ones. Rejecting it would need a shape model per peak,
   which this detector avoids.

**Output.** Same fields as the default path, plus `z`. For an added peak,
`prominence` is the persistence prominence of `I` at that index if it is a
local maximum and `0.0` otherwise (a shoulder has none). `sharpness` is the
Savitzky-Golay sharpness at the index; a value ≤ 0 (common for shoulders)
is replaced by the median of the positive sharpness values of the returned
peaks. `score` takes the CV of sharpness and falls back to CV = 0 when the
mean is ≤ 0; a near-zero substitute would inflate the CV and demote exactly
the peaks the hybrid adds, and the median keeps them neutral.

Without `σ`, `findpeaks(q, I; method = :hybrid)` estimates
`σ² = g·max(I, 0)`, clamped from below at a noise scale `s` taken from second
differences (`1.4826·MAD/√6`, then their standard deviation over `√6` for
quantised counts, then `mean(|I|)`, then 1). The gain `g` is the median of
`(d²/6)/I` over points with `I > 3s`, divided by 0.455. All of `g`, `s` and `I`
scale with the intensity units. Both noise scales are estimated from the trace
itself, so only the relative shape of `σ` matters.

**Why this use of σ is not the rejected "σ-based SNR".** Nothing is fitted
per peak, and `σ` is never trusted as an absolute noise level. It only
weights the filter. Both `z` and `zraw` are rescaled by noise scales measured on
the trace's own data: the clipped spread of each width's responses, and of the
second differences.

**Evidence.** On a set of curated real traces it recovered several times
more user-added peaks than the default path while detecting fewer
user-rejected ones, and the phase calls from the resulting peak lists improved. This is not enough to change the default; it is opt-in until
it has been tried on more data.

## References

[1] Huber, S. (2021). [Persistent Topology for Peak Detection](https://www.sthu.org/blog/13-perstopology-peakdetection/index.html). *sthu.org*. — Accessible explanation of the flooding intuition applied specifically to 1D peak detection.

[2] Edelsbrunner, H., Letscher, D., & Zomorodian, A. (2002). Topological Persistence and Simplification. *Discrete & Computational Geometry*, 28, 511–533. https://doi.org/10.1007/s00454-002-2885-2 — Foundational paper introducing persistent homology.

## Intentionally out of scope

- **Parametric shape fitting.** Brittle on real SAXS backgrounds; phase-ratio coherence is a stronger downstream filter.
- **Sub-pixel peak positions.** Peaks are returned at grid positions. Parabolic interpolation is a straightforward follow-up if needed.
- **Background subtraction in the default path.** SNIP is used only by the opt-in hybrid detector.
- **Multi-trace batch statistics.** Thresholds are per-trace; no pooling across scans.
- **q-banded kneedle** (independent thresholds per q-band). Considered when designing the no-peak-trace gating; rejected because every band of a noise-only trace is itself noise, and per-band kneedle would still find spurious knees.
- **σ(q)-based noise floor** using the third column of `_tot.dat`. Considered for the same purpose; rejected because parametric noise modelling at the candidate-selection stage re-opens the design dimension this document already closed (see "Why not σ-based SNR" above). The relative-prominence floor uses the data's own prominence distribution and avoids this.
- **Low-q trim.** The empirical failure mode is high-q only. A symmetric trim risks clipping Lamellar 1st-order peaks; we wait for a fixture exhibiting low-q noise contamination before adding one.
