# Peak finding — design notes

> Narrative reference for the `findpeaks` algorithm. For the formal
> specification, see
> [specs/2026-04-22-persistence-peakfinding-design.md](superpowers/specs/2026-04-22-persistence-peakfinding-design.md).
> For history, see the superseded
> [v1 spec](superpowers/specs/2026-04-22-peakfinding-rewrite-design.md).

## What `findpeaks` does (v2)

`findpeaks(q, I, σ)` is a shape-agnostic detector with two independently
thresholded criteria:

1. **Topological persistence** (prominence) — "does the peak stand above
   the local floor?" Implemented as a thin wrapper over `Peaks.peakproms`.
2. **Curvature / sharpness** via Savitzky-Golay 2nd derivative (default)
   or multi-scale Ricker CWT — "is the peak sharp, not a broad hump?"

Both thresholds are set *adaptively per trace* via the kneedle
algorithm on the data's own distribution. No per-trace tuning.

## Non-obvious defaults

- **`normalize_by_σ = false`.** For Poisson-dominated SAXS data, σ ≈ √I,
  so `I/σ ≈ √I` — normalising *compresses* the dynamic range we need
  to detect peaks. Real Bragg peaks lose their prominence advantage over
  noise when normalised. The kwarg exists for homoscedastic data (where
  σ is approximately constant), but default to `false` for SAXS.
- **`sharpness_method = :savgol`.** Faster, single-scale. Use `:cwt`
  when peaks vary significantly in width within one trace.

## Intentionally out of scope

- **Parametric shape fitting** (Lorentzian / Gaussian / Voigt). v1 tried
  this and capped at ~60% recall on real traces. v2 drops it.
- **σ-based SNR scoring.** Phase-ratio coherence (via `indexpeaks`)
  is a stronger filter than per-peak σ tests — a single 5σ excess can
  be noise, but four peaks at Pn3m's √2 : √3 : √4 : √6 ratios cannot.
- **Sub-pixel peak positions.** Peaks are returned at grid positions.
  Parabolic interpolation would be an easy follow-up if needed.
- **Background subtraction** (SNIP, asymmetric least squares). Noted as
  a potential follow-up if steep-background peaks become important.

## Version history

- `v0.5.0` (on `peakfinding-rewrite` branch, not yet merged) — v2
  persistence + sharpness detector.
- `v0.4.x` (current `main`) — v1 prominence-threshold + Savitzky-Golay
  detector requiring per-trace θ tuning.
