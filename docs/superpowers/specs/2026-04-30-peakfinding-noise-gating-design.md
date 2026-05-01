# Peak-finding noise gating: high-q trim + relative-prominence floor

**Status:** spec
**Date:** 2026-04-30
**Touches:** `src/peakfinding.jl`, `test/peakfinding_real.jl`, `docs/peak-finding.md`

## Problem

`findpeaks` returns spurious peaks on traces dominated by a form-factor envelope and pure shot noise — i.e., traces with no real Bragg peaks. The regression fixture `test/data/form-factor_tot.dat` is exactly such a trace; the existing test allows up to 2 spurious peaks (`SPURIOUS_CEILING["form-factor_tot.dat"] = 2`) and we currently produce exactly 2.

Empirically, both spurious peaks sit at the extreme high-q tail (`qfrac ≥ 0.99`) and have prominences (1.4 and 3.6) that are 10–25× the candidate-prominence median. By contrast, every real peak in the `example` and `cubic` fixtures has prominence at least 52× the candidate median.

The root cause is that the kneedle elbow finder, run on the sorted prominence distribution of a trace with no real peaks, picks an arbitrary knee in what is effectively a noise decay curve. Kneedle has no way to say "this trace has no real peaks." A trace-level absence-of-signal failure produces per-peak false positives.

## Goal

Drop `SPURIOUS_CEILING["form-factor_tot.dat"]` from 2 to 0 without regressing `RECALL_FLOOR` on `example` (7/10) or `cubic` (8/23).

## Non-goals

- q-banded kneedle (independent thresholds per q-band). Considered; rejected because every band of `form-factor` is noise and per-band kneedle would still find spurious knees.
- σ(q)-based noise floor (using the third column of `_tot.dat`). Considered; rejected because parametric noise modelling re-opens the design dimension that `docs/peak-finding.md` explicitly closed.
- Low-q trim. The empirical failure mode is high-q only. A low-q trim risks clipping Lamellar 1st-order peaks that we have no fixture for.
- Sharpness relative-floor. Sharpness shows the same separation pattern empirically, but adding a second adaptive gate widens the tuning surface without an observed failure that requires it.
- Sub-pixel positions, background subtraction, parametric shape fitting — already out of scope per `docs/peak-finding.md`.

## Design

Two changes to `findpeaks`, both adaptive (no fixed instrument-dependent constants), both layered on top of the existing prominence + sharpness AND-gate.

### A. High-q trim (post-filter)

Discard returned peaks whose `q` exceeds `qmin + (1 - q_trim_high) * (qmax - qmin)`.

- New kwarg: `q_trim_high::Real = 0.05`.
- `q_trim_high = 0.0` disables the trim (full backward compatibility for callers that pass it explicitly).
- The trim happens **after** kneedle, as a post-filter on the kept peaks. This is load-bearing: trimming the trace before kneedle removes the high-q noise candidates that anchor the lower part of the sorted-prominence curve, which shifts the knee upward and rejects real peaks. Empirically we observed `example_tot.dat` recall drop from 7 to 6 with pre-kneedle trimming; the post-filter ordering preserves recall while still discarding spurious high-q peaks.

**Defaults rationale.** Real peaks across our fixtures sit at `qfrac ∈ [0.10, 0.34]`. `0.05` clears the observed real-peak band by a factor of >10 in margin while covering the observed spurious peaks (`qfrac ≥ 0.99`).

### B. Relative-prominence floor

After kneedle but before returning, additionally require:

```
prominence ≥ prom_ratio_floor * median(candidate_prominence)
```

where `candidate_prominence` is the prominence of *all* candidates produced by `persistence` on the (trimmed) trace, not just the kneedle-survivors.

- New kwarg: `prom_ratio_floor::Real = 15.0`.
- `prom_ratio_floor = 0.0` disables the floor.
- The effective lower bound on prominence is `max(resolved_prom_floor, prom_ratio_floor * median_cand_prom)`, where `resolved_prom_floor` is whatever `something(prom_floor, knee(...))` produces today. The ratio floor is therefore always combined with — not bypassed by — a manual `prom_floor`. A user who wants to opt out passes `prom_ratio_floor = 0.0` explicitly.
- The floor is **skipped** when the number of candidates is below `RATIO_FLOOR_MIN_CANDIDATES` (currently 20). At low candidate counts (e.g., synthetic single-peak traces) the median is dominated by the peaks themselves, so the ratio is ~N× higher than the peak being measured against and would suppress real signal. The minimum count is high enough that all real-data fixtures (≥ 141 candidates) trigger the gate, and low enough that synthetic single-peak unit tests do not.

**Defaults rationale.** On our fixtures:
- `example`: minimum kept-peak prominence ratio = 52
- `cubic`: minimum kept-peak prominence ratio = 243
- `form-factor`: maximum kept-peak prominence ratio = 24.8

The default was originally `30.0` (chosen to sit between form-factor's worst spurious at 24.8 and example's best real peak at 52). It was lowered to `15.0` after a real SSRL Pn3m trace was found where genuine Bragg peaks had prominence ratios in the 17–31 range — `30.0` was killing real peaks. With `15.0` the backstop now sits *below* the form-factor noise band (24.8), so on that fixture the high-q trim (`q_trim_high = 0.05`) carries the spurious-rejection load alone; the backstop fires only on noise traces whose candidates concentrate in the lower-q region the trim doesn't touch. Real-peak fixtures remain a no-op (kneedle threshold dominates).

## Failure modes to consider

- **A trace with very few candidates total (< 20).** The relative-prominence floor is skipped entirely in this regime, deferring to kneedle alone. Synthetic test cases with one or two clean peaks land here; on real-data traces we have never seen fewer than 141 candidates, so this carve-out doesn't affect production behaviour.
- **A real peak between qfrac 0.95 and 1.0.** Trim removes it. We have no fixture exhibiting this. If a real-world case appears, the `q_trim_high` kwarg is a one-line override at the call site (or a config knob in `experiment.toml` later, but that change is out of scope here).
- **`prom_ratio_floor` interacting with manual `prom_floor`.** Both are lower bounds; `findpeaks` takes their max. Manual `prom_floor` cannot drop below the ratio floor unless the caller also sets `prom_ratio_floor = 0.0`. This is intentional — manual `prom_floor` is a *trust the data* override; `prom_ratio_floor` is a *don't trust pure-noise traces* backstop. They serve different purposes and should compose, not cancel.

## Test plan

`test/peakfinding_real.jl` is updated:

1. `SPURIOUS_CEILING["form-factor_tot.dat"]` drops from 2 to 0.
2. `RECALL_FLOOR["example_tot.dat"]` stays at 7.
3. `RECALL_FLOOR["cubic_tot.dat"]` stays at 8.
4. Comment updates: drop the "current: 2 spurious" annotation.

Additional unit-level tests in a new `test/peakfinding_gating.jl`:

- `q_trim_high = 0.0` reproduces pre-change behaviour bit-for-bit.
- `q_trim_high = 0.5` returns no peaks at qfrac ≥ 0.5 in any fixture.
- `prom_ratio_floor = 0.0` reproduces pre-change behaviour bit-for-bit.
- `prom_ratio_floor = 1e6` returns zero peaks on all fixtures (extreme floor sanity check).
- Manual `prom_floor` and `prom_ratio_floor` compose as a max: setting `prom_floor = 0.0` alone does not reproduce pre-change behaviour; setting both to `0.0` does.

## Documentation updates

`docs/peak-finding.md`:

- Document both kwargs in the defaults section.
- Add a short subsection on "no-peak traces" explaining why a relative-prominence floor exists alongside kneedle.
- Update the "Intentionally out of scope" list: q-banded kneedle and σ-based floor get added with one-line rationales.

## Backwards compatibility

Both kwargs default to non-zero values, so existing callers will see behaviour change. The only known semantic change is on traces resembling `form-factor_tot.dat` (no real peaks, candidates concentrated in noise). Traces with real peaks are unaffected — the kneedle threshold is ~50× higher than `prom_ratio_floor * median` on every real-peak fixture we have.

A caller wanting strictly pre-change behaviour can pass `q_trim_high = 0.0, prom_ratio_floor = 0.0`.

## Implementation notes

- The trim happens in `findpeaks`, not in `persistence`/`sharpness`/`kneedle`. Those primitives stay general-purpose.
- `median(candidate_prominence)` runs in `persistence`-output ordering, not sorted. `Statistics.median` is already a transitive dep via `Peaks.jl`, no new dep.
- A regression test scaffold for "kept_prom / median_cand_prom ratio" can be added later if we want to monitor that quantity over time, but it's not on the critical path.
