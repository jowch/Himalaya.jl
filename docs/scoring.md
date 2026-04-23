# Index Scoring

This document explains why `score(index)` is defined the way it is and what physical intuitions it encodes.

## Why scoring matters

`score` is load-bearing in two places:

1. **`auto_group`** (pipeline) — greedily selects the non-overlapping set of indices to present as the auto-analysis result by sorting candidates descending by score and claiming peaks first-come first-served. The highest-scoring index gets to claim its peaks before lower-scoring alternatives, so a wrong ordering can lock out the correct assignment.

2. **`remove_subsets`** — when two indices share the same basis but one's peaks are a subset of the other's, the lower-scoring one is removed. Score breaks the tie.

Getting the ordering right matters more than the absolute values. False negatives (correct phase ranked low) are the dominant failure mode — the user has to manually promote the right index from the Alternatives list.

## The formula

```
score(index) = coverage(index) × consistency(index)
```

Both factors are in `(0, 1]`. Their product is in `(0, 1]`.

### Coverage

```
coverage = Σ(1/rank_i  for each found peak i)
         / Σ(1/rank_j  for all expected peaks j)
```

`rank` is the 1-based ordinal position in the phase's ratio series — rank 1 is the fundamental reflection, rank 2 is the second, and so on.

**Why harmonic weighting?** Peaks appear in order of increasing q: lower-order reflections are stronger and are observed first as signal improves. If you see the 3rd peak you expect to also see the 1st and 2nd. Harmonic weights (`1/rank`) encode this: missing the fundamental costs `1/1 = 1.0` unit out of the denominator, while missing the 8th peak costs only `1/8 = 0.125`. An index that matches only the first three peaks of a 14-peak phase scores higher than one that matches three scattered high-order peaks.

The denominator is fixed per phase (sum of `1/rank` over all expected peaks), so coverage is always normalised to `[0, 1]` regardless of how many peaks a phase predicts.

### Consistency

```
cv = std(sharpnesses) / mean(sharpnesses)   # coefficient of variation
consistency = 1 / (1 + cv)
```

where `sharpnesses` are the per-peak sharpness values at the assigned ratio positions. When only one peak is assigned (CV undefined), `consistency = 1.0`.

**Why sharpness CV?** Peaks from the same liquid-crystalline phase have the same underlying origin and therefore similar line shapes: wide peaks come with wide peaks, sharp peaks with sharp peaks. An assignment that mixes a very sharp peak with a very broad one is suspicious — it may be combining peaks from different phases or misassigning noise. CV is unit-free and ranges from 0 (identical sharpness) upward. The `1/(1+CV)` transform maps this to `(0, 1]`: CV = 0 gives consistency 1.0 (perfect), CV = 1 gives 0.5, and large CV approaches 0.

**Why not `clamp(1 - CV, 0, 1)`?** That formula collapses to 0 whenever CV > 1, which is common (sharpness can vary by an order of magnitude across peaks in a real trace). A score of exactly 0 interacts badly with `remove_subsets`: if two related indices both score 0, neither is identified as a subset of the other, so both survive and pollute the candidate list. `1/(1+CV)` degrades gracefully to near-zero for very heterogeneous assignments without ever reaching it.

**NaN guard.** If all sharpness values in an index are exactly 0, `std/mean = 0/0 = NaN`. The implementation guards this:

```julia
cv = if length(sharps) > 1 && mean(sharps) > 0
    std(sharps) / mean(sharps)
else
    0.0
end
```

Zero-sharpness falls back to `cv = 0` (treated as maximally consistent), degrading to coverage-only scoring.

## R² — not in the score

R² from the least-squares lattice fit is stored per index in the database but is **not part of `score`**. In practice, any assignment Himalaya produces has very high R² — the fitting step is constrained enough that poor fits are rare. R² does not discriminate between competing good assignments.

Instead, R² acts as a **hard gate in the UI**: alternatives with `r_squared < 0.98` are visually dimmed in PhasePanel with a "low R²" label. Users can still manually promote them. This keeps the score formula clean and ensures R² doesn't dominate over the physically-motivated coverage and consistency signals.

## Design constraints

Several design choices are non-obvious enough to be worth stating explicitly:

- **No unit mixing.** Coverage and consistency are both dimensionless ratios in `(0, 1]`. Their product is always interpretable.
- **No gap count.** Gaps are already captured by coverage: a missing rank costs `1/rank` units of coverage, weighted by position. An explicit gap-penalty term would double-count this and introduce a free hyperparameter.
- **Score never reaches exactly 0.** Both `1/(1+CV)` and harmonic coverage approach 0 asymptotically. This matters for `remove_subsets`, which relies on score ordering to break ties — a zero score would make two bad indices indistinguishable.

## Further reading

- [docs/peak-finding.md](peak-finding.md) — how sharpness is computed by `findpeaks`.
