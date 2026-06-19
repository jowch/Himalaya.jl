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

Instead, R² is surfaced **informationally** in the UI: the comb/residual chart renders a `fit R²` readout (`src/print/comb/ResidualChart.tsx`) so users can see the lattice-fit quality, but it is not used as a hard gate. There is no `r_squared < 0.98` dimming threshold and no "low R²" label — the candidate/assignment UI (`AssignmentRail` / `CandidateRow`) ranks on `score` alone. This keeps the score formula clean and ensures R² doesn't dominate over the physically-motivated coverage and consistency signals.

## Gauss–Bonnet coexistence flag — not in the score

When a bicontinuous cubic (Pn3m / Im3m / Ia3d) is in an exposure's assignment, the `bonnet` field on the indices response (`GET /api/exposures/{id}/indices`) flags candidates whose measured lattice matches what the Gauss–Bonnet ratio predicts for a *coexisting* cubic of a different phase (`a_Pn3m : a_Im3m : a_Ia3d = 1.000 : 1.279 : 1.576`; kernel in `src/bonnet.jl`). It is a **display-and-ranking affordance computed from the live assignment**, recomputed per request and **never persisted, never folded into `score`**. Folding it in would corrupt the `auto_group` / `remove_subsets` ordering, which relies on `score` being coverage×consistency alone — so the Bonnet match surfaces only as the ⭙ badge (and may sort up within the candidate-list view), never as a score change.

## Peak-to-ratio assignment is greedy, not optimal

Before `score` ever runs, `indexpeaks` has to decide *which* observed peak fills
*which* slot in a phase's ratio series. When several observed peaks fall within
`tol` of the same ratio — or one peak sits within tolerance of two adjacent
ratios — there is a conflict, and the conflict is resolved greedily.

The implementation (`src/index.jl:168-180`) collects every `(ratio, peak)`
match whose residual is under tolerance, sorts the matches by ascending
residual, and walks them best-first. A match is committed only if **both** its
ratio slot and its peak are still free; once either is claimed it is locked:

```julia
rs, ps, δs = findnz(residuals)
used_ratios = Set()
assign = view(assignments, :, i)

# go through matches in order of increasing error
# only do assignment if ratio and peak haven't already been assigned
for j = sortperm(δs)
    r, p = rs[j], ps[j]
    if r ∉ used_ratios && assign[p] == 0
        assign[p] = r
        push!(used_ratios, r)
    end
end
```

This produces a 1:1 matching (no ratio is filled twice, no peak is used twice),
but it is **greedy, not optimal** — it is not the Hungarian/Kuhn–Munkres
assignment that minimises total residual. Greedy and optimal can disagree: a
slightly-offset strong peak can win a slot by having the single smallest
residual, even when the geometrically-correct peak would have given a better
*overall* matching once the knock-on assignments are accounted for. The
displaced peak then either lands in a worse slot or finds none free.

The failure this can cause is quiet. After assignment, bases that do not clear
`minpeaks(P)` are dropped (`src/index.jl:184-185`); if a greedy mis-assignment
costs an index even one slot, an otherwise-valid index can fall below the
minimum and be **silently discarded** — it never reaches `score`, never appears
in the candidate list, and there is no diagnostic. This is rare in practice
(tolerances are tight and real phase ratios are well-separated), but it is the
mechanism to suspect when a phase you expect is simply absent rather than
ranked low. Switching to an optimal assignment is the obvious remedy if a
fixture ever exhibits it.

## Design constraints

Several design choices are non-obvious enough to be worth stating explicitly:

- **No unit mixing.** Coverage and consistency are both dimensionless ratios in `(0, 1]`. Their product is always interpretable.
- **No gap count.** Gaps are already captured by coverage: a missing rank costs `1/rank` units of coverage, weighted by position. An explicit gap-penalty term would double-count this and introduce a free hyperparameter.
- **Score never reaches exactly 0.** Both `1/(1+CV)` and harmonic coverage approach 0 asymptotically. This matters for `remove_subsets`, which relies on score ordering to break ties — a zero score would make two bad indices indistinguishable.

## Further reading

- [docs/peak-finding.md](peak-finding.md) — how sharpness is computed by `findpeaks`.
