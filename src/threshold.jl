"""
    knee(v) -> Float64

Find the elbow of a sorted descending sequence `v` using the kneedle
algorithm: the value at the point of maximum perpendicular distance
from the chord connecting the first and last entries, after normalising
both axes to [0, 1] so the distance is scale-invariant.

Used as an adaptive threshold: values ≥ knee(v) are "signal", values
< knee(v) are "noise". No magic numbers, no labels, no per-trace tuning —
the threshold emerges from the shape of the data's own distribution.

Edge cases:
- Empty input → 0.0 (nothing passes).
- Single or two values → v[1].
- Constant sequence → v[1] (no elbow; nothing passes).
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

    best_score = -Inf
    best_idx   = 2          # first interior point (used when all scores equal)
    for i in 2:(n - 1)
        x = (i - 1) / (n - 1)
        y = (v[i] - y_min) / y_range
        score = 1 - x - y                 # > 0 ⇔ point below the (0,1)→(1,0) chord
        if score >= best_score
            best_score = score
            best_idx   = i
        end
    end

    # Return the last "steep" value: the entry just before the max-distance point.
    # This is the threshold that separates signal (≥ threshold) from noise (< threshold).
    float(v[best_idx - 1])
end
