using Peaks: findmaxima, peakproms

"""
    persistence(y) -> (; indices, prominence)

Compute the topological persistence — equivalently, the prominence — of
every local maximum in the 1D signal `y`.

Intuition: imagine flooding the signal from above. Each local max is
born when the water level passes it and dies when its hill merges with
a taller neighbour (or hits the array boundary). The vertical distance
between birth and death is its prominence.

This is the cleanest mathematical formalisation of criterion (A): a
peak is real if it stands clearly above the local floor. Internally
delegates to `Peaks.peakproms`, which implements the same computation
with careful handling of edges and plateaus.
"""
function persistence(y)
    maxima = findmaxima(y)
    with_proms = peakproms(maxima)
    (indices = with_proms.indices, prominence = with_proms.proms)
end
