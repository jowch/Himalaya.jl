"""
    findpeaks(q, I, σ; normalize_by_σ = false,
                        sharpness_method = :savgol,
                        prom_floor  = nothing,
                        sharp_floor = nothing) -> NamedTuple
    findpeaks(q, I;    kwargs...) -> NamedTuple

Detect peaks in a 1D signal using topological prominence (criterion A) and
curvature/sharpness (criterion B), each adaptively thresholded via the
kneedle elbow finder. Shape-agnostic; no per-trace tuning.

# Arguments
- `q`, `I`, `σ`: equal-length vectors. `σ` is per-point intensity
  uncertainty (e.g., the third column of a SAXS `_tot.dat` file).
- `normalize_by_σ`: if `true`, work on `I ./ σ` so prominence is in
  noise-relative units. Default `false`. For Poisson-dominated SAXS data
  where σ ≈ √I, normalising by σ compresses the very dynamic range we
  are trying to detect — peaks lose their prominence advantage over
  noise. The kwarg is kept for data with approximately homoscedastic
  noise (e.g., detector electronic noise dominating).
- `sharpness_method`: `:savgol` (default, single-scale 2nd derivative) or
  `:cwt` (multi-scale Ricker max response).
- `prom_floor`, `sharp_floor`: optional manual thresholds. When `nothing`
  (default), the kneedle algorithm chooses each from the data's own
  distribution.

# Returns
`(; indices, q, prominence, sharpness)` — four equal-length vectors,
sorted by ascending q.
"""
function findpeaks(q, I, σ; normalize_by_σ    = false,
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

    # Sort surviving candidates by ascending q (== ascending index).
    idx  = cands.indices[keep]
    perm = sortperm(idx)

    (indices    = idx[perm],
     q          = q[idx[perm]],
     prominence = cands.prominence[keep][perm],
     sharpness  = sharps_at_peaks[keep][perm])
end

findpeaks(q, I; kwargs...) = findpeaks(q, I, ones(length(I)); kwargs...)
