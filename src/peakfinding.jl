using Statistics: median

"""
    findpeaks(q, I, σ; normalize_by_σ = false,
                        sharpness_method = :savgol,
                        prom_floor       = nothing,
                        sharp_floor      = nothing,
                        q_trim_high      = 0.05,
                        prom_ratio_floor = 15.0) -> NamedTuple
    findpeaks(q, I;    kwargs...) -> NamedTuple

Detect peaks in a 1D signal using topological prominence (criterion A) and
curvature/sharpness (criterion B), each adaptively thresholded via the
kneedle elbow finder, plus two backstops against pure-noise traces.
Shape-agnostic; no per-trace tuning.

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
  distribution. Manual values compose with `prom_ratio_floor` as a max.
- `q_trim_high`: fraction of the q-range to discard from the high-q end
  as a post-filter (applied after kneedle). Default `0.05`. Set to `0.0`
  to disable. Real Bragg peaks rarely land in the upper few percent of
  the swept q-range on typical SAXS configurations; the high-q tail is
  dominated by shot noise on a near-zero baseline and produces spurious
  candidates. The post-filter ordering is load-bearing — truncating
  before kneedle removes the noise candidates that anchor the
  sorted-prominence curve and shifts the knee upward, regressing recall
  on real fixtures.
- `prom_ratio_floor`: relative prominence floor as a multiple of the
  candidate-prominence median. Default `15.0`. Set to `0.0` to disable.
  Combined with the kneedle threshold (and any manual `prom_floor`) as a
  max, this acts as a backstop on traces where kneedle's elbow falls
  inside a noise-only distribution. Skipped when fewer than
  `RATIO_FLOOR_MIN_CANDIDATES` (currently 20) candidates are present, so
  synthetic single-peak traces aren't suppressed by their own median.

# Returns
`(; indices, q, prominence, sharpness)` — four equal-length vectors,
sorted by ascending q. `indices` reference positions in the original
input vectors (the trim is invisible to callers).
"""
const RATIO_FLOOR_MIN_CANDIDATES = 20

function findpeaks(q, I, σ; normalize_by_σ    = false,
                              sharpness_method = :savgol,
                              prom_floor       = nothing,
                              sharp_floor      = nothing,
                              q_trim_high      = 0.05,
                              prom_ratio_floor = 15.0)
    @assert length(q) == length(I) == length(σ) "q, I, σ must be equal length"
    @assert 0 <= q_trim_high < 1 "q_trim_high must be in [0, 1)"
    @assert prom_ratio_floor >= 0 "prom_ratio_floor must be ≥ 0"

    y = normalize_by_σ ? I ./ σ : I

    cands           = persistence(y)
    sharps_full     = sharpness(y; method = sharpness_method)
    sharps_at_peaks = sharps_full[cands.indices]

    pf_kneedle = something(prom_floor,  knee(sort(cands.prominence; rev = true)))
    sf         = something(sharp_floor, knee(sort(sharps_at_peaks;  rev = true)))

    # Relative-prominence backstop. Only applied when the candidate
    # population is large enough for `median` to track the noise floor
    # rather than the peaks themselves — otherwise (e.g. on synthetic
    # single-peak traces) the floor can suppress the very peak it should
    # be measuring against.
    pf_ratio = (prom_ratio_floor > 0 &&
                length(cands.prominence) >= RATIO_FLOOR_MIN_CANDIDATES) ?
               prom_ratio_floor * median(cands.prominence) : 0.0
    pf = max(pf_kneedle, pf_ratio)

    keep = (cands.prominence .>= pf) .& (sharps_at_peaks .>= sf)

    # Sort surviving candidates by ascending q (== ascending index).
    idx  = cands.indices[keep]
    perm = sortperm(idx)

    # Apply the high-q trim as a post-filter, so kneedle continues to see
    # the full candidate distribution (truncating the high-q noise tail
    # before kneedle would shift the knee upward and reject real peaks).
    if q_trim_high > 0 && !isempty(idx)
        qmin, qmax = extrema(q)
        qcut = qmin + (1 - q_trim_high) * (qmax - qmin)
        in_window = q[idx[perm]] .<= qcut
        perm = perm[in_window]
    end

    (indices    = idx[perm],
     q          = q[idx[perm]],
     prominence = cands.prominence[keep][perm],
     sharpness  = sharps_at_peaks[keep][perm])
end

findpeaks(q, I; kwargs...) = findpeaks(q, I, ones(length(I)); kwargs...)
