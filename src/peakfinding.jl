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
- `method`: `:persistence` (default, described here) or `:hybrid`, which adds
  a matched-filter detector (see "Hybrid" below).
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
  as a post-filter (applied after kneedle). Default `0.05` (`0.07` for
  `method = :hybrid`). Set to `0.0`
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

# Hybrid (`method = :hybrid`)
Opt-in. `k` and `kb` apply only to `:hybrid` and are ignored otherwise;
`normalize_by_σ` affects only the persistence part (the matched filter
always works on `I`). Runs the persistence detector with `q_trim_high = 0.07`,
keeps its peaks whose matched-filter z (max within ±3 points) is at least `kb`
(default 5; a peak where z is undefined because even the narrowest kernel runs
off the trace is kept), then adds local maxima of z with `z ≥ k` (default 12),
more than 9 points apart and more than 9 points from a kept peak, at or below
the same q cut as the persistence path. Added peaks get `prominence` = the
persistence prominence of `I` there if it is a local maximum, else `0.0`. A
returned peak with `sharpness` ≤ 0 gets the median of the positive sharpness
of the returned peaks (fallback: the smallest positive value in
`sharpness_full`, then `eps()`). `score` falls back to CV = 0 when the mean
sharpness is ≤ 0, and a near-zero value would inflate the CV and demote the
peak; the median keeps it neutral. The result has an extra field `z`
(max-over-width matched-filter z, `NaN` where undefined). If the valid range
of the trace is shorter than 150 points the persistence result (with the
hybrid's `q_trim_high`) is returned and `z` is all `NaN`. Without `σ`
(`findpeaks(q, I; method = :hybrid)`) σ is estimated as
`sqrt(max(I, 0)·s)` clamped from below at the data noise scale `s`; both
terms scale with `I`, and the per-width MAD normalisation removes the overall
scale, so only the relative shape matters. See `docs/peak-finding.md`.

# Returns
`(; indices, q, prominence, sharpness, sharpness_full)` —
`indices`, `q`, `prominence`, `sharpness` are four equal-length vectors,
sorted by ascending q. `indices` reference positions in the original
input vectors (the trim is invisible to callers). `sharpness_full` is
the per-sample sharpness vector over the *entire* trace
(`length(sharpness_full) == length(I)`); callers can reuse it to avoid
re-running the Savitzky–Golay pass.
"""
const RATIO_FLOOR_MIN_CANDIDATES = 20

function _persistence_findpeaks(q, I, σ; normalize_by_σ, sharpness_method,
                                prom_floor, sharp_floor, q_trim_high, prom_ratio_floor)
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

    result = (indices        = idx[perm],
              q              = q[idx[perm]],
              prominence     = cands.prominence[keep][perm],
              sharpness      = sharps_at_peaks[keep][perm],
              sharpness_full = sharps_full)
    result, cands
end

# Hybrid detector: persistence detector (gated by matched-filter z) plus
# matched-filter peaks. See src/matchedfilter.jl and docs/peak-finding.md.
function _hybrid_findpeaks(q, I, σ; k, kb, q_trim_high, sep = 9, kwargs...)
    base, cands = _persistence_findpeaks(q, I, σ; q_trim_high, kwargs...)
    n = length(I)
    z = matched_filter_z(I, σ)
    # Too short a valid range: the matched filter is undefined, so return the
    # persistence result (with the hybrid's q_trim_high) unchanged.
    all(isnan, z) && return (; base..., z = z)
    qmin, qmax = extrema(q)
    qcut = qmin + (1 - q_trim_high) * (qmax - qmin)
    allowed = q .<= qcut
    idx = hybrid_merge(base.indices, z, k, kb, sep, allowed)

    sharps_full = base.sharpness_full
    sharp = sharps_full[idx]
    # `score` takes the CV of sharpness; a near-zero substitute inflates it and
    # demotes the peak, so use the median of the positive values instead.
    pos = filter(>(0), sharp)
    fill_value = if !isempty(pos)
        median(pos)
    else
        tiny = Inf
        for s in sharps_full
            s > 0 && s < tiny && (tiny = s)
        end
        isfinite(tiny) ? tiny : eps()
    end
    for m in eachindex(sharp)
        sharp[m] > 0 || (sharp[m] = fill_value)
    end
    pmap = Dict(zip(cands.indices, cands.prominence))
    prom = [get(pmap, i, 0.0) for i in idx]
    (indices = idx, q = q[idx], prominence = prom, sharpness = sharp,
     sharpness_full = sharps_full, z = z)
end

function findpeaks(q, I, σ; method = :persistence,
                              normalize_by_σ    = false,
                              sharpness_method = :savgol,
                              prom_floor       = nothing,
                              sharp_floor      = nothing,
                              q_trim_high      = nothing,
                              prom_ratio_floor = 15.0,
                              k                = 12.0,
                              kb               = 5.0)
    if method === :persistence
        return first(_persistence_findpeaks(q, I, σ; normalize_by_σ, sharpness_method,
            prom_floor, sharp_floor, q_trim_high = something(q_trim_high, 0.05),
            prom_ratio_floor))
    elseif method === :hybrid
        return _hybrid_findpeaks(q, I, σ; k, kb,
            q_trim_high = something(q_trim_high, 0.07), normalize_by_σ, sharpness_method,
            prom_floor, sharp_floor, prom_ratio_floor)
    end
    throw(ArgumentError("Unknown method: $(repr(method)) (use :persistence or :hybrid)"))
end

# Without σ, the default path uses a constant σ. The hybrid path needs the
# relative noise shape, so it uses the Poisson shape sqrt(max(I, 0) · s) with
# the data noise scale s, clamped from below at s. Both terms scale with I, so
# the estimate is independent of intensity units, and it does not understate
# noise where I is near zero or negative. The per-width MAD normalisation removes the
# overall scale of σ.
function _estimate_sigma(I)
    s = noise_scale(I, Float64[])
    [s * max(isfinite(x) && x > 0 ? sqrt(x / s) : 0.0, 1.0) for x in I]
end

findpeaks(q, I; kwargs...) =
    findpeaks(q, I, get(kwargs, :method, :persistence) === :hybrid ?
                    _estimate_sigma(I) : ones(length(I)); kwargs...)
