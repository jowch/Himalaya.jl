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
always works on `I`). Runs the persistence detector with `q_trim_high = 0.07`
and keeps its peaks unless the matched-filter response judged against the
second-difference noise scale (max over widths within ±3 points) is defined at
the peak, below `kb` (default 5), and not cut short by a masked run. Then adds
local maxima of the clipped matched-filter z with `z ≥ k` (default 20), more
than 9 points apart and from a kept peak, at or below the same q cut, that are
maxima at their own width and not valleys at a finer width. See
`src/matchedfilter.jl` and `docs/peak-finding.md` for the masking, the two noise
scales and the checks. Added peaks get `prominence` = the persistence
prominence of `I` there if it is a local maximum, else `0.0`. A returned peak
with `sharpness` ≤ 0 gets the median of the positive sharpness of the returned
peaks (fallback: the smallest positive value in `sharpness_full`, then
`eps()`). `score` falls back to CV = 0 when the mean sharpness is ≤ 0, and a
near-zero value would inflate the CV and demote the peak; the median keeps it
neutral. The result has an extra field `z` (max-over-width clipped z, `NaN`
where undefined or missing). If the trace has fewer than 150 points the
persistence result (with the hybrid's `q_trim_high`) is returned and `z` is
all `NaN`. Without `σ` (`findpeaks(q, I; method = :hybrid)`) σ is estimated as
`sqrt(g·max(I, 0))` clamped from below at a second-difference noise scale, with
the gain `g` from the same differences; everything scales with `I`, and both
noise scales are measured on the trace, so only the relative shape matters.

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
    # Too short for the widest kernel to give a usable scale: return the
    # persistence result (with the hybrid's q_trim_high) unchanged.
    n < MF_MIN_LENGTH && return (; base..., z = fill(NaN, n))
    Z, zraw, blocked = mf_width_z(I, σ)
    qmin, qmax = extrema(q)
    qcut = qmin + (1 - q_trim_high) * (qmax - qmin)
    allowed = q .<= qcut
    idx, z = hybrid_merge(base.indices, Z, zraw, blocked, k, kb, sep, allowed)

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
                              k                = 20.0,
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
# relative noise shape, so it estimates σ² = g · max(I, 0), clamped from below
# at the second-difference noise scale s. The gain g comes from the same second
# differences: the median of (d²/6) / I over points with I > 3s, divided by 0.455
# (the median of a χ² variable with one degree of freedom). g, s and I all scale
# with the intensity units, so the estimate does too. With no such point the
# estimate is the constant s.
function _estimate_sigma(I)
    s = second_difference_scale(I)
    n = length(I)
    ratios = Float64[]
    for i in 2:(n - 1)
        x = I[i - 1] - 2I[i] + I[i + 1]
        isfinite(x) && isfinite(I[i]) && I[i] > 3s && push!(ratios, x^2 / 6 / I[i])
    end
    g = isempty(ratios) ? 0.0 : median(ratios) / 0.455
    [max(isfinite(x) && x > 0 ? sqrt(g * x) : 0.0, s) for x in I]
end

findpeaks(q, I; kwargs...) =
    findpeaks(q, I, get(kwargs, :method, :persistence) === :hybrid ?
                    _estimate_sigma(I) : ones(length(I)); kwargs...)
