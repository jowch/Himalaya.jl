# Candidate generation and per-candidate features for one trace.
# A candidate is a grid index; features are computed from the trace alone, so
# `predict_peaks` needs nothing but (q, I, σ).

const MF_Z_FLOOR   = 3.0   # candidates are local maxima of the matched-filter z above this
const CAND_SEP     = 3     # points; closer candidates collapse to the stronger one
const Q_TRIM_HIGH  = 0.07  # same top-of-range trim as the hybrid detector
const MIN_POINTS   = Himalaya.MF_MIN_LENGTH

const WIDTH_NAMES = ["z_w$k" for k in 1:length(Himalaya.MF_WIDTHS)]
const FEATURE_NAMES = vcat(WIDTH_NAMES, ["zraw", "z_best", "best_width", "n_widths",
    "rel_height", "rel_height_wide", "persist_prom", "sharp_rel", "curv_rel", "log_bg_slope",
    "q", "pos_frac", "log_I", "sigma_rel", "z_over_max", "z_rank_frac", "n_strong",
    "near_dist", "near_z", "in_persistence", "in_hybrid", "blocked"])

# Features that echo which detector found the candidate. They line up with the
# labels when the stored auto peaks came from the same detector, so the
# ablation drops them.
const DETECTOR_FEATURES = ["in_persistence", "in_hybrid"]

"""
    candidate_indices(q, I, σ) -> (idx, Z, zraw, blocked, base, hyb)

Union of the persistence detector, the hybrid detector, and local maxima of the
matched-filter z above `MF_Z_FLOOR`. Empty `idx` for traces shorter than
`MIN_POINTS`.
"""
function candidate_indices(q, I, σ)
    n = length(I)
    n < MIN_POINTS && return (idx = Int[], Z = zeros(0, 0), zraw = Float64[],
                              blocked = falses(0), base = Int[], hyb = Int[])
    base = Himalaya.findpeaks(q, I, σ; method = :persistence, q_trim_high = Q_TRIM_HIGH)
    hyb  = Himalaya.findpeaks(q, I, σ; method = :hybrid)
    Z, zraw, blocked = Himalaya.mf_width_z(I, σ)
    z  = Himalaya.mf_best(Z)
    wb = Himalaya.mf_best_width(Z)
    qmin, qmax = extrema(q)
    allowed = q .<= qmin + (1 - Q_TRIM_HIGH) * (qmax - qmin)
    mf = Himalaya.mf_local_maxima(z, Z, wb, MF_Z_FLOOR, CAND_SEP, allowed)
    pool = unique(vcat(base.indices, hyb.indices, mf))
    zz(i) = isnan(z[i]) ? -Inf : z[i]
    sort!(pool; by = i -> (i in base.indices ? 0 : 1, -zz(i)))
    keep = Int[]
    for i in pool
        all(j -> abs(i - j) > CAND_SEP, keep) && push!(keep, i)
    end
    (idx = sort!(keep), Z = Z, zraw = zraw, blocked = blocked,
     base = base.indices, hyb = hyb.indices)
end

# Least-squares quadratic over ±h points around i; returns (slope, curvature).
# With `base`, fits `(y .- base) ./ scale` over the window only.
function local_quadratic(y, i, h, base = nothing, scale = 1.0)
    a = max(1, i - h); b = min(length(y), i + h)
    b - a < 4 && return (0.0, 0.0)
    x = collect(a:b) .- i
    v = base === nothing ? y[a:b] : (y[a:b] .- base[a:b]) ./ scale
    all(isfinite, v) || return (0.0, 0.0)
    c = [ones(length(x)) x x .^ 2] \ v
    (c[2], c[3])
end

nz(x) = isfinite(x) ? Float64(x) : 0.0

"""
    candidate_features(q, I, σ) -> (idx, X, hyb)

`hyb` is the hybrid detector's peak indices, from the same candidate pass.
`X` is `length(idx) × length(FEATURE_NAMES)`. Matched-filter responses are
scale-free (noise-normalised), and the remaining features are ratios or
logs, so the values do not depend on intensity units.
"""
function candidate_features(q, I, σ)
    c = candidate_indices(q, I, σ)
    idx = c.idx
    isempty(idx) && return idx, zeros(0, length(FEATURE_NAMES)), c.hyb
    n = length(I)
    s = let p = filter(x -> isfinite(x) && x > 0, σ); isempty(p) ? 1.0 : median(p) end
    B  = Himalaya.snip_background(I, s)
    B2 = Himalaya.snip_background(I, s, 200)
    Bs = max.(B, s)
    B2s = max.(B2, s)
    z  = Himalaya.mf_best(c.Z)
    wb = Himalaya.mf_best_width(c.Z)
    zr = [nz(x) for x in c.zraw]
    pers = Himalaya.persistence((I .- B) ./ Bs)
    pmap = Dict(zip(pers.indices, pers.prominence))
    sh = Himalaya.sharpness(I)
    lq = log.(q); lB = log.(Bs)
    zc = [nz(z[i]) for i in idx]
    zmax = max(maximum(zc), 1e-9)
    order = sortperm(-zc)
    rank = zeros(Int, length(idx)); rank[order] = 1:length(idx)
    nstrong = count(>(5.0), zc)
    X = zeros(length(idx), length(FEATURE_NAMES))
    for (r, i) in enumerate(idx)
        w = max(1, i - 2):min(n, i + 2)
        zw = [maximum((nz(c.Z[t, k]) for t in w); init = 0.0) for k in axes(c.Z, 2)]
        ndef = count(k -> !isnan(c.Z[i, k]), axes(c.Z, 2))
        rel = maximum((nz((I[t] - B[t]) / Bs[t]) for t in w))
        rel2 = maximum((nz((I[t] - B2[t]) / B2s[t]) for t in w))
        _, cv = local_quadratic(I, i, 3, B, Bs[i])
        sl, _ = local_quadratic(lB, i, 10)
        dq = lq[min(n, i + 1)] - lq[max(1, i - 1)]
        slope = dq > 0 ? nz(sl * 2 / dq) : 0.0
        others = [j for j in 1:length(idx) if j != r]
        nd = isempty(others) ? n : minimum(abs(idx[j] - i) for j in others)
        nj = isempty(others) ? 0.0 : zc[others[argmin([abs(idx[j] - i) for j in others])]]
        X[r, :] = vcat(zw, [zr[i], nz(z[i]), wb[i], ndef,
            rel, rel2, nz(get(pmap, i, 0.0)), nz(nz(sh[i]) / Bs[i]), cv,
            sign(slope) * log1p(abs(slope)),
            q[i], i / n, log(max(nz(I[i]), 1e-12)), nz(σ[i] / Bs[i]), zc[r] / zmax, rank[r] / length(idx),
            nstrong, nd, nj,
            i in c.base ? 1.0 : 0.0, i in c.hyb ? 1.0 : 0.0, any(view(c.blocked, w)) ? 1.0 : 0.0])
    end
    idx, X, c.hyb
end
