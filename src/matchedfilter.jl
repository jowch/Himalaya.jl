using Statistics: median

# Matched-filter significance used by the opt-in `:hybrid` detector in
# `findpeaks` (see docs/peak-finding.md, "Hybrid detector").

const MF_WIDTHS      = (2.0, 3.0, 4.5, 7.0, 11.0, 16.0)
const MF_SNIP_WINDOW = 60

"""
    second_difference_scale(I) -> Float64

Noise scale of `I` from second differences, for use without a provided `σ`.
For independent noise of standard deviation `σ`, `I[i-1] - 2I[i] + I[i+1]` has
standard deviation `σ√6`, and second differences remove a locally linear
background. The estimate is `1.4826·MAD/√6`; if that is 0 (quantised counts)
it is the standard deviation of the differences over `√6`, and if that is also
0 it is `mean(|I|)`, then 1. It is proportional to `I`.
"""
function second_difference_scale(I)
    d = Float64[]
    # Differences below rounding error are zero: scaled integer counts must give
    # the same exact zeros at every scale.
    tol = 1e-9 * maximum((abs(x) for x in I if isfinite(x)); init = 0.0)
    for i in 2:(length(I) - 1)
        x = I[i - 1] - 2I[i] + I[i + 1]
        isfinite(x) && push!(d, abs(x) <= tol ? 0.0 : x)
    end
    if !isempty(d)
        m = median(d)
        mad = 1.4826 * median(abs.(d .- m))
        mad > 0 && return mad / sqrt(6)
        mu = sum(d) / length(d)
        sd = sqrt(sum(abs2, d .- mu) / length(d))
        sd > 0 && return sd / sqrt(6)
    end
    f = [abs(x) for x in I if isfinite(x)]
    isempty(f) && return 1.0
    mu = sum(f) / length(f)
    mu > 0 ? mu : 1.0
end

"""
    snip_background(I, s, m) -> Vector{Float64}

SNIP baseline of `I` with half-window `m` points, using the LLS operator on
data scaled by the noise scale `s`: `y = max(I/s, 0)`,
`v = log(log(sqrt(y + 1) + 1) + 1)`, clip `v[i] = min(v[i], (v[i-w] + v[i+w]) / 2)`
for `w = m … 1`, then invert and multiply by `s`. The operator compresses the
dynamic range without a hole at zero, so zero and negative points need no
special handling. Nonfinite points are replaced by linear interpolation in
`v` between the nearest finite neighbours (nearest value at the ends); with
no finite point the baseline is zero.
"""
function snip_background(I, s::Real, m::Integer = MF_SNIP_WINDOW)
    n = length(I)
    v = Vector{Float64}(undef, n)
    u = Vector{Float64}(undef, n)
    prev = 0
    @inbounds for i in 1:n
        if isfinite(I[i])
            y = max(I[i] / s, 0.0)
            v[i] = log(log(sqrt(y + 1) + 1) + 1)
            if prev != i - 1
                if prev == 0
                    for j in 1:(i - 1); v[j] = v[i]; end
                else
                    for j in (prev + 1):(i - 1)
                        v[j] = v[prev] + (v[i] - v[prev]) * (j - prev) / (i - prev)
                    end
                end
            end
            prev = i
        end
    end
    prev == 0 && return zeros(n)
    @inbounds for j in (prev + 1):n; v[j] = v[prev]; end
    for w in m:-1:1
        copyto!(u, v)
        @inbounds for i in (w + 1):(n - w)
            u[i] = min(v[i], (v[i - w] + v[i + w]) / 2)
        end
        v, u = u, v
    end
    @inbounds for i in 1:n
        v[i] = s * ((exp(exp(v[i]) - 1) - 1)^2 - 1)
    end
    v
end

# Unit-norm, zero-mean Ricker kernel of half-width ceil(4s).
function ricker_kernel(s::Real)
    h = ceil(Int, 4s)
    k = Vector{Float64}(undef, 2h + 1)
    @inbounds for j in -h:h
        x2 = (j / s)^2
        k[j + h + 1] = (1 - x2) * exp(-x2 / 2)
    end
    k .-= sum(k) / length(k)
    k ./= sqrt(sum(abs2, k))
    k
end

const MF_KERNELS = map(ricker_kernel, MF_WIDTHS)

const MF_MIN_LENGTH = 150
const MF_MASK_NEIGHBOURS = 7

# Maximal runs of points that are exactly 0.0 or nonfinite and whose
# neighbourhood is clearly above the noise: the median of the nearest
# `MF_MASK_NEIGHBOURS` finite, non-run points on every available side exceeds
# 3 times their rms σ. Such a run is a masked region (beamstop, detector gap),
# not data. Zero runs in a near-zero region (a true zero baseline, a negative
# tail) fail the test and stay data. Returns (masked, missing): `masked` marks
# the runs that pass; `missing` additionally marks every nonfinite point.
function mf_masked_runs(I, σ)
    n = length(I)
    inrun = [!isfinite(I[i]) || I[i] == 0.0 for i in 1:n]
    masked = falses(n)
    i = 1
    while i <= n
        if !inrun[i]
            i += 1
            continue
        end
        j = i
        while j < n && inrun[j + 1]; j += 1; end
        sides = Vector{Vector{Int}}()
        left = Int[]
        t = i - 1
        while t >= 1 && length(left) < MF_MASK_NEIGHBOURS
            inrun[t] || push!(left, t)
            t -= 1
        end
        right = Int[]
        t = j + 1
        while t <= n && length(right) < MF_MASK_NEIGHBOURS
            inrun[t] || push!(right, t)
            t += 1
        end
        isempty(left) || push!(sides, left)
        isempty(right) || push!(sides, right)
        ok = !isempty(sides) && all(sides) do side
            median(I[side]) > 3 * sqrt(sum(σ[t]^2 for t in side) / length(side))
        end
        ok && (masked[i:j] .= true)
        i = j + 1
    end
    masked, [!isfinite(I[i]) for i in 1:n] .| masked
end

# Replace the points marked in `fill` by linear interpolation between the
# nearest unmarked points (nearest value at the ends).
function mf_interpolate(I, fill)
    n = length(I)
    out = Float64.(I)
    i = 1
    while i <= n
        if !fill[i]
            i += 1
            continue
        end
        j = i
        while j < n && fill[j + 1]; j += 1; end
        lo, hi = i - 1, j + 1
        for t in i:j
            out[t] = if lo >= 1 && hi <= n
                out[lo] + (out[hi] - out[lo]) * (t - lo) / (hi - lo)
            elseif lo >= 1
                out[lo]
            elseif hi <= n
                out[hi]
            else
                0.0
            end
        end
        i = j + 1
    end
    out
end

const MF_CLIP = 3.5
const MF_MIN_NOISE_POINTS = 30

# Robust centre and scale of `x` over the points marked in `use`, with points
# far from the centre (signal) removed iteratively: points with
# |x - centre| > MF_CLIP·scale, together with `dilate` points on each side, are
# dropped and the estimate repeated until the set stops changing. Returns
# `nothing` when fewer than MF_MIN_NOISE_POINTS points (or a fifth of the
# starting set) remain or the scale is zero. After 20 rounds that do not
# settle, the estimate on the last set is returned.
function mf_noise_scale(x, use, dilate)
    n = length(x)
    minpts = max(MF_MIN_NOISE_POINTS, count(use) ÷ 5)
    cur = copy(use)
    buf = Float64[]
    for round in 1:21
        empty!(buf)
        for i in 1:n
            cur[i] && push!(buf, x[i])
        end
        length(buf) >= minpts || return nothing
        med = median!(buf)
        buf .= abs.(buf .- med)
        mad = 1.4826 * median!(buf)
        mad > 0 || return nothing
        round == 21 && return (med, mad)
        next = copy(use)
        for i in 1:n
            if use[i] && abs(x[i] - med) > MF_CLIP * mad
                next[max(1, i - dilate):min(n, i + dilate)] .= false
            end
        end
        next .&= use
        next == cur && return (med, mad)
        cur = next
    end
end

# Matched-filter responses, NaN where undefined. Returns
# - `Z`: n × length(MF_WIDTHS), each width on its own clipped scale;
# - `zraw`: per point, the max over widths of the response on the
#   second-difference scale (`mf_white_scale`);
# - `blocked`: points where some width's window covers a missing point, or
#   would cross a span edge that borders a missing run (a masked beamstop).
# See `matched_filter_z`.
function mf_width_z(I, σ; snip_window::Integer = MF_SNIP_WINDOW)
    n = length(I)
    length(σ) == n || throw(DimensionMismatch("I and σ must be equal length"))
    Z = fill(NaN, n, length(MF_KERNELS))
    zraw = fill(NaN, n)
    blocked = falses(n)
    pos = [Float64(x) for x in σ if isfinite(x) && x > 0]
    σmin = isempty(pos) ? 1.0 : minimum(pos)
    s = isempty(pos) ? 1.0 : median(pos)
    σc = [isfinite(x) && x > 0 ? Float64(x) : σmin for x in σ]
    _, missing = mf_masked_runs(I, σc)
    a = findfirst(!, missing)
    a === nothing && return Z, zraw, blocked
    b = findlast(!, missing)
    seg = a:b
    Iw = mf_interpolate(view(I, seg), view(missing, seg))
    r = Iw .- snip_background(Iw, s, snip_window)
    σ2 = view(σc, seg) .^ 2
    m = length(seg)
    # cmiss[t] = number of missing points among seg[1:t-1]
    cmiss = zeros(Int, m + 1)
    for t in 1:m
        cmiss[t + 1] = cmiss[t] + missing[seg[t]]
    end
    cw = mf_white_scale(Iw, view(σc, seg), view(missing, seg))
    zw = fill(NaN, m)
    use = falses(m)
    for (w, k) in enumerate(MF_KERNELS)
        h = (length(k) - 1) ÷ 2
        fill!(zw, NaN)
        fill!(use, false)
        # A span edge next to a missing run blocks like an interior run does.
        a > 1 && (blocked[seg[1:min(h, m)]] .= true)
        b < n && (blocked[seg[max(1, m - h + 1):m]] .= true)
        @inbounds for i in (h + 1):(m - h)
            if cmiss[i + h + 1] - cmiss[i - h] > 0
                blocked[seg[i]] = true
                continue
            end
            num = 0.0
            var = 0.0
            for j in -h:h
                kj   = k[j + h + 1]
                num += kj * r[i + j]
                var += kj * kj * σ2[i + j]
            end
            zw[i] = num / sqrt(var)
            use[i] = true
        end
        @inbounds for i in 1:m
            if use[i]
                v = zw[i] / cw
                zraw[seg[i]] = isnan(zraw[seg[i]]) ? v : max(zraw[seg[i]], v)
            end
        end
        est = mf_noise_scale(zw, use, ceil(Int, MF_WIDTHS[w]))
        est === nothing && continue
        med, mad = est
        @inbounds for i in 1:m
            use[i] && (Z[seg[i], w] = (zw[i] - med) / mad)
        end
    end
    Z, zraw, blocked
end

# Noise scale of I relative to σ from second differences, which a smooth
# background barely affects: the clipped scale (`mf_noise_scale`) of
# (I[i-1] - 2I[i] + I[i+1]) / sqrt(σ[i-1]² + 4σ[i]² + σ[i+1]²) over triples with
# no missing point, so the curvature of sharp peaks is removed as signal. When
# that scale is undefined (mostly zero differences, as in sparse counts) it is
# their standard deviation, as in `second_difference_scale`; 1 if that is also
# zero or there are no triples.
function mf_white_scale(I, σ, missing)
    m = length(I)
    d = zeros(m)
    use = falses(m)
    # Differences below rounding error are zero, as in `second_difference_scale`.
    tol = 1e-9 * maximum((abs(x) for x in I if isfinite(x)); init = 0.0)
    for i in 2:(m - 1)
        (missing[i - 1] || missing[i] || missing[i + 1]) && continue
        x = I[i - 1] - 2I[i] + I[i + 1]
        d[i] = abs(x) <= tol ? 0.0 : x / sqrt(σ[i - 1]^2 + 4σ[i]^2 + σ[i + 1]^2)
        use[i] = true
    end
    est = mf_noise_scale(d, use, 2)
    est === nothing || return est[2]
    v = d[use]
    isempty(v) && return 1.0
    mu = sum(v) / length(v)
    sd = sqrt(sum(abs2, v .- mu) / length(v))
    sd > 0 ? sd : 1.0
end

"""
    matched_filter_z(I, σ; snip_window = MF_SNIP_WINDOW) -> Vector{Float64}

Per-point matched-filter significance: the maximum over the Ricker widths
`MF_WIDTHS` of the SNIP-subtracted residual correlated with a unit-norm,
zero-mean kernel, divided by the kernel-weighted noise `sqrt(Σ k² σ²)`.

Masked runs (see `mf_masked_runs`) and nonfinite points are missing. The trace
is cut to the span between the first and last present points; inside it,
missing points are filled by linear interpolation for the background only. The
background is the SNIP baseline (`snip_background`) of that span, scaled by `s`,
the median positive `σ`, so the result does not depend on intensity units.

A width is defined at a point when its kernel fits inside the span and covers
no missing point. Each width is normalised by a noise centre and scale
estimated from its own values with signal removed (`mf_noise_scale`), so dense
or strong peaks do not inflate the scale. A width with too few signal-free
points, or zero spread among them, is skipped. The maximum at a point runs over the widths defined there;
where none is, the value is `NaN`. Nonpositive or nonfinite `σ` are replaced by
the smallest positive `σ` (or 1 if there is none).
"""
function matched_filter_z(I, σ; snip_window::Integer = MF_SNIP_WINDOW)
    mf_best(first(mf_width_z(I, σ; snip_window)))
end

function mf_best(Z)
    n = size(Z, 1)
    z = fill(NaN, n)
    @inbounds for w in axes(Z, 2), i in 1:n
        v = Z[i, w]
        isnan(v) || (z[i] = isnan(z[i]) ? v : max(z[i], v))
    end
    z
end

# Width index with the largest z at each point (0 where none is defined).
function mf_best_width(Z)
    n = size(Z, 1)
    wb = zeros(Int, n)
    @inbounds for i in 1:n
        best = -Inf
        for w in axes(Z, 2)
            v = Z[i, w]
            if !isnan(v) && v > best
                best = v
                wb[i] = w
            end
        end
    end
    wb
end

# An added peak at `i` must not sit in a valley at a finer scale: at every
# width narrower than its best width `w`, z within ±1 point is above
# -MF_VALLEY_Z. Between two close peaks a wide kernel sees one blob and peaks
# at the midpoint, while narrower kernels sit in the valley and respond
# strongly negatively; a broad peak gives a weak or positive narrow response.
const MF_VALLEY_Z = 3.0
function mf_narrower_agrees(Z, i, w)
    n = size(Z, 1)
    for v in 1:(w - 1)
        top = -Inf
        for t in max(1, i - 1):min(n, i + 1)
            isnan(Z[t, v]) || (top = max(top, Z[t, v]))
        end
        top == -Inf && continue
        top > -MF_VALLEY_Z || return false
    end
    true
end

# Local maxima of z with z ≥ k at positions where `allowed[i]`, greedily
# thinned so chosen peaks are more than `sep` points apart (strongest first).
# A point must also be a local maximum at its own best width `wb[i]`, with both
# neighbours defined there: where a wide width becomes defined (at the trace
# ends or next to a masked run) the maximum over widths steps up, which is not
# a peak.
function mf_local_maxima(z, Z, wb, k, sep, allowed)
    n = length(z)
    cand = Int[]
    @inbounds for i in 2:(n - 1)
        (z[i] >= z[i - 1] && z[i] > z[i + 1] && z[i] >= k && allowed[i]) || continue
        w = wb[i]
        a, b, c = Z[i - 1, w], Z[i, w], Z[i + 1, w]
        (isnan(a) || isnan(c) || b < a || b <= c) && continue
        push!(cand, i)
    end
    sort!(cand; by = i -> -z[i])
    keep = Int[]
    for i in cand
        all(j -> abs(i - j) > sep, keep) && push!(keep, i)
    end
    sort!(keep)
end

# Combine persistence-detector peaks with matched-filter peaks. A persistence
# peak is dropped only if `zraw` is defined at it, stays below `kb` within ±3
# points, and no kernel there is `blocked` by missing data: too close to the
# trace ends no kernel is centred on the peak, and next to missing data (a
# masked run, including a beamstop at either end, or a single nonfinite point)
# the wide kernels cannot see a broad peak, so it is kept. `zraw` is judged against the
# second-difference noise scale, not the per-width clipped scale, so a trace
# that is mostly signal cannot lower it. An added peak is a local maximum of the
# clipped z with z ≥ k, more than `sep` points from every kept peak, that passes
# `mf_narrower_agrees`.
function hybrid_merge(base_idx, Z, zraw, blocked, k, kb, sep, allowed)
    z = mf_best(Z)
    wb = mf_best_width(Z)
    n = length(z)
    kept = Int[]
    for j in base_idx
        r = max(1, j - 3):min(n, j + 3)
        rmax = maximum((zraw[t] for t in r if !isnan(zraw[t])); init = -Inf)
        (isnan(zraw[j]) || rmax >= kb || any(view(blocked, r))) && push!(kept, j)
    end
    added = [i for i in mf_local_maxima(z, Z, wb, k, sep, allowed)
             if !any(j -> abs(i - j) <= sep, kept) && mf_narrower_agrees(Z, i, wb[i])]
    sort!(vcat(kept, added)), z
end
