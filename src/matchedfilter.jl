using Statistics: median

# Matched-filter significance used by the opt-in `:hybrid` detector in
# `findpeaks` (see docs/peak-finding.md, "Hybrid detector").

const MF_WIDTHS      = (2.0, 3.0, 4.5, 7.0, 11.0, 16.0)
const MF_SNIP_WINDOW = 60

"""
    noise_scale(I, σ) -> Float64

Robust point-to-point noise scale of `I` from the data itself:
`1.4826·median(|d − median(d)|)/√2` with `d = diff(I)` over finite pairs.
Falls back to the median positive `σ`, then to 1. It scales with `I`, so
everything built on it is independent of intensity units.
"""
function noise_scale(I, σ)
    d = Float64[]
    for i in 2:length(I)
        x = I[i] - I[i - 1]
        isfinite(x) && push!(d, x)
    end
    if !isempty(d)
        m = median!(d)
        for i in eachindex(d)
            d[i] = abs(d[i] - m)
        end
        s = 1.4826 * median!(d) / sqrt(2)
        s > 0 && isfinite(s) && return s
    end
    pos = [x for x in σ if isfinite(x) && x > 0]
    isempty(pos) ? 1.0 : median!(pos)
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

"""
    matched_filter_core(I, σ, s, snip_window) -> Vector{Float64}

Per-point matched-filter significance: the maximum over the Ricker widths
`MF_WIDTHS` of the SNIP-subtracted residual correlated with a unit-norm,
zero-mean kernel, divided by the kernel-weighted noise `sqrt(Σ k² σ²)`. The
SNIP background uses the LLS operator on data scaled by `noise_scale`, so the
result does not depend on intensity units. Nonfinite intensities have zero
residual. Each width is normalised by its median and `1.4826·MAD` over the
points where that kernel fits inside the trace (so edge points do not bias
the scale), and a width with zero MAD is skipped. The maximum at a point runs
over the widths that fit there; where none fits the value is `NaN`
(undefined). Nonpositive or nonfinite `σ` are replaced by the smallest
positive `σ` (or 1 if there is none).
"""
function matched_filter_core(I, σ, s::Real, snip_window::Integer)
    n = length(I)
    length(σ) == n || throw(DimensionMismatch("I and σ must be equal length"))
    σmin = Inf
    for x in σ
        isfinite(x) && x > 0 && x < σmin && (σmin = x)
    end
    isfinite(σmin) || (σmin = 1.0)
    σ2 = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        x = isfinite(σ[i]) && σ[i] > 0 ? σ[i] : σmin
        σ2[i] = x * x
    end
    B = snip_background(I, s, snip_window)
    r = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        r[i] = isfinite(I[i]) ? I[i] - B[i] : 0.0
    end
    zbest = fill(NaN, n)
    zw    = Vector{Float64}(undef, n)
    tmp   = Vector{Float64}(undef, n)
    for k in MF_KERNELS
        h = (length(k) - 1) ÷ 2
        n > 2h || continue
        @inbounds for i in (h + 1):(n - h)
            num = 0.0
            var = 0.0
            for j in -h:h
                kj   = k[j + h + 1]
                num += kj * r[i + j]
                var += kj * kj * σ2[i + j]
            end
            zw[i] = num / sqrt(var)
        end
        m   = n - 2h
        med = median!(copyto!(view(tmp, 1:m), view(zw, (h + 1):(n - h))))
        @inbounds for i in 1:m
            tmp[i] = abs(zw[h + i] - med)
        end
        mad = 1.4826 * median!(view(tmp, 1:m))
        mad > 0 || continue
        @inbounds for i in (h + 1):(n - h)
            z = (zw[i] - med) / mad
            zbest[i] = isnan(zbest[i]) ? z : max(zbest[i], z)
        end
    end
    zbest
end

const MF_MIN_RANGE = 150
const MF_RANGE_WINDOW = 13

# Replace nonpositive points that sit in a clearly positive neighbourhood (median
# of the ±5 neighbours above 3 noise units) by linear interpolation between the
# nearest points not so replaced. Nonpositive points in a near-zero region (a
# negative tail, a masked block) are left alone.
function fill_isolated_nonpositive(I, s)
    n = length(I)
    out = Float64.(I)
    bad = falses(n)
    nb = Float64[]
    for i in 1:n
        (isfinite(I[i]) && I[i] <= 0) || continue
        empty!(nb)
        for j in max(1, i - 5):min(n, i + 5)
            j != i && isfinite(I[j]) && push!(nb, I[j])
        end
        !isempty(nb) && median!(nb) > 3s && (bad[i] = true)
    end
    i = 1
    while i <= n
        if !bad[i]
            i += 1
            continue
        end
        j = i
        while j <= n && bad[j]; j += 1; end
        lo, hi = i - 1, j
        for t in i:(j - 1)
            if lo >= 1 && hi <= n
                out[t] = out[lo] + (out[hi] - out[lo]) * (t - lo) / (hi - lo)
            else
                out[t] = lo >= 1 ? out[lo] : out[hi]
            end
        end
        i = j
    end
    out
end

# First and last index whose `MF_RANGE_WINDOW`-point mean exceeds 3 standard
# errors, 3·rms(σ)/√c over the window's c finite points; (1, 0) if none.
# Outside this range the trace carries no signal. The test uses the local σ,
# so a steep low-q region does not set the threshold for the high-q end.
function mf_valid_range(I, σ)
    n = length(I)
    h = MF_RANGE_WINDOW ÷ 2
    lo, hi = 0, 0
    for i in 1:n
        t = 0.0; v = 0.0; c = 0
        for j in max(1, i - h):min(n, i + h)
            if isfinite(I[j])
                t += I[j]; v += σ[j]^2; c += 1
            end
        end
        if c > 0 && t > 3 * sqrt(v)
            lo == 0 && (lo = i)
            hi = i
        end
    end
    lo == 0 ? (1, 0) : (lo, hi)
end

"""
    matched_filter_z(I, σ; snip_window = MF_SNIP_WINDOW) -> Vector{Float64}

Matched-filter significance `z` (see `matched_filter_core`) computed only on
the valid range of the trace: from the first to the last point whose
13-point mean exceeds 3 standard errors (3·rms(σ)/√13, local σ), so a masked or empty
block at either end does not enter the background or the normalisation.
Nonpositive points inside a clearly positive region are first interpolated
over. `z` is `NaN` outside the valid range, and everywhere if the range is
shorter than `MF_MIN_RANGE` (150 points; the widest kernel needs 129 points to
fit, so a shorter range leaves too few points for the per-width scale).
"""
function matched_filter_z(I, σ; snip_window::Integer = MF_SNIP_WINDOW)
    n = length(I)
    length(σ) == n || throw(DimensionMismatch("I and σ must be equal length"))
    s  = noise_scale(I, σ)
    Iw = fill_isolated_nonpositive(I, s)
    σmin = Inf
    for x in σ
        isfinite(x) && x > 0 && x < σmin && (σmin = x)
    end
    isfinite(σmin) || (σmin = 1.0)
    σc = [isfinite(x) && x > 0 ? Float64(x) : σmin for x in σ]
    a, b = mf_valid_range(Iw, σc)
    z = fill(NaN, n)
    b - a + 1 >= MF_MIN_RANGE || return z
    z[a:b] = matched_filter_core(view(Iw, a:b), view(σ, a:b), s, snip_window)
    z
end

# Local maxima of z with z ≥ k at positions where `allowed[i]`, greedily
# thinned so chosen peaks are more than `sep` points apart (strongest first).
function mf_local_maxima(z, k, sep, allowed)
    n = length(z)
    cand = Int[]
    @inbounds for i in 2:(n - 1)
        z[i] >= z[i - 1] && z[i] > z[i + 1] && z[i] >= k && allowed[i] && push!(cand, i)
    end
    sort!(cand; by = i -> -z[i])
    keep = Int[]
    for i in cand
        all(j -> abs(i - j) > sep, keep) && push!(keep, i)
    end
    sort!(keep)
end

# Combine persistence-detector peaks with matched-filter peaks. A persistence
# peak is dropped only if z is defined at it and stays below `kb` within ±3
# points. An added peak must be more than `sep` points from every kept peak.
function hybrid_merge(base_idx, z, k, kb, sep, allowed)
    n = length(z)
    kept = Int[]
    for j in base_idx
        if isnan(z[j])
            push!(kept, j)
            continue
        end
        zmax = -Inf
        for t in max(1, j - 3):min(n, j + 3)
            isnan(z[t]) || (zmax = max(zmax, z[t]))
        end
        zmax >= kb && push!(kept, j)
    end
    added = [i for i in mf_local_maxima(z, k, sep, allowed) if !any(j -> abs(i - j) <= sep, kept)]
    sort!(vcat(kept, added))
end
