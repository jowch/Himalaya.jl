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

"""
    matched_filter_z(I, σ; snip_window = MF_SNIP_WINDOW) -> Vector{Float64}

Per-point matched-filter significance: the maximum over the Ricker widths
`MF_WIDTHS` of the SNIP-subtracted residual correlated with a unit-norm,
zero-mean kernel, divided by the kernel-weighted noise `sqrt(Σ k² σ²)`.

All statistics run over the whole trace. Masked runs (see `mf_masked_runs`)
and nonfinite points are first replaced by linear interpolation of `I`, which
is used for both the background and the residual; `z` is `NaN` on masked
points. The background is the SNIP baseline (`snip_background`) of `I` scaled by
`s`, the median positive `σ`, so the result does not depend on intensity
units. Each width is normalised by its median and `1.4826·MAD` over the points
where that kernel fits inside the trace and that are not masked; a width with
zero MAD is skipped. The maximum at a point runs over the widths that fit
there; where none fits the value is `NaN` (undefined). Nonpositive or
nonfinite `σ` are replaced by the smallest positive `σ` (or 1 if there is
none).
"""
function matched_filter_z(I, σ; snip_window::Integer = MF_SNIP_WINDOW)
    n = length(I)
    length(σ) == n || throw(DimensionMismatch("I and σ must be equal length"))
    pos = [Float64(x) for x in σ if isfinite(x) && x > 0]
    σmin = isempty(pos) ? 1.0 : minimum(pos)
    s = isempty(pos) ? 1.0 : median(pos)
    σc = [isfinite(x) && x > 0 ? Float64(x) : σmin for x in σ]
    masked, missing = mf_masked_runs(I, σc)
    Iw = mf_interpolate(I, missing)
    B = snip_background(Iw, s, snip_window)
    r = Iw .- B
    σ2 = σc .^ 2
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
        m = 0
        @inbounds for i in (h + 1):(n - h)
            masked[i] || (tmp[m += 1] = zw[i])
        end
        m > 0 || continue
        med = median!(view(tmp, 1:m))
        c = 0
        @inbounds for i in (h + 1):(n - h)
            masked[i] || (tmp[c += 1] = abs(zw[i] - med))
        end
        mad = 1.4826 * median!(view(tmp, 1:m))
        mad > 0 || continue
        @inbounds for i in (h + 1):(n - h)
            masked[i] && continue
            z = (zw[i] - med) / mad
            zbest[i] = isnan(zbest[i]) ? z : max(zbest[i], z)
        end
    end
    zbest
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
