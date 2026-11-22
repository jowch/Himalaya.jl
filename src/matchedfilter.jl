using Statistics: median

# Matched-filter significance used by the opt-in `:hybrid` detector in
# `findpeaks` (see docs/peak-finding.md, "Hybrid detector").

const MF_WIDTHS      = (2.0, 3.0, 4.5, 7.0, 11.0, 16.0)
const MF_SNIP_WINDOW = 60

"""
    snip_background(I, m) -> Vector{Float64}

SNIP baseline of `I`, computed on `log(I)` with half-window `m` points: for
`w = m … 1`, `v[i] = min(v[i], (v[i-w] + v[i+w]) / 2)`. Points of `I` that are
nonpositive or nonfinite are replaced in log space by linear interpolation
between the nearest valid neighbours, so the result does not depend on
intensity scale. A leading or trailing block of such points carries no
information and is set to `Inf` in log space: it never clips its neighbours
(extending the nearest valid value would put a flat step under a steep
background and create a false peak at the block edge). The baseline is `Inf`
there. With no valid point the baseline is all zeros.
"""
function snip_background(I, m::Integer = MF_SNIP_WINDOW)
    n = length(I)
    v = Vector{Float64}(undef, n)
    u = Vector{Float64}(undef, n)
    valid(x) = isfinite(x) && x > 0
    prev = 0
    @inbounds for i in 1:n
        if valid(I[i])
            v[i] = log(I[i])
            if prev != i - 1
                if prev == 0
                    for j in 1:(i - 1); v[j] = Inf; end
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
    @inbounds for j in (prev + 1):n; v[j] = Inf; end
    for w in m:-1:1
        copyto!(u, v)
        @inbounds for i in (w + 1):(n - w)
            u[i] = min(v[i], (v[i - w] + v[i + w]) / 2)
        end
        v, u = u, v
    end
    @inbounds for i in 1:n
        v[i] = exp(v[i])
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
    matched_filter_z(I, σ; snip_window = MF_SNIP_WINDOW) -> Vector{Float64}

Per-point matched-filter significance: the maximum over the Ricker widths
`MF_WIDTHS` of the SNIP-subtracted residual correlated with a unit-norm,
zero-mean kernel, divided by the kernel-weighted noise `sqrt(Σ k² σ²)`. The
residual uses the original `I` at valid points and is 0 at nonpositive or nonfinite ones (a masked point is missing data, and a hole in the residual would ring through the kernels into false peaks). Each width is normalised by its median and
`1.4826·MAD`, both taken only over points where that kernel fits inside the
trace, so widths are comparable and edge points do not bias the scale. The
maximum at a point runs over the widths that fit there; where none fits the
value is `NaN` (undefined). A width with zero MAD, or one that does not fit
anywhere, is skipped. Nonpositive or nonfinite `σ` are replaced by the
smallest positive `σ` (or 1 if there is none).
"""
function matched_filter_z(I, σ; snip_window::Integer = MF_SNIP_WINDOW)
    n = length(I)
    length(σ) == n || throw(DimensionMismatch("I and σ must be equal length"))
    B = snip_background(I, snip_window)
    σmin = Inf
    for s in σ
        isfinite(s) && s > 0 && s < σmin && (σmin = s)
    end
    isfinite(σmin) || (σmin = 1.0)
    r  = Vector{Float64}(undef, n)
    σ2 = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        r[i]  = isfinite(I[i]) && I[i] > 0 ? I[i] - B[i] : 0.0
        s     = isfinite(σ[i]) && σ[i] > 0 ? σ[i] : σmin
        σ2[i] = s * s
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
