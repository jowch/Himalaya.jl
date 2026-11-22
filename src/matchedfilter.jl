using Statistics: median

# Matched-filter significance used by the opt-in `:hybrid` detector in
# `findpeaks` (see docs/peak-finding.md, "Hybrid detector").

const MF_WIDTHS      = (2.0, 3.0, 4.5, 7.0, 11.0, 16.0)
const MF_SNIP_WINDOW = 60

"""
    snip_background(I, m) -> Vector{Float64}

SNIP baseline of `I`, computed on `log(max(I, tiny))` with half-window `m`
points: for `w = m … 1`, `v[i] = min(v[i], (v[i-w] + v[i+w]) / 2)`.
"""
function snip_background(I, m::Integer = MF_SNIP_WINDOW)
    n = length(I)
    v = Vector{Float64}(undef, n)
    u = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        v[i] = log(max(I[i], 1e-6))
    end
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

Per-point matched-filter significance: the maximum over the Ricker widths `MF_WIDTHS` of
the SNIP-subtracted residual correlated with a unit-norm, zero-mean kernel,
divided by the kernel-weighted noise `sqrt(Σ k² σ²)`. Each width is
normalised by its trace median and `1.4826·MAD` before the maximum is taken,
so widths are comparable. Points whose kernel runs off the trace score 0
before normalisation. Nonpositive `σ` are replaced by the smallest positive
`σ` (or 1 if there is none).
"""
function matched_filter_z(I, σ; snip_window::Integer = MF_SNIP_WINDOW)
    n = length(I)
    length(σ) == n || throw(DimensionMismatch("I and σ must be equal length"))
    B = snip_background(I, snip_window)
    σmin = Inf
    for s in σ
        s > 0 && s < σmin && (σmin = s)
    end
    isfinite(σmin) || (σmin = 1.0)
    r  = Vector{Float64}(undef, n)
    σ2 = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        r[i]  = I[i] - B[i]
        s     = σ[i] > 0 ? σ[i] : σmin
        σ2[i] = s * s
    end
    zbest = fill(-Inf, n)
    zw    = Vector{Float64}(undef, n)
    tmp   = Vector{Float64}(undef, n)
    for k in MF_KERNELS
        h = (length(k) - 1) ÷ 2
        fill!(zw, 0.0)
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
        copyto!(tmp, zw)
        med = median!(tmp)
        @inbounds for i in 1:n
            tmp[i] = abs(zw[i] - med)
        end
        mad = 1.4826 * median!(tmp)
        if mad > 0
            @inbounds for i in 1:n
                zbest[i] = max(zbest[i], (zw[i] - med) / mad)
            end
        else
            @inbounds for i in 1:n
                zbest[i] = max(zbest[i], 0.0)
            end
        end
    end
    zbest
end

# Local maxima of z with z ≥ k, in the first `n - ntrim` points, greedily
# thinned so chosen peaks are more than `sep` points apart (strongest first).
function mf_local_maxima(z, k, sep, ntrim)
    n = length(z)
    cand = Int[]
    @inbounds for i in 2:(n - 1)
        z[i] >= z[i - 1] && z[i] > z[i + 1] && z[i] >= k && i <= n - ntrim && push!(cand, i)
    end
    sort!(cand; by = i -> -z[i])
    keep = Int[]
    for i in cand
        all(j -> abs(i - j) > sep, keep) && push!(keep, i)
    end
    sort!(keep)
end

# Combine persistence-detector peaks with matched-filter peaks. Returns the
# merged sorted index vector and a mask telling which entries came from the
# persistence detector.
function hybrid_merge(base_idx, z, k, kb, sep, ntrim, guard)
    n = length(z)
    kept = Int[]
    for j in base_idx
        a = max(1, j - 3); b = min(n, j + 3)
        maximum(view(z, a:b)) >= kb && push!(kept, j)
    end
    added = [i for i in mf_local_maxima(z, k, sep, ntrim) if !any(j -> abs(i - j) <= guard, kept)]
    sort!(vcat(kept, added)), kept
end
