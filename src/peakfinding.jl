using LsqFit
using Statistics: median

# findpeaks is being rebuilt — see docs/superpowers/specs/2026-04-22-peakfinding-rewrite-design.md
# Tasks 4-9 fill in helpers; Task 10 implements the public function.

"""
    ricker(t, a)

Ricker (Mexican-hat) wavelet of width `a`, evaluated at `t`. Equal to the
negative normalized second derivative of a Gaussian.

    ψ_a(t) = (2 / (√(3a) π^(1/4))) · (1 − t²/a²) · exp(−t² / 2a²)
"""
function ricker(t, a)
    norm_const = 2 / (sqrt(3a) * pi^0.25)
    norm_const * (1 - (t/a)^2) * exp(-(t^2) / (2a^2))
end

"""
    cwt(y, scales)

Continuous wavelet transform of `y` using Ricker wavelets at each width in
`scales`. Returns a `length(y) × length(scales)` matrix of coefficients.

Convolution is direct (no FFT). At array edges, the signal is reflected.
"""
function cwt(y, scales)
    n = length(y)
    coeffs = zeros(n, length(scales))
    for (j, a) in enumerate(scales)
        # truncate the kernel beyond ±5a (wavelet ≈ 0)
        m = ceil(Int, 5a)
        kernel = [ricker(t, a) for t in -m:m]
        for i in 1:n
            s = 0.0
            for (k, w) in enumerate(kernel)
                idx = i + (k - m - 1)
                # reflect at edges
                if idx < 1
                    idx = 2 - idx
                elseif idx > n
                    idx = 2n - idx
                end
                s += w * y[idx]
            end
            # normalise by scale so coefficients are comparable across scales
            coeffs[i, j] = s / a
        end
    end
    coeffs
end

"""
    local_maxima(v)

Return indices `i` where `v[i] > v[i-1]` and `v[i] > v[i+1]` and `v[i] > 0`.
Endpoints are never returned.
"""
function local_maxima(v)
    [i for i in 2:(length(v) - 1) if v[i] > 0 && v[i] > v[i-1] && v[i] > v[i+1]]
end

"""
    find_ridges(coeffs, scales; min_ridge_length = 3, max_gap = 1)

Greedy ridge linking across scales. Walks from largest scale to smallest,
extending a ridge if the next-smaller scale has a local maximum within
`±max(1, scale/4)` points of the current ridge head. Allows up to `max_gap`
missing scales before terminating a ridge.

Returns a vector of NamedTuples `(; index, scale, length, max_coeff)`:
- `index`: position of the ridge at the scale where its CWT coefficient is
  maximal (the "best-matched scale")
- `scale`: that best-matched scale value
- `length`: number of scales the ridge spans
- `max_coeff`: maximum CWT coefficient along the ridge

Ridges shorter than `min_ridge_length` scales are discarded.
"""
function find_ridges(coeffs, scales; min_ridge_length = 3, max_gap = 1)
    nscales = length(scales)
    # Per-scale local maxima
    maxima = [local_maxima(coeffs[:, j]) for j in 1:nscales]

    # Each ridge: vector of (scale_idx, position, coeff). Walk largest → smallest.
    ridges = Vector{Vector{Tuple{Int,Int,Float64}}}()
    # Track which maxima at each scale have been consumed
    consumed = [falses(length(maxima[j])) for j in 1:nscales]

    # Seed ridges at the largest scale
    for (k, p) in enumerate(maxima[nscales])
        push!(ridges, [(nscales, p, coeffs[p, nscales])])
        consumed[nscales][k] = true
    end

    # Extend each ridge downward through scales
    for j in (nscales - 1):-1:1
        tol = max(1, ceil(Int, scales[j] / 4))
        for ridge in ridges
            head_scale_idx, head_pos, _ = ridge[end]
            gap = head_scale_idx - j - 1
            if gap > max_gap
                continue   # ridge already terminated
            end
            # find an unconsumed maximum at scale j within tolerance
            best_k = 0
            best_dist = tol + 1
            for (k, p) in enumerate(maxima[j])
                if consumed[j][k]; continue; end
                d = abs(p - head_pos)
                if d <= tol && d < best_dist
                    best_dist = d
                    best_k = k
                end
            end
            if best_k > 0
                p = maxima[j][best_k]
                push!(ridge, (j, p, coeffs[p, j]))
                consumed[j][best_k] = true
            end
        end
        # Any unconsumed maxima at scale j start new ridges
        for (k, p) in enumerate(maxima[j])
            if !consumed[j][k]
                push!(ridges, [(j, p, coeffs[p, j])])
                consumed[j][k] = true
            end
        end
    end

    # Filter and summarize
    out = NamedTuple{(:index, :scale, :length, :max_coeff),
                     Tuple{Int,Float64,Int,Float64}}[]
    for ridge in ridges
        if length(ridge) < min_ridge_length; continue; end
        # Best-matched scale = where coeff is maximal
        _, idx_in_ridge = findmax(map(t -> t[3], ridge))
        scale_idx, pos, coeff = ridge[idx_in_ridge]
        push!(out, (index = pos, scale = scales[scale_idx],
                    length = length(ridge), max_coeff = coeff))
    end
    out
end

"""
    fit_peak(q, I, σ, centre_idx, width_pts)

Locally fit a Gaussian + linear baseline to the trace `(q, I)` with
per-point uncertainties `σ`, in a window of `±2.5 · width_pts` points
around `centre_idx`. Uses weighted Levenberg-Marquardt (`LsqFit.curve_fit`)
with weights `1/σᵢ²`.

Returns a NamedTuple `(; A, q0, σw, snr)` on success:
- `A`: fitted peak amplitude (above baseline)
- `q0`: fitted peak centre in q-units
- `σw`: fitted Gaussian σ in q-units
- `snr`: `A / σ_A` where `σ_A` is the standard error on `A` from the
  covariance matrix

Returns `nothing` on convergence failure or non-physical fits (A ≤ 0,
σw ≤ 0, or window too small).
"""
function fit_peak(q, I, σ, centre_idx, width_pts)
    n = length(q)
    half = max(5, ceil(Int, 2.5 * width_pts))
    lo = max(1, centre_idx - half)
    hi = min(n, centre_idx + half)
    if hi - lo < 8
        return nothing
    end

    qw = q[lo:hi]
    Iw = I[lo:hi]
    σw_arr = σ[lo:hi]

    Δq = q[2] - q[1]
    A0 = I[centre_idx] - median(Iw)
    q0_init = q[centre_idx]
    σw_init = max(width_pts * Δq / 2.355, Δq)
    b0 = (Iw[1] + Iw[end]) / 2
    m0 = (Iw[end] - Iw[1]) / (qw[end] - qw[1])

    model(x, p) = p[1] .* exp.(-((x .- p[2]).^2) ./ (2 * p[3]^2)) .+ p[4] .+ p[5] .* x
    p0 = [A0, q0_init, σw_init, b0, m0]
    weights = 1.0 ./ (σw_arr .^ 2)

    fit = try
        curve_fit(model, qw, Iw, weights, p0)
    catch
        return nothing
    end

    fit.converged || return nothing
    A_fit, q0_fit, σw_fit = fit.param[1], fit.param[2], fit.param[3]
    if A_fit <= 0 || σw_fit <= 0
        return nothing
    end

    σ_A = try
        stderror(fit)[1]
    catch
        return nothing
    end
    if !isfinite(σ_A) || σ_A <= 0
        return nothing
    end

    (A = A_fit, q0 = q0_fit, σw = σw_fit, snr = A_fit / σ_A)
end
