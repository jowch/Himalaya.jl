using LinearAlgebra: I

const DEFAULT_SCALES = [2.0, 3.0, 5.0, 8.0, 13.0]

# ---------------------------------------------------------------------------
# Public: top-level sharpness dispatcher
# ---------------------------------------------------------------------------

"""
    sharpness(y; method = :savgol, scales = DEFAULT_SCALES, m = 5) -> Vector{Float64}

Compute a per-sample sharpness measure for the 1D signal `y`. High value =
sharply curved local feature. Two interchangeable methods:

- `:savgol` (default) — magnitude of the smoothed second derivative via a
  Savitzky-Golay filter (window `2m+1`, polynomial order 4). Single-scale,
  fast, O(n).
- `:cwt` — maximum Ricker (Mexican-hat) wavelet response across `scales`,
  normalised by 1/a so coefficients are scale-comparable. Multi-scale,
  width-aware. O(n × length(scales)).
"""
function sharpness(y; method = :savgol, scales = DEFAULT_SCALES, m = 5)
    method === :savgol && return sharpness_savgol(y, m)
    method === :cwt    && return sharpness_cwt(y, scales)
    throw(ArgumentError("Unknown sharpness method: $(repr(method)) (use :savgol or :cwt)"))
end

# ---------------------------------------------------------------------------
# Savitzky-Golay second derivative (sharpness_savgol)
# ---------------------------------------------------------------------------

"""
    sharpness_savgol(y, m) -> Vector{Float64}

Return `-d²y/dx²` computed via Savitzky-Golay smoothing (window `2m+1`,
polynomial order 4). A sharp local maximum has strongly negative d²y, so
we flip sign — large positive output = sharply peaked.
"""
sharpness_savgol(y, m) = -savitzky_golay(m, 4, y; order = 2)

# `savitzky_golay(m, n, y; order)` — preserved verbatim from the v0.4.5
# src/peakfinding.jl (commit fddd611). General-purpose SG: any window,
# polynomial order, derivative order.
function savitzky_golay(m, n, y; order = 0)
    num_y = length(y)
    z = -m:m
    J = zeros(2m + 1, n + 1)

    for i = 0:n
        @inbounds J[:, i + 1] .= z .^ i
    end

    # The convolution term matrix
    C = J' \ I(n .+ 1)[:, order .+ 1]   # = pinv(J) picking out the requested order(s)
    Y = zeros(num_y, length(order))

    for i in 1:num_y
        if i <= m
            window_indices = abs.(z .+ i) .+ 1
        elseif i > num_y - m
            window_indices = -abs.(z .+ i .- num_y) .+ num_y
        else
            window_indices = z .+ i
        end

        for j in eachindex(order)
            @inbounds Y[i, j] = C[:, j]' * y[window_indices]
        end
    end

    # SG extracts polynomial coefficients c_k; the k-th derivative is k! * c_k.
    fac = factorial.(order)
    for j in eachindex(order)
        Y[:, j] .*= fac[j]
    end

    length(order) == 1 ? Y[:, 1] : Y
end

# ---------------------------------------------------------------------------
# CWT-based sharpness (sharpness_cwt)
# ---------------------------------------------------------------------------

"""
    sharpness_cwt(y, scales) -> Vector{Float64}

For each sample of `y`, return the maximum Ricker-wavelet response across
`scales`, normalised by 1/a so coefficients are scale-comparable. Peak-like
features produce their maximum response at the scale matching their width.
"""
function sharpness_cwt(y, scales)
    n = length(y)
    out = fill(-Inf, n)
    for a in scales
        m = ceil(Int, 5a)
        kernel = [ricker(t, a) for t in -m:m]
        response = convolve_reflect(y, kernel) ./ a
        out .= max.(out, response)
    end
    out
end

"""
    ricker(t, a)

Ricker (Mexican-hat) wavelet of width `a`, evaluated at `t`. Equal to the
negative normalized second derivative of a Gaussian:

    ψ_a(t) = (2 / (√(3a) π^(1/4))) · (1 − t²/a²) · exp(−t² / (2a²))
"""
ricker(t, a) = (2 / (sqrt(3a) * π^0.25)) * (1 - (t/a)^2) * exp(-(t^2) / (2a^2))

# ---------------------------------------------------------------------------
# Shared: direct convolution with edge reflection
# ---------------------------------------------------------------------------

"""
    convolve_reflect(y, kernel) -> Vector{Float64}

Direct 1D convolution of `y` with a symmetric `kernel` of length `2m+1`,
reflecting at edges (no wrap-around, no zero-padding artefacts).
"""
function convolve_reflect(y, kernel)
    n = length(y)
    m = (length(kernel) - 1) ÷ 2
    out = zeros(n)
    for i in 1:n, k in 1:length(kernel)
        idx = i + (k - m - 1)
        idx = idx < 1 ? 2 - idx : idx > n ? 2n - idx : idx
        out[i] += kernel[k] * y[idx]
    end
    out
end
