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

# `savitzky_golay(m, n, y; order)` — general-purpose SG: any window, polynomial
# order, derivative order. Originally from v0.4.5 src/peakfinding.jl (commit
# fddd611); the per-sample convolution was rewritten allocation-free for issue
# #128.
#
# Edges use the "interp" boundary scheme (as in scipy.signal.savgol_filter):
# interior samples convolve a *centred* window, but the first/last `m` samples
# reuse the nearest full boundary window and evaluate the fitted degree-`n`
# polynomial's derivative at the sample's true offset within it. That is exact
# for polynomials up to degree `n` at *every* sample — unlike mirror padding,
# which folds the signal at the edge and injects spurious curvature ∝ the true
# boundary slope. (Earlier versions reflected; that was both inexact at edges
# and, on the left, off-by-one — see git history.)
function savitzky_golay(m, n, y; order = 0)
    num_y = length(y)
    # The boundary windows span y[1:2m+1] and y[num_y-2m:num_y]; both fit only
    # when the window fits inside the trace. Guarding here keeps every inner
    # loop `@inbounds`.
    num_y >= 2m + 1 ||
        throw(ArgumentError("trace length $num_y is shorter than the SG window $(2m + 1)"))
    z = -m:m
    J = zeros(2m + 1, n + 1)
    for i = 0:n
        @inbounds J[:, i + 1] .= z .^ i
    end

    # P = pinv(J)'  (size (2m+1)×(n+1)). For a window `d`, the least-squares
    # polynomial coefficients are a = pinv(J)·d, so the order-th derivative at
    # window position `t` is gᵀ·a = (P·g)ᵀ·d, where `g = _sg_deriv_covector`.
    # Thus `P·g(t)` is the length-(2m+1) convolution weight vector for that
    # position. Built once; reused for every sample.
    P = J' \ Matrix{Float64}(I, n + 1, n + 1)

    orders = order isa AbstractVector ? collect(order) : [order]
    Y = zeros(num_y, length(orders))
    for (j, o) in enumerate(orders)
        # Interior: centred window (t=0). `wc` is a concrete Vector, so the
        # allocation-free convolution runs behind a function barrier (#128).
        wc = P * _sg_deriv_covector(n, o, 0)
        _sg_interior!(view(Y, :, j), wc, y, m, num_y)
        # Left boundary: fixed window y[1:2m+1]; sample i sits at offset i-(m+1).
        @inbounds for i in 1:m
            w = P * _sg_deriv_covector(n, o, i - m - 1)
            acc = 0.0
            for k in 1:(2m + 1)
                acc += w[k] * y[k]
            end
            Y[i, j] = acc
        end
        # Right boundary: fixed window y[num_y-2m:num_y]; offset i-(num_y-m).
        base = num_y - 2m - 1
        @inbounds for i in (num_y - m + 1):num_y
            w = P * _sg_deriv_covector(n, o, i - (num_y - m))
            acc = 0.0
            for k in 1:(2m + 1)
                acc += w[k] * y[base + k]
            end
            Y[i, j] = acc
        end
    end
    length(orders) == 1 ? Y[:, 1] : Y
end

# Covector `g` such that `g·a` is the `order`-th derivative, evaluated at window
# position `t`, of the polynomial p(z)=Σ aₖ zᵏ. The k-th term contributes
# dᵒʳᵈᵉʳ/dtᵒʳᵈᵉʳ (tᵏ) = k!/(k-order)! · t^(k-order); `prod` avoids factorial
# overflow and yields 1 for the empty range (order=0). This already folds in the
# derivative factor, so no separate k!·cₖ scaling is needed downstream.
function _sg_deriv_covector(n, order, t)
    g = zeros(n + 1)
    for k in order:n
        g[k + 1] = prod((k - order + 1):k) * float(t)^(k - order)
    end
    g
end

# Interior convolution: each output sample is a centred dot product of `wc`
# against `y`. Function-barriered so it specialises on the concrete element
# types and keeps the inner loop allocation-free.
function _sg_interior!(col, wc, y, m, num_y)
    win = 2m + 1
    @inbounds for i in (m + 1):(num_y - m)
        acc = 0.0
        for k in 1:win
            acc += wc[k] * y[i - m - 1 + k]
        end
        col[i] = acc
    end
    col
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
