using LsqFit

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
