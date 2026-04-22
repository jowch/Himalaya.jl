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
    fit_peak(q, I, σ, centre_idx, width_pts; shape = :lorentzian)

Locally fit a peak + linear baseline to the trace `(q, I)` with
per-point uncertainties `σ`, in a window of `±2.5 · width_pts` points
around `centre_idx`. Uses weighted Levenberg-Marquardt (`LsqFit.curve_fit`)
with weights `1/σᵢ²`.

`shape` selects the peak model:
- `:lorentzian` (default): `A / (1 + ((x - q0) / w)²) + b + m·x`
- `:gaussian`: `A · exp(-(x - q0)² / (2w²)) + b + m·x`

Returns a NamedTuple `(; A, q0, fwhm, snr)` on success:
- `A`: fitted peak amplitude (above baseline)
- `q0`: fitted peak centre in q-units
- `fwhm`: full-width-at-half-maximum in q-units (shape-aware)
- `snr`: `A / σ_A` where `σ_A` is the standard error on `A` from the
  covariance matrix

Returns `nothing` on convergence failure or non-physical fits (A ≤ 0,
w ≤ 0, or window too small). Throws `ArgumentError` for unknown `shape`.
"""
function fit_peak(q, I, σ, centre_idx, width_pts; shape = :lorentzian)
    if shape !== :gaussian && shape !== :lorentzian
        throw(ArgumentError("Unknown peak shape $(repr(shape)). Must be :gaussian or :lorentzian."))
    end

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
    w_init = if shape === :gaussian
        max(width_pts * Δq / 2.355, Δq)
    else
        max(width_pts * Δq / 2, Δq)
    end
    b0 = (Iw[1] + Iw[end]) / 2
    m0 = (Iw[end] - Iw[1]) / (qw[end] - qw[1])

    model = if shape === :gaussian
        (x, p) -> p[1] .* exp.(-((x .- p[2]).^2) ./ (2 * p[3]^2)) .+ p[4] .+ p[5] .* x
    else
        (x, p) -> p[1] ./ (1 .+ ((x .- p[2]) ./ p[3]).^2) .+ p[4] .+ p[5] .* x
    end
    p0 = [A0, q0_init, w_init, b0, m0]
    weights = 1.0 ./ (σw_arr .^ 2)

    fit = try
        curve_fit(model, qw, Iw, weights, p0)
    catch
        return nothing
    end

    fit.converged || return nothing
    A_fit, q0_fit, w_fit = fit.param[1], fit.param[2], fit.param[3]
    if A_fit <= 0 || w_fit <= 0
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

    fwhm = shape === :gaussian ? 2.355 * w_fit : 2 * w_fit
    (A = A_fit, q0 = q0_fit, fwhm = fwhm, snr = A_fit / σ_A)
end

const DEFAULT_SCALES = [2.0, 2.6, 3.5, 4.5, 6.0, 8.0, 10.0, 13.0, 17.0, 22.0, 28.0, 35.0]

"""
    findpeaks(q, I, σ; nσ = 5, scales = nothing, min_ridge_length = 3, shape = :lorentzian)

Detect Bragg peaks in a SAXS integration trace.

# Arguments
- `q`, `I`, `σ`: equal-length vectors from a `_tot.dat` file. `σ` is the
  per-point intensity uncertainty (third column).
- `nσ`: SNR threshold in sigmas. A candidate is accepted iff `A / σ_A ≥ nσ`,
  where `A` is the fitted peak amplitude and `σ_A` its standard error.
  Default 5.
- `scales`: vector of CWT widths in *points*. `nothing` ⇒ a default
  geometric series from 2 to 35 points covering typical Bragg widths.
- `min_ridge_length`: a CWT ridge must persist across this many adjacent
  scales to be considered. Default 3.
- `shape`: peak shape model used in Stage 2 fit, `:lorentzian` (default) or
  `:gaussian`. Lipid mesophase Bragg peaks are typically Lorentzian; use
  `:gaussian` for samples dominated by strain or instrumental broadening.
- `reject_edge_ridges`: when `true` (default), reject candidates whose CWT
  ridge matches best at the largest scale (a sign the feature is wider than
  the searched range — typically broad form-factor background, not a Bragg
  peak). A summary `@warn` is emitted listing how many were rejected so the
  signal isn't silent. Set `false` to keep edge-matched ridges.

# Returns
A NamedTuple `(; indices, q, snr, width)` of equal-length vectors:
- `indices`: nearest grid indices into the input arrays
- `q`: sub-pixel peak positions from the local fit
- `snr`: physical SNR (`A / σ_A`)
- `width`: fitted FWHM in q-units

Peaks are returned sorted by ascending `q`.
"""
function findpeaks(q, I, σ;
                   nσ = 5, scales = nothing, min_ridge_length = 3,
                   shape = :lorentzian, reject_edge_ridges = true)
    @assert length(q) == length(I) == length(σ) "q, I, σ must be the same length"
    sc = scales === nothing ? DEFAULT_SCALES : collect(float.(scales))

    # Stage 1: CWT on log10(I) (compress dynamic range so high-q peaks register)
    y = log10.(max.(I, eps()))
    coeffs = cwt(y, sc)
    ridges = find_ridges(coeffs, sc; min_ridge_length = min_ridge_length)

    # Stage 2: local fit + σ-based SNR validation
    indices = Int[]
    qs = Float64[]
    snrs = Float64[]
    widths = Float64[]
    n_edge_rejected = 0

    for r in ridges
        # Scale-edge rejection: a ridge that matches best at the largest
        # available scale has likely "railed the meter" — its natural width
        # exceeds what we searched for, so it's a broad background feature,
        # not a Bragg peak.
        if reject_edge_ridges && r.scale == sc[end]
            n_edge_rejected += 1
            continue
        end

        width_pts = r.scale
        fit = fit_peak(q, I, σ, r.index, width_pts; shape = shape)
        fit === nothing && continue
        fit.snr < nσ && continue

        # Snap fitted q0 back to nearest grid index
        idx = argmin(abs.(q .- fit.q0))
        push!(indices, idx)
        push!(qs, fit.q0)
        push!(snrs, fit.snr)
        push!(widths, fit.fwhm)
    end

    if n_edge_rejected > 0
        @warn "findpeaks: rejected $n_edge_rejected candidate(s) at the largest CWT scale ($(sc[end]) pts) — likely broad background features. Extend `scales` upward or pass `reject_edge_ridges = false` if these are genuine peaks."
    end

    # Sort by q
    perm = sortperm(qs)
    (indices = indices[perm], q = qs[perm], snr = snrs[perm], width = widths[perm])
end
