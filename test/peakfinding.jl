@testset "ricker wavelet" begin
    # Ricker (Mexican hat) is the negative normalized 2nd derivative of a Gaussian.
    # At t=0 it has its maximum; at t=±a it crosses zero; for |t|>a it's negative.
    a = 5.0
    @test Himalaya.ricker(0.0, a) > 0
    @test isapprox(Himalaya.ricker(a, a), 0.0; atol = 1e-12)
    @test isapprox(Himalaya.ricker(-a, a), 0.0; atol = 1e-12)
    @test Himalaya.ricker(2a, a) < 0
    # symmetry
    @test Himalaya.ricker(1.3, a) == Himalaya.ricker(-1.3, a)
    # integral over a wide window should be ~0 (zero-mean wavelet)
    ts = -10a:0.01:10a
    @test abs(sum(Himalaya.ricker.(ts, a)) * 0.01) < 1e-6
end

@testset "cwt" begin
    # A pure Gaussian peak should produce its largest CWT response at the
    # scale matching its width.
    n = 300
    σ_true = 4.0          # Gaussian std in points
    centre = 150
    y = [exp(-((i - centre)^2) / (2σ_true^2)) for i in 1:n]

    scales = [1.5, 2.5, 4.0, 6.5, 10.5]
    coeffs = Himalaya.cwt(y, scales)

    @test size(coeffs) == (n, length(scales))
    # the row at the peak centre should peak in the column nearest σ_true
    centre_row = coeffs[centre, :]
    best_scale_idx = argmax(centre_row)
    @test scales[best_scale_idx] == 4.0
end

@testset "local_maxima" begin
    # A vector with two clear peaks at indices 3 and 7
    v = [0.0, 1.0, 3.0, 1.0, 0.0, 2.0, 5.0, 2.0, 0.0]
    @test Himalaya.local_maxima(v) == [3, 7]

    # Negative values should not produce maxima (the CWT filtering wants positives only)
    v2 = [-1.0, -0.5, -1.0, 0.0, 1.0, 0.0]
    @test Himalaya.local_maxima(v2) == [5]

    # Endpoints don't count
    v3 = [3.0, 1.0, 2.0, 1.0, 3.0]
    @test Himalaya.local_maxima(v3) == [3]
end

@testset "find_ridges" begin
    # Synthesize a CWT coefficient matrix where one Gaussian peak at
    # column-position 50 produces ridges across scales.
    n = 100
    scales = [1.5, 2.5, 4.0, 6.5, 10.5]
    σ_true = 4.0
    y = [exp(-((i - 50)^2) / (2σ_true^2)) for i in 1:n]
    coeffs = Himalaya.cwt(y, scales)

    ridges = Himalaya.find_ridges(coeffs, scales; min_ridge_length = 3)

    # Exactly one ridge survives, centred near i=50
    @test length(ridges) == 1
    r = ridges[1]
    @test abs(r.index - 50) <= 1
    @test r.scale == 4.0     # best-matched scale
    @test r.length >= 3
end

@testset "find_ridges rejects single-scale noise" begin
    # A noisy CWT with maxima only at the smallest scale should produce no ridges.
    n = 100
    scales = [1.5, 2.5, 4.0, 6.5, 10.5]
    coeffs = zeros(n, length(scales))
    coeffs[20, 1] = 1.0   # one isolated max at smallest scale only
    coeffs[60, 1] = 1.0   # another isolated max
    ridges = Himalaya.find_ridges(coeffs, scales; min_ridge_length = 3)
    @test isempty(ridges)
end

using Statistics: median

@testset "fit_peak: gaussian" begin
    # Synthesize a clean Gaussian + linear baseline + Gaussian noise.
    # Then assert fit_peak recovers amplitude and a high SNR.
    n = 60
    q_range = range(0.0, 1.0; length = n)
    q = collect(q_range)
    dq = step(q_range)
    A_true = 100.0
    q0_true = 0.5
    σw_true = 0.05
    baseline = 5.0 .+ 2.0 .* q
    peak = A_true .* exp.(-((q .- q0_true).^2) ./ (2σw_true^2))
    σ_arr = fill(1.0, n)
    I = peak .+ baseline   # noiseless to make assertions tight
    centre_idx = argmin(abs.(q .- q0_true))

    fit = Himalaya.fit_peak(q, I, σ_arr, centre_idx, σw_true / dq; shape = :gaussian)

    @test fit !== nothing
    @test isapprox(fit.A, A_true; rtol = 0.01)
    @test isapprox(fit.q0, q0_true; atol = dq)
    @test isapprox(fit.fwhm, 2.355 * σw_true; rtol = 0.05)
    # SNR should be very high for a noiseless peak with σ=1
    @test fit.snr > 50
end

@testset "fit_peak: lorentzian" begin
    n = 60
    q_range = range(0.0, 1.0; length = n)
    q = collect(q_range)
    dq = step(q_range)
    A_true = 100.0
    q0_true = 0.5
    γ_true = 0.05
    baseline = 5.0 .+ 2.0 .* q
    peak = A_true ./ (1 .+ ((q .- q0_true) ./ γ_true).^2)
    σ_arr = fill(1.0, n)
    I = peak .+ baseline
    centre_idx = argmin(abs.(q .- q0_true))

    fit = Himalaya.fit_peak(q, I, σ_arr, centre_idx, γ_true / dq)  # default shape=:lorentzian

    @test fit !== nothing
    @test isapprox(fit.A, A_true; rtol = 0.01)
    @test isapprox(fit.q0, q0_true; atol = dq)
    @test isapprox(fit.fwhm, 2 * γ_true; rtol = 0.05)
    @test fit.snr > 50
end

@testset "fit_peak: invalid shape throws" begin
    n = 60
    q = collect(range(0.0, 1.0; length = n))
    I = collect(1.0:n)
    σ_arr = fill(1.0, n)
    @test_throws ArgumentError Himalaya.fit_peak(q, I, σ_arr, 30, 5.0; shape = :triangle)
end

@testset "fit_peak returns nothing on failure" begin
    # Pathological window: flat data, no peak. LsqFit may not converge or
    # may converge to A≈0 with huge σ_A. fit_peak should return nothing
    # in either case (we treat both as "no detectable peak here").
    n = 60
    q = collect(range(0.0, 1.0; length = n))
    I = fill(10.0, n)
    σ_arr = fill(1.0, n)
    fit = Himalaya.fit_peak(q, I, σ_arr, 30, 5.0)
    @test fit === nothing || fit.snr < 1
end

@testset "findpeaks: single isolated peak" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    A = 100.0
    q0 = 0.2
    γ = 0.003
    σ_arr = fill(5.0, n)   # ⇒ A/σ = 20 at the peak
    # Use a Lorentzian since that's the new default shape
    I = 50.0 .+ A ./ (1 .+ ((q .- q0) ./ γ).^2)

    pk = findpeaks(q, I, σ_arr)
    @test length(pk.q) == 1
    @test abs(pk.q[1] - q0) < step(range(0.005, 0.4; length = n))
    @test pk.snr[1] > 10
end

@testset "findpeaks: 3 peaks at Pn3m ratios on power-law background" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    # Pn3m first three normalized ratios: √2, √3, √4 → divide by √2 → 1, √(3/2), √2
    base_q = 0.05
    peak_qs = base_q .* [1.0, sqrt(3/2), sqrt(2)]
    γ = 0.003
    background = 1000.0 .* (q ./ q[1]).^(-2)   # power law
    I = copy(background)
    for q0 in peak_qs
        I .+= 100.0 ./ (1 .+ ((q .- q0) ./ γ).^2)
    end
    σ_arr = fill(10.0, n)   # A/σ ≈ 10 at peak

    pk = findpeaks(q, I, σ_arr)
    @test length(pk.q) == 3
    for q0 in peak_qs
        @test any(abs.(pk.q .- q0) .< step(range(0.005, 0.4; length = n)))
    end
end

@testset "findpeaks: below-threshold peak rejected" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    γ = 0.003
    σ_arr = fill(50.0, n)
    I = 100.0 .+ 100.0 ./ (1 .+ ((q .- 0.2) ./ γ).^2)   # A/σ ≈ 2

    pk = findpeaks(q, I, σ_arr; nσ = 5)
    @test isempty(pk.q)
end

@testset "findpeaks: just-above-threshold peak detected" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    γ = 0.003
    σ_arr = fill(15.0, n)
    I = 100.0 .+ 100.0 ./ (1 .+ ((q .- 0.2) ./ γ).^2)   # A/σ ≈ 6

    pk = findpeaks(q, I, σ_arr; nσ = 5)
    @test length(pk.q) == 1
end

@testset "findpeaks: smooth Lorentzian background produces no peaks" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    # Wide Lorentzian (width ~100 points in q-units)
    width_q = 100 * step(range(0.005, 0.4; length = n))
    background = 1.0 ./ (1 .+ ((q .- 0.2) ./ width_q).^2)
    I = 1000.0 .* background
    σ_arr = sqrt.(I)   # roughly Poisson-scaled

    pk = findpeaks(q, I, σ_arr)
    @test isempty(pk.q)
end

@testset "findpeaks: single-pixel spikes do not trigger" begin
    using Random
    Random.seed!(0)
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    I = fill(1000.0, n)
    σ_arr = fill(30.0, n)
    # Three single-pixel outliers, each ~10σ above baseline
    for idx in (200, 500, 800)
        I[idx] = 1300.0
    end

    pk = findpeaks(q, I, σ_arr)
    @test isempty(pk.q)
end

@testset "findpeaks: two close-but-resolved peaks" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    γ = 0.003
    # Lorentzian peaks need ≳2-3 HWHM separation to be resolved (heavier tails
    # than Gaussian fill the dip between close peaks). γ=0.003 ≈ 7 q-grid points,
    # so 20 grid points apart is ~3 HWHM — clearly resolvable.
    Δ = 20 * step(range(0.005, 0.4; length = n))
    q01, q02 = 0.2 - Δ/2, 0.2 + Δ/2
    σ_arr = fill(10.0, n)
    I = 100.0 .+ 100.0 ./ (1 .+ ((q .- q01) ./ γ).^2) .+
                  100.0 ./ (1 .+ ((q .- q02) ./ γ).^2)

    pk = findpeaks(q, I, σ_arr)
    @test length(pk.q) == 2
end
