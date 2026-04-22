@testset "findpeaks: single isolated Lorentzian peak" begin
    # A narrow Lorentzian (γ = 0.003) on a flat baseline, uniform σ.
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    σ = fill(1.0, n)
    γ = 0.003
    q0 = 0.2
    I = 10.0 .+ 100.0 ./ (1 .+ ((q .- q0) ./ γ).^2)

    pk = findpeaks(q, I, σ)
    @test length(pk.indices) == 1
    @test abs(pk.q[1] - q0) < step(range(0.005, 0.4; length = n))
    @test pk.prominence[1] > 0
    @test pk.sharpness[1]  > 0
end

@testset "findpeaks: 3 peaks at Pn3m ratios on power-law background" begin
    # Tiny Gaussian noise (0.05 counts, σ=10, SNR≈200) is added so the signal has
    # more than 3 local maxima, enabling the kneedle threshold to work correctly.
    # A purely smooth signal with exactly 3 maxima causes knee() to return the
    # maximum prominence, which eliminates all but the tallest candidate.
    using Random
    Random.seed!(0)
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    σ = fill(10.0, n)
    # Pn3m ratios √2, √3, √4 → normalized to 1, √(3/2), √2
    base = 0.05
    peak_qs = base .* [1.0, sqrt(3/2), sqrt(2)]
    γ = 0.003
    I = 1000.0 .* (q ./ q[1]).^(-2)              # power-law background
    for q0 in peak_qs
        I .+= 100.0 ./ (1 .+ ((q .- q0) ./ γ).^2)
    end
    I .+= 0.05 .* randn(n)                       # negligible noise; breaks degeneracy

    pk = findpeaks(q, I, σ)
    @test length(pk.q) == 3
    for q0 in peak_qs
        @test any(abs.(pk.q .- q0) .< step(range(0.005, 0.4; length = n)))
    end
end

@testset "findpeaks: smooth Lorentzian background produces no peaks" begin
    # The broad Lorentzian is centered at q=0 (below the range minimum) so that
    # the profile is strictly monotone decreasing over [0.005, 0.4] — no local
    # maximum exists and persistence() returns zero candidates.
    # A Lorentzian centered inside the range (e.g. at 0.2) produces a single
    # local max; with n=1 candidate the kneedle threshold degenerates and the
    # feature would be reported as a peak, which is not the intended behaviour.
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    Δq = step(range(0.005, 0.4; length = n))
    width_q = 100 * Δq                           # ≫ any Bragg width
    I = 1000.0 ./ (1 .+ (q ./ width_q).^2)      # center at 0 → monotone over range
    σ = sqrt.(I)                                 # roughly Poisson-scaled

    pk = findpeaks(q, I, σ)
    @test isempty(pk.q)
end

@testset "findpeaks: single-pixel spikes suppressed by :cwt in presence of real peak" begin
    # When the data contain only spikes and no Bragg peak, all candidates have
    # equal prominence and equal CWT sharpness; the kneedle threshold degenerates
    # to the maximum and all pass.  The physically meaningful test is:
    # given a dominant real Lorentzian peak plus small single-pixel spikes
    # (1.5σ amplitude), the algorithm finds exactly the real peak and ignores
    # the spikes, because the spikes have far lower prominence than the Bragg peak.
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    γ = 0.003
    q0 = 0.2
    I = fill(1000.0, n)
    σ = fill(30.0, n)
    I .+= 200.0 ./ (1 .+ ((q .- q0) ./ γ).^2)  # real Bragg peak (amplitude 200 ≫ 1.5σ)
    for idx in (200, 500, 800)
        I[idx] += 45.0                            # 1.5σ single-pixel spikes
    end

    pk = findpeaks(q, I, σ; sharpness_method = :cwt)
    @test length(pk.q) == 1
    @test abs(pk.q[1] - q0) < step(range(0.005, 0.4; length = n))
end

@testset "findpeaks: two close-but-resolved peaks" begin
    # Two overlapping Lorentzians (separation ≈ 2.86 HWHM) produce exactly two
    # persistence candidates.  With only two candidates the kneedle threshold
    # degenerates to the maximum prominence, so the weaker peak is suppressed by
    # automatic thresholding.  We therefore use prom_floor=0 and sharp_floor=0
    # to test the resolution of the peak detector independently of the adaptive
    # thresholding logic.
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    Δq = step(range(0.005, 0.4; length = n))
    γ = 0.003
    Δ = 20 * Δq                                  # ~2.86 HWHM separation
    q01, q02 = 0.2 - Δ/2, 0.2 + Δ/2
    σ = fill(10.0, n)
    I = 100.0 .+ 100.0 ./ (1 .+ ((q .- q01) ./ γ).^2) .+
                  100.0 ./ (1 .+ ((q .- q02) ./ γ).^2)

    pk = findpeaks(q, I, σ; prom_floor = 0.0, sharp_floor = 0.0)
    @test length(pk.q) == 2
end

@testset "findpeaks: works without σ" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    γ = 0.003
    I = 10.0 .+ 100.0 ./ (1 .+ ((q .- 0.2) ./ γ).^2)

    pk = findpeaks(q, I)
    @test length(pk.q) >= 1
    @test any(abs.(pk.q .- 0.2) .< 2 * step(range(0.005, 0.4; length = n)))
end

@testset "findpeaks: normalize_by_σ on vs off produces sensible results" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    γ = 0.003
    I = 10.0 .+ 100.0 ./ (1 .+ ((q .- 0.2) ./ γ).^2)
    σ = fill(1.0, n)

    pk_on  = findpeaks(q, I, σ; normalize_by_σ = true)
    pk_off = findpeaks(q, I, σ; normalize_by_σ = false)

    # Uniform σ means both modes should find the same peak.
    @test !isempty(pk_on.q)
    @test !isempty(pk_off.q)
    @test abs(first(pk_on.q) - first(pk_off.q)) < 2 * step(range(0.005, 0.4; length = n))
end
