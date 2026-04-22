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
