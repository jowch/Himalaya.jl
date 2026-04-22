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
