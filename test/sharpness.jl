@testset "ricker wavelet" begin
    a = 5.0
    @test Himalaya.ricker(0.0, a) > 0
    @test isapprox(Himalaya.ricker(a, a), 0.0; atol = 1e-12)
    @test Himalaya.ricker(2a, a) < 0
    @test Himalaya.ricker(1.3, a) == Himalaya.ricker(-1.3, a)
end

@testset "convolve_reflect" begin
    # Reflecting convolution of a delta with a known kernel reproduces the kernel.
    y = zeros(21); y[11] = 1.0
    k = [1.0, 2.0, 3.0, 2.0, 1.0]
    out = Himalaya.convolve_reflect(y, k)
    @test out[11] == 3.0
    @test out[10] == 2.0
    @test out[12] == 2.0
    @test out[9]  == 1.0
    @test out[13] == 1.0
end

@testset "savitzky_golay reproduces a polynomial's derivative" begin
    xs = collect(-10.0:10.0)
    y = xs.^3
    d2 = Himalaya.savitzky_golay(5, 4, y; order = 2)
    # d²/dx² x³ = 6x; check interior points (away from edges).
    for i in 6:16
        @test isapprox(d2[i], 6 * xs[i]; atol = 1e-6)
    end
end

@testset "sharpness: :savgol picks up a Gaussian peak" begin
    xs = 1:201
    y = [exp(-((x - 101)^2) / (2 * 5^2)) for x in xs]
    s = Himalaya.sharpness(y; method = :savgol, m = 5)
    @test length(s) == length(y)
    @test argmax(s) == 101
    @test s[101] > 0
    @test abs(s[1])   < 0.01
    @test abs(s[end]) < 0.01
end

@testset "sharpness: :cwt picks up a Gaussian peak" begin
    xs = 1:201
    y = [exp(-((x - 101)^2) / (2 * 5^2)) for x in xs]
    s = Himalaya.sharpness(y; method = :cwt)
    @test length(s) == length(y)
    @test argmax(s) == 101
    @test s[101] > 0
end

@testset "sharpness: invalid method throws" begin
    y = rand(100)
    @test_throws ArgumentError Himalaya.sharpness(y; method = :bogus)
end
