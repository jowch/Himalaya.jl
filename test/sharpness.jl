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

@testset "savitzky_golay: allocation-free inner loop" begin
    # Regression guard for issue #128. The pre-fix implementation allocated a
    # length-(2m+1) index vector *and* a gathered window *per inner-loop
    # iteration*, so total allocation scaled with num_y × window (~280 B per
    # sample). After the fix the inner loop allocates nothing — the only
    # per-sample cost is the O(num_y) output vector itself.
    #
    # Measuring the *marginal* allocation between two trace lengths cancels
    # the constant one-time setup (the Vandermonde matrix + linear solve) and
    # isolates exactly the per-sample behaviour the issue is about.
    sg = Himalaya.savitzky_golay
    y1 = collect(range(0.0, 1.0; length = 1000)) .^ 2
    y2 = collect(range(0.0, 1.0; length = 2000)) .^ 2
    sg(5, 4, y1; order = 2); sg(5, 4, y2; order = 2)  # warm up / compile
    a1 = @allocated sg(5, 4, y1; order = 2)
    a2 = @allocated sg(5, 4, y2; order = 2)
    per_sample = (a2 - a1) / 1000
    @test per_sample < 64
end

@testset "savitzky_golay: rejects a trace shorter than the window" begin
    # The edge-reflection index math is only in-bounds when num_y >= 2m+1.
    # A shorter trace must fail loudly rather than reach the @inbounds loop.
    @test_throws ArgumentError Himalaya.savitzky_golay(5, 4, rand(8); order = 2)
    @test Himalaya.savitzky_golay(5, 4, rand(11); order = 2) isa Vector  # 2m+1 ok
end

@testset "savitzky_golay: numerical output pinned on a fixture trace" begin
    A = readdlm(joinpath(@__DIR__, "data", "example_tot.dat"))
    intensity = A[:, 2]
    s = Himalaya.sharpness(intensity; method = :savgol, m = 5)
    @test length(s) == 922
    # Pinned values spanning boundary-low, interior, and boundary-high regions.
    pinned = Dict(
        1   => 3579.4866724945023,
        2   => 954.0417132869541,
        6   => -792.756275058146,
        11  => -209.51359498828796,
        100 => -0.4906289510479224,
        461 => 0.026168164336535318,
        800 => 0.10356053613089689,
        917 => -0.2534995337991986,
        921 => -0.009534947552105598,
        922 => -0.14693493006959005,
    )
    for (idx, val) in pinned
        @test isapprox(s[idx], val; rtol = 1e-12)
    end
    @test isapprox(sum(s), -3113.2448241053053; rtol = 1e-12)
end
