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
    # d²/dx² x³ = 6x. The cubic has degree ≤ n=4, so the "interp" boundary
    # scheme is exact at EVERY sample, edges included (mirror padding was not).
    for i in eachindex(xs)
        @test isapprox(d2[i], 6 * xs[i]; atol = 1e-6)
    end
end

@testset "savitzky_golay: interp edges are exact for degree ≤ n polynomials" begin
    # Acceptance test for the "interp" boundary scheme: for a polynomial of
    # degree ≤ n, the order-d derivative estimate must be exact at every sample
    # including the first/last m — there is no reflection assumption to violate.
    x = collect(1.0:1.0:30.0)
    c = [1.0, 2.0, -0.5, 0.05, 0.001]                     # a₀..a₄, degree 4 = n
    y  = [sum(c[k + 1] * xi^k for k in 0:4) for xi in x]
    d0 = y                                                 # order 0 ⇒ smoothed == input
    d1 = [sum(k * c[k + 1] * xi^(k - 1) for k in 1:4) for xi in x]
    d2 = [sum(k * (k - 1) * c[k + 1] * xi^(k - 2) for k in 2:4) for xi in x]
    for (o, truth) in ((0, d0), (1, d1), (2, d2))
        est = Himalaya.savitzky_golay(5, 4, y; order = o)
        for i in eachindex(x)                              # edges included
            @test isapprox(est[i], truth[i]; atol = 1e-6)
        end
    end
end

@testset "savitzky_golay: edge handling preserves input symmetry" begin
    # A signal symmetric about its midpoint must produce a symmetric second
    # derivative — a convention-agnostic correctness property, not a pinned
    # magic number. Regression guard against asymmetric edge handling: an
    # earlier mirror-padding version mis-centred the left window (folding
    # half-sample about 0.5) and shifted the first m outputs one sample left,
    # breaking this symmetry. The interp scheme treats both edges identically.
    xs = collect(-15.0:1.0:15.0)   # 31 points, even count flanking the centre
    y = 3 .* xs .^ 2               # even function ⇒ samples symmetric about i=16
    d2 = Himalaya.savitzky_golay(5, 4, y; order = 2)
    n = length(d2)
    for i in 1:n
        @test isapprox(d2[i], d2[n + 1 - i]; atol = 1e-9)
    end
end

@testset "savitzky_golay: left edge stays centred on its own sample" begin
    # Regression guard: an earlier off-by-one made output[m] == output[m+1]
    # exactly (the first m samples were the interior values shifted one index
    # left). On a strictly-curved monotone signal the seam between the boundary
    # samples (i ≤ m) and the interior (i > m) must not duplicate a sample.
    m = 5
    y = collect(1.0:1.0:40.0) .^ 2
    d2 = Himalaya.savitzky_golay(m, 4, y; order = 2)
    @test d2[m] != d2[m + 1]
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
    # The boundary pins (1, 2, 921, 922) are the "interp" boundary values; the
    # interior pins (6 ≤ idx ≤ 917) are convention-independent and unchanged.
    pinned = Dict(
        1   => -2122.6039906770857,
        2   => -1821.240713287396,
        6   => -792.756275058146,
        11  => -209.51359498828796,
        100 => -0.4906289510479224,
        461 => 0.026168164336535318,
        800 => 0.10356053613089689,
        917 => -0.2534995337991986,
        921 => 1.0763154778542567,
        922 => 1.7567832167810487,
    )
    for (idx, val) in pinned
        @test isapprox(s[idx], val; rtol = 1e-12)
    end
    @test isapprox(sum(s), -11117.585612599207; rtol = 1e-12)
end
