@testset "persistence" begin
    # Single Gaussian peak on a flat baseline (baseline = 1, peak height = 10).
    n = 101
    xs = 1:n
    y = [1.0 + 9.0 * exp(-((x - 50)^2) / (2 * 5^2)) for x in xs]

    p = Himalaya.persistence(y)

    @test length(p.indices) == 1
    @test p.indices[1] == 50
    # Prominence = peak height above flat baseline = 9.0, within numerical tolerance.
    @test isapprox(only(p.prominence), 9.0; atol = 1e-6)

    # Two peaks: heights 10 and 6 at indices 30 and 70, with gradual descent.
    y2 = [1.0 + 0.01*i for i in 1:100]
    y2[30] = 10.0; y2[29] = 6.0; y2[31] = 6.0
    y2[70] =  6.0; y2[69] = 3.0; y2[71] = 3.0

    p2 = Himalaya.persistence(y2)
    @test length(p2.indices) == 2
    @test 30 in p2.indices
    @test 70 in p2.indices

    # Monotonic signal has no local maxima (endpoints excluded by convention).
    p3 = Himalaya.persistence(collect(1.0:100.0))
    @test isempty(p3.indices)
    @test isempty(p3.prominence)
end
