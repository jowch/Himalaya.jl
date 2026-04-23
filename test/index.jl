
@testset "Indexing" begin
    @test all(hasfield.(Index, [:basis, :peaks, :sharpness]))
    @test !hasfield(Index, :prom)

    test_peaks     = [1.0, √3, √3 + eps(), 2.0]
    test_sharpness = ones(length(test_peaks))
    indices = indexpeaks(test_peaks, test_sharpness, 0:0.1:2; gaps = false)
    @test length(indices) == 1

    index = only(indices)
    @test phase(index) == Hexagonal
    @test basis(index) == 1
    @test peaks(index) == test_peaks[[1, 2, 4]]
    @test numpeaks(index) == 3

    @test predictpeaks(index) == phaseratios(Hexagonal)

    (; d, R²) = fit(index)
    @test round(d; digits = 3) == 7.255
    @test R² == 1

    # score: coverage × consistency
    # Hexagonal has 14 ratios; peaks at ranks [1,2,3], uniform sharpness → consistency=1
    # coverage = (1/1 + 1/2 + 1/3) / sum(1/r for r in 1:14) ≈ 1.8333 / 3.2515 ≈ 0.564
    @test isapprox(score(index), 0.564; atol = 0.01)
    @test 0.0 ≤ score(index) ≤ 1.0

    # sharpness consistency: heterogeneous peaks score lower
    test_peaks2     = [1.0, √3, 2.0]
    test_sharpness2 = [1.0, 10.0, 1.0]   # mixed wide/sharp
    indices2 = indexpeaks(test_peaks2, test_sharpness2, 0:0.1:2; gaps = false)
    @test length(indices2) == 2  # both Lamellar and Hexagonal survive (not subsets of each other)
    hex_idx = filter(idx -> phase(idx) == Hexagonal, indices2)[1]
    @test score(hex_idx) < score(index)  # Hexagonal with mixed sharpness scores lower than uniform

    test_peaks3 = 1:5
    indices3 = indexpeaks(test_peaks3)
    @test length(indices3) == 2

    a, b = indices3
    @test a != b
    @test issubset(a, a)
    @test !issubset(a, b)

    # End-to-end: top-scoring index is a supported phase, score in [0,1]
    A = readdlm(joinpath(@__DIR__, "data", "example_tot.dat"))
    q, I, σ = A[:, 1], A[:, 2], A[:, 3]
    pk = findpeaks(q, I, σ)
    @test !isempty(pk.q)
    indices4 = indexpeaks(pk.q, pk.sharpness)
    @test !isempty(indices4)
    top = first(sort(indices4; by = score, rev = true))
    @test phase(top) in (Pn3m, Im3m, Ia3d, Lamellar, Hexagonal, Square, Fm3m, Fd3m)
    @test 0.0 ≤ score(top) ≤ 1.0
end
