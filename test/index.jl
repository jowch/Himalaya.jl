
@testset "Indexing" begin
    @test all(hasfield.(Index, [:basis, :peaks])) 

    test_peaks = [1, √3, √3 + eps(), 2]
    indices = indexpeaks(test_peaks, ones(length(test_peaks)), 0:0.1:2; gaps = false)
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
    
    @test isapprox(score(index), 3; atol = 1e10)

    test_peaks = 1:5
    indices = indexpeaks(test_peaks)
    @test length(indices) == 2

    a, b = indices
    @test a != b
    @test issubset(a, a)
    @test !issubset(a, b)

    # End-to-end: load a real trace, run findpeaks (v2), run indexpeaks,
    # assert the top-scoring index is in the supported phase set.
    A = readdlm(joinpath(@__DIR__, "data", "example_tot.dat"))
    q, I, σ = A[:, 1], A[:, 2], A[:, 3]
    pk = findpeaks(q, I, σ)
    @test !isempty(pk.q)
    indices = indexpeaks(pk.q, pk.prominence)
    @test !isempty(indices)
    top = first(sort(indices; by = score, rev = true))
    @test phase(top) in (Pn3m, Im3m, Ia3d, Lamellar, Hexagonal, Square, Fm3m, Fd3m)

end
