
@testset "Indexing" begin
    @test all(hasfield.(Index, [:basis, :peaks])) 

    test_peaks = [1, √3, 2]
    indices = indexpeaks(test_peaks, ones(length(test_peaks)), 0:0.1:2; gaps = false)
    @test length(indices) == 1

    index = only(indices)
    @test phase(index) == Hexagonal
    @test basis(index) == 1
    @test peaks(index) == test_peaks
    @test numpeaks(index) == 3

    @test predictpeaks(index) == phaseratios(Hexagonal)

    d, R² = fit(index)
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


    test_peaks = [√2, √3]
    indices = indexpeaks(test_peaks, ones(length(test_peaks)), test_peaks)
    @test isempty(indices)
    indices = indexpeaks(test_peaks, ones(length(test_peaks)), test_peaks; req_min = false)
    @test !isempty(indices)
end
