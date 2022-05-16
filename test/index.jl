
@testset "Indexing" begin
    @test all(hasfield.(Index, [:basis, :peaks])) 

    test_peaks = [1, √3, 2]
    indices = indexpeaks(0:0.1:2, test_peaks)
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
    
    @test score(index) == 3

    test_peaks = 1:5
    indices = indexpeaks(0:0.1:2, test_peaks)
    @test length(indices) == 2

    a, b = indices
    @test a != b
    @test issubset(a, a)
    @test !issubset(a, b)
end
