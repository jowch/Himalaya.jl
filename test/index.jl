
@testset "Indexing" begin
    @test all(hasfield.(Index, [:basis, :peaks, :observed])) 

    peaks = 1:5
    indices = indexpeaks(1:5, peaks)
end
