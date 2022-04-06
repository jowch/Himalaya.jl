
@testset "Indexing" begin
    @test all(hasfield.(Index, [:phase, :basis, :peaks])) 

    peaks = 1:5
    indices = indexpeaks(peaks)

end
