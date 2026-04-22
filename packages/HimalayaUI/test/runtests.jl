using Test

@testset "HimalayaUI" begin
    include("test_db.jl")
    include("test_datfile.jl")
    include("test_manifest.jl")
end
