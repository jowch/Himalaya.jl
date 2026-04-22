using Himalaya
using Test

@testset "Himalaya" begin
    include("index.jl")
    include("threshold.jl")
    include("persistence.jl")
    include("sharpness.jl")
end
