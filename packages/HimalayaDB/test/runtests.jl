using Test
using HimalayaDB
include("fixture.jl")

@testset "HimalayaDB" begin
    include("test_connect.jl")
end
