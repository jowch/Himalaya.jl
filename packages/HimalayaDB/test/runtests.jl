using Test
using HimalayaDB
include("fixture.jl")

@testset "HimalayaDB" begin
    include("test_connect.jl")
    include("test_queries.jl")
end
