using Test
using HimalayaDB
include("fixture.jl")

@testset "HimalayaDB" begin
    include("test_connect.jl")
    include("test_queries.jl")
    include("test_reconstruct.jl")
    include("test_trace.jl")
    include("test_dataframes.jl")
end
