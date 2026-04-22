using Test
using HimalayaUI: load_dat

const EXAMPLE_DAT = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")

@testset "load_dat" begin
    q, I, σ = load_dat(EXAMPLE_DAT)

    @test length(q) == 922
    @test length(I) == 922
    @test length(σ) == 922

    @test q[1] ≈ 6.258469e-03 atol=1e-8
    @test I[1] ≈ 4.527260e+04 atol=1.0
    @test σ[1] ≈ 2.127736e+02 atol=1e-2

    @test all(q .> 0)
    @test all(I .> 0)
    @test all(σ .> 0)
    @test issorted(q)
end
