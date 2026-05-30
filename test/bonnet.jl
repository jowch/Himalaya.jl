using Test
using Himalaya
using Himalaya: Pn3m, Im3m, Ia3d, Lamellar, bonnet_lattice

@testset "bonnet_lattice predicts the coexisting cubic" begin
    # a_Im3m ≈ 1.279 · a_Pn3m  (the canonical ratio).
    @test bonnet_lattice(Pn3m, 100.0, Im3m) ≈ 127.9 atol=1e-6
    @test bonnet_lattice(Im3m, 127.9, Pn3m) ≈ 100.0 atol=1e-3

    # Ia3d uses the 1.576 reference (confirm with domain expert before shipping).
    @test bonnet_lattice(Pn3m, 100.0, Ia3d) ≈ 157.6 atol=1e-6

    # Non-bicontinuous phase → nothing (no Bonnet relation defined).
    @test bonnet_lattice(Pn3m, 100.0, Lamellar) === nothing
    @test bonnet_lattice(Lamellar, 100.0, Im3m) === nothing

    # Same phase → identity.
    @test bonnet_lattice(Pn3m, 100.0, Pn3m) ≈ 100.0 atol=1e-9
end

@testset "bonnet_consistent within relative tolerance" begin
    # Observed Im3m at 128 vs predicted 127.9 from a Pn3m at 100 → consistent.
    @test Himalaya.bonnet_consistent(Pn3m, 100.0, Im3m, 128.0; reltol=0.02)
    # Observed Im3m at 140 → off by ~9% → not consistent at 2%.
    @test !Himalaya.bonnet_consistent(Pn3m, 100.0, Im3m, 140.0; reltol=0.02)
    # Undefined pair → false (no relation to be consistent with).
    @test !Himalaya.bonnet_consistent(Pn3m, 100.0, Lamellar, 100.0; reltol=0.02)
end
