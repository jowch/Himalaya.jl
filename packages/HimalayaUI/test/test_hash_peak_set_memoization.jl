using Test
using HimalayaUI: hash_peak_set
import HimalayaUI

@testset "hash_peak_set: cache hit returns same hash without recomputing" begin
    eff = (q = [0.1, 0.2, 0.3], sharpness = [1.0, 0.8, 0.6])
    h1 = hash_peak_set(eff)
    h2 = hash_peak_set(eff)
    @test h1 == h2

    # Mutating inputs invalidates: same vector identity but different content.
    eff_mut = (q = copy(eff.q), sharpness = copy(eff.sharpness))
    push!(eff_mut.q, 0.4)
    push!(eff_mut.sharpness, 0.5)
    h3 = hash_peak_set(eff_mut)
    @test h3 != h1
end

@testset "hash_peak_set: cache works across separate NamedTuple constructions" begin
    HimalayaUI._clear_peak_set_hash_cache!()
    h1 = hash_peak_set((q = [0.11, 0.22], sharpness = [1.0, 0.9]))
    n_after_first = length(HimalayaUI._PEAK_SET_HASH_CACHE)
    h2 = hash_peak_set((q = [0.11, 0.22], sharpness = [1.0, 0.9]))
    n_after_second = length(HimalayaUI._PEAK_SET_HASH_CACHE)
    @test h1 == h2
    @test n_after_second == n_after_first  # second call hit the cache, didn't insert
end

@testset "hash_peak_set: cache cap" begin
    HimalayaUI._clear_peak_set_hash_cache!()
    # Pump > cap entries to confirm cap holds.
    cap = HimalayaUI.PEAK_SET_HASH_CACHE_CAP
    for i in 1:(cap + 50)
        hash_peak_set((q = [float(i)], sharpness = [1.0]))
    end
    @test length(HimalayaUI._PEAK_SET_HASH_CACHE) <= cap
end
