using Test
using Himalaya
using Himalaya: indexpeaks, Index, peaks, score
using HimalayaUI: auto_group
using SparseArrays

@testset "auto_group claims by peak_id, not float equality" begin
    qs        = [0.1000, 0.1414, 0.1732, 0.2000, 0.2236, 0.2449]
    sharpness = [1.0, 0.9, 0.85, 0.8, 0.7, 0.65]
    peak_ids  = [101, 102, 103, 104, 105, 106]
    eff = (q = qs, sharpness = sharpness, peak_id = peak_ids,
           peak_kind = fill(:auto, length(qs)))

    candidates = indexpeaks(eff.q, eff.sharpness)
    @test !isempty(candidates)

    group = auto_group(candidates, eff)

    # Convert each idx's q-values back to peak_ids via eff.
    q_to_id = Dict(eff.q[i] => eff.peak_id[i] for i in eachindex(eff.q))
    id_sets = [Set(q_to_id[q] for q in peaks(idx)) for idx in group]
    for i in eachindex(id_sets), j in eachindex(id_sets)
        i == j && continue
        @test isempty(intersect(id_sets[i], id_sets[j]))
    end
end

@testset "auto_group(indices) single-arg form remains q-based" begin
    # The legacy single-arg signature MUST keep working for fixtures that
    # don't materialize an eff tuple (e.g. test_pipeline.jl::auto_group).
    qs        = [0.1000, 0.1414, 0.2000]
    sharpness = [1.0, 0.8, 0.6]
    candidates = indexpeaks(qs, sharpness)
    group = auto_group(candidates)
    if !isempty(candidates)
        @test !isempty(group)
    end
end

@testset "auto_group members never share a peak" begin
    qs        = [0.1, 0.1414, 0.1732, 0.2, 0.2449]
    sharpness = ones(length(qs))
    candidates = indexpeaks(qs, sharpness)
    @test !isempty(candidates)
    eff = (q = qs, sharpness = sharpness,
           peak_id = collect(1:length(qs)),
           peak_kind = fill(:auto, length(qs)))
    for group in (auto_group(candidates, eff), auto_group(candidates))
        @test !isempty(group)
        all_q = reduce(vcat, [collect(peaks(i)) for i in group])
        @test length(all_q) == length(unique(all_q))
    end
end

const HEX_Q  = 0.09789 .* [1, sqrt(3), 2, sqrt(7), 3, sqrt(12), sqrt(13)]
const HEX_SH = [4.57, 0.747, 0.383, 0.111, 0.0543, 0.0205, 0.0397]

_eff(q, sh; ids = collect(1:length(q)), kinds = fill(:auto, length(q))) =
    (q = q, sharpness = sh, peak_id = ids, peak_kind = kinds)
_phases(g) = sort([string(nameof(Himalaya.phase(x))) for x in g])
_both(cands, eff) = (auto_group(cands, eff), auto_group(cands))

@testset "auto_group picks the index that contains the top index and a sharper peak" begin
    cands = collect(indexpeaks(HEX_Q, HEX_SH))
    top = cands[argmax(score.(cands))]
    @test top isa Himalaya.Index{Himalaya.Ia3d}
    for group in _both(cands, _eff(HEX_Q, HEX_SH; ids = collect(11:17)))
        @test first(group) isa Himalaya.Index{Himalaya.Hexagonal}
        @test !any(idx -> idx isa Himalaya.Index{Himalaya.Ia3d}, group)
    end
end

@testset "auto_group: hexagonal and lamellar coexistence" begin
    q  = [0.1009, 0.1048, 0.1748, 0.2017, 0.2095, 0.2668, 0.302, 0.3141, 0.3489, 0.3632]
    sh = [11.1, 12.0, 2.05, 1.17, 1.5, 0.35, 0.195, 0.333, 0.0858, 0.0649]
    cands = collect(indexpeaks(q, sh))
    top = cands[argmax(score.(cands))]
    @test top isa Himalaya.Index{Himalaya.Ia3d}
    for group in _both(cands, _eff(q, sh))
        @test _phases(group) == ["Hexagonal", "Lamellar"]
    end
end

@testset "auto_group: one spurious sharp peak does not change the call" begin
    order = sortperm(vcat(HEX_Q, 0.0421))
    q  = vcat(HEX_Q, 0.0421)[order]
    sh = vcat(HEX_SH, 10.0)[order]
    cands = collect(indexpeaks(q, sh))
    clean = collect(indexpeaks(HEX_Q, HEX_SH))
    # No candidate contains the spurious peak together with the top index's peaks, so
    # the spurious peak cannot win on its own; the call matches the clean trace.
    for (group, ref) in zip(_both(cands, _eff(q, sh)), _both(clean, _eff(HEX_Q, HEX_SH)))
        @test !(length(group) == 1 && first(group) isa Himalaya.Index{Himalaya.Lamellar})
        @test _phases(group) == _phases(ref) == ["Hexagonal"]
    end
end

@testset "auto_group: colliding peak_ids across kinds" begin
    base = _eff(HEX_Q, HEX_SH; ids = collect(11:17))
    ids = collect(11:17); ids[7] = ids[1]
    kinds = fill(:auto, 7); kinds[7] = :curation
    coll = _eff(HEX_Q, HEX_SH; ids = ids, kinds = kinds)
    cands = collect(indexpeaks(HEX_Q, HEX_SH))
    @test _phases(auto_group(cands, coll)) == _phases(auto_group(cands, base))
    @test first(auto_group(cands, coll)) isa Himalaya.Index{Himalaya.Hexagonal}
end

