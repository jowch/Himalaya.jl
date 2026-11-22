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

@testset "auto_group is score-ordered descending" begin
    qs        = [0.1, 0.1414, 0.1732, 0.2, 0.2449]
    sharpness = ones(length(qs))
    candidates = indexpeaks(qs, sharpness)
    @test !isempty(candidates)
    eff = (q = qs, sharpness = sharpness,
           peak_id = collect(1:length(qs)),
           peak_kind = fill(:auto, length(qs)))
    group = auto_group(candidates, eff)
    scores = score.(group)
    @test issorted(scores; rev = true)
end

@testset "auto_group first pick contains the sharpest peak" begin
    q  = 0.09789 .* [1, sqrt(3), 2, sqrt(7), 3, sqrt(12), sqrt(13)]
    sh = [4.57, 0.747, 0.383, 0.111, 0.0543, 0.0205, 0.0397]
    candidates = collect(indexpeaks(q, sh))
    @test !isempty(candidates)

    # Precondition: without the rule, an index that skips the sharpest peak wins on score.
    top = candidates[argmax(score.(candidates))]
    @test !(q[1] in peaks(top))

    eff = (q = q, sharpness = sh, peak_id = collect(11:17),
           peak_kind = fill(:auto, length(q)))
    for group in (auto_group(candidates, eff), auto_group(candidates))
        @test !isempty(group)
        @test q[1] in peaks(first(group))
        @test first(group) isa Himalaya.Index{Himalaya.Hexagonal}
        @test !any(idx -> idx isa Himalaya.Index{Himalaya.Ia3d}, group)
    end
end

@testset "auto_group unchanged when no candidate contains the sharpest peak" begin
    q  = 0.09789 .* [1, sqrt(3), 2, sqrt(7), 3, sqrt(12), sqrt(13)]
    sh = [4.57, 0.747, 0.383, 0.111, 0.0543, 0.0205, 0.0397]
    cands = [i for i in indexpeaks(q, sh) if !(q[1] in peaks(i))]
    @test !isempty(cands)
    eff = (q = q, sharpness = sh, peak_id = collect(1:7), peak_kind = fill(:auto, 7))
    expected = first(sort(cands; by = score, rev = true))
    @test first(auto_group(cands, eff)) === expected
    @test first(auto_group(cands)) === expected
end
