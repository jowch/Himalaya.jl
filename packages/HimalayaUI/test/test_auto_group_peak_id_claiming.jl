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
    isempty(candidates) && return

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
    isempty(candidates) && return
    eff = (q = qs, sharpness = sharpness,
           peak_id = collect(1:length(qs)),
           peak_kind = fill(:auto, length(qs)))
    group = auto_group(candidates, eff)
    scores = score.(group)
    @test issorted(scores; rev = true)
end
