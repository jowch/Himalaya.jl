using Test
using Himalaya: indexpeaks, Index, peaks, score
using HimalayaUI: auto_group

@testset "auto_group" begin
    qs    = [0.1000, 0.1414, 0.2000]
    proms = [1.0, 0.8, 0.6]

    candidates = indexpeaks(qs, proms)

    group = auto_group(candidates)

    if !isempty(candidates)
        @test !isempty(group)
        peak_sets = [Set(peaks(idx)) for idx in group]
        for i in eachindex(peak_sets), j in eachindex(peak_sets)
            i == j && continue
            @test isempty(intersect(peak_sets[i], peak_sets[j]))
        end
    end
end
