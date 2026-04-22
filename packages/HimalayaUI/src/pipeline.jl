using Himalaya

"""
    auto_group(indices) -> Vector{Index}

Greedily select a non-overlapping set of indices by descending score.
An index is added to the group only if none of its peaks are already
claimed by a previously selected index.
"""
function auto_group(indices::Vector{<:Himalaya.Index})::Vector{<:Himalaya.Index}
    isempty(indices) && return indices
    sorted  = sort(indices; by = Himalaya.score, rev = true)
    claimed = Set{Float64}()
    group   = eltype(indices)[]

    for idx in sorted
        idx_peaks = Set(Himalaya.peaks(idx))
        isempty(intersect(idx_peaks, claimed)) || continue
        push!(group, idx)
        union!(claimed, idx_peaks)
    end
    group
end
