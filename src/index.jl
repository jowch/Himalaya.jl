using Statistics


PHASES = Dict(
	:Ia3d		=> [√6, √8, √14, √16, √20, √22, √24, √26],
	:Pn3m		=> [√2, √3, √4, √6, √8, √9, √10, √11],
	:Im3m		=> [√2, √4, √6, √8, √10, √12, √14, √16, √18],
	:Lamellar	=> [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
	:Hexagonal  => [1, √3, √4, √7, √9, √11, √12],
	:Fd3m		=> [√3, √8, √11, √12, √16, √19, √24, √27, √32, √35, √36]
)

minpeaks = Dict(
	:Ia3d		=> 3,
	:Pn3m		=> 3,
	:Im3m		=> 3,
	:Lamellar	=> 2,
	:Hexagonal  => 2,
	:Fd3m		=> 2
)

struct Index{T}
    phase::Symbol
    basis::T
    peaks::Vector{Union{Missing, T}}

    Index(phase, basis, peaks) = begin
        if phase ∉ keys(PHASES)
            throw(ArgumentError("Invalid phase $phase"))
        end
        new{typeof(basis)}(Symbol(phase), basis, peaks)
    end
end

function Index(phase, basis, peaks)
    Index(
        phase,
        basis,
        convert(Array{Union{Missing, typeof(basis)}}, peaks),
    )
end

function ==(a::Index, b::Index)
    a.phase == b.phase && a.basis == b.basis
end

function npeak(index::Index)
    observed_idx = findall(!ismissing, index.peaks)
    length(index.peaks[observed_idx])
end

function predictpeaks(index::Index)
    expected_ratios = PHASES[index.phase]
    expected_ratios .* index.basis
end

function indexpeaks(peaks; tol=0.005, all=false)
    indices = []

    for (i, basis) in enumerate(peaks)
        for phase in keys(PHASES)
            index = Index(phase, basis, [basis])
            indexpeaks!(index, peaks[i+1:end], tol=tol)
            push!(indices, index)
        end
    end

    if !all
        filter!(indices) do index
            all(!ismissing.(index.peaks[1:minpeaks[idx.phase]]))
        end
    end

    indices
end

"""
At any given step of the search process, there are two options, either the next
peak is missing, or that 
"""
function indexpeaks!(index::Index{T}, peaks; tol=0.005) where T
    if isempty(peaks)
        return index
    end

    predicted_peaks = predictpeaks(index)

    if length(index.peaks) >= length(predicted_peaks)
        return index
    end

    expected_peak = predicted_peaks[length(index.peaks) + 1]
    candidates = peaks[isapprox(expected_peak; atol=tol).(peaks)]

    if isempty(candidates)
        remaining_peaks = peaks[peaks .> expected_peak]

        push!(index.peaks, missing)
    else
        selected_peak = candidates[argmin(candidates .- expected_peak)]
        remaining_peaks = peaks[(peaks .> expected_peak) .& (peaks .> last(candidates))]

        push!(index.peaks, selected_peak)
    end

    indexpeaks!(index, remaining_peaks)
end


function score(index::Index)
    # find indices of non-missing peaks
    nonmissing_idx = findall(!ismissing, index.peaks)
    nonmissing_peaks = index.peaks[nonmissing_idx]
    expected_ratios = PHASES[index.phase][nonmissing_idx]

    # use least squares fit to compute lattice parameter
    m = (expected_ratios' * expected_ratios) \ (expected_ratios' * nonmissing_peaks)

    # compute R² value for fit
    RSS = sum((nonmissing_peaks .- (m * expected_ratios)).^2)
    TSS = sum((nonmissing_peaks .- mean(nonmissing_peaks)).^2)

    R² = 1 - (RSS/TSS)

    R² * npeak(index)
end

function fit(index::Index)
    # find indices of non-missing peaks
    nonmissing_idx = findall(!ismissing, index.peaks)
    nonmissing_peaks = index.peaks[nonmissing_idx]
    expected_ratios = PHASES[index.phase][nonmissing_idx]

    # use least squares fit to compute lattice parameter
    m = (expected_ratios' * expected_ratios) \ (expected_ratios' * nonmissing_peaks)

    # compute R² value for fit
    RSS = sum((nonmissing_peaks .- (m * expected_ratios)).^2)
    TSS = sum((nonmissing_peaks .- mean(nonmissing_peaks)).^2)

    R² = 1 - (RSS/TSS)

    if index.phase == "hexagonal"
        λ = 2 / √3
    else
        λ = 1
    end
    
    # a = 2pi / m
    2π * λ / m, R²
end
