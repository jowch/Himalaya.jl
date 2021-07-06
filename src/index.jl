using Statistics


phases = Dict(
	"Ia3d" 		=> [√6, √8, √14, √16, √20, √22, √24, √26],
	"Pn3m" 		=> [√2, √3, √4, √6, √8, √9, √10, √11],
	"Im3m" 		=> [√2, √4, √6, √8, √10, √12, √14, √16, √18],
	"lamellar" 	=> [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
	"hexagonal" => [1, √3, √4, √7, √9, √11, √12],
	"Fd3m" 		=> [√3, √8, √11, √12, √16]
)

minpeaks = Dict(
	"Ia3d" 		=> 3,
	"Pn3m" 		=> 3,
	"Im3m" 		=> 3,
	"lamellar" 	=> 2,
	"hexagonal" => 2,
	"Fd3m" 		=> 2
)

struct Index{T}
    phase::String
    basis::T
    peaks::Vector{Union{Missing, T}}
    predicted_peaks::Vector{T}
end

function Index(phase, basis, peaks)
    basis_position = findfirst(!ismissing, peaks)
    expected_ratios = phases[phase]
    predicted_peaks = expected_ratios ./ expected_ratios[basis_position] .* basis

    Index(
        phase,
        basis,
        convert(Array{Union{Missing, typeof(basis)}}, peaks),
        predicted_peaks
    )
end

function index_peaks(peaks; tol=0.01)
    indices = []

    for phase in keys(phases)
        for offset in 0:(length(phases[phase])-2)
            for i in 1:length(peaks)
                index = Index(
                    phase,
                    peaks[i],
                    [Array{Missing}(missing, offset)..., peaks[i]],
                )

                index_peaks!(index, peaks[peaks .> peaks[i]], tol=tol)
 
                push!(indices, index)
            end
        end
    end

    filter!(idx -> count(!ismissing, idx.peaks) >= minpeaks[idx.phase], indices)

    subsets = any([check_subset(i, j) for i in indices, j in indices]; dims=2)[:]
    deleteat!(indices, subsets)

    indices[sortperm(score.(indices), rev=true)]
end

"""
At any given step of the search process, there are two options, either the next
peak is missing, or that 
"""
function index_peaks!(index::Index{T}, peaks; tol=0.005) where T
    if isempty(peaks)
        return index
    end

    if length(index.peaks) >= length(index.predicted_peaks)
        return index
    end

    expected_peak = index.predicted_peaks[length(index.peaks) + 1]
    candidates = peaks[isapprox(expected_peak; rtol=tol).(peaks)]

    if isempty(candidates)
        remaining_peaks = peaks[peaks .> expected_peak]

        push!(index.peaks, missing)
    else
        selected_peak = candidates[argmin(candidates .- expected_peak)]
        remaining_peaks = peaks[(peaks .> expected_peak) .& (peaks .> last(candidates))]

        push!(index.peaks, selected_peak)
    end

    index_peaks!(index, remaining_peaks)
end

function check_subset(i::Index, j::Index)
    x = i.peaks[findfirst(!ismissing, i.peaks):end]
    !isequal(x, j.peaks) && issubset(replace.([x, j.peaks], missing => 0)...)
end

function score(index::Index)
    observed_idx = findall(!ismissing, index.peaks)

    observed_peaks = index.peaks[observed_idx]
    expected_peaks = index.predicted_peaks[observed_idx]

    last_observed = findlast(!ismissing, index.peaks)
    num_obs_before_trailing = length(index.predicted_peaks[1:last_observed])

    completeness = length(observed_idx) / num_obs_before_trailing
    rmse = sqrt(mean((observed_peaks .- expected_peaks).^2))

    completeness - rmse
end

function fit_cubic(index::Index)
    # find indices of non-missing peaks
    nonmissing_idx = findall(!ismissing, index.peaks)
    nonmissing_peaks = index.peaks[nonmissing_idx]
    expected_ratios = phases[index.phase][nonmissing_idx]

    # use least squares fit to compute lattice parameter
    m = (expected_ratios' * expected_ratios) \ (expected_ratios' * nonmissing_peaks)

    # compute R² value for fit
    RSS = sum((nonmissing_peaks .- (m * expected_ratios)).^2)
    TSS = sum((nonmissing_peaks .- mean(nonmissing_peaks)).^2)

    R² = 1 - (RSS/TSS)
    
    # a = 2pi / m
    2pi / m, R²
end

function fit(index::Index)
    if (index.phase !== "lamellar")
        return fit_cubic(index)

    else
        nonmissing_idx = findall(!ismissing, index.peaks)
        nonmissing_peaks = index.peaks[nonmissing_idx]
        expected_peaks = index.predicted_peaks[nonmissing_idx]

        deviation = sum(abs.(expected_peaks .- nonmissing_peaks))

        return inv(mean(diff(nonmissing_peaks))), deviation
    end
end
