using Statistics


phases = Dict(
	"Ia3d" 		=> [√6, √8, √14, √16, √20, √22, √24, √26],
	"Pn3m" 		=> [√2, √3, √4, √6, √8, √9, √10, √11],
	"Im3m" 		=> [√2, √4, √6, √8, √10, √12, √14, √16, √18],
	"lamellar" 	=> [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
	"hexagonal" => [1, √3, √4, √7, √9, √11, √12],
	"Fd3m" 		=> [√3, √8, √11, √12, √16, √19, √24, √27, √32, √35, √36]
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

function index_peaks(peaks, n_offset=0; tol=0.005, all=false)
    indices = []

    for phase in keys(phases)
        for offset in 0:n_offset
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

    if !all
        filter!(idx -> count(!ismissing, idx.peaks) >= minpeaks[idx.phase], indices)
    end

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
    candidates = peaks[isapprox(expected_peak; atol=tol).(peaks)]

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

"""

- more peaks is better than fewer
- fewer missing peaks is better
- less error is better
- smaller basis is better
"""
function score(index::Index)
    observed_idx = findall(!ismissing, index.peaks)

    observed_peaks = index.peaks[observed_idx]
    expected_peaks = index.predicted_peaks[observed_idx]

    first_observed = findfirst(!ismissing, index.peaks)
    last_observed = findlast(!ismissing, index.peaks)
    connectedness = (count(!ismissing, index.peaks[first_observed:last_observed])
                     / (last_observed - first_observed + 1))


    # completeness = length(observed_idx) / num_obs_before_trailing
    rmse = sqrt(mean((observed_peaks .- expected_peaks).^2))

    length(observed_idx) * connectedness - rmse - index.basis
end

function fit(index::Index)
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

    if index.phase == "hexagonal"
        λ = 2 / √3
    else
        λ = 1
    end
    
    # a = 2pi / m
    2π * λ / m, R²
end
