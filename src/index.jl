using Statistics


const PHASES = Dict(
	:Ia3d		=> [√6, √8, √14, √16, √20, √22, √24, √26],
	:Pn3m		=> [√2, √3, √4, √6, √8, √9, √10, √11],
	:Im3m		=> [√2, √4, √6, √8, √10, √12, √14, √16, √18],
	:Lamellar	=> [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
	:Hexagonal  => [1, √3, √4, √7, √9, √11, √12],
	:Fd3m		=> [√3, √8, √11, √12, √16, √19, √24, √27, √32, √35, √36]
)

const MINPEAKS = Dict(
	:Ia3d		=> 3,
	:Pn3m		=> 4,
	:Im3m		=> 3,
	:Lamellar	=> 2, 
	:Hexagonal  => 3,
	:Fd3m		=> 3 
)

struct Index{T}
    phase::Symbol
    basis::T
    peaks::Vector{Union{Missing, T}}

    Index(phase, basis, peaks) = begin
        if phase ∉ keys(PHASES)
            throw(ArgumentError("Invalid phase $phase"))
        end
        
        new{typeof(basis)}(
            Symbol(phase),
            basis,
            convert(Array{Union{Missing, typeof(basis)}}, peaks)
        )
    end
end

function ==(a::Index, b::Index)
    a.phase == b.phase && a.basis == b.basis
end

"""
    npeak(index::Index)

Returns the number of non-missing peaks for a given `index`
"""
function npeak(index::Index)
    count(!ismissing, index.peaks)
end

"""
    predictpeaks(phase, basis)

Predicts the expected peaks for a given `phase` and `basis` peak
"""
function predictpeaks(phase::Symbol, basis::T) where T
    expected_ratios = PHASES[phase]
    expected_ratios ./ first(expected_ratios) .* basis
end

"""
    predictpeaks(index::Index)

Predicts the expected peaks for a given `index` using its phase and basis
"""
function predictpeaks(index::Index)
    predictpeaks(index.phase, index.basis)
end

"""
    indexpeaks(peaks)

Compute valid indexings for a collection of `peaks` by considering each phase
in `PHASES` and each choice of basis peak.

Given a collection of peaks, we can select each peak as a basis for an
indexing and then evaluate different phases with this choice of basis. This is
done by computing the expected peak positions based on the phase and basis.
For each expected peak, we use `isapprox` to look for observed `peaks` no more
than `tol` away. If not found, that peak is marked as `missing` in the resulting
indexing.

An indexing is only considered to be valid if it meets the required number of
peaks for its phase.
"""
function indexpeaks(peaks; tol=0.0025)
    indices = []

    for basis in peaks
        for phase in keys(PHASES)
            valid_peaks = view(peaks, peaks .> basis)

            if length(valid_peaks) < MINPEAKS[phase]
                continue
            end

            expected_peaks = filter(predictpeaks(phase, basis)) do expected_peak
                basis < expected_peak <= maximum(valid_peaks)
            end

            if length(expected_peaks) < MINPEAKS[phase]
                continue
            end

            index_peaks = Array{Union{Missing, typeof(basis)}}(missing, 1 + length(expected_peaks))
            index_peaks[1] = basis

            for (j, expected_peak) in enumerate(expected_peaks)
                candidates = view(
                    valid_peaks,
                    isapprox(expected_peak; atol = tol).(valid_peaks)
                )

                if isempty(candidates)
                    valid_peaks = view(peaks, peaks .> expected_peak)
                else
                    index_peaks[j + 1] = candidates[argmin(candidates .- expected_peak)]                  
                    valid_peaks = view(peaks, peaks .> last(candidates))
                end
            end

            # ensure that none of the first `minpeaks` number of peaks are missing.
            # ie. require that the first `number_required` are present
            if all(.!ismissing.(index_peaks[1:MINPEAKS[phase]]))
                push!(indices, Index(phase, basis, index_peaks))
            end
        end
    end

    indices
end

"""
    R²(ŷ, y)

Compute the coefficient of determination (R²) given some predicted values `ŷ`
and the target values `y`.
"""
function R²(ŷ, y)
    RSS = sum((y .- ŷ).^2)
    TSS = sum((y .- mean(y)).^2)

    1 - (RSS / TSS)
end

"""
    fit(index::Index)

Fit and return the lattice constant for a given `index` along with the
associated R² for the fit.
"""
function fit(index::Index)
    # find indices of non-missing peaks
    nonmissing_idx = findall(!ismissing, index.peaks)
    nonmissing_peaks = index.peaks[nonmissing_idx]
    expected_ratios = PHASES[index.phase][nonmissing_idx]

    # use least squares fit to compute lattice parameter
    d = (expected_ratios' * expected_ratios) \ (expected_ratios' * nonmissing_peaks)

    if index.phase == "hexagonal"
        λ = 2 / √3
    else
        λ = 1
    end
    
    # a = 2pi / d
    2π * λ / d, R²(d * expected_ratios, nonmissing_peaks)
end

"""
    score(index::Index)

Score the validity of a given `index`. This is simply the percentage of
reasonably expected peaks actually found.
"""
function score(index::Index)
    count(!ismissing, index.peaks) / length(index.peaks)
end
