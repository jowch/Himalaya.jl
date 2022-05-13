using Statistics

"""
    Index{::Phase}

Represents an index assignment. Stores the `basis` and the observed `peaks`.
"""
struct Index{P<:Phase}
    basis::Real
    peaks::Vector{<:Real}
    observed::BitVector
end

# getters
phase(::Index{P}) where P = P
basis(index::Index) = index.basis
peaks(index::Index) = index.peaks
numpeaks(index::Index) = length(index.peaks)

"""
    predictpeaks(index)

Predicts the expected peaks for the given `index`.
"""
function predictpeaks(index::Index{P}) where P
    phaseratios(P) * basis(index)
end

"""
    missingpeaks(index)

Computes the unobserved, yet expected peaks for the given `index`.
"""
function missingpeaks(index::Index)
    predictpeaks(index)[.!index.observed]
end

# operations
==(a::Index{P}, b::Index{Q}) where {P,Q} = P <: Q && basis(a) == basis(b)
issubset(a::Index, b::Index) = basis(a) == basis(b) && issubset(peaks(a), peaks(b))

"""
    indexpeaks(domain, peaks; tol = 0.005)

Compute valid indexings for a collection of `peaks` by considering each phase
in `PHASES` and all bases in `domain`.

A `Phase`'s peaks positions are determined by its basis peak position and its
ratios. Given a `Phase` and `domain` of possible basis values, we can evaluate
each phase and basis to see whether the expected peak positions are present in
the provided observed `peaks`.

More specifically, we can divide each observed peak value by a given basis to
find their hypothetical ratios and then compare these ratios to the defining
ratios for a `Phase`. Only the "correct" choices of phase and basis will match
the hypothetical ratios. We can set a tolerance `tol` to be the maximum
acceptable deviation of hypothetical ratios from the expected ratios. The peaks
that result in residuals within tolerance are `candidates` for an `Index` of
this phase and basis.

Furthermore, each `Phase` has a minimum number of candidates to be considered
reasonable.

See also `Phase`, `minpeaks`.
"""
function indexpeaks(domain, peaks; tol = 0.005)
    [index for phase in (Lamellar, Hexagonal, Pn3m, Im3m, Ia3d, Fd3m)
           for index in indexpeaks(phase, domain, peaks, tol)]
end

function indexpeaks(::Type{P}, domain, peaks, tol) where P
    observed_ratios = let
        # consider all elements in `domain` as a potential basis value
        B = reshape(domain, length(domain), 1)
        X = reshape(peaks, 1, length(peaks))

        # compute the ratios of observed peaks given each candidate basis
        1 ./ B * X
    end
    ratios = phaseratios(P)
    min_peaks = minpeaks(P)

    # find subset of rows that have hypothetical ratios comparable to phase ratios
    valid_idx, comparable_mask, comparable_ratios = let
        criteria = 1 .<= observed_ratios .<= maximum(ratios)
        valid = vec(count(criteria; dims = 2) .>= min_peaks)

        valid, criteria[valid, :], view(observed_ratios, valid, :)
    end

    if !any(valid_idx)
        return []
    end

    # compute the residuals between each observed ratio and phase ratio
    residuals, matched_ratios = let
        d = zeros(size(comparable_mask)..., length(ratios))

        for (i, ratio) in enumerate(ratios)
            @inbounds d[:,:,i] = abs.(comparable_ratios .- ratio)
        end

        argminima = dropdims(argmin(abs.(d); dims = 3); dims = 3)

        d[argminima, 1], getindex.(argminima, 3) .* comparable_mask 
    end

    # candidate peaks for each basis need to have residuals within tolerance `tol`
    candidate_idx = any(residuals .<= tol; dims = 3) .&& comparable_mask 
    matches = matched_ratios .* candidate_idx

    for i = 1:size(matches, 1)
        ms = view(matches, i, :)
        rs = view(residuals, i, :)
        idxs = unique(ms[ms .!= 0])

        for idx in idxs
            target_idx = ms .== idx
            num_idx = count(target_idx)

            if num_idx .> 1
                suppressed = zeros(Int64, num_idx)
                suppressed[argmin(rs[target_idx])] = idx

                @inbounds matches[i, target_idx] .= suppressed
                @inbounds candidate_idx[i, target_idx] .= suppressed .> 0
            end
        end
    end

    # only keep bases with the minimum number of candidates
    num_candidates = count(candidate_idx, dims = 2)
    valid_bases = vec(num_candidates .>= min_peaks)

    if !any(valid_bases)
        return []
    end

    # compute the total error for each basis; optimal bases will minimize error
    basis_error = sum(abs, residuals[valid_bases, :], dims = 2) |> vec

    [Index{P}(domain[valid_idx][valid_bases][idx],
              peaks[candidate_idx[valid_bases, :][idx, :]],
              BitVector(i ∈ matches[valid_bases, :][idx, :] for i in 1:length(ratios)))
     for idx in sortperm(basis_error)]
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
    fit(index)

Fit and return the lattice constant for a given `index` along with the
associated R² for the fit.
"""
function fit(index::Index{P}) where P
    observed_ratios = phaseratios(P)[index.observed]
    observed_peaks = peaks(index)

    # use least squares fit to compute lattice parameter
    d = (observed_ratios' * observed_ratios) \ (observed_ratios' * observed_peaks)

    if P <: Hexagonal
        λ = 2 / √3
    else
        λ = 1
    end
    
    # a = 2pi / d
    2π * λ / d, R²(d * observed_ratios, observed_peaks)
end

"""
    score(index::Index)

Score the validity of a given `index`. This is simply the number of peaks
indexed times the quality of the index's `fit`.
"""
function score(index::Index)
    _, rsquared = fit(index)
    missing_first = first(index.observed) == 0
    has_gaps = any(index.observed[findfirst(index.observed):findlast(index.observed)] .== 0)

    numpeaks(index) * rsquared - missing_first - has_gaps
end

"""
    remove_duplicates(indices::Vector{Index})

Remove any indices that are strict subsets of another index. This is the case
when two indices share the same basis but the first has a subset of the second's
peaks.

See also `issubset`.
"""
function remove_subsets(indices::Vector{Index})
    subsets = falses(length(indices), length(indices))

    for i = 1:length(indices)
        for j = i:length(indices)
            if i != j
                @inbounds subsets[i, j] = issubset(indices[i], indices[j])
            end
        end
    end

    indices[.!any(subsets; dims = 2) |> vec]
end
