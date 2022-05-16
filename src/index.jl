"""
    Index{::Phase}

Represents an index assignment. Stores the `basis` and the observed `peaks`.
"""
struct Index{P<:Phase}
    basis::Real
    peaks::SparseVector{<:Real}
end

function show(io::IO, index::Index{P}) where P
    idx, xs = findnz(index.peaks)
    peak_str = fill("⋅", length(index.peaks))
    peak_str[idx] = string.(round.(xs; digits = 5))

    print(io, "Index(::$P, $(round(basis(index); digits = 5)), [$(join(peak_str, ' ')...)]")
end

# getters
phase(::Index{P}) where P = P
basis(index::Index) = index.basis
peaks(index::Index) = nonzeros(index.peaks)
numpeaks(index::Index) = nnz(index.peaks)

"""
    predictpeaks(index)

Predicts the expected peaks for the given `index`.
"""
predictpeaks(index::Index{P}) where P = basis(index) * phaseratios(P)

# operations
==(a::Index{P}, b::Index{Q}) where {P,Q} = P <: Q && basis(a) == basis(b)
issubset(a::Index, b::Index) = basis(a) == basis(b) && issubset(peaks(a), peaks(b))

"""
    indexpeaks(peaks, [domain]; tol = 0.005)

Compute valid indexings for a collection of `peaks` by considering each phase
in `PHASES` and all bases in `domain`. If domain is not provided, then `peaks`
will be used as the domain (equivalent to requiring an observed basis).

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

If `gaps`, then allow gaps between observed peaks.

See also `Phase`, `minpeaks`.
"""
function indexpeaks(domain, peaks; tol = 0.005, gaps = true)
    indices = remove_subsets([
        index for phase in (Lamellar, Hexagonal, Pn3m, Im3m, Ia3d, Fd3m)
              for index in indexpeaks(phase, domain, peaks, tol)
    ])

    # apply filter
    if !gaps
        filter!(indices) do index
            peak_idx, _ = findnz(index.peaks)
            !any(==(0), view(index.peaks, first(peak_idx):last(peak_idx)))
        end
    end

    indices
end

indexpeaks(peaks; kwargs...) = indexpeaks(peaks, peaks; kwargs...)

function indexpeaks(::Type{P}, domain, peaks, tol) where {P<:Phase}
    indices = Index[]
    ratios = phaseratios(P)
    observed_ratios = let
        # consider all elements in `domain` as a potential basis value
        B = reshape(domain, length(domain), 1)
        X = reshape(peaks, 1, length(peaks))

        # compute the ratios of observed peaks given each candidate basis
        1 ./ B * X
    end
    candidate_mask = 1 .<= observed_ratios .<= maximum(ratios)

    if !any(candidate_mask)
        return indices
    end

    # match the observed ratios to phase ratios by computing residuals and
    # finding pairs within tolerance
    matched_ratios = let
        comparable_ratios = view(observed_ratios, candidate_mask)
        residuals = zeros(length(comparable_ratios), length(ratios))

        for (i, ratio) in enumerate(ratios)
            @inbounds residuals[:, i] = abs.(comparable_ratios .- ratio)
        end

        # find argminima residuals
        argmin_index = vec(argmin(residuals; dims = 2))
        within_tol = residuals[argmin_index] .< tol

        # update mask with positions that are within tolerance
        candidate_mask[candidate_mask, :] .&= within_tol

        # pull out ratio indices that are within tolerance
        matched_ratios = spzeros(UInt8, size(candidate_mask)...)
        matched_ratios[candidate_mask] .= getindex.(argmin_index[within_tol], 2)

        matched_ratios
    end

    # only keep bases with the minimum number of candidates
    num_candidates = count(candidate_mask, dims = 2)
    candidate_mask[vec(num_candidates .< minpeaks(P)), :] .= 0

    if !any(candidate_mask)
        return indices
    end

    candidate_bases = getindex.(findall(any(candidate_mask; dims = 2)), 1)

    for i in candidate_bases
        peak_idx, ratio_idx = findnz(matched_ratios[i, :])

        push!(
            indices,
            Index{P}(
                domain[i],
                SparseVector{eltype(peaks),UInt8}(
                    length(ratios), ratio_idx, peaks[peak_idx]
                )
            )
        )
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
    fit(index)

Fit and return the lattice constant for a given `index` along with the
associated R² for the fit.
"""
function fit(index::Index{P}) where P
    observed_idx, observed_peaks = findnz(index.peaks)
    observed_ratios = phaseratios(P)[observed_idx]

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
    missing_first = ismissing(first(index.peaks))
    has_gaps = any(ismissing, index.peaks[
        findfirst(!ismissing, index.peaks):findlast(!ismissing, index.peaks)
    ])

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
    subsets = [a != b && issubset(a, b) for a = indices, b = indices]
    indices[.!any(subsets; dims = 2) |> vec]
end
