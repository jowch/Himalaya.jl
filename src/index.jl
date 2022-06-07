"""
    Index{::Phase}

Represents an index assignment. Stores the `basis` and the observed `peaks`.
"""
struct Index{P<:Phase}
    basis::Real
    peaks::SparseVector{<:Real}
    prom::Real
end

function show(io::IO, index::Index{P}) where P
    idx, xs = findnz(index.peaks)
    peak_str = fill("⋅", length(index.peaks))
    peak_str[idx] = string.(round.(xs; digits = 5))

    print(io, "Index(::$P, $(round(basis(index); digits = 5)), [$(join(peak_str, ' ')...)])")
end

# getters
phase(::Index{P}) where P = P
basis(index::Index) = index.basis
peaks(index::Index) = nonzeros(index.peaks)
numpeaks(index::Index) = nnz(index.peaks)

"""
    totalprom(index)

The total prominence of an index is the sum of the prominences of its peaks.
"""
totalprom(index::Index) = index.prom

"""
    predictpeaks(index)

Predicts the expected peaks for the given `index`.
"""
predictpeaks(index::Index{P}) where P = basis(index) * phaseratios(P; normalize = true)

"""
    missingpeaks(index)

Computes the predicted, missing peaks for a given `index`.
"""
function missingpeaks(index::Index{P}) where P
    peak_idx, _ = findnz(index.peaks)
    missing_idx = setdiff(1:length(phaseratios(P)), peak_idx)
    predictpeaks(index)[missing_idx]
end

# operations
"""
    ==(a::Index, b::Index)

An index `a` is equals (`==`) to an index `b` iff they have the same phase and
the same basis.
"""
==(a::Index{P}, b::Index{Q}) where {P,Q} = P <: Q && basis(a) == basis(b)

"""
    issubset(a::Index, b::Index)

An index `a` is a subset of index `b` iff they have the same basis and the peaks
of `a` are a subset of the peaks of `b`.
"""
issubset(a::Index, b::Index) = basis(a) == basis(b) && issubset(peaks(a), peaks(b))

"""
    indexpeaks([phase], peaks, [proms], [domain];
               tol = 0.0025, gaps = true, requiremin = true)

Compute possible index assingments for a collection of `peaks` by considering
a specific `phase` or all available phases and all bases in `domain`. If a
`domain` is not provided, then `peaks` will be used as the domain, which is
equivalent to requiring bases to be observed peaks.

A `Phase`'s peaks positions are determined by its basis peak position and its
ratios. This is computed by multiplying the phase's normalized ratios (first
ratio is 1) by the basis. Thus, given a `Phase` and `domain` of possible basis
values, we can evaluate each phase and basis to see whether the expected peak
positions are present in the provided observed `peaks`.

More specifically, we can divide each observed peak value by a given basis to
find their hypothetical ratios and then compare these ratios to the defining
ratios for a `Phase`. Only the "correct" choices of phase and basis will match
the hypothetical ratios. 

Tolerance `tol` is the maximum acceptable deviation of candidate peaks from
observed peaks. The peaks that result in residuals (deviation) within this
tolerance are candidates for an `Index` of this phase and basis.

If `gaps`, peaks are allowed to be missing between the first and last observed
peaks.

If `requiremin`, require the minimum number of peaks for each phase as defined
by `minpeaks`.

See also `Phase`, `minpeaks`.
"""
function indexpeaks(peaks, proms, domain; kwargs...)
    indices = Index[]

    for phase in (Lamellar, Hexagonal, Pn3m, Im3m, Ia3d, Fd3m)
        push!(indices, indexpeaks(phase, peaks, proms, domain; kwargs...)...)
    end

    remove_subsets(indices)
end

function indexpeaks(::Type{P}, peaks, proms, domain; gaps = true, tol = 0.0025) where {P<:Phase}
    indices = indexpeaks(P, peaks, proms, domain, tol)

    # apply filter
    if !gaps
        filter!(indices) do index
            peak_idx, _ = findnz(index.peaks)
            !any(==(0), view(index.peaks, first(peak_idx):last(peak_idx)))
        end
    end

    indices
end

indexpeaks(peaks; kwargs...) = indexpeaks(peaks, ones(length(peaks)), peaks; kwargs...)
indexpeaks(peaks, proms; kwargs...) = indexpeaks(peaks, proms, peaks; kwargs...)
indexpeaks(phase::Type{<:Phase}, peaks; kwargs...) = indexpeaks(phase, peaks, ones(length(peaks)), peaks; kwargs...)
indexpeaks(phase::Type{<:Phase}, peaks, proms; kwargs...) = indexpeaks(phase, peaks, proms, peaks; kwargs...)

function indexpeaks(::Type{P}, peaks, proms, domain, tol) where {P<:Phase}
    indices = Index[]
    ratios = phaseratios(P; normalize = true)

    # compute the ratios of observed peaks given each candidate basis
    # also compute the adjusted tolerable peak error as ρ = (x + ε) / b
    observed_ratios, tols = let
        # consider all elements in `domain` as a potential basis value
        B = reshape(domain, 1, length(domain))
        X = reshape(peaks, length(peaks), 1)

        observed_ratios = X * (1 ./ B)

        # compute mask of valid ratios
        # TODO: work out addition tol on max side
        mask = 1 - tol .<= observed_ratios .<= maximum(ratios)
        
        # these are (# peaks X # domain)
        sparse(observed_ratios .* mask), tol ./ domain
    end

    if nnz(observed_ratios) == 0
        return indices
    end

    assignments = spzeros(UInt8, size(observed_ratios))

    for (i, θ) in enumerate(tols)
        # for each basis
        residuals = spzeros(length(ratios), length(peaks))
        nz, vals = findnz(observed_ratios[:, i])

        for (j, ratio) in enumerate(ratios)
            for (k, δ) in zip(nz, abs.(vals .- ratio))
                if δ < θ
                    @inbounds residuals[j, k] = δ + eps()
                end
            end
        end

        rs, ps, δs = findnz(residuals)
        used_ratios = Set()
        assign = view(assignments, :, i)

        # go through matches in order of increasing error
        # only do assignment if ratio and peak haven't already been assigned
        for j = sortperm(δs)
            r, p = rs[j], ps[j]
            if r ∉ used_ratios && assign[p] == 0
                @inbounds assign[p] = r
                push!(used_ratios, r)
            end
        end
    end

    # only keep bases with the minimum number of candidates
    num_candidates = count(>(0), assignments; dims = 1)
    candidate_idx = vec(num_candidates .>= minpeaks(P))
    # zero-out non-candidates
    assignments[:, .!candidate_idx] .= 0
    dropzeros!(assignments)

    if nnz(observed_ratios) == 0
        return indices
    end

    for i in findall(candidate_idx)
        peak_idx, ratio_idx = findnz(assignments[:, i])

        push!(
            indices,
            Index{P}(
                domain[i],
                SparseVector{eltype(peaks),UInt8}(
                    length(ratios), ratio_idx, peaks[peak_idx]
                ),
                sum(proms[peak_idx])
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

Score the validity of a given `index`. This score attempts to quantify the
likeliness of a given `index`, where higher values are better. It captures the
number of peaks, the prominence of those peaks, the fit of the index, and the
number of missing peaks.

See also `fit` and `totalprom`.
"""
function score(index::Index)
    if numpeaks(index) > 1
        _, rsquared = fit(index)
    else
        rsquared = 1
    end

    num_gaps = let
        peak_idx, _ = findnz(index.peaks)
        count(==(0), view(index.peaks, first(peak_idx):last(peak_idx)))
    end
    (numpeaks(index) + totalprom(index)) * (1 - 0.25 * num_gaps) * rsquared
end

"""
    remove_duplicates(indices::Vector{Index})

Remove any indices that are strict subsets of another index. This is the case
when two indices share the same basis but the first has a subset of the second's
peaks.

See also `issubset` and `score`.
"""
function remove_subsets(indices::Vector{<:Index})
    subsets = [
        a != b && issubset(a, b) && score(a) < score(b)
        for a = indices, b = indices
    ]
    indices[.!any(subsets; dims = 2) |> vec]
end
