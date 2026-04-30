using Himalaya
using SparseArrays
using SQLite, DBInterface, Tables

"""
Tolerance (relative to predicted q) used when snapping observed peaks to
predicted ratio positions during speculative-index construction. Matches
`Himalaya.indexpeaks`'s default `tol` so the snap UI agrees with the
phase-search engine on what counts as "close enough".
"""
const SNAP_TOL = 0.0025

"""
    resolve_phase(name) -> Union{Type{<:Himalaya.Phase}, Nothing}

Resolve a phase name string (e.g. `"Pn3m"`, `"Himalaya.Pn3m"`) to the phase
type. Returns `nothing` if the name isn't a known `Himalaya.Phase` subtype.
"""
function resolve_phase(name::AbstractString)
    bare = last(split(String(name), '.'))
    P = try
        getfield(Himalaya, Symbol(bare))
    catch
        return nothing
    end
    (P isa Type && P <: Himalaya.Phase) ? P : nothing
end

"""
    compute_snap(peak_rows, phase, anchor_q, anchor_ratio; tol = SNAP_TOL)
        -> Vector{NamedTuple}

For each ratio position of `phase`, return the predicted q (basis × ratio)
and the closest non-anchor peak whose relative deviation is within `tol`.
`peak_rows` is any iterable of NamedTuples with `id::Int` and `q::Real`
fields (e.g. the result of `Tables.rowtable` over a `peaks` query).

Returns one row per ratio position (1-indexed) with fields:
- `ratio_position::Int`
- `predicted_q::Float64`
- `suggested_peak_id::Union{Int, Nothing}`
- `suggested_q::Union{Float64, Nothing}`
- `suggested_residual::Union{Float64, Nothing}`
"""
function compute_snap(peak_rows, phase::Type{P}, anchor_q::Real, anchor_ratio::Int;
                     tol::Real = SNAP_TOL) where {P<:Himalaya.Phase}
    ratios = Himalaya.phaseratios(P; normalize = true)
    1 <= anchor_ratio <= length(ratios) || error("anchor_ratio $anchor_ratio out of range for $P (1..$(length(ratios)))")
    basis = Float64(anchor_q) / ratios[anchor_ratio]

    out = NamedTuple[]
    for (rpos, ratio) in enumerate(ratios)
        predicted_q = basis * ratio
        best_id      = nothing
        best_q       = nothing
        best_resid   = nothing
        best_relresid = Inf
        for pr in peak_rows
            q  = Float64(pr.q)
            relresid = abs(q - predicted_q) / predicted_q
            if relresid <= tol && relresid < best_relresid
                best_relresid = relresid
                best_id    = Int(pr.id)
                best_q     = q
                best_resid = abs(q - predicted_q)
            end
        end
        push!(out, (
            ratio_position     = rpos,
            predicted_q        = predicted_q,
            suggested_peak_id  = best_id,
            suggested_q        = best_q,
            suggested_residual = best_resid,
        ))
    end
    out
end

"""
    build_speculative_index(peak_rows, phase, ratio_to_peak_id::Dict{Int,Int})
        -> NamedTuple

Construct the numeric fields a speculative `indices` row needs, given a
ratio-position → peak-id assignment. Reuses `Himalaya.fit` and
`Himalaya.score` so speculative indices stay numerically consistent with
auto indices.

`peak_rows` must contain rows with `id`, `q`, and `sharpness` fields.

Returns `(; basis, score, r_squared, lattice_d, residuals)` where
`residuals` is `Dict{Int, Float64}` mapping ratio_position → |q - ideal|.
"""
function build_speculative_index(peak_rows, phase::Type{P},
                                  ratio_to_peak_id::Dict{Int,Int}) where {P<:Himalaya.Phase}
    isempty(ratio_to_peak_id) && error("at least one ratio assignment required")

    ratios = Himalaya.phaseratios(P; normalize = true)

    by_id = Dict{Int, eltype(peak_rows)}()
    for pr in peak_rows
        by_id[Int(pr.id)] = pr
    end

    rpos_sorted = sort(collect(keys(ratio_to_peak_id)))
    qvals      = Float64[]
    sharpvals  = Float64[]
    for rpos in rpos_sorted
        pid = ratio_to_peak_id[rpos]
        haskey(by_id, pid) || error("peak id $pid not in supplied peak_rows")
        pr = by_id[pid]
        push!(qvals,     Float64(pr.q))
        push!(sharpvals, pr.sharpness === nothing || ismissing(pr.sharpness) ? 0.0 : Float64(pr.sharpness))
    end

    n = length(ratios)
    peaks_sv     = SparseVector{Float64, Int}(n, rpos_sorted, qvals)
    sharpness_sv = SparseVector{Float64, Int}(n, rpos_sorted, sharpvals)

    # Least-squares fit through assigned (ratio, q) pairs (intercept fixed at 0).
    # Mirrors `Himalaya.fit` exactly — extracts (idx, q) from the sparse vector.
    observed_ratios_full = Himalaya.phaseratios(P)  # un-normalized
    observed_ratios_used = observed_ratios_full[rpos_sorted]
    basis_unnorm = observed_ratios_used \ qvals  # 1/d in fit's terms

    idx = Himalaya.Index{P}(basis_unnorm, peaks_sv, sharpness_sv)
    fit_result = Himalaya.fit(idx)
    s          = Himalaya.score(idx)

    residuals = Dict{Int, Float64}()
    ratios_unnorm = observed_ratios_full
    for (rpos, qv) in zip(rpos_sorted, qvals)
        ideal = ratios_unnorm[rpos] * basis_unnorm
        residuals[rpos] = abs(qv - ideal)
    end

    (; basis = basis_unnorm,
       score = s,
       r_squared = fit_result.R²,
       lattice_d = fit_result.d,
       residuals = residuals)
end

"""
    insert_speculative_index!(db, exposure_id, phase, ratio_to_peak_id)
        -> Int (new index id)

Builds the speculative index from supplied peak assignments and inserts
it into `indices` with `kind='speculative'`, plus the per-peak rows in
`index_peaks`. Returns the new id. Does **not** add the index to any
group — the caller is responsible for active-group membership.
"""
function insert_speculative_index!(db::SQLite.DB, exposure_id::Int,
                                    phase::Type{P},
                                    ratio_to_peak_id::Dict{Int,Int}) where {P<:Himalaya.Phase}
    peak_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, q, sharpness FROM peaks WHERE exposure_id = ? AND excluded = 0",
        [exposure_id]))

    built = build_speculative_index(peak_rows, P, ratio_to_peak_id)

    res = DBInterface.execute(db,
        """INSERT INTO indices
             (exposure_id, phase, basis, score, r_squared, lattice_d, status, kind)
           VALUES (?, ?, ?, ?, ?, ?, 'candidate', 'speculative')""",
        [exposure_id, string(nameof(P)), built.basis,
         built.score, built.r_squared, built.lattice_d])
    new_id = Int(DBInterface.lastrowid(res))

    for (rpos, peak_id) in ratio_to_peak_id
        DBInterface.execute(db,
            """INSERT INTO index_peaks (index_id, peak_id, ratio_position, residual)
               VALUES (?, ?, ?, ?)""",
            [new_id, peak_id, rpos, built.residuals[rpos]])
    end
    new_id
end
