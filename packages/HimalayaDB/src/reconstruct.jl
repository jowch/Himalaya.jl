# Local copy of HimalayaUI's phase-name→type resolver (speculative.jl:19).
# Tiny and stable; reimplemented rather than depending on a HimalayaUI internal.
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
    reconstruct_index(db, index_id) -> Himalaya.Index{P}

Rebuild a core `Himalaya.Index{P}` from an `indices` row and its `index_peaks`
supporting peaks. Peaks/sharpness are sparse vectors indexed by ratio position.
Throws `ArgumentError` if the index is missing or its phase name is unknown.

**Fidelity note:** a supporting peak backed by an add-curation (`peak_kind =
'curation'`) reconstructs with `sharpness = 0`, because `peak_curations` stores
no sharpness column — the pipeline samples sharpness from the trace only to
compute the stored `indices.score`, and does not persist that sample. As a
result, `Himalaya.score` of a reconstructed index that includes an add-peak can
diverge from the authoritative stored score (observed: ~0.434 reconstructed vs
~0.469 stored, for a Pn3m index with one add-peak among its supporting peaks).
`fit`, `predictpeaks`, `missingpeaks`, `==`, and `issubset` are unaffected by
this gap — only sharpness-weighted scoring is. When you need the authoritative
score, read it from `indices.score` (via `index_candidates`) rather than
recomputing `Himalaya.score` on a reconstructed index.
"""
function reconstruct_index(db::SQLite.DB, index_id::Integer)
    irows = Tables.rowtable(DBInterface.execute(db,
        "SELECT phase, basis FROM indices WHERE id = ?", [Int(index_id)]))
    isempty(irows) && throw(ArgumentError("reconstruct_index: no index $index_id"))
    irow = irows[1]

    P = resolve_phase(irow.phase)
    P === nothing && throw(ArgumentError(
        "reconstruct_index: unknown phase $(irow.phase) for index $index_id"))

    prows = Tables.rowtable(DBInterface.execute(db, """
        SELECT ip.ratio_position,
               COALESCE(ap.q, pc.q) AS q_observed,
               ap.sharpness AS sharpness
        FROM index_peaks ip
        LEFT JOIN auto_peaks ap     ON ap.id = ip.peak_id AND ip.peak_kind = 'auto'
        LEFT JOIN peak_curations pc ON pc.id = ip.peak_id AND ip.peak_kind = 'curation'
        WHERE ip.index_id = ?
        ORDER BY ip.ratio_position
    """, [Int(index_id)]))

    n = length(Himalaya.phaseratios(P))
    peaks = spzeros(Float64, n)
    sharp = spzeros(Float64, n)
    for r in prows
        pos = Int(r.ratio_position)
        # Persisted positions can outlive a phase's ratio series — src/phase.jl
        # trims and extends series tails (#304 removed Hexagonal's √11, taking
        # it 14 → 13). Skip out-of-range rows rather than BoundsError-ing this
        # read-only API on a DB whose migration hasn't run yet. Same guard as
        # the intent re-resolve loop in HimalayaUI's pipeline.jl.
        1 <= pos <= n || continue
        peaks[pos] = Float64(r.q_observed)
        sharp[pos] = r.sharpness === missing ? 0.0 : Float64(r.sharpness)
    end

    return Himalaya.Index{P}(Float64(irow.basis), peaks, sharp)
end
