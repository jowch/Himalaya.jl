using Himalaya
using SparseArrays
using Tables

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

"""
Tolerance (relative to the snapshotted basis) used to re-attach custom-group
members across reanalysis: a freshly inserted candidate of the same phase
counts as the "same" indexing if `|Δbasis| ≤ MEMBER_REATTACH_RELTOL · basis`.
"""
const MEMBER_REATTACH_RELTOL = 0.05

function persist_analysis!(db::SQLite.DB, exposure_id::Int,
                            q_full::Vector{Float64},
                            I_full::Vector{Float64},
                            peaks_result::NamedTuple,
                            candidates::Vector{<:Himalaya.Index},
                            group_indices::Vector{<:Himalaya.Index})

    # Snapshot the custom group's members by *semantic identity* (phase + basis)
    # BEFORE we delete the indices rows. Without this, the deletion below
    # invalidates every PK in `index_group_members`, leaving the active set
    # empty after reanalysis.
    custom_member_identities = Tables.rowtable(DBInterface.execute(db, """
        SELECT g.id AS group_id, i.phase, i.basis
        FROM index_groups g
        JOIN index_group_members m ON m.group_id = g.id
        JOIN indices i ON i.id = m.index_id
        WHERE g.exposure_id = ? AND g.kind = 'custom'
        """, [exposure_id]))

    # Remove prior auto peaks, indices, and auto groups for this exposure.
    # Manual peaks (source='manual') are preserved.
    DBInterface.execute(db,
        "DELETE FROM index_group_members WHERE group_id IN
         (SELECT id FROM index_groups WHERE exposure_id = ? AND kind = 'auto')",
        [exposure_id])
    DBInterface.execute(db,
        "DELETE FROM index_groups WHERE exposure_id = ? AND kind = 'auto'",
        [exposure_id])
    # Clear stale custom-group memberships too — we'll re-attach by semantic
    # identity below. Orphan rows would otherwise sit alongside the new ones
    # and confuse the JOIN in `_group_with_members`.
    DBInterface.execute(db,
        "DELETE FROM index_group_members WHERE group_id IN
         (SELECT id FROM index_groups WHERE exposure_id = ? AND kind = 'custom')",
        [exposure_id])
    DBInterface.execute(db,
        "DELETE FROM index_peaks WHERE index_id IN
         (SELECT id FROM indices WHERE exposure_id = ?)", [exposure_id])
    DBInterface.execute(db,
        "DELETE FROM indices WHERE exposure_id = ?", [exposure_id])

    # Snapshot the q-values of any auto peaks the user explicitly excluded so
    # we can re-apply that override after re-detection. Without this, every
    # reanalysis silently undoes the user's curation work.
    excluded_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT q FROM peaks WHERE exposure_id = ? AND source = 'auto' AND excluded = 1",
        [exposure_id]))
    excluded_qs = Float64[Float64(r.q) for r in excluded_rows]

    DBInterface.execute(db,
        "DELETE FROM peaks WHERE exposure_id = ? AND source = 'auto'", [exposure_id])

    # Persist auto peaks. If the q-value matches (within a small tolerance) one
    # the user previously excluded, carry the `excluded = 1` flag forward.
    EXCLUDE_TOL = 1e-6
    function was_excluded(q::Float64)::Bool
        for eq in excluded_qs
            if abs(eq - q) <= max(EXCLUDE_TOL, abs(q) * 0.001)
                return true
            end
        end
        false
    end

    q_to_peak_id = Dict{Float64, Int}()
    for i in eachindex(peaks_result.q)
        qval       = peaks_result.q[i]
        full_idx   = peaks_result.indices[i]
        intensity  = I_full[full_idx]
        res = DBInterface.execute(db,
            "INSERT INTO peaks (exposure_id, q, intensity, prominence, sharpness, source, excluded)
             VALUES (?, ?, ?, ?, ?, 'auto', ?)",
            [exposure_id, qval, intensity,
             peaks_result.prominence[i], peaks_result.sharpness[i],
             Int(was_excluded(qval))])
        q_to_peak_id[qval] = Int(DBInterface.lastrowid(res))
    end

    # Persist candidate indices
    candidate_to_db_id = Dict{Int, Int}()
    for (ci, idx) in enumerate(candidates)
        P          = Himalaya.phase(idx)
        fit_result = Himalaya.fit(idx)
        s          = Himalaya.score(idx)
        res = DBInterface.execute(db,
            "INSERT INTO indices (exposure_id, phase, basis, score, r_squared, lattice_d, status)
             VALUES (?, ?, ?, ?, ?, ?, 'candidate')",
            [exposure_id, string(nameof(P)), Himalaya.basis(idx),
             s, fit_result.R², fit_result.d])
        db_id = Int(DBInterface.lastrowid(res))
        candidate_to_db_id[ci] = db_id

        # Persist index_peaks: for each assigned (ratio_position, q-value) pair
        ratio_positions, peak_qvals = SparseArrays.findnz(idx.peaks)
        ratios_normed = Himalaya.phaseratios(P; normalize=true)
        for (rpos, qval) in zip(ratio_positions, peak_qvals)
            peak_id = get(q_to_peak_id, qval, nothing)
            peak_id === nothing && continue
            ideal = ratios_normed[rpos] * Himalaya.basis(idx)
            resid = abs(qval - ideal)
            DBInterface.execute(db,
                "INSERT OR IGNORE INTO index_peaks (index_id, peak_id, ratio_position, residual)
                 VALUES (?, ?, ?, ?)",
                [db_id, peak_id, rpos, resid])
        end
    end

    # Persist auto group
    res = DBInterface.execute(db,
        "INSERT INTO index_groups (exposure_id, kind, active) VALUES (?, 'auto', 1)",
        [exposure_id])
    group_db_id = Int(DBInterface.lastrowid(res))

    group_set = Set(group_indices)
    for (ci, idx) in enumerate(candidates)
        idx in group_set || continue
        db_id = candidate_to_db_id[ci]
        DBInterface.execute(db,
            "INSERT INTO index_group_members (group_id, index_id) VALUES (?, ?)",
            [group_db_id, db_id])
    end

    # ── Re-attach custom-group members by semantic identity ────────────────
    # For each (group_id, phase, basis) we snapshotted before the delete,
    # find the freshly-inserted candidate of the same phase whose basis is
    # closest to the snapshot. If the closest match is within
    # MEMBER_REATTACH_RELTOL of the snapshotted basis, attach it; otherwise
    # silently drop the member (the change in peaks invalidated that
    # indexing — that's the honest outcome).
    if !isempty(custom_member_identities)
        new_by_phase = Dict{String, Vector{Tuple{Float64, Int}}}()
        for (ci, idx) in enumerate(candidates)
            phase_name = string(nameof(Himalaya.phase(idx)))
            push!(get!(new_by_phase, phase_name, Tuple{Float64, Int}[]),
                  (Himalaya.basis(idx), candidate_to_db_id[ci]))
        end

        for r in custom_member_identities
            phase = String(r.phase)
            prev_basis = Float64(r.basis)
            cands = get(new_by_phase, phase, nothing)
            cands === nothing && continue
            best_id, best_delta = 0, Inf
            for (b, id) in cands
                d = abs(b - prev_basis)
                if d < best_delta
                    best_delta, best_id = d, id
                end
            end
            tol = max(MEMBER_REATTACH_RELTOL * prev_basis, 1e-9)
            if best_id != 0 && best_delta <= tol
                DBInterface.execute(db,
                    "INSERT OR IGNORE INTO index_group_members (group_id, index_id)
                     VALUES (?, ?)", [Int(r.group_id), best_id])
            end
        end

        # If any custom group survived with at least one re-attached member,
        # it remains the active set — demote the freshly-created auto group.
        # Mirrors the post-curate semantics of `ensure_custom_group!`.
        custom_still_populated = Tables.rowtable(DBInterface.execute(db,
            "SELECT g.id FROM index_groups g
             WHERE g.exposure_id = ? AND g.kind = 'custom'
               AND EXISTS (SELECT 1 FROM index_group_members m WHERE m.group_id = g.id)",
            [exposure_id]))
        if !isempty(custom_still_populated)
            DBInterface.execute(db,
                "UPDATE index_groups SET active = 0 WHERE id = ?", [group_db_id])
        end
    end
end

function get_peaks_for_exposure(db::SQLite.DB, exposure_id::Int)
    Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM peaks WHERE exposure_id = ? ORDER BY q", [exposure_id]))
end

function get_indices_for_exposure(db::SQLite.DB, exposure_id::Int)
    Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM indices WHERE exposure_id = ? ORDER BY score DESC", [exposure_id]))
end

function get_groups_for_exposure(db::SQLite.DB, exposure_id::Int)
    Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM index_groups WHERE exposure_id = ? ORDER BY id", [exposure_id]))
end

"""
    init_experiment!(db; kwargs...) -> experiment_id

Create the experiment row. Thin wrapper over `create_experiment!`.
"""
function init_experiment!(db::SQLite.DB; kwargs...)
    create_experiment!(db; kwargs...)
end

"""
    analyze_exposure!(db, exposure_id, analysis_dir)

Load the .dat file for `exposure_id`, run findpeaks + indexpeaks,
auto-group results, and persist everything to the DB.
The .dat filename is taken from `exposures.filename` with `.dat` appended.
"""
function analyze_exposure!(db::SQLite.DB, exposure_id::Int, analysis_dir::String)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT filename FROM exposures WHERE id = ?", [exposure_id]))
    isempty(rows) && error("exposure $exposure_id not found")
    filename = rows[1].filename
    dat_path = joinpath(analysis_dir, filename * ".dat")
    isfile(dat_path) || error("dat file not found: $dat_path")

    q, I, σ      = load_dat(dat_path)
    peaks_result = Himalaya.findpeaks(q, I, σ)
    candidates   = Himalaya.indexpeaks(peaks_result.q, peaks_result.sharpness)
    group        = auto_group(candidates)

    persist_analysis!(db, exposure_id, q, I, peaks_result, candidates, group)
end
