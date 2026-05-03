using Himalaya
using SparseArrays
using Tables
using JSON3

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

"""
    effective_peaks(db, exposure_id, q_grid, I) -> NamedTuple

Compute the effective peak set for analysis. Returns
`(q::Vector{Float64}, sharpness::Vector{Float64},
  peak_id::Vector{Int}, peak_kind::Vector{Symbol})` sorted by q.
Auto peaks whose q matches an `exclude` curation are dropped;
`add` curations are unioned in with sharpness sampled from the
current trace.
"""
function effective_peaks(db::SQLite.DB, exposure_id::Int,
                          q_grid::Vector{Float64}, I::Vector{Float64})
    auto = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, q, sharpness FROM auto_peaks WHERE exposure_id = ?", [exposure_id]))
    excludes = Tables.rowtable(DBInterface.execute(db,
        "SELECT q FROM peak_curations WHERE exposure_id = ? AND kind = 'exclude'",
        [exposure_id]))
    adds = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, q FROM peak_curations WHERE exposure_id = ? AND kind = 'add'",
        [exposure_id]))

    tol(q) = max(1e-6, abs(q) * 0.001)
    is_excluded(qv) = any(e -> abs(Float64(e.q) - qv) <= tol(qv), excludes)

    qs        = Float64[]
    shs       = Float64[]
    peak_id   = Int[]
    peak_kind = Symbol[]
    for r in auto
        qv = Float64(r.q)
        is_excluded(qv) && continue
        push!(qs, qv)
        push!(shs, ismissing(r.sharpness) ? 0.0 : Float64(r.sharpness))
        push!(peak_id,   Int(r.id))
        push!(peak_kind, :auto)
    end
    sharp_full = isempty(adds) ? Float64[] : Himalaya.sharpness(I)
    for r in adds
        qv = Float64(r.q)
        push!(qs, qv)
        push!(shs, isempty(sharp_full) ? 0.0 : sharp_full[argmin(abs.(q_grid .- qv))])
        push!(peak_id,   Int(r.id))
        push!(peak_kind, :curation)
    end

    perm = sortperm(qs)
    (q = qs[perm], sharpness = shs[perm],
     peak_id = peak_id[perm], peak_kind = peak_kind[perm])
end

"""
Match new findpeaks output against existing auto_peaks rows by q-value within
tolerance, UPDATE matched rows in place, INSERT unmatched, DELETE orphans.
Populates `findpeaks_index` from `peaks_result.indices[i]` on every INSERT/UPDATE.
"""
function diff_update_auto_peaks!(db::SQLite.DB, exposure_id::Int,
                                  peaks_result::NamedTuple,
                                  I_full::Vector{Float64})
    EXCLUDE_TOL = 1e-6
    tol(q) = max(EXCLUDE_TOL, abs(q) * 0.001)

    existing = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, q FROM auto_peaks WHERE exposure_id = ?", [exposure_id]))
    remaining = Set{Int}(Int(r.id) for r in existing)

    q_to_id = Dict{Float64, Int}()
    for i in eachindex(peaks_result.q)
        qval      = peaks_result.q[i]
        full_idx  = peaks_result.indices[i]
        intensity = I_full[full_idx]
        prom      = peaks_result.prominence[i]
        sharp     = peaks_result.sharpness[i]

        # Find closest existing auto peak within tolerance.
        best_id, best_d = 0, Inf
        for r in existing
            Int(r.id) in remaining || continue
            d = abs(Float64(r.q) - qval)
            if d < best_d && d <= tol(qval)
                best_d, best_id = d, Int(r.id)
            end
        end

        if best_id != 0
            DBInterface.execute(db,
                """UPDATE auto_peaks SET q = ?, intensity = ?, prominence = ?,
                                         sharpness = ?, findpeaks_index = ?
                   WHERE id = ?""",
                [qval, intensity, prom, sharp, full_idx, best_id])
            delete!(remaining, best_id)
            q_to_id[qval] = best_id
        else
            res = DBInterface.execute(db,
                """INSERT INTO auto_peaks (exposure_id, q, intensity, prominence,
                                           sharpness, findpeaks_index)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                [exposure_id, qval, intensity, prom, sharp, full_idx])
            q_to_id[qval] = Int(DBInterface.lastrowid(res))
        end
    end

    # Drop auto peaks that no longer correspond to any new detection.
    # Also drop the index_peaks rows that referenced them (kind='auto').
    for orphan_id in remaining
        DBInterface.execute(db,
            "DELETE FROM index_peaks WHERE peak_id = ? AND peak_kind = 'auto'", [orphan_id])
        DBInterface.execute(db, "DELETE FROM auto_peaks WHERE id = ?", [orphan_id])
    end

    q_to_id
end

function persist_analysis!(db::SQLite.DB, exposure_id::Int,
                            q_full::Vector{Float64},
                            I_full::Vector{Float64},
                            peaks_result::NamedTuple,
                            candidates::Vector{<:Himalaya.Index},
                            group_indices::Vector{<:Himalaya.Index},
                            eff::NamedTuple)
    SQLite.transaction(db) do
    _persist_analysis_inner!(db, exposure_id, q_full, I_full, peaks_result,
                             candidates, group_indices, eff)
    end
end

function _persist_analysis_inner!(db::SQLite.DB, exposure_id::Int,
                                   q_full::Vector{Float64},
                                   I_full::Vector{Float64},
                                   peaks_result::NamedTuple,
                                   candidates::Vector{<:Himalaya.Index},
                                   group_indices::Vector{<:Himalaya.Index},
                                   eff::NamedTuple)

    EXCLUDE_TOL = 1e-6   # tolerance for q-value matching in speculative re-resolution

    # Snapshot the custom group's members by *semantic identity* (phase + basis)
    # BEFORE we delete the indices rows. Without this, the deletion below
    # invalidates every PK in `index_group_members`, leaving the active set
    # empty after reanalysis. Restricted to kind='auto' members because
    # speculative members keep stable ids across the wipe (see below).
    custom_member_identities = Tables.rowtable(DBInterface.execute(db, """
        SELECT g.id AS group_id, i.phase, i.basis
        FROM index_groups g
        JOIN index_group_members m ON m.group_id = g.id
        JOIN indices i ON i.id = m.index_id
        WHERE g.exposure_id = ? AND g.kind = 'custom' AND i.kind = 'auto'
        """, [exposure_id]))

    # Snapshot speculative indices' index_peaks rows by q-value. Auto peaks are
    # about to be re-detected with new ids, so the FK in index_peaks would
    # dangle without this remap.
    speculative_assignments = Tables.rowtable(DBInterface.execute(db, """
        SELECT i.id AS index_id, ip.ratio_position, ip.peak_kind,
               ip.peak_id AS old_peak_id,
               COALESCE(ap.q, pc.q) AS q_value
        FROM indices i
        JOIN index_peaks ip ON ip.index_id = i.id
        LEFT JOIN auto_peaks ap     ON ap.id = ip.peak_id     AND ip.peak_kind = 'auto'
        LEFT JOIN peak_curations pc ON pc.id = ip.peak_id     AND ip.peak_kind = 'curation'
        WHERE i.exposure_id = ? AND i.kind = 'speculative'
        """, [exposure_id]))
    speculative_index_ids = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, phase, basis FROM indices WHERE exposure_id = ? AND kind = 'speculative'",
        [exposure_id]))

    # Remove prior auto peaks, auto indices, and auto groups for this exposure.
    # Curation rows (peak_curations) are preserved. Speculative indices are
    # preserved but their index_peaks are wiped and re-resolved below.
    DBInterface.execute(db,
        "DELETE FROM index_group_members WHERE group_id IN
         (SELECT id FROM index_groups WHERE exposure_id = ? AND kind = 'auto')",
        [exposure_id])
    DBInterface.execute(db,
        "DELETE FROM index_groups WHERE exposure_id = ? AND kind = 'auto'",
        [exposure_id])
    # Clear custom-group memberships of *auto* indices only — we'll re-attach
    # by semantic identity below. Speculative members keep their membership
    # rows because their index ids are stable.
    DBInterface.execute(db, """
        DELETE FROM index_group_members
        WHERE index_id IN (
          SELECT id FROM indices WHERE exposure_id = ? AND kind = 'auto'
        )""", [exposure_id])
    DBInterface.execute(db, """
        DELETE FROM index_peaks
        WHERE index_id IN (SELECT id FROM indices WHERE exposure_id = ?)
        """, [exposure_id])
    DBInterface.execute(db,
        "DELETE FROM indices WHERE exposure_id = ? AND kind = 'auto'", [exposure_id])

    # Build O(1) q→(peak_id, peak_kind) lookup from effective peaks.
    eff_lookup = Dict{Float64, Tuple{Int, Symbol}}()
    for i in eachindex(eff.q)
        eff_lookup[eff.q[i]] = (eff.peak_id[i], eff.peak_kind[i])
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
            peak_info = get(eff_lookup, qval, nothing)
            peak_info === nothing && continue
            (peak_id, peak_kind) = peak_info
            ideal = ratios_normed[rpos] * Himalaya.basis(idx)
            resid = abs(qval - ideal)
            DBInterface.execute(db,
                """INSERT OR IGNORE INTO index_peaks
                     (index_id, peak_id, peak_kind, ratio_position, residual)
                   VALUES (?, ?, ?, ?, ?)""",
                [db_id, peak_id, String(peak_kind), rpos, resid])
        end
    end

    # ── Re-resolve speculative indices' peak references ─────────────────────
    # Auto peaks were just re-detected with new ids, so each speculative
    # index's `index_peaks` rows need fresh peak_ids. We match by q-value
    # (auto peaks: lookup in eff_lookup; curation-add peaks: stable id).
    if !isempty(speculative_index_ids)
        # Build a lookup: peak_id → (q, sharpness, peak_kind) from eff
        eff_by_id = Dict{Int, NamedTuple}()
        for i in eachindex(eff.q)
            eff_by_id[eff.peak_id[i]] = (q = eff.q[i], sharpness = eff.sharpness[i],
                                          peak_kind = eff.peak_kind[i])
        end

        function _resolve_peak_id(snap_row)
            snap_kind = String(snap_row.peak_kind)
            if snap_kind == "curation"
                # Curation-add peaks keep their id across reanalysis.
                pid = Int(snap_row.old_peak_id)
                return haskey(eff_by_id, pid) ? pid : nothing
            else
                # Auto peak — match by q-value within tolerance.
                qv = Float64(snap_row.q_value)
                # Exact match via eff_lookup
                info = get(eff_lookup, qv, nothing)
                info !== nothing && info[2] == :auto && return info[1]
                # Fuzzy match (peak shifted slightly under new detection)
                tol_val = max(EXCLUDE_TOL, abs(qv) * 0.001)
                best, best_delta = nothing, Inf
                for i in eachindex(eff.q)
                    eff.peak_kind[i] == :auto || continue
                    d = abs(eff.q[i] - qv)
                    if d < best_delta && d <= tol_val
                        best_delta = d
                        best = eff.peak_id[i]
                    end
                end
                return best
            end
        end

        # Group snapshot rows by index_id.
        by_index = Dict{Int, Vector{NamedTuple}}()
        for r in speculative_assignments
            push!(get!(by_index, Int(r.index_id), NamedTuple[]), r)
        end

        for ix_row in speculative_index_ids
            ix_id      = Int(ix_row.id)
            phase_name = String(ix_row.phase)
            P          = resolve_phase(phase_name)
            P === nothing && continue  # malformed phase string — leave as-is

            ratios_unnorm = Himalaya.phaseratios(P)
            ratios_normed = Himalaya.phaseratios(P; normalize = true)
            n             = length(ratios_normed)

            snaps = get(by_index, ix_id, NamedTuple[])
            ratio_to_peak = Dict{Int, Tuple{Int, Float64, Float64}}()  # rpos → (peak_id, q, sharpness)
            for s in snaps
                rpos = Int(s.ratio_position)
                pid  = _resolve_peak_id(s)
                pid === nothing && continue
                pr = get(eff_by_id, pid, nothing)
                pr === nothing && continue
                sharp = Float64(pr.sharpness)
                ratio_to_peak[rpos] = (pid, Float64(pr.q), sharp)
            end

            # Determine a working basis we can use to discover *new* peaks that
            # fit the speculative's predicted ratio positions. Two paths:
            #  - Survived snapshot has ≥ 2 peaks: refit basis from those.
            #  - Stale-recovery path (< 2 survived): use the original q-values
            #    from the snapshot — they encode the user's hypothesis even if
            #    the underlying peaks no longer exist.
            basis_for_snap = if length(ratio_to_peak) >= 2
                # Nominal path: use currently-resolved peaks to refit basis.
                rpos_seed = sort(collect(keys(ratio_to_peak)))
                qvals_seed = [ratio_to_peak[r][2] for r in rpos_seed]
                ratios_unnorm[rpos_seed] \ qvals_seed
            else
                # Stale-recovery path: use stored snapshot q-values, but skip
                # any whose underlying peak has been deleted (q_value = NULL/missing).
                valid_snaps = filter(s -> !ismissing(s.q_value), snaps)
                if length(valid_snaps) >= 2
                    snap_pairs = sort([(Int(s.ratio_position), Float64(s.q_value)) for s in valid_snaps])
                    rpos_seed  = [first(p) for p in snap_pairs]
                    qvals_seed = [last(p)  for p in snap_pairs]
                    ratios_unnorm[rpos_seed] \ qvals_seed
                else
                    # Last resort: use the persisted basis on the indices row.
                    if ismissing(ix_row.basis)
                        nothing
                    else
                        stored = Float64(ix_row.basis)
                        stored > 0 ? stored : nothing
                    end
                end
            end

            # Auto-discovery pass: pull in any current peak that fits an unfilled
            # ratio position within snap tolerance and isn't already claimed
            # by another ratio of this index.
            if basis_for_snap !== nothing && basis_for_snap > 0
                claimed_pids = Set{Int}(p[1] for p in values(ratio_to_peak))
                for rpos in 1:n
                    haskey(ratio_to_peak, rpos) && continue
                    predicted_q = basis_for_snap * ratios_normed[rpos]
                    best_pid::Union{Int, Nothing} = nothing
                    best_q     = 0.0
                    best_sharp = 0.0
                    best_relresid = SNAP_TOL
                    for i in eachindex(eff.q)
                        pid = eff.peak_id[i]
                        pid in claimed_pids && continue
                        qv = eff.q[i]
                        relresid = abs(qv - predicted_q) / predicted_q
                        if relresid <= best_relresid
                            best_relresid = relresid
                            best_pid      = pid
                            best_q        = qv
                            best_sharp    = eff.sharpness[i]
                        end
                    end
                    if best_pid !== nothing
                        ratio_to_peak[rpos] = (best_pid, best_q, best_sharp)
                        push!(claimed_pids, best_pid)
                    end
                end
            end

            if length(ratio_to_peak) < 2
                # Still not enough peaks to keep the index live — mark stale.
                DBInterface.execute(db,
                    "UPDATE indices SET status = 'stale' WHERE id = ?", [ix_id])
                continue
            end

            rpos_sorted = sort(collect(keys(ratio_to_peak)))
            qvals      = [ratio_to_peak[r][2] for r in rpos_sorted]
            sharpvals  = [ratio_to_peak[r][3] for r in rpos_sorted]

            # Recompute basis/r²/d/score using new peak values
            observed_ratios_used = ratios_unnorm[rpos_sorted]
            new_basis = observed_ratios_used \ qvals

            peaks_sv     = SparseArrays.SparseVector{Float64, Int}(n, rpos_sorted, qvals)
            sharpness_sv = SparseArrays.SparseVector{Float64, Int}(n, rpos_sorted, sharpvals)
            new_idx = Himalaya.Index{P}(new_basis, peaks_sv, sharpness_sv)
            fit_result = Himalaya.fit(new_idx)
            new_score = Himalaya.score(new_idx)

            DBInterface.execute(db,
                """UPDATE indices SET basis = ?, score = ?, r_squared = ?, lattice_d = ?,
                                       status = 'candidate'
                   WHERE id = ?""",
                [new_basis, new_score, fit_result.R², fit_result.d, ix_id])

            for (rpos, (pid, qv, _)) in ratio_to_peak
                # Determine peak_kind from eff_by_id
                pk_info = get(eff_by_id, pid, nothing)
                pk_str = if pk_info !== nothing
                    String(pk_info.peak_kind)
                else
                    "auto"
                end
                ideal = ratios_unnorm[rpos] * new_basis
                resid = abs(qv - ideal)
                DBInterface.execute(db,
                    """INSERT OR IGNORE INTO index_peaks
                         (index_id, peak_id, peak_kind, ratio_position, residual)
                       VALUES (?, ?, ?, ?, ?)""",
                    [ix_id, pid, pk_str, rpos, resid])
            end
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

    end

    # If any custom group survived with at least one member (auto-reattached
    # OR speculative), it remains the active set — demote the freshly-created
    # auto group.
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

function get_peaks_for_exposure(db::SQLite.DB, exposure_id::Int)
    Tables.rowtable(DBInterface.execute(db, """
        SELECT a.id, a.exposure_id, a.q, a.intensity, a.prominence, a.sharpness,
               'auto' AS source,
               CASE WHEN c.q IS NOT NULL THEN 1 ELSE 0 END AS excluded
        FROM auto_peaks a
        LEFT JOIN peak_curations c
            ON c.exposure_id = a.exposure_id
           AND c.kind = 'exclude'
           AND ABS(c.q - a.q) <= MAX(1e-6, ABS(a.q) * 0.001)
        WHERE a.exposure_id = ?
        UNION ALL
        SELECT id, exposure_id, q, NULL AS intensity, NULL AS prominence, NULL AS sharpness,
               'manual' AS source, 0 AS excluded
        FROM peak_curations
        WHERE exposure_id = ? AND kind = 'add'
        ORDER BY q
    """, [exposure_id, exposure_id]))
end

"""
    _serialized_indices_bytes(db, exposure_id) -> Int

Compute the serialized size (in bytes) of the indices for an exposure. Used
by the slow-path `analyze_run` event to emit `post_state_size_bytes` for
observability.

Intentionally distinct from `_serialized_indices_for_broadcast` below: this
helper returns just the bytes of a thin metadata-only projection (no joined
peaks, no predicted_q, no ngc), which is all the telemetry payload needs.
The broadcast helper builds the full enriched cache-replay payload used by
curation routes — it joins index_peaks and computes predicted_q per index.
Consolidating the two would force `analyze_run` (which fires on every
exposure analysis) to pay the heavier query just to measure size. The
duplication here is structural, not redundant.
"""
function _serialized_indices_bytes(db::SQLite.DB, exposure_id::Int)::Int
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, phase, basis, score, r_squared, lattice_d, status, kind, inputs_hash FROM indices WHERE exposure_id = ?",
        [exposure_id]))
    arr = [Dict(
        :id          => Int(r.id),
        :phase       => ismissing(r.phase) ? nothing : String(r.phase),
        :basis       => ismissing(r.basis) ? nothing : Float64(r.basis),
        :score       => ismissing(r.score) ? nothing : Float64(r.score),
        :r_squared   => ismissing(r.r_squared) ? nothing : Float64(r.r_squared),
        :lattice_d   => ismissing(r.lattice_d) ? nothing : Float64(r.lattice_d),
        :status      => ismissing(r.status) ? nothing : String(r.status),
        :kind        => ismissing(r.kind) ? nothing : String(r.kind),
        :inputs_hash => ismissing(r.inputs_hash) ? nothing : String(r.inputs_hash),
    ) for r in rows]
    sizeof(JSON3.write(arr))
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
    synthesize_peaks_result(db, exposure_id, q, I) -> NamedTuple

Reconstruct a peaks_result NamedTuple from the persisted `auto_peaks` rows,
matching the output shape of `Himalaya.findpeaks`. Used when findpeaks is
skipped due to unchanged trace_hash but indexpeaks still needs to run
(e.g. peak set changed via curation).

Field `findpeaks_index` from legacy R2 migration rows may be NULL; falls back
to `argmin` nearest-grid-index for those rows (heals on next diff_update).
"""
function synthesize_peaks_result(db::SQLite.DB, exposure_id::Int,
                                  q::Vector{Float64}, I::Vector{Float64})
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT q, prominence, sharpness, findpeaks_index
           FROM auto_peaks WHERE exposure_id = ? ORDER BY q""", [exposure_id]))
    qs    = Float64[Float64(r.q)                                           for r in rows]
    proms = Float64[ismissing(r.prominence) ? 0.0 : Float64(r.prominence)  for r in rows]
    shs   = Float64[ismissing(r.sharpness)  ? 0.0 : Float64(r.sharpness)   for r in rows]
    idxs  = Int[
        # Prefer the persisted findpeaks_index (exact local-maximum sample).
        # Fall back to nearest-grid-index for legacy rows from R2 migration
        # where findpeaks_index is NULL — these heal on next diff_update.
        ismissing(r.findpeaks_index) ?
            argmin(abs.(q .- Float64(r.q))) :
            Int(r.findpeaks_index)
        for r in rows
    ]
    (q = qs, indices = idxs, prominence = proms, sharpness = shs)
end

# ──────────────────────────────────────────────────────────────────────────────
# Lightweight DB-only readers used by analyze_exposure!'s fast-skip path.
# Kept separate from the body so they can be reused (and so the fast path is
# obviously file-I/O-free).
# ──────────────────────────────────────────────────────────────────────────────

function read_trace_hash(db::SQLite.DB, exposure_id::Int)::Union{String, Nothing}
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT trace_hash FROM exposures WHERE id = ?", [exposure_id]))
    isempty(rows) && return nothing
    ismissing(rows[1].trace_hash) ? nothing : String(rows[1].trace_hash)
end

function read_inputs_hash(db::SQLite.DB, exposure_id::Int)::Union{String, Nothing}
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT analysis_inputs_hash FROM exposures WHERE id = ?", [exposure_id]))
    isempty(rows) && return nothing
    ismissing(rows[1].analysis_inputs_hash) ? nothing : String(rows[1].analysis_inputs_hash)
end

function count_auto_peaks(db::SQLite.DB, exposure_id::Int)::Int
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM auto_peaks WHERE exposure_id = ?", [exposure_id]))
    Int(rows[1].c)
end

function count_indices(db::SQLite.DB, exposure_id::Int)::Int
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM indices WHERE exposure_id = ?", [exposure_id]))
    Int(rows[1].c)
end

function any_add_curations(db::SQLite.DB, exposure_id::Int)::Bool
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 AS x FROM peak_curations WHERE exposure_id = ? AND kind = 'add' LIMIT 1",
        [exposure_id]))
    !isempty(rows)
end

"""
    hash_peak_set_from_db(db, exposure_id) -> String

Compute the analysis_inputs_hash from `auto_peaks` + `exclude` curations alone,
without loading the trace. Equivalent to `hash_peak_set(effective_peaks(...))`
when no `add` curations exist for the exposure. Caller MUST verify
`!any_add_curations(...)` first — `add` rows need sharpness sampled from the
trace, which requires file I/O.

Mirrors the `hash_peak_set` byte encoding exactly: sort by q, write 16 bytes
per peak (8 q + 8 sharpness, native order), SHA-256.
"""
function hash_peak_set_from_db(db::SQLite.DB, exposure_id::Int)::String
    auto = Tables.rowtable(DBInterface.execute(db,
        "SELECT q, sharpness FROM auto_peaks WHERE exposure_id = ?", [exposure_id]))
    excludes = Tables.rowtable(DBInterface.execute(db,
        "SELECT q FROM peak_curations WHERE exposure_id = ? AND kind = 'exclude'",
        [exposure_id]))
    tol(q) = max(1e-6, abs(q) * 0.001)
    is_excluded(qv) = any(e -> abs(Float64(e.q) - qv) <= tol(qv), excludes)

    qs  = Float64[]
    shs = Float64[]
    for r in auto
        qv = Float64(r.q)
        is_excluded(qv) && continue
        push!(qs, qv)
        push!(shs, ismissing(r.sharpness) ? 0.0 : Float64(r.sharpness))
    end

    n = length(qs)
    buf = Vector{UInt8}(undef, 16n)
    perm = sortperm(qs)
    for (i, k) in enumerate(perm)
        q_bytes  = reinterpret(UInt8, [qs[k]])
        sh_bytes = reinterpret(UInt8, [shs[k]])
        copyto!(buf, 16(i-1) + 1, q_bytes)
        copyto!(buf, 16(i-1) + 9, sh_bytes)
    end
    bytes2hex(SHA.sha256(buf))
end

"""
    analyze_exposure!(db, exposure_id, analysis_dir; trace_known_unchanged=false)

Load the .dat file for `exposure_id`, run findpeaks + indexpeaks,
auto-group results, and persist everything to the DB.

Hash-guarded: findpeaks is skipped when `trace_hash` matches the persisted
value AND auto_peaks already exist. indexpeaks is skipped when
`analysis_inputs_hash` matches the persisted value AND indices already exist.

When `trace_known_unchanged=true`, the caller asserts that the .dat file has
not changed since the last `trace_hash` was computed (e.g. inside a curation
route handler that does not touch the trace). In that case, if no `add`
curations exist for the exposure, the function performs ZERO file I/O —
inputs hash is computed from the DB alone and matched against the stored
value; if they match and indices already exist, the call is effectively free.

The .dat filename is constructed from `exposures.filename` and the integration
pattern stored in the experiment's config (defaults to `{name}.dat` for
experiments without an explicit config).
"""
function analyze_exposure!(db::SQLite.DB, exposure_id::Int, analysis_dir::String;
                            trace_known_unchanged::Bool=false,
                            defer_broadcast::Bool=false)
    t0 = time()

    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT e.filename, x.id AS experiment_id
           FROM exposures e
           JOIN samples s ON s.id = e.sample_id
           JOIN experiments x ON x.id = s.experiment_id
           WHERE e.id = ?""", [exposure_id]))
    isempty(rows) && error("exposure $exposure_id not found")
    filename      = rows[1].filename
    experiment_id = rows[1].experiment_id

    cfg = config_from_db(db, experiment_id)
    pattern_filename = replace(cfg.integration_pattern, "{name}" => filename)
    dat_path = joinpath(analysis_dir, pattern_filename)

    # DB-only reads first so the fast path can short-circuit before any file I/O.
    stored_trace_hash  = read_trace_hash(db, exposure_id)
    stored_inputs_hash = read_inputs_hash(db, exposure_id)
    autopeaks_count    = count_auto_peaks(db, exposure_id)
    indices_count      = count_indices(db, exposure_id)

    new_trace_hash = nothing
    findpeaks_skipped = false
    if trace_known_unchanged
        new_trace_hash = stored_trace_hash
        findpeaks_skipped = autopeaks_count > 0
    end

    # Fast path: no file I/O at all when caller guarantees trace unchanged AND
    # no add curations (which require sharpness from the trace).
    if findpeaks_skipped && !any_add_curations(db, exposure_id)
        new_inputs_hash = hash_peak_set_from_db(db, exposure_id)
        indexpeaks_skipped = (stored_inputs_hash == new_inputs_hash) && (indices_count > 0)
        if indexpeaks_skipped
            # Fast-path no-op: skip the durable analyze_run row (M0.4 already drops the SSE frame). The hashes prove no-op-ness; durable counting offers no load-bearing value here.
            return
        end
    end

    # Slow path: from here on, behavior mirrors the pre-refactor implementation.
    isfile(dat_path) || error("dat file not found: $dat_path")

    if new_trace_hash === nothing
        new_trace_hash = hash_trace_file(dat_path)
        findpeaks_skipped = (stored_trace_hash == new_trace_hash) && (autopeaks_count > 0)
    end

    q, I, σ = load_dat(dat_path)
    fresh_peaks_result = nothing
    if !findpeaks_skipped
        fresh_peaks_result = Himalaya.findpeaks(q, I, σ)
        diff_update_auto_peaks!(db, exposure_id, fresh_peaks_result, I)
        DBInterface.execute(db,
            "UPDATE exposures SET trace_hash = ? WHERE id = ?",
            [new_trace_hash, exposure_id])
    end

    eff = effective_peaks(db, exposure_id, q, I)
    new_inputs_hash = hash_peak_set(eff)

    indexpeaks_skipped = (stored_inputs_hash == new_inputs_hash) && (indices_count > 0)

    if !indexpeaks_skipped
        peaks_result_for_persist = fresh_peaks_result === nothing ?
            synthesize_peaks_result(db, exposure_id, q, I) :
            fresh_peaks_result
        candidates = Himalaya.indexpeaks(eff.q, eff.sharpness)
        group = auto_group(candidates)
        persist_analysis!(db, exposure_id, q, I, peaks_result_for_persist,
                          candidates, group, eff)
        DBInterface.execute(db,
            "UPDATE exposures SET analysis_inputs_hash = ? WHERE id = ?",
            [new_inputs_hash, exposure_id])
        DBInterface.execute(db,
            "UPDATE indices SET inputs_hash = ? WHERE exposure_id = ?",
            [new_inputs_hash, exposure_id])
    end

    duration_ms = round(Int, (time() - t0) * 1000)
    post_state_size_bytes = _serialized_indices_bytes(db, exposure_id)
    apply_event!(db, _system_request();
        kind        = "analyze_run",
        entity_type = "exposure",
        entity_id   = exposure_id,
        payload     = Dict(
            :trace_hash_before     => stored_trace_hash,
            :trace_hash_after      => new_trace_hash,
            :inputs_hash_before    => stored_inputs_hash,
            :inputs_hash_after     => new_inputs_hash,
            :findpeaks_skipped     => findpeaks_skipped,
            :indexpeaks_skipped    => indexpeaks_skipped,
            :duration_ms           => duration_ms,
            :effective_peaks_count => length(eff.q),
            :post_state_size_bytes => post_state_size_bytes,
        ),
        defer_broadcast = defer_broadcast)
end

"""
    _serialized_indices_for_broadcast(db, exposure_id) -> Vector{Dict}

Build the indices array for the SSE `post_state` payload — same shape as the
`IndexEntry[]` consumed by the frontend: `id`, `exposure_id`, `phase`, `basis`,
`score`, `r_squared`, `lattice_d`, `status`, `kind`, `inputs_hash`, plus the
joined `peaks` (each `peak_id`, `ratio_position`, `residual`, `q_observed`)
and the `predicted_q` array. `ngc` mirrors what `GET /api/exposures/:id/indices`
emits and is computed lazily by the frontend in some places, but we include it
here too so the cache can stay authoritative.
"""
function _serialized_indices_for_broadcast(db::SQLite.DB, exposure_id::Int)
    indices = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM indices WHERE exposure_id = ? ORDER BY score DESC", [exposure_id]))
    map(indices) do ix
        peak_rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT ip.peak_id, ip.ratio_position, ip.residual,
                      COALESCE(ap.q, pc.q) AS q_observed
               FROM index_peaks ip
               LEFT JOIN auto_peaks ap     ON ap.id = ip.peak_id AND ip.peak_kind = 'auto'
               LEFT JOIN peak_curations pc ON pc.id = ip.peak_id AND ip.peak_kind = 'curation'
               WHERE ip.index_id = ? ORDER BY ip.ratio_position""",
            [Int(ix.id)]))
        d = row_to_json(ix)
        d[:peaks]       = rows_to_json(peak_rows)
        d[:predicted_q] = predicted_q_for_phase(String(ix.phase), Float64(ix.basis))
        d[:ngc]         = _ngc_for_phase(String(ix.phase), ix.lattice_d)
        d
    end
end
