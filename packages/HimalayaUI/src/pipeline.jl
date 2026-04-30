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
    SQLite.transaction(db) do
    _persist_analysis_inner!(db, exposure_id, q_full, I_full, peaks_result,
                             candidates, group_indices)
    end
end

function _persist_analysis_inner!(db::SQLite.DB, exposure_id::Int,
                                   q_full::Vector{Float64},
                                   I_full::Vector{Float64},
                                   peaks_result::NamedTuple,
                                   candidates::Vector{<:Himalaya.Index},
                                   group_indices::Vector{<:Himalaya.Index})

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
        SELECT i.id AS index_id, ip.ratio_position, p.q AS q_value, p.source AS peak_source, p.id AS old_peak_id
        FROM indices i
        JOIN index_peaks ip ON ip.index_id = i.id
        JOIN peaks p ON p.id = ip.peak_id
        WHERE i.exposure_id = ? AND i.kind = 'speculative'
        """, [exposure_id]))
    speculative_index_ids = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, phase, basis FROM indices WHERE exposure_id = ? AND kind = 'speculative'",
        [exposure_id]))

    # Remove prior auto peaks, auto indices, and auto groups for this exposure.
    # Manual peaks (source='manual') are preserved. Speculative indices are
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

    # Manual peaks survive the auto-peak wipe; map their stable ids by q so
    # candidate indices that incorporate them get proper `index_peaks` rows.
    #
    # Float64 equality is load-bearing here: `analyze_exposure!` reads these
    # same q-values out of SQLite and `vcat`s them into the indexpeaks input,
    # so the q's stored in `IndexEntry.peaks` are bit-identical copies of the
    # `r.q` values queried below — no float drift on the round-trip. If
    # `Himalaya.indexpeaks` ever starts snapping or transforming its input
    # qs, this lookup will silently miss and the manual-peak `index_peaks`
    # rows won't be written. The "incorporates manual peaks" testset in
    # test_pipeline.jl will catch that regression.
    manual_peak_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, q FROM peaks WHERE exposure_id = ? AND source = 'manual'",
        [exposure_id]))
    for r in manual_peak_rows
        q_to_peak_id[Float64(r.q)] = Int(r.id)
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

    # ── Re-resolve speculative indices' peak references ─────────────────────
    # Auto peaks were just re-detected with new ids, so each speculative
    # index's `index_peaks` rows need fresh peak_ids. We match by q-value
    # (auto peaks: lookup in q_to_peak_id; manual peaks: by stable id).
    # An index loses ratio assignments whose q-values no longer correspond to
    # any current peak. If ≥ 2 assignments survive, we recompute its
    # basis/r²/d/score; if not, we mark status='stale' but keep the row so
    # the user can decide what to do.
    if !isempty(speculative_index_ids)
        # Build a fresh q→peak_id lookup that includes auto + manual peaks.
        all_peak_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, q, sharpness, source, excluded FROM peaks WHERE exposure_id = ?",
            [exposure_id]))
        q_to_id_full   = Dict{Int, NamedTuple}()  # peak_id → full row
        for pr in all_peak_rows
            q_to_id_full[Int(pr.id)] = pr
        end
        manual_q_to_id = Dict{Int, Int}()
        for pr in all_peak_rows
            String(pr.source) == "manual" && (manual_q_to_id[Int(pr.id)] = Int(pr.id))
        end

        function _resolve_peak_id(snap_row)
            if String(snap_row.peak_source) == "manual"
                # Manual peaks keep their id across reanalysis.
                pid = Int(snap_row.old_peak_id)
                return haskey(q_to_id_full, pid) ? pid : nothing
            else
                # Auto peak — match by q-value within tolerance.
                qv = Float64(snap_row.q_value)
                pid = get(q_to_peak_id, qv, nothing)
                pid !== nothing && return pid
                # Fall back to fuzzy match (peak shifted slightly under new
                # detection): find closest current auto peak within EXCLUDE_TOL.
                # Note this is tighter than `SNAP_TOL` (≈0.25%) — peaks that
                # drift more than 0.1% won't bind here, but the auto-discovery
                # loop further down uses `SNAP_TOL` as a safety net so they're
                # still picked up under their predicted ratio position.
                tol = max(EXCLUDE_TOL, abs(qv) * 0.001)
                best, best_delta = nothing, Inf
                for pr in all_peak_rows
                    String(pr.source) == "auto" || continue
                    d = abs(Float64(pr.q) - qv)
                    if d < best_delta && d <= tol
                        best_delta = d
                        best = Int(pr.id)
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
                pr = q_to_id_full[pid]
                # Treat newly-excluded peaks as if they don't exist.
                (pr.excluded == 1 || pr.excluded === true) && continue
                sharp = pr.sharpness === nothing || ismissing(pr.sharpness) ? 0.0 : Float64(pr.sharpness)
                ratio_to_peak[rpos] = (pid, Float64(pr.q), sharp)
            end

            # Determine a working basis we can use to discover *new* peaks that
            # fit the speculative's predicted ratio positions. Two paths:
            #  - Survived snapshot has ≥ 2 peaks: refit basis from those.
            #  - Stale-recovery path (< 2 survived): use the original q-values
            #    from the snapshot — they encode the user's hypothesis even if
            #    the underlying peaks no longer exist. This is what makes
            #    "speculative goes stale → user adds a manual peak → reanalyze
            #    pulls it in" actually work.
            basis_for_snap = if length(ratio_to_peak) >= 2
                rpos_seed = sort(collect(keys(ratio_to_peak)))
                qvals_seed = [ratio_to_peak[r][2] for r in rpos_seed]
                ratios_unnorm[rpos_seed] \ qvals_seed
            elseif length(snaps) >= 2
                snap_pairs = sort([(Int(s.ratio_position), Float64(s.q_value)) for s in snaps])
                rpos_seed  = [first(p) for p in snap_pairs]
                qvals_seed = [last(p)  for p in snap_pairs]
                ratios_unnorm[rpos_seed] \ qvals_seed
            else
                # Once a speculative goes stale, its index_peaks (and therefore
                # the snapshot) is empty. The `basis` column on the indices row
                # itself was last set when the index was built, so use it as
                # the persisted record of user intent for auto-discovery.
                # `phaseratios` is normalized to ratio[1] = 1, so the ratio[1]
                # equivalent of `basis` is just `basis`. The DDL doesn't
                # enforce NOT NULL on `basis` (yet), so guard against it.
                if ismissing(ix_row.basis) || ix_row.basis === nothing
                    nothing
                else
                    stored = Float64(ix_row.basis)
                    stored > 0 ? stored : nothing
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
                    for pr in all_peak_rows
                        (pr.excluded == 1 || pr.excluded === true) && continue
                        pid = Int(pr.id)
                        pid in claimed_pids && continue
                        qv = Float64(pr.q)
                        relresid = abs(qv - predicted_q) / predicted_q
                        if relresid <= best_relresid
                            best_relresid = relresid
                            best_pid      = pid
                            best_q        = qv
                            best_sharp    = pr.sharpness === nothing || ismissing(pr.sharpness) ? 0.0 : Float64(pr.sharpness)
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
                ideal = ratios_unnorm[rpos] * new_basis
                resid = abs(qv - ideal)
                DBInterface.execute(db,
                    """INSERT OR IGNORE INTO index_peaks (index_id, peak_id, ratio_position, residual)
                       VALUES (?, ?, ?, ?)""",
                    [ix_id, pid, rpos, resid])
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

    end

    # If any custom group survived with at least one member (auto-reattached
    # OR speculative), it remains the active set — demote the freshly-created
    # auto group. Mirrors the post-curate semantics of `ensure_custom_group!`.
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

The .dat filename is constructed from `exposures.filename` and the integration
pattern stored in the experiment's config (defaults to `{name}.dat` for
experiments without an explicit config).
"""
function analyze_exposure!(db::SQLite.DB, exposure_id::Int, analysis_dir::String)
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT e.filename, x.id AS experiment_id
           FROM exposures e
           JOIN samples s ON s.id = e.sample_id
           JOIN experiments x ON x.id = s.experiment_id
           WHERE e.id = ?""", [exposure_id]))
    isempty(rows) && error("exposure $exposure_id not found")
    filename      = rows[1].filename
    experiment_id = rows[1].experiment_id

    cfg              = config_from_db(db, experiment_id)
    pattern_filename = replace(cfg.integration_pattern, "{name}" => filename)
    dat_path         = joinpath(analysis_dir, pattern_filename)
    isfile(dat_path) || error("dat file not found: $dat_path")

    q, I, σ      = load_dat(dat_path)
    peaks_result = Himalaya.findpeaks(q, I, σ)

    # Curation snapshot: previously-excluded auto q's, plus user-added manual q's.
    # Both adjust what the indexer sees relative to raw `findpeaks` output:
    #   - excluded auto peaks must not contribute to candidate scoring
    #     (otherwise the user's "this is noise" verdict gets ignored every
    #     reanalysis), and
    #   - manual peaks must be included so a peak the user marked at a
    #     predicted ratio position can land in `IndexEntry.peaks`.
    excluded_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT q FROM peaks WHERE exposure_id = ? AND source = 'auto' AND excluded = 1",
        [exposure_id]))
    excluded_qs = Float64[Float64(r.q) for r in excluded_rows]
    manual_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT q FROM peaks WHERE exposure_id = ? AND source = 'manual'",
        [exposure_id]))
    manual_q = Float64[Float64(r.q) for r in manual_rows]

    is_excluded(qv::Float64) = any(
        eq -> abs(eq - qv) <= max(1e-6, abs(qv) * 0.001), excluded_qs)

    keep = [i for i in eachindex(peaks_result.q) if !is_excluded(peaks_result.q[i])]
    auto_q_kept     = peaks_result.q[keep]
    auto_sharp_kept = peaks_result.sharpness[keep]

    candidates = if isempty(manual_q)
        Himalaya.indexpeaks(auto_q_kept, auto_sharp_kept)
    else
        # Sharpness for a manual peak: nearest-grid lookup against the full
        # SG-second-derivative trace. One O(N) sharpness pass regardless of
        # how many manual peaks exist, and consistent with the unit `findpeaks`
        # uses for its own peaks so they mix correctly in `consistency`.
        sharp_full   = Himalaya.sharpness(I)
        manual_sharp = Float64[sharp_full[argmin(abs.(q .- mq))] for mq in manual_q]
        all_q        = vcat(auto_q_kept, manual_q)
        all_sharp    = vcat(auto_sharp_kept, manual_sharp)
        perm         = sortperm(all_q)
        Himalaya.indexpeaks(all_q[perm], all_sharp[perm])
    end

    group = auto_group(candidates)

    persist_analysis!(db, exposure_id, q, I, peaks_result, candidates, group)
end
