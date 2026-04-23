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

function persist_analysis!(db::SQLite.DB, exposure_id::Int,
                            q_full::Vector{Float64},
                            I_full::Vector{Float64},
                            peaks_result::NamedTuple,
                            candidates::Vector{<:Himalaya.Index},
                            group_indices::Vector{<:Himalaya.Index})

    # Remove prior auto peaks, indices, and auto groups for this exposure.
    # Manual peaks (source='manual') are preserved.
    DBInterface.execute(db,
        "DELETE FROM index_group_members WHERE group_id IN
         (SELECT id FROM index_groups WHERE exposure_id = ? AND kind = 'auto')",
        [exposure_id])
    DBInterface.execute(db,
        "DELETE FROM index_groups WHERE exposure_id = ? AND kind = 'auto'",
        [exposure_id])
    DBInterface.execute(db,
        "DELETE FROM index_peaks WHERE index_id IN
         (SELECT id FROM indices WHERE exposure_id = ?)", [exposure_id])
    DBInterface.execute(db,
        "DELETE FROM indices WHERE exposure_id = ?", [exposure_id])
    DBInterface.execute(db,
        "DELETE FROM peaks WHERE exposure_id = ? AND source = 'auto'", [exposure_id])

    # Persist auto peaks
    q_to_peak_id = Dict{Float64, Int}()
    for i in eachindex(peaks_result.q)
        qval       = peaks_result.q[i]
        full_idx   = peaks_result.indices[i]
        intensity  = I_full[full_idx]
        res = DBInterface.execute(db,
            "INSERT INTO peaks (exposure_id, q, intensity, prominence, sharpness, source)
             VALUES (?, ?, ?, ?, ?, 'auto')",
            [exposure_id, qval, intensity,
             peaks_result.prominence[i], peaks_result.sharpness[i]])
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
    candidates   = Himalaya.indexpeaks(peaks_result.q, peaks_result.prominence)
    group        = auto_group(candidates)

    persist_analysis!(db, exposure_id, q, I, peaks_result, candidates, group)
end
