const PEAK_MATCH_POINTS = 6

function peak_match_tol(q::Real, q_step)
    base = max(1e-6, 0.001 * abs(Float64(q)))
    (q_step === nothing || q_step === missing) && return base
    max(base, PEAK_MATCH_POINTS * Float64(q_step))
end

function exposure_q_step(db::SQLite.DB, exposure_id::Integer)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT q_step FROM exposures WHERE id = ?", [Int(exposure_id)]))
    (isempty(rows) || ismissing(rows[1].q_step)) && return nothing
    s = Float64(rows[1].q_step)
    (isfinite(s) && s > 0) ? s : nothing
end

function grid_step(q::AbstractVector{<:Real})
    length(q) >= 2 || return nothing
    d = sort!(Float64.(diff(q)))
    m = isodd(length(d)) ? d[(end + 1) ÷ 2] : (d[end ÷ 2] + d[end ÷ 2 + 1]) / 2
    (isfinite(m) && m > 0) ? m : nothing
end

function _assign_nearest(auto_q, cur_q, q_step, blocked)
    pairs = Tuple{Float64, Int, Int}[]
    for (ci, cq) in enumerate(cur_q), (ai, aq) in enumerate(auto_q)
        blocked[ai] && continue
        d = abs(Float64(cq) - Float64(aq))
        d <= peak_match_tol(aq, q_step) && push!(pairs, (d, ai, ci))
    end
    sort!(pairs)
    owner = zeros(Int, length(auto_q))
    used  = falses(length(cur_q))
    for (_, ai, ci) in pairs
        (owner[ai] != 0 || used[ci]) && continue
        owner[ai] = ci
        used[ci]  = true
    end
    owner
end

_claim_nearest(auto_q, cur_q, q_step, blocked) =
    BitVector(_assign_nearest(auto_q, cur_q, q_step, blocked) .!= 0)

function match_curations(auto_q::AbstractVector, exclude_q::AbstractVector,
                         add_q::AbstractVector, q_step)
    n = length(auto_q)
    if q_step === nothing || q_step === missing
        excluded = BitVector([any(e -> abs(Float64(e) - Float64(a)) <= peak_match_tol(a, nothing),
                                  exclude_q) for a in auto_q])
        return (excluded = excluded, superseded = falses(n))
    end
    excluded   = _claim_nearest(auto_q, exclude_q, q_step, falses(n))
    superseded = _claim_nearest(auto_q, add_q, q_step, excluded)
    (excluded = excluded, superseded = superseded)
end

function curation_flags(db::SQLite.DB, exposure_id::Integer, auto_q::AbstractVector)
    ex = Tables.rowtable(DBInterface.execute(db,
        "SELECT q FROM peak_curations WHERE exposure_id = ? AND kind = 'exclude'",
        [Int(exposure_id)]))
    adds = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, q FROM peak_curations WHERE exposure_id = ? AND kind = 'add' ORDER BY id",
        [Int(exposure_id)]))
    m = match_curations(auto_q, [Float64(e.q) for e in ex],
                        [Float64(a.q) for a in adds], exposure_q_step(db, exposure_id))
    (excluded = m.excluded, superseded = m.superseded, adds = adds)
end

function peak_listing(db::SQLite.DB, exposure_id::Integer; include_superseded::Bool = false)
    auto = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, exposure_id, q, intensity, prominence, sharpness FROM auto_peaks WHERE exposure_id = ? ORDER BY id",
        [Int(exposure_id)]))
    fl = curation_flags(db, exposure_id, [Float64(r.q) for r in auto])
    rows = NamedTuple[]
    for (i, r) in enumerate(auto)
        (fl.superseded[i] && !include_superseded) && continue
        push!(rows, (id = r.id, exposure_id = r.exposure_id, q = r.q,
                     intensity = r.intensity, prominence = r.prominence,
                     sharpness = r.sharpness, source = "auto",
                     excluded = Int(fl.excluded[i])))
    end
    for a in fl.adds
        push!(rows, (id = a.id, exposure_id = Int(exposure_id), q = a.q,
                     intensity = missing, prominence = missing, sharpness = missing,
                     source = "manual", excluded = 0))
    end
    sort!(rows; by = r -> Float64(r.q))
    rows
end

function exclude_ids_for_peak(db::SQLite.DB, exposure_id::Integer, q::Real)
    q0   = Float64(q)
    step = exposure_q_step(db, exposure_id)
    ex = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, q FROM peak_curations WHERE exposure_id = ? AND kind = 'exclude' ORDER BY id",
        [Int(exposure_id)]))
    base = peak_match_tol(q0, nothing)
    legacy = [Int(e.id) for e in ex if abs(Float64(e.q) - q0) <= base]
    step === nothing && return legacy

    auto = Tables.rowtable(DBInterface.execute(db,
        "SELECT q FROM auto_peaks WHERE exposure_id = ? ORDER BY id", [Int(exposure_id)]))
    aq = [Float64(r.q) for r in auto]
    exq = [Float64(e.q) for e in ex]
    ti = isempty(aq) ? 0 : argmin(abs.(aq .- q0))
    if ti != 0 && abs(aq[ti] - q0) <= base
        owner = _assign_nearest(aq, exq, step, falses(length(aq)))
        return owner[ti] == 0 ? Int[] : [Int(ex[owner[ti]].id)]
    end
    tol = peak_match_tol(q0, step)
    near = [i for i in eachindex(exq) if abs(exq[i] - q0) <= tol]
    ids = Set(legacy)
    isempty(near) || push!(ids, Int(ex[near[argmin([abs(exq[i] - q0) for i in near])]].id))
    collect(ids)
end
