# Matching user curations (`exclude` / `add` rows in peak_curations) to auto
# peaks. Every Julia and SQL site that attaches a curation to an auto peak goes
# through this file; the tolerance constant lives only here.

"""Grid points of |Δq| within which a curation attaches to an auto peak."""
const PEAK_MATCH_POINTS = 6

"""
    peak_match_tol(q, q_step) -> Float64

Largest |Δq| at which a curation at `q` attaches to an auto peak. `q_step` is
the trace's grid step (`exposures.q_step`); `nothing`/`missing` means unknown
and gives the original rule, `max(1e-6, 0.001·|q|)`.
"""
function peak_match_tol(q::Real, q_step)
    base = max(1e-6, 0.001 * abs(Float64(q)))
    (q_step === nothing || q_step === missing) && return base
    max(base, PEAK_MATCH_POINTS * Float64(q_step))
end

"""
    exposure_q_step(db, exposure_id) -> Union{Float64, Nothing}

The stored grid step, or `nothing` when unknown (NULL, non-positive or non-finite).
"""
function exposure_q_step(db::SQLite.DB, exposure_id::Integer)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT q_step FROM exposures WHERE id = ?", [Int(exposure_id)]))
    (isempty(rows) || ismissing(rows[1].q_step)) && return nothing
    s = Float64(rows[1].q_step)
    (isfinite(s) && s > 0) ? s : nothing
end

"""
    grid_step(q) -> Union{Float64, Nothing}

Median spacing of the q grid; `nothing` when it cannot be measured.
"""
function grid_step(q::AbstractVector{<:Real})
    length(q) >= 2 || return nothing
    d = sort!(Float64.(diff(q)))
    m = isodd(length(d)) ? d[(end + 1) ÷ 2] : (d[end ÷ 2] + d[end ÷ 2 + 1]) / 2
    (isfinite(m) && m > 0) ? m : nothing
end

# Greedy one-to-one assignment, closest pairs first. `blocked` autos are not
# available. Returns the mask of claimed auto peaks.
function _claim_nearest(auto_q, cur_q, q_step, blocked)
    pairs = Tuple{Float64, Int, Int}[]
    for (ci, cq) in enumerate(cur_q), (ai, aq) in enumerate(auto_q)
        blocked[ai] && continue
        d = abs(Float64(cq) - Float64(aq))
        d <= peak_match_tol(aq, q_step) && push!(pairs, (d, ai, ci))
    end
    sort!(pairs)
    claimed = falses(length(auto_q))
    used    = falses(length(cur_q))
    for (_, ai, ci) in pairs
        (claimed[ai] || used[ci]) && continue
        claimed[ai] = true
        used[ci]    = true
    end
    claimed
end

"""
    match_curations(auto_q, exclude_q, add_q, q_step) -> (excluded, superseded)

Boolean masks over `auto_q`.

With a known `q_step`, matching is one-to-one and nearest-first: an `exclude`
hides the one auto peak nearest to it within tolerance, and an `add` replaces
the nearest remaining auto peak within tolerance (the user's peak wins, so no
doublet). Excludes claim first, so an add cannot take an excluded peak.

With `q_step === nothing` the original rule applies: an exclude hides every
auto peak within `max(1e-6, 0.001·|q|)` and an add never replaces anything.
"""
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

"""
    curation_flags(db, exposure_id, auto_q) -> (excluded, superseded, adds)

Reads the exposure's curations and grid step, matches them against `auto_q`.
`adds` is the rowtable `(id, q)` of `add` curations.
"""
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

"""
    peak_listing(db, exposure_id; include_superseded = false) -> Vector{NamedTuple}

The rows the peaks API returns: every auto peak (with `excluded` 0/1) plus the
`add` curations (`source = "manual"`), sorted by q. Auto peaks replaced by an
`add` are left out unless `include_superseded`.
"""
function peak_listing(db::SQLite.DB, exposure_id::Integer; include_superseded::Bool = false)
    auto = Tables.rowtable(DBInterface.execute(db,
        """SELECT id, exposure_id, q, intensity, prominence, sharpness
           FROM auto_peaks WHERE exposure_id = ? ORDER BY id""", [Int(exposure_id)]))
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
