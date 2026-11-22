import HimalayaUI

"""
One exposure with the user's curation reduced to q values.

- `adds`: peaks the user added by hand.
- `kept`: stored auto peaks the user did not exclude (and no add replaced).
- `excluded`: stored auto peaks the user excluded.
- `curated`: has at least one add or exclude. `kept` is only trusted as
  "the user looked and accepted" on curated exposures.
- `form_factor`: the sample was called form factor (no Bragg peaks expected).
- `group`: cross-validation fold key, an experiment id by default.
"""
struct LabeledExposure
    id::Int
    group::Int
    q::Vector{Float64}
    I::Vector{Float64}
    σ::Vector{Float64}
    adds::Vector{Float64}
    kept::Vector{Float64}
    excluded::Vector{Float64}
    curated::Bool
    form_factor::Bool
end

match_tol(q, q_step) = HimalayaUI.peak_match_tol(q, q_step)

# Candidate labels: 1 positive, 0 strong negative, -1 weak negative, -2 dropped.
const POS, NEG, WEAK, DROP = 1, 0, -1, -2

"""
    label_candidates(q, idx, adds, kept, excluded, curated, form_factor) -> Vector{Int}

A candidate is positive when it is the nearest candidate to an add or a kept
auto peak within the match tolerance. Other candidates inside the tolerance of
a positive are dropped (duplicates of a real peak, resolved by non-maximum
suppression at prediction time). A candidate within tolerance of an excluded
auto peak, or on a form-factor exposure, is a strong negative. Any other
candidate on a curated exposure is a weak negative: the user looked at the
exposure and did not add it, but users do not add every faint peak. Candidates
on exposures without curations are dropped.
"""
function label_candidates(q, idx, adds, kept, excluded, curated::Bool, form_factor::Bool)
    q_step = HimalayaUI.grid_step(q)
    lab = fill(DROP, length(idx))
    pos = form_factor ? collect(Float64, adds) : vcat(adds, kept)
    qi = q[idx]
    near_pos = Set{Int}()
    for p in pos
        isempty(qi) && break
        d = abs.(qi .- p)
        m = argmin(d)
        d[m] <= match_tol(p, q_step) && push!(near_pos, m)
    end
    for (m, x) in enumerate(qi)
        if m in near_pos
            lab[m] = POS
        elseif any(p -> abs(x - p) <= match_tol(p, q_step), pos)
            lab[m] = DROP
        elseif form_factor || any(e -> abs(x - e) <= match_tol(e, q_step), excluded)
            lab[m] = NEG
        elseif curated
            lab[m] = WEAK
        end
    end
    lab
end

"""
Candidate-level table for a set of exposures. Rows are candidates; `exposure`
maps each row to its index in the input, `group` to its fold key.
"""
struct Dataset
    X::Matrix{Float64}
    label::Vector{Int}
    exposure::Vector{Int}
    group::Vector{Int}
    cand::Vector{Vector{Int}}     # candidate grid indices per exposure
    hybrid::Vector{Vector{Int}}   # hybrid detector peaks per exposure
end

function build_dataset(exps::AbstractVector{LabeledExposure}; verbose::Bool = false)
    parts = Vector{Any}(undef, length(exps))
    Threads.@threads for k in eachindex(exps)
        e = exps[k]
        c = candidate_indices(e.q, e.I, e.σ)
        idx, X = candidate_features(e.q, e.I, e.σ)
        lab = label_candidates(e.q, idx, e.adds, e.kept, e.excluded, e.curated, e.form_factor)
        parts[k] = (idx, X, lab, c.hyb)
    end
    Xs = [p[2] for p in parts]
    exposure = reduce(vcat, [fill(k, length(parts[k][1])) for k in eachindex(parts)]; init = Int[])
    group = [exps[k].group for k in exposure]
    Dataset(isempty(Xs) ? zeros(0, length(FEATURE_NAMES)) : reduce(vcat, Xs),
            reduce(vcat, [p[3] for p in parts]; init = Int[]),
            exposure, group, [p[1] for p in parts], [p[4] for p in parts])
end
