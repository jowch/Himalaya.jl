"""
One labelled exposure with the user's curation reduced to q values.

- `adds`: peaks the user added by hand.
- `kept`: peaks the user saw, left alone and the hybrid detector found
  (`found_by` hybrid or both), from the latest review snapshot minus excluded
  and superseded ones.
- `excluded`: q of every exclude the user recorded, whether or not a stored
  auto peak still matches it.
- `labelled`: the user acted on this exposure's own peaks or indexing (a direct
  review snapshot). `kept` is trusted as "the user looked and accepted" on
  labelled exposures.
- `form_factor`: the sample was called form factor (no Bragg peaks expected).
- `group`: cross-validation fold key, an experiment id by default.
- `kept_model`: like `kept`, for peaks only a model found. They are not
  evidence for the hybrid, so they are never positives; neighbours of one are
  dropped rather than labelled.
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
    labelled::Bool
    form_factor::Bool
    kept_model::Vector{Float64}
end

LabeledExposure(id, group, q, I, σ, adds, kept, excluded, labelled, form_factor) =
    LabeledExposure(id, group, q, I, σ, adds, kept, excluded, labelled, form_factor, Float64[])

match_tol(q, q_step) = peak_match_tol(q, q_step)

# Candidate labels: 1 positive, 0 strong negative, -1 weak negative, -2 dropped.
const POS, NEG, WEAK, DROP = 1, 0, -1, -2

"""
    label_candidates(q, idx, adds, kept, excluded, labelled, form_factor; kept_model = Float64[])
        -> Vector{Int}

Rules apply in order, first match wins:

1. the nearest candidate to an add or a kept peak, within the match tolerance,
   is positive;
2. other candidates within the tolerance of a positive are dropped (duplicates
   of a real peak, resolved by non-maximum suppression at prediction time);
3. candidates within the tolerance of a `kept_model` peak are dropped;
4. candidates within the tolerance of an exclude are strong negatives, whether
   or not a stored auto peak matches the exclude;
5. remaining candidates on a form-factor exposure are strong negatives;
6. remaining candidates on a labelled exposure are weak negatives: the user
   looked at the exposure and did not add them, but users do not add every
   faint peak.

Other candidates are dropped. On a form-factor exposure only adds count as
positives, and `kept_model` is ignored.
"""
function label_candidates(q, idx, adds, kept, excluded, labelled::Bool, form_factor::Bool;
                          kept_model = Float64[])
    q_step = grid_step(q)
    lab = fill(DROP, length(idx))
    pos = form_factor ? collect(Float64, adds) : vcat(adds, kept)
    model_only = form_factor ? Float64[] : collect(Float64, kept_model)
    qi = q[idx]
    near(x, set) = any(p -> abs(x - p) <= match_tol(p, q_step), set)
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
        elseif near(x, pos) || near(x, model_only)
            lab[m] = DROP
        elseif form_factor || near(x, excluded)
            lab[m] = NEG
        elseif labelled
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
    feature_names::Vector{String} # columns of X
end

"""
    build_dataset(exps; features = DefaultFeatures()) -> Dataset

Candidates, features and labels for each exposure. A recipe can train on a
dataset whose columns include all of its feature names. With a `cache`
(`DiskFeatureCache`), an exposure's candidates, features and hybrid peaks are
read from its file when there is one and written after computing otherwise.
`float32` rounds features to Float32 (the cache's precision) and defaults to on with a cache.
"""
function build_dataset(exps::AbstractVector{LabeledExposure}; features::FeatureSet = DefaultFeatures(),
                       verbose::Bool = false, cache::Union{DiskFeatureCache, Nothing} = nothing,
                       float32::Bool = cache !== nothing)
    parts = Vector{Any}(undef, length(exps))
    ncol = length(feature_names(features))
    Threads.@threads for k in eachindex(exps)
        e = exps[k]
        path = cache === nothing ? nothing : cache_path(cache, e)
        hit = path === nothing ? nothing : read_features(path, ncol)
        if hit === nothing
            idx, X, hyb = candidates_and_features(features, e.q, e.I, e.σ)
            hyb_idx = hyb === nothing ? Int[] : hyb.indices
            float32 && (X = Float64.(Float32.(X)))
            if cache !== nothing
                path === nothing || try_write_features(path, idx, X, hyb_idx)
            end
        else
            idx, X, hyb_idx = hit
        end
        lab = label_candidates(e.q, idx, e.adds, e.kept, e.excluded, e.labelled, e.form_factor;
                              kept_model = e.kept_model)
        parts[k] = (idx, X, lab, hyb_idx)
    end
    Xs = [p[2] for p in parts]
    exposure = reduce(vcat, [fill(k, length(parts[k][1])) for k in eachindex(parts)]; init = Int[])
    group = [exps[k].group for k in exposure]
    names = feature_names(features)
    Dataset(isempty(Xs) ? zeros(0, length(names)) : reduce(vcat, Xs),
            reduce(vcat, [p[3] for p in parts]; init = Int[]),
            exposure, group, [p[1] for p in parts], [p[4] for p in parts], names)
end
