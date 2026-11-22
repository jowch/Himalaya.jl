"""
Trained peak classifier. `threshold` and `nms` (points) are applied at
prediction time; they are stored with the model so a saved file predicts the
way it was validated.
"""
struct PeakModel
    forest::RandomForestClassifier
    feature_names::Vector{String}
    threshold::Float64
    nms::Int
    n_train::NamedTuple{(:positive, :negative), Tuple{Int, Int}}
end

"""
    train(ds::Dataset; rows = nothing, n_trees = 200, min_samples_leaf = 5,
          threshold = 0.5, nms = 8, weak_keep = 0.5, seed = 1) -> PeakModel

Trains on positives and negatives. Weak negatives are subsampled with
probability `weak_keep` (DecisionTree.jl has no sample weights). `rows`
restricts training to a subset of candidate rows.
"""
function train(ds::Dataset; rows = nothing, n_trees::Int = 200, min_samples_leaf::Int = 5,
               threshold::Real = 0.5, nms::Int = 8, weak_keep::Real = 0.5, seed::Int = 1)
    rng = MersenneTwister(seed)
    rows === nothing && (rows = 1:length(ds.label))
    sel = Int[]
    for r in rows
        l = ds.label[r]
        (l == POS || l == NEG || (l == WEAK && rand(rng) < weak_keep)) && push!(sel, r)
    end
    y = ds.label[sel] .== POS
    (any(y) && !all(y)) || error("training needs both positive and negative candidates ($(count(y)) positive of $(length(y)))")
    forest = RandomForestClassifier(n_trees = n_trees, min_samples_leaf = min_samples_leaf,
                                    rng = rng)
    DecisionTree.fit!(forest, ds.X[sel, :], Int.(y))
    PeakModel(forest, FEATURE_NAMES, Float64(threshold), nms,
              (positive = count(y), negative = count(!, y)))
end

function proba(m::PeakModel, X::AbstractMatrix)
    isempty(X) && return Float64[]
    P = DecisionTree.predict_proba(m.forest, X)
    j = findfirst(==(1), m.forest.classes)
    j === nothing ? zeros(size(X, 1)) : P[:, j]
end

"""
Indices (into `idx`) of candidates with probability ≥ threshold, thinned so
kept grid points are more than `nms` apart, strongest first.
"""
function select_peaks(idx, p, threshold, nms)
    order = sort([k for k in eachindex(idx) if p[k] >= threshold]; by = k -> -p[k])
    keep = Int[]
    for k in order
        all(j -> abs(idx[k] - idx[j]) > nms, keep) && push!(keep, k)
    end
    sort!(keep)
end

"""
    predict_peaks(model, q, I, σ = nothing) -> Vector{Float64}

q values of the peaks the model accepts. Without `σ` it is estimated from the
trace as the hybrid detector does.
"""
function predict_peaks(m::PeakModel, q, I, σ = nothing)
    σ === nothing && (σ = Himalaya._estimate_sigma(I))
    idx, X = candidate_features(q, I, σ)
    isempty(idx) && return Float64[]
    p = proba(m, X)
    Float64.(q[idx[select_peaks(idx, p, m.threshold, m.nms)]])
end

function save_model(path::AbstractString, m::PeakModel)
    mkpath(dirname(abspath(path)))
    serialize(path, (version = 1, julia = string(VERSION), model = m))
    path
end

function load_model(path::AbstractString)
    d = deserialize(path)
    d.version == 1 || error("unsupported model file version $(d.version)")
    d.model.feature_names == FEATURE_NAMES ||
        error("model was trained with a different feature set; retrain it")
    d.model
end
