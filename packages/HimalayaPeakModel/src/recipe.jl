"""
Everything that defines how a model is trained and applied: feature set,
learner, label options and the decision rule. Two models trained from the same
recipe on the same data are the same model; `recipe_id` names the recipe.

- `weak_keep`: probability that a weak negative is kept for training.
- `threshold`, `nms`: a candidate is accepted at probability ≥ `threshold`,
  then thinned so kept grid points are more than `nms` apart.
"""
struct Recipe
    features::FeatureSet
    learner::Learner
    weak_keep::Float64
    threshold::Float64
    nms::Int
end

Recipe(features::FeatureSet, learner::Learner; weak_keep::Real = 0.5, threshold::Real = 0.5, nms::Int = 8) =
    Recipe(features, learner, Float64(weak_keep), Float64(threshold), nms)

"""
    default_recipe(; n_trees = 200, min_samples_leaf = 5, threshold = 0.5, nms = 8,
                   weak_keep = 0.5, seed = 1, detector_features = true) -> Recipe

`DefaultFeatures` with a `RandomForestLearner`.
"""
default_recipe(; n_trees::Int = 200, min_samples_leaf::Int = 5, threshold::Real = 0.5, nms::Int = 8,
               weak_keep::Real = 0.5, seed::Int = 1, detector_features::Bool = true) =
    Recipe(DefaultFeatures(detector_features),
           RandomForestLearner(; n_trees, min_samples_leaf, seed); weak_keep, threshold, nms)

to_dict(r::Recipe) = Dict{String, Any}("version" => 1, "features" => to_dict(r.features),
    "learner" => to_dict(r.learner), "weak_keep" => r.weak_keep, "threshold" => r.threshold, "nms" => r.nms)

function from_dict(::Type{Recipe}, d)
    check_version(d, 1, "recipe")
    Recipe(from_dict(FeatureSet, fld(d, "features")), from_dict(Learner, fld(d, "learner")),
           Float64(fld(d, "weak_keep")), Float64(fld(d, "threshold")), Int(fld(d, "nms")))
end

# Canonical text: sorted keys, shortest round-trip floats; independent of Dict order.
function canonical(x)
    x isa AbstractDict ? "{" * join(("$(repr(String(k))):$(canonical(x[k]))" for k in sort!(String.(collect(keys(x))))), ",") * "}" :
    x isa AbstractVector ? "[" * join(canonical.(x), ",") * "]" :
    x isa AbstractString ? repr(String(x)) :
    x isa Integer || x isa Bool ? string(Int(x)) :
    x isa Real ? repr(Float64(x)) : error("cannot hash $(typeof(x))")
end

"""Stable 16-hex-digit identifier of a recipe (SHA-1 of its canonical dict)."""
recipe_id(r::Recipe) = bytes2hex(sha1(canonical(to_dict(r))))[1:16]

"""
A fitted model with the recipe that produced it. `feature_names` are the
columns `fitted` expects. `trained_on` counts the training candidates.
"""
struct TrainedModel
    recipe::Recipe
    fitted::Fitted
    threshold::Float64
    nms::Int
    feature_names::Vector{String}
    trained_on::NamedTuple{(:positive, :negative), Tuple{Int, Int}}
end

"""
    train(recipe, ds::Dataset; rows = nothing) -> TrainedModel

Trains on positives and negatives. Weak negatives are subsampled with
probability `recipe.weak_keep` (DecisionTree.jl has no sample weights). `rows`
restricts training to a subset of candidate rows. `ds` must have been built
with a feature set that has every column the recipe's feature set names.
"""
function train(recipe::Recipe, ds::Dataset; rows = nothing)
    rng = MersenneTwister(seed(recipe.learner))
    rows === nothing && (rows = 1:length(ds.label))
    sel = Int[]
    for r in rows
        l = ds.label[r]
        (l == POS || l == NEG || (l == WEAK && rand(rng) < recipe.weak_keep)) && push!(sel, r)
    end
    y = ds.label[sel] .== POS
    (any(y) && !all(y)) || error("training needs both positive and negative candidates ($(count(y)) positive of $(length(y)))")
    names = feature_names(recipe.features)
    cols = indexin(names, ds.feature_names)
    any(isnothing, cols) && error("dataset lacks features $(names[isnothing.(cols)]); build it with the recipe's feature set")
    fitted = fit(recipe.learner, ds.X[sel, cols], y; rng)
    TrainedModel(recipe, fitted, recipe.threshold, recipe.nms, names,
                 (positive = count(y), negative = count(!, y)))
end

"""
    train(ds::Dataset; rows = nothing, kwargs...) -> TrainedModel

Trains `default_recipe(; kwargs...)`.
"""
train(ds::Dataset; rows = nothing, kwargs...) = train(default_recipe(; kwargs...), ds; rows)

"""Probabilities for the `rows` of `ds` (all columns of `ds` available, the model takes its own)."""
function proba(m::TrainedModel, ds::Dataset, rows = 1:length(ds.label))
    isempty(rows) && return Float64[]
    predict_proba(m.fitted, ds.X[rows, indexin(m.feature_names, ds.feature_names)])
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
    predict_peaks(model, q, I, σ = nothing; hybrid = nothing) -> Vector{Float64}

q values of the peaks the model accepts. Pass the file's `σ` when you have it:
training uses it. Without `σ` it is estimated from the trace as the hybrid
detector does, and the σ-dependent features then differ from training, so
predictions can shift. Pass `hybrid` (a `findpeaks(...; method = :hybrid)`
result for this trace) to skip running the detector again.
"""
function predict_peaks(m::TrainedModel, q, I, σ = nothing; hybrid = nothing)
    σ === nothing && (σ = Himalaya._estimate_sigma(I))
    idx, X, _ = candidates_and_features(m.recipe.features, q, I, σ; hybrid)
    isempty(idx) && return Float64[]
    p = predict_proba(m.fitted, X)
    Float64.(q[idx[select_peaks(idx, p, m.threshold, m.nms)]])
end

to_dict(m::TrainedModel) = Dict{String, Any}("version" => 1, "recipe" => to_dict(m.recipe),
    "recipe_id" => recipe_id(m.recipe), "feature_version" => feature_version(m.recipe.features),
    "feature_names" => m.feature_names, "threshold" => m.threshold, "nms" => m.nms,
    "trained_on" => Dict{String, Any}("positive" => m.trained_on.positive, "negative" => m.trained_on.negative),
    "fitted" => to_dict(m.fitted))

function from_dict(::Type{TrainedModel}, d)
    check_version(d, 1, "model")
    recipe = from_dict(Recipe, fld(d, "recipe"))
    names = String[String(x) for x in fld(d, "feature_names")]
    names == feature_names(recipe.features) ||
        error("model was trained with different features than its recipe names; retrain it")
    t = fld(d, "trained_on")
    TrainedModel(recipe, from_dict(Fitted, fld(d, "fitted")), Float64(fld(d, "threshold")),
                 Int(fld(d, "nms")), names,
                 (positive = Int(fld(t, "positive")), negative = Int(fld(t, "negative"))))
end
