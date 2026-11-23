function fit end

"""Probability of a real peak for each row of `X`."""
function predict_proba end

"""Seed for the training RNG. Defaults to the learner's `seed` field; override for a learner without one."""
seed(l::Learner) = l.seed

"""
    RandomForestLearner(; n_trees = 200, min_samples_leaf = 5, seed = 1,
                        n_subfeatures = -1, partial_sampling = 0.7, max_depth = -1,
                        min_samples_split = 2)

DecisionTree.jl random forest.
"""
Base.@kwdef struct RandomForestLearner <: Learner
    n_trees::Int = 200
    min_samples_leaf::Int = 5
    seed::Int = 1
    n_subfeatures::Int = -1
    partial_sampling::Float64 = 0.7
    max_depth::Int = -1
    min_samples_split::Int = 2
end

struct ForestFitted <: Fitted
    forest::RandomForestClassifier
end

function fit(l::RandomForestLearner, X::AbstractMatrix, y::AbstractVector{Bool}; rng = MersenneTwister(l.seed))
    forest = RandomForestClassifier(n_trees = l.n_trees, min_samples_leaf = l.min_samples_leaf,
        n_subfeatures = l.n_subfeatures, partial_sampling = l.partial_sampling,
        max_depth = l.max_depth, min_samples_split = l.min_samples_split, rng = rng)
    DecisionTree.fit!(forest, X, Int.(y))
    ForestFitted(forest)
end

function predict_proba(f::ForestFitted, X::AbstractMatrix)
    isempty(X) && return Float64[]
    P = DecisionTree.predict_proba(f.forest, X)
    j = findfirst(==(1), f.forest.classes)
    j === nothing ? zeros(size(X, 1)) : Float64.(P[:, j])
end

# A leaf stores its majority class and the training counts per class. A forest
# prediction uses only the majority (one vote per tree in DecisionTree 0.12);
# the counts rebuild `Leaf.values`, so per-tree class fractions match too.
function tree_dict(l::DecisionTree.Leaf)
    n = Dict{String, Any}("0" => count(==(0), l.values), "1" => count(==(1), l.values))
    Dict{String, Any}("c" => Int(l.majority), "n" => n)
end
tree_dict(n::DecisionTree.Node) = Dict{String, Any}("f" => n.featid, "t" => Float64(n.featval),
    "l" => tree_dict(n.left), "r" => tree_dict(n.right))

function tree_from_dict(d)
    if haskey(d, :c) || haskey(d, "c")
        n = fld(d, "n")
        values = Int[fill(0, Int(fld(n, "0"))); fill(1, Int(fld(n, "1")))]
        DecisionTree.Leaf{Int}(Int(fld(d, "c")), values)
    else
        DecisionTree.Node{Float64, Int}(Int(fld(d, "f")), Float64(fld(d, "t")),
            tree_from_dict(fld(d, "l")), tree_from_dict(fld(d, "r")))
    end
end

learner_params(l::RandomForestLearner) = Dict{String, Any}("n_trees" => l.n_trees,
    "min_samples_leaf" => l.min_samples_leaf, "seed" => l.seed, "n_subfeatures" => l.n_subfeatures,
    "partial_sampling" => l.partial_sampling, "max_depth" => l.max_depth,
    "min_samples_split" => l.min_samples_split)

to_dict(l::RandomForestLearner) = merge(Dict{String, Any}("kind" => "random_forest", "version" => 1),
                                        learner_params(l))

function rf_learner_from_dict(d)
    check_version(d, 1, "random_forest learner")
    RandomForestLearner(n_trees = Int(fld(d, "n_trees")), min_samples_leaf = Int(fld(d, "min_samples_leaf")),
        seed = Int(fld(d, "seed")), n_subfeatures = Int(fld(d, "n_subfeatures")),
        partial_sampling = Float64(fld(d, "partial_sampling")), max_depth = Int(fld(d, "max_depth")),
        min_samples_split = Int(fld(d, "min_samples_split")))
end

function to_dict(f::ForestFitted)
    ens = f.forest.ensemble
    Dict{String, Any}("kind" => "random_forest", "version" => 2, "n_features" => ens.n_feat,
        "classes" => Int.(f.forest.classes), "trees" => Any[tree_dict(t) for t in ens.trees])
end

function rf_fitted_from_dict(d)
    check_version(d, 2, "random_forest model")
    n = Int(fld(d, "n_features"))
    trees = DecisionTree.LeafOrNode{Float64, Int}[tree_from_dict(t) for t in fld(d, "trees")]
    classes = Int[Int(c) for c in fld(d, "classes")]
    ens = DecisionTree.Ensemble{Float64, Int}(trees, n, zeros(n))
    ForestFitted(RandomForestClassifier(n_trees = length(trees), ensemble = ens, classes = classes))
end

"""
    LogisticLearner(; steps = 300, lr = 0.5, l2 = 1e-3, seed = 1)

L2-regularised logistic regression on standardised features, fitted by
full-batch gradient descent. Deterministic; `seed` is unused.
"""
Base.@kwdef struct LogisticLearner <: Learner
    steps::Int = 300
    lr::Float64 = 0.5
    l2::Float64 = 1e-3
    seed::Int = 1
end

struct LogisticFitted <: Fitted
    mean::Vector{Float64}
    scale::Vector{Float64}
    w::Vector{Float64}
    b::Float64
end

sigmoid(z) = 1 / (1 + exp(-clamp(z, -35.0, 35.0)))

function fit(l::LogisticLearner, X::AbstractMatrix, y::AbstractVector{Bool}; rng = nothing)
    n, p = size(X)
    μ = vec(mean(X; dims = 1))
    s = [(v = std(view(X, :, j))) > 0 && isfinite(v) ? v : 1.0 for j in 1:p]
    Z = (X .- μ') ./ s'
    w = zeros(p); b = log(max(count(y), 1) / max(count(!, y), 1))
    t = Float64.(y)
    for _ in 1:l.steps
        r = sigmoid.(Z * w .+ b) .- t
        w .-= l.lr .* (Z' * r ./ n .+ l.l2 .* w)
        b -= l.lr * sum(r) / n
    end
    (all(isfinite, w) && isfinite(b)) || error("logistic fit diverged (non-finite weights); lower lr")
    LogisticFitted(μ, s, w, b)
end

function predict_proba(f::LogisticFitted, X::AbstractMatrix)
    isempty(X) && return Float64[]
    sigmoid.(((X .- f.mean') ./ f.scale') * f.w .+ f.b)
end

to_dict(l::LogisticLearner) = Dict{String, Any}("kind" => "logistic", "version" => 1,
    "steps" => l.steps, "lr" => l.lr, "l2" => l.l2, "seed" => l.seed)

function logistic_learner_from_dict(d)
    check_version(d, 1, "logistic learner")
    LogisticLearner(steps = Int(fld(d, "steps")), lr = Float64(fld(d, "lr")),
                    l2 = Float64(fld(d, "l2")), seed = Int(fld(d, "seed")))
end

to_dict(f::LogisticFitted) = Dict{String, Any}("kind" => "logistic", "version" => 1,
    "mean" => f.mean, "scale" => f.scale, "w" => f.w, "b" => f.b)

function logistic_fitted_from_dict(d)
    check_version(d, 1, "logistic model")
    LogisticFitted(floats(fld(d, "mean")), floats(fld(d, "scale")), floats(fld(d, "w")), Float64(fld(d, "b")))
end

register_learner!("random_forest", rf_learner_from_dict, rf_fitted_from_dict)
register_learner!("logistic", logistic_learner_from_dict, logistic_fitted_from_dict)
