"""
A learning algorithm and its hyperparameters. Implement
`fit(::MyLearner, X, y; rng) -> Fitted` and `to_dict`, give it an `Int` field
`seed` (or define `HimalayaPeakModel.seed(::MyLearner)`), and make it loadable with
`register_learner!`.
"""
abstract type Learner end

"""A trained classifier. Implement `predict_proba(f, X)` and `to_dict`."""
abstract type Fitted end

"""
Candidate generation plus per-candidate features. Implement `feature_names`,
`feature_version`, `candidates_and_features`, `to_dict` and register a loader
with `register_features!`.
"""
abstract type FeatureSet end

"""
    to_dict(x) -> Dict{String, Any}

JSON-safe form of a `Learner`, `Fitted`, `FeatureSet`, `Recipe` or
`TrainedModel`. Learners, fitted models and feature sets carry `"kind"` and
`"version"`.
"""
function to_dict end

"""
    from_dict(T, d)

Inverse of `to_dict` for `T` in `Learner`, `Fitted`, `FeatureSet`, `Recipe`
and `TrainedModel`. `d` may have `String` or `Symbol` keys (so `JSON.parse`
output works). An unregistered `kind` is an error.
"""
function from_dict end

fld(d, k::AbstractString) = haskey(d, Symbol(k)) ? d[Symbol(k)] :
    haskey(d, k) ? d[k] : error("model data has no field \"$k\"")

const LEARNERS = Dict{String, Function}()
const FITTED = Dict{String, Function}()
const FEATURESETS = Dict{String, Function}()

"""
    register_learner!(kind, learner_loader, fitted_loader; replace = false)

Makes a learner selectable by recipes and loadable from stored models. Each
loader takes the dictionary and returns the object. Loaders must not execute
anything taken from the data. Registering a `kind` twice is an error unless
`replace = true`.

An external package must call this from its `__init__`: a top-level call runs
during that package's precompilation, and the change to this module's registry
is not saved.
"""
function register_learner!(kind::AbstractString, learner_loader, fitted_loader; replace::Bool = false)
    replace || !(haskey(LEARNERS, kind) || haskey(FITTED, kind)) || error("learner kind \"$kind\" is already registered")
    LEARNERS[kind] = learner_loader
    FITTED[kind] = fitted_loader
    nothing
end

"""`register_features!(kind, loader; replace = false)`; same rules as `register_learner!`."""
function register_features!(kind::AbstractString, loader; replace::Bool = false)
    replace || !haskey(FEATURESETS, kind) || error("feature set kind \"$kind\" is already registered")
    FEATURESETS[kind] = loader
    nothing
end

function load_kind(registry, what, d)
    kind = String(fld(d, "kind"))
    haskey(registry, kind) || error("unknown $what kind \"$kind\"; known: $(join(sort!(collect(keys(registry))), ", "))")
    registry[kind](d)
end

from_dict(::Type{Learner}, d) = load_kind(LEARNERS, "learner", d)
from_dict(::Type{Fitted}, d) = load_kind(FITTED, "fitted model", d)
from_dict(::Type{FeatureSet}, d) = load_kind(FEATURESETS, "feature set", d)

function check_version(d, v, what)
    got = Int(fld(d, "version"))
    got == v || error("unsupported $what version $got (expected $v); retrain it")
end

floats(v) = Float64[Float64(x) for x in v]
