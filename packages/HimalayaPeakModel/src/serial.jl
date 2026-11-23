# Models are stored as plain dictionaries of JSON types. `from_dict` reads
# only those values and dispatches on a registered `kind`; it never evaluates
# or deserialises code.

"""
A learning algorithm and its hyperparameters. Implement
`fit(::MyLearner, X, y; rng) -> Fitted`, `to_dict`, `from_dict(::Type{Learner}, d)`
via `register_learner!`, and `seed`.
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
and `TrainedModel`. `d` may have `String` or `Symbol` keys (so `JSON3.read`
output works). An unregistered `kind` is an error.
"""
function from_dict end

fld(d, k::AbstractString) = haskey(d, Symbol(k)) ? d[Symbol(k)] :
    haskey(d, k) ? d[k] : error("model data has no field \"$k\"")

const LEARNERS = Dict{String, Function}()
const FITTED = Dict{String, Function}()
const FEATURESETS = Dict{String, Function}()

"""
    register_learner!(kind, learner_loader, fitted_loader)

Makes a learner selectable by recipes and loadable from stored models. Each
loader takes the dictionary and returns the object. Loaders must not execute
anything taken from the data.
"""
function register_learner!(kind::AbstractString, learner_loader, fitted_loader)
    LEARNERS[kind] = learner_loader
    FITTED[kind] = fitted_loader
    nothing
end

register_features!(kind::AbstractString, loader) = (FEATURESETS[kind] = loader; nothing)

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
