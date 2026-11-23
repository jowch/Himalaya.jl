function feature_names end

"""
Identifies the feature definitions and the code they depend on. Cached features
and trained models are valid only for an equal version string.
"""
function feature_version end

"""
    candidates_and_features(fs, q, I, σ; hybrid = nothing) -> (idx, X, hybrid)

Candidate grid indices, the `length(idx) × length(feature_names(fs))` matrix,
and the hybrid detector's result. Pass a precomputed `hybrid` result to avoid
running the detector again.
"""
function candidates_and_features end

"""
    DefaultFeatures(detector_features = true)

The matched-filter, persistence and context features of `candidate_features`.
`detector_features = false` drops the columns in `DETECTOR_FEATURES`.
"""
struct DefaultFeatures <: FeatureSet
    detector_features::Bool
end
DefaultFeatures() = DefaultFeatures(true)

feature_names(fs::DefaultFeatures) =
    fs.detector_features ? FEATURE_NAMES : filter(!in(DETECTOR_FEATURES), FEATURE_NAMES)

feature_version(fs::DefaultFeatures) =
    "default-1;detector_features=$(fs.detector_features);himalaya=$(pkgversion(Himalaya))"

function candidates_and_features(fs::DefaultFeatures, q, I, σ; hybrid = nothing)
    idx, X, hyb = candidate_features(q, I, σ; hybrid)
    fs.detector_features ? (idx, X, hyb) : (idx, X[:, indexin(feature_names(fs), FEATURE_NAMES)], hyb)
end

to_dict(fs::DefaultFeatures) = Dict{String, Any}("kind" => "default", "version" => 1,
                                                 "detector_features" => fs.detector_features)

function default_features_from_dict(d)
    check_version(d, 1, "default feature set")
    DefaultFeatures(Bool(fld(d, "detector_features")))
end

register_features!("default", default_features_from_dict)
from_dict(::Type{DefaultFeatures}, d) = default_features_from_dict(d)
