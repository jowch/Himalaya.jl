module HimalayaPeakModel

using Himalaya, Statistics, Random, Printf, SHA
import DecisionTree
using DecisionTree: RandomForestClassifier

include("matching.jl")
include("features.jl")
include("serial.jl")
include("featureset.jl")
include("learner.jl")
include("cache.jl")
include("dataset.jl")
include("recipe.jl")
include("cv.jl")

export LabeledExposure, Dataset, build_dataset, DiskFeatureCache, oof_proba, train, predict_peaks,
       cross_validate, print_cv, FEATURE_NAMES, DETECTOR_FEATURES,
       peak_match_tol, grid_step,
       Learner, Fitted, fit, predict_proba, RandomForestLearner, LogisticLearner,
       register_learner!, register_features!, to_dict, from_dict,
       FeatureSet, DefaultFeatures, feature_names, feature_version, candidates_and_features,
       Recipe, TrainedModel, default_recipe, recipe_id

end
