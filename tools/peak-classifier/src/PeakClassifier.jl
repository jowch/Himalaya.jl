module PeakClassifier

using Himalaya, Statistics, Random, Printf, Serialization
using DecisionTree

include("features.jl")
include("dataset.jl")
include("model.jl")
include("cv.jl")
include("db.jl")

export LabeledExposure, build_dataset, train, predict_peaks, save_model, load_model,
       cross_validate, load_labeled_exposures, FEATURE_NAMES

end
