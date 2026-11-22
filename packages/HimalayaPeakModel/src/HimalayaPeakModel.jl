module HimalayaPeakModel

using Himalaya, Statistics, Random, Printf
using DecisionTree

include("matching.jl")
include("features.jl")
include("dataset.jl")
include("model.jl")
include("cv.jl")

export LabeledExposure, Dataset, PeakModel, build_dataset, train, predict_peaks,
       cross_validate, print_cv, FEATURE_NAMES, DETECTOR_FEATURES,
       peak_match_tol, grid_step

end
