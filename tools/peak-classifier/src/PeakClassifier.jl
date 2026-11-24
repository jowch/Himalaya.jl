module PeakClassifier

import JSON
using HimalayaPeakModel
import HimalayaUI
import Himalaya
import HimalayaPeakModel as HPM

include("modelio.jl")
include("db.jl")
include("evaluate.jl")

export open_readonly, load_labeled_exposures, save_model, load_model,
       LabeledExposure, build_dataset, train, predict_peaks, cross_validate, print_cv, FEATURE_NAMES,
       Recipe, TrainedModel, DefaultFeatures, RandomForestLearner, LogisticLearner, default_recipe,
       evaluate_database, print_evaluation, write_evaluation, phase_table

end
