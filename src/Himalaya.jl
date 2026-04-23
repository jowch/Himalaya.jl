module Himalaya

using LinearAlgebra
using Peaks
using SparseArrays
using Statistics

import Base: ==, issubset, show

export

# phases
Phase, Lamellar, Hexagonal, Square, Pn3m, Im3m, Ia3d, Fm3m, Fd3m, phaseratios,

# indexing
Index, phase, basis, peaks, numpeaks, predictpeaks, missingpeaks,
==, issubset, show, indexpeaks, fit, score,

# peaks
findpeaks, persistence, knee

include("util.jl")
include("phase.jl")
include("threshold.jl")
include("persistence.jl")
include("sharpness.jl")
include("peakfinding.jl")
include("index.jl")

end # module
