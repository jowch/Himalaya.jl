module Himalaya

using Peaks
using Statistics

import Base: ==, issubset, show

export 

# phases
Phase, Lamellar, Hexagonal, Pn3m, Im3m, Ia3d, Fd3m, phaseratios,

# indexing
Index, phase, basis, peaks, numpeaks, predictpeaks, ==, issubset, show,
indexpeaks, fit, score,

# peaks
findpeaks

# include("util.jl")
include("phase.jl")
include("peakfinding.jl")
include("index.jl")

end # module
