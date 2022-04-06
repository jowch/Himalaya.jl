module Himalaya

using Peaks
using Statistics

import Base: ==

export 

# indexing
Index, indexpeaks, predictpeaks, score, fit, npeak, ==,

# peaks
findpeaks

# include("util.jl")
include("peakfinding.jl")
include("index.jl")

end # module
