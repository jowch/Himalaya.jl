module Himalaya

using Peaks
using Statistics
using Images, TiffImages
using DelimitedFiles

import Base: ==

export 

# indexing
Index, indexpeaks, score, fit, npeak, ==,

# peaks
findpeaks,

# images
load_image,

# traces
load_trace

# include("util.jl")
include("peakfinding.jl")
include("index.jl")
include("image.jl")
include("trace.jl")

end # module
