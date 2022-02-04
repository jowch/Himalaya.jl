module Himalaya

using Peaks
using Statistics
using Images, TiffImages

export 

# indexing
Index, indexpeaks, score, fit, npeak,

# peaks
findpeaks,

# images
load_image


# include("util.jl")
include("peakfinding.jl")
include("index.jl")
include("image.jl")

end # module
