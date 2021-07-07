module Himalaya

using Peaks
using Statistics
using Images, TiffImages

export 

# indexing
Index, index_peaks, score, fit,

# peaks
find_peaks,

# images
load_image


include("util.jl")
include("peakfinding.jl")
include("index.jl")
include("image.jl")

end # module
