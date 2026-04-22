# Example utilities for loading SAXS detector images and integration files.
#
# These helpers are NOT part of the Himalaya package — they require extra
# dependencies that are not in Himalaya's Project.toml:
#   - DelimitedFiles  (stdlib, free)
#   - Images, TiffImages
#
# To use, install those packages in your own environment and `include` this file.

using DelimitedFiles
using Images, TiffImages

function load_image(path; nbins = 256, kwargs...)
    img = reinterpret(Gray{N0f32}, TiffImages.load(path))
    
    # take log of intensities
    a = log.(Array{Float32}(img) .+ eps())
    # normalize log intensities to 0-1 range
    a_unit = (a .- minimum(a)) ./ (maximum(a) - minimum(a))
    # improve contrast and make visualization easier to understand
    # a_eq = adjust_histogram(a_unit, Images.Equalization())
    a_eq = clahe(a_unit, nbins; kwargs...)
    
    # clamp noisy float digits and return
    Gray.(convert.(N0f32, clamp.(a_eq, 0, 1)))
end

function load_trace(path)
    readdlm(path, ' ', Float64, '\n')
end
