
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
