
function load_image(path)
    img = reinterpret(Gray{N0f32}, TiffImages.load(target_path))
    
    # take log of intensities
    a = log.(Array{Float32}(img) .+ eps())
    # normalize log intensities to 0-1 range
    a_unit = (a .- minimum(a)) ./ (maximum(a) - minimum(a))
    # improve contrast and make visualization easier to understand
    a_eq = adjust_histogram(a_unit, Images.Equalization())
    
    # clamp noisy float digits and return
    Gray.(convert.(N0f32, clamp.(a_eq, 0, 1)))
end


function load_image(image_dir, fname)
    image_paths = filter(contains(".tif"), readdir(image_dir, join=true))
    target_path = image_paths[findfirst(contains(fname), image_paths)]

    load_image(target_path)
end
