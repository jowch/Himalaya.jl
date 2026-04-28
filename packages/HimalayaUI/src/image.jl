using FileIO, TiffImages, ImageIO, ImageTransformations, ImageCore, ColorTypes

"""
    load_and_lognormalize(path) -> Matrix{Gray{Float32}}

Load a TIFF, convert to grayscale, apply log1p normalization to [0,1].
"""
function load_and_lognormalize(path::String)
    raw = load(path)

    # The TIFF stores photon counts as Q0f31 fixed-point (Int32 / 2^31).
    # Extracting raw Int32 counts gives the true photon numbers needed for
    # meaningful log compression. Zero and negative pixels are detector gaps,
    # beamstop shadow, or noise — they map to 0.0 (background).
    counts = Float32.(reinterpret.(Int32, channelview(raw)))

    lv = log1p.(max.(counts, 0f0))  # negatives → 0 before log

    # Clip at p99 of *positive* pixels only, so the direct beam / hot pixels
    # don't crush contrast in the diffraction rings.
    pos = filter(>(0f0), vec(lv))
    if isempty(pos)
        return colorview(Gray, lv)
    end
    hi = sort(pos)[min(end, round(Int, 0.99 * length(pos)))]

    normed = hi > 0f0 ? clamp.(lv ./ hi, 0f0, 1f0) : lv
    colorview(Gray, normed)
end

"""
    encode_png(img) -> Vector{UInt8}

Encode any grayscale image as PNG bytes.
"""
function encode_png(img)
    buf = IOBuffer()
    FileIO.save(FileIO.Stream{FileIO.format"PNG"}(buf), img)
    take!(buf)
end

"""
    resize_to_fit(img, max_px) -> image

Downscale so the longest side is ≤ max_px; no-op if already smaller.
"""
function resize_to_fit(img, max_px::Int)
    h, w = size(img)
    max(h, w) <= max_px && return img
    scale = max_px / max(h, w)
    imresize(img; ratio=scale)
end
