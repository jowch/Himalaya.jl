using FileIO, TiffImages, ImageIO, ImageTransformations, ImageCore, ColorTypes

"""
Bump this string whenever `load_and_lognormalize` semantics change. It feeds
into the per-image cache-busting token; bumping forces all browsers to re-fetch.
"""
const IMAGE_PROCESSING_VERSION = "v1"

"""
    image_version_token(path) -> String

Stable token combining `IMAGE_PROCESSING_VERSION` with the source TIFF's mtime.
The token is appended as `?v=<token>` to image URLs so the browser can cache
aggressively while still picking up changes when the source file is rewritten
or our image-processing code is bumped.

Returns `""` if `path` is `nothing`/`missing` or the file is missing.
"""
function image_version_token(path)
    (path === nothing || path isa Missing) && return ""
    p = String(path)
    isfile(p) || return ""
    string(IMAGE_PROCESSING_VERSION, "-", Int(round(mtime(p))))
end

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
