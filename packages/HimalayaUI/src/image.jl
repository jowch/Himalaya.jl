using FileIO, TiffImages, ImageIO, ImageTransformations, ImageCore, ColorTypes
using SQLite, DBInterface, Tables
using Base.Threads: @threads, Atomic, atomic_add!

const IMAGE_PROCESSING_VERSION = "v4"

const LOW_CLIP_PERCENTILE = 0.05f0
const HIGH_CLIP_PERCENTILE = 0.99f0

function image_version_token(path)
    (path === nothing || path isa Missing) && return ""
    p = String(path)
    isfile(p) || return ""
    string(IMAGE_PROCESSING_VERSION, "-", Int(round(mtime(p))))
end

# physics: the TIFF stores photon counts as Q0f31 fixed point (Int32 / 2^31),
# physics: so reinterpreting as Int32 recovers the counts. Zero and negative
# physics: pixels are detector gaps, beamstop shadow or noise.
_photon_counts(path::String) = max.(Float32.(reinterpret.(Int32, channelview(load(path)))), 0f0)

load_and_lognormalize(path::String) = colorview(Gray, _clip_normalize(log1p.(_photon_counts(path))))

# physics: clipping at p99 of positive pixels stops the direct beam and hot
# physics: pixels from crushing ring contrast. Clipping at p5 sends the faint
# physics: background floor to the same black as dead pixels.
function _clip_normalize(lv::AbstractMatrix{Float32})
    pos = sort!(filter(>(0f0), vec(lv)))
    isempty(pos) && return lv
    lo = pos[max(1, round(Int, LOW_CLIP_PERCENTILE * length(pos)))]
    hi = pos[min(end, round(Int, HIGH_CLIP_PERCENTILE * length(pos)))]
    # physics: on a beam-off frame almost every positive pixel is the same dark
    # physics: count, so the clip points coincide. Everything at or below that
    # physics: count renders black, and only brighter pixels render white.
    hi > lo || return Float32.(lv .> lo)
    clamp.((lv .- lo) ./ (hi - lo), 0f0, 1f0)
end

function encode_png(img)
    buf = IOBuffer()
    FileIO.save(FileIO.Stream{FileIO.format"PNG"}(buf), img)
    take!(buf)
end

function resize_to_fit(img, max_px::Int)
    h, w = size(img)
    max(h, w) <= max_px && return img
    scale = max_px / max(h, w)
    imresize(img; ratio=scale)
end

function load_and_lognormalize_thumb(path::String, max_px::Int = 128)
    # physics: the clip points come from the downscaled raster, so a faint narrow
    # physics: Bragg ring can under-render. Fine for a 128 px preview, never for
    # physics: the full-size view.
    img_small = resize_to_fit(colorview(Gray, _photon_counts(path)), max_px)
    lv = log1p.(Float32.(channelview(img_small)))
    colorview(Gray, _clip_normalize(lv))
end

function thumb_cache_dir(db::SQLite.DB)
    f = db.file
    (f === nothing || isempty(f) || f == ":memory:") && return nothing
    joinpath(dirname(f), "cache", "thumb-128")
end

function thumb_cache_path(db::SQLite.DB, exposure_id::Integer, token::AbstractString)
    dir = thumb_cache_dir(db)
    dir === nothing && return nothing
    joinpath(dir, string(exposure_id, "-", token, ".png"))
end

function _write_thumb_atomic(cache_path::String, bytes::Vector{UInt8})
    mkpath(dirname(cache_path))
    tmp, io = mktemp(dirname(cache_path))
    try
        write(io, bytes)
    finally
        close(io)
    end
    mv(tmp, cache_path; force = true)
    nothing
end

function ensure_thumb_cached(db::SQLite.DB, exposure_id::Integer, path::String;
                             token::AbstractString = image_version_token(path))
    cache_path = thumb_cache_path(db, exposure_id, token)

    if cache_path !== nothing
        hit = _read_cached(cache_path)
        hit === nothing || return hit
    end

    bytes = encode_png(load_and_lognormalize_thumb(path))
    cache_path === nothing || _write_thumb_atomic(cache_path, bytes)
    bytes
end

function full_cache_path(db::SQLite.DB, exposure_id::Integer, token::AbstractString)
    dir = thumb_cache_dir(db)
    dir === nothing && return nothing
    joinpath(dirname(dir), "full-1536", string(exposure_id), string(token, ".png"))
end

_full_size_path(cache_path::String) = first(splitext(cache_path)) * ".size"

function _read_cached(cache_path::String)
    f = try
        Base.Filesystem.open(cache_path, Base.Filesystem.JL_O_RDONLY)
    catch err
        err isa Base.IOError || err isa SystemError || rethrow()
        return nothing
    end
    try
        try
            touch(f)
        catch err
            err isa Base.IOError || err isa SystemError || rethrow()
        end
        read(f)
    finally
        close(f)
    end
end

function _read_full_cached(cache_path::String)
    try
        dims = tryparse.(Int, split(read(_full_size_path(cache_path), String)))
        (length(dims) == 2 && !any(isnothing, dims)) || return nothing
        bytes = _read_cached(cache_path)
        bytes === nothing && return nothing
        (bytes = bytes, raw_h = dims[1], raw_w = dims[2])
    catch err
        err isa SystemError || err isa Base.IOError || rethrow()
        nothing
    end
end

function _remove_stale_full!(cache_path::String)
    keep = (basename(cache_path), basename(_full_size_path(cache_path)))
    removed = 0
    for f in readdir(dirname(cache_path))
        (f in keep || !(endswith(f, ".png") || endswith(f, ".size"))) && continue
        removed += filesize(joinpath(dirname(cache_path), f))
        rm(joinpath(dirname(cache_path), f); force = true)
    end
    removed
end

const IMAGE_CACHE_CAP_BYTES = 1 << 30
const _IMAGE_CACHE_LOCK = ReentrantLock()
const _CACHE_BYTES = Dict{String, Int}()
const _ORPHAN_TEMP_AGE_S = 3600

function _cache_usage(root::String)
    total = 0
    entries = Tuple{Float64, Int, String}[]
    for (dir, _, files) in walkdir(root; onerror = _ -> nothing), f in files
        path = joinpath(dir, f)
        st = stat(path)
        if endswith(f, ".png")
            push!(entries, (st.mtime, st.size, path))
        elseif !endswith(f, ".size") && time() - st.mtime > _ORPHAN_TEMP_AGE_S
            rm(path; force = true)
            continue
        end
        total += st.size
    end
    total, entries
end

function _evict_cache!(root::String, keep::String, target::Integer)
    total, entries = _cache_usage(root)
    sort!(entries)
    for (_, bytes, png) in entries
        total <= target && break
        png == keep && continue
        size_path = _full_size_path(png)
        total -= bytes + filesize(size_path)
        rm(png; force = true)
        rm(size_path; force = true)
        basename(dirname(dirname(png))) == "full-1536" && isempty(readdir(dirname(png))) &&
            rm(dirname(png); force = true)
    end
    total
end

function _store_full!(cache_path::String, bytes::Vector{UInt8}, raw_h::Integer, raw_w::Integer, cap::Integer)
    root = dirname(dirname(dirname(cache_path)))
    size_bytes = Vector{UInt8}(string(raw_h, " ", raw_w))
    lock(_IMAGE_CACHE_LOCK) do
        total = get!(() -> first(_cache_usage(root)), _CACHE_BYTES, root)
        total -= filesize(cache_path) + filesize(_full_size_path(cache_path))
        _write_thumb_atomic(_full_size_path(cache_path), size_bytes)
        _write_thumb_atomic(cache_path, bytes)
        total += length(bytes) + length(size_bytes) - _remove_stale_full!(cache_path)
        total > cap && (total = _evict_cache!(root, cache_path, 9 * cap ÷ 10))
        _CACHE_BYTES[root] = total
    end
end

function ensure_full_cached(db::SQLite.DB, exposure_id::Integer, path::String;
                            token::AbstractString = image_version_token(path),
                            cap::Integer = IMAGE_CACHE_CAP_BYTES)
    cache_path = full_cache_path(db, exposure_id, token)
    if cache_path !== nothing
        hit = _read_full_cached(cache_path)
        hit === nothing || return hit
    end
    img = load_and_lognormalize(path)
    raw_h, raw_w = size(img)
    bytes = encode_png(resize_to_fit(img, 1536))
    cache_path === nothing || _store_full!(cache_path, bytes, raw_h, raw_w, cap)
    (bytes = bytes, raw_h = raw_h, raw_w = raw_w)
end

function prewarm_thumbnails!(db::SQLite.DB; threads::Bool = true,
                             experiment_id::Union{Integer, Nothing} = nothing)
    rows = experiment_id === nothing ?
        sqlrows(DBInterface.execute(db,
            "SELECT id, image_path FROM exposures WHERE image_path IS NOT NULL")) :
        sqlrows(DBInterface.execute(db,
            "SELECT id, image_path FROM exposures WHERE image_path IS NOT NULL AND experiment_id = ?",
            [Int(experiment_id)]))

    n_warmed  = Atomic{Int}(0)
    n_skipped = Atomic{Int}(0)

    work = function (row)
        path = String(row.image_path)
        if !isfile(path)
            atomic_add!(n_skipped, 1)
            @info "thumb prewarm skip" exposure_id=row.id reason="TIFF not on disk" path
            return
        end
        token = image_version_token(path)
        cached = thumb_cache_path(db, Int(row.id), token)
        if cached !== nothing && isfile(cached)
            atomic_add!(n_warmed, 1)
            return
        end
        try
            ensure_thumb_cached(db, Int(row.id), path; token)
        catch e
            atomic_add!(n_skipped, 1)
            @info "thumb prewarm skip" exposure_id=row.id reason="render failed" path exception=e
            return
        end
        atomic_add!(n_warmed, 1)
    end

    if threads
        @threads for row in rows
            work(row)
        end
    else
        for row in rows
            work(row)
        end
    end

    @info "thumb prewarm complete" warmed=n_warmed[] skipped=n_skipped[]
    (warmed = n_warmed[], skipped = n_skipped[])
end
