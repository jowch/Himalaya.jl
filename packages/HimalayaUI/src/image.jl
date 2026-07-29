using FileIO, TiffImages, ImageIO, ImageTransformations, ImageCore, ColorTypes
using SQLite, DBInterface, Tables
using Base.Threads: @threads, Atomic, atomic_add!

"""
Bump this string whenever `load_and_lognormalize` semantics change. It feeds
into the per-image cache-busting token; bumping forces all browsers to re-fetch.
"""
const IMAGE_PROCESSING_VERSION = "v3"

"""
Percentile clips applied to the log-counts of *positive* pixels in
`load_and_lognormalize`. Bump `IMAGE_PROCESSING_VERSION` whenever these
change so cached PNGs are invalidated.
"""
const LOW_CLIP_PERCENTILE = 0.05f0
const HIGH_CLIP_PERCENTILE = 0.99f0

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
    # don't crush contrast in the diffraction rings. Also clip the low end at
    # p5 of positives so the faint background floor of real signal collapses
    # to the same pure black as detector dead pixels — otherwise the whole
    # frame sits on a slightly-lifted gray haze that hides ring contrast.
    pos = sort!(filter(>(0f0), vec(lv)))
    if isempty(pos)
        return colorview(Gray, lv)
    end
    lo = pos[max(1, round(Int, LOW_CLIP_PERCENTILE * length(pos)))]
    hi = pos[min(end, round(Int, HIGH_CLIP_PERCENTILE * length(pos)))]

    normed = hi > lo ? clamp.((lv .- lo) ./ (hi - lo), 0f0, 1f0) : lv
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

"""
    load_and_lognormalize_thumb(path, max_px = 128) -> Matrix{Gray{Float32}}

Triage-preview variant of [`load_and_lognormalize`](@ref) for thumbnails: the
loaded TIFF is downscaled to `max_px` max-side *before* the log1p +
percentile-clip + normalize runs, so the percentile math sees ~16K pixels
instead of the full ~4M (issue #261, U-5 — an order-of-magnitude per-request
speedup on every cold thumb fetch).

Because the percentile clip points are computed on the *downsampled* raster,
they differ slightly from the full-resolution version: this is a triage
preview where a faint, narrow Bragg ring may under-render relative to the
full-quality image. For a 128px contact-sheet thumbnail that trade is
perceptually invisible; never use this for the science-quality view (the
`size="full"` route path keeps `load_and_lognormalize`).
"""
function load_and_lognormalize_thumb(path::String, max_px::Int = 128)
    raw = load(path)

    # Recover raw Q0f31 photon counts exactly as load_and_lognormalize does;
    # negatives (beamstop/dead pixels) clip to 0.
    counts = Float32.(reinterpret.(Int32, channelview(raw)))
    img_full = colorview(Gray, max.(counts, 0f0))

    # Downscale FIRST — the whole point of the thumb variant.
    img_small = resize_to_fit(img_full, max_px)

    lv = log1p.(Float32.(channelview(img_small)))
    pos = sort!(filter(>(0f0), vec(lv)))
    isempty(pos) && return colorview(Gray, lv)
    lo = pos[max(1, round(Int, LOW_CLIP_PERCENTILE * length(pos)))]
    hi = pos[min(end, round(Int, HIGH_CLIP_PERCENTILE * length(pos)))]
    normed = hi > lo ? clamp.((lv .- lo) ./ (hi - lo), 0f0, 1f0) : lv
    colorview(Gray, normed)
end

"""
    thumb_cache_dir(db) -> Union{String, Nothing}

The on-disk thumbnail cache directory for the corpus backing `db`:
`<db_dir>/cache/thumb-128/`, co-located with the SQLite file so it travels with
the corpus. Returns `nothing` for an in-memory (`:memory:`) DB — those have no
stable directory, so callers fall back to rendering fresh every request.
"""
function thumb_cache_dir(db::SQLite.DB)
    f = db.file
    (f === nothing || isempty(f) || f == ":memory:") && return nothing
    joinpath(dirname(f), "cache", "thumb-128")
end

"""
    thumb_cache_path(db, exposure_id, token) -> Union{String, Nothing}

Cache file path for one exposure's thumbnail. The filename embeds the
`image_version_token` (which already folds in `IMAGE_PROCESSING_VERSION` + the
TIFF mtime), so a TIFF rewrite or a processing-code bump yields a new key and
the stale entry is simply never read again. `nothing` for `:memory:` DBs.
"""
function thumb_cache_path(db::SQLite.DB, exposure_id::Integer, token::AbstractString)
    dir = thumb_cache_dir(db)
    dir === nothing && return nothing
    joinpath(dir, string(exposure_id, "-", token, ".png"))
end

"""
    _write_thumb_atomic(cache_path, bytes)

Write `bytes` to `cache_path` via a temp file + `mv` so a concurrent reader (or
a parallel prewarm thread) never observes a half-written PNG. `force=true` makes
the rename overwrite an existing entry — used by the unconditional reingest
re-warm (issue #261, must-fix: whole-second mtime granularity would otherwise
leave a same-second re-ingest serving a stale cached thumb).
"""
function _write_thumb_atomic(cache_path::String, bytes::Vector{UInt8})
    mkpath(dirname(cache_path))
    tmp = string(cache_path, ".", getpid(), ".", Threads.threadid(), ".tmp")
    write(tmp, bytes)
    mv(tmp, cache_path; force = true)
    nothing
end

"""
    ensure_thumb_cached(db, exposure_id, path; token = image_version_token(path), overwrite = false) -> Vector{UInt8}

Return the thumbnail PNG bytes for `exposure_id`, reading the on-disk cache when
present and rendering + persisting on a miss. On a `:memory:` DB (no cache dir)
the thumbnail is rendered fresh and returned without a write.

The cache key is `token` — by default `image_version_token(path)` (folds in
`IMAGE_PROCESSING_VERSION` + the source TIFF's mtime). A cache HIT returns the
stored bytes via a single `read`, NEVER re-touching the source TIFF — so the
thumbnail still serves after the source is moved/deleted. Note callers MUST pass
the token (or call while the source exists), because `image_version_token`
returns `""` for a missing file: recomputing it lazily on a hit-after-delete
would compute the WRONG key and miss the cached entry. The route and prewarm
both compute the token while the guaranteed-present source is still on disk.

`overwrite = true` forces a re-render + atomic replace even on a cache hit; the
reingest prewarm uses it to close the same-second-mtime stale hole (the new
bytes are produced regardless, so the rewrite is free).
"""
function ensure_thumb_cached(db::SQLite.DB, exposure_id::Integer, path::String;
                             token::AbstractString = image_version_token(path),
                             overwrite::Bool = false)
    cache_path = thumb_cache_path(db, exposure_id, token)

    if cache_path !== nothing && !overwrite && isfile(cache_path)
        return read(cache_path)
    end

    bytes = encode_png(load_and_lognormalize_thumb(path))
    cache_path === nothing || _write_thumb_atomic(cache_path, bytes)
    bytes
end

"""
    prewarm_thumbnails!(db; threads = true, overwrite = false, experiment_id = nothing,
                        on_progress = nothing) -> NamedTuple

Populate the on-disk thumbnail cache for every exposure with a non-NULL
`image_path`. Run from `init` and `reingest` (issue #261) so the very first
contact-sheet visit on a freshly-prepared corpus is fast (no cold-load stagger).

The exposure list is materialized via `Tables.rowtable` *before* the parallel
loop so the threaded workers touch only the filesystem — never the shared DB
connection (the residual #122 reader-vs-writer `stmt_wrappers` race). Workers are
FS-only by construction.

Unreadable source TIFFs are skipped-with-`@info` (matching the route's graceful
404), so a partially-present corpus still warms what it can — both a file missing
from disk AND one that is present but fails to decode. The decode guard matters
because `work` runs under `@threads`: an uncaught throw there aborts the rest of
the batch, silently leaving every later exposure unwarmed. Returns
`(warmed, skipped)` counts.

`experiment_id` restricts the warm to one experiment. `scan_and_group!` passes
its own id: the whole-DB default walks EVERY exposure in EVERY experiment on
every scan, which is O(corpus) per ingest and a long silent tail after the last
progress tick — the bar sits at 100% while an unrelated experiment's thumbnails
regenerate.

`overwrite = true` re-renders unconditionally, defeating whole-second mtime token
granularity (a TIFF rewritten in the same second as its render keeps the same
token, so the cache would serve a stale PNG). NOTE the ingest path deliberately
passes `overwrite = false`: closing that same-second hole costs a full re-render
of every exposure on every scan, which is what made the unscoped call so
expensive. With `false`, an exposure whose TIFF is rewritten in a *different*
second gets a fresh token → fresh cache path → re-renders normally; only the
same-second rewrite is missed. On a `:memory:` DB this is a no-op pass (rendering
with nowhere to write).
"""
function prewarm_thumbnails!(db::SQLite.DB; threads::Bool = true,
                             overwrite::Bool = false,
                             experiment_id::Union{Integer, Nothing} = nothing,
                             on_progress::Union{Function, Nothing} = nothing)
    rows = experiment_id === nothing ?
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, image_path FROM exposures WHERE image_path IS NOT NULL")) :
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, image_path FROM exposures WHERE image_path IS NOT NULL AND experiment_id = ?",
            [Int(experiment_id)]))

    n_warmed  = Atomic{Int}(0)
    n_skipped = Atomic{Int}(0)
    # Progress counts rows RETIRED (warmed or skipped), so the segment fills even
    # on a corpus that is mostly cache hits or mostly missing files. Atomic
    # because `work` runs under `@threads` — and the tick is emitted from the
    # worker, so `on_progress` must tolerate being called off the main task
    # (broadcast_progress! is lock-guarded, so it does).
    n_done   = Atomic{Int}(0)
    n_total  = length(rows)
    step     = max(1, n_total ÷ 100)
    tick = function ()
        d = atomic_add!(n_done, 1) + 1     # atomic_add! returns the OLD value
        if on_progress !== nothing && (d % step == 0 || d == n_total)
            on_progress(d, n_total)
        end
    end

    work = function (row)
        path = String(row.image_path)
        if !isfile(path)
            atomic_add!(n_skipped, 1)
            @info "thumb prewarm skip" exposure_id=row.id reason="TIFF not on disk" path
            tick()
            return
        end
        # A present-but-undecodable TIFF must not abort the @threads batch and
        # take every later exposure down with it — skip it like a missing file.
        try
            ensure_thumb_cached(db, Int(row.id), path; overwrite = overwrite)
        catch e
            atomic_add!(n_skipped, 1)
            @info "thumb prewarm skip" exposure_id=row.id reason="render failed" path exception=e
            tick()
            return
        end
        atomic_add!(n_warmed, 1)
        tick()
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

    # Close the segment explicitly, mirroring the analyze loop's terminal tick.
    # Two reasons the in-loop ticks cannot be trusted to land it full: a zero-row
    # prewarm never ticks at all (segment stuck at 0), and under `threads` the
    # per-worker ticks can be broadcast out of order, so the LAST frame observed is
    # not necessarily the highest. This fires after @threads has joined.
    on_progress === nothing || on_progress(n_total, n_total)

    @info "thumb prewarm complete" warmed=n_warmed[] skipped=n_skipped[]
    (warmed = n_warmed[], skipped = n_skipped[])
end
