using Test, FileIO, TiffImages, ImageCore, ColorTypes
using ImageCore.FixedPointNumbers: Q0f31

@testset "load_and_lognormalize" begin
    # Build a 3x3 with three regions:
    #   (1,1) dead pixel  — negative raw count, must clip to pure black
    #   count=50 cells   — faint background, must collapse to black via p5 low clip
    #   count=100 cells  — mid-range, must fall strictly between black and white
    #   count=1e6 cells  — direct beam / bright ring, must clip to white via p99
    # The Q0f31 bit pattern survives the TIFF round trip; the production path
    # recovers raw counts via `reinterpret(Int32, channelview(raw))`.
    counts = Int32[
        -1     50      100;
         100   100     1_000_000;
         50    100     1_000_000;
    ]
    img = Gray.(reinterpret.(Q0f31, counts))
    path = tempname() * ".tiff"
    save(path, img)
    try
        out = HimalayaUI.load_and_lognormalize(path)
        arr = Float32.(channelview(out))

        @test arr[1, 1] == 0f0           # dead pixel → black (same as background)
        @test arr[1, 2] == 0f0           # faint background (count=50) → black
        @test arr[3, 1] == 0f0           # faint background (count=50) → black
        @test arr[2, 3] == 1f0           # bright (count=1e6) → white (clipped at hi)
        @test arr[3, 3] == 1f0           # bright (count=1e6) → white
        @test 0f0 < arr[1, 3] < 1f0      # mid (count=100) sits strictly between
    finally
        rm(path; force=true)
    end
end

# ---------------------------------------------------------------------------
# Issue #261 (H) — thumb downscale-before-lognormalize + disk cache + prewarm
# ---------------------------------------------------------------------------

using SQLite, DBInterface

"""
Build a synthetic detector TIFF: a `dim`×`dim` faint background with a BROAD
bright square block (issue #261 saxs test-pin 2 — a broad feature, NOT a 1px
ring, which can legitimately wash out at 128px). Returns the file path.
"""
function _make_detector_tiff(dim::Int)
    counts = fill(Int32(50), dim, dim)             # faint background
    lo = dim ÷ 4; hi = 3 * dim ÷ 4
    counts[lo:hi, lo:hi] .= Int32(1_000_000)       # broad bright block
    img = Gray.(reinterpret.(Q0f31, counts))
    path = tempname() * ".tiff"
    save(path, img)
    path
end

@testset "load_and_lognormalize_thumb — downscale-first, broad feature survives" begin
    path = _make_detector_tiff(512)
    try
        out = HimalayaUI.load_and_lognormalize_thumb(path, 128)
        @test maximum(size(out)) <= 128          # actually downscaled
        arr = Float32.(channelview(out))
        @test all(0f0 .<= arr .<= 1f0)           # normalized range
        # A BROAD bright block must still read near-white at 128px; the faint
        # background floor must collapse toward black. Floor-style assertions.
        @test maximum(arr) > 0.9f0
        @test minimum(arr) < 0.1f0
        # The bright block (center) is brighter than a background corner.
        c = size(arr) .÷ 2
        @test arr[c[1], c[2]] > arr[1, 1]
    finally
        rm(path; force=true)
    end
end

@testset "thumb_cache_dir / thumb_cache_path — :memory: guard" begin
    mem = SQLite.DB()
    @test HimalayaUI.thumb_cache_dir(mem) === nothing
    @test HimalayaUI.thumb_cache_path(mem, 1, "v3-123") === nothing

    mktempdir() do tmp
        disk = SQLite.DB(joinpath(tmp, "himalaya.db"))
        dir = HimalayaUI.thumb_cache_dir(disk)
        @test dir == joinpath(tmp, "cache", "thumb-128")
        @test HimalayaUI.thumb_cache_path(disk, 7, "v3-123") ==
              joinpath(tmp, "cache", "thumb-128", "7-v3-123.png")
    end
end

@testset "ensure_thumb_cached — miss writes, hit reads from disk" begin
    mktempdir() do tmp
        disk = SQLite.DB(joinpath(tmp, "himalaya.db"))
        path = _make_detector_tiff(256)
        try
            token = HimalayaUI.image_version_token(path)
            cpath = HimalayaUI.thumb_cache_path(disk, 42, token)
            @test !isfile(cpath)

            # Miss: renders + writes the cache file.
            b1 = HimalayaUI.ensure_thumb_cached(disk, 42, path)
            @test isfile(cpath)
            @test length(b1) > 100

            # Hit: serve from disk even after the SOURCE TIFF is gone — proves the
            # second call read the cached file rather than re-rendering.
            rm(path; force=true)
            b2 = HimalayaUI.ensure_thumb_cached(disk, 42, path)
            @test b2 == b1
        finally
            rm(path; force=true)
        end
    end
end

@testset "ensure_thumb_cached — token change writes a NEW file, leaves the stale one" begin
    mktempdir() do tmp
        disk = SQLite.DB(joinpath(tmp, "himalaya.db"))
        path = _make_detector_tiff(256)
        try
            HimalayaUI.ensure_thumb_cached(disk, 9, path)
            token1 = HimalayaUI.image_version_token(path)
            stale  = HimalayaUI.thumb_cache_path(disk, 9, token1)
            @test isfile(stale)

            # Bump the source mtime into the future → a new token → a new key.
            touch(path; mtime = time() + 5)
            token2 = HimalayaUI.image_version_token(path)
            @test token2 != token1
            HimalayaUI.ensure_thumb_cached(disk, 9, path)
            fresh = HimalayaUI.thumb_cache_path(disk, 9, token2)
            @test isfile(fresh)
            @test isfile(stale)            # stale entry left in place (bounded)
        finally
            rm(path; force=true)
        end
    end
end

@testset "ensure_thumb_cached — :memory: renders without writing" begin
    mem = SQLite.DB()
    path = _make_detector_tiff(256)
    try
        b = HimalayaUI.ensure_thumb_cached(mem, 1, path)
        @test length(b) > 100             # rendered fresh, no cache dir, no throw
    finally
        rm(path; force=true)
    end
end

@testset "prewarm_thumbnails! — populates cache, skips missing TIFF, ingest-shaped" begin
    mktempdir() do tmp
        # db_dir is the tmp dir; deliberately NOT an experiment dir, so the
        # cache writes never collide with the read-only-experiment-dir invariant.
        db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
        present = _make_detector_tiff(256)
        gone    = tempname() * ".tiff"     # path set, file absent
        try
            exp = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
            samp = HimalayaUI.create_sample!(db; experiment_id=exp)
            e1 = HimalayaUI.create_exposure!(db; sample_id=samp, image_path=present)
            HimalayaUI.create_exposure!(db; sample_id=samp, image_path=gone)
            HimalayaUI.create_exposure!(db; sample_id=samp)   # NULL image_path

            res = HimalayaUI.prewarm_thumbnails!(db; threads=false)
            @test res.warmed == 1
            @test res.skipped == 1          # the missing TIFF, skipped not fatal

            token = HimalayaUI.image_version_token(present)
            @test isfile(HimalayaUI.thumb_cache_path(db, e1, token))

            # threads=true smoke run must not throw.
            res2 = HimalayaUI.prewarm_thumbnails!(db; threads=true)
            @test res2.warmed == 1

            # overwrite=true re-renders the present one (the reingest contract).
            res3 = HimalayaUI.prewarm_thumbnails!(db; threads=false, overwrite=true)
            @test res3.warmed == 1
        finally
            rm(present; force=true)
        end
    end
end
