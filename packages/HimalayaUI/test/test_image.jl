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

@testset "a flat dark-floor frame (beam off) renders and encodes" begin
    counts = fill(Int32(3), 64, 64)
    counts[1, 1] = Int32(-1)
    counts[2, 2] = Int32(40)
    img = Gray.(reinterpret.(Q0f31, counts))
    path = tempname() * ".tiff"
    save(path, img)
    try
        for out in (HimalayaUI.load_and_lognormalize(path), HimalayaUI.load_and_lognormalize_thumb(path, 16))
            arr = Float32.(channelview(out))
            @test all(0f0 .<= arr .<= 1f0)
            @test !isempty(HimalayaUI.encode_png(out))
        end
        arr = Float32.(channelview(HimalayaUI.load_and_lognormalize(path)))
        @test arr[5, 5] == 0f0
        @test arr[2, 2] == 1f0
    finally
        rm(path; force=true)
    end
end

# ---------------------------------------------------------------------------
# Issue #261 (H) — thumb downscale-before-lognormalize + disk cache + prewarm
# ---------------------------------------------------------------------------

using SQLite, DBInterface, Tables

import JSON

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
            b1 = HimalayaUI.ensure_thumb_cached(disk, 42, path; token=token)
            @test isfile(cpath)
            @test length(b1.bytes) > 100
            @test b1.raw_h == 256 && b1.raw_w == 256

            # Hit: serve from disk even after the SOURCE TIFF is gone — proves the
            # second call read the cached file rather than re-rendering. The route
            # computes the token while the source is present, so we pass the same
            # `token` here (image_version_token() of a deleted file is "", which
            # would compute the wrong key — exactly the resilience contract).
            rm(path; force=true)
            b2 = HimalayaUI.ensure_thumb_cached(disk, 42, path; token=token)
            @test b2 == b1
        finally
            rm(path; force=true)
        end
    end
end

@testset "ensure_thumb_cached: a thumbnail cached before the size was stored gets its size once, then never reparses" begin
    mktempdir() do tmp
        disk = SQLite.DB(joinpath(tmp, "himalaya.db"))
        rows = Int32[50 for _ in 1:96, _ in 1:160]
        path = tempname() * ".tiff"
        save(path, Gray.(reinterpret.(Q0f31, rows)))
        try
            token = HimalayaUI.image_version_token(path)
            cpath = HimalayaUI.thumb_cache_path(disk, 5, token)
            old = HimalayaUI.ensure_thumb_cached(disk, 5, path; token=token)
            rm(HimalayaUI._full_size_path(cpath))

            healed = HimalayaUI.ensure_thumb_cached(disk, 5, path; token=token)
            @test healed.bytes == old.bytes
            @test healed.raw_h == 96 && healed.raw_w == 160
            @test isfile(HimalayaUI._full_size_path(cpath))

            rm(path; force=true)
            again = HimalayaUI.ensure_thumb_cached(disk, 5, path; token=token)
            @test again.raw_h == 96 && again.raw_w == 160
        finally
            rm(path; force=true)
        end
    end
end

@testset "ensure_full_cached: a miss renders and stores, a hit keeps the raw size without the source" begin
    mktempdir() do tmp
        disk = SQLite.DB(joinpath(tmp, "himalaya.db"))
        path = _make_detector_tiff(256)
        try
            token = HimalayaUI.image_version_token(path)
            first = HimalayaUI.ensure_full_cached(disk, 42, path; token=token)
            @test first.raw_h == 256 && first.raw_w == 256
            @test first.bytes == HimalayaUI.encode_png(HimalayaUI.load_and_lognormalize(path))
            @test isfile(HimalayaUI.full_cache_path(disk, 42, token))
            @test HimalayaUI.full_cache_path(disk, 42, token) != HimalayaUI.thumb_cache_path(disk, 42, token)

            rm(path; force=true)
            again = HimalayaUI.ensure_full_cached(disk, 42, path; token=token)
            @test again == first
        finally
            rm(path; force=true)
        end
    end
    mem = SQLite.DB()
    path = _make_detector_tiff(64)
    try
        @test HimalayaUI.full_cache_path(mem, 1, "v3-1") === nothing
        @test HimalayaUI.ensure_full_cached(mem, 1, path).raw_w == 64
    finally
        rm(path; force=true)
    end
end

@testset "image cache: the cache stays under its byte cap by evicting the least recently viewed frames" begin
    mktempdir() do tmp
        disk = SQLite.DB(joinpath(tmp, "himalaya.db"))
        path = _make_detector_tiff(64)
        try
            token = HimalayaUI.image_version_token(path)
            entry(id) = HimalayaUI.full_cache_path(disk, id, token)
            HimalayaUI.ensure_full_cached(disk, 1, path; token=token)
            one = filesize(entry(1)) + filesize(HimalayaUI._full_size_path(entry(1)))
            cap = 3one - 1
            viewed(id) = (sleep(0.05); HimalayaUI.ensure_full_cached(disk, id, path; token=token, cap=cap))

            viewed(2)
            viewed(1)
            viewed(3)
            @test isfile(entry(1)) && isfile(entry(3))
            @test !isfile(entry(2)) && !isfile(HimalayaUI._full_size_path(entry(2)))
            @test !isdir(dirname(entry(2)))

            sleep(0.05)
            HimalayaUI.ensure_full_cached(disk, 4, path; token=token, cap=one ÷ 2)
            @test isfile(entry(4)) && !isfile(entry(1)) && !isfile(entry(3))
        finally
            rm(path; force=true)
        end
    end
end

@testset "image cache: old thumbnails make room for a new full frame, and a restart reseeds the total" begin
    mktempdir() do tmp
        disk = SQLite.DB(joinpath(tmp, "himalaya.db"))
        path = _make_detector_tiff(64)
        try
            token = HimalayaUI.image_version_token(path)
            thumbs = [HimalayaUI.ensure_thumb_cached(disk, id, path; token=token) for id in 1:20]
            entry(id) = HimalayaUI.full_cache_path(disk, id, token)
            sleep(0.05)
            HimalayaUI.ensure_full_cached(disk, 100, path; token=token)
            one = filesize(entry(100)) + filesize(HimalayaUI._full_size_path(entry(100)))
            cap = 3one
            @test sum(t -> length(t.bytes), thumbs) > cap

            empty!(HimalayaUI._CACHE_BYTES)
            sleep(0.05)
            HimalayaUI.ensure_full_cached(disk, 101, path; token=token, cap=cap)
            @test isfile(entry(100)) && isfile(entry(101))
            @test count(id -> isfile(HimalayaUI.thumb_cache_path(disk, id, token)), 1:20) < 20
            root = joinpath(tmp, "cache")
            @test HimalayaUI._CACHE_BYTES[root] == first(HimalayaUI._cache_usage(root))
            @test HimalayaUI._CACHE_BYTES[root] <= cap
        finally
            rm(path; force=true)
        end
    end
end

@testset "ensure_full_cached: a new token replaces the old entry, and a broken entry renders again" begin
    mktempdir() do tmp
        disk = SQLite.DB(joinpath(tmp, "himalaya.db"))
        path = _make_detector_tiff(64)
        try
            old = HimalayaUI.full_cache_path(disk, 5, "v0-1")
            HimalayaUI.ensure_full_cached(disk, 5, path; token="v0-1")
            HimalayaUI.ensure_full_cached(disk, 6, path; token="v0-1")
            @test isfile(old)
            fresh = HimalayaUI.full_cache_path(disk, 5, "v0-2")
            HimalayaUI.ensure_full_cached(disk, 5, path; token="v0-2")
            @test sort(readdir(dirname(fresh))) == sort(basename.([fresh, HimalayaUI._full_size_path(fresh)]))
            @test isfile(HimalayaUI.full_cache_path(disk, 6, "v0-1"))

            rm(fresh)
            @test HimalayaUI.ensure_full_cached(disk, 5, path; token="v0-2").raw_w == 64
            write(HimalayaUI._full_size_path(fresh), "garbage")
            @test HimalayaUI.ensure_full_cached(disk, 5, path; token="v0-2").raw_h == 64
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
            token1 = HimalayaUI.image_version_token(path)
            HimalayaUI.ensure_thumb_cached(disk, 9, path; token=token1)
            stale  = HimalayaUI.thumb_cache_path(disk, 9, token1)
            @test isfile(stale)

            # Age the source mtime into the future so image_version_token yields a
            # different token (the route's natural invalidator on a TIFF rewrite).
            # Julia's `touch` has no mtime kwarg; shell out to `touch -t`.
            future = Libc.strftime("%Y%m%d%H%M", time() + 120)
            run(`touch -t $future $path`)
            token2 = HimalayaUI.image_version_token(path)
            @test token2 != token1
            HimalayaUI.ensure_thumb_cached(disk, 9, path; token=token2)
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
        @test length(b.bytes) > 100
        @test b.raw_h == 256 && b.raw_w == 256
    finally
        rm(path; force=true)
    end
end

@testset "prewarm_thumbnails! — populates cache, skips missing TIFF, ingest-shaped" begin
    mktempdir() do tmp
        # db_dir is the tmp dir; deliberately NOT an experiment dir, so the
        # cache writes never collide with the read-only-experiment-dir invariant.
        db = open_prepared_clone(tmp)
        present = _make_detector_tiff(256)
        gone    = tempname() * ".tiff"     # path set, file absent
        try
            exp = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
            samp = HimalayaUI.create_sample!(db; experiment_id=exp, name="S")
            e1 = HimalayaUI.create_exposure!(db; experiment_id=exp, sample_id=samp, image_path=present)
            HimalayaUI.create_exposure!(db; experiment_id=exp, sample_id=samp, image_path=gone)
            HimalayaUI.create_exposure!(db; experiment_id=exp, sample_id=samp)   # NULL image_path

            res = HimalayaUI.prewarm_thumbnails!(db; threads=false)
            @test res.warmed == 1
            @test res.skipped == 1          # the missing TIFF, skipped not fatal

            token = HimalayaUI.image_version_token(present)
            @test isfile(HimalayaUI.thumb_cache_path(db, e1, token))
            @test isfile(HimalayaUI._full_size_path(HimalayaUI.thumb_cache_path(db, e1, token)))

            # threads=true smoke run must not throw.
            res2 = HimalayaUI.prewarm_thumbnails!(db; threads=true)
            @test res2.warmed == 1
        finally
            rm(present; force=true)
        end
    end

end

@testset "prewarm_thumbnails! experiment_id scopes the warm to one experiment" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        mine  = _make_detector_tiff(256)
        other = _make_detector_tiff(256)
        try
            exp_a = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
            exp_b = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
            samp_a = HimalayaUI.create_sample!(db; experiment_id=exp_a, name="A")
            samp_b = HimalayaUI.create_sample!(db; experiment_id=exp_b, name="B")
            e_a = HimalayaUI.create_exposure!(db; experiment_id=exp_a, sample_id=samp_a, image_path=mine)
            e_b = HimalayaUI.create_exposure!(db; experiment_id=exp_b, sample_id=samp_b, image_path=other)

            # Scoped warm touches ONLY exp_a. The unscoped default warmed every
            # exposure in every experiment on every scan — O(corpus) per ingest.
            res = HimalayaUI.prewarm_thumbnails!(db; threads=false, experiment_id=exp_a)
            @test res.warmed == 1
            @test isfile(HimalayaUI.thumb_cache_path(db, e_a, HimalayaUI.image_version_token(mine)))
            @test !isfile(HimalayaUI.thumb_cache_path(db, e_b, HimalayaUI.image_version_token(other)))

            # An experiment with no exposures warms nothing and does not throw.
            exp_c = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
            @test HimalayaUI.prewarm_thumbnails!(db; threads=false, experiment_id=exp_c).warmed == 0

            # The scoped pass is the REPAIR path: re-running it warms whatever is
            # still missing. A new-exposures-only warm would never revisit exp_b.
            res_b = HimalayaUI.prewarm_thumbnails!(db; threads=false, experiment_id=exp_b)
            @test res_b.warmed == 1
            @test isfile(HimalayaUI.thumb_cache_path(db, e_b, HimalayaUI.image_version_token(other)))
        finally
            rm(mine;  force=true)
            rm(other; force=true)
        end
    end
end

@testset "prewarm_thumbnails! skips an undecodable TIFF without aborting the batch" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        good = _make_detector_tiff(256)
        # Present on disk, but not a decodable TIFF — the !isfile guard does not
        # catch this, and an uncaught throw inside the @threads loop used to take
        # every later exposure in the batch down with it.
        bad = tempname() * ".tiff"
        write(bad, "this is not a tiff")
        try
            exp = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
            samp = HimalayaUI.create_sample!(db; experiment_id=exp, name="S")
            # `bad` is inserted BEFORE `good` on purpose: pre-guard the throw took
            # every later exposure in the batch down with it.
            HimalayaUI.create_exposure!(db; experiment_id=exp, sample_id=samp, image_path=bad)
            e_good = HimalayaUI.create_exposure!(db; experiment_id=exp, sample_id=samp, image_path=good)

            res = HimalayaUI.prewarm_thumbnails!(db; threads=false, experiment_id=exp)
            @test res.skipped == 1
            @test res.warmed == 1
            # The good one after the bad one still got warmed.
            @test isfile(HimalayaUI.thumb_cache_path(db, e_good, HimalayaUI.image_version_token(good)))
        finally
            rm(good; force=true)
            rm(bad;  force=true)
        end
    end
end

@testset "an ingest pass completes before its thumbnails are warmed" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        data = joinpath(tmp, "data"); mkpath(data)
        tiff = _make_detector_tiff(256)
        cp(tiff, joinpath(data, "A_S0001_0_001.tif"))
        rm(tiff; force = true)
        exp = HimalayaUI.create_experiment!(db; name = "thumbs", path = tmp, data_dir = data, analysis_dir = data)
        pending = Channel{String}(256)
        sub = (pending = pending,)
        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end
        gate = Channel{Nothing}(1)
        HimalayaUI._THUMB_PREWARM_GATE[] = gate
        try
            @test HimalayaUI.run_ingest!(db, exp) === :ok
            stages = String[]
            kinds = String[]
            while isready(pending)
                frame = take!(pending)
                lines = filter(l -> startswith(l, "data: "), split(frame, '\n'))
                isempty(lines) && continue
                obj = JSON.parse(replace(first(lines), r"^data: " => ""))
                haskey(obj, :kind) && push!(kinds, String(obj.kind))
                haskey(obj, :payload) && haskey(obj.payload, :stage) && push!(stages, String(obj.payload.stage))
            end
            @test "ingest_complete" in kinds
            @test "analyzing" in stages
            @test !("thumbnails" in stages)
            @test first(Tables.rowtable(DBInterface.execute(db,
                "SELECT ingest_status FROM experiments WHERE id = ?", [exp]))).ingest_status == "complete"

            e = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT id, image_path FROM exposures WHERE experiment_id = ?", [exp])))
            thumb = HimalayaUI.thumb_cache_path(db, Int(e.id), HimalayaUI.image_version_token(String(e.image_path)))
            @test !isfile(thumb)

            put!(gate, nothing)
            HimalayaUI.wait_prewarm(exp)
            @test isfile(thumb)
        finally
            HimalayaUI._THUMB_PREWARM_GATE[] = nothing
            isready(gate) || put!(gate, nothing)
            HimalayaUI.wait_prewarm(exp)
            lock(HimalayaUI.SSE_LOCK) do
                filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
            end
            close(pending)
        end
    end
end
