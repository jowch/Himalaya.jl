using Test, HTTP, JSON3, SQLite, DBInterface, FileIO, ImageCore, TiffImages

@testset "GET /api/exposures/:id/image" begin
    tiff_path = tempname() * ".tiff"
    test_img = Gray.(rand(Float32, 512, 384))  # larger than 128px so thumb is actually smaller
    save(tiff_path, test_img)

    # image_path points at a file that does not exist on disk (e.g. the source
    # TIFF was moved/deleted after ingest). The route must 404 gracefully, like
    # the trace route does for a missing .dat — not throw an unhandled
    # ArgumentError into a 500 (finding BE-1 / L-2, issue #233).
    missing_path = tempname() * ".tiff"

    db      = SQLite.DB()
    HimalayaUI.create_schema!(db)
    HimalayaUI.migrate_schema!(db)
    exp_id  = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
    samp_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
    eid     = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=samp_id, image_path=tiff_path)
    eid_noi = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=samp_id)  # no image
    eid_gone = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=samp_id, image_path=missing_path)  # path set, file missing

    with_test_server(db) do port, base
        # Full image
        r = HTTP.get("$base/api/exposures/$eid/image")
        @test r.status == 200
        @test Dict(r.headers)["Content-Type"] == "image/png"
        @test length(r.body) > 100

        # Thumb variant is smaller
        rt = HTTP.get("$base/api/exposures/$eid/image?thumb=1")
        @test rt.status == 200
        @test length(rt.body) < length(r.body)

        # null image_path → 404
        r404 = HTTP.get("$base/api/exposures/$eid_noi/image"; status_exception=false)
        @test r404.status == 404

        # nonexistent exposure → 404
        r404b = HTTP.get("$base/api/exposures/9999/image"; status_exception=false)
        @test r404b.status == 404

        # image_path set but source file missing → graceful 404, not 500
        r404c = HTTP.get("$base/api/exposures/$eid_gone/image"; status_exception=false)
        @test r404c.status == 404
        # thumb variant of a missing source must also 404, not 500
        r404d = HTTP.get("$base/api/exposures/$eid_gone/image?thumb=1"; status_exception=false)
        @test r404d.status == 404
    end

    rm(tiff_path; force=true)
end

@testset "GET /api/exposures/:id/image — full path capped at 1536px (issue #260, G)" begin
    # A >1536px detector: the full path must cap the RENDERED raster at 1536
    # max-side (G) while the percentile-clip math still runs at full resolution.
    big_path = tempname() * ".tiff"
    save(big_path, Gray.(rand(Float32, 2048, 2048)))
    # A small detector: resize_to_fit no-ops, dimensions preserved.
    small_path = tempname() * ".tiff"
    save(small_path, Gray.(rand(Float32, 600, 400)))

    db = SQLite.DB()
    HimalayaUI.create_schema!(db); HimalayaUI.migrate_schema!(db)
    exp  = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
    samp = HimalayaUI.create_sample!(db; experiment_id=exp, name="S")
    big  = HimalayaUI.create_exposure!(db; experiment_id=exp, sample_id=samp, image_path=big_path)
    small = HimalayaUI.create_exposure!(db; experiment_id=exp, sample_id=samp, image_path=small_path)

    decode_dims(body) = size(FileIO.load(FileIO.Stream{FileIO.format"PNG"}(IOBuffer(body))))

    with_test_server(db) do port, base
        rb = HTTP.get("$base/api/exposures/$big/image")
        @test rb.status == 200
        @test maximum(decode_dims(rb.body)) <= 1536      # capped

        rs = HTTP.get("$base/api/exposures/$small/image")
        @test rs.status == 200
        @test maximum(decode_dims(rs.body)) == 600       # untouched (<=1536)
    end

    rm(big_path; force=true); rm(small_path; force=true)
end

@testset "full image emits raw X-Image-Width/Height; thumb does not" begin
    big_path = tempname() * ".tiff"
    save(big_path, Gray.(rand(Float32, 2048, 2048)))      # > 1536 → resized for display
    rect_path = tempname() * ".tiff"
    save(rect_path, Gray.(rand(Float32, 400, 600)))        # (rows=400, cols=600) → W=600, H=400

    db = SQLite.DB()
    HimalayaUI.create_schema!(db); HimalayaUI.migrate_schema!(db)
    exp  = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
    samp = HimalayaUI.create_sample!(db; experiment_id=exp, name="S")
    big  = HimalayaUI.create_exposure!(db; experiment_id=exp, sample_id=samp, image_path=big_path)
    rect = HimalayaUI.create_exposure!(db; experiment_id=exp, sample_id=samp, image_path=rect_path)

    with_test_server(db) do port, base
        rb = HTTP.get("$base/api/exposures/$big/image")
        h  = Dict(rb.headers)
        @test h["X-Image-Width"]  == "2048"   # RAW, not the <=1536 displayed size
        @test h["X-Image-Height"] == "2048"

        rr = HTTP.get("$base/api/exposures/$rect/image")
        hr = Dict(rr.headers)
        @test hr["X-Image-Width"]  == "600"   # cols
        @test hr["X-Image-Height"] == "400"   # rows

        rt = HTTP.get("$base/api/exposures/$big/image?thumb=1")
        @test !haskey(Dict(rt.headers), "X-Image-Width")
    end
    rm(big_path; force=true); rm(rect_path; force=true)
end
