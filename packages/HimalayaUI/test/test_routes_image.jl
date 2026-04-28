using Test, HTTP, JSON3, SQLite, DBInterface, FileIO, ImageCore, TiffImages

@testset "GET /api/exposures/:id/image" begin
    tiff_path = tempname() * ".tiff"
    test_img = Gray.(rand(Float32, 512, 384))  # larger than 128px so thumb is actually smaller
    save(tiff_path, test_img)

    db      = SQLite.DB()
    HimalayaUI.create_schema!(db)
    HimalayaUI.migrate_schema!(db)
    exp_id  = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
    samp_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
    eid     = HimalayaUI.create_exposure!(db; sample_id=samp_id, image_path=tiff_path)
    eid_noi = HimalayaUI.create_exposure!(db; sample_id=samp_id)  # no image

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
    end

    rm(tiff_path; force=true)
end
