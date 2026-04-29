using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "mention lookup routes" begin
    tmp = mktempdir()
    db  = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.create_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="A", name="sampleA")
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="run001")

    # Insert a peak manually
    res = DBInterface.execute(db,
        "INSERT INTO peaks (exposure_id, q, source, excluded) VALUES (?, ?, 'auto', 0)",
        [e_id, 1.223])
    pk_id = Int(DBInterface.lastrowid(res))

    with_test_server(db) do port, base
        @testset "GET /api/peaks/:id" begin
            r = HTTP.get("$base/api/peaks/$pk_id")
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test body.id == pk_id
            @test body.q ≈ 1.223
            @test body.source == "auto"
            @test body.excluded === false

            r404 = HTTP.get("$base/api/peaks/999999"; status_exception=false)
            @test r404.status == 404
        end

        @testset "GET /api/exposures/:id" begin
            r = HTTP.get("$base/api/exposures/$e_id")
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test body.id == e_id
            @test body.filename == "run001"

            r404 = HTTP.get("$base/api/exposures/999999"; status_exception=false)
            @test r404.status == 404
        end

        @testset "GET /api/samples/:id" begin
            r = HTTP.get("$base/api/samples/$s_id")
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test body.id == s_id
            @test body.name == "sampleA"

            r404 = HTTP.get("$base/api/samples/999999"; status_exception=false)
            @test r404.status == 404
        end

        @testset "GET /api/indices/:id" begin
            # No indices until analysis runs — just test 404 path
            r404 = HTTP.get("$base/api/indices/999999"; status_exception=false)
            @test r404.status == 404
        end
    end
end
