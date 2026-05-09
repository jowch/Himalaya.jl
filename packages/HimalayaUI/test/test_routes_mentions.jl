using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "mention lookup routes" begin
    tmp = mktempdir()
    db  = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.create_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="A", display_name="sampleA")
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="run001")

    # Insert a peak manually (use auto_peaks — the new schema table)
    res = DBInterface.execute(db,
        "INSERT INTO auto_peaks (exposure_id, q) VALUES (?, ?)",
        [e_id, 1.223])
    pk_id = Int(DBInterface.lastrowid(res))

    # Insert an index + index_peaks row manually so we can assert the 200 path
    # without standing up a full analysis pipeline (the worktree may not have
    # the matching Himalaya core for findpeaks). Phase "Pn3m" with a basis
    # gives the route enough to compute predicted_q via predicted_q_for_phase.
    res = DBInterface.execute(db,
        """
        INSERT INTO indices (exposure_id, phase, basis, score, r_squared, lattice_d, status)
        VALUES (?, 'Pn3m', 1.223, 0.91, 0.998, 5.14, 'candidate')
        """,
        [e_id])
    ix_id = Int(DBInterface.lastrowid(res))
    DBInterface.execute(db,
        "INSERT INTO index_peaks (index_id, peak_id, peak_kind, ratio_position, residual) VALUES (?, ?, 'auto', 1, 0.0)",
        [ix_id, pk_id])

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
            @test body.name == "A"

            r404 = HTTP.get("$base/api/samples/999999"; status_exception=false)
            @test r404.status == 404
        end

        @testset "GET /api/indices/:id" begin
            r = HTTP.get("$base/api/indices/$ix_id")
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test body.id == ix_id
            @test body.phase == "Pn3m"
            @test body.basis ≈ 1.223
            @test body.score ≈ 0.91
            @test body.r_squared ≈ 0.998
            @test body.lattice_d ≈ 5.14
            @test body.status == "candidate"
            @test length(body.peaks) == 1
            @test body.peaks[1].peak_id == pk_id
            @test body.peaks[1].ratio_position == 1
            @test body.peaks[1].q_observed ≈ 1.223
            # predicted_q is derived from phase-ratios × basis
            @test length(body.predicted_q) > 0
            @test body.predicted_q[1] ≈ 1.223

            r404 = HTTP.get("$base/api/indices/999999"; status_exception=false)
            @test r404.status == 404
        end
    end
end
