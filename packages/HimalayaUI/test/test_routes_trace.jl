using Test, HTTP, JSON3

@testset "trace route" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot",
        transmitted_counts=7365.0, exposure_time=15.0, counting_time=15.8)
    HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot_0_002",
        transmitted_counts=1.0)

    with_inproc_routes(db) do call
        r = call("GET", "/api/exposures/$e_id/trace")
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test haskey(body, :q) && haskey(body, :I) && haskey(body, :sigma)
        @test length(body.frames) == 1
        @test body.frames[1].exposure_id == e_id
        @test body.frames[1].transmitted_counts ≈ 7365.0
        @test body.frames[1].exposure_time ≈ 15.0
        @test body.frames[1].counting_time ≈ 15.8
        @test length(body.q) == length(body.I) == length(body.sigma)
        @test length(body.q) > 100
        @test all(q -> q > 0, body.q)

        # 404 for unknown exposure
        r = call("GET", "/api/exposures/99999/trace")
        @test r.status == 404
    end
end

@testset "trace route honors integration_pattern" begin
    # Regression: filename stem stored in DB lacks the integration suffix.
    # With integration = "{name}_tot.dat", the actual file is
    # `<stem>_tot.dat`. The route must consult the experiment's config
    # rather than hardcoding ".dat".
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "JC_D01_1_S2449_tot.dat"))


    db     = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir,
        integration_pattern="{name}_tot.dat")
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="JC_D01_1_S2449")

    with_inproc_routes(db) do call
        r = call("GET", "/api/exposures/$e_id/trace")
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test length(body.q) > 100
    end
end

@testset "trace route resolves migrated full-stem filename to per-acquisition _tot.dat" begin
    # Migration scenario: filename is the full PER-FRAME stem (rewritten for rescan
    # dedup), but the trace total is named by the ACQUISITION stem (frame suffix
    # dropped). The route must fall back via the shared resolver.
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "JC_D01_1_S2449_tot.dat"))   # acquisition-stem name

    db     = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir, integration_pattern="{name}_tot.dat")
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    # Full per-frame stem.
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id,
        filename="JC_D01_1_S2449_0_001", transmitted_counts=7000.0, exposure_time=15.0, counting_time=15.8)
    e2_id  = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id,
        filename="JC_D01_1_S2449_0_002", transmitted_counts=6000.0, exposure_time=15.0, counting_time=15.8)
    HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id,
        filename="JC_D01_1_S24490_0_001", transmitted_counts=1.0)

    with_inproc_routes(db) do call
        r = call("GET", "/api/exposures/$e_id/trace")
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test length(body.q) > 100
        @test [f.exposure_id for f in body.frames] == [e_id, e2_id]
        @test [f.transmitted_counts for f in body.frames] ≈ [7000.0, 6000.0]
        r2 = call("GET", "/api/exposures/$e2_id/trace")
        @test [f.exposure_id for f in JSON3.read(String(r2.body)).frames] == [e_id, e2_id]
    end
end
