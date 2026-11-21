using Test, HTTP, JSON3

@testset "trace route" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "cubic_tot.dat"),
       joinpath(analysis_dir, "cubic_tot.dat"))
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "form-factor_tot.dat"),
       joinpath(analysis_dir, "form-factor_tot.dat"))
    db     = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    sib_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="cubic_tot",
        transmitted_counts=1.0, exposure_time=5.0, counting_time=5.4)
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot",
        transmitted_counts=7365.0, exposure_time=15.0, counting_time=15.8)
    bare_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="form-factor_tot")

    with_inproc_routes(db) do call
        r = call("GET", "/api/exposures/$e_id/trace")
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test haskey(body, :q) && haskey(body, :I) && haskey(body, :sigma)
        @test body.monitors.transmitted_counts ≈ 7365.0
        @test body.monitors.exposure_time ≈ 15.0
        @test body.monitors.counting_time ≈ 15.8
        @test !haskey(body, :frames)
        @test length(body.q) == length(body.I) == length(body.sigma)
        @test length(body.q) > 100
        @test all(q -> q > 0, body.q)

        sib = JSON3.read(String(call("GET", "/api/exposures/$sib_id/trace").body))
        @test sib.monitors.transmitted_counts ≈ 1.0
        @test sib.monitors.exposure_time ≈ 5.0
        @test sib.monitors.counting_time ≈ 5.4
        @test sib.I != body.I

        bare = JSON3.read(String(call("GET", "/api/exposures/$bare_id/trace").body))
        @test haskey(bare.monitors, :transmitted_counts) && bare.monitors.transmitted_counts === nothing
        @test haskey(bare.monitors, :exposure_time) && bare.monitors.exposure_time === nothing
        @test haskey(bare.monitors, :counting_time) && bare.monitors.counting_time === nothing

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

    config_blob = """
    [experiment]
    name = "PatternTest"
    [layout]
    data_dir = "data"
    analysis_dir = "analysis/automatic_analysis"
    exposure_type = "simple"
    [files]
    integration = "{name}_tot.dat"
    image = "{name}_0_001.tiff"
    """

    db     = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir,
        config=config_blob)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="JC_D01_1_S2449")

    with_inproc_routes(db) do call
        r = call("GET", "/api/exposures/$e_id/trace")
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test length(body.q) > 100
    end
end

@testset "trace route resolves a _0_001 filename to its own exposure's _tot.dat" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "JC_D01_1_S2449_tot.dat"))
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "cubic_tot.dat"),
       joinpath(analysis_dir, "JC_D01_1_S2448_tot.dat"))

    config_blob = """
    [experiment]
    name = "MigratedTrace"
    [layout]
    data_dir = "data"
    analysis_dir = "analysis/automatic_analysis"
    exposure_type = "simple"
    [files]
    integration = "{name}_tot.dat"
    image = "{name}_0_001.tiff"
    """
    db     = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir, config=config_blob)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    sib_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id,
        filename="JC_D01_1_S2448_0_001", transmitted_counts=6000.0, exposure_time=15.0, counting_time=15.8)
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id,
        filename="JC_D01_1_S2449_0_001", transmitted_counts=7000.0, exposure_time=15.0, counting_time=15.8)

    with_inproc_routes(db) do call
        r = call("GET", "/api/exposures/$e_id/trace")
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test length(body.q) > 100
        @test body.monitors.transmitted_counts ≈ 7000.0

        sib = JSON3.read(String(call("GET", "/api/exposures/$sib_id/trace").body))
        @test sib.monitors.transmitted_counts ≈ 6000.0
        @test sib.I != body.I
    end
end
