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
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot")

    with_inproc_routes(db) do call
        r = call("GET", "/api/exposures/$e_id/trace")
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test haskey(body, :q) && haskey(body, :I) && haskey(body, :sigma)
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
