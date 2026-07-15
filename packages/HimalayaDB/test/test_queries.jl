@testset "navigation" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    exps = experiments(db)
    @test length(exps) == 1
    @test exps[1].id == ids.experiment_id
    @test exps[1].name == "fixture"

    smps = samples(db; experiment=ids.experiment_id)
    @test length(smps) == 1
    @test smps[1].id == ids.sample_id

    exps_all = samples(db)  # no filter
    @test length(exps_all) == 1

    exs = exposures(db; sample=ids.sample_id)
    @test length(exs) == 2                 # exposure_id + the uncurated exposure2_id
    @test exs[1].id == ids.exposure_id
    @test exs[1].filename == "s1"

    close(db)
end
