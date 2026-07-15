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

@testset "curated_peaks" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    peaks = curated_peaks(db, ids.exposure_id)
    # 3 auto + 1 manual add = 4 rows
    @test length(peaks) == 4
    autos = filter(p -> p.source == "auto", peaks)
    manuals = filter(p -> p.source == "manual", peaks)
    @test length(autos) == 3
    @test length(manuals) == 1
    @test manuals[1].q == 0.20
    # the middle auto peak (0.1414) is excluded
    excluded = filter(p -> p.excluded == 1, peaks)
    @test length(excluded) == 1
    @test excluded[1].q == 0.1414
    # sorted by q
    @test issorted([p.q for p in peaks])

    close(db)
end

@testset "indices" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    cands = index_candidates(db, ids.exposure_id)
    @test length(cands) == 1
    @test cands[1].id == ids.index_id
    @test cands[1].phase == "Pn3m"
    @test cands[1].kind == "auto"

    conf = confirmed_indices(db, ids.exposure_id)
    @test length(conf) == 1
    @test conf[1].id == ids.index_id
    @test conf[1].phase == "Pn3m"

    # uncurated exposure2: has a candidate index in an ACTIVE auto group, but no
    # custom group -> confirmed_indices must be empty (filter is kind='custom',
    # not active=1). index_candidates still sees the index.
    @test length(index_candidates(db, ids.exposure2_id)) == 1
    @test isempty(confirmed_indices(db, ids.exposure2_id))

    close(db)
end
