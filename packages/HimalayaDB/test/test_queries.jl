@testset "navigation" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    exps = experiments(db)
    @test length(exps) == 1
    @test exps[1].id == ids.experiment_id
    @test exps[1].name == "fixture"
    @test exps[1].path == dir

    sts = setups(db; experiment=ids.experiment_id)
    @test length(sts) == 1
    @test sts[1].id == ids.setup_id
    @test sts[1].analysis_dir == dir
    @test sts[1].energy_kev == 11.0
    @test sts[1].flight_path_m == 1.1
    @test sts[1].integration_pattern == "{name}_tot.dat"
    @test sts[1].trace_normalisation == "i2"
    @test sts[1].normalisation_base == 2.5e7
    @test sts[1].mode == "capillary"
    @test sts[1].mode_source == "user"
    @test length(setups(db)) == 1
    @test isempty(setups(db; experiment=999999))

    smps = samples(db; experiment=ids.experiment_id)
    @test [s.id for s in smps] == [ids.sample_id, ids.sample2_id, ids.empty_id]
    @test smps[1].setup_id == ids.setup_id
    @test length(samples(db)) == 3
    @test isempty(exposures(db; sample=ids.empty_id))

    exs = exposures(db; sample=ids.sample_id)
    @test length(exs) == 2
    @test exs[1].id == ids.exposure_id
    @test exs[1].filename == "s1"
    @test exs[1].setup_id == ids.setup_id
    @test exs[1].experiment_id == ids.experiment_id

    close(db)
end

@testset "curated_peaks" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    peaks = curated_peaks(db, ids.exposure_id)
    @test length(peaks) == 4
    autos = filter(p -> p.source == "auto", peaks)
    manuals = filter(p -> p.source == "manual", peaks)
    @test length(autos) == 3
    @test length(manuals) == 1
    @test manuals[1].q == 0.20
    excluded = filter(p -> p.excluded == 1, peaks)
    @test length(excluded) == 1
    @test excluded[1].q == 0.1414
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

    @test length(index_candidates(db, ids.exposure2_id)) == 1
    @test isempty(confirmed_indices(db, ids.exposure2_id))

    close(db)
end

@testset "phase_call" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    call = phase_call(db, ids.sample_id)
    @test call.state == "indexed"
    @test call.phase == "Pn3m"
    @test call.call_phases == ["Im3m", "Pn3m"]
    @test [m.index_id for m in call.members] == [ids.index_id, ids.index2_id]
    @test [m.exposure_id for m in call.members] == [ids.exposure_id, ids.exposure2_id]

    shared = phase_call(db, ids.sample2_id)
    @test shared.state == "indexed"
    @test shared.phase == "Pn3m"
    @test shared.call_phases == ["Im3m", "Pn3m"]
    @test [m.exposure_id for m in shared.members] == [ids.exposure_id, ids.exposure2_id]

    close(db)
end

@testset "sample annotations" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    ls = sample_labels(db, ids.sample_id)
    @test [(l.label, l.state, l.source) for l in ls] ==
          [("Calibrant", "suggested", "system"), ("Hot", "confirmed", "user")]
    @test [l.label for l in sample_labels(db, ids.sample2_id)] == ["Buffer"]

    @test subtraction_background(db, ids.sample_id) == ids.sample2_id
    @test subtraction_background(db, ids.sample2_id) === nothing

    @test bridged_samples(db, ids.sample_id) == [ids.sample2_id]
    @test bridged_samples(db, ids.sample2_id) == [ids.sample_id]
    @test isempty(bridged_samples(db, 999999))

    close(db)
end
