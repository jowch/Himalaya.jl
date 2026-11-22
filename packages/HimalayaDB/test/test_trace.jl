@testset "load_trace" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    open(joinpath(dir, "s1_tot.dat"), "w") do io
        println(io, "# q I sigma")
        println(io, "0.10 100.0 1.0")
        println(io, "0.20 50.0 0.5")
    end

    db = HimalayaDB.connect(dbpath)
    q, I, σ = load_trace(db, ids.exposure_id)
    @test q == [0.10, 0.20]
    @test I == [100.0, 50.0]
    @test σ == [1.0, 0.5]

    open(joinpath(dir, "custom_s1.dat"), "w") do io
        println(io, "0.30 10.0 0.1")
    end
    qo, Io, σo = load_trace(db, ids.exposure_id; pattern="custom_{name}.dat")
    @test qo == [0.30]

    @test_throws ArgumentError load_trace(db, 999999)
    close(db)
end

@testset "load_trace reads the exposure's own setup and strips the _0_001-style suffix" begin
    dir = mktempdir()
    other = mkpath(joinpath(dir, "setup2", "analysis"))
    dbpath = joinpath(dir, "himalaya.db")
    db0 = HimalayaUI.open_db(dbpath)
    eid = HimalayaUI.create_experiment!(db0; name="acq", path=dir, data_dir=dir, analysis_dir=dir)
    sid2 = HimalayaUI.create_setup!(db0; experiment_id=eid, data_dir=other, analysis_dir=other,
        position=1)
    xid = HimalayaUI.create_exposure!(db0; experiment_id=eid, setup_id=sid2,
        filename="acq_1_2", kind="file", status="accepted")
    orphan = Int(DBInterface.lastrowid(DBInterface.execute(db0,
        "INSERT INTO exposures (experiment_id, filename, kind, status) VALUES (?, 'lost', 'file', 'accepted')", [eid])))
    close(db0)
    open(joinpath(dir, "acq.dat"), "w") do io; println(io, "9.99 9.0 0.2"); end
    open(joinpath(other, "acq.dat"), "w") do io; println(io, "0.11 9.0 0.2"); end
    db = HimalayaDB.connect(dbpath)
    q, I, sig = load_trace(db, xid)
    @test q == [0.11]
    @test_throws ArgumentError load_trace(db, orphan)
    close(db)
end
