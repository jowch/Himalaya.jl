@testset "load_trace" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    # fixture filename is "s1"; default pattern -> "s1.dat" in analysis_dir (=dir)
    open(joinpath(dir, "s1.dat"), "w") do io
        println(io, "# q I sigma")
        println(io, "0.10 100.0 1.0")
        println(io, "0.20 50.0 0.5")
    end

    db = HimalayaDB.connect(dbpath)
    q, I, σ = load_trace(db, ids.exposure_id)
    @test q == [0.10, 0.20]
    @test I == [100.0, 50.0]
    @test σ == [1.0, 0.5]

    # explicit pattern= override
    open(joinpath(dir, "custom_s1.dat"), "w") do io
        println(io, "0.30 10.0 0.1")
    end
    qo, Io, σo = load_trace(db, ids.exposure_id; pattern="custom_{name}.dat")
    @test qo == [0.30]

    # missing exposure errors
    @test_throws ArgumentError load_trace(db, 999999)
    close(db)
end

@testset "load_trace resolves via experiment_id + strips per-frame suffix" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    db0 = HimalayaUI.open_db(dbpath)
    eid = HimalayaUI.create_experiment!(db0; name="acq", path=dir, data_dir=dir, analysis_dir=dir)
    # NULL sample_id (transient/ungrouped) + per-frame filename; only the shared
    # per-acquisition (suffix-stripped) .dat exists on disk.
    r = DBInterface.execute(db0,
        "INSERT INTO exposures (experiment_id, filename, kind, status) VALUES (?, 'acq_1_2', 'file', 'accepted')", [eid])
    xid = Int(DBInterface.lastrowid(r))
    close(db0)
    open(joinpath(dir, "acq.dat"), "w") do io; println(io, "0.11 9.0 0.2"); end
    db = HimalayaDB.connect(dbpath)
    q, I, sig = load_trace(db, xid)
    @test q == [0.11]
    close(db)
end
