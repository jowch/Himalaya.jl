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
