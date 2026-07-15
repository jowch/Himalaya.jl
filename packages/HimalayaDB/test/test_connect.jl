@testset "connect" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    build_fixture(dbpath, dir)

    db = HimalayaDB.connect(dbpath)
    @test db isa SQLite.DB
    # read works
    @test !isempty(Tables.rowtable(DBInterface.execute(db, "SELECT id FROM experiments")))
    # writes are blocked by query_only=ON
    @test_throws Exception DBInterface.execute(db, "CREATE TABLE t (x INTEGER)")
    close(db)

    # missing file errors, does not create it
    missing_path = joinpath(dir, "nope.db")
    @test_throws ArgumentError HimalayaDB.connect(missing_path)
    @test !isfile(missing_path)
end
