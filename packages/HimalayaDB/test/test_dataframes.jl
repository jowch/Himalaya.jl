using DataFrames

@testset "dataframe extension" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    df = dataframe(curated_peaks(db, ids.exposure_id))
    @test df isa DataFrame
    @test nrow(df) == 4
    @test "source" in names(df)
    close(db)
end
