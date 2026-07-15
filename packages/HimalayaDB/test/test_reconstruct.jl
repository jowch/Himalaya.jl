import Himalaya

@testset "reconstruct_index" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)
    db = HimalayaDB.connect(dbpath)

    idx = reconstruct_index(db, ids.index_id)
    @test idx isa Himalaya.Index
    @test Himalaya.phase(idx) === Himalaya.Pn3m
    @test Himalaya.basis(idx) == 0.10
    # two supporting peaks at ratio positions 1 and 2
    @test Himalaya.numpeaks(idx) == 2

    @test_throws ArgumentError reconstruct_index(db, 999999)
    close(db)
end
