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
    # three supporting peaks: auto at positions 1,2 + curation-backed add at 3
    @test Himalaya.numpeaks(idx) == 3

    # curation-supported peak reconstructs with sharpness 0 (peak_curations stores
    # none) — this is why Himalaya.score of a reconstructed index with an add-peak
    # diverges from the stored score. Pinned so the fidelity boundary is conscious.
    @test idx.peaks[3] == 0.20
    @test idx.sharpness[3] == 0.0
    @test idx.sharpness[1] == 1.0   # auto peak keeps its sampled sharpness

    @test_throws ArgumentError reconstruct_index(db, 999999)
    close(db)
end
