import HimalayaUI

@testset "contract vs HimalayaUI" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)

    # HimalayaUI's own readers (read-write open) on the same DB
    uidb = HimalayaUI.open_db(dbpath)
    ui_peaks = HimalayaUI.get_peaks_for_exposure(uidb, ids.exposure_id)
    ui_indices = HimalayaUI.get_indices_for_exposure(uidb, ids.exposure_id)
    close(uidb)

    db = HimalayaDB.connect(dbpath)
    db_peaks = curated_peaks(db, ids.exposure_id)
    db_indices = index_candidates(db, ids.exposure_id)
    close(db)

    # Same rows, same order, same effective-peaks logic.
    @test [(p.q, p.source, p.excluded) for p in db_peaks] ==
          [(p.q, p.source, p.excluded) for p in ui_peaks]
    @test [(i.id, i.phase, i.score) for i in db_indices] ==
          [(i.id, i.phase, i.score) for i in ui_indices]
end
