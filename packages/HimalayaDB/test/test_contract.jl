import HimalayaUI

@testset "contract vs HimalayaUI" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    ids = build_fixture(dbpath, dir)

    uidb = HimalayaUI.open_db(dbpath)
    ui_peaks = HimalayaUI.get_peaks_for_exposure(uidb, ids.exposure_id)
    ui_indices = HimalayaUI.get_indices_for_exposure(uidb, ids.exposure_id)
    ui_calls = Dict(sid => HimalayaUI._sample_call_body(uidb, sid) for sid in (ids.sample_id, ids.sample2_id))
    ui_labels = HimalayaUI.labels_for_experiment(uidb, ids.experiment_id)
    ui_setup = HimalayaUI.exposure_setup(uidb, ids.exposure_id)
    ui_pattern = HimalayaUI.integration_pattern_for(ui_setup)
    close(uidb)

    db = HimalayaDB.connect(dbpath)
    db_peaks = curated_peaks(db, ids.exposure_id)
    db_indices = index_candidates(db, ids.exposure_id)
    db_calls = Dict(sid => phase_call(db, sid) for sid in (ids.sample_id, ids.sample2_id))
    db_labels = [l for sid in (ids.sample_id, ids.sample2_id) for l in sample_labels(db, sid)]
    db_setup = only(filter(s -> s.id == ui_setup.id, setups(db)))
    close(db)

    @test [(p.q, p.source, p.excluded) for p in db_peaks] ==
          [(p.q, p.source, p.excluded) for p in ui_peaks]
    @test [(i.id, i.phase, i.score) for i in db_indices] ==
          [(i.id, i.phase, i.score) for i in ui_indices]
    for sid in (ids.sample_id, ids.sample2_id)
        ui, ours = ui_calls[sid], db_calls[sid]
        @test ours.state == ui[:state]
        @test ours.phase == ui[:phase]
        @test ours.call_phases == ui[:call_phases]
        @test [(m.index_id, m.exposure_id, m.phase) for m in ours.members] ==
              [(m[:index_id], m[:exposure_id], m[:phase]) for m in ui[:members]]
    end
    @test sort([(l.id, l.label, l.state, l.source) for l in db_labels]) ==
          sort([(l[:id], l[:label], l[:state], l[:source]) for l in ui_labels])
    @test db_setup.analysis_dir == ui_setup.analysis_dir
    @test db_setup.integration_pattern == ui_pattern
end
