import HimalayaUI
using DBInterface

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
    ui_subs = HimalayaUI.compares_for_experiment(uidb, ids.experiment_id)[:subtractions]
    ui_bridges = HimalayaUI.bridge_groups_for_experiment(uidb, ids.experiment_id)
    ui_sample_ids = [Int(s.id) for s in HimalayaUI.get_samples(uidb, ids.experiment_id)]
    close(uidb)

    db = HimalayaDB.connect(dbpath)
    db_sample_ids = [s.id for s in samples(db; experiment=ids.experiment_id)]
    db_peaks = curated_peaks(db, ids.exposure_id)
    db_indices = index_candidates(db, ids.exposure_id)
    db_calls = Dict(sid => phase_call(db, sid) for sid in (ids.sample_id, ids.sample2_id))
    db_labels = [l for sid in (ids.sample_id, ids.sample2_id) for l in sample_labels(db, sid)]
    db_setup = only(filter(s -> s.id == ui_setup.id, setups(db)))
    all_ids = (ids.sample_id, ids.sample2_id)
    db_subs = [(sid, subtraction_background(db, sid)) for sid in all_ids]
    db_bridges = Dict(sid => bridged_samples(db, sid) for sid in all_ids)
    db_all_labels = [l for sid in all_ids for l in sample_labels(db, sid)]
    close(db)

    @test db_sample_ids == ui_sample_ids
    @test ids.empty_id in db_sample_ids
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
    @test db_setup.mode == ui_setup.mode == "capillary"
    @test db_setup.mode_source == ui_setup.mode_source == "user"
    @test HimalayaDB.integration_pattern_for(db_setup.integration_pattern) == ui_pattern

    @test sort([(sid, ref) for (sid, ref) in db_subs if ref !== nothing]) ==
          sort([(s[:sample_id], s[:ref_sample_id]) for s in ui_subs])
    ui_group = Dict(b[:sample_id] => b[:group_id] for b in ui_bridges)
    for sid in all_ids
        expected = haskey(ui_group, sid) ?
            sort([b[:sample_id] for b in ui_bridges if b[:group_id] == ui_group[sid] && b[:sample_id] != sid]) : Int[]
        @test db_bridges[sid] == expected
    end
    @test sort([l.id for l in db_all_labels]) == sort([l[:id] for l in ui_labels])

    @test HimalayaDB.DEFAULT_INTEGRATION_PATTERN == HimalayaUI.DEFAULT_INTEGRATION_PATTERN
    for column in (missing, nothing, "", "frames.dat", "{name}_tot.dat")
        @test HimalayaDB.integration_pattern_for(column) ==
              HimalayaUI.integration_pattern_for((; integration_pattern = column))
    end
end

@testset "load_trace falls back to HimalayaUI's default pattern" begin
    dir = mktempdir()
    dbpath = joinpath(dir, "himalaya.db")
    uidb = HimalayaUI.open_db(dbpath)
    eid = HimalayaUI.create_experiment!(uidb; name="nopattern", path=dir, data_dir=dir, analysis_dir=dir)
    xid = HimalayaUI.create_exposure!(uidb; experiment_id=eid, filename="f1", kind="file", status="accepted")
    DBInterface.execute(uidb, "UPDATE setups SET integration_pattern = 'no-placeholder.dat'")
    expected = joinpath(dir, replace(HimalayaUI.integration_pattern_for_exposure(uidb, xid), "{name}" => "f1"))
    close(uidb)
    open(expected, "w") do io; println(io, "0.12 3.0 0.1"); end
    db = HimalayaDB.connect(dbpath)
    @test first(load_trace(db, xid)) == [0.12]
    close(db)
end
