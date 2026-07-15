using Test, SQLite, DBInterface, Tables
using Himalaya

# Seed helper: one experiment/sample/exposure and one speculative index row
# written RAW (bypassing insert_speculative_index!) so we can simulate
# legacy-convention rows.
function _mig_fixture()
    tmp = mktempdir()
    db  = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="x")
    (; db, experiment_id = exp_id, exposure_id = e_id)
end

# Re-arm the one-shot migration so it can be driven directly against
# hand-seeded legacy rows (open_db already ran it on the fresh DB).
function _rearm!(db)
    DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?",
        [HimalayaUI.MIGRATION_SPECULATIVE_PEAK_DURABILITY])
end

@testset "migration: cubic speculative basis rescaled once" begin
    fx = _mig_fixture()
    # Legacy Pn3m row: un-normalized basis 0.1/√2 for an anchor at q=0.1.
    old_basis = 0.1 / first(Himalaya.phaseratios(Himalaya.Pn3m))
    res = DBInterface.execute(fx.db, """
        INSERT INTO indices (exposure_id, phase, basis, status, kind)
        VALUES (?, 'Pn3m', ?, 'candidate', 'speculative')""",
        [fx.exposure_id, old_basis])
    ix_id = Int(DBInterface.lastrowid(res))

    _rearm!(fx.db)
    HimalayaUI.migrate_speculative_peak_durability!(fx.db)

    b = Float64(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT basis FROM indices WHERE id = ?", [ix_id]))[1].basis)
    @test b ≈ 0.1 atol=1e-12

    # Sentinel: a second direct call must be a no-op (no double rescale).
    HimalayaUI.migrate_speculative_peak_durability!(fx.db)
    b2 = Float64(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT basis FROM indices WHERE id = ?", [ix_id]))[1].basis)
    @test b2 ≈ 0.1 atol=1e-12
end

@testset "migration: intents backfilled from surviving index_peaks" begin
    fx = _mig_fixture()
    res = DBInterface.execute(fx.db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.05)",
        [fx.exposure_id])
    p1 = Int(DBInterface.lastrowid(res))
    spec_id = HimalayaUI.insert_speculative_index!(fx.db, fx.exposure_id,
        Himalaya.Lamellar, Dict{Int,Int}(1 => p1))
    # Simulate a pre-fix DB: no intents rows yet (Task 4 writes them at
    # creation; delete to model the legacy state).
    DBInterface.execute(fx.db,
        "DELETE FROM speculative_peak_intents WHERE index_id = ?", [spec_id])

    _rearm!(fx.db)
    HimalayaUI.migrate_speculative_peak_durability!(fx.db)

    intents = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT ratio_position, q FROM speculative_peak_intents WHERE index_id = ?",
        [spec_id]))
    @test length(intents) == 1
    @test Int(intents[1].ratio_position) == 1
    @test Float64(intents[1].q) ≈ 0.05 atol=1e-12
end

@testset "migration: analysis_inputs_hash NULLed only for peak-less-speculative exposures" begin
    fx = _mig_fixture()
    # Exposure A: peak-less speculative (prod shape) → hash must be NULLed.
    res = DBInterface.execute(fx.db, """
        INSERT INTO indices (exposure_id, phase, basis, status, kind)
        VALUES (?, 'Square', 0.17, 'candidate', 'speculative')""", [fx.exposure_id])
    DBInterface.execute(fx.db,
        "UPDATE exposures SET analysis_inputs_hash = 'aaaa' WHERE id = ?",
        [fx.exposure_id])

    # Exposure B (same DB): healthy speculative → hash untouched.
    s2 = HimalayaUI.create_sample!(fx.db; experiment_id=fx.experiment_id, name="D2")
    e2 = HimalayaUI.create_exposure!(fx.db; sample_id=s2, filename="y")
    res = DBInterface.execute(fx.db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.05)", [e2])
    p1 = Int(DBInterface.lastrowid(res))
    HimalayaUI.insert_speculative_index!(fx.db, e2, Himalaya.Lamellar,
        Dict{Int,Int}(1 => p1))
    DBInterface.execute(fx.db,
        "UPDATE exposures SET analysis_inputs_hash = 'bbbb' WHERE id = ?", [e2])

    _rearm!(fx.db)
    HimalayaUI.migrate_speculative_peak_durability!(fx.db)

    ha = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT analysis_inputs_hash AS h FROM exposures WHERE id = ?", [fx.exposure_id]))[1].h
    hb = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT analysis_inputs_hash AS h FROM exposures WHERE id = ?", [e2]))[1].h
    @test ismissing(ha)          # gate reopened
    @test String(hb) == "bbbb"   # healthy exposure untouched
end
