using Test, SQLite, DBInterface, Tables

function _ratio_fixture()
    tmp = mktempdir()
    path = joinpath(tmp, "himalaya.db")
    db = HimalayaUI.open_db(path)
    exp_id = HimalayaUI.init_experiment!(db; path = tmp,
        data_dir = joinpath(tmp, "data"), analysis_dir = joinpath(tmp, "analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "R1")
    e_id = HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = s_id, filename = "x")
    ix = Int(DBInterface.lastrowid(DBInterface.execute(db, """
        INSERT INTO indices (exposure_id, phase, basis, status, kind)
        VALUES (?, 'Pn3m', 0.1, 'candidate', 'auto')""", [e_id])))
    (; db, path, exposure_id = e_id, index_id = ix)
end

function _seed_peak!(db, exposure_id, index_id, pos)
    pk = Int(DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO auto_peaks (exposure_id, q, intensity, prominence, sharpness) VALUES (?, ?, 1.0, 1.0, 1.0)",
        [exposure_id, 0.01 * Int(pos)])))
    DBInterface.execute(db,
        "INSERT INTO index_peaks (index_id, peak_id, ratio_position) VALUES (?, ?, ?)",
        Any[index_id, pk, pos])
end

_ratio_rows(db) = [(r.c, Int(r.v)) for r in Tables.rowtable(DBInterface.execute(db,
    "SELECT typeof(ratio_position) AS c, ratio_position AS v FROM index_peaks ORDER BY ratio_position"))]

_rearm_ratio!(db) = DBInterface.execute(db,
    "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_RATIO_POSITION_INT])

@testset "ratio_position migration: serialized UInt8 blobs become integers on open" begin
    fx = _ratio_fixture()
    _seed_peak!(fx.db, fx.exposure_id, fx.index_id, UInt8(3))
    _seed_peak!(fx.db, fx.exposure_id, fx.index_id, UInt8(200))
    _seed_peak!(fx.db, fx.exposure_id, fx.index_id, 5)
    @test sort(first.(_ratio_rows(fx.db))) == ["blob", "blob", "integer"]

    _rearm_ratio!(fx.db)
    close(fx.db)
    db = HimalayaUI.open_db(fx.path)
    @test _ratio_rows(db) == [("integer", 3), ("integer", 5), ("integer", 200)]
    @test HimalayaUI._migrated(db, HimalayaUI.MIGRATION_RATIO_POSITION_INT)
end

@testset "ratio_position migration: a blob it cannot decode is left alone and the sentinel stays open" begin
    fx = _ratio_fixture()
    pk = Int(DBInterface.lastrowid(DBInterface.execute(fx.db,
        "INSERT INTO auto_peaks (exposure_id, q, intensity, prominence, sharpness) VALUES (?, 0.02, 1.0, 1.0, 1.0)",
        [fx.exposure_id])))
    DBInterface.execute(fx.db,
        "INSERT INTO index_peaks (index_id, peak_id, ratio_position) VALUES (?, ?, X'0102')",
        [fx.index_id, pk])
    _seed_peak!(fx.db, fx.exposure_id, fx.index_id, UInt8(4))

    _rearm_ratio!(fx.db)
    @test_logs (:warn,) match_mode = :any HimalayaUI.migrate_index_peaks_ratio_int!(fx.db)
    classes = [r.c for r in Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT typeof(ratio_position) AS c FROM index_peaks ORDER BY rowid"))]
    @test classes == ["blob", "integer"]
    @test !HimalayaUI._migrated(fx.db, HimalayaUI.MIGRATION_RATIO_POSITION_INT)
end
