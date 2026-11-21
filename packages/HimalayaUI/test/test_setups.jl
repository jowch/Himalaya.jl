using Test, Dates, SQLite, DBInterface, Tables, HimalayaUI
isdefined(@__MODULE__, :write_prp) || include("test_ingestion_core.jl")

function _write_setup!(root; stems, pipe_length_mm = 1700, h0 = 58.9)
    data = mkpath(joinpath(root, "data"))
    analysis = mkpath(joinpath(root, "analysis"))
    t0 = DateTime(2026, 4, 26, 23, 14, 8)
    for (i, stem) in enumerate(stems)
        write_prp(joinpath(data, "$stem.prp");
            timestamp = Dates.format(t0 + Second(19 * (i - 1)), "dd u yyyy HH:MM:SS"),
            pipe_length_mm, horizontal_position_mm = h0)
        write(joinpath(data, "$stem.tif"), "fake tif")
    end
    (data, analysis)
end

_count(db, sql, args...) = Int(only(Tables.rowtable(DBInterface.execute(db, sql, collect(args)))).c)

@testset "setups" begin
    @testset "one experiment scans each of its setups on its own" begin
        bt = mktempdir()
        stems = ["MO_S0001_0_001", "MO_S0002_0_001"]
        d1, a1 = _write_setup!(joinpath(bt, "setup1"); stems)
        d2, a2 = _write_setup!(joinpath(bt, "setup2"); stems, pipe_length_mm = 700)
        db = fresh_db()
        eid = HimalayaUI.create_experiment!(db; name = "bt", path = bt, data_dir = d1, analysis_dir = a1)
        s2 = HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = d2, analysis_dir = a2, position = 1)

        r = HimalayaUI.scan_and_group!(db, eid; analyze = false)
        @test r.status == :ok
        @test _count(db, "SELECT COUNT(*) AS c FROM exposures WHERE experiment_id = ?", eid) == 4
        for st in HimalayaUI.setups_of(db, eid)
            @test _count(db, "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ?", st.id) == 2
            @test _count(db, "SELECT COUNT(*) AS c FROM loads WHERE setup_id = ? AND load_index = 1", st.id) == 1
            @test _count(db, """SELECT COUNT(*) AS c FROM samples s
                                 WHERE s.setup_id = ? AND s.load_id NOT IN (SELECT id FROM loads WHERE setup_id = ?)""",
                         st.id, st.id) == 0
        end
        @test only(Tables.rowtable(DBInterface.execute(db,
            "SELECT folder FROM setups WHERE id = ?", [s2]))).folder == "setup2"

        again = HimalayaUI.scan_and_group!(db, eid; analyze = false)
        @test again.added_exposures == 0
        @test isempty(HimalayaUI._insertable_stems(db, eid))
    end

    @testset "an excluded setup is not scanned" begin
        bt = mktempdir()
        d1, a1 = _write_setup!(joinpath(bt, "setup1"); stems = ["A_S0001_0_001"])
        d2, a2 = _write_setup!(joinpath(bt, "setup2"); stems = ["B_S0001_0_001"])
        db = fresh_db()
        eid = HimalayaUI.create_experiment!(db; name = "bt", path = bt, data_dir = d1, analysis_dir = a1)
        s2 = HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = d2, analysis_dir = a2,
            position = 1, included = false)
        HimalayaUI.scan_and_group!(db, eid; analyze = false)
        @test _count(db, "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ?", s2) == 0
        @test _count(db, "SELECT COUNT(*) AS c FROM exposures WHERE experiment_id = ?", eid) == 1
        @test isempty(HimalayaUI._insertable_stems(db, eid))
    end

    @testset "geometry fills on the setup it was derived from" begin
        bt = mktempdir()
        d1, a1 = _write_setup!(joinpath(bt, "setup1"); stems = ["A_S0001_0_001"], pipe_length_mm = 1700)
        d2, a2 = _write_setup!(joinpath(bt, "setup2"); stems = ["B_S0001_0_001"], pipe_length_mm = 700)
        db = fresh_db()
        eid = HimalayaUI.create_experiment!(db; name = "bt", path = bt, data_dir = d1, analysis_dir = a1)
        HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = d2, analysis_dir = a2, position = 1)
        HimalayaUI.scan_and_group!(db, eid; analyze = false)
        fps = [st.flight_path_m for st in HimalayaUI.setups_of(db, eid)]
        @test fps ≈ [1.7, 0.7]
    end

    @testset "a trace is read through its exposure's own setup" begin
        bt = mktempdir()
        d1, a1 = _write_setup!(joinpath(bt, "setup1"); stems = String[])
        d2, a2 = _write_setup!(joinpath(bt, "setup2"); stems = String[])
        db = fresh_db()
        eid = HimalayaUI.create_experiment!(db; name = "bt", path = bt, data_dir = d1, analysis_dir = a1,
            integration_pattern = "{name}.dat")
        s2 = HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = d2, analysis_dir = a2,
            position = 1, integration_pattern = "{name}_tot.dat")
        x = HimalayaUI.create_exposure!(db; experiment_id = eid, setup_id = s2, filename = "f")
        @test HimalayaUI._resolve_analysis_dir(db, x) == a2
        @test HimalayaUI.integration_pattern_for_exposure(db, x) == "{name}_tot.dat"
    end

    @testset "migration gives each legacy experiment one setup and keeps every id" begin
        dir = mktempdir(); path = joinpath(dir, "legacy.db")
        legacy = SQLite.DB(path)
        HimalayaUI.create_schema!(legacy)
        add_legacy_experiment_columns!(legacy)
        for (c, decl) in (("integration_pattern", "TEXT"), ("flight_path_m", "REAL"),
                          ("flight_path_m_source", "TEXT"), ("q_units", "TEXT"))
            DBInterface.execute(legacy, "ALTER TABLE experiments ADD COLUMN $c $decl")
        end
        DBInterface.execute(legacy, """INSERT INTO experiments
            (id, name, path, data_dir, analysis_dir, integration_pattern, flight_path_m, flight_path_m_source, q_units)
            VALUES (7, '20261103/setup1', '/bt/setup1/data', '/bt/setup1/data', '/bt/setup1/analysis',
                    '{name}_tot.dat', 1.8095, 'setup', 'nm^-1')""")
        DBInterface.execute(legacy, """INSERT INTO experiments (id, name, path, data_dir, analysis_dir)
            VALUES (8, 'other', '/x', '/x', '/x/analysis')""")
        DBInterface.execute(legacy, "INSERT INTO loads (id, experiment_id, load_index) VALUES (70, 7, 1)")
        DBInterface.execute(legacy, "INSERT INTO samples (id, experiment_id, name, load_id) VALUES (700, 7, 'MO (S01P01)', 70)")
        DBInterface.execute(legacy, "INSERT INTO samples (id, experiment_id, name) VALUES (800, 8, 'B')")
        DBInterface.execute(legacy, """INSERT INTO exposures (id, experiment_id, sample_id, filename, load_id)
            VALUES (7000, 7, 700, 'MO_S0001_0_001', 70)""")
        SQLite.close(legacy)

        db = HimalayaUI.open_db(path)
        st7 = only(HimalayaUI.setups_of(db, 7))
        @test st7.folder == "setup1"
        @test st7.root_dir == "/bt/setup1"
        @test st7.data_dir == "/bt/setup1/data"
        @test st7.integration_pattern == "{name}_tot.dat"
        @test st7.flight_path_m == 1.8095 && st7.flight_path_m_source == "setup"
        @test st7.q_units == "nm^-1"
        @test only(HimalayaUI.setups_of(db, 8)).folder == "x"
        for (t, id) in (("loads", 70), ("samples", 700), ("exposures", 7000))
            @test only(Tables.rowtable(DBInterface.execute(db, "SELECT setup_id FROM $t WHERE id = ?", [id]))).setup_id == st7.id
        end
        @test only(Tables.rowtable(DBInterface.execute(db, "SELECT name, path FROM experiments WHERE id = 7"))) ==
              (name = "20261103/setup1", path = "/bt/setup1")
        for c in HimalayaUI.RETIRED_EXPERIMENT_COLUMNS
            @test !(c in HimalayaUI.cols_of(db, "experiments"))
        end
        @test HimalayaUI.integration_pattern_for_exposure(db, 7000) == "{name}_tot.dat"
        SQLite.close(db)

        db2 = HimalayaUI.open_db(path)
        @test length(HimalayaUI.setups_of(db2, 7)) == 1
        @test _count(db2, "SELECT COUNT(*) AS c FROM setups") == 2
        SQLite.close(db2)
    end
end
