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
                    '{name}_tot.dat', 1.81, 'setup', 'nm^-1')""")
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
        @test st7.flight_path_m == 1.81 && st7.flight_path_m_source == "setup"
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
    @testset "upgrade from the pre-setups schema" begin
        path = joinpath(mktempdir(), "pre.db")
        old = SQLite.DB(path)
        for stmt in split(read(joinpath(@__DIR__, "fixtures", "pre_setups_schema.sql"), String), ";\n")
            isempty(strip(stmt)) || DBInterface.execute(old, stmt)
        end
        DBInterface.execute(old, """INSERT INTO experiments
            (id, name, path, data_dir, analysis_dir, integration_pattern,
             flight_path_m, flight_path_m_source, energy_kev, energy_kev_source)
            VALUES (1, 'one', '/bt/setup1/data', '/bt/setup1/data', '/bt/setup1/analysis', '{name}_tot.dat',
                    2.5, 'user', 9.0, NULL),
                   (2, 'two', '/elsewhere', '/bt/setup2/data', '/bt/setup2/analysis', NULL,
                    NULL, 'default', NULL, 'default')""")
        for (eid, lid) in ((1, 10), (2, 20))
            DBInterface.execute(old, "INSERT INTO loads (id, experiment_id, load_index) VALUES (?, ?, 1)", [lid, eid])
            DBInterface.execute(old, "INSERT INTO samples (id, experiment_id, name, load_id, slot_index) VALUES (?, ?, 'MO (S01P01)', ?, 1)",
                [lid + 1, eid, lid])
            DBInterface.execute(old, "INSERT INTO exposures (id, experiment_id, sample_id, load_id, filename) VALUES (?, ?, ?, ?, 'MO_S0001')",
                [lid + 2, eid, lid + 1, lid])
        end
        DBInterface.execute(old, "INSERT INTO sample_calls (sample_id, state) VALUES (11, 'form_factor')")
        close(old)

        db = HimalayaUI.open_db(path)
        @test isempty(Tables.rowtable(DBInterface.execute(db, "PRAGMA foreign_key_check")))
        @test only(Tables.rowtable(DBInterface.execute(db, "PRAGMA integrity_check"))).integrity_check == "ok"
        @test HimalayaUI.cols_of(db, "experiments") ==
              ["id", "name", "path", "created_at", "last_scanned_at", "ingest_status", "ingest_error", "description"]
        sts = Tables.rowtable(DBInterface.execute(db, "SELECT * FROM setups ORDER BY experiment_id"))
        @test [(r.experiment_id, r.folder, r.root_dir, r.data_dir) for r in sts] ==
              [(1, "setup1", "/bt/setup1", "/bt/setup1/data"), (2, "setup2", "/bt/setup2", "/bt/setup2/data")]
        @test (sts[1].flight_path_m, sts[1].flight_path_m_source) == (2.5, "user")
        @test (sts[1].energy_kev, sts[1].energy_kev_source) == (9.0, "default")
        @test sts[1].integration_pattern == "{name}_tot.dat"
        @test [r.path for r in Tables.rowtable(DBInterface.execute(db, "SELECT path FROM experiments ORDER BY id"))] ==
              ["/bt/setup1", "/elsewhere"]
        for t in ("loads", "samples", "exposures")
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT x.setup_id, st.experiment_id AS sx, x.experiment_id AS ex FROM $t x LEFT JOIN setups st ON st.id = x.setup_id"))
            @test length(rows) == 2
            @test all(r -> !ismissing(r.setup_id) && r.sx == r.ex, rows)
        end
        idx = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE name = 'exposures_unique_filename'"))).sql
        @test occursin("(setup_id, filename)", idx)
        @test _count(db, "SELECT COUNT(*) AS c FROM sample_calls WHERE sample_id = 11 AND state = 'form_factor'") == 1
        @test HimalayaUI._migrated(db, HimalayaUI.MIGRATION_SETUP_IDS)

        DBInterface.execute(db, "INSERT INTO exposures (experiment_id, filename) VALUES (1, 'orphan')")
        close(db)
        db = HimalayaUI.open_db(path)
        @test ismissing(only(Tables.rowtable(DBInterface.execute(db,
            "SELECT setup_id FROM exposures WHERE filename = 'orphan'"))).setup_id)
        @test _count(db, "SELECT COUNT(*) AS c FROM setups") == 2
    end
    @testset "a scan stores each frame's pipe length and an upgraded DB rereads it" begin
        bt = mktempdir()
        d1, a1 = _write_setup!(joinpath(bt, "setup1"); stems = ["MO_S0001_0_001"], pipe_length_mm = 700)
        db = fresh_db()
        eid = HimalayaUI.create_experiment!(db; name = "bt", path = bt, data_dir = d1, analysis_dir = a1)
        HimalayaUI.scan_and_group!(db, eid; analyze = false)
        HimalayaUI.suggest_labels_for_experiment!(db, eid)
        @test only(Tables.rowtable(DBInterface.execute(db, "SELECT pipe_length_m FROM exposures"))).pipe_length_m ≈ 0.7
        @test any(r -> r.label == "0.7 m", Tables.rowtable(DBInterface.execute(db, "SELECT label FROM sample_labels")))

        DBInterface.execute(db, "UPDATE exposures SET pipe_length_m = NULL")
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_EXPOSURES_PIPE_LENGTH])
        HimalayaUI.migrate_exposures_pipe_length!(db)
        @test only(Tables.rowtable(DBInterface.execute(db, "SELECT monitors_read FROM exposures"))).monitors_read == 0
        @test HimalayaUI.backfill_monitors!(db, eid) == 1
        @test only(Tables.rowtable(DBInterface.execute(db, "SELECT pipe_length_m FROM exposures"))).pipe_length_m ≈ 0.7
    end
end

@testset "trace normalisation is read from the setup's integ.mpp" begin
    dir = mktempdir()
    mpp(text) = (p = joinpath(dir, "integ.mpp"); write(p, text); p)
    @test HimalayaUI.parse_integ_normalisation(mpp("-o x\n-i yes i2 2.5e6\n")) == (by = "i2", base = 2.5e6)
    @test HimalayaUI.parse_integ_normalisation(mpp("  -i yes I1 100\n")) == (by = "i1", base = 100.0)
    @test HimalayaUI.parse_integ_normalisation(mpp("-o x\n-s 1\n")) == (by = "none", base = nothing)
    @test HimalayaUI.parse_integ_normalisation(mpp("-i no\n")) == (by = "none", base = nothing)
    @test HimalayaUI.parse_integ_normalisation(mpp("-i yes ion3 1\n")).by == "unknown"
    @test HimalayaUI.parse_integ_normalisation(mpp("-i yes\n")).by == "unknown"

    root = mktempdir()
    analysis = joinpath(root, "setup1", "analysis", "automatic_analysis")
    mkpath(analysis)
    @test HimalayaUI.trace_normalisation(analysis).by == "unknown"
    write(joinpath(root, "setup1", "analysis", "integ.mpp"), "-i yes i2 1e6\n")
    @test HimalayaUI.trace_normalisation(analysis) == (by = "i2", base = 1e6)

    db = HimalayaUI.open_db(joinpath(mktempdir(), "h.db"))
    other = joinpath(root, "setup2", "analysis")
    mkpath(other)
    eid = HimalayaUI.create_experiment!(db; name = "norm", path = root,
        data_dir = joinpath(root, "setup1", "data"), analysis_dir = analysis)
    sid2 = HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = joinpath(root, "setup2", "data"),
        analysis_dir = other, position = 1)
    stored() = Dict(Int(r.id) => (r.trace_normalisation, r.normalisation_base)
                    for r in Tables.rowtable(DBInterface.execute(db, "SELECT id, trace_normalisation, normalisation_base FROM setups")))
    @test HimalayaUI.read_trace_normalisations!(db; experiment_id = eid) == 1
    @test isequal(stored()[sid2], ("unknown", missing))
    @test any(v -> v == ("i2", 1e6), values(stored()))

    write(joinpath(other, "integ.mpp"), "-o x\n")
    @test HimalayaUI.read_trace_normalisations!(db) == 1
    @test isequal(stored()[sid2], ("none", missing))
    write(joinpath(other, "integ.mpp"), "-i yes i2 5\n")
    @test HimalayaUI.read_trace_normalisations!(db) == 0
    @test isequal(stored()[sid2], ("none", missing))
end
