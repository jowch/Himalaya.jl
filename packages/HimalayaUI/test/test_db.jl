using Test, SQLite, DBInterface, Tables, Logging
using HimalayaUI: create_schema!, migrate_schema!, create_experiment!, create_sample!,
                  create_exposure!, get_experiment, get_samples, get_exposures,
                  migrate_r2_split_peaks!

@testset "db schema" begin
    db = SQLite.DB()  # in-memory
    create_schema!(db)

    tables = Set(r[1] for r in DBInterface.execute(db,
        "SELECT name FROM sqlite_master WHERE type='table'"))

    for t in ["users", "experiments", "samples",
              "exposures", "exposure_sources",
              "indices", "index_peaks",
              "auto_peaks", "peak_curations",
              "index_groups", "index_group_members", "user_actions"]
        @test t in tables
    end
    # peaks is removed from SCHEMA in R2.2 — fresh DBs must NOT have it
    @test "peaks" ∉ tables
end

@testset "open_db sets journal_mode = WAL" begin
    # WAL lets concurrent readers proceed alongside one writer, which is
    # load-bearing for the parallel-request perf fix (#115). Without WAL,
    # SQLite's default rollback journal serializes all readers behind any
    # in-flight writer.
    mktempdir() do dir
        db_path = joinpath(dir, "h.db")
        db = HimalayaUI.open_db(db_path)
        rows = Tables.rowtable(DBInterface.execute(db, "PRAGMA journal_mode"))
        @test lowercase(String(rows[1].journal_mode)) == "wal"
    end
end

@testset "exposures schema migration" begin
    db = SQLite.DB()
    create_schema!(db)
    migrate_schema!(db)
    # Verify columns exist by inserting and reading back
    exp_id  = create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
    samp_id = create_sample!(db; experiment_id=exp_id, name = "S")
    DBInterface.execute(db,
        "INSERT INTO exposures (sample_id, filename, kind, status, image_path)
         VALUES (?, 'test.dat', 'file', 'accepted', '/tmp/test.tiff')", [samp_id])
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT status, image_path FROM exposures"))
    @test length(rows) == 1
    @test rows[1].status == "accepted"
    @test rows[1].image_path == "/tmp/test.tiff"
    # Calling migrate_schema! again is idempotent
    migrate_schema!(db)
end

@testset "db CRUD" begin
    db = SQLite.DB()
    create_schema!(db)

    exp_id = create_experiment!(db;
        name        = "TestRun",
        path        = "/data/exp1",
        data_dir    = "/data/exp1/data",
        analysis_dir = "/data/exp1/analysis/automatic_analysis")
    @test exp_id == 1

    exp = get_experiment(db, exp_id)
    @test exp.name == "TestRun"
    @test exp.path == "/data/exp1"

    s_id = create_sample!(db; experiment_id = exp_id, name = "D1")
    @test s_id == 1

    samples = get_samples(db, exp_id)
    @test length(samples) == 1
    @test first(samples).name == "D1"

    e_id = create_exposure!(db; experiment_id = exp_id, sample_id = s_id, filename = "JC001", kind = "file")
    @test e_id == 1

    exposures = get_exposures(db, s_id)
    @test length(exposures) == 1
    @test first(exposures).filename == "JC001"
end

@testset "setup configuration lives on setups, not experiments" begin
    db = SQLite.DB()
    create_schema!(db)
    exp_cols = HimalayaUI.cols_of(db, "experiments")
    for c in ("data_dir", "analysis_dir", "integration_pattern", "energy_kev", "flight_path_m", "config")
        @test !(c in exp_cols)
    end
    setup_cols = HimalayaUI.cols_of(db, "setups")
    for c in ("folder", "display_name", "included", "root_dir", "data_dir", "analysis_dir",
              "integration_pattern", "energy_kev", "flight_path_m", "q_units")
        @test c in setup_cols
    end
end

@testset "create_experiment! stores its setup" begin
    db = SQLite.DB()
    create_schema!(db)
    exp_id = create_experiment!(db;
        name = "ConfigTest",
        path = "/tmp/x",
        data_dir = "/tmp/x/data",
        analysis_dir = "/tmp/x/analysis",
        energy_kev = 12.0,
        flight_path_m = 2.5)
    setups = HimalayaUI.setups_of(db, exp_id)
    @test length(setups) == 1
    st = only(setups)
    @test st.folder == "x"
    @test st.root_dir == "/tmp/x"
    @test st.data_dir == "/tmp/x/data"
    @test st.energy_kev == 12.0
    @test st.flight_path_m == 2.5
    @test ismissing(st.display_name)
    @test st.included == 1
end

@testset "migrate_schema! moves a legacy experiment's folders into one setup" begin
    db = SQLite.DB()
    DBInterface.execute(db, """
        CREATE TABLE experiments (
            id INTEGER PRIMARY KEY,
            name TEXT,
            path TEXT NOT NULL,
            data_dir TEXT NOT NULL,
            analysis_dir TEXT NOT NULL,
            manifest_path TEXT,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP
        )""")
    DBInterface.execute(db,
        "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'old', '/bt/setup1/data', '/bt/setup1/data', '/bt/setup1/analysis')")
    create_schema!(db)
    migrate_schema!(db)
    @test !("data_dir" in HimalayaUI.cols_of(db, "experiments"))
    st = only(HimalayaUI.setups_of(db, 1))
    @test st.folder == "setup1"
    @test st.root_dir == "/bt/setup1"
    @test st.data_dir == "/bt/setup1/data"
    @test st.analysis_dir == "/bt/setup1/analysis"
    @test only(Tables.rowtable(DBInterface.execute(db, "SELECT path FROM experiments WHERE id = 1"))).path == "/bt/setup1"
end

@testset "index_groups partial unique constraint on custom" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

        DBInterface.execute(db,
            "INSERT INTO index_groups (exposure_id, kind) VALUES (?, 'custom')", [e_id])
        @test_throws SQLite.SQLiteException DBInterface.execute(db,
            "INSERT INTO index_groups (exposure_id, kind) VALUES (?, 'custom')", [e_id])

        # Auto groups can multiply (only 'custom' is unique per exposure).
        DBInterface.execute(db,
            "INSERT INTO index_groups (exposure_id, kind) VALUES (?, 'auto')", [e_id])
        DBInterface.execute(db,
            "INSERT INTO index_groups (exposure_id, kind) VALUES (?, 'auto')", [e_id])
    end
end

@testset "migrate_r2_split_peaks! on fresh R2.2 DB is no-op" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        # On a fresh R2.2 DB, `peaks` does not exist in the schema. The migration
        # sentinel (peaks_exists check) should return early immediately.
        auto_before = Tables.rowtable(DBInterface.execute(db, "SELECT COUNT(*) AS n FROM auto_peaks"))[1].n
        HimalayaUI.migrate_r2_split_peaks!(db)
        auto_after = Tables.rowtable(DBInterface.execute(db, "SELECT COUNT(*) AS n FROM auto_peaks"))[1].n
        # Nothing changed (still a no-op).
        @test Int(auto_before) == Int(auto_after)
        # peaks does not exist on R2.2+ fresh DBs.
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM sqlite_master WHERE name = 'peaks'")))
    end
end

@testset "migrate_r2_split_peaks! on legacy DB backfills" begin
    mktempdir() do dir
        db_path = joinpath(dir, "h.db")
        db = HimalayaUI.open_db(db_path)
        # Simulate a pre-R2.1 legacy DB by: dropping the new R2.1 tables,
        # and creating a legacy `peaks` table with data (data makes sqlite_sequence
        # aware of it, which is how the migration sentinel detects "this is a
        # legacy DB that was actually used, not a fresh R2.1 DB").
        # Note: after R2.2, `peaks` is not in SCHEMA so it doesn't exist on
        # fresh DBs — no need to drop it; just drop the R2.1 destination tables.
        DBInterface.execute(db, "DROP TABLE auto_peaks")
        DBInterface.execute(db, "DROP TABLE peak_curations")
        DBInterface.execute(db, """
            CREATE TABLE peaks (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                exposure_id INTEGER, q REAL, intensity REAL, prominence REAL,
                sharpness REAL, source TEXT DEFAULT 'auto', excluded INTEGER DEFAULT 0
            )
        """)
        # Set up an exposure
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)
        # Three peaks: one auto kept, one auto excluded, one manual.
        # The INSERTs populate sqlite_sequence["peaks"], marking this as a
        # "used" DB (not a fresh R2.1 DB waiting for first analysis).
        DBInterface.execute(db, "INSERT INTO peaks (exposure_id, q, source, excluded) VALUES (?, 0.10, 'auto', 0)", [e_id])
        DBInterface.execute(db, "INSERT INTO peaks (exposure_id, q, source, excluded) VALUES (?, 0.15, 'auto', 1)", [e_id])
        DBInterface.execute(db, "INSERT INTO peaks (exposure_id, q, source, excluded) VALUES (?, 0.20, 'manual', 0)", [e_id])

        # Run migration directly (simulates what open_db → migrate_schema! does
        # on a legacy DB: create_schema! creates the destination tables, then
        # migrate_r2_split_peaks! backfills them).
        HimalayaUI.create_schema!(db)
        HimalayaUI.migrate_r2_split_peaks!(db)

        autos = Tables.rowtable(DBInterface.execute(db, "SELECT * FROM auto_peaks"))
        @test length(autos) == 2  # both auto peaks (excluded ones still in auto_peaks)
        curs  = Tables.rowtable(DBInterface.execute(db, "SELECT * FROM peak_curations"))
        @test length(curs) == 2  # one exclude, one add
        @test Set([String(c.kind) for c in curs]) == Set(["exclude", "add"])
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM sqlite_master WHERE name = 'peaks'")))
    end
end

@testset "migrate_r2_split_peaks! is idempotent" begin
    mktempdir() do dir
        path = joinpath(dir, "h.db")
        db = HimalayaUI.open_db(path)
        # First open: record table list.
        tables1 = Set(String(r.name) for r in Tables.rowtable(
            DBInterface.execute(db, "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")))
        # Second open — must not error or duplicate any table.
        @test_nowarn HimalayaUI.open_db(path)
        db2 = HimalayaUI.open_db(path)
        tables2 = Set(String(r.name) for r in Tables.rowtable(
            DBInterface.execute(db2, "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")))
        @test tables1 == tables2
    end
end

@testset "migrate_r2_split_peaks! preserves index_peaks for speculatives anchored on manual peaks" begin
    # The bug this guards against: a user-built speculative index that
    # anchored on a manual peak loses its anchor when migration drops the
    # manual peak. The repoint step in the migration prevents this.
    mktempdir() do dir
        db_path = joinpath(dir, "h.db")
        db = HimalayaUI.open_db(db_path)
        # Construct legacy schema by hand and seed with a speculative
        # referencing a manual peak. Drop the R2.1 destination tables and
        # create a legacy `peaks` table (after R2.2, peaks is not in SCHEMA
        # so it doesn't exist on fresh DBs).
        DBInterface.execute(db, "DROP TABLE auto_peaks")
        DBInterface.execute(db, "DROP TABLE peak_curations")
        DBInterface.execute(db, """
            CREATE TABLE peaks (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                exposure_id INTEGER, q REAL, intensity REAL, prominence REAL,
                sharpness REAL, source TEXT DEFAULT 'auto', excluded INTEGER DEFAULT 0
            )
        """)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)
        # Manual peak that the speculative will anchor on.
        res = DBInterface.execute(db,
            "INSERT INTO peaks (exposure_id, q, source) VALUES (?, 0.20, 'manual')", [e_id])
        manual_peak_id = Int(DBInterface.lastrowid(res))
        # Pseudo-speculative index referencing the manual peak.
        res = DBInterface.execute(db,
            "INSERT INTO indices (exposure_id, phase, basis, kind) VALUES (?, 'Pn3m', 0.20, 'speculative')",
            [e_id])
        ix_id = Int(DBInterface.lastrowid(res))
        DBInterface.execute(db,
            "INSERT INTO index_peaks (index_id, peak_id, ratio_position, residual) VALUES (?, ?, 1, 0.0)",
            [ix_id, manual_peak_id])

        # Run migration.
        HimalayaUI.create_schema!(db)
        HimalayaUI.migrate_r2_split_peaks!(db)

        # Assert the speculative still has its anchor — repointed at the
        # new peak_curations.id, not dangling or deleted.
        ip = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM index_peaks WHERE index_id = ?", [ix_id]))
        @test length(ip) == 1
        @test String(ip[1].peak_kind) == "curation"
        # The new peak_id should be a valid peak_curations row of the right shape.
        curations = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, kind, q FROM peak_curations WHERE id = ?", [Int(ip[1].peak_id)]))
        @test length(curations) == 1
        @test String(curations[1].kind) == "add"
        @test Float64(curations[1].q) ≈ 0.20
    end
end

@testset "open_db rejects pre-existing duplicate custom index_groups" begin
    mktempdir() do dir
        db_path = joinpath(dir, "h.db")
        # Build a legacy DB with the partial unique index missing AND a
        # pre-existing duplicate-custom-group row (the multiplayer-era TOCTOU
        # outcome the index now prevents).
        legacy = SQLite.DB(db_path)
        DBInterface.execute(legacy, "PRAGMA foreign_keys = ON")
        # Minimal subset of the schema needed to seed the duplicate.
        for stmt in split(HimalayaUI.SCHEMA, ";")
            s = strip(stmt)
            isempty(s) && continue
            DBInterface.execute(legacy, s)
        end
        # Drop the partial unique index that open_db will try to add.
        DBInterface.execute(legacy,
            "DROP INDEX IF EXISTS idx_one_custom_group_per_exposure")
        exp_id = HimalayaUI.create_experiment!(legacy; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(legacy; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(legacy; experiment_id=exp_id, sample_id=s_id)
        DBInterface.execute(legacy,
            "INSERT INTO index_groups (exposure_id, kind) VALUES (?, 'custom')", [e_id])
        DBInterface.execute(legacy,
            "INSERT INTO index_groups (exposure_id, kind) VALUES (?, 'custom')", [e_id])
        SQLite.DBInterface.close!(legacy)

        @test_throws ErrorException HimalayaUI.open_db(db_path)
    end
end

@testset "open_db heals corrupted _migrate_old_* FK references on populated tables" begin
    mktempdir() do tmp
        db_path = joinpath(tmp, "corrupted.db")
        db = SQLite.DB(db_path)
        DBInterface.execute(db, """
            CREATE TABLE samples (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                experiment_id INTEGER,
                label TEXT, name TEXT, notes TEXT
            )""")
        DBInterface.execute(db, """
            CREATE TABLE users (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                username TEXT UNIQUE NOT NULL,
                first_name TEXT, last_name TEXT
            )""")
        DBInterface.execute(db, """
            CREATE TABLE sample_bridges (
                sample_id INTEGER PRIMARY KEY REFERENCES "_migrate_old_samples"(id) ON DELETE CASCADE,
                group_id  INTEGER NOT NULL
            )""")
        DBInterface.execute(db, "INSERT INTO samples (label) VALUES ('S1')")
        SQLite.close(db)

        verify = SQLite.DB(db_path)
        broken_before = Tables.rowtable(DBInterface.execute(verify,
            "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE '%_migrate_old_%'"))
        @test length(broken_before) == 1
        SQLite.close(verify)

        healed = HimalayaUI.open_db(db_path)
        broken_after = Tables.rowtable(DBInterface.execute(healed,
            "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE '%_migrate_old_%'"))
        @test isempty(broken_after)

        DBInterface.execute(healed,
            "INSERT INTO sample_bridges (sample_id, group_id) VALUES (?, ?)", [1, 1])
        @test first(Tables.rowtable(DBInterface.execute(healed,
            "SELECT COUNT(*) AS c FROM sample_bridges"))).c == 1
    end
end

@testset "_fix_fk_references_after_autoincrement_migration! heals bare-form (no quotes)" begin
    # The smoke regression covered the QUOTED form ("_migrate_old_samples").
    # The other branch is the BARE form (_migrate_old_samples without quotes)
    # — SQLite accepts unquoted table names in CREATE statements, and that's
    # the more interesting branch since it'd be missed by a quoted-only fix
    # (review issue #23). Tested via the heal helper directly to avoid
    # tangling with create_schema! migrations.
    mktempdir() do tmp
        db_path = joinpath(tmp, "corrupted-bare.db")
        db = SQLite.DB(db_path)
        DBInterface.execute(db,
            "CREATE TABLE exposures (id INTEGER PRIMARY KEY AUTOINCREMENT, x TEXT)")
        # Bare form — no quotes around the FK target.
        DBInterface.execute(db, """
            CREATE TABLE refs_bare_exposure (
                id INTEGER PRIMARY KEY,
                exposure_id INTEGER REFERENCES _migrate_old_exposures(id)
            )""")

        broken_before = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE '%_migrate_old_%'"))
        @test length(broken_before) == 1

        HimalayaUI._fix_fk_references_after_autoincrement_migration!(db)

        broken_after = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE '%_migrate_old_%'"))
        @test isempty(broken_after)
        SQLite.close(db)
    end
end

@testset "_fix_fk_references_after_autoincrement_migration! heals multiple entities" begin
    # Smoke regression covered only `samples`. Real corruption has multiple
    # entities at once. This test exercises the heal loop directly (no
    # open_db / create_schema! interference) so we can include corrupted
    # FKs targeting `samples`, `exposures`, AND `experiments` in one DB
    # without colliding with production schema migrations (review issue #23).
    mktempdir() do tmp
        db_path = joinpath(tmp, "corrupted-multi.db")
        db = SQLite.DB(db_path)
        # Three entity targets with distinct corrupted-FK references to each.
        # Use synthetic referent table names so create_schema! is irrelevant.
        DBInterface.execute(db,
            "CREATE TABLE samples (id INTEGER PRIMARY KEY AUTOINCREMENT, x TEXT)")
        DBInterface.execute(db,
            "CREATE TABLE exposures (id INTEGER PRIMARY KEY AUTOINCREMENT, x TEXT)")
        DBInterface.execute(db,
            "CREATE TABLE experiments (id INTEGER PRIMARY KEY AUTOINCREMENT, x TEXT)")
        DBInterface.execute(db, """
            CREATE TABLE refs_to_samples (
                id INTEGER PRIMARY KEY,
                sample_id INTEGER REFERENCES "_migrate_old_samples"(id)
            )""")
        DBInterface.execute(db, """
            CREATE TABLE refs_to_exposures (
                id INTEGER PRIMARY KEY,
                exposure_id INTEGER REFERENCES _migrate_old_exposures(id)
            )""")
        DBInterface.execute(db, """
            CREATE TABLE refs_to_experiments (
                id INTEGER PRIMARY KEY,
                experiment_id INTEGER REFERENCES "_migrate_old_experiments"(id)
            )""")

        broken_before = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE '%_migrate_old_%'"))
        @test length(broken_before) == 3

        # Call the heal helper directly.
        HimalayaUI._fix_fk_references_after_autoincrement_migration!(db)

        broken_after = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE '%_migrate_old_%'"))
        @test isempty(broken_after)
        SQLite.close(db)
    end
end

@testset "_fix_fk_references_after_autoincrement_migration! is idempotent on healed DB" begin
    # Once the DB is healed, calling the heal helper again must not error or
    # re-rewrite. Detector LIKE clause must yield empty so the heal loop
    # short-circuits at the early return (review issue #23).
    mktempdir() do tmp
        db_path = joinpath(tmp, "corrupted-idem.db")
        db = SQLite.DB(db_path)
        DBInterface.execute(db,
            "CREATE TABLE samples (id INTEGER PRIMARY KEY AUTOINCREMENT, x TEXT)")
        DBInterface.execute(db, """
            CREATE TABLE refs_to_samples (
                id INTEGER PRIMARY KEY,
                sample_id INTEGER REFERENCES "_migrate_old_samples"(id)
            )""")

        # First call heals.
        HimalayaUI._fix_fk_references_after_autoincrement_migration!(db)
        broken1 = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE '%_migrate_old_%'"))
        @test isempty(broken1)

        # Second call must short-circuit without throwing.
        HimalayaUI._fix_fk_references_after_autoincrement_migration!(db)
        broken2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE '%_migrate_old_%'"))
        @test isempty(broken2)
        SQLite.close(db)
    end
end

@testset "user_actions.client_id column exists on fresh DB" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        cols = Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(user_actions)"))
        @test any(c -> c.name == "client_id", cols)
        SQLite.close(db)
    end
end

@testset "open_db adds client_id to legacy user_actions schema" begin
    mktempdir() do tmp
        path = joinpath(tmp, "test.db")
        # Build a legacy user_actions table without client_id. Mirror the
        # CURRENT create_schema! DDL ([db.jl:146-156]) verbatim minus the
        # new column — column names, defaults, and FKs must match exactly,
        # otherwise the test exercises a schema that never shipped.
        db = SQLite.DB(path)
        DBInterface.execute(db, """
            CREATE TABLE user_actions (
                id              INTEGER PRIMARY KEY,
                user_id         INTEGER REFERENCES users(id) ON DELETE SET NULL,
                timestamp       DATETIME DEFAULT CURRENT_TIMESTAMP,
                action          TEXT,
                entity_type     TEXT,
                entity_id       INTEGER,
                note            TEXT,
                payload         TEXT,
                undoes_event_id INTEGER REFERENCES user_actions(id)
            )
        """)
        SQLite.close(db)

        # Re-open via open_db: migrate_schema!'s ALTER TABLE adds the column
        db = HimalayaUI.open_db(path)
        cols = Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(user_actions)"))
        @test any(c -> c.name == "client_id", cols)
        SQLite.close(db)
    end
end

@testset "user_actions.client_op_id column exists on fresh DB" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        cols = Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(user_actions)"))
        @test any(c -> c.name == "client_op_id", cols)
        SQLite.close(db)
    end
end

@testset "idempotent_responses table exists on fresh DB" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        tables = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table' AND name='idempotent_responses'"))
        @test length(tables) == 1
        cols = Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(idempotent_responses)"))
        @test Set(c.name for c in cols) == Set(["client_op_id", "status_code", "body", "created_at"])
        @test any(c -> c.name == "client_op_id" && c.pk == 1, cols)
        SQLite.close(db)
    end
end

@testset "open_db adds client_op_id to legacy user_actions schema" begin
    mktempdir() do tmp
        path = joinpath(tmp, "test.db")
        db = SQLite.DB(path)
        DBInterface.execute(db, """
            CREATE TABLE user_actions (
                id              INTEGER PRIMARY KEY,
                user_id         INTEGER REFERENCES users(id) ON DELETE SET NULL,
                timestamp       DATETIME DEFAULT CURRENT_TIMESTAMP,
                action          TEXT,
                entity_type     TEXT,
                entity_id       INTEGER,
                note            TEXT,
                payload         TEXT,
                undoes_event_id INTEGER REFERENCES user_actions(id),
                client_id       TEXT
            )
        """)
        SQLite.close(db)
        db = HimalayaUI.open_db(path)
        cols = Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(user_actions)"))
        @test any(c -> c.name == "client_op_id", cols)
        SQLite.close(db)
    end
end

@testset "open_db creates idempotent_responses on legacy DB" begin
    mktempdir() do tmp
        path = joinpath(tmp, "test.db")
        db = SQLite.DB(path)
        DBInterface.execute(db, "CREATE TABLE foo (x INTEGER)")
        SQLite.close(db)
        db = HimalayaUI.open_db(path)
        tables = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table' AND name='idempotent_responses'"))
        @test length(tables) == 1
        SQLite.close(db)
    end
end

@testset "client_op_id partial index present" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        idx = Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='index' AND name='idx_events_by_client_op_id'"))
        @test length(idx) == 1
        @test occursin("WHERE client_op_id IS NOT NULL", String(idx[1].sql))
        SQLite.close(db)
    end
end

@testset "I2: partial unique index on user_actions(client_op_id, action, entity_id)" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        idx = Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='index' AND name='idx_events_unique_op'"))
        @test length(idx) == 1
        sql = String(idx[1].sql)
        @test occursin("UNIQUE", uppercase(sql))
        @test occursin("WHERE", sql)
        @test occursin("client_op_id IS NOT NULL", sql)
        SQLite.close(db)
    end
end

@testset "I2: partial unique index also installed on legacy DB via migrate_schema!" begin
    mktempdir() do tmp
        path = joinpath(tmp, "test.db")
        # Simulate a legacy DB (any DB that didn't already have this index).
        db = SQLite.DB(path)
        DBInterface.execute(db, "CREATE TABLE foo (x INTEGER)")
        SQLite.close(db)
        db = HimalayaUI.open_db(path)
        idx = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='index' AND name='idx_events_unique_op'"))
        @test length(idx) == 1
        SQLite.close(db)
    end
end

@testset "open_db: legacy user_actions populated with NULL client_op_id rows survives upgrade (issue #15)" begin
    # Reviewer-flagged gap: prior tests created legacy user_actions tables but
    # inserted ZERO rows before calling open_db, so they didn't exercise the
    # "rows exist before ALTER + index install" path. Pin the two properties
    # the partial unique index design depends on:
    #   (a) Pre-existing NULL-client_op_id rows survive ALTER ADD COLUMN +
    #       CREATE UNIQUE INDEX (the partial WHERE excludes them).
    #   (b) Two such rows with same (action, entity_id) coexist (partial
    #       WHERE keeps NULL pairs out of the unique constraint).
    #   (c) After upgrade, a non-NULL duplicate IS rejected by the live
    #       constraint.
    mktempdir() do tmp
        path = joinpath(tmp, "test.db")
        db = SQLite.DB(path)
        # Legacy schema mirroring create_schema! pre-Plan-8.
        DBInterface.execute(db, """
            CREATE TABLE users (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                username TEXT UNIQUE NOT NULL
            )""")
        DBInterface.execute(db, """
            CREATE TABLE user_actions (
                id              INTEGER PRIMARY KEY,
                user_id         INTEGER REFERENCES users(id) ON DELETE SET NULL,
                timestamp       DATETIME DEFAULT CURRENT_TIMESTAMP,
                action          TEXT,
                entity_type     TEXT,
                entity_id       INTEGER,
                note            TEXT
            )""")
        DBInterface.execute(db, "INSERT INTO users (username) VALUES ('alice')")
        # Seed two rows with same (action, entity_id), both NULL on the
        # not-yet-existing client_op_id column.
        DBInterface.execute(db,
            "INSERT INTO user_actions (user_id, action, entity_type, entity_id) VALUES (?, ?, ?, ?)",
            [1, "peak_added", "exposure", 42])
        DBInterface.execute(db,
            "INSERT INTO user_actions (user_id, action, entity_type, entity_id) VALUES (?, ?, ?, ?)",
            [1, "peak_added", "exposure", 42])
        SQLite.close(db)

        # Upgrade — open_db calls migrate_schema!. Must not throw despite the
        # partial-legacy schema (no exposures/samples/etc tables) AND must
        # preserve all rows.
        db = HimalayaUI.open_db(path)

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, action, client_op_id FROM user_actions ORDER BY id"))
        @test length(rows) == 2                            # (a) rows survive
        @test all(r -> ismissing(r.client_op_id), rows)
        @test rows[1].action == "peak_added"
        # Same (action, entity_id) appears twice with NULL op_id without
        # tripping the unique index — (b).
        same = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions WHERE action = ? AND entity_id = ?",
            ["peak_added", 42]))
        @test same[1].c == 2

        # The partial unique index installed.
        idx = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='index' AND name='idx_events_unique_op'"))
        @test length(idx) == 1

        # (c) A non-NULL op_id duplicate is rejected post-upgrade.
        DBInterface.execute(db, """
            INSERT INTO user_actions (action, entity_type, entity_id, client_op_id)
            VALUES ('peak_added', 'exposure', 42, 'op-fresh')""")
        @test_throws Exception DBInterface.execute(db, """
            INSERT INTO user_actions (action, entity_type, entity_id, client_op_id)
            VALUES ('peak_added', 'exposure', 42, 'op-fresh')""")

        SQLite.close(db)
    end
end

@testset "I2: duplicate (client_op_id, action, entity_id) rejected at DB level" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        # Seed FK targets.
        DBInterface.execute(db,
            "INSERT INTO experiments (name, path) VALUES ('e', '/p')")
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, name) VALUES (1, 'A1')")
        res = DBInterface.execute(db,
            "INSERT INTO exposures (sample_id, filename) VALUES (1, 'f')")
        exp_id = Int(DBInterface.lastrowid(res))

        DBInterface.execute(db, """
            INSERT INTO user_actions (action, entity_type, entity_id, client_op_id)
            VALUES ('peak_added', 'exposure', ?, 'op-x')""", [exp_id])
        @test_throws Exception DBInterface.execute(db, """
            INSERT INTO user_actions (action, entity_type, entity_id, client_op_id)
            VALUES ('peak_added', 'exposure', ?, 'op-x')""", [exp_id])

        # Multiple events under one op_id with different actions are still allowed.
        DBInterface.execute(db, """
            INSERT INTO user_actions (action, entity_type, entity_id, client_op_id)
            VALUES ('index_confirmed', 'exposure', ?, 'op-x')""", [exp_id])

        # NULL client_op_id rows are not constrained — partial WHERE clause.
        for _ in 1:3
            DBInterface.execute(db, """
                INSERT INTO user_actions (action, entity_type, entity_id, client_op_id)
                VALUES ('peak_added', 'exposure', ?, NULL)""", [exp_id])
        end

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM user_actions WHERE action = 'peak_added' AND entity_id = ?", [exp_id]))
        @test length(rows) == 4  # 1 with op-x + 3 with NULL
        SQLite.close(db)
    end
end

@testset "no INSERT in routes_*.jl uses empty-string placeholder pattern (#67 guardrail)" begin
    src_dir = joinpath(dirname(@__FILE__), "..", "src")
    routes = filter(f -> startswith(basename(f), "routes_") && endswith(f, ".jl"),
                    readdir(src_dir, join=true))
    @test !isempty(routes)
    for f in routes
        text = read(f, String)
        # Reject any `VALUES (...)` clause containing two or more `''`
        # literals — the empty-string-in-NOT-NULL placeholder pattern from
        # #67. Single `''` (e.g. as a deliberate sentinel for one column
        # later overwritten unconditionally) is permitted.
        m = match(r"VALUES\s*\([^)]*''[^)]*''[^)]*\)", text)
        @test m === nothing
    end
end

# Retired (ingestion redesign, Phase A): migrate_samples_naming! no longer suffixes
# duplicate sample names. Sample labels may legitimately repeat across loads; identity is
# (load_id, slot_index), not the label. The "duplicate labels survive naming+collapse"
# testset in test_ingestion_schema.jl covers the new contract. (Two dup-suffix testsets
# removed here: "duplicate names suffixed by ascending id" and "dup suffix avoids
# user-named -N collision".) See db.jl migrate_samples_name_collapse! and spec §10.

@testset "migrate_samples_naming! — idempotent on second run" begin
    mktempdir() do tmp
        db = SQLite.DB(joinpath(tmp, "h.db"))
        # Build canonical post-migration shape directly.
        DBInterface.execute(db, """CREATE TABLE samples (
            id INTEGER PRIMARY KEY AUTOINCREMENT, experiment_id INTEGER,
            name TEXT, display_name TEXT, notes TEXT)""")
        DBInterface.execute(db,
            "CREATE UNIQUE INDEX samples_unique_name ON samples(experiment_id, name)")
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (1, 1, 'JC001', 'DOPC')")
        # Migration should be a no-op (sentinel triggers).
        HimalayaUI.migrate_samples_naming!(db)
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name, display_name FROM samples"))
        @test length(rows) == 1
        @test rows[1].name == "JC001"
        @test rows[1].display_name == "DOPC"
    end
end

@testset "migrate_samples_naming! — legacy (label, name) → (name, display_name)" begin
    mktempdir() do tmp
        dbpath = joinpath(tmp, "h.db")
        db     = SQLite.DB(dbpath)
        # Synthetic LEGACY shape (pre-rename), bypassing open_db's migrations.
        DBInterface.execute(db, "PRAGMA foreign_keys = ON")
        DBInterface.execute(db, """CREATE TABLE experiments (
            id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, path TEXT,
            data_dir TEXT, analysis_dir TEXT, manifest_path TEXT, config TEXT,
            experiment_type TEXT, energy_kev REAL, flight_path_m REAL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP)""")
        DBInterface.execute(db, """CREATE TABLE samples (
            id INTEGER PRIMARY KEY AUTOINCREMENT, experiment_id INTEGER REFERENCES experiments(id),
            label TEXT, name TEXT, notes TEXT)""")
        DBInterface.execute(db,
            "INSERT INTO experiments (id, name, path) VALUES (?, ?, ?)", [1, "exp", "/tmp"])
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, label, name, notes) VALUES (?, ?, ?, ?)",
            [1, "JC001", "DOPC + chol", "first run"])
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, label, name, notes) VALUES (?, ?, ?, ?)",
            [1, "", "fallback only", nothing])
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, label, name, notes) VALUES (?, ?, ?, ?)",
            [1, "JC002", "second run", nothing])

        HimalayaUI.migrate_samples_naming!(db)

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name, display_name, notes FROM samples ORDER BY id"))
        @test rows[1].name == "JC001";          @test rows[1].display_name == "DOPC + chol"
        @test rows[2].name == "fallback only";  @test rows[2].display_name == "fallback only"
        @test rows[3].name == "JC002";          @test rows[3].display_name == "second run"

        # Old `label` column is gone.
        cols = [r.name for r in Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info('samples')"))]
        @test "label" ∉ cols
        @test "display_name" ∈ cols

        # Unique index exists.
        idxs = [r.name for r in Tables.rowtable(DBInterface.execute(db,
            "PRAGMA index_list('samples')"))]
        @test "samples_unique_name" ∈ idxs
    end
end

@testset "open_db: pre-Plan-7 legacy (label, name) DB collapses to the single friendly name" begin
    mktempdir() do tmp
        dbpath = joinpath(tmp, "h.db")
        # Build the worst-case fixture: legacy (label, name) shape WITHOUT AUTOINCREMENT.
        # This emulates a pre-Plan-7 deployment that never ran any migrations.
        db = SQLite.DB(dbpath)
        DBInterface.execute(db, "PRAGMA foreign_keys = ON")
        DBInterface.execute(db, """CREATE TABLE experiments (
            id INTEGER PRIMARY KEY, name TEXT, path TEXT,
            data_dir TEXT, analysis_dir TEXT, manifest_path TEXT, config TEXT,
            experiment_type TEXT, energy_kev REAL, flight_path_m REAL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP)""")
        DBInterface.execute(db, """CREATE TABLE samples (
            id INTEGER PRIMARY KEY,
            experiment_id INTEGER REFERENCES experiments(id),
            label TEXT, name TEXT, notes TEXT)""")
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'exp', '/tmp', '/tmp', '/tmp')")
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, label, name) VALUES (?, ?, ?, ?)",
            [1, 1, "JC001", "DOPC + chol"])
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, label, name) VALUES (?, ?, ?, ?)",
            [2, 1, "JC002", "POPC"])
        # Close and reopen via open_db to trigger the FULL migration chain.
        SQLite.close(db)

        db2 = HimalayaUI.open_db(dbpath)
        scols = lowercase.(String.(getproperty.(Tables.rowtable(
            DBInterface.execute(db2, "PRAGMA table_info(samples)")), :name)))
        @test "name" in scols && !("display_name" in scols)   # collapsed to one label
        rows = Tables.rowtable(DBInterface.execute(db2,
            "SELECT id, name FROM samples ORDER BY id"))
        @test length(rows) == 2
        # Post-collapse the single `name` is the friendly text (old `name` column, routed
        # name→display_name→name by the naming+collapse migrations). The old machine `label`
        # (JC001/JC002) is dropped — identity is now (load_id, slot_index).
        @test rows[1].name == "DOPC + chol"
        @test rows[2].name == "POPC"
    end
end

@testset "migrate_exposures_unique_filename!" begin
    @testset "clean DB: index exists, no warnings" begin
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
            try
                # open_db already runs the migration; assert the index is in place.
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT name FROM sqlite_master WHERE type='index' AND name='exposures_unique_filename'"))
                @test length(rows) == 1
            finally
                close(db)
            end
        end
    end

    @testset "idempotent re-run: no-op" begin
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
            try
                # Calling the helper directly a second time should be a no-op (no error,
                # no warnings, index still present).
                HimalayaUI.migrate_exposures_unique_filename!(db)
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT name FROM sqlite_master WHERE type='index' AND name='exposures_unique_filename'"))
                @test length(rows) == 1
            finally
                close(db)
            end
        end
    end

    @testset "synthetic duplicates: rename + warn + index created" begin
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
            try
                # Drop the index so we can synthesize duplicates and re-run the helper.
                DBInterface.execute(db, "DROP INDEX IF EXISTS exposures_unique_filename")

                # Seed a sample + two duplicate exposures via raw SQL (bypassing the upsert).
                eid = let res = DBInterface.execute(db,
                          "INSERT INTO experiments (name, path) VALUES ('e', '/p')")
                    Int(DBInterface.lastrowid(res))
                end
                sid = let res = DBInterface.execute(db,
                          "INSERT INTO samples (experiment_id, name) VALUES (?, 'S1')", [eid])
                    Int(DBInterface.lastrowid(res))
                end
                # Two rows with the same (sample_id, filename).
                DBInterface.execute(db,
                    "INSERT INTO exposures (sample_id, filename, kind) VALUES (?, 'JC001-007', 'simple')", [sid])
                DBInterface.execute(db,
                    "INSERT INTO exposures (sample_id, filename, kind) VALUES (?, 'JC001-007', 'simple')", [sid])

                # Run the helper directly (FK-heal pattern).
                @test_logs (:warn, r"Renamed duplicate exposure"i) min_level=Logging.Warn match_mode=:any begin
                    HimalayaUI.migrate_exposures_unique_filename!(db)
                end

                # Oldest id keeps the bare name; second is suffixed.
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, filename FROM exposures WHERE sample_id = ? ORDER BY id ASC", [sid]))
                @test length(rows) == 2
                @test rows[1].filename == "JC001-007"
                @test rows[2].filename == "JC001-007-2"

                # Index now present.
                idx = Tables.rowtable(DBInterface.execute(db,
                    "SELECT name FROM sqlite_master WHERE type='index' AND name='exposures_unique_filename'"))
                @test length(idx) == 1
            finally
                close(db)
            end
        end
    end
end

@testset "open_db chmods owned DB file to 0664" begin
    # SQLite creates files at 0644 (hardcoded in os_unix.c, ignores umask).
    # open_db chmods to 0664 so other curators in the shared group can
    # write the same DB on multi-user deploys. Regression for #127: the
    # chmod path must only skip the cross-user case (another user owns
    # the file). On a file WE own — which is what this test exercises —
    # chmod must succeed; any failure (read-only mount, immutable bit,
    # etc.) must propagate, not silently leave 0644.
    #
    # WAL sidecars (-wal / -shm) are not asserted: -wal is checkpointed
    # away on close, and -shm doesn't exist at the chmod call site on a
    # fresh DB (chmod runs immediately after WAL PRAGMA, before any
    # write). The sidecar branch shares the same code path as the main
    # file — covered by inspection.
    Sys.isunix() || return
    mktempdir() do tmp
        path = joinpath(tmp, "h.db")
        db = HimalayaUI.open_db(path)
        close(db)
        @test (stat(path).mode & 0o777) == 0o664
    end
end

@testset "open_db drops the retired series, comparison and sample chat tables, keeping their log rows" begin
    mktempdir() do tmp
        path = joinpath(tmp, "legacy.db")
        db = HimalayaUI.open_db(path)
        exp_id = create_experiment!(db; name = "E", path = tmp,
            data_dir = joinpath(tmp, "data"), analysis_dir = joinpath(tmp, "analysis"))
        s_id = create_sample!(db; experiment_id = exp_id, name = "S")
        e_id = create_exposure!(db; experiment_id = exp_id, sample_id = s_id, filename = "F")
        u_id = Int(DBInterface.lastrowid(DBInterface.execute(db,
            "INSERT INTO users (username) VALUES ('alice')")))
        for stmt in (
            """CREATE TABLE comparisons (id INTEGER PRIMARY KEY AUTOINCREMENT, title TEXT,
                   created_by INTEGER REFERENCES users(id) ON DELETE SET NULL,
                   forked_from_id INTEGER REFERENCES comparisons(id) ON DELETE SET NULL)""",
            """CREATE TABLE comparison_members (id INTEGER PRIMARY KEY,
                   comparison_id INTEGER NOT NULL REFERENCES comparisons(id) ON DELETE CASCADE,
                   exposure_id INTEGER REFERENCES exposures(id) ON DELETE SET NULL)""",
            """CREATE TABLE comparison_messages (id INTEGER PRIMARY KEY,
                   comparison_id INTEGER NOT NULL REFERENCES comparisons(id) ON DELETE CASCADE,
                   body TEXT NOT NULL)""",
            """CREATE TABLE comparison_pins (
                   user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
                   comparison_id INTEGER NOT NULL REFERENCES comparisons(id) ON DELETE CASCADE,
                   PRIMARY KEY (user_id, comparison_id))""",
            """CREATE TABLE series (id INTEGER PRIMARY KEY AUTOINCREMENT, title TEXT,
                   forked_from_id INTEGER REFERENCES series(id) ON DELETE SET NULL)""",
            """CREATE TABLE series_members (id INTEGER PRIMARY KEY,
                   series_id INTEGER NOT NULL REFERENCES series(id) ON DELETE CASCADE,
                   exposure_id INTEGER REFERENCES exposures(id) ON DELETE SET NULL)""",
            """CREATE TABLE series_samples (id INTEGER PRIMARY KEY,
                   series_id INTEGER NOT NULL REFERENCES series(id) ON DELETE CASCADE,
                   sample_id INTEGER NOT NULL REFERENCES samples(id) ON DELETE CASCADE,
                   position INTEGER NOT NULL)""",
            """CREATE TABLE series_messages (id INTEGER PRIMARY KEY,
                   series_id INTEGER NOT NULL REFERENCES series(id) ON DELETE CASCADE,
                   body TEXT NOT NULL)""",
            """CREATE TABLE series_pins (
                   user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
                   series_id INTEGER NOT NULL REFERENCES series(id) ON DELETE CASCADE,
                   PRIMARY KEY (user_id, series_id))""",
            """CREATE TABLE sample_messages (id INTEGER PRIMARY KEY,
                   sample_id INTEGER REFERENCES samples(id),
                   author_id INTEGER REFERENCES users(id) ON DELETE SET NULL,
                   body TEXT NOT NULL)""",
        )
            DBInterface.execute(db, stmt)
        end
        DBInterface.execute(db, "INSERT INTO comparisons (id, title, created_by) VALUES (1, 'c', ?)", [u_id])
        DBInterface.execute(db, "INSERT INTO comparisons (id, title, forked_from_id) VALUES (2, 'fork', 1)")
        DBInterface.execute(db, "INSERT INTO comparison_members (comparison_id, exposure_id) VALUES (1, ?)", [e_id])
        DBInterface.execute(db, "INSERT INTO comparison_messages (comparison_id, body) VALUES (1, 'hi')")
        DBInterface.execute(db, "INSERT INTO comparison_pins (user_id, comparison_id) VALUES (?, 1)", [u_id])
        DBInterface.execute(db, "INSERT INTO series (id, title) VALUES (1, 's')")
        DBInterface.execute(db, "INSERT INTO series (id, title, forked_from_id) VALUES (2, 'fork', 1)")
        DBInterface.execute(db, "INSERT INTO series_members (series_id, exposure_id) VALUES (1, ?)", [e_id])
        DBInterface.execute(db, "INSERT INTO series_samples (series_id, sample_id, position) VALUES (1, ?, 0)", [s_id])
        DBInterface.execute(db, "INSERT INTO series_messages (series_id, body) VALUES (1, 'hi')")
        DBInterface.execute(db, "INSERT INTO series_pins (user_id, series_id) VALUES (?, 1)", [u_id])
        DBInterface.execute(db, "INSERT INTO sample_messages (sample_id, author_id, body) VALUES (?, ?, 'hi')", [s_id, u_id])
        DBInterface.execute(db,
            """INSERT INTO user_actions (user_id, action, entity_type, entity_id, payload)
               VALUES (?, 'series_created', 'series', 1, '{}')""", [u_id])
        close(db)

        db = HimalayaUI.open_db(path)
        tables = Set(String(r.name) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type = 'table'")))
        for t in HimalayaUI.RETIRED_TABLES
            @test t ∉ tables
        end
        @test isempty(Tables.rowtable(DBInterface.execute(db, "PRAGMA foreign_key_check")))
        @test length(get_exposures(db, s_id)) == 1
        @test length(get_samples(db, exp_id)) == 1
        @test only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM user_actions WHERE action = 'series_created'"))).n == 1
        close(db)

        db = HimalayaUI.open_db(path)
        @test "series" ∉ Set(String(r.name) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type = 'table'")))
        close(db)
    end
end

@testset "open_db drops the legacy tag tables, and their old log rows still replay" begin
    mktempdir() do tmp
        path = joinpath(tmp, "legacy.db")
        db = HimalayaUI.open_db(path)
        exp_id = create_experiment!(db; name = "E", path = tmp,
            data_dir = joinpath(tmp, "data"), analysis_dir = joinpath(tmp, "analysis"))
        s_id = create_sample!(db; experiment_id = exp_id, name = "S")
        e_id = create_exposure!(db; experiment_id = exp_id, sample_id = s_id, filename = "F")
        DBInterface.execute(db, """CREATE TABLE IF NOT EXISTS sample_tags (id INTEGER PRIMARY KEY,
            sample_id INTEGER REFERENCES samples(id), key TEXT NOT NULL, value TEXT NOT NULL,
            source TEXT DEFAULT 'manual')""")
        DBInterface.execute(db, """CREATE TABLE IF NOT EXISTS exposure_tags (id INTEGER PRIMARY KEY,
            exposure_id INTEGER REFERENCES exposures(id), key TEXT NOT NULL, value TEXT NOT NULL,
            source TEXT DEFAULT 'manual')""")
        DBInterface.execute(db, "INSERT INTO sample_tags (sample_id, key, value) VALUES (?, 'k', 'v')", [s_id])
        DBInterface.execute(db,
            "INSERT INTO exposure_tags (exposure_id, key, value) VALUES (?, 'rejection_reason', 'r')", [e_id])
        for (kind, ent, eid) in (("add_tag", "sample", s_id), ("add_tag", "exposure", e_id),
                                 ("edit_tag", "exposure", e_id), ("remove_tag", "exposure", e_id))
            DBInterface.execute(db,
                """INSERT INTO user_actions (user_id, action, entity_type, entity_id, payload)
                   VALUES (NULL, ?, ?, ?, '{"key":"rejection_reason","value":"r"}')""",
                [kind, ent, eid])
        end
        close(db)

        db = HimalayaUI.open_db(path)
        tables = Set(String(r.name) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type = 'table'")))
        @test "sample_tags" ∉ tables
        @test "exposure_tags" ∉ tables
        @test isempty(Tables.rowtable(DBInterface.execute(db, "PRAGMA foreign_key_check")))
        @test only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM user_actions WHERE action LIKE '%_tag'"))).n == 4
        @test_nowarn HimalayaUI.rebuild_views_from_log!(db, e_id)
        @test_nowarn HimalayaUI.rebuild_views_from_log!(db, s_id; entity_type = "sample")
        @test length(get_exposures(db, s_id)) == 1
        close(db)
    end
end

@testset "a fresh DB has schema_migrations and none of the retired tables" begin
    db = SQLite.DB()
    create_schema!(db)
    migrate_schema!(db)
    tables = Set(String(r.name) for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT name FROM sqlite_master WHERE type = 'table'")))
    @test "schema_migrations" in tables
    for t in HimalayaUI.RETIRED_TABLES
        @test t ∉ tables
    end
end

@testset "sample_subtractions can hold an explicit none, and an older database gains that without losing rows" begin
    q(db, sql) = Tables.rowtable(DBInterface.execute(db, sql))
    info(db) = [(r.name, r.type, r.notnull, r.pk) for r in q(db, "PRAGMA table_info(sample_subtractions)")]
    fresh = HimalayaUI.open_db(joinpath(mktempdir(), "fresh.db"))
    path = joinpath(mktempdir(), "old.db")
    db = HimalayaUI.open_db(path)
    exp_id = create_experiment!(db; path = "/d", data_dir = "/d", analysis_dir = "/a")
    a = create_sample!(db; experiment_id = exp_id, name = "lipidA 01A (S0101)")
    b = create_sample!(db; experiment_id = exp_id, name = "lipidA 01B (B0100)")
    DBInterface.execute(db, "DROP TABLE sample_subtractions")
    DBInterface.execute(db, """
        CREATE TABLE sample_subtractions (
            sample_id     INTEGER PRIMARY KEY REFERENCES samples(id) ON DELETE CASCADE,
            ref_sample_id INTEGER NOT NULL REFERENCES samples(id) ON DELETE CASCADE
        )""")
    DBInterface.execute(db, "INSERT INTO sample_subtractions (sample_id, ref_sample_id) VALUES (?, ?)", [a, b])
    DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_SUBTRACTION_NONE])
    SQLite.close(db)
    migrated = HimalayaUI.open_db(path)
    @test info(migrated) == info(fresh)
    @test [(Int(r.sample_id), Int(r.ref_sample_id)) for r in q(migrated, "SELECT * FROM sample_subtractions")] == [(a, b)]
    DBInterface.execute(migrated, "INSERT INTO sample_subtractions (sample_id, ref_sample_id) VALUES (?, NULL)", [b])
    @test ismissing(only(q(migrated, "SELECT ref_sample_id FROM sample_subtractions WHERE sample_id = $b")).ref_sample_id)
    DBInterface.execute(migrated, "DELETE FROM samples WHERE id = ?", [a])
    @test [Int(r.sample_id) for r in q(migrated, "SELECT * FROM sample_subtractions")] == [b]
end
