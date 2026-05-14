using Test, SQLite, DBInterface, Tables, Logging
using HimalayaUI: create_schema!, migrate_schema!, create_experiment!, create_sample!,
                  create_exposure!, get_experiment, get_samples, get_exposures,
                  migrate_r2_split_peaks!

@testset "db schema" begin
    db = SQLite.DB()  # in-memory
    create_schema!(db)

    tables = Set(r[1] for r in DBInterface.execute(db,
        "SELECT name FROM sqlite_master WHERE type='table'"))

    for t in ["users", "experiments", "samples", "sample_tags",
              "exposures", "exposure_sources", "exposure_tags",
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
    samp_id = create_sample!(db; experiment_id=exp_id)
    DBInterface.execute(db,
        "INSERT INTO exposures (sample_id, filename, kind, selected, status, image_path)
         VALUES (?, 'test.dat', 'file', 0, 'accepted', '/tmp/test.tiff')", [samp_id])
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

    s_id = create_sample!(db; experiment_id = exp_id, name = "D1", display_name = "UX1")
    @test s_id == 1

    samples = get_samples(db, exp_id)
    @test length(samples) == 1
    @test first(samples).name == "D1"

    e_id = create_exposure!(db; sample_id = s_id, filename = "JC001", kind = "file")
    @test e_id == 1

    exposures = get_exposures(db, s_id)
    @test length(exposures) == 1
    @test first(exposures).filename == "JC001"
end

@testset "experiments table has config columns" begin
    db = SQLite.DB()
    create_schema!(db)
    cols = [r.name for r in Tables.rowtable(DBInterface.execute(db,
        "PRAGMA table_info(experiments)"))]
    @test "config"          in cols
    @test "experiment_type" in cols
    @test "energy_kev"      in cols
    @test "flight_path_m"   in cols
end

@testset "create_experiment! stores config columns" begin
    db = SQLite.DB()
    create_schema!(db)
    exp_id = create_experiment!(db;
        name = "ConfigTest",
        path = "/tmp/x",
        data_dir = "/tmp/x/data",
        analysis_dir = "/tmp/x/analysis",
        config = "[experiment]\nname=\"X\"\n",
        experiment_type = "simple",
        energy_kev = 12.0,
        flight_path_m = 2.5)

    row = Tables.rowtable(DBInterface.execute(db,
        "SELECT config, experiment_type, energy_kev, flight_path_m FROM experiments WHERE id=?",
        [exp_id]))[1]

    @test contains(row.config, "[experiment]")
    @test row.experiment_type == "simple"
    @test row.energy_kev == 12.0
    @test row.flight_path_m == 2.5
end

@testset "migrate_schema! adds config columns to legacy DB" begin
    db = SQLite.DB()
    # Simulate legacy schema without new columns
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
    migrate_schema!(db)
    cols = [r.name for r in Tables.rowtable(DBInterface.execute(db,
        "PRAGMA table_info(experiments)"))]
    @test "config" in cols
    @test "experiment_type" in cols
    @test "energy_kev" in cols
    @test "flight_path_m" in cols
end

@testset "index_groups partial unique constraint on custom" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

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
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
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
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        # Open again — should not error or duplicate rows.
        db2 = HimalayaUI.open_db(joinpath(dir, "h.db"))
        @test true  # if no exception thrown, we're good
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
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
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
        s_id   = HimalayaUI.create_sample!(legacy; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(legacy; sample_id=s_id)
        DBInterface.execute(legacy,
            "INSERT INTO index_groups (exposure_id, kind) VALUES (?, 'custom')", [e_id])
        DBInterface.execute(legacy,
            "INSERT INTO index_groups (exposure_id, kind) VALUES (?, 'custom')", [e_id])
        SQLite.DBInterface.close!(legacy)

        @test_throws ErrorException HimalayaUI.open_db(db_path)
    end
end

@testset "open_db heals corrupted _migrate_old_* FK references on populated tables" begin
    # Real prod DB at /tmp/himalaya-dev.db exhibited 6 tables with FKs that
    # referenced the dead `_migrate_old_<entity>` temp-table name from a
    # prior `migrate_pk_to_autoincrement!` run. The previous fix function
    # only handled `auto_peaks`/`peak_curations`; the populated tables
    # (`sample_messages`, `sample_tags`, `exposure_tags`, `exposure_sources`,
    # `index_groups`, `index_group_members`) silently kept the dead refs
    # and any subsequent INSERT on a table with such an FK trips
    # SQLiteException("no such table: main._migrate_old_*"). Caught only by
    # running the actual server (smoke checklist), not by tests using
    # freshly-built DBs.
    #
    # This regression test simulates the post-corruption state directly:
    # builds a DB whose `sample_messages` references `_migrate_old_samples`,
    # then runs `open_db` and asserts (a) the FK is rewritten back to
    # `samples`, (b) an INSERT succeeds without 500.
    mktempdir() do tmp
        db_path = joinpath(tmp, "corrupted.db")
        db = SQLite.DB(db_path)
        # Build a minimal corrupted state: samples + sample_messages with
        # a dead FK reference. Mirrors what migrate_pk_to_autoincrement!
        # leaves behind on a DB whose ALTER TABLE RENAME tracking corrupted
        # references inside the rename transaction.
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
            CREATE TABLE sample_messages (
                id         INTEGER PRIMARY KEY,
                sample_id  INTEGER REFERENCES "_migrate_old_samples"(id),
                author_id  INTEGER REFERENCES users(id) ON DELETE SET NULL,
                body       TEXT NOT NULL,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP
            )""")
        DBInterface.execute(db, "INSERT INTO samples (label) VALUES ('S1')")
        SQLite.close(db)

        # Pre-condition: corrupted reference present.
        verify = SQLite.DB(db_path)
        broken_before = Tables.rowtable(DBInterface.execute(verify,
            "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE '%_migrate_old_%'"))
        @test length(broken_before) == 1
        SQLite.close(verify)

        # Heal via open_db (calls migrate_schema! which calls the fix helper).
        healed = HimalayaUI.open_db(db_path)
        broken_after = Tables.rowtable(DBInterface.execute(healed,
            "SELECT name FROM sqlite_master WHERE type='table' AND sql LIKE '%_migrate_old_%'"))
        @test isempty(broken_after)

        # The healed DB must accept an INSERT into sample_messages — the
        # operation that would 500 on the corrupted DB.
        u_id = HimalayaUI.get_or_create_user!(healed, "alice")
        DBInterface.execute(healed,
            "INSERT INTO sample_messages (sample_id, author_id, body) VALUES (?, ?, ?)",
            [1, u_id, "smoke-regression"])
        @test first(Tables.rowtable(DBInterface.execute(healed,
            "SELECT COUNT(*) AS c FROM sample_messages"))).c == 1
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
            "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e', '/p', '/d', '/a')")
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

# ─────────────────────────────────────────────────────────────────────────────
# Compare page (Plan §Phase 1, Task 1.1): comparisons / comparison_members /
# comparison_messages schema. See docs/superpowers/specs/2026-05-02-compare-page-design.md.
# ─────────────────────────────────────────────────────────────────────────────

@testset "Compare schema: tables and indexes exist after open_db" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        tables = Set(String[String(r.name) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table'"))])
        @test "comparisons" in tables
        @test "comparison_members" in tables
        @test "comparison_messages" in tables

        indexes = Set(String[String(r.name) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='index'"))])
        @test "idx_comparisons_forked_from" in indexes
        @test "idx_comparison_members_by_comparison" in indexes
        @test "idx_comparison_messages_comparison" in indexes
    end
end

@testset "Compare schema: comparisons has expected columns" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        cols = Set(String[String(c.name) for c in Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(comparisons)"))])
        for c in ("id", "title", "description", "content_hash",
                  "created_by", "created_at", "updated_at",
                  "forked_from_id", "forked_at_hash")
            @test c in cols
        end
        # comparisons.id must use AUTOINCREMENT (mention-target rule).
        sql = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='table' AND name='comparisons'"))).sql
        @test occursin("AUTOINCREMENT", String(sql))
    end
end

@testset "Compare schema: comparison_members columns and JSON CHECKs" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        cols = Set(String[String(c.name) for c in Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(comparison_members)"))])
        for c in ("id", "comparison_id", "exposure_id", "display_order",
                  "band_height", "y_offset", "normalization", "color_override",
                  "label_override", "q_window_min", "q_window_max",
                  "peak_display", "snapshot", "created_by", "created_at")
            @test c in cols
        end
        # Plain INTEGER PRIMARY KEY (not AUTOINCREMENT — members are not @-mentioned).
        sql = String(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='table' AND name='comparison_members'"))).sql)
        @test !occursin("AUTOINCREMENT", sql)
        # snapshot CHECK (json_valid(...))
        @test occursin("json_valid(snapshot)", sql)
        @test occursin("json_valid(peak_display)", sql)
    end
end

@testset "Compare schema: comparison_messages plain PK" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        sql = String(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='table' AND name='comparison_messages'"))).sql)
        @test !occursin("AUTOINCREMENT", sql)
        cols = Set(String[String(c.name) for c in Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(comparison_messages)"))])
        for c in ("id", "comparison_id", "author_id", "body", "created_at")
            @test c in cols
        end
    end
end

@testset "Compare schema: idempotent — second open_db is a no-op" begin
    mktempdir() do tmp
        path = joinpath(tmp, "h.db")
        db1 = HimalayaUI.open_db(path)
        SQLite.close(db1)
        # Second open must not error and must not duplicate the schema.
        db2 = HimalayaUI.open_db(path)
        tables = [String(r.name) for r in Tables.rowtable(DBInterface.execute(db2,
            "SELECT name FROM sqlite_master WHERE type='table' AND name='comparisons'"))]
        @test length(tables) == 1
        SQLite.close(db2)
    end
end

@testset "Compare schema: comparisons.id AUTOINCREMENT does not reuse freed ids" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        DBInterface.execute(db, """
            INSERT INTO comparisons (title, content_hash, created_at, updated_at)
            VALUES ('a', 'h1', '2026-01-01', '2026-01-01')""")
        DBInterface.execute(db, """
            INSERT INTO comparisons (title, content_hash, created_at, updated_at)
            VALUES ('b', 'h2', '2026-01-01', '2026-01-01')""")
        first_max = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT MAX(id) AS m FROM comparisons"))).m
        DBInterface.execute(db, "DELETE FROM comparisons")
        DBInterface.execute(db, """
            INSERT INTO comparisons (title, content_hash, created_at, updated_at)
            VALUES ('c', 'h3', '2026-01-01', '2026-01-01')""")
        new_id = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparisons"))).id
        @test Int(new_id) > Int(first_max)  # AUTOINCREMENT didn't reuse
    end
end

@testset "Compare schema: FK on comparison_members.comparison_id rejects orphans" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        @test_throws SQLite.SQLiteException DBInterface.execute(db, """
            INSERT INTO comparison_members
              (comparison_id, display_order, snapshot, created_at)
            VALUES (9999, 0, '{}', '2026-01-01')""")
    end
end

@testset "Compare schema: ON DELETE SET NULL on comparison_members.exposure_id" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db, """
            INSERT INTO comparisons (title, content_hash, created_at, updated_at)
            VALUES ('t', 'h', '2026-01-01', '2026-01-01')""")
        cmp_id = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparisons"))).id)
        DBInterface.execute(db, """
            INSERT INTO comparison_members
              (comparison_id, exposure_id, display_order, snapshot, created_at)
            VALUES (?, ?, 0, '{}', '2026-01-01')""", [cmp_id, e_id])
        DBInterface.execute(db, "DELETE FROM exposures WHERE id = ?", [e_id])
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id FROM comparison_members WHERE comparison_id = ?", [cmp_id]))
        @test length(rows) == 1
        @test ismissing(rows[1].exposure_id)
    end
end

@testset "Compare schema: ON DELETE CASCADE on comparison_members.comparison_id" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        DBInterface.execute(db, """
            INSERT INTO comparisons (title, content_hash, created_at, updated_at)
            VALUES ('t', 'h', '2026-01-01', '2026-01-01')""")
        cmp_id = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparisons"))).id)
        DBInterface.execute(db, """
            INSERT INTO comparison_members
              (comparison_id, display_order, snapshot, created_at)
            VALUES (?, 0, '{}', '2026-01-01')""", [cmp_id])
        DBInterface.execute(db, """
            INSERT INTO comparison_messages (comparison_id, body)
            VALUES (?, 'hi')""", [cmp_id])
        DBInterface.execute(db, "DELETE FROM comparisons WHERE id = ?", [cmp_id])
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparison_members WHERE comparison_id = ?", [cmp_id])))
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparison_messages WHERE comparison_id = ?", [cmp_id])))
    end
end

@testset "Compare schema: ON DELETE SET NULL on forked_from_id" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        DBInterface.execute(db, """
            INSERT INTO comparisons (title, content_hash, created_at, updated_at)
            VALUES ('parent', 'hp', '2026-01-01', '2026-01-01')""")
        parent_id = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparisons"))).id)
        DBInterface.execute(db, """
            INSERT INTO comparisons
              (title, content_hash, created_at, updated_at, forked_from_id, forked_at_hash)
            VALUES ('fork', 'hf', '2026-01-01', '2026-01-01', ?, ?)""",
            [parent_id, "hp"])
        fork_id = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparisons WHERE forked_from_id = ?", [parent_id]))).id)
        DBInterface.execute(db, "DELETE FROM comparisons WHERE id = ?", [parent_id])
        # Fork survives, forked_from_id is now NULL, forked_at_hash preserved.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT forked_from_id, forked_at_hash FROM comparisons WHERE id = ?", [fork_id]))
        @test length(rows) == 1
        @test ismissing(rows[1].forked_from_id)
        @test String(rows[1].forked_at_hash) == "hp"
    end
end

# Issue #67 — eliminate empty-string-in-NOT-NULL placeholder pattern.
# After PR #66 patched the symptom (NULLIF in the dispatcher COALESCE),
# this set of tests pins the structural fix: the four columns the route's
# mint-the-id INSERT used to seed with '' (`title`, `content_hash`,
# `created_at`, `updated_at`) are now nullable, so the route can use a
# clean NULL placeholder (`INSERT INTO comparisons DEFAULT VALUES`) and
# the dispatcher's plain `COALESCE(col, ?)` Just Works.

# Build a comparisons schema with the OLD strict NOT NULL columns so the
# `migrate_compare_relax_nullability!` heal path can be exercised. Mirrors
# the migrate_compare! shape pre-#67 (the four mint-the-id columns are
# `TEXT NOT NULL` rather than nullable).
function _build_legacy_strict_compare_db(path::String)
    db = SQLite.DB(path)
    DBInterface.execute(db, "PRAGMA foreign_keys = ON")
    DBInterface.execute(db, """
        CREATE TABLE users (id INTEGER PRIMARY KEY, username TEXT)
    """)
    DBInterface.execute(db, """
        CREATE TABLE comparisons (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            title           TEXT NOT NULL,
            description     TEXT,
            content_hash    TEXT NOT NULL,
            created_by      INTEGER REFERENCES users(id) ON DELETE SET NULL,
            created_at      TEXT NOT NULL,
            updated_at      TEXT NOT NULL,
            forked_from_id  INTEGER REFERENCES comparisons(id) ON DELETE SET NULL,
            forked_at_hash  TEXT
        )
    """)
    DBInterface.execute(db, """
        CREATE TABLE comparison_members (
            id              INTEGER PRIMARY KEY,
            comparison_id   INTEGER NOT NULL REFERENCES comparisons(id) ON DELETE CASCADE,
            exposure_id     INTEGER,
            display_order   INTEGER NOT NULL,
            snapshot        TEXT    NOT NULL CHECK (json_valid(snapshot)),
            created_at      TEXT    NOT NULL
        )
    """)
    DBInterface.execute(db, """
        CREATE TABLE comparison_pins (
            user_id        INTEGER NOT NULL REFERENCES users(id)       ON DELETE CASCADE,
            comparison_id  INTEGER NOT NULL REFERENCES comparisons(id) ON DELETE CASCADE,
            pinned_at      TEXT    NOT NULL DEFAULT CURRENT_TIMESTAMP,
            PRIMARY KEY (user_id, comparison_id)
        )
    """)
    db
end

@testset "Compare schema: comparisons NOT NULL columns relaxed (#67)" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        sql = String(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='table' AND name='comparisons'"))).sql)
        # All four mint-the-id placeholder columns must NOT carry NOT NULL on
        # a fresh DB. The route's `INSERT … DEFAULT VALUES` depends on this.
        for col in ("title", "content_hash", "created_at", "updated_at")
            @test !occursin(Regex("\\b$col\\s+TEXT\\s+NOT\\s+NULL\\b", "i"), sql)
        end
    end
end

@testset "Compare schema: relaxed schema accepts INSERT DEFAULT VALUES (#67)" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        # Proves the route can use a clean NULL placeholder. If this fires
        # `NOT NULL constraint failed`, the relax migration didn't run.
        DBInterface.execute(db, "INSERT INTO comparisons DEFAULT VALUES")
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, title, content_hash, created_at, updated_at FROM comparisons"))
        @test length(rows) == 1
        @test rows[1].id isa Integer
        @test ismissing(rows[1].title)
        @test ismissing(rows[1].content_hash)
        @test ismissing(rows[1].created_at)
        @test ismissing(rows[1].updated_at)
    end
end

@testset "migrate_compare_relax_nullability!: legacy NOT NULL schema rebuilds (#67)" begin
    mktempdir() do tmp
        path = joinpath(tmp, "legacy.db")
        db = _build_legacy_strict_compare_db(path)
        # Seed: a row with the #54-bug `created_at = ''` plus a normal row.
        DBInterface.execute(db, """
            INSERT INTO comparisons (title, content_hash, created_at, updated_at)
            VALUES ('broken', 'h1', '', '')""")
        DBInterface.execute(db, """
            INSERT INTO comparisons (title, content_hash, created_at, updated_at)
            VALUES ('normal', 'h2', '2026-05-01T00:00:00', '2026-05-01T00:00:00')""")
        broken_id = 1
        normal_id = 2
        # FK-related rows in BOTH comparison_members AND comparison_pins so
        # the heal walks both referrers (review suggestion #2 on PR #72 —
        # second referrer matters because _fix_fk_references_after_autoincrement_migration!
        # had a prod incident where only the first referrer was rebuilt).
        DBInterface.execute(db, """
            INSERT INTO comparison_members
              (comparison_id, display_order, snapshot, created_at)
            VALUES (?, 0, '{}', '2026-05-01T00:00:00')""",
            [broken_id])
        DBInterface.execute(db, """
            INSERT INTO users (id, username) VALUES (1, 'alice')""")
        DBInterface.execute(db, """
            INSERT INTO comparison_pins (user_id, comparison_id, pinned_at)
            VALUES (1, ?, '2026-05-01T00:00:00')""", [broken_id])

        HimalayaUI.migrate_compare_relax_nullability!(db)

        # Schema is relaxed
        sql = String(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='table' AND name='comparisons'"))).sql)
        for col in ("title", "content_hash", "created_at", "updated_at")
            @test !occursin(Regex("\\b$col\\s+TEXT\\s+NOT\\s+NULL\\b", "i"), sql)
        end

        # Existing rows survive
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, title, created_at FROM comparisons ORDER BY id"))
        @test length(rows) == 2
        @test String(rows[1].title) == "broken"
        @test String(rows[2].title) == "normal"
        @test String(rows[2].created_at) == "2026-05-01T00:00:00"

        # Issue body's edge case: the stale '' created_at gets healed.
        # Format MUST match the fresh-row format (`yyyy-mm-ddTHH:MM:SS.sssZ`
        # from `comparison_now_iso()`) so ComparisonSidebar's lexical sort
        # on `updated_at` keeps healed and fresh rows in monotonic order
        # (review suggestion #1 on PR #72: SQLite's CURRENT_TIMESTAMP uses
        # `YYYY-MM-DD HH:MM:SS` and space < T, so it would lexically precede
        # fresh rows of the same instant).
        @test rows[1].created_at != ""
        @test occursin(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}", String(rows[1].created_at))
        @test endswith(String(rows[1].created_at), "Z")

        # FK from comparison_members still resolves (heal worked).
        DBInterface.execute(db, """
            INSERT INTO comparison_members
              (comparison_id, display_order, snapshot, created_at)
            VALUES (?, 1, '{}', '2026-05-01T00:00:00')""", [broken_id])
        # FK from comparison_pins (the SECOND referrer) also resolves —
        # this is the suggestion #2 coverage: prove the heal walked both
        # referrers, not just the first one.
        DBInterface.execute(db, """
            INSERT INTO comparison_pins (user_id, comparison_id, pinned_at)
            VALUES (1, ?, '2026-05-02T00:00:00')""", [normal_id])
        # And the orphan-rejection path: inserting a member pointing at a
        # nonexistent comparison must still 19=>SQLITE_CONSTRAINT.
        @test_throws SQLite.SQLiteException DBInterface.execute(db, """
            INSERT INTO comparison_members
              (comparison_id, display_order, snapshot, created_at)
            VALUES (9999, 0, '{}', '2026-05-01T00:00:00')""")
        # Same orphan rejection on comparison_pins (proves the FK
        # constraint, not just the column, survived the heal).
        @test_throws SQLite.SQLiteException DBInterface.execute(db, """
            INSERT INTO comparison_pins (user_id, comparison_id, pinned_at)
            VALUES (1, 9999, '2026-05-01T00:00:00')""")

        # ON DELETE CASCADE still fires for BOTH referrers: deleting a
        # comparison drops its members AND its pins (proves the FK action
        # survived the rebuild for both edges).
        m_before = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM comparison_members WHERE comparison_id = ?",
            [broken_id]))).n)
        p_before = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM comparison_pins WHERE comparison_id = ?",
            [broken_id]))).n)
        @test m_before > 0
        @test p_before > 0
        DBInterface.execute(db, "DELETE FROM comparisons WHERE id = ?", [broken_id])
        m_after = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM comparison_members WHERE comparison_id = ?",
            [broken_id]))).n)
        p_after = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM comparison_pins WHERE comparison_id = ?",
            [broken_id]))).n)
        @test m_after == 0
        @test p_after == 0

        # Now the route's clean placeholder works
        DBInterface.execute(db, "INSERT INTO comparisons DEFAULT VALUES")

        SQLite.close(db)
    end
end

@testset "migrate_compare_relax_nullability! is idempotent (#67)" begin
    mktempdir() do tmp
        path = joinpath(tmp, "legacy.db")
        db = _build_legacy_strict_compare_db(path)
        HimalayaUI.migrate_compare_relax_nullability!(db)
        sql_first = String(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='table' AND name='comparisons'"))).sql)
        # Second call must not change the schema.
        HimalayaUI.migrate_compare_relax_nullability!(db)
        sql_second = String(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='table' AND name='comparisons'"))).sql)
        @test sql_first == sql_second
        # Third call from within open_db (which always runs migrate_schema!)
        # must also be a no-op — the bigger guarantee.
        SQLite.close(db)
        db2 = HimalayaUI.open_db(path)
        sql_third = String(first(Tables.rowtable(DBInterface.execute(db2,
            "SELECT sql FROM sqlite_master WHERE type='table' AND name='comparisons'"))).sql)
        @test sql_first == sql_third
        SQLite.close(db2)
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

@testset "migrate_samples_naming! — duplicate names suffixed by ascending id" begin
    mktempdir() do tmp
        db = SQLite.DB(joinpath(tmp, "h.db"))
        DBInterface.execute(db, "PRAGMA foreign_keys = ON")
        DBInterface.execute(db, """CREATE TABLE experiments (
            id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, path TEXT,
            data_dir TEXT, analysis_dir TEXT, manifest_path TEXT, config TEXT,
            experiment_type TEXT, energy_kev REAL, flight_path_m REAL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP)""")
        DBInterface.execute(db, """CREATE TABLE samples (
            id INTEGER PRIMARY KEY AUTOINCREMENT, experiment_id INTEGER REFERENCES experiments(id),
            label TEXT, name TEXT, notes TEXT)""")
        DBInterface.execute(db, "INSERT INTO experiments (id, name) VALUES (1, 'exp')")
        # Three rows that collide post-COALESCE on (1, "JC001").
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, label, name) VALUES (?, ?, ?, ?)",
            [10, 1, "JC001", "v1"])
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, label, name) VALUES (?, ?, ?, ?)",
            [11, 1, "JC001", "v2"])
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, label, name) VALUES (?, ?, ?, ?)",
            [12, 1, "JC001", "v3"])

        HimalayaUI.migrate_samples_naming!(db)

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name FROM samples ORDER BY id"))
        @test rows[1].name == "JC001"
        @test rows[2].name == "JC001-2"
        @test rows[3].name == "JC001-3"
    end
end

@testset "migrate_samples_naming! — dup suffix avoids user-named -N collision" begin
    mktempdir() do tmp
        db = SQLite.DB(joinpath(tmp, "h.db"))
        DBInterface.execute(db, "PRAGMA foreign_keys = ON")
        DBInterface.execute(db, """CREATE TABLE experiments (
            id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, path TEXT,
            data_dir TEXT, analysis_dir TEXT, manifest_path TEXT, config TEXT,
            experiment_type TEXT, energy_kev REAL, flight_path_m REAL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP)""")
        DBInterface.execute(db, """CREATE TABLE samples (
            id INTEGER PRIMARY KEY AUTOINCREMENT, experiment_id INTEGER REFERENCES experiments(id),
            label TEXT, name TEXT, notes TEXT)""")
        DBInterface.execute(db, "INSERT INTO experiments (id, name) VALUES (1, 'exp')")
        # Two duplicates of JC001 PLUS a user-named JC001-2 already in the table.
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, label, name) VALUES (?, ?, ?, ?)",
            [10, 1, "JC001", "v1"])
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, label, name) VALUES (?, ?, ?, ?)",
            [11, 1, "JC001", "v2"])
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, label, name) VALUES (?, ?, ?, ?)",
            [12, 1, "JC001-2", "user-named"])

        HimalayaUI.migrate_samples_naming!(db)

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name FROM samples ORDER BY id"))
        @test rows[1].name == "JC001"        # oldest of the dup pair keeps bare name
        @test rows[2].name == "JC001-3"      # second dup gets suffix-3 (skipping -2)
        @test rows[3].name == "JC001-2"      # user-named sample preserved
    end
end

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

@testset "migrate_experiment_config_label_to_name! — legacy [manifest].label blob rewritten in place" begin
    mktempdir() do tmp
        dbpath = joinpath(tmp, "h.db")
        db = SQLite.DB(dbpath)
        DBInterface.execute(db, "PRAGMA foreign_keys = ON")
        DBInterface.execute(db, """CREATE TABLE experiments (
            id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, path TEXT,
            data_dir TEXT, analysis_dir TEXT, manifest_path TEXT, config TEXT)""")
        legacy = """
        [experiment]
        name = "Legacy Exp"

        [manifest]
        delimiter = ","
        sample_id = 1
        label     = 2
        name      = 3
        filenames = 9

        [files]
        integration = "{name}_tot.dat"
        """
        DBInterface.execute(db,
            "INSERT INTO experiments (id, name, path, config) VALUES (?, ?, ?, ?)",
            [1, "exp", "/tmp", legacy])
        # Pre-condition: _build_config rejects the legacy blob.
        @test_throws ErrorException HimalayaUI.config_from_db(db, 1)

        HimalayaUI.migrate_experiment_config_label_to_name!(db)

        # Post: blob no longer has `[manifest].label`; has `display_name`; `name`
        # is still present (was column 3 before, now column 2 with the rewrite).
        blob = String(only(Tables.rowtable(DBInterface.execute(db,
            "SELECT config FROM experiments WHERE id = 1"))).config)
        section = ""
        seen_label = false; seen_display_name = false; seen_name_in_manifest = false
        for line in eachline(IOBuffer(blob))
            m = match(r"^\s*\[([A-Za-z0-9_]+)\]\s*$", line)
            if m !== nothing; section = m.captures[1]; continue; end
            if section == "manifest"
                occursin(r"^\s*label\s*=", line)        && (seen_label = true)
                occursin(r"^\s*display_name\s*=", line) && (seen_display_name = true)
                occursin(r"^\s*name\s*=", line)         && (seen_name_in_manifest = true)
            end
        end
        @test !seen_label
        @test seen_display_name
        @test seen_name_in_manifest

        # And config_from_db succeeds, with column indices matching the rewrite:
        # original `label = 2` ⇒ `name = 2` (stable identifier),
        # original `name = 3`  ⇒ `display_name = 3` (editable).
        cfg = HimalayaUI.config_from_db(db, 1)
        @test cfg.col_name         == 2
        @test cfg.col_display_name == 3
    end
end

@testset "migrate_experiment_config_label_to_name! — idempotent on already-migrated blob" begin
    mktempdir() do tmp
        db = SQLite.DB(joinpath(tmp, "h.db"))
        DBInterface.execute(db, """CREATE TABLE experiments (
            id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, path TEXT, config TEXT)""")
        canonical = "[manifest]\nname = 2\ndisplay_name = 3\n"
        DBInterface.execute(db,
            "INSERT INTO experiments (id, name, path, config) VALUES (?, ?, ?, ?)",
            [1, "e", "/", canonical])
        HimalayaUI.migrate_experiment_config_label_to_name!(db)
        # Bytes unchanged.
        blob = String(only(Tables.rowtable(DBInterface.execute(db,
            "SELECT config FROM experiments WHERE id = 1"))).config)
        @test blob == canonical
        # Second run still no-op.
        HimalayaUI.migrate_experiment_config_label_to_name!(db)
        @test String(only(Tables.rowtable(DBInterface.execute(db,
            "SELECT config FROM experiments WHERE id = 1"))).config) == canonical
    end
end

@testset "migrate_experiment_config_label_to_name! — NULL/missing config rows tolerated" begin
    mktempdir() do tmp
        db = SQLite.DB(joinpath(tmp, "h.db"))
        DBInterface.execute(db, """CREATE TABLE experiments (
            id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, path TEXT, config TEXT)""")
        DBInterface.execute(db,
            "INSERT INTO experiments (id, name, path, config) VALUES (?, ?, ?, NULL)",
            [1, "legacy-no-config", "/"])
        # Should not raise; config column NULL is filtered by the WHERE clause.
        @test (HimalayaUI.migrate_experiment_config_label_to_name!(db); true)
    end
end

@testset "migrate_experiment_config_label_to_name! — corrupt 'both label AND display_name' raises" begin
    mktempdir() do tmp
        db = SQLite.DB(joinpath(tmp, "h.db"))
        DBInterface.execute(db, """CREATE TABLE experiments (
            id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, path TEXT, config TEXT)""")
        corrupt = "[manifest]\nlabel = 2\nname = 3\ndisplay_name = 4\n"
        DBInterface.execute(db,
            "INSERT INTO experiments (id, name, path, config) VALUES (?, ?, ?, ?)",
            [1, "e", "/", corrupt])
        @test_throws ErrorException HimalayaUI.migrate_experiment_config_label_to_name!(db)
    end
end

@testset "open_db: prod-shape DB with legacy experiments.config blob heals end-to-end" begin
    mktempdir() do tmp
        dbpath = joinpath(tmp, "h.db")
        # Build a DB via open_db FIRST so the schema is current, then synthesize
        # the legacy blob. This mirrors the actual prod state: schema is recent
        # (samples migrated, AUTOINCREMENT etc.) but the experiments.config blob
        # is stale because PR #107 didn't include an in-DB migration.
        HimalayaUI.open_db(dbpath)  # initial open establishes schema

        # Now poison the experiments.config blob and re-open.
        db = SQLite.DB(dbpath)
        legacy = "[manifest]\nlabel = 2\nname = 3\nfilenames = 9\n"
        DBInterface.execute(db,
            "INSERT INTO experiments (id, name, path, data_dir, analysis_dir, config) VALUES (?, ?, ?, ?, ?, ?)",
            [1, "e", "/", "/d", "/a", legacy])
        close(db)

        # Re-open should run migrate_experiment_config_label_to_name! and heal it.
        db2 = HimalayaUI.open_db(dbpath)
        cfg = HimalayaUI.config_from_db(db2, 1)
        @test cfg.col_name         == 2
        @test cfg.col_display_name == 3
        close(db2)
    end
end

@testset "open_db: pre-Plan-7 legacy DB (non-AUTOINCREMENT + (label, name)) preserves identifiers" begin
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
        rows = Tables.rowtable(DBInterface.execute(db2,
            "SELECT id, name, display_name FROM samples ORDER BY id"))
        @test length(rows) == 2
        # Stable identifier (was label) preserved as name:
        @test rows[1].name == "JC001"
        @test rows[2].name == "JC002"
        # Friendly text (was name) preserved as display_name:
        @test rows[1].display_name == "DOPC + chol"
        @test rows[2].display_name == "POPC"
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
                          "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e','/p','/d','/a')")
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
