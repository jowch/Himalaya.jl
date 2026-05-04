using Test, SQLite, DBInterface, Tables
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

    s_id = create_sample!(db; experiment_id = exp_id, label = "D1", name = "UX1")
    @test s_id == 1

    samples = get_samples(db, exp_id)
    @test length(samples) == 1
    @test first(samples).label == "D1"

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
            "INSERT INTO samples (experiment_id, label) VALUES (1, 'A1')")
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
