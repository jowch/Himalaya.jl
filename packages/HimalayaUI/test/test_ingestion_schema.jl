# packages/HimalayaUI/test/test_ingestion_schema.jl
using Test
using SQLite, DBInterface, Tables
using HimalayaUI

"Build a fresh DB at the CURRENT (post-migration) schema and return its path."
function fresh_db()
    dir = mktempdir()
    path = joinpath(dir, "test.db")
    db = SQLite.DB(path)
    HimalayaUI.create_schema!(db)    # canonical tables
    HimalayaUI.migrate_schema!(db)   # NOTE: takes an SQLite.DB, not a path (db.jl:247)
    SQLite.close(db)
    return path
end

"Insert a minimal experiment honouring its NOT NULL columns (name/path/data_dir/analysis_dir, db.jl:21-23). Returns the id."
function seed_experiment(db; name="e")
    DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES (?, ?, ?, ?)",
        [name, "/exp/$name", "/exp/$name/data", "/exp/$name/analysis"]))
end

"Open a path, run f(db), always close."
function with_db(f, path)
    db = SQLite.DB(path)
    try
        return f(db)
    finally
        SQLite.close(db)
    end
end

"Column names of a table, lowercased."
cols(db, table) = lowercase.(getproperty.(Tables.rowtable(
    DBInterface.execute(db, "PRAGMA table_info($table)")), :name))

"Names of indexes on a table."
indexes(db, table) = String.(getproperty.(Tables.rowtable(
    DBInterface.execute(db, "PRAGMA index_list($table)")), :name))

@testset "ingestion schema (Phase A)" begin
    # tasks append @testset blocks below

    @testset "migrations recorded" begin
        path = fresh_db()
        with_db(path) do db
            applied = Set(String.(getproperty.(Tables.rowtable(
                DBInterface.execute(db, "SELECT name FROM schema_migrations")), :name)))
            @test HimalayaUI.MIGRATION_LOADS_TABLE in applied
            @test HimalayaUI.MIGRATION_EXPOSURES_EXPERIMENT_ID in applied
            @test HimalayaUI.MIGRATION_EXPERIMENTS_GEOMETRY in applied
            @test HimalayaUI.MIGRATION_SAMPLES_NAME_COLLAPSE in applied
        end
    end

    @testset "loads table" begin
        path = fresh_db()
        with_db(path) do db
            @test "loads" in lowercase.(String.(getproperty.(Tables.rowtable(
                DBInterface.execute(db, "SELECT name FROM sqlite_master WHERE type='table'")), :name)))
            c = cols(db, "loads")
            for col in ["id","experiment_id","load_index","session_id","start_time","end_time","frame_count","note"]
                @test col in c
            end
            # FK to experiments enforced
            DBInterface.execute(db, "PRAGMA foreign_keys=ON")
            @test_throws Exception DBInterface.execute(db,
                "INSERT INTO loads (experiment_id, load_index) VALUES (999999, 0)")
        end
    end

    @testset "exposures experiment_id + dedup" begin
        path = fresh_db()
        with_db(path) do db
            c = cols(db, "exposures")
            for col in ["experiment_id","prp_path","timestamp","exposure_time",
                        "horizontal_position","scan_id","frame_no","load_id","content_fingerprint"]
                @test col in c
            end
            # new dedup index present, old one gone
            idx = indexes(db, "exposures")
            @test any(i -> occursin("exposures_unique_filename", i), idx)
            # the new index is on (experiment_id, filename): same filename, different sample, same experiment → rejected
            DBInterface.execute(db, "PRAGMA foreign_keys=ON")
            eid = seed_experiment(db)   # honours NOT NULL name/path/data_dir/analysis_dir
            s1 = DBInterface.lastrowid(DBInterface.execute(db,
                "INSERT INTO samples (experiment_id, name) VALUES (?, 'A')", [eid]))
            s2 = DBInterface.lastrowid(DBInterface.execute(db,
                "INSERT INTO samples (experiment_id, name) VALUES (?, 'B')", [eid]))
            DBInterface.execute(db,
                "INSERT INTO exposures (experiment_id, sample_id, filename) VALUES (?, ?, 'f.tif')", [eid, s1])
            @test_throws Exception DBInterface.execute(db,
                "INSERT INTO exposures (experiment_id, sample_id, filename) VALUES (?, ?, 'f.tif')", [eid, s2])
        end
    end

    # Revert the exposures-denorm migration on an already-migrated DB so a
    # re-run of migrate_schema! exercises the real backfill/preflight branches
    # against hand-seeded legacy rows. Drops the sentinel + new (experiment_id,
    # filename) index, restores the OLD (sample_id, filename) index, and clears
    # experiment_id so the backfill fires again.
    function _revert_exposures_denorm!(db)
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?",
            [HimalayaUI.MIGRATION_EXPOSURES_EXPERIMENT_ID])
        DBInterface.execute(db, "DROP INDEX IF EXISTS exposures_unique_filename")
        DBInterface.execute(db, "DROP INDEX IF EXISTS exposures_experiment_idx")
        DBInterface.execute(db, "DROP INDEX IF EXISTS exposures_load_idx")
        DBInterface.execute(db,
            "CREATE UNIQUE INDEX exposures_unique_filename ON exposures(sample_id, filename)")
        DBInterface.execute(db, "UPDATE exposures SET experiment_id = NULL")
    end

    @testset "migration fails on cross-sample (experiment_id, filename) collision" begin
        path = fresh_db()
        with_db(path) do db
            _revert_exposures_denorm!(db)
            eid = seed_experiment(db)
            s1 = DBInterface.lastrowid(DBInterface.execute(db,
                "INSERT INTO samples (experiment_id, name) VALUES (?, 'A')", [eid]))
            s2 = DBInterface.lastrowid(DBInterface.execute(db,
                "INSERT INTO samples (experiment_id, name) VALUES (?, 'B')", [eid]))
            # Two samples, same experiment, same filename — legal under the old key,
            # a collision under the new (experiment_id, filename) key.
            DBInterface.execute(db,
                "INSERT INTO exposures (sample_id, filename) VALUES (?, 'dup.tif')", [s1])
            DBInterface.execute(db,
                "INSERT INTO exposures (sample_id, filename) VALUES (?, 'dup.tif')", [s2])
            err = nothing
            @test_throws ErrorException try
                HimalayaUI.migrate_schema!(db)
            catch e
                err = e
                rethrow()
            end
            @test occursin("dup.tif", sprint(showerror, err))
        end
    end

    @testset "migration fails on orphaned exposure (no derivable experiment_id)" begin
        path = fresh_db()
        with_db(path) do db
            _revert_exposures_denorm!(db)
            # sample_id NULL → nothing to backfill experiment_id from.
            DBInterface.execute(db,
                "INSERT INTO exposures (sample_id, filename) VALUES (NULL, 'orphan.tif')")
            @test_throws ErrorException HimalayaUI.migrate_schema!(db)
        end
    end

    @testset "experiments geometry+scan columns" begin
        path = fresh_db()
        with_db(path) do db
            c = cols(db, "experiments")
            for col in ["beam_center_x","beam_center_y","pixel_size_um","q_units",
                        "energy_kev_source","flight_path_m_source","beam_center_x_source",
                        "beam_center_y_source","pixel_size_um_source","q_units_source",
                        "last_scanned_at","scan_signature","ingest_status",
                        "last_scan_tier","consecutive_empty_ticks"]
                @test col in c
            end
            @test !("detector_distance_mm" in c)   # spec §10: no separate column
        end
    end

    @testset "samples name collapse + grouping cols" begin
        path = fresh_db()
        with_db(path) do db
            c = cols(db, "samples")
            @test "name" in c
            @test !("display_name" in c)            # collapsed away
            for col in ["load_id","slot_index","grouping_source","name_source"]
                @test col in c
            end
            # label is NOT unique: two samples in one experiment may share a name
            DBInterface.execute(db, "PRAGMA foreign_keys=ON")
            eid = seed_experiment(db)   # honours NOT NULL name/path/data_dir/analysis_dir
            DBInterface.execute(db, "INSERT INTO samples (experiment_id, name) VALUES (?, 'HA85 (S01P15)')", [eid])
            @test_nowarn DBInterface.execute(db,
                "INSERT INTO samples (experiment_id, name) VALUES (?, 'HA85 (S01P15)')", [eid])
            @test !any(i -> occursin("samples_unique_name", i), indexes(db, "samples"))
        end
    end

    @testset "create_experiment! geometry kwargs" begin
        path = fresh_db()
        with_db(path) do db
            id = HimalayaUI.create_experiment!(db; name="geo",
                path="/exp/geo", data_dir="/exp/geo/data", analysis_dir="/exp/geo/analysis",
                energy_kev=9.0, energy_kev_source="prp",
                flight_path_m=1.8095, flight_path_m_source="setup",
                beam_center_x=421.409, beam_center_y=836.946, beam_center_x_source="setup",
                pixel_size_um=172.0, q_units="A^-1")
            row = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT flight_path_m, flight_path_m_source, beam_center_x FROM experiments WHERE id=?", [id])))
            @test row.flight_path_m ≈ 1.8095
            @test row.flight_path_m_source == "setup"
            @test row.beam_center_x ≈ 421.409
        end
    end

    @testset "create_sample! grouping fields, no display_name" begin
        path = fresh_db()
        with_db(path) do db
            eid = seed_experiment(db)   # honours NOT NULL name/path/data_dir/analysis_dir
            lid = DBInterface.lastrowid(DBInterface.execute(db,
                "INSERT INTO loads (experiment_id, load_index) VALUES (?, 1)", [eid]))
            sid = HimalayaUI.create_sample!(db; experiment_id=eid, name="HA85 (S01P15)",
                load_id=lid, slot_index=15, grouping_source="auto_position", name_source="auto")
            row = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT name, load_id, slot_index, grouping_source, name_source FROM samples WHERE id=?", [sid])))
            @test row.name == "HA85 (S01P15)"
            @test row.load_id == lid
            @test row.slot_index == 15
            @test row.name_source == "auto"
        end
    end

    @testset "duplicate labels survive naming+collapse" begin
        # Build the pre-redesign DB with create_schema! (the canonical samples shape on this
        # branch is still name+display_name+samples_unique_name — i.e. the pre-collapse legacy
        # shape). A hand-built minimal fixture is NOT a faithful legacy DB: it lacks the tables
        # the historical AUTOINCREMENT rebuild (migrate_pk_to_autoincrement!) rebuilds and copies
        # under PRAGMA foreign_keys=OFF, and a NULL experiments.data_dir trips its NOT NULL rebuild.
        # Real legacy DBs were always made by an earlier create_schema!, so this is the honest fixture.
        dir = mktempdir(); path = joinpath(dir, "duplabel.db")
        db = SQLite.DB(path)
        HimalayaUI.create_schema!(db)
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1,'e','/exp/e','/exp/e/data','/exp/e/analysis')")
        # distinct machine names (unique), SAME display_name (the future label)
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (1,1,'m1','HA85 (S01P15)')")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (2,1,'m2','HA85 (S01P15)')")
        SQLite.close(db)

        dbm = SQLite.DB(path); HimalayaUI.migrate_schema!(dbm); SQLite.close(dbm)

        with_db(path) do db2
            names = String.(getproperty.(Tables.rowtable(DBInterface.execute(db2,
                "SELECT name FROM samples ORDER BY id")), :name))
            @test names == ["HA85 (S01P15)", "HA85 (S01P15)"]   # NOT suffixed to '...-2'
            @test !("display_name" in cols(db2, "samples"))
            @test !any(i -> occursin("samples_unique_name", i), indexes(db2, "samples"))
        end
    end

    @testset "create_exposure! experiment_id + prp fields, nullable sample_id" begin
        path = fresh_db()
        with_db(path) do db
            eid = seed_experiment(db)   # honours NOT NULL name/path/data_dir/analysis_dir
            # sample_id omitted (transient pre-group state must be allowed)
            xid = HimalayaUI.create_exposure!(db; experiment_id=eid, filename="HA_85_001.tif",
                prp_path="/d/HA_85_001.prp", timestamp="2026-04-12T10:02:00",
                exposure_time=2.0, horizontal_position=12.4, scan_id=2404, frame_no=1)
            row = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT experiment_id, sample_id, horizontal_position, frame_no FROM exposures WHERE id=?", [xid])))
            @test row.experiment_id == eid
            @test row.sample_id === missing
            @test row.horizontal_position ≈ 12.4
            @test row.frame_no == 1
        end
    end
end
