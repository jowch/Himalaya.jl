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
end
