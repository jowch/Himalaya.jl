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
end
