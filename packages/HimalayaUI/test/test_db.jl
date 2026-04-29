using Test, SQLite, DBInterface, Tables
using HimalayaUI: create_schema!, migrate_schema!, create_experiment!, create_sample!,
                  create_exposure!, get_experiment, get_samples, get_exposures

@testset "db schema" begin
    db = SQLite.DB()  # in-memory
    create_schema!(db)

    tables = Set(r[1] for r in DBInterface.execute(db,
        "SELECT name FROM sqlite_master WHERE type='table'"))

    for t in ["users", "experiments", "samples", "sample_tags",
              "exposures", "exposure_sources", "exposure_tags",
              "peaks", "indices", "index_peaks",
              "index_groups", "index_group_members", "user_actions"]
        @test t in tables
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
