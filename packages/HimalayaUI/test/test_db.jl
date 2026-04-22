using Test, SQLite, DBInterface
using HimalayaUI: create_schema!, create_experiment!, create_sample!,
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
