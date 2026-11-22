using Test, HimalayaUI, SQLite

function _pattern_db(; kwargs...)
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.create_experiment!(db;
        name = "X", path = "/tmp/x", data_dir = "/tmp/x/data",
        analysis_dir = "/tmp/x/analysis", kwargs...)
    db, exp_id
end

_exposure_pattern(db, id) = HimalayaUI.integration_pattern_for_exposure(db,
    HimalayaUI.create_exposure!(db; experiment_id = id, filename = "f"))

@testset "integration_pattern_for" begin
    @testset "defaults to {name}.dat when the column is unset" begin
        db, id = _pattern_db()
        @test _exposure_pattern(db, id) == "{name}.dat"
    end

    @testset "uses the setup's column" begin
        db, id = _pattern_db(integration_pattern = "{name}_tot.dat")
        @test _exposure_pattern(db, id) == "{name}_tot.dat"
    end

    @testset "refuses a new value without {name} and ignores a stored one" begin
        @test_throws ArgumentError _pattern_db(integration_pattern = "trace.dat")
        db, id = _pattern_db()
        SQLite.DBInterface.execute(db, "UPDATE setups SET integration_pattern = 'trace.dat'")
        @test _exposure_pattern(db, id) == "{name}.dat"
    end

    @testset "errors on an exposure without a setup" begin
        db, _ = _pattern_db()
        @test_throws ErrorException HimalayaUI.integration_pattern_for_exposure(db, 999)
    end
end

@testset "default_db_path respects HIMALAYA_DB_PATH" begin
    old = get(ENV, "HIMALAYA_DB_PATH", nothing)
    try
        ENV["HIMALAYA_DB_PATH"] = "/tmp/himalaya-test/central.db"
        @test HimalayaUI.default_db_path() == "/tmp/himalaya-test/central.db"
    finally
        old === nothing ? delete!(ENV, "HIMALAYA_DB_PATH") :
                          (ENV["HIMALAYA_DB_PATH"] = old)
    end
end

@testset "open_db creates parent directory" begin
    mktempdir() do dir
        nested = joinpath(dir, "deeply", "nested", "himalaya.db")
        db = HimalayaUI.open_db(nested)
        @test isfile(nested)
        close(db)
    end
end
