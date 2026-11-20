using Test, HimalayaUI, SQLite

function _pattern_db(; kwargs...)
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.create_experiment!(db;
        name = "X", path = "/tmp/x", data_dir = "/tmp/x/data",
        analysis_dir = "/tmp/x/analysis", kwargs...)
    db, exp_id
end

@testset "integration_pattern_for" begin
    @testset "defaults to {name}.dat when the column and the blob are unset" begin
        db, id = _pattern_db()
        @test HimalayaUI.integration_pattern_for(db, id) == "{name}.dat"
    end

    @testset "uses the experiment's column" begin
        db, id = _pattern_db(integration_pattern = "{name}_tot.dat")
        @test HimalayaUI.integration_pattern_for(db, id) == "{name}_tot.dat"
    end

    @testset "falls back to [files].integration in a legacy config blob" begin
        db, id = _pattern_db(config = "[files]\nintegration = \"integrated/{name}.dat\"\n")
        @test HimalayaUI.integration_pattern_for(db, id) == "integrated/{name}.dat"
    end

    @testset "the column wins over the legacy blob" begin
        db, id = _pattern_db(integration_pattern = "{name}_tot.dat",
                             config = "[files]\nintegration = \"{name}.dat\"\n")
        @test HimalayaUI.integration_pattern_for(db, id) == "{name}_tot.dat"
    end

    @testset "a column without {name} falls through to the legacy blob" begin
        db, id = _pattern_db(integration_pattern = "trace.dat",
                             config = "[files]\nintegration = \"integrated/{name}.dat\"\n")
        @test HimalayaUI.integration_pattern_for(db, id) == "integrated/{name}.dat"
    end

    @testset "ignores values without {name}" begin
        db, id = _pattern_db(integration_pattern = "trace.dat",
                             config = "[files]\nintegration = \"also.dat\"\n")
        @test HimalayaUI.integration_pattern_for(db, id) == "{name}.dat"
    end

    @testset "tolerates a malformed legacy blob" begin
        db, id = _pattern_db(config = "not [ valid toml")
        @test HimalayaUI.integration_pattern_for(db, id) == "{name}.dat"
    end

    @testset "errors on an unknown experiment" begin
        db, _ = _pattern_db()
        @test_throws ErrorException HimalayaUI.integration_pattern_for(db, 999)
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
