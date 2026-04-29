using Test, HimalayaUI

@testset "load_config" begin
    mktempdir() do dir
        toml = joinpath(dir, "experiment.toml")
        write(toml, """
        [experiment]
        name        = "TestRun/ExpA"
        description = "test"
        manifest    = "manifest.csv"

        [beamline]
        energy_kev    = 12.0
        flight_path_m = 2.5

        [manifest]
        delimiter      = "\\t"
        skip_rows      = 1
        header_row     = 0
        sample_id      = 1
        label          = 2
        name           = 3
        filenames      = 9
        notes_sample   = 10
        notes_exposure = 11

        [layout]
        data_dir      = "data"
        analysis_dir  = "analysis/automatic_analysis"
        exposure_type = "simple"

        [files]
        integration = "{name}.dat"
        image       = "{name}.tiff"
        """)
        cfg = HimalayaUI.load_config(toml)
        @test cfg.name == "TestRun/ExpA"
        @test cfg.energy_kev == 12.0
        @test cfg.flight_path_m == 2.5
        @test cfg.manifest_file == "manifest.csv"
        @test cfg.delimiter == "\t"
        @test cfg.skip_rows == 1
        @test cfg.header_row == 0
        @test cfg.col_sample_id == 1
        @test cfg.col_label == 2
        @test cfg.col_name == 3
        @test cfg.col_filenames == 9
        @test cfg.col_notes_sample == 10
        @test cfg.col_notes_exposure == 11
        @test cfg.data_dir == "data"
        @test cfg.analysis_dir == "analysis/automatic_analysis"
        @test cfg.exposure_type == "simple"
        @test cfg.integration_pattern == "{name}.dat"
        @test cfg.image_pattern == "{name}.tiff"
    end
end

@testset "load_config validates path traversal" begin
    mktempdir() do dir
        toml = joinpath(dir, "experiment.toml")
        write(toml, """
        [experiment]
        name = "X"
        description = ""
        manifest = "manifest.csv"
        [beamline]
        energy_kev = 0.0
        flight_path_m = 0.0
        [manifest]
        delimiter = "\\t"
        skip_rows = 0
        header_row = 0
        sample_id = 1
        label = 2
        name = 3
        filenames = 9
        notes_sample = 10
        notes_exposure = 11
        [layout]
        data_dir = "data"
        analysis_dir = "analysis/automatic_analysis"
        exposure_type = "simple"
        [files]
        integration = "../{name}.dat"
        image = "{name}.tiff"
        """)
        @test_throws ErrorException HimalayaUI.load_config(toml)
    end
end
