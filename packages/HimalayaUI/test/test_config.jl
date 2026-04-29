using Test, HimalayaUI

# Helper: build a minimal valid TOML body, allowing overrides for the
# `[files]` and `[manifest]` sections so individual test cases can poke at
# specific validator branches.
function _make_toml(; integration="{name}.dat", image="{name}.tiff", sample_id="1")
    """
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
    sample_id = $sample_id
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
    integration = "$integration"
    image = "$image"
    """
end

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

@testset "load_config validates patterns" begin
    @testset "rejects upward traversal (..)" begin
        mktempdir() do dir
            toml = joinpath(dir, "experiment.toml")
            write(toml, _make_toml(integration="../{name}.dat"))
            @test_throws ErrorException HimalayaUI.load_config(toml)
        end
    end

    @testset "rejects absolute path" begin
        mktempdir() do dir
            toml = joinpath(dir, "experiment.toml")
            write(toml, _make_toml(integration="/abs/{name}.dat"))
            @test_throws ErrorException HimalayaUI.load_config(toml)
        end
    end

    @testset "rejects pattern missing {name}" begin
        mktempdir() do dir
            toml = joinpath(dir, "experiment.toml")
            write(toml, _make_toml(integration="fixed.dat"))
            @test_throws ErrorException HimalayaUI.load_config(toml)
        end
    end

    @testset "rejects missing file" begin
        @test_throws ErrorException HimalayaUI.load_config("/nonexistent/experiment.toml")
    end
end

@testset "load_config rejects non-int/non-string column" begin
    mktempdir() do dir
        toml = joinpath(dir, "experiment.toml")
        # `sample_id = 1.5` is neither Integer nor String — should hit _coerce_col's error.
        write(toml, _make_toml(sample_id="1.5"))
        @test_throws ErrorException HimalayaUI.load_config(toml)
    end
end

@testset "load_config accepts string column header" begin
    mktempdir() do dir
        toml = joinpath(dir, "experiment.toml")
        write(toml, _make_toml(sample_id="\"SampleID\""))
        cfg = HimalayaUI.load_config(toml)
        @test cfg.col_sample_id == "SampleID"
    end
end

@testset "load_builtin_config simple" begin
    cfg = HimalayaUI.load_builtin_config("simple")
    @test cfg.delimiter == "\t"
    @test cfg.col_sample_id == 1
    @test cfg.data_dir == "data"
    @test cfg.analysis_dir == "analysis/automatic_analysis"
    @test cfg.integration_pattern == "{name}.dat"
    @test cfg.image_pattern == "{name}.tiff"
    @test cfg.exposure_type == "simple"
end

@testset "list_config_types includes simple" begin
    types = HimalayaUI.list_config_types()
    @test "simple" in types
end

@testset "load_builtin_config rejects unknown type" begin
    @test_throws ErrorException HimalayaUI.load_builtin_config("nonexistent_type_xyz")
end

@testset "resolve_files" begin
    mktempdir() do dir
        for name in ["JC001", "JC002", "JC003", "AB001"]
            write(joinpath(dir, name * ".dat"), "")
        end

        cfg = HimalayaUI.load_builtin_config("simple")

        # Single stem prefix: JC001 → finds JC001.dat only (only file starting with JC001)
        @test HimalayaUI.resolve_files(cfg, dir, "JC001", cfg.integration_pattern) == ["JC001"]

        # Range prefix: JC002
        @test HimalayaUI.resolve_files(cfg, dir, "JC002", cfg.integration_pattern) == ["JC002"]

        # Broad prefix: JC → finds all JC*.dat sorted
        @test HimalayaUI.resolve_files(cfg, dir, "JC", cfg.integration_pattern) == ["JC001", "JC002", "JC003"]

        # Non-existent prefix → empty
        @test HimalayaUI.resolve_files(cfg, dir, "ZZ", cfg.integration_pattern) == String[]

        # Non-existent directory → empty
        @test HimalayaUI.resolve_files(cfg, "/no/such/dir", "JC", cfg.integration_pattern) == String[]
    end
end

@testset "resolve_files with subdir pattern" begin
    mktempdir() do dir
        subdir = joinpath(dir, "integrated")
        mkpath(subdir)
        write(joinpath(subdir, "SA001.dat"), "")
        write(joinpath(subdir, "SA002.dat"), "")

        cfg = HimalayaUI.load_builtin_config("simple")
        results = HimalayaUI.resolve_files(cfg, dir, "SA", "integrated/{name}.dat")
        @test results == ["SA001", "SA002"]
    end
end

@testset "resolve_files filters by suffix" begin
    mktempdir() do dir
        write(joinpath(dir, "X001.dat"), "")
        write(joinpath(dir, "X001.tiff"), "")  # different suffix
        write(joinpath(dir, "X002.dat"), "")

        cfg = HimalayaUI.load_builtin_config("simple")
        @test HimalayaUI.resolve_files(cfg, dir, "X", cfg.integration_pattern) == ["X001", "X002"]
    end
end

@testset "parse_manifest with config: positional columns" begin
    csv = join([
        "header row to skip",
        "1\tD1\tUX1\tT\tt\t\t\t\tJC001-003\tnote_s\tnote_e",
        "2\tD2\tUX2\tT\tt\t\t\t\tJC004\t\t",
        "not_a_number\tskip\tme\tT\tt\t\t\t\tJC999\t\t",
    ], "\n")
    cfg = HimalayaUI.load_builtin_config("simple")
    samples = HimalayaUI.parse_manifest(cfg, IOBuffer(csv))

    @test length(samples) == 2
    @test samples[1].label == "D1"
    @test samples[1].name  == "UX1"
    @test samples[1].filenames == ["JC001", "JC002", "JC003"]
    @test samples[1].notes_sample   == "note_s"
    @test samples[1].notes_exposure == "note_e"
    @test samples[2].filenames == ["JC004"]
    @test samples[2].notes_sample   == ""
    @test samples[2].notes_exposure == ""
end

@testset "parse_manifest with config: named columns" begin
    csv = join([
        "#,Sample,Name,Type,Time,x,y,z,Filename(s),Notes (Sample),Notes (Exposure)",
        "1,D1,UX1,,,,,,AB001-002,s_note,e_note",
        "2,D2,UX2,,,,,,AB003,,",
    ], "\n")

    mktempdir() do dir
        toml_path = joinpath(dir, "experiment.toml")
        write(toml_path, """
        [experiment]
        name = "X"
        description = ""
        manifest = "manifest.csv"
        [beamline]
        energy_kev = 0.0
        flight_path_m = 0.0
        [manifest]
        delimiter = ","
        skip_rows = 0
        header_row = 1
        sample_id = "#"
        label = "Sample"
        name = "Name"
        filenames = "Filename(s)"
        notes_sample = "Notes (Sample)"
        notes_exposure = "Notes (Exposure)"
        [layout]
        data_dir = "data"
        analysis_dir = "analysis/automatic_analysis"
        exposure_type = "simple"
        [files]
        integration = "{name}.dat"
        image = "{name}.tiff"
        """)
        cfg = HimalayaUI.load_config(toml_path)

        samples = HimalayaUI.parse_manifest(cfg, IOBuffer(csv))
        @test length(samples) == 2
        @test samples[1].label    == "D1"
        @test samples[1].filenames == ["AB001", "AB002"]
        @test samples[1].notes_sample   == "s_note"
        @test samples[1].notes_exposure == "e_note"
        @test samples[2].filenames == ["AB003"]
    end
end

@testset "parse_manifest backward-compat wrapper still works" begin
    csv = "header\n1\tD1\tUX1\tT\tt\t\t\t\tJC001-002\t\t"
    samples_old = HimalayaUI.parse_manifest(IOBuffer(csv))
    samples_new = HimalayaUI.parse_manifest(HimalayaUI.load_builtin_config("simple"), IOBuffer(csv))
    @test length(samples_old) == length(samples_new) == 1
    @test samples_old[1].filenames == samples_new[1].filenames
end
