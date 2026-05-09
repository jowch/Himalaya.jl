using Test
using HimalayaUI: cli_migrate_toml, migrate_manifest_toml_text

const _OLD_TOML = """
[experiment]
name        = "demo"

[beamline]
# axis label units (ASCII)
q_units = "A-1"

[manifest]
delimiter = "\\t"
skip_rows = 1
sample_id = 1
label     = 2  # column for sample identifier
name      = 3
filenames = 9

[files]
integration = "{name}.dat"
image       = "{name}.tiff"
"""

const _NEW_TOML = """
[experiment]
name        = "demo"

[beamline]
# axis label units (ASCII)
q_units = "A-1"

[manifest]
delimiter = "\\t"
skip_rows = 1
sample_id = 1
name      = 2  # column for sample identifier
display_name = 3
filenames = 9

[files]
integration = "{name}.dat"
image       = "{name}.tiff"
"""

@testset "migrate-toml — happy path" begin
    mktempdir() do dir
        path = joinpath(dir, "experiment.toml")
        write(path, _OLD_TOML)
        cli_migrate_toml([dir])
        out = read(path, String)
        @test occursin(r"^name\s*=\s*2\s+#"m, out)        # comment preserved
        @test occursin(r"^display_name\s*=\s*3"m, out)
        @test !occursin(r"^\s*label\s*="m, out)            # old key gone
        # The beamline comment "axis label units" must not be misfired:
        @test occursin("# axis label units", out)
    end
end

@testset "migrate-toml — idempotent on already-migrated file" begin
    mktempdir() do dir
        path = joinpath(dir, "experiment.toml")
        write(path, _NEW_TOML)
        @test_logs (:info, r"already migrated") cli_migrate_toml([dir])
        @test read(path, String) == _NEW_TOML  # unchanged
    end
end

@testset "migrate-toml — errors on file with both old and new keys" begin
    mktempdir() do dir
        path = joinpath(dir, "experiment.toml")
        write(path, """
[manifest]
sample_id    = 1
label        = 2
display_name = 3
filenames    = 9
""")
        @test_throws ErrorException cli_migrate_toml([dir])
    end
end

@testset "migrate-toml — errors on missing experiment.toml" begin
    mktempdir() do dir
        @test_throws ErrorException cli_migrate_toml([dir])
    end
end

@testset "migrate_manifest_toml_text — newline preservation matches source byte-for-byte" begin
    # A blob with NO trailing newline must come out with NO trailing newline
    # (rewritten line follows the same rule). This pins the behavior change
    # vs the old cli_migrate_toml that always appended `\n`. Important for
    # the in-DB blob path where we round-trip the bytes through SQLite.
    no_newline = "[manifest]\nlabel = 2\nname = 3"
    out, changed = migrate_manifest_toml_text(no_newline)
    @test changed
    @test out == "[manifest]\nname = 2\ndisplay_name = 3"
    @test !endswith(out, '\n')

    # And a blob WITH trailing newline keeps it.
    with_newline = "[manifest]\nlabel = 2\nname = 3\n"
    out2, changed2 = migrate_manifest_toml_text(with_newline)
    @test changed2
    @test out2 == "[manifest]\nname = 2\ndisplay_name = 3\n"
end
