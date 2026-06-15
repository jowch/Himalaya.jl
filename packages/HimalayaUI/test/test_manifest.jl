using Test
using HimalayaUI: expand_filename_range, expand_filename_field, parse_manifest,
    ManifestSample

@testset "expand_filename_range" begin
    @test expand_filename_range("JC001-004") == ["JC001", "JC002", "JC003", "JC004"]
    @test expand_filename_range("JC013-JC016") == ["JC013", "JC014", "JC015", "JC016"]
    @test expand_filename_range("JC001") == ["JC001"]
end

@testset "expand_filename_field: multi-range" begin
    # Single range / single filename still works
    @test expand_filename_field("JC001-004") == ["JC001", "JC002", "JC003", "JC004"]
    @test expand_filename_field("JC001") == ["JC001"]

    # Semicolon-separated multiple ranges (no CSV quoting needed)
    @test expand_filename_field("JC037-040;JC153-156;JC161-164") ==
        ["JC037", "JC038", "JC039", "JC040",
         "JC153", "JC154", "JC155", "JC156",
         "JC161", "JC162", "JC163", "JC164"]

    # Comma-separated multiple ranges (cell already CSV-unquoted by parser)
    @test expand_filename_field("JC037-JC040,JC153-156,JC161-164") ==
        ["JC037", "JC038", "JC039", "JC040",
         "JC153", "JC154", "JC155", "JC156",
         "JC161", "JC162", "JC163", "JC164"]

    # Mixed range styles (with/without prefix on second bound) and stray spaces
    @test expand_filename_field(" JC001-002 ; JC005-JC006 ") ==
        ["JC001", "JC002", "JC005", "JC006"]

    # Empty segments are skipped
    @test expand_filename_field(";JC001;;JC003;") == ["JC001", "JC003"]

    # Edge cases: empty cell and all-delimiter cell both yield no filenames
    @test expand_filename_field("") == String[]
    @test expand_filename_field(" , ; ") == String[]
end

const MANIFEST_CSV = """
#\tSample\tName\tType\tTime(s)\t\t#\t\tFilename(s)\tNotes (Sample)\tNotes (Exposure)
\tFlight Path: DNA, 0.7 m, Capillaries\t\t\t\t\t\t\t\t\t
1\tD1\tUX1\tControl\t20\t\t\t\tJC001-004\tclear\t
2\tD2\tUX2\tControl\t20\t\t\t\tJC005-008\tclear\t
3\tD3\tUL1\tControl\t20\t\t\t\tJC009-JC012\tclear\t
4\tD4\tUL2\tSample\t20\t\t\t\tJC013-JC016\tcondensed\tsq
"""

@testset "parse_manifest" begin
    samples = parse_manifest(IOBuffer(MANIFEST_CSV))

    @test length(samples) == 4

    s1 = samples[1]
    @test s1.name         == "D1"
    @test s1.display_name == "UX1"
    @test s1.notes_sample == "clear"
    @test s1.filenames == ["JC001", "JC002", "JC003", "JC004"]

    s3 = samples[3]
    @test s3.filenames == ["JC009", "JC010", "JC011", "JC012"]

    s4 = samples[4]
    @test s4.name           == "D4"
    @test s4.notes_sample   == "condensed"
    @test s4.notes_exposure == "sq"
    @test s4.filenames == ["JC013", "JC014", "JC015", "JC016"]
end

# Both fixtures below mirror the user's real /data/ssrl/2026_01/0p7/ files:
# the same logical sample with multiple non-contiguous filename ranges, expressed
# either as `;`-separated ranges (no quoting needed) or as `,`-separated ranges
# with the cell wrapped in CSV double quotes.
const MANIFEST_MULTIRANGE_SEMI = """
#,Sample,Name,Type,Time(s),,#,,Filename(s),Notes (Sample),Notes (Exposure)
1,D1,XX,Control,,,,,JC001-004,,
10,D10,B1,Sample,,,,,JC037-040;JC153-156;JC161-164,,Try to get stronger peaks; spin
"""

const MANIFEST_MULTIRANGE_QUOTED = """
#,Sample,Name,Type,Time(s),,#,,Filename(s),Notes (Sample),Notes (Exposure)
1,D1,XX,Control,,,,,JC001-JC004,,
10,D10,B1,Sample,,,,,"JC037-JC040,JC153-156,JC161-164",,Try to get stronger peaks; spin
"""

@testset "parse_manifest: multi-range filename cells" begin
    cfg = HimalayaUI.ExperimentConfig(
        "x", "", "manifest.csv",
        nothing, nothing, "A-1",
        nothing, nothing, nothing,   # beam_center_x/y, pixel_size_um
        ",", 1, 0,
        1, 2, 3, 9, 10, 11,
        "data", "analysis", "simple",
        "{name}.dat", "{name}.tiff",
    )

    for csv in (MANIFEST_MULTIRANGE_SEMI, MANIFEST_MULTIRANGE_QUOTED)
        samples = parse_manifest(cfg, IOBuffer(csv))
        @test length(samples) == 2

        @test samples[1].name == "D1"
        @test samples[1].filenames == ["JC001", "JC002", "JC003", "JC004"]

        @test samples[2].name         == "D10"
        @test samples[2].display_name == "B1"
        @test samples[2].filenames == [
            "JC037", "JC038", "JC039", "JC040",
            "JC153", "JC154", "JC155", "JC156",
            "JC161", "JC162", "JC163", "JC164",
        ]
        @test samples[2].notes_exposure == "Try to get stronger peaks; spin"
    end
end

@testset "parse_manifest: multi-range with named-header columns" begin
    # Same multi-range payload, but resolved via header_row=1 + String column
    # names instead of positional Int indices. Exercises the CSV.jl named-
    # column code path with the new field-expander.
    cfg = HimalayaUI.ExperimentConfig(
        "x", "", "manifest.csv",
        nothing, nothing, "A-1",
        nothing, nothing, nothing,   # beam_center_x/y, pixel_size_um
        ",", 0, 1,
        "#", "Sample", "Name", "Filename(s)",
        "Notes (Sample)", "Notes (Exposure)",
        "data", "analysis", "simple",
        "{name}.dat", "{name}.tiff",
    )
    samples = parse_manifest(cfg, IOBuffer(MANIFEST_MULTIRANGE_QUOTED))
    @test length(samples) == 2
    @test samples[2].name == "D10"
    @test samples[2].filenames == [
        "JC037", "JC038", "JC039", "JC040",
        "JC153", "JC154", "JC155", "JC156",
        "JC161", "JC162", "JC163", "JC164",
    ]
end

@testset "parse_manifest — column 2 is identifier (name), column 3 is display_name" begin
    cfg = HimalayaUI.load_builtin_config("simple")
    csv = """sample_id\tname\tdisplay_name\tcol4\tcol5\tcol6\tcol7\tcol8\tfilenames\tnotes_sample\tnotes_exposure
1\tJC001\tDOPC + chol\t\t\t\t\t\tJC001\t\t
2\tJC002\tPOPC\t\t\t\t\t\tJC002\t\t
"""
        # skip_rows=1 with header_row=0 means line 1 is skipped, parsing starts at line 2.
    samples = HimalayaUI.parse_manifest(cfg, IOBuffer(csv))
    @test length(samples) == 2
    @test samples[1].name == "JC001"
    @test samples[1].display_name == "DOPC + chol"
    @test samples[2].name == "JC002"
    @test samples[2].display_name == "POPC"
end
