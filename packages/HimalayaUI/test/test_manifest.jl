using Test
using HimalayaUI: expand_filename_range, parse_manifest, ManifestSample

@testset "expand_filename_range" begin
    @test expand_filename_range("JC001-004") == ["JC001", "JC002", "JC003", "JC004"]
    @test expand_filename_range("JC013-JC016") == ["JC013", "JC014", "JC015", "JC016"]
    @test expand_filename_range("JC001") == ["JC001"]
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
    @test s1.label == "D1"
    @test s1.name  == "UX1"
    @test s1.notes_sample == "clear"
    @test s1.filenames == ["JC001", "JC002", "JC003", "JC004"]

    s3 = samples[3]
    @test s3.filenames == ["JC009", "JC010", "JC011", "JC012"]

    s4 = samples[4]
    @test s4.label == "D4"
    @test s4.notes_sample   == "condensed"
    @test s4.notes_exposure == "sq"
    @test s4.filenames == ["JC013", "JC014", "JC015", "JC016"]
end
