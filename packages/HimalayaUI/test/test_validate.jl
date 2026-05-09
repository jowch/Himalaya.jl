using Test
using HimalayaUI: ManifestSample, ManifestViolation, validate_manifest

@testset "validate_manifest — empty name" begin
    samples = [ManifestSample("", "DOPC", "", "", ["JC001"])]
    vs = validate_manifest(samples)
    @test length(vs) == 1
    @test vs[1].kind == :empty_name
    @test vs[1].sample_index == 1
end

@testset "validate_manifest — bad name characters" begin
    samples = [ManifestSample("JC 001", "DOPC", "", "", ["JC001"])]  # space
    vs = validate_manifest(samples)
    @test length(vs) == 1
    @test vs[1].kind == :bad_name_chars
    @test occursin("JC 001", vs[1].detail)
end

@testset "validate_manifest — duplicate name" begin
    samples = [
        ManifestSample("JC001", "first",  "", "", ["JC001"]),
        ManifestSample("JC001", "second", "", "", ["JC002"]),
    ]
    vs = validate_manifest(samples)
    dup = filter(v -> v.kind == :duplicate_name, vs)
    @test length(dup) == 1
end

@testset "validate_manifest — duplicate filenames within sample" begin
    samples = [ManifestSample("JC001", "x", "", "", ["JC001", "JC001"])]
    vs = validate_manifest(samples)
    @test any(v -> v.kind == :duplicate_filename_in_sample, vs)
end

@testset "validate_manifest — overlapping filenames across samples" begin
    samples = [
        ManifestSample("JC001", "x", "", "", ["JC001", "JC002"]),
        ManifestSample("JC002", "y", "", "", ["JC002", "JC003"]),
    ]
    vs = validate_manifest(samples)
    @test any(v -> v.kind == :overlapping_filenames, vs)
end

@testset "validate_manifest — clean manifest produces no violations" begin
    samples = [
        ManifestSample("JC001", "first",  "", "", ["JC001"]),
        ManifestSample("JC002", "second", "", "", ["JC002"]),
    ]
    @test isempty(validate_manifest(samples))
end

@testset "validate_manifest — collects all violations, no fail-fast" begin
    samples = [
        ManifestSample("",       "x", "", "", ["JC001"]),         # empty
        ManifestSample("JC 002", "y", "", "", ["JC001", "JC001"]) # bad chars + dup filenames
    ]
    vs = validate_manifest(samples)
    @test length(vs) >= 3   # at least: empty_name, bad_name_chars, duplicate_filename_in_sample
end
