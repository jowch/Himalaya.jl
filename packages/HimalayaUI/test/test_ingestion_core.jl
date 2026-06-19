# packages/HimalayaUI/test/test_ingestion_core.jl
using Test
using SQLite, DBInterface, Tables
using Dates
using HimalayaUI

# ---------------------------------------------------------------------------
# Shared fixture helpers
# ---------------------------------------------------------------------------

"""
Write a synthetic PRP file to `path` with the given field overrides.
Defaults match the real HA_85_422_S2404_0_001.prp from the SSRL 2026-04 run.

Note: scan id and frame number are NOT PRP fields — they live in the *filename
stem* (`_S<digits>_…_<frame>`) and are recovered by `_parse_scan_frame` (Task 7),
so `write_prp` takes no `scan_id`/`frame_no` kwargs. Tests exercise scan/frame
parsing by naming the fixture file accordingly (e.g. `HA_85_422_S2404_0_001.prp`).
"""
function write_prp(path::AbstractString;
        timestamp = "26 Apr 2026 18:20:27",
        beam_energy_ev = 9000.027604502573,
        pipe_length_mm = 1700,
        detector = "Pilatus 1M",
        exposure_time = 15,
        horizontal_position_mm = 58.9)
    open(path, "w") do io
        println(io, "Image file name: $(basename(path)[1:end-4]).tif")
        println(io, "Time this file was written: $timestamp")
        println(io, "Exposure time=$exposure_time")
        println(io, "Counting time=$(exposure_time + 0.8)")
        println(io, "Beam energy=$beam_energy_ev eV")
        println(io, "Pipe length=$pipe_length_mm mm")
        println(io, "Scan motor=motor")
        println(io, "Scan range=delta mm")
        println(io, "Horizontal position=$horizontal_position_mm mm")
        println(io, "Vertical position=26.55 mm")
        println(io, "dispx position=9.8 mm")
        println(io, "dispy position=35.9 mm")
        println(io, "Detector=$detector")
        println(io, "Detector mode=Pilatus; 0 for normal and 1 for dezingered")
        println(io, "Phi position=157")
        println(io, "Attenuator=0 um Al")
    end
end

"""
Write a synthetic setup_info file to `path`.
Defaults match the real setup_info_20260425_181705.txt from SSRL 2026-04.
"""
function write_setup_info(path::AbstractString;
        beam_center_x = 421.409,
        beam_center_y = 836.946,
        mean_distance_mm = 1809.5,
        ring_distances = [1806.0, 1811.1, 1811.5])
    open(path, "w") do io
        println(io, "File used for calibration:")
        println(io, "../data/agbe_1p7m_S1958_0_001.tif")
        println(io, "")
        println(io, "Beam center is at ($beam_center_x, $beam_center_y)")
        println(io, "")
        println(io, "Calculated detector distances are:")
        for (i, d) in enumerate(ring_distances)
            println(io, "Ring $i:\t pixel $(247.82 * i):\t $d mm")
        end
        println(io, "----------------------------------------------------------")
        println(io, "Mean distance:\t\t $mean_distance_mm mm")
    end
end

"Open a fresh DB with the full migrated schema."
function fresh_db()
    path = joinpath(mktempdir(), "h.db")
    HimalayaUI.open_db(path)
end

@testset "ingestion core (Phase B)" begin
    @testset "detector pitch lookup" begin
        # Known detectors
        @test HimalayaUI.detector_pixel_size_um("Pilatus 1M") ≈ 172.0
        @test HimalayaUI.detector_pixel_size_um("Pilatus 2M") ≈ 172.0
        @test HimalayaUI.detector_pixel_size_um("Pilatus 300K") ≈ 172.0
        @test HimalayaUI.detector_pixel_size_um("Eiger 1M")  ≈ 75.0
        @test HimalayaUI.detector_pixel_size_um("Eiger 4M")  ≈ 75.0
        # Unknown detector → missing
        @test HimalayaUI.detector_pixel_size_um("Unknown XRD") === missing
    end

    @testset "parse_setup_info" begin
        dir = mktempdir()
        setup_path = joinpath(dir, "setup_info_20260425_181705.txt")
        write_setup_info(setup_path;
            beam_center_x = 421.409,
            beam_center_y = 836.946,
            mean_distance_mm = 1809.5)
        info = HimalayaUI.parse_setup_info(setup_path)
        @test info.beam_center_x ≈ 421.409
        @test info.beam_center_y ≈ 836.946
        @test info.mean_distance_m ≈ 1.8095

        # Missing Mean distance line → missing
        bad_path = joinpath(dir, "setup_bad.txt")
        write(bad_path, "Beam center is at (100.0, 200.0)\n")
        bad_info = HimalayaUI.parse_setup_info(bad_path)
        @test bad_info.beam_center_x ≈ 100.0
        @test bad_info.mean_distance_m === missing

        # Completely empty → all missing
        empty_path = joinpath(dir, "setup_empty.txt")
        write(empty_path, "")
        empty_info = HimalayaUI.parse_setup_info(empty_path)
        @test empty_info.beam_center_x === missing
        @test empty_info.mean_distance_m === missing
    end

    @testset "derive_geometry" begin
        dir = mktempdir()
        data_dir = joinpath(dir, "data")
        mkpath(data_dir)
        analysis_dir = joinpath(dir, "analysis")
        mkpath(analysis_dir)

        # Two PRP files with consistent geometry
        for (name, hpos) in [("HA_001", 58.9), ("HA_002", 63.1)]
            write_prp(joinpath(data_dir, "$name.prp");
                beam_energy_ev = 9000.027604502573,
                pipe_length_mm = 1700,
                detector = "Pilatus 1M",
                exposure_time = 15.0,
                horizontal_position_mm = hpos)
        end

        # One setup file with calibrated distance
        write_setup_info(joinpath(analysis_dir, "setup_info_20260425_181705.txt");
            beam_center_x = 421.409, beam_center_y = 836.946, mean_distance_mm = 1809.5)

        prp_paths   = [joinpath(data_dir, "HA_001.prp"), joinpath(data_dir, "HA_002.prp")]
        setup_files = [joinpath(analysis_dir, "setup_info_20260425_181705.txt")]

        geo, discrepancies = HimalayaUI.derive_geometry(prp_paths, setup_files)

        # Energy from PRP
        @test geo.energy_kev ≈ 9.000027604502573
        @test geo.energy_kev_source == "prp"

        # flight_path_m from calibrated setup file (authority per spec §6)
        @test geo.flight_path_m ≈ 1.8095
        @test geo.flight_path_m_source == "setup"

        # Beam center from setup file
        @test geo.beam_center_x ≈ 421.409
        @test geo.beam_center_y ≈ 836.946
        @test geo.beam_center_x_source == "setup"

        # Pixel size from detector→pitch lookup
        @test geo.pixel_size_um ≈ 172.0
        @test geo.pixel_size_um_source == "prp"

        # No discrepancies when fields are consistent
        @test isempty(discrepancies)

        # Discrepancy when PRP beam energies disagree
        write_prp(joinpath(data_dir, "HA_003.prp");
            beam_energy_ev = 8900.0,    # different!
            pipe_length_mm = 1700, detector = "Pilatus 1M",
            horizontal_position_mm = 68.0)
        geo2, disc2 = HimalayaUI.derive_geometry(
            vcat(prp_paths, [joinpath(data_dir, "HA_003.prp")]), setup_files)
        @test any(d -> d.field == "beam_energy_ev", disc2)

        # No setup file → fall back to PRP pipe length
        geo3, _ = HimalayaUI.derive_geometry(prp_paths, String[])
        @test geo3.flight_path_m ≈ 1.700
        @test geo3.flight_path_m_source == "prp"

        # Unknown detector → pixel_size_um missing, discrepancy flagged
        write_prp(joinpath(data_dir, "HA_004.prp");
            beam_energy_ev = 9000.0, pipe_length_mm = 1700,
            detector = "SuperDetector X",
            horizontal_position_mm = 72.0)
        geo4, disc4 = HimalayaUI.derive_geometry([joinpath(data_dir, "HA_004.prp")], String[])
        @test geo4.pixel_size_um === missing
        @test any(d -> d.field == "pixel_size_um", disc4)
    end

    @testset "scan_directory" begin
        dir = mktempdir()
        data_dir = joinpath(dir, "data")
        analysis_dir = joinpath(dir, "analysis")
        mkpath(data_dir); mkpath(analysis_dir)

        # Write 6 PRP files with 3 distinct horizontal positions (2 frames per slot)
        stems = ["HA_01_S001_0_001", "HA_01_S002_0_001",   # slot 1 (H ≈ 58.9)
                 "HA_02_S003_0_001", "HA_02_S004_0_001",   # slot 2 (H ≈ 63.1)
                 "HA_03_S005_0_001", "HA_03_S006_0_001"]   # slot 3 (H ≈ 67.3)
        h_positions = [58.9, 58.91, 63.1, 63.09, 67.3, 67.31]
        for (stem, hpos) in zip(stems, h_positions)
            write_prp(joinpath(data_dir, "$stem.prp");
                horizontal_position_mm = hpos)
            write(joinpath(data_dir, "$stem.tif"), "fake tif")
            write(joinpath(analysis_dir, "$stem.dat"), "fake dat")
        end

        # Write one setup_info file
        write_setup_info(joinpath(analysis_dir, "setup_info_20260425_181705.txt"))

        metas = HimalayaUI.scan_directory(data_dir, analysis_dir)

        # All 6 stems found
        @test length(metas) == 6
        # Each has a prp_path, tif_path, dat_path, and parsed prp fields
        m = first(metas)
        @test haskey(m, :stem)
        @test haskey(m, :prp_path)
        @test haskey(m, :tif_path)
        @test haskey(m, :prp)  # the parse_prp result
        # Stems returned sorted
        @test [m.stem for m in metas] == sort([m.stem for m in metas])
        # A missing PRP (has tif but no prp) → prp is missing
        write(joinpath(data_dir, "HA_04_S007_0_001.tif"), "fake tif")
        metas2 = HimalayaUI.scan_directory(data_dir, analysis_dir)
        @test length(metas2) == 7
        tif_only = first(filter(m -> m.stem == "HA_04_S007_0_001", metas2))
        @test tif_only.prp_path === nothing
    end

    @testset "load segmentation" begin
        # Build 8 ExposureMeta with synthetic timestamps:
        #   4 in Load 1 at T+0, T+19, T+38, T+57
        #   2578s gap (reload)
        #   4 in Load 2 at T+2635, T+2654, T+2673, T+2692
        t0 = DateTime(2026, 4, 26, 23, 14, 8)
        offsets = [0, 19, 38, 57, 2635, 2654, 2673, 2692]  # seconds
        metas = [HimalayaUI.ExposureMeta(
            "HA_$(lpad(i,2,'0'))_S$(lpad(i,4,'0'))_0_001",
            nothing, nothing, nothing,
            (timestamp = t0 + Second(s), horizontal_position_mm = 58.9 + Float64(i)*4.0,
             beam_energy_ev = 9000.0, energy_kev = 9.0, pipe_length_m = 1.7,
             detector = "Pilatus 1M", exposure_time_s = 15.0),
        ) for (i, s) in enumerate(offsets)]

        loads = HimalayaUI._segment_loads(metas)
        @test length(loads) == 2
        @test length(loads[1]) == 4
        @test length(loads[2]) == 4
        @test loads[1][1].stem == "HA_01_S0001_0_001"
        @test loads[2][1].stem == "HA_05_S0005_0_001"

        # Unimodal fallback: all gaps equal → one load, returned with flag
        uniform_metas = [HimalayaUI.ExposureMeta(
            "HA_$(lpad(i,2,'0'))_S$(lpad(i,4,'0'))_0_001",
            nothing, nothing, nothing,
            (timestamp = t0 + Second(i * 20), horizontal_position_mm = 58.9 + Float64(i)*4.0,
             beam_energy_ev = 9000.0, energy_kev = 9.0, pipe_length_m = 1.7,
             detector = "Pilatus 1M", exposure_time_s = 15.0),
        ) for i in 1:4]
        uni_loads, uni_flag = HimalayaUI._segment_loads_with_flag(uniform_metas)
        @test length(uni_loads) == 1
        @test uni_flag == :unimodal_fallback

        # Single exposure → one load, no-crash
        single = [first(metas)]
        sl, sf = HimalayaUI._segment_loads_with_flag(single)
        @test length(sl) == 1
        @test sf == :ok
    end

    @testset "parse_prp" begin
        dir = mktempdir()
        prp_path = joinpath(dir, "HA_85_422_S2404_0_001.prp")
        write_prp(prp_path;
            timestamp = "26 Apr 2026 18:20:27",
            beam_energy_ev = 9000.027604502573,
            pipe_length_mm = 1700,
            detector = "Pilatus 1M",
            exposure_time = 15.0,
            horizontal_position_mm = 58.9)
        prp = HimalayaUI.parse_prp(prp_path)

        # Energy: strip "eV", convert to keV
        @test prp.beam_energy_ev ≈ 9000.027604502573
        @test prp.energy_kev ≈ 9.000027604502573

        # Pipe length: strip "mm", store as metres
        @test prp.pipe_length_m ≈ 1.700

        # Detector model string (drives pitch lookup in geometry.jl)
        @test prp.detector == "Pilatus 1M"

        # Exposure time (seconds)
        @test prp.exposure_time_s ≈ 15.0

        # Horizontal position (mm)
        @test prp.horizontal_position_mm ≈ 58.9

        # Timestamp parses to a DateTime
        @test prp.timestamp isa DateTime
        @test year(prp.timestamp) == 2026
        @test month(prp.timestamp) == 4

        # Graceful missing: a truncated PRP returns missing for absent fields
        trunc_path = joinpath(dir, "trunc.prp")
        write(trunc_path, "Image file name: trunc.tif\nBeam energy=9000 eV\n")
        trunc_prp = HimalayaUI.parse_prp(trunc_path)
        @test trunc_prp.horizontal_position_mm === missing
        @test trunc_prp.detector === missing
    end

    @testset "slot clustering" begin
        # Simulate one load: 3 slots × 4 frames each, H positions stepped ~4 mm apart,
        # within-burst jitter ≈ 0.3 mm. Mirrors the real HA data pattern from spec §5.
        t0 = DateTime(2026, 4, 26, 23, 14, 8)
        slot_h = [74.80, 70.85, 67.22]  # three slot centers
        jitter = [0.0, 0.13, -0.24, 0.05]  # within-burst jitter, 4 frames
        frames = HimalayaUI.ExposureMeta[]
        for (si, h_center) in enumerate(slot_h)
            for (fi, j) in enumerate(jitter)
                offset = (si - 1) * 4 * 19 + (fi - 1) * 19
                stem = "HA_$(si)_$(fi)_S$(lpad((si-1)*4+fi, 4,'0'))_0_001"
                push!(frames, HimalayaUI.ExposureMeta(stem, nothing, nothing, nothing,
                    (timestamp = t0 + Second(offset),
                     horizontal_position_mm = h_center + j,
                     beam_energy_ev = 9000.0, energy_kev = 9.0, pipe_length_m = 1.7,
                     detector = "Pilatus 1M", exposure_time_s = 15.0)))
            end
        end

        slots = HimalayaUI._cluster_slots(frames)
        @test length(slots) == 3
        @test all(length(s) == 4 for s in slots)

        # All 4 frames of slot 1 should have H ≈ 74.80
        slot1_h = [m.prp.horizontal_position_mm for m in slots[1]]
        @test all(abs(h - 74.80) < 0.5 for h in slot1_h)

        # Median-delta-near-zero fallback (the `med_delta < 1e-6` branch).
        #
        # This exercises the documented algorithm intent: within-slot revisits collapse the
        # MEDIAN consecutive-frame delta to ~0 (most consecutive frames sit at the *same*
        # position), so the plain `med_delta × slot_k` tolerance is degenerate (0) and would
        # split on every frame. The fallback instead learns the slot-spacing tolerance from
        # the non-zero deltas (the burst→burst jumps). Fixture: 3 slots, each a 3-frame burst
        # at a fixed position (consecutive in-burst deltas == 0), with ~4 mm jumps between
        # bursts. The 6 in-burst deltas are 0 and the 2 between-burst deltas are ~4 mm, so
        # median(all 8) == 0 → fallback branch → tolerance = median(nonzero)≈3.8 mm / 5 ≈
        # 0.76 mm, which sits below the ~4 mm burst jumps (so they split) and above the
        # 0 mm within-burst jitter (so the bursts stay intact). 3 bursts → 3 slots.
        burst_frames = HimalayaUI.ExposureMeta[]
        for (si, h_center) in enumerate([74.80, 70.85, 67.22])
            for _ in 1:3
                push!(burst_frames, HimalayaUI.ExposureMeta("b_$(si)_$(length(burst_frames))",
                    nothing, nothing, nothing,
                    (timestamp = t0 + Second(length(burst_frames) * 19),
                     horizontal_position_mm = h_center,   # identical within the burst → zero deltas
                     beam_energy_ev = 9000.0, energy_kev = 9.0, pipe_length_m = 1.7,
                     detector = "Pilatus 1M", exposure_time_s = 15.0)))
            end
        end

        slots_bf = HimalayaUI._cluster_slots(burst_frames)
        @test length(slots_bf) == 3                       # one slot per burst position
        @test all(length(s) == 3 for s in slots_bf)       # 3 frames each
    end

    @testset "group_into_samples" begin
        # Two-load synthetic directory: 2 loads × 3 slots × 2 frames
        t0 = DateTime(2026, 4, 26, 23, 14, 8)
        slot_h_l1 = [74.80, 70.85, 67.22]
        slot_h_l2 = [74.75, 70.90, 67.18]  # same positions in load 2 (reload)
        frames = HimalayaUI.ExposureMeta[]
        global_s = 1
        for (li, slot_h) in enumerate([slot_h_l1, slot_h_l2])
            load_start = li == 1 ? 0 : 600  # 600 s reload gap
            for (si, h) in enumerate(slot_h)
                for fi in 1:2
                    offset = load_start + (si - 1) * 2 * 19 + (fi - 1) * 19
                    stem = "HA_$(li)_$(si)_$(fi)_S$(lpad(global_s, 4,'0'))_0_001"
                    global_s += 1
                    push!(frames, HimalayaUI.ExposureMeta(stem, "/d/$stem.tif", "/d/$stem.prp", "/a/$stem.dat",
                        (timestamp = t0 + Second(offset),
                         horizontal_position_mm = h + (fi == 2 ? 0.1 : 0.0),
                         beam_energy_ev = 9000.0, energy_kev = 9.0, pipe_length_m = 1.7,
                         detector = "Pilatus 1M", exposure_time_s = 15.0)))
                end
            end
        end

        result = HimalayaUI.group_into_samples(frames)

        # Structure: 2 loads
        @test length(result.loads) == 2
        # 3 slots per load
        @test all(length(l.samples) == 3 for l in result.loads)
        # 2 frames per sample
        @test all(length(s.exposures) == 2 for l in result.loads for s in l.samples)

        # Auto-naming: first sample of load 1, slot 1 → "HA ... (S01P01)"
        name11 = result.loads[1].samples[1].name
        @test occursin("S01P01", name11)
        # Second load: S02P01
        name21 = result.loads[2].samples[1].name
        @test occursin("S02P01", name21)

        # Load indices 1-based
        @test result.loads[1].load_index == 1
        @test result.loads[2].load_index == 2

        # Frame counts
        @test result.loads[1].frame_count == 6  # 3 slots × 2 frames
        @test result.loads[2].frame_count == 6

        # Load time range
        @test result.loads[1].start_time isa DateTime
        @test result.loads[1].end_time   isa DateTime
    end
end
