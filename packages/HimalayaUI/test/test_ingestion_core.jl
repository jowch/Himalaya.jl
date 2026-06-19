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

        # Two PRP files with consistent geometry. The nominal pipe length (1700 mm)
        # differs from the setup calibrated distance (1809.5 mm) by ~6.4%, but that
        # expected calibration gap is NOT flagged as a discrepancy (option A).
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

        # No setup file → fall back to PRP pipe length (HA_001/HA_002 have pipe_length=1700mm)
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
        @test hasproperty(m, :stem)
        @test hasproperty(m, :prp_path)
        @test hasproperty(m, :tif_path)
        @test hasproperty(m, :prp)  # the parse_prp result
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

    @testset "grouping regression floor — SSRL 2026-04 Load 1 fixture" begin
        using Dates

        # Mirror the real 1p7m run structure (Load 1 of 13):
        # 4 slots, 4 frames each, H steps ~3.95 mm, within-burst jitter ≈ 0.3 mm,
        # 19 s between frames, 15 s exposure time (from PRP).
        # Confirmed parameters from real data sweep 2026-06-18.
        dir = mktempdir()
        data_dir     = joinpath(dir, "data")
        analysis_dir = joinpath(dir, "analysis")
        mkpath(data_dir); mkpath(analysis_dir)

        SLOT_CENTERS_MM = [74.80, 70.85, 67.22, 63.49]  # from real HA data
        WITHIN_BURST_JITTER = [0.0, 0.13, -0.24, 0.05]
        FRAME_GAP_S = 19
        t0 = DateTime(2026, 4, 25, 23, 14, 8)

        scan_id = 2001
        stems = String[]
        for (si, h_center) in enumerate(SLOT_CENTERS_MM)
            for (fi, jitter) in enumerate(WITHIN_BURST_JITTER)
                offset = (si - 1) * length(WITHIN_BURST_JITTER) * FRAME_GAP_S + (fi - 1) * FRAME_GAP_S
                stem = "HA_$(si)_$(lpad(scan_id, 4,'0'))_S$(scan_id)_0_001"
                scan_id += 1
                push!(stems, stem)
                write_prp(joinpath(data_dir, "$stem.prp");
                    timestamp        = Dates.format(t0 + Second(offset), "dd u yyyy HH:MM:SS"),
                    beam_energy_ev   = 9000.027604502573,
                    pipe_length_mm   = 1700,
                    detector         = "Pilatus 1M",
                    exposure_time    = 15.0,
                    horizontal_position_mm = h_center + jitter)
                write(joinpath(data_dir, "$stem.tif"), "fake tif")
            end
        end
        write_setup_info(joinpath(analysis_dir, "setup_info_20260425_181705.txt");
            beam_center_x = 421.409, beam_center_y = 836.946, mean_distance_mm = 1809.5)

        metas = HimalayaUI.scan_directory(data_dir, analysis_dir)

        # Regression floors for scan_directory
        @test length(metas) >= 16          # ≥ all 16 exposures found
        @test length(metas) <= 16          # ≤ only the 16 we wrote (no phantom)

        result = HimalayaUI.group_into_samples(metas)

        # One load (19 s gaps are well within one load)
        @test length(result.loads) >= 1
        # At least 4 slots detected
        total_samples = sum(length(l.samples) for l in result.loads)
        @test total_samples >= 4           # floor: found at least 4 slots
        @test total_samples <= 6           # ceiling: shouldn't massively over-split
        # Frame count matches
        total_frames = sum(l.frame_count for l in result.loads)
        @test total_frames == 16

        # Geometry derivation
        prp_paths   = filter(!isnothing, [m.prp_path for m in metas])
        setup_files = [joinpath(analysis_dir, "setup_info_20260425_181705.txt")]
        geo, disc   = HimalayaUI.derive_geometry(String.(prp_paths), String.(setup_files))

        @test geo.energy_kev ≈ 9.000027604502573
        @test geo.flight_path_m ≈ 1.8095            # calibrated, not PRP nominal
        @test geo.flight_path_m_source == "setup"
        @test geo.beam_center_x ≈ 421.409
        @test geo.pixel_size_um ≈ 172.0

        # The ~6.4% nominal-vs-calibrated flight-path gap is EXPECTED calibration
        # physics and is deliberately NOT flagged as a discrepancy (option A).
        @test !any(d -> occursin("flight_path_m_nominal_vs_calibrated", d.field), disc)

        # Auto-naming: every sample name contains the coordinate anchor
        for load in result.loads
            for (si, samp) in enumerate(load.samples)
                @test occursin("P$(lpad(si, 2,'0'))", samp.name)
            end
        end
    end

    @testset "scan_and_group! inserts loads/samples/exposures" begin
        dir = mktempdir()
        data_dir     = joinpath(dir, "data")
        analysis_dir = joinpath(dir, "analysis")
        mkpath(data_dir); mkpath(analysis_dir)

        # 2 slots × 2 frames in one load
        slot_h = [58.9, 63.1]
        t0 = DateTime(2026, 4, 26, 23, 14, 8)
        all_stems = String[]
        for (si, h) in enumerate(slot_h)
            for fi in 1:2
                offset = (si - 1) * 2 * 19 + (fi - 1) * 19
                stem = "HA_$(si)_$(fi)_S$(lpad((si-1)*2+fi, 4,'0'))_0_001"
                push!(all_stems, stem)
                write_prp(joinpath(data_dir, "$stem.prp");
                    timestamp = Dates.format(t0 + Second(offset), "dd u yyyy HH:MM:SS"),
                    beam_energy_ev = 9000.027604502573,
                    pipe_length_mm = 1700, detector = "Pilatus 1M",
                    exposure_time = 15.0, horizontal_position_mm = h)
                write(joinpath(data_dir, "$stem.tif"), "fake tif")
                # No .dat files: analyze_exposure! will fail to read them, which is OK
                # for this insert-path test (we set analyze=false)
            end
        end
        write_setup_info(joinpath(analysis_dir, "setup_info_20260425_181705.txt"))

        db = fresh_db()
        # Create a bare experiment (Phase-A create_experiment! signature)
        exp_id = HimalayaUI.create_experiment!(db;
            name = "test-ingest", path = dir,
            data_dir = data_dir, analysis_dir = analysis_dir)

        HimalayaUI.scan_and_group!(db, exp_id; analyze = false)

        # Exposures in DB
        exps = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, experiment_id, sample_id, filename, image_path, scan_id, frame_no FROM exposures WHERE experiment_id = ?",
            [exp_id]))
        @test length(exps) == 4
        @test all(e.experiment_id == exp_id for e in exps)
        # After grouping, sample_id is populated
        @test all(!ismissing(e.sample_id) for e in exps)
        # image_path persisted (non-NULL, equals the written .tif path) — the image route
        # and prewarm_thumbnails! both filter `WHERE image_path IS NOT NULL`.
        @test all(!ismissing(e.image_path) for e in exps)
        @test all(endswith(String(e.image_path), ".tif") for e in exps)
        # scan_id / frame_no parsed from the filename stem (`_S<digits>_…_<frame>`) and persisted (spec §4/§10)
        @test all(!ismissing(e.scan_id) for e in exps)
        @test all(e.frame_no == 1 for e in exps)   # every fixture stem ends `_001`

        # Samples: 2 (one per slot)
        samps = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name, load_id, slot_index FROM samples WHERE experiment_id = ?", [exp_id]))
        @test length(samps) == 2
        @test all(s.slot_index in [1, 2] for s in samps)
        @test any(occursin("S01P01", s.name) for s in samps)

        # Loads: 1
        loads = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, load_index, frame_count FROM loads WHERE experiment_id = ?", [exp_id]))
        @test length(loads) == 1
        @test loads[1].load_index == 1
        @test loads[1].frame_count == 4

        # Geometry written to experiments row
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT energy_kev, flight_path_m, flight_path_m_source, beam_center_x FROM experiments WHERE id = ?",
            [exp_id])))
        @test row.energy_kev ≈ 9.000027604502573
        @test row.flight_path_m ≈ 1.8095   # from setup file
        @test row.flight_path_m_source == "setup"
        @test row.beam_center_x ≈ 421.409

        # Idempotent: a second scan of the same directory is a true no-op. Insert-only
        # discipline applies to loads AND samples AND exposures — re-running must not mint
        # a fresh load_id (which would re-key the sample dedup and duplicate samples).
        HimalayaUI.scan_and_group!(db, exp_id; analyze = false)
        exps2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM exposures WHERE experiment_id = ?", [exp_id]))
        @test length(exps2) == 4  # exposures unchanged
        loads2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM loads WHERE experiment_id = ?", [exp_id]))
        @test length(loads2) == 1  # loads unchanged (no orphan duplicate load row)
        samps2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM samples WHERE experiment_id = ?", [exp_id]))
        @test length(samps2) == 2  # samples unchanged (sample dedup keyed on a stable load_id)
    end

    @testset "derive_sample_flags" begin
        # Helper to build a rollup-shaped exposure NamedTuple (mirrors get_loads_rollup).
        exp(id, fname, hpos) = (id = id, filename = fname,
            horizontal_position = hpos, timestamp = nothing)
        smp(sid, name, slot, exps) = (sample_id = sid, name = name, slot_index = slot,
            grouping_source = "auto_position", name_source = "auto",
            merged_into_id = nothing, exposures = exps)
        load(lid, lidx, smps) = (load_id = lid, load_index = lidx,
            start_time = nothing, end_time = nothing,
            frame_count = sum(length(s.exposures) for s in smps), samples = smps)

        # --- Cross-load MERGE candidate: label "HA_85" recurs across two loads ---
        # Load 1 slot 1: HA_85 burst (clean, ~0.1 mm jitter) → merge candidate (recurs in L2)
        # Load 2 slot 1: HA_85 burst (clean)                  → merge candidate (recurs in L1)
        # Load 1 slot 2: HA_90 burst (no recurrence)          → no flag
        s1 = smp(101, "HA85 (S01P01)", 1,
                 [exp(1, "HA_85_S2001_0_001", 74.80), exp(2, "HA_85_S2002_0_002", 74.83)])
        s2 = smp(102, "HA90 (S01P02)", 2,
                 [exp(3, "HA_90_S2003_0_001", 70.85), exp(4, "HA_90_S2004_0_002", 70.88)])
        s3 = smp(201, "HA85 (S02P01)", 1,
                 [exp(5, "HA_85_S2101_0_001", 74.75), exp(6, "HA_85_S2102_0_002", 74.79)])

        rollup_merge = [load(1, 1, [s1, s2]), load(2, 1, [s3])]
        flags = HimalayaUI.derive_sample_flags(rollup_merge)

        # s1 and s3 are merge candidates of each other (shared label "HA_85")
        @test haskey(flags, 101)
        @test flags[101] isa HimalayaUI.MergeFlag
        @test flags[101].merge_with_sample_id == 201
        @test flags[101].merge_with_label == "HA_85"
        @test haskey(flags, 201)
        @test flags[201].merge_with_sample_id == 101
        # s2 (label HA_90, no recurrence) is NOT flagged → absent from the Dict (contract null)
        @test !haskey(flags, 102)

        # --- Intra-sample SPLIT: one sample spans a ~4 mm position jump ---
        # Load 1 has only this sample, so there is no cross-load recurrence;
        # within the sample the position jumps 74.80 → 67.20 (≫ local jitter).
        split_sample = smp(301, "HA77 (S01P01)", 1,
            [exp(10, "HA_77_S3001_0_001", 74.80),
             exp(11, "HA_77_S3002_0_002", 74.85),   # still ~slot center
             exp(12, "HA_77_S3003_0_003", 67.20),   # JUMP here (index 3)
             exp(13, "HA_77_S3004_0_004", 67.25)])
        rollup_split = [load(3, 1, [split_sample])]
        sflags = HimalayaUI.derive_sample_flags(rollup_split)
        @test haskey(sflags, 301)
        @test sflags[301] isa HimalayaUI.SplitFlag
        @test sflags[301].split_at_index == 3            # the exposure index where the jump occurs
        @test sflags[301].jump_from ≈ 74.85
        @test sflags[301].jump_to   ≈ 67.20

        # --- Pure single-frame round-robin (the prior Plan-B deferral, now a concrete split) ---
        # Each "exposure" is a different slot position visited in one round; _cluster_slots
        # under-split it to ONE sample (KNOWN LIMITATION). derive_sample_flags flags the
        # first position jump as a split suggestion (spec §5: "single-frame … else flag").
        rr_sample = smp(401, "JC (S01P01)", 1,
            [exp(20, "JC_C01_S4001_0_001", 39.5),
             exp(21, "JC_C02_S4002_0_001", 71.0),    # jump → split at index 2
             exp(22, "JC_C03_S4003_0_001", 87.0)])
        rr_flags = HimalayaUI.derive_sample_flags([load(4, 1, [rr_sample])])
        @test haskey(rr_flags, 401)
        @test rr_flags[401] isa HimalayaUI.SplitFlag
        @test rr_flags[401].split_at_index == 2

        # --- No flags for a clean directory (one load, one clean burst) ---
        clean = smp(501, "HA60 (S01P01)", 1,
            [exp(30, "HA_60_S5001_0_001", 58.9), exp(31, "HA_60_S5002_0_002", 58.95)])
        @test isempty(HimalayaUI.derive_sample_flags([load(5, 1, [clean])]))
    end

    @testset "cheap_change_check" begin
        dir = mktempdir()
        data_dir     = joinpath(dir, "data")
        analysis_dir = joinpath(dir, "analysis")
        mkpath(data_dir); mkpath(analysis_dir)

        t0 = DateTime(2026, 4, 26, 23, 14, 8)
        # Write 2 TIF/PRP pairs
        for (i, h) in enumerate([58.9, 63.1])
            stem = "HA_$(i)_S$(lpad(i, 4,'0'))_0_001"
            write_prp(joinpath(data_dir, "$stem.prp");
                timestamp = Dates.format(t0 + Second((i-1)*19), "dd u yyyy HH:MM:SS"),
                beam_energy_ev = 9000.0, pipe_length_mm = 1700,
                detector = "Pilatus 1M", exposure_time = 15.0,
                horizontal_position_mm = h)
            write(joinpath(data_dir, "$stem.tif"), "fake tif")
        end
        write_setup_info(joinpath(analysis_dir, "setup_info_20260425_181705.txt"))

        db = fresh_db()
        exp_id = HimalayaUI.create_experiment!(db;
            name = "cheap-check", path = dir,
            data_dir = data_dir, analysis_dir = analysis_dir)

        # Before any ingest: 2 files on disk, 0 exposures in DB → changed
        @test HimalayaUI.cheap_change_check(db, exp_id) == true

        # After ingest: 2 files, 2 exposures → unchanged (true no-op tick)
        HimalayaUI.scan_and_group!(db, exp_id; analyze = false)
        @test HimalayaUI.cheap_change_check(db, exp_id) == false

        # Add a new TIF/PRP pair on disk (not yet ingested) → changed again
        new_stem = "HA_3_S0003_0_001"
        write_prp(joinpath(data_dir, "$new_stem.prp");
            timestamp = Dates.format(t0 + Second(600), "dd u yyyy HH:MM:SS"),
            beam_energy_ev = 9000.0, pipe_length_mm = 1700,
            detector = "Pilatus 1M", exposure_time = 15.0,
            horizontal_position_mm = 67.3)
        write(joinpath(data_dir, "$new_stem.tif"), "fake tif")
        @test HimalayaUI.cheap_change_check(db, exp_id) == true

        # Re-scan picks up the new file; check goes quiet again
        HimalayaUI.scan_and_group!(db, exp_id; analyze = false)
        @test HimalayaUI.cheap_change_check(db, exp_id) == false

        # Defensive: a non-existent data dir is treated as "no change" (nothing to ingest),
        # never an error (the scheduler tick must not crash on a vanished volume).
        bad_db = fresh_db()
        bad_id = HimalayaUI.create_experiment!(bad_db;
            name = "missing-dir", path = dir,
            data_dir = joinpath(dir, "does_not_exist"),
            analysis_dir = analysis_dir)
        @test HimalayaUI.cheap_change_check(bad_db, bad_id) == false
    end
end
