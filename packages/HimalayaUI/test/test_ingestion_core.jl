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
Defaults are invented values of the size an SSRL .prp holds.

Note: scan id and frame number are NOT PRP fields — they live in the *filename
stem* (`_S<digits>_…_<frame>`) and are recovered by `_parse_scan_frame` (Task 7),
so `write_prp` takes no `scan_id`/`frame_no` kwargs. Tests exercise scan/frame
parsing by naming the fixture file accordingly (e.g. `SY_12_034_S0567_0_001.prp`).
"""
function write_prp(path::AbstractString;
        timestamp = "14 Mar 2026 14:05:11",
        beam_energy_ev = 9000.0,
        pipe_length_mm = 1700,
        detector = "Pilatus 1M",
        exposure_time = 15,
        horizontal_position_mm = 58.9,
        i_signal = nothing,
        i_0 = nothing)
    open(path, "w") do io
        println(io, "Image file name: $(basename(path)[1:end-4]).tif")
        println(io, "Time this file was written: $timestamp")
        if i_signal !== nothing
            println(io, "I_Signal (i.e. I1 in header file) is from: i_apd")
            println(io, "I_signal=$i_signal (corrected for dark counts)")
        end
        i_0 === nothing || println(io, "I_0(daqIon3.0)=$i_0")
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
Defaults match the real setup_info_20260301_120000.txt from SSRL 2026-04.
"""
function write_setup_info(path::AbstractString;
        beam_center_x = 420.5,
        beam_center_y = 835.5,
        mean_distance_mm = 1810.0,
        ring_distances = [1807.0, 1810.5, 1812.5])
    open(path, "w") do io
        println(io, "File used for calibration:")
        println(io, "../data/agbe_1p7m_S1000_0_001.tif")
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
        setup_path = joinpath(dir, "setup_info_20260301_120000.txt")
        write_setup_info(setup_path;
            beam_center_x = 420.5,
            beam_center_y = 835.5,
            mean_distance_mm = 1810.0)
        info = HimalayaUI.parse_setup_info(setup_path)
        @test info.beam_center_x ≈ 420.5
        @test info.beam_center_y ≈ 835.5
        @test info.mean_distance_m ≈ 1.81

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

    @testset "_find_setup_files ascends from a tot_files analysis_dir" begin
        # Real SSRL layout: setup_info lives in the analysis ROOT, but the
        # experiment's analysis_dir is the .../automatic_analysis/tot_files leaf.
        root = mktempdir()
        analysis = joinpath(root, "analysis")
        tot_files = joinpath(analysis, "automatic_analysis", "tot_files")
        mkpath(tot_files)
        setup = joinpath(analysis, "setup_info_20260301_120000.txt")
        write_setup_info(setup; beam_center_x = 420.5, beam_center_y = 835.5,
                         mean_distance_mm = 1810.0)

        # From the leaf, it must ascend to the analysis root and find the file.
        found = HimalayaUI._find_setup_files(tot_files)
        @test found == [setup]
        # When analysis_dir directly holds the file (the mini layout), found there.
        @test HimalayaUI._find_setup_files(analysis) == [setup]
        # No setup file anywhere → empty (not an error).
        @test HimalayaUI._find_setup_files(mktempdir()) == String[]
        # A non-existent dir is tolerated (walks up; returns empty).
        @test HimalayaUI._find_setup_files(joinpath(root, "does", "not", "exist")) == String[]
    end

    @testset "derive_geometry" begin
        dir = mktempdir()
        data_dir = joinpath(dir, "data")
        mkpath(data_dir)
        analysis_dir = joinpath(dir, "analysis")
        mkpath(analysis_dir)

        # Two PRP files with consistent geometry. The nominal pipe length (1700 mm)
        # differs from the setup calibrated distance (1810.0 mm) by ~6.4%, but that
        # expected calibration gap is NOT flagged as a discrepancy (option A).
        for (name, hpos) in [("SY_001", 58.9), ("SY_002", 63.1)]
            write_prp(joinpath(data_dir, "$name.prp");
                beam_energy_ev = 9000.0,
                pipe_length_mm = 1700,
                detector = "Pilatus 1M",
                exposure_time = 15.0,
                horizontal_position_mm = hpos)
        end

        # One setup file with calibrated distance
        write_setup_info(joinpath(analysis_dir, "setup_info_20260301_120000.txt");
            beam_center_x = 420.5, beam_center_y = 835.5, mean_distance_mm = 1810.0)

        prp_paths   = [joinpath(data_dir, "SY_001.prp"), joinpath(data_dir, "SY_002.prp")]
        setup_files = [joinpath(analysis_dir, "setup_info_20260301_120000.txt")]

        geo, discrepancies = HimalayaUI.derive_geometry(prp_paths, setup_files)

        # Energy from PRP
        @test geo.energy_kev ≈ 9.0
        @test geo.energy_kev_source == "prp"

        # flight_path_m from calibrated setup file (authority per spec §6)
        @test geo.flight_path_m ≈ 1.81
        @test geo.flight_path_m_source == "setup"

        # Beam center from setup file
        @test geo.beam_center_x ≈ 420.5
        @test geo.beam_center_y ≈ 835.5
        @test geo.beam_center_x_source == "setup"

        # Pixel size from detector→pitch lookup
        @test geo.pixel_size_um ≈ 172.0
        @test geo.pixel_size_um_source == "prp"

        # No discrepancies when fields are consistent
        @test isempty(discrepancies)

        # Discrepancy when PRP beam energies disagree
        write_prp(joinpath(data_dir, "SY_003.prp");
            beam_energy_ev = 8900.0,    # different!
            pipe_length_mm = 1700, detector = "Pilatus 1M",
            horizontal_position_mm = 68.0)
        geo2, disc2 = HimalayaUI.derive_geometry(
            vcat(prp_paths, [joinpath(data_dir, "SY_003.prp")]), setup_files)
        @test any(d -> d.field == "beam_energy_ev", disc2)

        # No setup file → fall back to PRP pipe length (SY_001/SY_002 have pipe_length=1700mm)
        geo3, _ = HimalayaUI.derive_geometry(prp_paths, String[])
        @test geo3.flight_path_m ≈ 1.700
        @test geo3.flight_path_m_source == "prp"

        # EMPTY prp_paths must still dispatch and return an all-"default" geometry.
        # The path-based method delegates via a comprehension, and `parse_prp`'s
        # inferred return is a non-concrete NamedTuple, so a bare `[...]` collects
        # to Vector{Any} and matches NEITHER method -> MethodError. That is an
        # invisible break: the element-type annotation guarding it is the kind of
        # thing a later "simplification" removes with no other CI signal, and the
        # MethodError aborts the whole ingest. Reachable from routes_fs.jl (a
        # directory with TIFs and no PRPs).
        geo_empty, disc_empty = HimalayaUI.derive_geometry(String[], String[])
        @test geo_empty.energy_kev === missing
        @test geo_empty.energy_kev_source == "default"
        @test geo_empty.pixel_size_um_source == "default"
        @test isempty(disc_empty)
        # The pre-parsed method must accept an empty vector too.
        geo_empty2, _ = HimalayaUI.derive_geometry(NamedTuple[], String[])
        @test geo_empty2.energy_kev_source == "default"

        # Unknown detector → pixel_size_um missing, discrepancy flagged
        write_prp(joinpath(data_dir, "SY_004.prp");
            beam_energy_ev = 9000.0, pipe_length_mm = 1700,
            detector = "SuperDetector X",
            horizontal_position_mm = 72.0)
        geo4, disc4 = HimalayaUI.derive_geometry([joinpath(data_dir, "SY_004.prp")], String[])
        @test geo4.pixel_size_um === missing
        @test any(d -> d.field == "pixel_size_um", disc4)
    end

    @testset "scan_directory" begin
        dir = mktempdir()
        data_dir = joinpath(dir, "data")
        analysis_dir = joinpath(dir, "analysis")
        mkpath(data_dir); mkpath(analysis_dir)

        # Write 6 PRP files with 3 distinct horizontal positions (2 frames per slot)
        stems = ["SY_01_S001_0_001", "SY_01_S002_0_001",   # slot 1 (H ≈ 58.9)
                 "SY_02_S003_0_001", "SY_02_S004_0_001",   # slot 2 (H ≈ 63.1)
                 "SY_03_S005_0_001", "SY_03_S006_0_001"]   # slot 3 (H ≈ 67.3)
        h_positions = [58.9, 58.91, 63.1, 63.09, 67.3, 67.31]
        for (stem, hpos) in zip(stems, h_positions)
            write_prp(joinpath(data_dir, "$stem.prp");
                horizontal_position_mm = hpos)
            write(joinpath(data_dir, "$stem.tif"), "fake tif")
            write(joinpath(analysis_dir, "$stem.dat"), "fake dat")
        end

        # Write one setup_info file
        write_setup_info(joinpath(analysis_dir, "setup_info_20260301_120000.txt"))

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
        write(joinpath(data_dir, "SY_04_S007_0_001.tif"), "fake tif")
        metas2 = HimalayaUI.scan_directory(data_dir, analysis_dir)
        @test length(metas2) == 7
        tif_only = first(filter(m -> m.stem == "SY_04_S007_0_001", metas2))
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
            "SY_$(lpad(i,2,'0'))_S$(lpad(i,4,'0'))_0_001",
            nothing, nothing, nothing,
            (timestamp = t0 + Second(s), horizontal_position_mm = 58.9 + Float64(i)*4.0,
             beam_energy_ev = 9000.0, energy_kev = 9.0, pipe_length_m = 1.7,
             detector = "Pilatus 1M", exposure_time_s = 15.0),
        ) for (i, s) in enumerate(offsets)]

        loads, flag = HimalayaUI._segment_loads_with_flag(metas)
        @test flag == :ok
        @test length(loads) == 2
        @test length(loads[1]) == 4
        @test length(loads[2]) == 4
        @test loads[1][1].stem == "SY_01_S0001_0_001"
        @test loads[2][1].stem == "SY_05_S0005_0_001"

        # Unimodal fallback: all gaps equal → one load, returned with flag
        uniform_metas = [HimalayaUI.ExposureMeta(
            "SY_$(lpad(i,2,'0'))_S$(lpad(i,4,'0'))_0_001",
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
        prp_path = joinpath(dir, "SY_12_034_S0567_0_001.prp")
        write_prp(prp_path;
            timestamp = "14 Mar 2026 14:05:11",
            beam_energy_ev = 9000.0,
            pipe_length_mm = 1700,
            detector = "Pilatus 1M",
            exposure_time = 15.0,
            horizontal_position_mm = 58.9)
        prp = HimalayaUI.parse_prp(prp_path)

        # Energy: strip "eV", convert to keV
        @test prp.beam_energy_ev ≈ 9000.0
        @test prp.energy_kev ≈ 9.0

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
        @test month(prp.timestamp) == 3

        # Graceful missing: a truncated PRP returns missing for absent fields
        trunc_path = joinpath(dir, "trunc.prp")
        write(trunc_path, "Image file name: trunc.tif\nBeam energy=9000 eV\n")
        trunc_prp = HimalayaUI.parse_prp(trunc_path)
        @test trunc_prp.horizontal_position_mm === missing
        @test trunc_prp.detector === missing
        @test trunc_prp.transmitted_counts === missing
        @test trunc_prp.incident_monitor === missing
        @test trunc_prp.counting_time_s === missing
    end

    @testset "parse_prp beam monitors" begin
        dir = mktempdir()
        p = joinpath(dir, "SYN_01_001_S100_0_001.prp")
        write(p, """
            Image file name: SYN_01_001_S100_0_001.tif
            Time this file was written: 03 Mar 2026 10:15:00
            I_Signal (i.e. I1 in header file) is from: i_apd
            I_signal=4100.0 (corrected for dark counts)
            I_0(daqIon3.0)=1.2345
            I_1(daqCounter1.0)=304300.0
            I_2(daqCounter1.1)=31000000.5 (corrected for dark counts)
            I_apd(daqCounter1.0)=304300.0
            Exposure time=15
            Counting time=15.8
            i1 dark counts=19000.0 per sec
            Beam energy=9000.0 eV
            Beam Intensity Reference: i1 i2 acquisition_time
            empty(air)=300.0 12000000.0 20260101_120000
            water=40000.0 8000000 20260101_1300
            """)
        prp = HimalayaUI.parse_prp(p)
        @test prp.transmitted_counts ≈ 4100.0
        @test prp.incident_monitor ≈ 1.2345
        @test prp.counting_time_s ≈ 15.8
        @test prp.exposure_time_s ≈ 15.0

        mixed = joinpath(dir, "mixed_case.prp")
        write(mixed, "I_Signal=812.5 (corrected for dark counts)\n")
        @test HimalayaUI.parse_prp(mixed).transmitted_counts ≈ 812.5

        q = joinpath(dir, "no_signal.prp")
        write(q, "I_1(daqCounter1.0)=304300.0\nExposure time=15\n")
        @test HimalayaUI.parse_prp(q).transmitted_counts === missing
    end

    @testset "slot clustering" begin
        # Simulate one load: 3 slots × 4 frames each, H positions stepped ~4 mm apart,
        # within-burst jitter ≈ 0.3 mm. Mirrors the SSRL slot pattern from spec §5.
        t0 = DateTime(2026, 4, 26, 23, 14, 8)
        slot_h = [74.80, 70.85, 67.22]  # three slot centers
        jitter = [0.0, 0.13, -0.24, 0.05]  # within-burst jitter, 4 frames
        frames = HimalayaUI.ExposureMeta[]
        for (si, h_center) in enumerate(slot_h)
            for (fi, j) in enumerate(jitter)
                offset = (si - 1) * 4 * 19 + (fi - 1) * 19
                stem = "SY_$(si)_$(fi)_S$(lpad((si-1)*4+fi, 4,'0'))_0_001"
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

        # Single-frame position scan (real-data load-20 case): one frame per slot,
        # each ~3.5 mm from the last, no bursts. The median consecutive Δ IS the
        # slot step (~3.5 mm), so the old `median × slot_k` (~17 mm) tolerance
        # exceeded every gap and collapsed all 6 into ONE slot. The physical-floor
        # regime (median Δ > min_slot_separation_mm → split at the 1 mm floor)
        # recovers one slot per scan position.
        scan_frames = HimalayaUI.ExposureMeta[]
        for (i, pos) in enumerate([110.0, 106.5, 103.0, 99.5, 96.0, 92.5])
            push!(scan_frames, HimalayaUI.ExposureMeta("scan_$(i)_S$(lpad(i, 4, '0'))_0_001",
                nothing, nothing, nothing,
                (timestamp = t0 + Second((i - 1) * 19),
                 horizontal_position_mm = pos,
                 beam_energy_ev = 9000.0, energy_kev = 9.0, pipe_length_m = 1.7,
                 detector = "Pilatus 1M", exposure_time_s = 15.0)))
        end
        slots_scan = HimalayaUI._cluster_slots(scan_frames)
        @test length(slots_scan) == 6                     # one slot per scan position
        @test all(length(s) == 1 for s in slots_scan)
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
                    stem = "SY_$(li)_$(si)_$(fi)_S$(lpad(global_s, 4,'0'))_0_001"
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

        # Auto-naming: first sample of load 1, slot 1 → "SY ... (S01P01)"
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

        SLOT_CENTERS_MM = [74.80, 70.85, 67.22, 63.49]  # typical SSRL slot spacing
        WITHIN_BURST_JITTER = [0.0, 0.13, -0.24, 0.05]
        FRAME_GAP_S = 19
        t0 = DateTime(2026, 4, 25, 23, 14, 8)

        scan_id = 2001
        stems = String[]
        for (si, h_center) in enumerate(SLOT_CENTERS_MM)
            for (fi, jitter) in enumerate(WITHIN_BURST_JITTER)
                offset = (si - 1) * length(WITHIN_BURST_JITTER) * FRAME_GAP_S + (fi - 1) * FRAME_GAP_S
                stem = "SY_$(si)_$(lpad(scan_id, 4,'0'))_S$(scan_id)_0_001"
                scan_id += 1
                push!(stems, stem)
                write_prp(joinpath(data_dir, "$stem.prp");
                    timestamp        = Dates.format(t0 + Second(offset), "dd u yyyy HH:MM:SS"),
                    beam_energy_ev   = 9000.0,
                    pipe_length_mm   = 1700,
                    detector         = "Pilatus 1M",
                    exposure_time    = 15.0,
                    horizontal_position_mm = h_center + jitter)
                write(joinpath(data_dir, "$stem.tif"), "fake tif")
            end
        end
        write_setup_info(joinpath(analysis_dir, "setup_info_20260301_120000.txt");
            beam_center_x = 420.5, beam_center_y = 835.5, mean_distance_mm = 1810.0)

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
        setup_files = [joinpath(analysis_dir, "setup_info_20260301_120000.txt")]
        geo, disc   = HimalayaUI.derive_geometry(String.(prp_paths), String.(setup_files))

        @test geo.energy_kev ≈ 9.0
        @test geo.flight_path_m ≈ 1.81            # calibrated, not PRP nominal
        @test geo.flight_path_m_source == "setup"
        @test geo.beam_center_x ≈ 420.5
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
                stem = "SY_$(si)_$(fi)_S$(lpad((si-1)*2+fi, 4,'0'))_0_001"
                push!(all_stems, stem)
                write_prp(joinpath(data_dir, "$stem.prp");
                    timestamp = Dates.format(t0 + Second(offset), "dd u yyyy HH:MM:SS"),
                    beam_energy_ev = 9000.0,
                    pipe_length_mm = 1700, detector = "Pilatus 1M",
                    exposure_time = 15.0, horizontal_position_mm = h)
                write(joinpath(data_dir, "$stem.tif"), "fake tif")
                # No .dat files: analyze_exposure! will fail to read them, which is OK
                # for this insert-path test (we set analyze=false)
            end
        end
        write_setup_info(joinpath(analysis_dir, "setup_info_20260301_120000.txt"))

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
            "SELECT energy_kev, flight_path_m, flight_path_m_source, beam_center_x FROM setups WHERE experiment_id = ?",
            [exp_id])))
        @test row.energy_kev ≈ 9.0
        @test row.flight_path_m ≈ 1.81   # from setup file
        @test row.flight_path_m_source == "setup"
        @test row.beam_center_x ≈ 420.5

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

    @testset "scan_and_group! stores beam monitors; backfill_monitors! fills old rows once" begin
        dir = mktempdir()
        data_dir     = joinpath(dir, "data")
        analysis_dir = joinpath(dir, "analysis")
        mkpath(data_dir); mkpath(analysis_dir)
        t0 = DateTime(2026, 4, 26, 23, 14, 8)
        stems = ["SY_1_$(fi)_S$(lpad(fi, 4, '0'))_0_001" for fi in 1:2]
        for (fi, stem) in enumerate(stems)
            write_prp(joinpath(data_dir, "$stem.prp");
                timestamp = Dates.format(t0 + Second(19 * fi), "dd u yyyy HH:MM:SS"),
                exposure_time = 15.0, i_signal = 7000.0 + fi, i_0 = 1.5)
            write(joinpath(data_dir, "$stem.tif"), "fake tif")
        end
        write_setup_info(joinpath(analysis_dir, "setup_info_20260301_120000.txt"))
        db = fresh_db()
        exp_id = HimalayaUI.create_experiment!(db;
            name = "test-monitors", path = dir,
            data_dir = data_dir, analysis_dir = analysis_dir)
        HimalayaUI.scan_and_group!(db, exp_id; analyze = false)

        monitors() = Tables.rowtable(DBInterface.execute(db,
            "SELECT filename, transmitted_counts, incident_monitor, counting_time FROM exposures
              WHERE experiment_id = ? ORDER BY filename", [exp_id]))
        rows = monitors()
        @test [r.transmitted_counts for r in rows] ≈ [7001.0, 7002.0]
        @test all(r.incident_monitor ≈ 1.5 for r in rows)
        @test all(r.counting_time ≈ 15.8 for r in rows)

        @test all(r -> r.monitors_read == 1, Tables.rowtable(DBInterface.execute(db,
            "SELECT monitors_read FROM exposures WHERE experiment_id = ?", [exp_id])))

        DBInterface.execute(db, """UPDATE exposures SET transmitted_counts = NULL,
            incident_monitor = NULL, counting_time = NULL, monitors_read = 0""")
        DBInterface.execute(db,
            "UPDATE exposures SET transmitted_counts = 1.0 WHERE filename = ?", [stems[2]])
        HimalayaUI.scan_and_group!(db, exp_id; analyze = false)
        @test all(r -> ismissing(r.counting_time), monitors())

        @test HimalayaUI.backfill_monitors!(db, exp_id) == 2
        rows = monitors()
        @test rows[1].transmitted_counts ≈ 7001.0
        @test rows[2].transmitted_counts ≈ 1.0
        @test all(r.counting_time ≈ 15.8 for r in rows)
        @test HimalayaUI.backfill_monitors!(db, exp_id) == 0

        write_prp(joinpath(data_dir, "$(stems[1]).prp"); exposure_time = 15.0)
        DBInterface.execute(db, """UPDATE exposures SET transmitted_counts = NULL,
            incident_monitor = NULL, monitors_read = 0 WHERE filename = ?""", [stems[1]])
        @test HimalayaUI.backfill_monitors!(db, exp_id) == 1
        @test ismissing(monitors()[1].transmitted_counts)
        @test HimalayaUI.backfill_monitors!(db, exp_id) == 0
    end

    @testset "server boot backfills beam monitors on an upgraded experiment" begin
        dir = mktempdir()
        data_dir     = joinpath(dir, "data")
        analysis_dir = joinpath(dir, "analysis")
        mkpath(data_dir); mkpath(analysis_dir)
        t0 = DateTime(2026, 4, 26, 23, 14, 8)
        for fi in 1:2
            stem = "SY_1_$(fi)_S$(lpad(fi, 4, '0'))_0_001"
            write_prp(joinpath(data_dir, "$stem.prp");
                timestamp = Dates.format(t0 + Second(19 * fi), "dd u yyyy HH:MM:SS"),
                exposure_time = 15.0, i_signal = 7000.0 + fi, i_0 = 1.5)
            write(joinpath(data_dir, "$stem.tif"), "fake tif")
        end
        write_setup_info(joinpath(analysis_dir, "setup_info_20260301_120000.txt"))
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db;
            name = "upgraded", path = dir, data_dir = data_dir, analysis_dir = analysis_dir)
        HimalayaUI.scan_and_group!(db, exp_id; analyze = false)
        DBInterface.execute(db, """UPDATE exposures SET transmitted_counts = NULL,
            incident_monitor = NULL, counting_time = NULL, monitors_read = 0""")
        DBInterface.execute(db, """UPDATE experiments SET ingest_status = 'complete' WHERE id = ?""", [exp_id])
        try
            task = HimalayaUI.boot!(db)
            fetch(task)
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT monitors_read, transmitted_counts FROM exposures WHERE experiment_id = ?
                  ORDER BY filename", [exp_id]))
            @test all(r -> r.monitors_read == 1, rows)
            @test isequal([r.transmitted_counts for r in rows], [7001.0, 7002.0])
        finally
            HimalayaUI.stop_gc_timer!()
            HimalayaUI.stop_peak_training!()
            HimalayaUI.Oxygen.resetstate()
            HimalayaUI._DB_REF[] = nothing
            SQLite.close(db)
        end
    end

    @testset "migrate_exposures_beam_monitors! adds the columns to an old exposures table" begin
        db = SQLite.DB()
        DBInterface.execute(db, "CREATE TABLE schema_migrations (name TEXT PRIMARY KEY, applied_at TEXT)")
        DBInterface.execute(db, "CREATE TABLE exposures (id INTEGER PRIMARY KEY, filename TEXT)")
        DBInterface.execute(db, "INSERT INTO exposures (filename) VALUES ('a')")
        HimalayaUI.migrate_exposures_beam_monitors!(db)
        HimalayaUI.migrate_exposures_beam_monitors!(db)
        cols = HimalayaUI.cols_of(db, "exposures")
        @test all(c in cols for c in ("transmitted_counts", "incident_monitor", "counting_time", "monitors_read"))
        r = first(Tables.rowtable(DBInterface.execute(db, "SELECT * FROM exposures")))
        @test ismissing(r.transmitted_counts)
        @test r.monitors_read == 0
    end

    @testset "ingest_has_work sees new files" begin
        dir = mktempdir()
        data_dir     = joinpath(dir, "data")
        analysis_dir = joinpath(dir, "analysis")
        mkpath(data_dir); mkpath(analysis_dir)

        t0 = DateTime(2026, 4, 26, 23, 14, 8)
        # Write 2 TIF/PRP pairs
        for (i, h) in enumerate([58.9, 63.1])
            stem = "SY_$(i)_S$(lpad(i, 4,'0'))_0_001"
            write_prp(joinpath(data_dir, "$stem.prp");
                timestamp = Dates.format(t0 + Second((i-1)*19), "dd u yyyy HH:MM:SS"),
                beam_energy_ev = 9000.0, pipe_length_mm = 1700,
                detector = "Pilatus 1M", exposure_time = 15.0,
                horizontal_position_mm = h)
            write(joinpath(data_dir, "$stem.tif"), "fake tif")
        end
        write_setup_info(joinpath(analysis_dir, "setup_info_20260301_120000.txt"))

        db = fresh_db()
        exp_id = HimalayaUI.create_experiment!(db;
            name = "cheap-check", path = dir,
            data_dir = data_dir, analysis_dir = analysis_dir)

        # Before any ingest: 2 files on disk, 0 exposures in DB → changed
        @test HimalayaUI.ingest_has_work(db, exp_id) == true

        HimalayaUI.scan_and_group!(db, exp_id)
        @test HimalayaUI.ingest_has_work(db, exp_id) == false

        # Add a new TIF/PRP pair on disk (not yet ingested) → changed again
        new_stem = "SY_3_S0003_0_001"
        write_prp(joinpath(data_dir, "$new_stem.prp");
            timestamp = Dates.format(t0 + Second(600), "dd u yyyy HH:MM:SS"),
            beam_energy_ev = 9000.0, pipe_length_mm = 1700,
            detector = "Pilatus 1M", exposure_time = 15.0,
            horizontal_position_mm = 67.3)
        write(joinpath(data_dir, "$new_stem.tif"), "fake tif")
        @test HimalayaUI.ingest_has_work(db, exp_id) == true

        HimalayaUI.scan_and_group!(db, exp_id)
        @test HimalayaUI.ingest_has_work(db, exp_id) == false

        bad_db = fresh_db()
        bad_id = HimalayaUI.create_experiment!(bad_db;
            name = "missing-dir", path = joinpath(dir, "does_not_exist"),
            data_dir = joinpath(dir, "does_not_exist"),
            analysis_dir = analysis_dir)
        @test HimalayaUI.ingest_has_work(bad_db, bad_id) == false
    end

    @testset "ingest_has_work reads image_pattern from experiment row" begin
        dir = mktempdir()
        data_dir = joinpath(dir, "data")
        mkpath(data_dir)

        # Write 2 files with a custom suffix ".tiff" (not the default ".tif")
        write(joinpath(data_dir, "EX001.tiff"), "fake tif")
        write(joinpath(data_dir, "EX002.tiff"), "fake tif")

        db = fresh_db()
        exp_id = HimalayaUI.create_experiment!(db;
            name = "pattern-check", path = dir,
            data_dir = data_dir, analysis_dir = joinpath(dir, "analysis"),
            image_pattern = "{name}.tiff")

        # With the row's image_pattern="{name}.tiff": 2 files, 0 exposures → changed
        @test HimalayaUI.ingest_has_work(db, exp_id) == true

        db2 = fresh_db()
        exp_id2 = HimalayaUI.create_experiment!(db2;
            name = "pattern-check-tif", path = dir,
            data_dir = data_dir, analysis_dir = joinpath(dir, "analysis"))
        # image_pattern is NULL → default ".tif" → 0 .tif files → false (nothing to ingest)
        @test HimalayaUI.ingest_has_work(db2, exp_id2) == false
    end

    @testset "_assign_sessions" begin
        t0 = DateTime(2026, 4, 25, 10, 0, 0)

        # Two loads with a >3 h gap → two distinct session indices.
        starts_two = [t0, t0 + Dates.Hour(4)]
        ids_two = HimalayaUI._assign_sessions(starts_two)
        @test ids_two[1] == 1
        @test ids_two[2] == 2
        @test length(unique(ids_two)) == 2

        # All loads close together → all session index 1.
        starts_close = [t0, t0 + Dates.Hour(1), t0 + Dates.Hour(2)]
        ids_close = HimalayaUI._assign_sessions(starts_close)
        @test all(==(1), ids_close)

        # A missing start_time in the middle stays in the current session and
        # does not advance prev: [real, missing, real_close] → all session 1.
        starts_miss = Union{DateTime, Missing}[t0, missing, t0 + Dates.Hour(2)]
        ids_miss = HimalayaUI._assign_sessions(starts_miss)
        @test ids_miss[1] == 1
        @test ids_miss[2] == 1   # missing stays in current session
        @test ids_miss[3] == 1   # 2 h gap (below 3 h) → same session

        # A missing that sits between two real times with a >3 h outer gap: the
        # session only advances when we see the next real timestamp.
        starts_miss2 = Union{DateTime, Missing}[t0, missing, t0 + Dates.Hour(5)]
        ids_miss2 = HimalayaUI._assign_sessions(starts_miss2)
        @test ids_miss2[1] == 1
        @test ids_miss2[2] == 1   # missing: stays in session 1
        @test ids_miss2[3] == 2   # 5 h gap from prev (t0) → new session
    end

    @testset "scan_and_group! persists non-null session_id on loads" begin
        dir = mktempdir()
        data_dir     = joinpath(dir, "data")
        analysis_dir = joinpath(dir, "analysis")
        mkpath(data_dir); mkpath(analysis_dir)

        # One load (all frames close together), so we expect session_id = 1.
        t0 = DateTime(2026, 4, 26, 23, 14, 8)
        for (i, h) in enumerate([58.9, 63.1])
            stem = "SY_$(i)_S$(lpad(i, 4,'0'))_0_001"
            write_prp(joinpath(data_dir, "$stem.prp");
                timestamp = Dates.format(t0 + Second((i-1)*19), "dd u yyyy HH:MM:SS"),
                beam_energy_ev = 9000.0, pipe_length_mm = 1700,
                detector = "Pilatus 1M", exposure_time = 15.0,
                horizontal_position_mm = h)
            write(joinpath(data_dir, "$stem.tif"), "fake tif")
        end
        write_setup_info(joinpath(analysis_dir, "setup_info_20260301_120000.txt"))

        db = fresh_db()
        exp_id = HimalayaUI.create_experiment!(db;
            name = "session-test", path = dir,
            data_dir = data_dir, analysis_dir = analysis_dir)
        HimalayaUI.scan_and_group!(db, exp_id; analyze = false)

        loads = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, session_id FROM loads WHERE experiment_id = ?", [exp_id]))
        @test length(loads) >= 1
        @test all(!ismissing(l.session_id) && l.session_id !== nothing for l in loads)
        # All frames in one load → one session (session_id = 1).
        @test all(Int(l.session_id) == 1 for l in loads)
    end

    # -----------------------------------------------------------------------
    # on_progress is (processed, total, stage). Each stage is a SEGMENT of the
    # UI's progress bar and reports against its OWN denominator, so the bar
    # advances left-to-right through segments instead of resetting one shared
    # scale. Helper: collect the ticks for one stage.
    # -----------------------------------------------------------------------
    _stage(ticks, s) = [(p, t) for (p, t, st) in ticks if st == s]

    function _progress_fixture(stems)
        db = fresh_db()
        data_dir = mktempdir()
        for stem in stems
            touch(joinpath(data_dir, "$stem.tif"))
            write_prp(joinpath(data_dir, "$stem.prp"))      # default kwargs OK; see :20
        end
        exp_id = HimalayaUI.create_experiment!(db; name="t", path=data_dir,
                                              data_dir=data_dir, analysis_dir=data_dir)
        (db, exp_id)
    end

    @testset "scan_and_group! ingests a TIF-only directory (no PRPs at all)" begin
        # Regression for the parsed_prps element type. ExposureMeta.prp is
        # Union{NamedTuple,Nothing}, so with NO PRPs anywhere the filtered
        # comprehension collects to Vector{Any} and derive_geometry MethodErrors --
        # turning "geometry falls back to defaults" into "the ingest crashes". This
        # is the standard shape of a MISCONFIGURED prp_pattern, i.e. the exact
        # situation where a crash is least helpful.
        db = fresh_db()
        data_dir = mktempdir()
        for stem in ("e1", "e2")
            touch(joinpath(data_dir, "$stem.tif"))      # deliberately NO .prp
        end
        exp_id = HimalayaUI.create_experiment!(db; name="tifonly", path=data_dir,
                                              data_dir=data_dir, analysis_dir=data_dir)
        res = HimalayaUI.scan_and_group!(db, exp_id; analyze=false)
        @test res.status == :ok
        @test res.added_exposures == 2
        # Geometry stayed at defaults rather than throwing.
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT energy_kev_source, pixel_size_um_source FROM setups WHERE experiment_id = ?", [exp_id])))
        @test row.energy_kev_source == "default"
        @test row.pixel_size_um_source == "default"
    end

    @testset "scan_and_group! on_progress reports each stage on its own denominator" begin
        db, exp_id = _progress_fixture(("e1", "e2", "e3"))
        ticks = Tuple{Int,Int,String}[]
        HimalayaUI.scan_and_group!(db, exp_id; analyze=true,
                                   on_progress=(p, t, s) -> push!(ticks, (p, t, s)))

        # Every tick must name a stage — the frontend keys its segment off this.
        @test all(t -> !isempty(t[3]), ticks)

        disc = _stage(ticks, "discovery")
        anal = _stage(ticks, "analyzing")
        @test !isempty(disc)
        @test !isempty(anal)
        # Each segment fills 0→full against its own total, monotonically.
        @test last(disc) == (3, 3)
        @test last(anal) == (3, 3)
        @test issorted(first.(disc))
        @test issorted(first.(anal))
        # Discovery reports BEFORE analysis begins: the whole point is that the
        # slow silent half now has a moving bar.
        @test findfirst(t -> t[3] == "discovery", ticks) <
              findfirst(t -> t[3] == "analyzing",  ticks)
    end

    @testset "scan_and_group! on_progress reports discovery even with analyze=false" begin
        db, exp_id = _progress_fixture(("e1", "e2", "e3"))
        ticks = Tuple{Int,Int,String}[]
        # analyze=false skips the analyze + thumbnail stages entirely, leaving only
        # discovery. Regression: discovery used to emit NOTHING — the total was
        # first published inside the analyze loop, so the UI had no scale at all
        # during the slow scan_directory + geometry + persist half and sat on
        # "0 / ~0" throughout.
        HimalayaUI.scan_and_group!(db, exp_id; analyze=false,
                                   on_progress=(p, t, s) -> push!(ticks, (p, t, s)))
        @test !isempty(ticks)
        @test all(t -> t[3] == "discovery", ticks)
        @test last(ticks) == (3, 3, "discovery")
        # Progress arrives DURING the walk, not just as a closing summary.
        @test first(ticks) == (1, 3, "discovery")
    end

    @testset "scan_and_group! on_progress closes the analyzing segment on a clean rescan" begin
        db, exp_id = _progress_fixture(("e1", "e2", "e3"))
        HimalayaUI.scan_and_group!(db, exp_id; analyze=true)   # first scan ingests everything

        # Second scan: insert-only dedup means ZERO new exposures, so the analyze
        # loop body never runs. It used to report `i / length(new_exposure_ids)` and
        # therefore emitted nothing at all, leaving the bar stalled until
        # ingest_complete wiped it. Now the segment is explicitly closed as 0-of-0,
        # which the frontend renders as a complete segment rather than an empty one.
        ticks = Tuple{Int,Int,String}[]
        HimalayaUI.scan_and_group!(db, exp_id; analyze=true,
                                   on_progress=(p, t, s) -> push!(ticks, (p, t, s)))
        anal = _stage(ticks, "analyzing")
        @test !isempty(anal)
        @test last(anal) == (0, 0)             # nothing to do, segment closed
        # Discovery still walked all three files (they exist on disk regardless).
        @test last(_stage(ticks, "discovery")) == (3, 3)
    end
end
