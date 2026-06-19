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
    # task @testsets are appended below
end
