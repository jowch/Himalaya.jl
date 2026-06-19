# packages/HimalayaUI/src/geometry.jl
#
# Geometry derivation (spec §6):
#   1. detector_pixel_size_um(model) — static lookup table.
#   2. parse_setup_info(path) — parse one analysis/setup_info_*.txt file.
#   3. derive_geometry(prp_paths, setup_files) — orchestrate (Task 3).
#
# Setup-file format confirmed against real SSRL data (2026-06-18):
#   analysis/setup_info_20260425_181705.txt (11 lines):
#     "Beam center is at (421.409, 836.946)"
#     "Mean distance:         1809.5 mm"
#
# Geometry authority (spec §6, confirmed Jonathan 2026-06-18):
#   flight_path_m ← setup file Mean distance (AgBe calibration) when present,
#                   else PRP Pipe length.  Source tagged accordingly.

# ---------------------------------------------------------------------------
# Detector → pixel pitch lookup
# ---------------------------------------------------------------------------

"""Detector model string → pixel pitch in µm, or `missing` for unknown models."""
const _DETECTOR_PITCH_UM = Dict{String, Float64}(
    # Pilatus family (Dectris): all 172 µm pixels
    "Pilatus 1M"   => 172.0,
    "Pilatus 2M"   => 172.0,
    "Pilatus 6M"   => 172.0,
    "Pilatus 300K" => 172.0,
    "Pilatus 200K" => 172.0,
    # Eiger family (Dectris): all 75 µm pixels
    "Eiger 1M"     => 75.0,
    "Eiger 4M"     => 75.0,
    "Eiger 9M"     => 75.0,
    "Eiger 16M"    => 75.0,
)

"""
    detector_pixel_size_um(model) -> Union{Float64, Missing}

Look up the pixel pitch (µm) for a detector model string.
Returns `missing` for unknown models; the caller should flag for manual entry
rather than guessing (spec §6).
"""
function detector_pixel_size_um(model::AbstractString)::Union{Float64, Missing}
    # Exact match first
    haskey(_DETECTOR_PITCH_UM, model) && return _DETECTOR_PITCH_UM[model]
    # Prefix-based match for strings like "Pilatus 1M (dezingered)"
    for (k, v) in _DETECTOR_PITCH_UM
        startswith(model, k) && return v
    end
    return missing
end

# ---------------------------------------------------------------------------
# setup_info_*.txt parser
# ---------------------------------------------------------------------------

"""
    parse_setup_info(path) -> NamedTuple

Parse one `analysis/setup_info_<YYYYMMDD_HHMMSS>.txt` file and extract the
beam center and calibrated detector distance. Returns `missing` for absent fields.

Confirmed format (SSRL 2026-04, setup_info_20260425_181705.txt):
  "Beam center is at (421.409, 836.946)"
  "Mean distance:         1809.5 mm"
"""
function parse_setup_info(path::AbstractString)
    isfile(path) || error("setup_info file not found: $path")

    beam_center_x    = missing
    beam_center_y    = missing
    mean_distance_m  = missing

    for line in eachline(path)
        # "Beam center is at (X, Y)"
        m_bc = match(r"Beam center is at \(\s*([0-9.]+)\s*,\s*([0-9.]+)\s*\)", line)
        if m_bc !== nothing
            try
                beam_center_x = parse(Float64, m_bc.captures[1])
                beam_center_y = parse(Float64, m_bc.captures[2])
            catch
                @warn "parse_setup_info: could not parse beam center" path line
            end
            continue
        end

        # "Mean distance:         1809.5 mm"
        m_md = match(r"Mean distance:\s*([0-9.]+)\s*mm", line)
        if m_md !== nothing
            try
                mean_distance_m = parse(Float64, m_md.captures[1]) / 1000.0
            catch
                @warn "parse_setup_info: could not parse Mean distance" path line
            end
            continue
        end
    end

    return (
        beam_center_x   = beam_center_x,
        beam_center_y   = beam_center_y,
        mean_distance_m = mean_distance_m,
    )
end
