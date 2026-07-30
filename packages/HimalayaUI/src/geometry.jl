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

# ---------------------------------------------------------------------------
# Geometry derivation orchestrator
# ---------------------------------------------------------------------------

"""
    GeometryDiscrepancy

A detected variance in a constant PRP field across multiple exposures,
or an unresolvable field (unknown detector). Surfaced as a Configuration-tab
banner (spec §6).
"""
struct GeometryDiscrepancy
    field   ::String
    message ::String
end

"""
    derive_geometry(prp_paths, setup_files) -> (geometry, discrepancies)

Sample the PRP files at `prp_paths` (reads all of them; at SSRL scale this is
~700 files × ~50 lines each, sub-second), derive per-field geometry with source
tags, and collect discrepancies.

`setup_files` is a vector of `analysis/setup_info_*.txt` paths; the latest by
filename sort is used for beam center + calibrated distance (filenames encode a
`YYYYMMDD_HHMMSS` timestamp so lexicographic sort is chronological).

Returns:
  `geometry` — a NamedTuple with fields:
      energy_kev, energy_kev_source,
      flight_path_m, flight_path_m_source,
      beam_center_x, beam_center_x_source,
      beam_center_y, beam_center_y_source,
      pixel_size_um, pixel_size_um_source
  `discrepancies` — a Vector{GeometryDiscrepancy}
"""
function derive_geometry(
    prp_paths::AbstractVector{<:AbstractString},
    setup_files::AbstractVector{<:AbstractString},
)
    # Parse all PRPs (small files, sequential read is fine), then delegate.
    # The NamedTuple[...] element type is REQUIRED, not decoration: parse_prp's
    # inferred return is a non-concrete NamedTuple{names,<:Tuple{...}}, so an
    # EMPTY prp_paths collects to Vector{Any} and matches neither method ->
    # MethodError, where this used to return an all-"default" geometry. Reachable
    # from routes_fs.jl (a directory with TIFs and no PRPs) and regroup_experiment!.
    derive_geometry(NamedTuple[parse_prp(p) for p in prp_paths], setup_files)
end

"""
    derive_geometry(parsed_prps, setup_files) -> (geometry, discrepancies)

Pre-parsed variant: takes `parse_prp` NamedTuples instead of paths, so a caller
that has already read the PRPs does not read them a second time.

`scan_directory` parses every PRP into `ExposureMeta.prp`, and the path-based
method above then re-read and re-parsed all of them — two full passes over every
PRP in the experiment per ingest. Over SMB each parse is a network round trip, so
on an ~1100-exposure experiment that was ~1100 avoidable file reads. The ingest
paths (`scan_and_group!`, `regroup_experiment!`) pass `m.prp` straight through;
the funnel preview (`routes_fs.jl`) has only paths and still uses the method
above.
"""
function derive_geometry(
    parsed::AbstractVector{<:NamedTuple},
    setup_files::AbstractVector{<:AbstractString},
)
    discrepancies = GeometryDiscrepancy[]

    # --- Energy (from PRP, should be constant) ---
    energy_vals = filter(!ismissing, [p.energy_kev for p in parsed])
    energy_kev = isempty(energy_vals) ? missing : first(energy_vals)
    if length(unique(round.(skipmissing([p.beam_energy_ev for p in parsed]); digits=2))) > 1
        push!(discrepancies, GeometryDiscrepancy("beam_energy_ev",
            "beam energy varies across PRPs: " *
            join(unique(round.(skipmissing([p.beam_energy_ev for p in parsed]); digits=2)), ", ") * " eV"))
    end

    # --- Pixel size from detector model (constant across run) ---
    detectors = unique(filter(!ismissing, [p.detector for p in parsed]))
    pixel_size_um    = missing
    pixel_size_source = "default"
    if length(detectors) == 1
        pitch = detector_pixel_size_um(detectors[1])
        if pitch === missing
            push!(discrepancies, GeometryDiscrepancy("pixel_size_um",
                "unknown detector model '$(detectors[1])': pixel pitch unknown, manual entry required"))
        else
            pixel_size_um    = pitch
            pixel_size_source = "prp"
        end
    elseif length(detectors) > 1
        push!(discrepancies, GeometryDiscrepancy("pixel_size_um",
            "multiple detector models found: $(join(detectors, ", "))"))
    end

    # --- Beam center + flight path from setup file ---
    beam_center_x        = missing; beam_center_x_source = "default"
    beam_center_y        = missing; beam_center_y_source = "default"
    flight_path_m        = missing; flight_path_m_source = "default"

    if !isempty(setup_files)
        # Filenames are "setup_info_YYYYMMDD_HHMMSS.txt"; lexicographic sort picks the latest.
        latest_setup = last(sort(setup_files))
        si = parse_setup_info(latest_setup)
        if !ismissing(si.beam_center_x)
            beam_center_x = si.beam_center_x; beam_center_x_source = "setup"
        end
        if !ismissing(si.beam_center_y)
            beam_center_y = si.beam_center_y; beam_center_y_source = "setup"
        end
        if !ismissing(si.mean_distance_m)
            flight_path_m = si.mean_distance_m; flight_path_m_source = "setup"
        end
    end

    # Fallback: use PRP Pipe length when setup file absent or mean distance unparseable
    if ismissing(flight_path_m)
        pipe_vals = filter(!ismissing, [p.pipe_length_m for p in parsed])
        if !isempty(pipe_vals)
            flight_path_m = first(pipe_vals); flight_path_m_source = "prp"
            # Flag if pipe lengths vary (multi-setup discrepancy)
            if length(unique(round.(pipe_vals; digits=4))) > 1
                push!(discrepancies, GeometryDiscrepancy("flight_path_m",
                    "PRP pipe lengths vary: " * join(unique(pipe_vals .* 1000), ", ") * " mm"))
            end
        end
    end

    # Note: a gap between the PRP nominal pipe length and the setup-file calibrated
    # Mean distance (≈6.4% on the real SSRL 2026-04 data) is EXPECTED — AgBe ring
    # calibration is precisely meant to refine the nominal mechanical value. We use
    # the calibrated value for flight_path_m (above) but do NOT emit a discrepancy
    # for the gap, since it would fire on essentially every real run (option A,
    # ratified 2026-06-19; design doc §6).

    geometry = (
        energy_kev             = energy_kev,
        energy_kev_source      = ismissing(energy_kev) ? "default" : "prp",
        flight_path_m          = flight_path_m,
        flight_path_m_source   = flight_path_m_source,
        beam_center_x          = beam_center_x,
        beam_center_x_source   = beam_center_x_source,
        beam_center_y          = beam_center_y,
        beam_center_y_source   = beam_center_y_source,
        pixel_size_um          = pixel_size_um,
        pixel_size_um_source   = pixel_size_source,
    )
    return (geometry, discrepancies)
end
