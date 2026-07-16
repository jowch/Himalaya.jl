# packages/HimalayaUI/src/prp.jl
#
# PRP file parser (§9.1 shared ingest core).
# Format: "Key=Value" lines with optional trailing unit tokens (mm, eV, etc.).
# Confirmed against the real SSRL 2026-04 PRP: HA_85_422_S2404_0_001.prp.
#
# Returned NamedTuple fields:
#   timestamp              :: Union{DateTime, Missing}   — "Time this file was written:"
#   beam_energy_ev         :: Union{Float64, Missing}    — raw eV value
#   energy_kev             :: Union{Float64, Missing}    — beam_energy_ev / 1000
#   pipe_length_m          :: Union{Float64, Missing}    — "Pipe length" converted m
#   detector               :: Union{String, Missing}     — e.g. "Pilatus 1M"
#   exposure_time_s        :: Union{Float64, Missing}    — "Exposure time"
#   horizontal_position_mm :: Union{Float64, Missing}    — "Horizontal position"
#
# Dropped fields (constant/noise per spec §4): Vertical position, dispx, dispy,
# detector_vert/horz, Phi, slit widths, attenuator, per-counter voltages, scan motor.

using Dates

"""
    parse_prp(path) -> NamedTuple

Parse one SSRL PRP file at `path` and return the fields the ingest system uses.
Unparseable or absent fields are `missing`. Never throws on malformed lines;
a `@warn` is emitted instead so batch ingestion can proceed.
"""
function parse_prp(path::AbstractString)
    isfile(path) || error("PRP file not found: $path")

    timestamp              = missing
    beam_energy_ev         = missing
    energy_kev             = missing
    pipe_length_m          = missing
    detector               = missing
    exposure_time_s        = missing
    horizontal_position_mm = missing

    for line in eachline(path)
        # Timestamp line has a different format: "Time this file was written: DD Mon YYYY HH:MM:SS"
        if startswith(line, "Time this file was written:")
            raw = strip(line[length("Time this file was written:") + 1:end])
            try
                timestamp = DateTime(raw, dateformat"dd u yyyy HH:MM:SS")
            catch
                @warn "parse_prp: could not parse timestamp" path line
            end
            continue
        end

        idx = findfirst('=', line)
        idx === nothing && continue
        key = strip(line[1:idx-1])
        val_raw = strip(line[idx+1:end])

        # Strip trailing unit token (first non-numeric, non-dot, non-sign char onward).
        # E.g. "9000.03 eV" → 9000.03; "1700 mm" → 1700; "15" → 15.
        val_str = let s = val_raw
            m = match(r"^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?", s)
            m === nothing ? nothing : m.match
        end

        if key == "Beam energy" && val_str !== nothing
            try
                beam_energy_ev = parse(Float64, val_str)
                energy_kev     = beam_energy_ev / 1000.0
            catch
                @warn "parse_prp: could not parse Beam energy" path val_raw
            end

        elseif key == "Pipe length" && val_str !== nothing
            try
                pipe_length_m = parse(Float64, val_str) / 1000.0
            catch
                @warn "parse_prp: could not parse Pipe length" path val_raw
            end

        elseif key == "Detector"
            # e.g. "Pilatus 1M" — string, no numeric conversion needed
            detector = val_raw

        elseif key == "Exposure time" && val_str !== nothing
            try
                exposure_time_s = parse(Float64, val_str)
            catch
                @warn "parse_prp: could not parse Exposure time" path val_raw
            end

        elseif key == "Horizontal position" && val_str !== nothing
            try
                horizontal_position_mm = parse(Float64, val_str)
            catch
                @warn "parse_prp: could not parse Horizontal position" path val_raw
            end
        end
    end

    return (
        timestamp              = timestamp,
        beam_energy_ev         = beam_energy_ev,
        energy_kev             = energy_kev,
        pipe_length_m          = pipe_length_m,
        detector               = detector,
        exposure_time_s        = exposure_time_s,
        horizontal_position_mm = horizontal_position_mm,
    )
end
