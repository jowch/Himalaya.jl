using Dates

function parse_prp(path::AbstractString)
    isfile(path) || error("PRP file not found: $path")

    timestamp              = missing
    beam_energy_ev         = missing
    energy_kev             = missing
    pipe_length_m          = missing
    detector               = missing
    exposure_time_s        = missing
    horizontal_position_mm = missing
    transmitted_counts     = missing
    incident_monitor       = missing
    counting_time_s        = missing

    for line in eachline(path)
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
            detector = val_raw

        elseif key == "Exposure time" && val_str !== nothing
            try
                exposure_time_s = parse(Float64, val_str)
            catch
                @warn "parse_prp: could not parse Exposure time" path val_raw
            end

        # physics: I_2 is the photodiode in the beamstop, the beam that got through the
        # physics: sample: incident flux times transmission. The beamline writes it with
        # physics: dark counts already removed, and they are a small part of the reading.
        # physics: I_signal is not used: it copies whichever counter the file's "is from:"
        # physics: line names, usually the APD, which reads close to the incident beam
        # physics: and is mostly dark counts.
        elseif startswith(key, "I_2(") && val_str !== nothing
            try
                transmitted_counts = parse(Float64, val_str)
            catch
                @warn "parse_prp: could not parse I_2" path val_raw
            end

        # physics: I_0 is the ion chamber upstream of the sample, the incident flux.
        elseif startswith(key, "I_0(") && val_str !== nothing
            try
                incident_monitor = parse(Float64, val_str)
            catch
                @warn "parse_prp: could not parse I_0" path val_raw
            end

        elseif key == "Counting time" && val_str !== nothing
            try
                counting_time_s = parse(Float64, val_str)
            catch
                @warn "parse_prp: could not parse Counting time" path val_raw
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
        transmitted_counts     = transmitted_counts,
        incident_monitor       = incident_monitor,
        counting_time_s        = counting_time_s,
    )
end
