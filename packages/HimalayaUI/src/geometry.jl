const _DETECTOR_PITCH_UM = Dict{String, Float64}(
    "Pilatus 1M"   => 172.0,
    "Pilatus 2M"   => 172.0,
    "Pilatus 6M"   => 172.0,
    "Pilatus 300K" => 172.0,
    "Pilatus 200K" => 172.0,
    "Eiger 1M"     => 75.0,
    "Eiger 4M"     => 75.0,
    "Eiger 9M"     => 75.0,
    "Eiger 16M"    => 75.0,
)

function detector_pixel_size_um(model::AbstractString)::Union{Float64, Missing}
    haskey(_DETECTOR_PITCH_UM, model) && return _DETECTOR_PITCH_UM[model]
    for (k, v) in _DETECTOR_PITCH_UM
        startswith(model, k) && return v
    end
    return missing
end

function parse_setup_info(path::AbstractString)
    isfile(path) || error("setup_info file not found: $path")

    beam_center_x    = missing
    beam_center_y    = missing
    mean_distance_m  = missing

    for line in eachline(path)
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

struct GeometryDiscrepancy
    field   ::String
    message ::String
end

function derive_geometry(
    prp_paths::AbstractVector{<:AbstractString},
    setup_files::AbstractVector{<:AbstractString},
)
    derive_geometry(NamedTuple[parse_prp(p) for p in prp_paths], setup_files)
end

function derive_geometry(
    parsed::AbstractVector{<:NamedTuple},
    setup_files::AbstractVector{<:AbstractString},
)
    discrepancies = GeometryDiscrepancy[]

    energy_vals = filter(!ismissing, [p.energy_kev for p in parsed])
    energy_kev = isempty(energy_vals) ? missing : first(energy_vals)
    if length(unique(round.(skipmissing([p.beam_energy_ev for p in parsed]); digits=2))) > 1
        push!(discrepancies, GeometryDiscrepancy("beam_energy_ev",
            "beam energy varies across PRPs: " *
            join(unique(round.(skipmissing([p.beam_energy_ev for p in parsed]); digits=2)), ", ") * " eV"))
    end

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

    beam_center_x        = missing; beam_center_x_source = "default"
    beam_center_y        = missing; beam_center_y_source = "default"
    flight_path_m        = missing; flight_path_m_source = "default"

    if !isempty(setup_files)
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

    # physics: the AgBe-calibrated setup distance beats the PRP's nominal pipe
    # physics: length. They differ by several percent on real runs, which is the point of
    # physics: calibrating, so the gap is not reported as a discrepancy.
    if ismissing(flight_path_m)
        pipe_vals = filter(!ismissing, [p.pipe_length_m for p in parsed])
        if !isempty(pipe_vals)
            flight_path_m = first(pipe_vals); flight_path_m_source = "prp"
            if length(unique(round.(pipe_vals; digits=4))) > 1
                push!(discrepancies, GeometryDiscrepancy("flight_path_m",
                    "PRP pipe lengths vary: " * join(unique(pipe_vals .* 1000), ", ") * " mm"))
            end
        end
    end

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
