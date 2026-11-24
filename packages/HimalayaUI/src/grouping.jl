struct ExposureMeta
    stem         ::String
    tif_path     ::Union{String, Nothing}
    prp_path     ::Union{String, Nothing}
    dat_path     ::Union{String, Nothing}
    prp          ::Union{NamedTuple, Nothing}
end

function tif_listing(data_dir::AbstractString, tif_pattern::AbstractString)
    tif_parts = split(tif_pattern, "{name}"; limit=2)
    length(tif_parts) == 2 || error("tif_pattern must contain {{name}}: $tif_pattern")
    tif_prefix_literal = String(tif_parts[1])
    tif_suffix = String(tif_parts[2])
    scan_subdir = dirname(tif_prefix_literal)
    tif_scan_dir = isempty(scan_subdir) ? String(data_dir) : joinpath(data_dir, scan_subdir)
    isdir(tif_scan_dir) || return nothing
    tif_files = sort!(filter(f -> startswith(f, tif_prefix_literal) && endswith(f, tif_suffix),
                             list_names(tif_scan_dir)))
    (tif_scan_dir, tif_files, tif_suffix)
end

function scan_directory(
    data_dir::AbstractString,
    analysis_dir::AbstractString;
    tif_pattern::String = "{name}.tif",
    prp_pattern::String = "{name}.prp",
    dat_pattern::String = "{name}.dat",
    on_progress::Union{Function, Nothing} = nothing,
)::Vector{ExposureMeta}
    sidecar(base, stem, pattern) =
        (p = joinpath(base, replace(pattern, "{name}" => stem)); isfile(p) ? p : nothing)

    listing = tif_listing(data_dir, tif_pattern)
    listing === nothing && return ExposureMeta[]
    tif_scan_dir, tif_files, tif_suffix = listing

    suffix_len = length(tif_suffix)
    metas = ExposureMeta[]
    n_total = length(tif_files)
    step    = max(1, n_total ÷ 100)
    for (i, fname) in enumerate(tif_files)
        stem = fname[1:end - suffix_len]
        tif_path  = joinpath(tif_scan_dir, fname)
        prp_path  = sidecar(data_dir, stem, prp_pattern)
        dat_path  = sidecar(analysis_dir, stem, dat_pattern)
        prp_parsed = prp_path !== nothing ? parse_prp(prp_path) : nothing
        push!(metas, ExposureMeta(stem, tif_path, prp_path, dat_path, prp_parsed))
        if on_progress !== nothing && (i % step == 0 || i == n_total)
            on_progress(i, n_total, "discovery")
        end
    end
    return metas
end

function _segment_loads_with_flag(metas::Vector{ExposureMeta}; gap_k::Float64 = 10.0)
    isempty(metas) && return (Vector{ExposureMeta}[], :ok)

    sorted = sort(metas; by = m -> begin
        ts = m.prp !== nothing ? m.prp.timestamp : missing
        ismissing(ts) ? DateTime(9999) : ts
    end)

    length(sorted) == 1 && return ([[sorted[1]]], :ok)

    timestamps = [begin
        ts = m.prp !== nothing ? m.prp.timestamp : missing
        ismissing(ts) ? missing : ts
    end for m in sorted]

    gaps = Float64[]
    for i in 2:length(timestamps)
        if !ismissing(timestamps[i]) && !ismissing(timestamps[i-1])
            push!(gaps, Float64(Dates.value(timestamps[i] - timestamps[i-1])) / 1000.0)
        end
    end

    if isempty(gaps)
        return ([sorted], :unimodal_fallback)
    end

    med = _median_inline(gaps)
    threshold = med * gap_k
    flag = :unimodal_fallback

    loads = Vector{ExposureMeta}[[sorted[1]]]
    for i in 2:length(sorted)
        ts_prev = sorted[i-1].prp !== nothing ? sorted[i-1].prp.timestamp : missing
        ts_curr = sorted[i].prp   !== nothing ? sorted[i].prp.timestamp   : missing
        gap_s = (!ismissing(ts_prev) && !ismissing(ts_curr)) ?
            Float64(Dates.value(ts_curr - ts_prev)) / 1000.0 : 0.0
        if gap_s > threshold
            push!(loads, ExposureMeta[])
            flag = :ok
        end
        push!(last(loads), sorted[i])
    end

    return (loads, flag)
end

function _median_inline(xs::Vector{Float64})::Float64
    v = sort(xs)
    n = length(v)
    iseven(n) ? (v[n÷2] + v[n÷2+1]) / 2.0 : v[(n+1)÷2]
end

function _cluster_slots(
    load_metas::Vector{ExposureMeta};
    slot_k::Float64 = 5.0,
    min_slot_separation_mm::Float64 = 1.0,
)::Vector{Vector{ExposureMeta}}
    isempty(load_metas) && return Vector{ExposureMeta}[]
    length(load_metas) == 1 && return [[load_metas[1]]]

    hpos = Union{Float64, Missing}[
        (m.prp !== nothing && !ismissing(m.prp.horizontal_position_mm)) ?
            m.prp.horizontal_position_mm : missing
        for m in load_metas
    ]

    if all(ismissing, hpos)
        return [load_metas]
    end

    deltas = Float64[]
    for i in 2:length(hpos)
        if !ismissing(hpos[i]) && !ismissing(hpos[i-1])
            push!(deltas, abs(hpos[i] - hpos[i-1]))
        end
    end

    if isempty(deltas)
        return [load_metas]
    end

    med_delta = _median_inline(deltas)
    local threshold::Float64
    # physics: capillary slots on the holder sit at least ~1 mm apart, so a larger
    # physics: step between consecutive frames is always a new slot. Below that,
    # physics: within-slot jitter dominates and the adaptive tolerance applies.
    if med_delta < 1e-6
        nonzero = filter(d -> d > 1e-6, deltas)
        threshold = isempty(nonzero) ? Inf : _median_inline(nonzero) / slot_k
    elseif med_delta > min_slot_separation_mm
        threshold = min_slot_separation_mm
    else
        threshold = med_delta * slot_k
    end

    slots = Vector{ExposureMeta}[[load_metas[1]]]
    for i in 2:length(load_metas)
        h_prev = hpos[i-1]
        h_curr = hpos[i]
        gap = (!ismissing(h_prev) && !ismissing(h_curr)) ? abs(h_curr - h_prev) : 0.0
        if gap > threshold
            push!(slots, ExposureMeta[])
        end
        push!(last(slots), load_metas[i])
    end
    return slots
end

struct GroupedExposure
    stem       ::String
    tif_path   ::Union{String, Nothing}
    prp_path   ::Union{String, Nothing}
    dat_path   ::Union{String, Nothing}
    timestamp  ::Union{DateTime, Missing}
    exposure_time_s        ::Union{Float64, Missing}
    horizontal_position_mm ::Union{Float64, Missing}
    transmitted_counts     ::Union{Float64, Missing}
    incident_monitor       ::Union{Float64, Missing}
    counting_time_s        ::Union{Float64, Missing}
    pipe_length_m          ::Union{Float64, Missing}
    scan_id    ::Union{Int, Missing}
    frame_no   ::Union{Int, Missing}
end

struct GroupedSample
    name           ::String
    name_source    ::String
    slot_index     ::Int
    grouping_source::String
    exposures      ::Vector{GroupedExposure}
end

struct GroupedLoad
    load_index    ::Int
    frame_count   ::Int
    start_time    ::Union{DateTime, Missing}
    end_time      ::Union{DateTime, Missing}
    samples       ::Vector{GroupedSample}
    flag          ::Symbol
    session_index ::Int
end

struct GroupingResult
    loads        ::Vector{GroupedLoad}
    discrepancies::Vector{String}
end

function _auto_name(label_hint::AbstractString, load_index::Int, slot_index::Int)
    coord = "S$(lpad(load_index, 2, '0'))P$(lpad(slot_index, 2, '0'))"
    isempty(strip(label_hint)) && return "($(coord))"
    return "$(strip(label_hint)) ($(coord))"
end

function _label_from_stem(stem::AbstractString)
    m = match(r"^(.+?)_S\d+_", stem)
    m === nothing && return String(stem)
    return String(m.captures[1])
end

function _parse_scan_frame(stem::AbstractString)
    m_scan  = match(r"_S(\d+)_", stem)
    m_frame = match(r"_(\d+)$", stem)
    scan_id  = m_scan  === nothing ? missing : parse(Int, m_scan.captures[1])
    frame_no = m_frame === nothing ? missing : parse(Int, m_frame.captures[1])
    return (scan_id, frame_no)
end

function _assign_sessions(start_times::AbstractVector; session_gap_hours::Float64 = 3.0)
    ids = Vector{Int}(undef, length(start_times))
    s = 1; prev = missing
    for (i, st) in enumerate(start_times)
        if !ismissing(st) && !ismissing(prev) && Dates.value(st - prev) / 3.6e6 > session_gap_hours
            s += 1
        end
        ids[i] = s
        ismissing(st) || (prev = st)
    end
    ids
end

_monitors(prp) = (
    transmitted_counts = prp === nothing ? missing : get(prp, :transmitted_counts, missing),
    incident_monitor   = prp === nothing ? missing : get(prp, :incident_monitor, missing),
    counting_time_s    = prp === nothing ? missing : get(prp, :counting_time_s, missing),
    pipe_length_m      = prp === nothing ? missing : get(prp, :pipe_length_m, missing),
)

function group_into_samples(metas::Vector{ExposureMeta})::GroupingResult
    isempty(metas) && return GroupingResult(GroupedLoad[], String[])

    load_groups, seg_flag = _segment_loads_with_flag(metas)
    discrepancies = String[]

    load_starts = map(load_groups) do lg
        ts = DateTime[m.prp.timestamp for m in lg if m.prp !== nothing && !ismissing(m.prp.timestamp)]
        isempty(ts) ? missing : minimum(ts)
    end
    session_ids = _assign_sessions(load_starts)

    if seg_flag == :unimodal_fallback
        push!(discrepancies,
            "load segmentation: no clear bimodal time-gap separation; treated the " *
            "entire directory as a single load — review load boundaries")
    end

    grouped_loads = GroupedLoad[]
    for (li, load_metas) in enumerate(load_groups)
        slot_groups = _cluster_slots(load_metas)
        grouped_samples = GroupedSample[]
        for (si, slot_metas) in enumerate(slot_groups)
            label_hint = _label_from_stem(first(slot_metas).stem)

            exps = [begin
                sid, fno = _parse_scan_frame(m.stem)
                ts  = (m.prp !== nothing && !ismissing(m.prp.timestamp)) ? m.prp.timestamp : missing
                ets = m.prp !== nothing ? m.prp.exposure_time_s : missing
                hpm = m.prp !== nothing ? m.prp.horizontal_position_mm : missing
                mon = _monitors(m.prp)
                GroupedExposure(
                    m.stem,
                    m.tif_path, m.prp_path, m.dat_path,
                    ts, ets, hpm,
                    mon.transmitted_counts, mon.incident_monitor, mon.counting_time_s,
                    mon.pipe_length_m,
                    sid, fno,
                )
            end for m in slot_metas]

            push!(grouped_samples, GroupedSample(
                _auto_name(label_hint, li, si),
                "auto",
                si,
                "auto_position",
                exps,
            ))
        end

        timestamps = DateTime[]
        for m in load_metas
            ts = m.prp !== nothing ? m.prp.timestamp : missing
            ismissing(ts) || push!(timestamps, ts)
        end
        start_time = isempty(timestamps) ? missing : minimum(timestamps)
        end_time   = isempty(timestamps) ? missing : maximum(timestamps)

        push!(grouped_loads, GroupedLoad(
            li,
            length(load_metas),
            start_time,
            end_time,
            grouped_samples,
            seg_flag,
            session_ids[li],
        ))
    end

    return GroupingResult(grouped_loads, discrepancies)
end
