const AUTOSAMPLER_STEM = r"^(.+)_(\d{2}[A-Z])_[SB]\d+$"
const AUTOSAMPLER_RUN_FILE = r"^AT_SolSAXS_.+_(\d+)\.txt$"

function autosampler_run_starts(data_dir::AbstractString)::Set{Int}
    names = try
        readdir(String(data_dir))
    catch
        String[]
    end
    Set{Int}(parse(Int, m.captures[1]) for m in (match(AUTOSAMPLER_RUN_FILE, f) for f in names) if m !== nothing)
end

function _acquisition_scan(stem::AbstractString)
    m = match(r"_[SB](\d+)$", stem)
    m === nothing ? missing : parse(Int, m.captures[1])
end

function _acquisition_label(stem::AbstractString)
    m = match(AUTOSAMPLER_STEM, stem)
    m === nothing ? nothing : (String(m.captures[1]), String(m.captures[2]))
end

function _acquisition_exposure(m::ExposureMeta)
    ts  = (m.prp !== nothing && !ismissing(m.prp.timestamp)) ? m.prp.timestamp : missing
    ets = m.prp !== nothing ? m.prp.exposure_time_s : missing
    hpm = m.prp !== nothing ? m.prp.horizontal_position_mm : missing
    mon = _monitors(m.prp)
    GroupedExposure(m.stem, m.tif_path, m.prp_path, m.dat_path, ts, ets, hpm,
        mon.transmitted_counts, mon.incident_monitor, mon.counting_time_s, mon.pipe_length_m,
        _acquisition_scan(m.stem), missing)
end

function _stable_hash(s::AbstractString)::Int
    h = 0x811c9dc5
    for b in codeunits(s)
        h = (h ⊻ UInt32(b)) * 0x01000193
    end
    Int(h)
end

function _acquisition_runs(sorted::Vector{ExposureMeta}, run_starts::Set{Int})
    runs = Vector{ExposureMeta}[]
    previous_label = nothing
    for m in sorted
        parsed = _acquisition_label(m.stem)
        label = parsed === nothing ? nothing : parsed[1]
        scan = _acquisition_scan(m.stem)
        new_run = isempty(runs) || label === nothing || label != previous_label || scan in run_starts
        new_run ? push!(runs, [m]) : push!(last(runs), m)
        previous_label = label
    end
    runs
end

function _acquisition_name(stem::AbstractString)
    m = match(r"^(.+)_(\d{2}[A-Z])_([SB])(\d+)$", stem)
    m === nothing ? String(stem) : "$(m.captures[1]) $(m.captures[2]) ($(m.captures[3])$(m.captures[4]))"
end

_acquisition_slot(m::ExposureMeta) = (s = _acquisition_scan(m.stem); ismissing(s) ? _stable_hash(m.stem) : s)

function group_autosampler(metas::Vector{ExposureMeta}, data_dir::AbstractString)::GroupingResult
    isempty(metas) && return GroupingResult(GroupedLoad[], String[])
    scanned = sort([m for m in metas if !ismissing(_acquisition_scan(m.stem))];
                   by = m -> (_acquisition_scan(m.stem), m.stem))
    scanless = sort([m for m in metas if ismissing(_acquisition_scan(m.stem))]; by = m -> m.stem)
    runs = _acquisition_runs(scanned, autosampler_run_starts(data_dir))
    indexed = Tuple{Int,Vector{ExposureMeta}}[(_acquisition_scan(first(run).stem), run) for run in runs]
    isempty(scanless) || push!(indexed, (0, scanless))
    stamps(run) = DateTime[m.prp.timestamp for m in run if m.prp !== nothing && !ismissing(m.prp.timestamp)]
    starts = [isempty(stamps(run)) ? missing : minimum(stamps(run)) for (_, run) in indexed]
    sessions = _assign_sessions(starts)
    loads = GroupedLoad[]
    for (i, (load_index, run)) in enumerate(indexed)
        samples = [GroupedSample(_acquisition_name(m.stem), "auto", _acquisition_slot(m), "auto_acquisition",
                                 [_acquisition_exposure(m)])
                   for m in run]
        ts = stamps(run)
        push!(loads, GroupedLoad(load_index, length(run),
            isempty(ts) ? missing : minimum(ts), isempty(ts) ? missing : maximum(ts),
            samples, :ok, sessions[i]))
    end
    GroupingResult(loads, String[])
end
