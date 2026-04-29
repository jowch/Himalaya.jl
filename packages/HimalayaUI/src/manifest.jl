struct ManifestSample
    label          ::String
    name           ::String
    notes_sample   ::String
    notes_exposure ::String
    filenames      ::Vector{String}
end

"""
    expand_filename_range(s) -> Vector{String}

Expand a filename range like "JC001-004" or "JC013-JC016" to individual
filenames. Returns a single-element vector for plain filenames.
"""
function expand_filename_range(s::AbstractString)::Vector{String}
    m = match(r"^([A-Za-z]+)(\d+)-(?:[A-Za-z]*)(\d+)$", s)
    m === nothing && return [s]
    prefix, start_s, stop_s = m[1], m[2], m[3]
    width  = length(start_s)
    start  = parse(Int, start_s)
    stop   = parse(Int, stop_s)
    [string(prefix, lpad(i, width, '0')) for i in start:stop]
end

"""
    parse_manifest(cfg::ExperimentConfig, source) -> Vector{ManifestSample}

Parse a manifest CSV/TSV using `cfg` to drive delimiter, skip_rows, header
discovery, and column resolution. `source` may be an IO or a file path.

For each `[manifest]` column field, if the config value is a `String` and a
header row is present (`header_row > 0`), the column index is looked up by
that header name; otherwise the value is treated as a 1-based positional index.

Rows whose sample_id column does not parse as `Int` are silently skipped
(handles lab-notebook section headers).
"""
function parse_manifest(cfg::ExperimentConfig, source)::Vector{ManifestSample}
    lines = readlines(source)
    isempty(lines) && return ManifestSample[]

    cfg.skip_rows >= length(lines) && return ManifestSample[]
    body = lines[cfg.skip_rows+1:end]

    header_map = Dict{String,Int}()
    data_start = 1
    if cfg.header_row > 0
        hdr_idx = cfg.header_row - cfg.skip_rows
        if 1 <= hdr_idx <= length(body)
            for (i, c) in enumerate(split(body[hdr_idx], cfg.delimiter))
                header_map[strip(String(c))] = i
            end
            data_start = hdr_idx + 1
        end
    end

    function resolve_col(col::Union{Int,String})::Int
        if col isa AbstractString
            idx = get(header_map, col, 0)
            idx == 0 && @warn "Manifest column '$col' not found in header — field will be empty"
            return idx
        end
        return col
    end

    idx_id        = resolve_col(cfg.col_sample_id)
    idx_label     = resolve_col(cfg.col_label)
    idx_name      = resolve_col(cfg.col_name)
    idx_filenames = resolve_col(cfg.col_filenames)
    idx_notes_s   = resolve_col(cfg.col_notes_sample)
    idx_notes_e   = resolve_col(cfg.col_notes_exposure)

    function safe_get(cols, idx::Int)::String
        idx == 0 && return ""
        idx <= length(cols) ? String(strip(cols[idx])) : ""
    end

    samples = ManifestSample[]
    for line in body[data_start:end]
        cols = split(line, cfg.delimiter)
        id_str = safe_get(cols, idx_id)
        tryparse(Int, id_str) === nothing && continue

        filenames_str = safe_get(cols, idx_filenames)
        isempty(filenames_str) && continue

        push!(samples, ManifestSample(
            safe_get(cols, idx_label),
            safe_get(cols, idx_name),
            safe_get(cols, idx_notes_s),
            safe_get(cols, idx_notes_e),
            expand_filename_range(filenames_str),
        ))
    end
    samples
end

"""
    parse_manifest(source) -> Vector{ManifestSample}

Backward-compatible wrapper: parses using the built-in `simple` config (current
hardcoded behavior — tab-separated, columns 1/2/3/9/10/11, skip 1 row, no header).
For new code, prefer `parse_manifest(cfg::ExperimentConfig, source)`.
"""
parse_manifest(source)::Vector{ManifestSample} =
    parse_manifest(load_builtin_config("simple"), source)
