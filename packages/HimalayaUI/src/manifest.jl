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
    parse_manifest(io_or_path) -> Vector{ManifestSample}

Parse a tab-separated manifest CSV exported from Google Sheets.
Skips section header rows (rows where the first column is non-numeric or empty).
"""
function parse_manifest(source)::Vector{ManifestSample}
    lines = readlines(source)
    samples = ManifestSample[]
    for line in lines[2:end]
        cols = split(line, '\t')
        length(cols) < 9 && continue
        num_str = strip(cols[1])
        tryparse(Int, num_str) === nothing && continue

        label          = String(strip(get(cols, 2,  "")))
        name           = String(strip(get(cols, 3,  "")))
        notes_sample   = String(strip(get(cols, 10, "")))
        notes_exposure = String(strip(get(cols, 11, "")))
        filename_str   = String(strip(get(cols, 9,  "")))

        isempty(filename_str) && continue
        filenames = expand_filename_range(filename_str)

        push!(samples, ManifestSample(label, name, notes_sample, notes_exposure, filenames))
    end
    samples
end
