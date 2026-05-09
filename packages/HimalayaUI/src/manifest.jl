using CSV
using Tables

struct ManifestSample
    name           ::String   # was: label — stable identifier (e.g. "JC001")
    display_name   ::String   # was: name  — user-friendly label (e.g. "DOPC + chol")
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
    expand_filename_field(s) -> Vector{String}

Expand a manifest filenames cell into individual filename stems. Supports
one or more ranges separated by `;` or `,` — e.g. `"JC037-040;JC153-156"`
or (CSV-quoted) `"JC037-JC040,JC153-156"`. Each segment may be a single
filename or a range like `JC001-004` / `JC013-JC016`. Surrounding
whitespace is stripped and empty segments are skipped.
"""
function expand_filename_field(s::AbstractString)::Vector{String}
    out = String[]
    for seg in split(s, r"[;,]")
        seg = strip(seg)
        isempty(seg) && continue
        append!(out, expand_filename_range(seg))
    end
    out
end

"""
    parse_manifest(cfg::ExperimentConfig, source) -> Vector{ManifestSample}

Parse a manifest CSV/TSV using `cfg` to drive delimiter, skip_rows, header
discovery, and column resolution. `source` may be an IO or a file path.

Row tokenization is delegated to CSV.jl, so RFC 4180-style quoted fields
(commas inside `"..."` cells) work for free. The filenames cell may
contain multiple ranges separated by `;` or `,` — see
[`expand_filename_field`](@ref).

Column resolution happens in two steps. First, CSV.jl produces a list of
column names: the actual header values when `header_row > 0`, or synthetic
`Column1`, `Column2`, … when `header_row == 0` (positional mode). Second,
each `[manifest]` column entry from `cfg` is resolved against that list:
a `String` value matches a column by name (warns and yields an empty field
on miss), and an `Int` value selects the i-th column (1-based). Cells are
then accessed via `getproperty(row, sym)` against the resolved Symbols.

Rows whose sample_id column does not parse as `Int` are silently skipped
(handles lab-notebook section headers and stray preamble rows).
"""
function parse_manifest(cfg::ExperimentConfig, source)::Vector{ManifestSample}
    cf = if cfg.header_row > 0
        # Named columns: row `header_row` is the header. CSV.jl skips
        # preamble rows above it automatically; `skip_rows` is implicit.
        CSV.File(source;
            delim = cfg.delimiter, header = cfg.header_row,
            types = String, stringtype = String,
            silencewarnings = true,
        )
    else
        # Positional columns: skip `skip_rows` lines, then treat the next
        # row as the start of data. CSV.jl assigns synthetic Column1,
        # Column2, ... names which we resolve by index.
        CSV.File(source;
            delim = cfg.delimiter, header = false,
            skipto = cfg.skip_rows + 1,
            types = String, stringtype = String,
            silencewarnings = true,
        )
    end

    isempty(cf) && return ManifestSample[]
    column_syms = collect(Tables.columnnames(cf))

    function resolve_col(col::Union{Int,String})::Symbol
        if col isa AbstractString
            sym = Symbol(col)
            sym in column_syms || (
                @warn "Manifest column '$col' not found in header — field will be empty";
                return Symbol("")
            )
            return sym
        end
        1 <= col <= length(column_syms) ? column_syms[col] : Symbol("")
    end

    sym_id               = resolve_col(cfg.col_sample_id)
    sym_name_col         = resolve_col(cfg.col_name)
    sym_display_name_col = resolve_col(cfg.col_display_name)
    sym_filenames        = resolve_col(cfg.col_filenames)
    sym_notes_s          = resolve_col(cfg.col_notes_sample)
    sym_notes_e          = resolve_col(cfg.col_notes_exposure)

    function safe_get(row, sym::Symbol)::String
        sym == Symbol("") && return ""
        v = getproperty(row, sym)
        v === missing ? "" : String(strip(v))
    end

    samples = ManifestSample[]
    for row in cf
        id_str = safe_get(row, sym_id)
        tryparse(Int, id_str) === nothing && continue

        filenames_str = safe_get(row, sym_filenames)
        isempty(filenames_str) && continue

        push!(samples, ManifestSample(
            safe_get(row, sym_name_col),         # column 2 → name field (stable identifier)
            safe_get(row, sym_display_name_col), # column 3 → display_name field (user-friendly label)
            safe_get(row, sym_notes_s),
            safe_get(row, sym_notes_e),
            expand_filename_field(filenames_str),
        ))
    end
    samples
end

"""
    parse_manifest(source) -> Vector{ManifestSample}

Backward-compatible wrapper: parses using the built-in `simple` config
(tab-separated, columns 1/2/3/9/10/11, skip 1 row, no header). New code
should prefer `parse_manifest(cfg::ExperimentConfig, source)`.
"""
parse_manifest(source)::Vector{ManifestSample} =
    parse_manifest(load_builtin_config("simple"), source)
