"""
    load_dat(path) -> (q, I, σ)

Parse a whitespace-delimited `.dat` trace file: columns q, I, σ. Blank lines and
`#` comments are skipped. Errors if a data row has fewer than 3 columns.
"""
function load_dat(path::AbstractString)
    q = Float64[]; I = Float64[]; σ = Float64[]
    for ln in eachline(path)
        s = strip(ln)
        (isempty(s) || startswith(s, '#')) && continue
        cols = split(s)
        length(cols) >= 3 || error("$path: expected ≥3 columns, got $(length(cols))")
        push!(q, parse(Float64, cols[1]))
        push!(I, parse(Float64, cols[2]))
        push!(σ, parse(Float64, cols[3]))
    end
    return (q, I, σ)
end

const DEFAULT_INTEGRATION_PATTERN = "{name}.dat"

_name_pattern(v) = (v isa AbstractString && occursin("{name}", v)) ? String(v) : nothing

"""
    load_trace(db, exposure_id; pattern=nothing) -> (q, I, σ)

Resolve an exposure's on-disk `.dat` path and parse it. The exposure's own
setup (`exposures.setup_id`) gives the directory and the file pattern, so each
setup of a beamtime resolves against its own `analysis_dir`. The candidate path
is `joinpath(setups.analysis_dir, replace(pat, "{name}" => exposures.filename))`,
where `pat` is the `pattern` keyword if given, else the setup's
`integration_pattern` when it contains `{name}`, else `"{name}.dat"`.

Mirrors `HimalayaUI.resolve_trace_path` (pipeline.jl): tries the exact name
first, then falls back to the name with the trailing `_0_001` suffix stripped.
Each exposure has its own `_tot.dat`, named by the stem without that suffix.

The file-not-found error names the path(s) it tried, which reveals a pattern
mismatch. Throws `ArgumentError` when the exposure is missing or has no setup.
"""
function load_trace(db::SQLite.DB, exposure_id::Integer;
                    pattern::Union{AbstractString,Nothing}=nothing)
    rows = Tables.rowtable(DBInterface.execute(db, """
        SELECT e.filename, st.analysis_dir, st.integration_pattern
        FROM exposures e
        LEFT JOIN setups st ON st.id = e.setup_id
        WHERE e.id = ?
    """, [Int(exposure_id)]))
    isempty(rows) && throw(ArgumentError("load_trace: no exposure $exposure_id"))
    row = rows[1]
    row.filename === missing && throw(ArgumentError(
        "load_trace: exposure $exposure_id has no filename"))
    row.analysis_dir === missing && throw(ArgumentError(
        "load_trace: exposure $exposure_id has no setup"))
    pat = pattern !== nothing ? String(pattern) :
          something(_name_pattern(coalesce(row.integration_pattern, nothing)), DEFAULT_INTEGRATION_PATTERN)
    dir = String(row.analysis_dir)
    fname = String(row.filename)

    exact = joinpath(dir, replace(pat, "{name}" => fname))
    isfile(exact) && return load_dat(exact)
    stripped = replace(fname, r"_\d+_\d+$" => "")
    if stripped != fname
        alt = joinpath(dir, replace(pat, "{name}" => stripped))
        isfile(alt) && return load_dat(alt)
    end
    throw(ArgumentError(
        "load_trace: no trace for exposure $exposure_id (tried $exact" *
        (stripped != fname ? " and its suffix-stripped name" : "") *
        "; wrong `pattern=`, or data dir not present?)"))
end
