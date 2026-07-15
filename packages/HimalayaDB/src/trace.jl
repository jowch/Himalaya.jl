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

"""
    load_trace(db, exposure_id; pattern=nothing) -> (q, I, σ)

Resolve an exposure's on-disk `.dat` path and parse it. `exposures.experiment_id`
is joined directly to `experiments` (an exposure can have a NULL `sample_id`, so
going through `samples` would miss it). The candidate path is
`joinpath(experiments.analysis_dir, replace(pat, "{name}" => exposures.filename))`,
where `pat` is the `pattern` keyword if given, else the experiment's
`integration_pattern` column, else `"{name}.dat"` when that column is NULL.

Mirrors `HimalayaUI.resolve_trace_path` (pipeline.jl): tries the exact per-frame
name first, then falls back to the per-acquisition name with the trailing
`_<rep>_<frame>` suffix stripped — the shared `_tot.dat` total integration is
named by the acquisition stem, not the per-frame stem.

Pass `pattern=` to override (e.g. `"{name}_tot.dat"`). The file-not-found error
names the path(s) it tried, which reveals a pattern mismatch.
"""
function load_trace(db::SQLite.DB, exposure_id::Integer;
                    pattern::Union{AbstractString,Nothing}=nothing)
    rows = Tables.rowtable(DBInterface.execute(db, """
        SELECT e.filename, x.analysis_dir, x.integration_pattern
        FROM exposures e
        JOIN experiments x ON x.id = e.experiment_id
        WHERE e.id = ?
    """, [Int(exposure_id)]))
    isempty(rows) && throw(ArgumentError("load_trace: no exposure $exposure_id"))
    row = rows[1]
    row.filename === missing && throw(ArgumentError(
        "load_trace: exposure $exposure_id has no filename"))
    pat = pattern !== nothing ? String(pattern) :
          (row.integration_pattern === missing ? "{name}.dat" : String(row.integration_pattern))
    dir = String(row.analysis_dir)
    fname = String(row.filename)

    # Mirror HimalayaUI.resolve_trace_path (pipeline.jl): exact per-frame name,
    # then the per-acquisition name with the trailing _<rep>_<frame> suffix
    # stripped — the shared _tot.dat is named by the acquisition stem.
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
