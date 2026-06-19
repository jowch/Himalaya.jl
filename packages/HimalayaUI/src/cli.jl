using ArgParse
using Printf
using DBInterface

"""
    _resolve_experiment(db, key) -> (id::Int, name::String, path::String)

Resolve an experiment from the central DB. `key` may be:
  - `nothing` — picks the sole experiment if there's exactly one; errors with
    a listing otherwise.
  - a numeric string — looked up by `experiments.id`.
  - a path-like string (starts with `/` or `.`, or contains a separator) —
    looked up by `experiments.path` after `abspath` (back-compat with the
    old positional-path CLI).
  - any other string — looked up by `experiments.name`.
"""
function _resolve_experiment(db::SQLite.DB, key::Union{Nothing,AbstractString})
    if key === nothing
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name, path FROM experiments ORDER BY id"))
        if length(rows) == 1
            return rows[1]
        elseif isempty(rows)
            error("No experiments registered. Ingest one via the HTTP scan API first.")
        else
            io = IOBuffer()
            println(io, "Multiple experiments registered — specify which:")
            for r in rows
                println(io, "  [$(r.id)]  $(r.name)  ($(r.path))")
            end
            error(String(take!(io)))
        end
    end

    # Heuristic: anything containing `/` is a path, not a name. Experiments
    # whose name happens to contain `/` would be misclassified — pass `-e`
    # with the numeric id to disambiguate in that (currently theoretical) case.
    looks_like_path = startswith(key, "/") || startswith(key, ".") || occursin('/', key)
    rows = if !isempty(key) && all(isdigit, key)
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name, path FROM experiments WHERE id = ?", [parse(Int, key)]))
    elseif looks_like_path
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name, path FROM experiments WHERE path = ?", [abspath(key)]))
    else
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name, path FROM experiments WHERE name = ?", [key]))
    end
    isempty(rows)     && error("No experiment matching '$key'.")
    length(rows) > 1  && error("Multiple experiments matching '$key' — disambiguate by id.")
    rows[1]
end

"""
    _analyze_experiment!(db, exp_id; sample_filter=nothing) -> nothing

Run peak-finding + indexing for every exposure of `exp_id`, optionally
restricted to a single sample label. Skips rejected exposures and auto-
selects the first accepted exposure when none is currently selected.
Errors per-exposure are caught and printed as `SKIP (...)` so one bad
.dat file doesn't abort the batch.

Invoked by `cli_analyze` (explicit user invocation). Ingestion is handled by
the HTTP scan path (`POST /api/experiments/{id}/scan`), not the CLI.
"""
function _analyze_experiment!(db::SQLite.DB, exp_id::Int; sample_filter=nothing)
    exp     = get_experiment(db, exp_id)
    samples = get_samples(db, exp_id)
    sample_filter !== nothing && filter!(sm -> sm.name == sample_filter, samples)

    for sample in samples
        exposures = get_exposures(db, Int(sample.id))

        has_selected = any(e -> Int(e.selected) == 1, exposures)
        if !has_selected
            first_accepted = findfirst(
                e -> !ismissing(e.status) && e.status == "accepted", exposures)
            if first_accepted !== nothing
                fallback_id = Int(exposures[first_accepted].id)
                DBInterface.execute(db,
                    "UPDATE exposures SET selected = 0 WHERE sample_id = ?",
                    [Int(sample.id)])
                DBInterface.execute(db,
                    "UPDATE exposures SET selected = 1 WHERE id = ?",
                    [fallback_id])
            end
        end

        for exp_row in exposures
            e_id = Int(exp_row.id)
            e_status = ismissing(exp_row.status) ? nothing : exp_row.status
            if e_status == "rejected"
                println("  Skipping $(sample.name) / $(exp_row.filename) (rejected)")
                continue
            end
            print("  Analyzing $(sample.name) / $(exp_row.filename) ... ")
            try
                analyze_exposure!(db, e_id)
                println("done")
            catch e
                msg = isa(e, ErrorException) ? e.msg : sprint(showerror, e)
                println("SKIP ($msg)")
            end
        end
    end
end

function cli_analyze(args)
    s = ArgParseSettings(prog = "himalaya analyze")
    @add_arg_table! s begin
        "--experiment", "-e"
            help     = "experiment id, name, or path (required)"
            required = true
        "--sample", "-s"
            help    = "sample name (stable identifier) (e.g. D1)"
            default = nothing
    end
    p = parse_args(args, s; as_symbols = true)

    db      = open_db()
    exp_row = _resolve_experiment(db, p[:experiment])
    _analyze_experiment!(db, Int(exp_row.id); sample_filter = p[:sample])
end

function cli_show(args)
    s = ArgParseSettings(prog = "himalaya show")
    @add_arg_table! s begin
        "--experiment", "-e"
            help    = "experiment id, name, or path (default: the sole registered experiment)"
            default = nothing
        "--sample", "-s"
            help     = "sample name (stable identifier)"
            required = true
    end
    p = parse_args(args, s; as_symbols = true)

    db      = open_db()
    exp_row = _resolve_experiment(db, p[:experiment])
    exp_id  = Int(exp_row.id)
    samples = get_samples(db, exp_id)
    idx     = findfirst(sm -> sm.name == p[:sample], samples)
    idx === nothing && error("sample $(p[:sample]) not found")
    sample_row = samples[idx]

    exposures = get_exposures(db, Int(sample_row.id))
    for exp_row in exposures
        e_id    = Int(exp_row.id)
        pks     = get_peaks_for_exposure(db, e_id)
        idxs    = get_indices_for_exposure(db, e_id)
        groups  = get_groups_for_exposure(db, e_id)

        println("\nExposure: $(exp_row.filename)")
        println("  Peaks ($(length(pks))):")
        for pk in pks
            @printf "    q=%.4f  prom=%.3f  sharp=%.3f  [%s]\n" pk.q pk.prominence pk.sharpness pk.source
        end
        println("  Indices ($(length(idxs))):")
        for ix in idxs
            @printf "    %-6s  basis=%.4f  score=%.3f  R²=%.4f  d=%.2f\n" ix.phase ix.basis ix.score ix.r_squared ix.lattice_d
        end
        active = findfirst(g -> g.active == 1, groups)
        if active !== nothing
            println("  Active group: $(groups[active].kind)")
        end
    end
end

"""
    cli_config_list(io=stdout)

Print the names of all built-in config templates, one per line.
"""
function cli_config_list(io::IO = stdout)
    types = list_config_types()
    if isempty(types)
        println(io, "(no built-in config types found)")
    else
        for t in types
            println(io, t)
        end
    end
end

"""
    cli_config_new(; type_name::String="simple", dir::String)

Copy the named built-in template to `<dir>/experiment.toml`. Errors if
the destination already exists with content (will not overwrite). An
empty or whitespace-only file is treated as a placeholder and gets
filled — this supports read-only experiment dirs where a sysadmin
pre-creates the file with the right permissions for a curator to
populate. This is the only documented operation that writes to an
experiment directory.
"""
function cli_config_new(; type_name::String = "simple", dir::String)
    isdir(dir) || error("Directory not found: $dir")
    dest = joinpath(dir, "experiment.toml")
    placeholder = isfile(dest)
    if placeholder && !isempty(strip(read(dest, String)))
        error("experiment.toml already exists at $dest with content — will not overwrite")
    end
    src = joinpath(configs_dir(), type_name * ".toml")
    isfile(src) || error("Unknown config type '$type_name'. Run 'himalaya config list' to see options.")
    # When filling a pre-existing placeholder, write into the existing file
    # rather than unlink+copy — the parent dir may be read-only (which is
    # the whole reason the placeholder exists). For new files we cp, so
    # mode/ownership stay sensible.
    if placeholder
        write(dest, read(src, String))
    else
        cp(src, dest)
    end
    println("Created $dest from template '$type_name'")
    println("Edit it to set your experiment name, beamline parameters, and manifest column mappings.")
    dest
end

function cli_config(args)
    if isempty(args) || first(args) in ("--help", "-h", "help")
        println("""
Usage: himalaya config <subcommand> [options]

Subcommands:
  list                                  List built-in config templates
  new --type <name> --dir <path>        Copy a template to <path>/experiment.toml
""")
        return
    end
    sub = popfirst!(args)
    if sub == "list"
        cli_config_list()
    elseif sub == "new"
        s = ArgParseSettings(prog = "himalaya config new")
        @add_arg_table! s begin
            "--type"
                default = "simple"
                help    = "built-in template name (see 'himalaya config list')"
            "--dir"
                required = true
                help     = "directory in which to write experiment.toml"
        end
        p = parse_args(args, s; as_symbols = true)
        cli_config_new(type_name = p[:type], dir = p[:dir])
    else
        println("Unknown config subcommand: $sub. Available: list, new")
    end
end

function cli_serve(args)
    s = ArgParseSettings(prog = "himalaya serve")
    @add_arg_table! s begin
        "--port"
            arg_type = Int
            default  = parse(Int, get(ENV, "HIMALAYA_PORT", "8080"))
        "--host"
            default  = get(ENV, "HIMALAYA_HOST", "127.0.0.1")
    end
    p = parse_args(args, s; as_symbols = true)

    db_path = default_db_path()
    isfile(db_path) || error("no database at $db_path — ingest an experiment via the HTTP scan API first")

    db = open_db(db_path)
    println("HimalayaUI serving DB at $db_path on http://$(p[:host]):$(p[:port])")
    serve(db; host = p[:host], port = p[:port])
end

"""
    cli_migrate_toml(args)

Rewrite `<dir>/experiment.toml` from the legacy `[manifest].label/name` shape
to the canonical `[manifest].name/display_name` shape. Section-aware: only
substitutes inside the `[manifest]` block so the "axis label units" comment
in `[beamline]` is not misfired. Idempotent. Atomic file write.
"""
function cli_migrate_toml(args::Vector{<:AbstractString})
    isempty(args) && error("Usage: himalaya migrate-toml <experiment-dir>")
    dir  = args[1]
    path = joinpath(dir, "experiment.toml")
    isfile(path) || error("experiment.toml not found in $dir")

    text = read(path, String)
    new_text, changed = try
        migrate_manifest_toml_text(text)
    catch err
        # Re-raise with a path-aware prefix; intentionally a fresh exception
        # rather than `rethrow(err)` because the helper's bare message lacks
        # the file path that operators need in this CLI context.
        throw(ErrorException("experiment.toml at $path: $(sprint(showerror, err))"))
    end
    if !changed
        @info "experiment.toml at $path already migrated (or no `[manifest].label` to migrate)"
        return nothing
    end

    tmp = path * ".tmp"
    open(tmp, "w") do io; print(io, new_text); end
    mv(tmp, path; force=true)
    @info "Migrated $path"
    nothing
end

const _USAGE = """
Usage: himalaya <command> [options]

Commands:
  config new --type <name> --dir <path>     Create experiment.toml from a template
  config list                               List built-in config templates
  analyze   -e <experiment> [-s <label>]    Run peak-finding + indexing
  show     [-e <experiment>] -s <label>     Print stored analysis for one sample
  serve    [--port N] [--host H]            Start the web server

Ingestion is performed via the HTTP scan API (`POST /api/experiments/{id}/scan`),
not the CLI. `-e <experiment>` accepts an id, name, or path. Required for
`analyze`; optional for the read-only `show` (defaults to the sole registered
experiment). Run `himalaya <command> --help` for full options.
Environment variables (HIMALAYA_DB_PATH, HIMALAYA_HOST, HIMALAYA_PORT, …) are
documented in packages/HimalayaUI/.env.example.
"""

function main(args = copy(ARGS))
    if isempty(args) || first(args) in ("--help", "-h", "help")
        println(_USAGE)
        return
    end
    cmd = popfirst!(args)

    if cmd == "analyze"
        cli_analyze(args)
    elseif cmd == "show"
        cli_show(args)
    elseif cmd == "serve"
        cli_serve(args)
    elseif cmd == "config"
        cli_config(args)
    elseif cmd == "migrate-toml"
        cli_migrate_toml(args)
    else
        println("Unknown command: $cmd")
        println()
        println(_USAGE)
    end
end
