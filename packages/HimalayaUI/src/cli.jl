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
    matched against `experiments.path` and `experiments.data_dir`, as typed
    (after `abspath`) and with symlinks resolved, ignoring a trailing slash.
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
        plain = rstrip(abspath(key), '/')
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name, path FROM experiments WHERE rtrim(path, '/') IN (?, ?) OR rtrim(data_dir, '/') IN (?, ?)",
            [plain, _canonical_dir(plain), plain, _canonical_dir(plain)]))
    else
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name, path FROM experiments WHERE name = ?", [key]))
    end
    isempty(rows)     && error("No experiment matching '$key'.")
    length(rows) > 1  && error("Multiple experiments matching '$key' — disambiguate by id.")
    rows[1]
end

"""
    _analyze_experiment!(db, exp_id; sample_filter=nothing) -> Int

Run peak-finding + indexing for every exposure of `exp_id`, optionally
restricted to a single sample label. Skips rejected exposures and leaves
`selected` alone: the representative is derived by `representative_exposure_id`.
Errors per-exposure are caught and printed as `SKIP (...)` so one bad
.dat file doesn't abort the batch.

Reports any assignment members the rebuild failed to re-attach (F-WIPE W1),
per-exposure and again as a closing summary — a withdrawn phase call is a human
decision disappearing, and this is the only surface that shows it on the CLI.

Invoked by `cli_analyze` (explicit user invocation). Ingestion is handled by
the HTTP scan path (`POST /api/experiments/{id}/scan`), not the CLI.
"""
function _analyze_experiment!(db::SQLite.DB, exp_id::Int; sample_filter=nothing)
    exp     = get_experiment(db, exp_id)
    samples = get_samples(db, exp_id)
    dropped_assignments = Tuple{String, String, String}[]
    skipped = 0
    sample_filter !== nothing && filter!(sm -> sm.name == sample_filter, samples)
    sample_filter !== nothing && isempty(samples) &&
        println("  (no sample named '$(sample_filter)')")

    for sample in samples
        exposures = get_exposures(db, Int(sample.id))

        for exp_row in exposures
            e_id = Int(exp_row.id)
            e_status = ismissing(exp_row.status) ? nothing : exp_row.status
            if e_status == "rejected"
                println("  Skipping $(sample.name) / $(exp_row.filename) (rejected)")
                continue
            end
            print("  Analyzing $(sample.name) / $(exp_row.filename) ... ")
            try
                res = analyze_exposure!(db, e_id)
                # A rebuild can fail to re-attach an assignment member (F-WIPE
                # W1) — the human's phase call is withdrawn. The routes thread
                # this into the SSE post_state; on the CLI it is the only signal
                # a curated call disappeared, so it must not be silent. Reached
                # in practice by the #304 migration: a Hexagonal call resting on
                # the removed √11 can drop below minpeaks and stop existing.
                if isempty(res.dropped_assignment_phases)
                    println("done")
                else
                    dropped = join(sort(res.dropped_assignment_phases), ", ")
                    println("done — DROPPED ASSIGNMENT: $dropped")
                    push!(dropped_assignments, (sample.name, exp_row.filename, dropped))
                end
            catch e
                msg = isa(e, ErrorException) ? e.msg : sprint(showerror, e)
                println("SKIP ($msg)")
                skipped += 1
            end
        end
    end

    if !isempty(dropped_assignments)
        println()
        println("$(length(dropped_assignments)) assignment(s) dropped — a curated phase call no longer reproduces:")
        for (sname, fname, phases) in dropped_assignments
            println("  $sname / $fname — $phases")
        end
    end
    return skipped
end

function cli_analyze(args)
    s = ArgParseSettings(prog = "himalaya analyze")
    @add_arg_table! s begin
        "--experiment", "-e"
            help    = "experiment id, name, or path"
            default = nothing
        "--sample", "-s"
            help    = "sample name (stable identifier) (e.g. D1)"
            default = nothing
        "--all"
            help   = "analyze every registered experiment"
            action = :store_true
    end
    p = parse_args(args, s; as_symbols = true)

    db = open_db()
    if p[:all]
        p[:experiment] === nothing ||
            error("--all and --experiment are mutually exclusive.")
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name FROM experiments ORDER BY id"))
        isempty(rows) && error("No experiments registered.")
        skipped = 0
        for r in rows
            println("Experiment $(r.name):")
            # Same catch-print-continue philosophy as the per-exposure loop:
            # one bad experiment must not abort the rest of the batch.
            try
                skipped += _analyze_experiment!(db, Int(r.id); sample_filter = p[:sample])
            catch e
                msg = isa(e, ErrorException) ? e.msg : sprint(showerror, e)
                println("  SKIP experiment ($msg)")
                skipped += 1
            end
        end
    else
        exp_row = _resolve_experiment(db, p[:experiment])
        skipped = _analyze_experiment!(db, Int(exp_row.id); sample_filter = p[:sample])
    end
    skipped > 0 && error("$skipped exposure(s) or experiment(s) skipped; see the SKIP lines above")
end

_fmt(x, spec) = ismissing(x) || x === nothing ? "n/a" : Printf.format(Printf.Format(spec), x)

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
            println("    q=", _fmt(pk.q, "%.4f"), "  prom=", _fmt(pk.prominence, "%.3f"),
                    "  sharp=", _fmt(pk.sharpness, "%.3f"), "  [", pk.source, "]")
        end
        println("  Indices ($(length(idxs))):")
        for ix in idxs
            println("    ", rpad(ix.phase, 6), "  basis=", _fmt(ix.basis, "%.4f"),
                    "  score=", _fmt(ix.score, "%.3f"), "  R²=", _fmt(ix.r_squared, "%.4f"),
                    "  d=", _fmt(ix.lattice_d, "%.2f"))
        end
        active = findfirst(g -> g.active == 1, groups)
        if active !== nothing
            println("  Active group: $(groups[active].kind)")
        end
    end
end

function parse_serve_args(args)
    s = ArgParseSettings(prog = "himalaya serve")
    @add_arg_table! s begin
        "--port"
            arg_type = Int
            default  = parse(Int, get(ENV, "HIMALAYA_PORT", "8080"))
        "--host"
            default  = get(ENV, "HIMALAYA_HOST", "127.0.0.1")
    end
    parse_args(args, s; as_symbols = true)
end

function cli_serve(args)
    p = parse_serve_args(args)

    db_path = default_db_path()
    # First-run bootstrap: serve an empty schema'd DB when none exists yet, so a
    # brand-new user lands on the empty experiments home and can create one via
    # the funnel (you can't reach the HTTP scan API without a running server, and
    # the server can't start without a DB — a chicken-and-egg the old hard error
    # left unresolvable). open_db creates + migrates. A loud notice keeps a typo'd
    # HIMALAYA_DB_PATH from silently masking existing data behind a fresh empty DB.
    fresh = !isfile(db_path)
    fresh && @warn "no database at $db_path — creating a new empty one. If you expected existing data, stop and check HIMALAYA_DB_PATH."

    db = open_db(db_path)
    println("HimalayaUI serving $(fresh ? "a new empty " : "")DB at $db_path on http://$(p[:host]):$(p[:port])")
    serve(db; host = p[:host], port = p[:port])
end

const _USAGE = """
Usage: himalaya <command> [options]

Commands:
  analyze   [-e <experiment> | --all] [-s <label>]  Run peak-finding + indexing
  show     [-e <experiment>] -s <label>     Print stored analysis for one sample
  serve    [--port N] [--host H]            Start the web server

Ingestion is performed via the HTTP scan API (`POST /api/experiments/{id}/scan`),
not the CLI. `-e <experiment>` accepts an id, name, or path; when omitted it
defaults to the sole registered experiment. `analyze --all` runs every
registered experiment. Run `himalaya <command> --help` for full options.
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
    else
        println("Unknown command: $cmd")
        println()
        println(_USAGE)
    end
end
