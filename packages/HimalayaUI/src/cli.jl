using ArgParse
using Printf
using DBInterface

function _resolve_experiment(db::SQLite.DB, key::Union{Nothing,AbstractString})
    if key === nothing
        rows = sqlrows(DBInterface.execute(db,
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

    looks_like_path = startswith(key, "/") || startswith(key, ".") || occursin('/', key)
    rows = if !isempty(key) && all(isdigit, key)
        sqlrows(DBInterface.execute(db,
            "SELECT id, name, path FROM experiments WHERE id = ?", [parse(Int, key)]))
    elseif looks_like_path
        plain = rstrip(abspath(key), '/')
        sqlrows(DBInterface.execute(db,
            """SELECT DISTINCT x.id, x.name, x.path FROM experiments x
                 LEFT JOIN setups st ON st.experiment_id = x.id
                WHERE rtrim(x.path, '/') IN (?1, ?2) OR rtrim(st.data_dir, '/') IN (?1, ?2)
                   OR rtrim(st.root_dir, '/') IN (?1, ?2)""",
            [plain, _canonical_dir(plain)]))
    else
        sqlrows(DBInterface.execute(db,
            "SELECT id, name, path FROM experiments WHERE name = ?", [key]))
    end
    isempty(rows)     && error("No experiment matching '$key'.")
    length(rows) > 1  && error("Multiple experiments matching '$key' — disambiguate by id.")
    rows[1]
end

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
        rows = sqlrows(DBInterface.execute(db,
            "SELECT id, name FROM experiments ORDER BY id"))
        isempty(rows) && error("No experiments registered.")
        skipped = 0
        for r in rows
            println("Experiment $(r.name):")
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
    skipped > 0 && error("$skipped image(s) or experiment(s) skipped; see the SKIP lines above")
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

    sid   = Int(sample_row.id)
    state = sample_call_state(db, sid)
    ms    = get(sample_call_members(db, [sid]), sid, [])
    phases = call_summary(state, ms).call_phases
    println("Call: ", export_call_state(state, ms), isempty(phases) ? "" : " ($(join(phases, ", ")))")

    exposures = get_exposures(db, sid)
    for exp_row in exposures
        e_id    = Int(exp_row.id)
        pks     = get_peaks_for_exposure(db, e_id)
        idxs    = get_indices_for_exposure(db, e_id)

        println("\nFrame: $(exp_row.filename)")
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
    fresh = !isfile(db_path)
    fresh && @warn "no database at $db_path — creating a new empty one. If you expected existing data, stop and check HIMALAYA_DB_PATH."

    db = open_db(db_path)
    println("HimalayaUI serving $(fresh ? "a new empty " : "")DB at $db_path on http://$(p[:host]):$(p[:port])")
    serve(db; host = p[:host], port = p[:port])
end

function parse_merge_args(args)
    s = ArgParseSettings(prog = "himalaya merge")
    @add_arg_table! s begin
        "experiments"
            nargs = '+'
            required = true
        "--name"
            default = nothing
        "--path"
            default = nothing
        "--apply"
            action = :store_true
    end
    parse_args(args, s; as_symbols = true)
end

function cli_merge(args; io::IO = stdout, db::Union{SQLite.DB,Nothing} = nothing)
    p = parse_merge_args(args)
    db = something(db, open_db(default_db_path()))
    ids = Int[Int(_resolve_experiment(db, k).id) for k in p[:experiments]]
    plan = merge_plan(db, ids; name = p[:name], path = p[:path])
    print_merge_plan(io, plan)
    if p[:apply]
        apply_merge!(db, plan)
        println(io, "Merged.")
    else
        println(io, "Dry run. Nothing was written. Run again with --apply to merge.")
    end
    plan
end

function parse_train_peak_model_args(args)
    s = ArgParseSettings(prog = "himalaya train-peak-model")
    @add_arg_table! s begin
        "--db"
            default = nothing
        "--trees"
            arg_type = Int
            default  = 200
        "--no-cache"
            dest_name = "no_cache"
            action = :store_true
    end
    parse_args(args, s; as_symbols = true)
end

function cli_train_peak_model(args; io::IO = stdout, log::IO = stderr)
    p = parse_train_peak_model_args(args)
    path = p[:db] === nothing ? default_db_path() : String(p[:db])
    isfile(path) || error("no database at $path")
    db = open_readonly_db(path)
    try
        result = run_peak_training(db; recipe = HimalayaPeakModel.default_recipe(n_trees = p[:trees]),
            cache_dir = p[:no_cache] ? nothing : peak_feature_cache_dir(path), io = log)
        JSON3.write(io, result)
        println(io)
    finally
        close(db)
    end
end

const _USAGE = """
Usage: himalaya <command> [options]

Commands:
  analyze   [-e <experiment> | --all] [-s <label>]  Run peak-finding + indexing
  show     [-e <experiment>] -s <label>     Print stored analysis for one sample
  serve    [--port N] [--host H]            Start the web server
  merge    <experiment>... [--name N] [--path P] [--apply]
                                            Merge per-setup experiments into one beamtime
  train-peak-model [--db PATH] [--trees N] [--no-cache]
                                            Train and evaluate a peak model read-only; prints one JSON object

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
    elseif cmd == "merge"
        cli_merge(args)
    elseif cmd == "train-peak-model"
        cli_train_peak_model(args)
    else
        println("Unknown command: $cmd")
        println()
        println(_USAGE)
    end
end
