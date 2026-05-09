using ArgParse
using Printf
using DBInterface

"""
    cli_init_with_db!(db, exp_dir; analyze=true) -> experiment_id

Read `experiment.toml` from `exp_dir`, register the experiment in `db`, parse
the manifest (if present), and create samples and exposures via filesystem
discovery using the config's integration pattern. When `analyze=true` (the
default), runs peak-finding + indexing on every newly-created exposure.

This function is read-only with respect to `exp_dir` — it does not create,
modify, or delete any file inside it. All writes go to `db`.

The auto-analyze step is *not* wrapped in an outer transaction (each
`persist_analysis!` is itself atomic). A Ctrl-C mid-init leaves a registered
experiment with partial analysis; recover with `himalaya analyze -e <id>`,
which is idempotent for the unanalyzed exposures and skips the analyzed ones
that already have peaks.
"""
function cli_init_with_db!(db::SQLite.DB, exp_dir::String; analyze::Bool = true)::Int
    exp_dir   = abspath(exp_dir)
    toml_path = joinpath(exp_dir, "experiment.toml")
    isfile(toml_path) || error("experiment.toml not found in $exp_dir. Run 'himalaya config new --dir $exp_dir' first.")

    existing = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, name FROM experiments WHERE path = ?", [exp_dir]))
    if !isempty(existing)
        ids = join((string(r.id) for r in existing), ", ")
        error("$exp_dir is already registered (experiment id=$ids). " *
              "Run `himalaya reingest $exp_dir` to update it instead.")
    end

    cfg  = load_config(toml_path)
    blob = config_to_toml(cfg)

    data_dir     = isabspath(cfg.data_dir)     ? cfg.data_dir     : joinpath(exp_dir, cfg.data_dir)
    analysis_dir = isabspath(cfg.analysis_dir) ? cfg.analysis_dir : joinpath(exp_dir, cfg.analysis_dir)
    manifest_path = isabspath(cfg.manifest_file) ? cfg.manifest_file :
                    joinpath(exp_dir, cfg.manifest_file)

    exp_name = isempty(cfg.name) ? basename(exp_dir) : cfg.name

    exp_id = create_experiment!(db;
        name            = exp_name,
        path            = exp_dir,
        data_dir        = data_dir,
        analysis_dir    = analysis_dir,
        manifest_path   = isfile(manifest_path) ? manifest_path : nothing,
        config          = blob,
        experiment_type = cfg.exposure_type,
        energy_kev      = cfg.energy_kev,
        flight_path_m   = cfg.flight_path_m,
    )

    if isfile(manifest_path)
        samples = parse_manifest(cfg, manifest_path)
        sample_count = 0
        exposure_count = 0
        for ms in samples
            s_id = create_sample!(db;
                experiment_id = exp_id,
                name          = ms.name,
                display_name  = ms.display_name,
                notes         = ms.notes_sample)
            sample_count += 1

            for prefix in ms.filenames
                stems = resolve_files(cfg, analysis_dir, prefix, cfg.integration_pattern)
                if isempty(stems)
                    @warn "No integration files found for prefix '$prefix' in $analysis_dir"
                end
                for stem in stems
                    image_path = resolve_file_path(cfg, data_dir, stem, cfg.image_pattern)
                    e_id = create_exposure!(db;
                        sample_id  = s_id,
                        filename   = stem,
                        image_path = image_path)
                    exposure_count += 1

                    if !isempty(ms.notes_exposure)
                        DBInterface.execute(db,
                            "INSERT INTO exposure_tags (exposure_id, key, value, source) VALUES (?, 'note', ?, 'manifest')",
                            [e_id, ms.notes_exposure])
                    end
                end
            end
        end
        println("Imported $sample_count samples and $exposure_count exposures from $(basename(manifest_path)).")
    else
        println("No manifest at $manifest_path — experiment registered without samples.")
    end

    println("Initialized experiment '$exp_name' (id=$exp_id) at $exp_dir")

    if analyze
        println("Running analysis (peak-finding + indexing)...")
        _analyze_experiment!(db, exp_id)
    end

    exp_id
end

"""
    reingest!(db, experiment_id, exp_dir)

Re-read `experiment.toml` and the manifest CSV from `exp_dir` and update the
experiment row in `db`. Inserts new samples and exposures discovered on disk;
preserves existing exposures that have a non-NULL status (curation: accepted/rejected)
or any peaks (manual or auto-analysis already run).

Read-only with respect to `exp_dir`.
"""
function reingest!(db::SQLite.DB, experiment_id::Int, exp_dir::String)
    exp_dir   = abspath(exp_dir)
    toml_path = joinpath(exp_dir, "experiment.toml")
    isfile(toml_path) || error("experiment.toml not found in $exp_dir")
    SQLite.transaction(db) do
        _reingest_inner!(db, experiment_id, exp_dir, toml_path)
    end
end

function _reingest_inner!(db::SQLite.DB, experiment_id::Int, exp_dir::String, toml_path::String)
    cfg  = load_config(toml_path)
    blob = config_to_toml(cfg)

    data_dir     = isabspath(cfg.data_dir)     ? cfg.data_dir     : joinpath(exp_dir, cfg.data_dir)
    analysis_dir = isabspath(cfg.analysis_dir) ? cfg.analysis_dir : joinpath(exp_dir, cfg.analysis_dir)
    manifest_path = isabspath(cfg.manifest_file) ? cfg.manifest_file :
                    joinpath(exp_dir, cfg.manifest_file)

    exp_name = isempty(cfg.name) ? basename(exp_dir) : cfg.name

    DBInterface.execute(db,
        """UPDATE experiments
              SET name = ?, config = ?, experiment_type = ?,
                  energy_kev = ?, flight_path_m = ?,
                  data_dir = ?, analysis_dir = ?, manifest_path = ?
            WHERE id = ?""",
        [exp_name, blob, cfg.exposure_type, cfg.energy_kev, cfg.flight_path_m,
         data_dir, analysis_dir,
         isfile(manifest_path) ? manifest_path : nothing,
         experiment_id])

    if !isfile(manifest_path)
        return (status = :no_manifest, added_samples = 0,
                added_exposures = 0, manifest_path = manifest_path)
    end

    samples = parse_manifest(cfg, manifest_path)
    inserted_samples = 0
    inserted_exposures = 0

    for ms in samples
        # Upsert sample (match by name within experiment)
        existing = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM samples WHERE experiment_id = ? AND name = ?",
            [experiment_id, ms.name]))
        s_id = if isempty(existing)
            inserted_samples += 1
            create_sample!(db;
                experiment_id = experiment_id,
                name          = ms.name,
                display_name  = ms.display_name,
                notes         = ms.notes_sample)
        else
            DBInterface.execute(db,
                "UPDATE samples SET display_name = ?, notes = ? WHERE id = ?",
                [ms.display_name, ms.notes_sample, existing[1].id])
            Int(existing[1].id)
        end

        for prefix in ms.filenames
            stems = resolve_files(cfg, analysis_dir, prefix, cfg.integration_pattern)
            for stem in stems
                existing_exp = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, image_path FROM exposures WHERE sample_id = ? AND filename = ?",
                    [s_id, stem]))
                if isempty(existing_exp)
                    image_path = resolve_file_path(cfg, data_dir, stem, cfg.image_pattern)
                    create_exposure!(db;
                        sample_id  = s_id,
                        filename   = stem,
                        image_path = image_path)
                    inserted_exposures += 1
                else
                    # Existing exposures keep all curation. Only backfill image_path
                    # when it's currently NULL — typically because a previous reingest
                    # ran with a wrong image pattern and the file wasn't findable.
                    # Never overwrite a non-NULL image_path.
                    if ismissing(existing_exp[1].image_path) || isnothing(existing_exp[1].image_path)
                        new_image = resolve_file_path(cfg, data_dir, stem, cfg.image_pattern)
                        if new_image !== nothing
                            DBInterface.execute(db,
                                "UPDATE exposures SET image_path = ? WHERE id = ?",
                                [new_image, Int(existing_exp[1].id)])
                        end
                    end
                end
            end
        end
    end

    return (status = :ok, added_samples = inserted_samples,
            added_exposures = inserted_exposures, manifest_path = manifest_path)
end

function cli_reingest(args)
    s = ArgParseSettings(prog = "himalaya reingest")
    @add_arg_table! s begin
        "--experiment", "-e"
            help     = "experiment id, name, or path (required)"
            required = true
    end
    p = parse_args(args, s; as_symbols = true)

    db      = open_db()
    exp_row = _resolve_experiment(db, p[:experiment])
    exp_id  = Int(exp_row.id)
    exp_dir = String(exp_row.path)
    res     = reingest!(db, exp_id, exp_dir)
    if res.status === :no_manifest
        println("Reingested experiment $exp_id: config updated; no manifest at $(res.manifest_path).")
    else
        println("Reingested experiment $exp_id: +$(res.added_samples) samples, +$(res.added_exposures) exposures.")
    end
    res
end

function cli_init(args)
    s = ArgParseSettings(prog = "himalaya init")
    @add_arg_table! s begin
        "experiment_path"
            help     = "path to experiment directory containing experiment.toml"
            required = true
        "--no-analyze"
            help     = "skip auto-analysis after registration (run `himalaya analyze` later)"
            action   = :store_true
    end
    p = parse_args(args, s; as_symbols = true)
    exp_dir = p[:experiment_path]
    db = open_db()
    cli_init_with_db!(db, exp_dir; analyze = !p[Symbol("no-analyze")])
end

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
            error("No experiments registered. Run `himalaya init <path>` first.")
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

Shared between `cli_analyze` (explicit user invocation) and
`cli_init_with_db!` (auto-analyze after registering a fresh experiment).
"""
function _analyze_experiment!(db::SQLite.DB, exp_id::Int; sample_filter=nothing)
    exp     = get_experiment(db, exp_id)
    samples = get_samples(db, exp_id)
    sample_filter !== nothing && filter!(sm -> sm.label == sample_filter, samples)

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
                println("  Skipping $(sample.label) / $(exp_row.filename) (rejected)")
                continue
            end
            print("  Analyzing $(sample.label) / $(exp_row.filename) ... ")
            try
                analyze_exposure!(db, e_id, exp.analysis_dir)
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
            help    = "analyze only this sample label (e.g. D1)"
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
            help     = "sample label"
            required = true
    end
    p = parse_args(args, s; as_symbols = true)

    db      = open_db()
    exp_row = _resolve_experiment(db, p[:experiment])
    exp_id  = Int(exp_row.id)
    samples = get_samples(db, exp_id)
    idx     = findfirst(sm -> sm.label == p[:sample], samples)
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
    isfile(db_path) || error("no database at $db_path — run `himalaya init` first")

    db = open_db(db_path)
    println("HimalayaUI serving DB at $db_path on http://$(p[:host]):$(p[:port])")
    serve(db; host = p[:host], port = p[:port])
end

const _USAGE = """
Usage: himalaya <command> [options]

Commands:
  config new --type <name> --dir <path>     Create experiment.toml from a template
  config list                               List built-in config templates
  init <experiment_path> [--no-analyze]     Register an experiment in the central DB
                                            (auto-analyzes unless --no-analyze)
  reingest  -e <experiment>                 Re-read experiment.toml + manifest, update DB
  analyze   -e <experiment> [-s <label>]    Run peak-finding + indexing
  show     [-e <experiment>] -s <label>     Print stored analysis for one sample
  serve    [--port N] [--host H]            Start the web server

`-e <experiment>` accepts an id, name, or path. Required for write commands
(`reingest`, `analyze`); optional for the read-only `show` (defaults to the
sole registered experiment). Run `himalaya <command> --help` for full options.
Environment variables (HIMALAYA_DB_PATH, HIMALAYA_HOST, HIMALAYA_PORT, …) are
documented in packages/HimalayaUI/.env.example.
"""

function main(args = copy(ARGS))
    if isempty(args) || first(args) in ("--help", "-h", "help")
        println(_USAGE)
        return
    end
    cmd = popfirst!(args)

    if cmd == "init"
        cli_init(args)
    elseif cmd == "analyze"
        cli_analyze(args)
    elseif cmd == "show"
        cli_show(args)
    elseif cmd == "serve"
        cli_serve(args)
    elseif cmd == "config"
        cli_config(args)
    elseif cmd == "reingest"
        cli_reingest(args)
    else
        println("Unknown command: $cmd")
        println()
        println(_USAGE)
    end
end
