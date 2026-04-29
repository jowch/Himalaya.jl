using ArgParse
using Printf
using DBInterface

"""
    cli_init_with_db!(db, exp_dir) -> experiment_id

Read `experiment.toml` from `exp_dir`, register the experiment in `db`, parse
the manifest (if present), and create samples and exposures via filesystem
discovery using the config's integration pattern.

This function is read-only with respect to `exp_dir` — it does not create,
modify, or delete any file inside it. All writes go to `db`.
"""
function cli_init_with_db!(db::SQLite.DB, exp_dir::String)::Int
    exp_dir   = abspath(exp_dir)
    toml_path = joinpath(exp_dir, "experiment.toml")
    isfile(toml_path) || error("experiment.toml not found in $exp_dir. Run 'himalaya config new --dir $exp_dir' first.")

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
                label         = ms.label,
                name          = ms.name,
                notes         = ms.notes_sample)
            sample_count += 1

            for prefix in ms.filenames
                stems = resolve_files(cfg, analysis_dir, prefix, cfg.integration_pattern)
                if isempty(stems)
                    @warn "No integration files found for prefix '$prefix' in $analysis_dir"
                end
                for stem in stems
                    image_rel = replace(cfg.image_pattern, "{name}" => stem)
                    image_full = joinpath(data_dir, image_rel)
                    image_path = isfile(image_full) ? image_full : nothing
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

    isfile(manifest_path) || (println("No manifest at $manifest_path — config updated only."); return)

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
                label         = ms.label,
                name          = ms.name,
                notes         = ms.notes_sample)
        else
            DBInterface.execute(db,
                "UPDATE samples SET label = ?, notes = ? WHERE id = ?",
                [ms.label, ms.notes_sample, existing[1].id])
            Int(existing[1].id)
        end

        for prefix in ms.filenames
            stems = resolve_files(cfg, analysis_dir, prefix, cfg.integration_pattern)
            for stem in stems
                existing_exp = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id FROM exposures WHERE sample_id = ? AND filename = ?",
                    [s_id, stem]))
                if isempty(existing_exp)
                    image_rel = replace(cfg.image_pattern, "{name}" => stem)
                    image_full = joinpath(data_dir, image_rel)
                    image_path = isfile(image_full) ? image_full : nothing
                    create_exposure!(db;
                        sample_id  = s_id,
                        filename   = stem,
                        image_path = image_path)
                    inserted_exposures += 1
                end
                # Existing exposures: do NOT modify — they may carry curation/peaks.
            end
        end
    end

    println("Reingested experiment $experiment_id: +$inserted_samples samples, +$inserted_exposures exposures.")
end

function cli_reingest(args)
    s = ArgParseSettings(prog = "himalaya reingest")
    @add_arg_table! s begin
        "experiment_path"
            help     = "experiment directory (must already be registered via 'himalaya init')"
            required = true
    end
    p = parse_args(args, s; as_symbols = true)
    exp_dir = abspath(p[:experiment_path])

    db = open_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM experiments WHERE path = ?", [exp_dir]))
    isempty(rows) && error("No experiment registered at $exp_dir. Run 'himalaya init' first.")
    reingest!(db, Int(rows[1].id), exp_dir)
end

function cli_init(args)
    s = ArgParseSettings(prog = "himalaya init")
    @add_arg_table! s begin
        "experiment_path"
            help     = "path to experiment directory containing experiment.toml"
            required = true
    end
    p = parse_args(args, s; as_symbols = true)
    exp_dir = p[:experiment_path]
    db = open_db()
    cli_init_with_db!(db, exp_dir)
end

function cli_analyze(args)
    s = ArgParseSettings(prog = "himalaya analyze")
    @add_arg_table! s begin
        "experiment_path"
            required = true
        "--sample", "-s"
            help    = "analyze only this sample label (e.g. D1)"
            default = nothing
    end
    p = parse_args(args, s; as_symbols = true)

    exp_dir = abspath(p[:experiment_path])
    db   = open_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM experiments WHERE path = ?", [exp_dir]))
    isempty(rows) && error("No experiment registered at $exp_dir. Run 'himalaya init' first.")
    exp_id        = Int(rows[1].id)
    exp           = get_experiment(db, exp_id)
    sample_filter = p[:sample]
    samples       = get_samples(db, exp_id)
    sample_filter !== nothing && filter!(sm -> sm.label == sample_filter, samples)

    for sample in samples
        exposures = get_exposures(db, Int(sample.id))

        # Auto-fallback: if no exposure is explicitly selected, use first accepted one
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
            # Mirror the auto-fallback's status guard above: rejected
            # exposures are explicitly out of the analysis set, so don't
            # waste compute or refresh peaks/indices for them.
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

function cli_show(args)
    s = ArgParseSettings(prog = "himalaya show")
    @add_arg_table! s begin
        "experiment_path"
            required = true
        "--sample", "-s"
            help     = "sample label"
            required = true
    end
    p = parse_args(args, s; as_symbols = true)

    exp_dir = abspath(p[:experiment_path])
    db   = open_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM experiments WHERE path = ?", [exp_dir]))
    isempty(rows) && error("No experiment registered at $exp_dir. Run 'himalaya init' first.")
    exp_id  = Int(rows[1].id)
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
the destination already exists (will not overwrite). This is the only
documented operation that writes to an experiment directory.
"""
function cli_config_new(; type_name::String = "simple", dir::String)
    isdir(dir) || error("Directory not found: $dir")
    dest = joinpath(dir, "experiment.toml")
    isfile(dest) && error("experiment.toml already exists at $dest — will not overwrite")
    src = joinpath(configs_dir(), type_name * ".toml")
    isfile(src) || error("Unknown config type '$type_name'. Run 'himalaya config list' to see options.")
    cp(src, dest)
    println("Created $dest from template '$type_name'")
    println("Edit it to set your experiment name, beamline parameters, and manifest column mappings.")
    dest
end

function cli_config(args)
    isempty(args) && (println("Usage: himalaya config <list|new> [options]"); return)
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
        "experiment_path"
            required = true
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

function main(args = copy(ARGS))
    isempty(args) && (println("Usage: himalaya <command> [args]"); return)
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
        println("Unknown command: $cmd. Available: init, analyze, show, serve, config, reingest")
    end
end
