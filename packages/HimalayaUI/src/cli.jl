using ArgParse
using Printf

function cli_init(args)
    s = ArgParseSettings(prog = "himalaya init")
    @add_arg_table! s begin
        "experiment_path"
            help     = "path to experiment directory"
            required = true
        "--manifest", "-m"
            help     = "path to manifest CSV"
            default  = nothing
        "--beamline"
            help     = "beamline profile name (default: default)"
            default  = "default"
        "--name"
            help     = "experiment name"
            default  = nothing
    end
    p = parse_args(args, s; as_symbols = true)

    exp_path      = p[:experiment_path]
    manifest_path = p[:manifest]

    data_dir     = joinpath(exp_path, "data")
    analysis_dir = joinpath(exp_path, "analysis", "automatic_analysis")

    db     = open_db(exp_path)
    exp_id = init_experiment!(db;
        name          = something(p[:name], basename(exp_path)),
        path          = exp_path,
        data_dir      = data_dir,
        analysis_dir  = analysis_dir,
        manifest_path = manifest_path)

    if manifest_path !== nothing && isfile(manifest_path)
        samples = parse_manifest(manifest_path)
        for ms in samples
            s_id = create_sample!(db;
                experiment_id = exp_id,
                label         = ms.label,
                name          = ms.name,
                notes         = ms.notes_sample)
            for filename in ms.filenames
                create_exposure!(db; sample_id = s_id, filename = filename)
            end
            if !isempty(ms.notes_exposure)
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id FROM exposures WHERE sample_id = ? LIMIT 1", [s_id]))
                if !isempty(rows)
                    e_id = Int(rows[1].id)
                    DBInterface.execute(db,
                        "INSERT INTO exposure_tags (exposure_id, key, value, source)
                         VALUES (?, 'note', ?, 'manifest')",
                        [e_id, ms.notes_exposure])
                end
            end
        end
        println("Imported $(length(samples)) samples from manifest.")
    end

    println("Initialized experiment #$exp_id at $exp_path")
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

    db            = open_db(p[:experiment_path])
    exp           = get_experiment(db, 1)
    sample_filter = p[:sample]
    samples       = get_samples(db, 1)
    sample_filter !== nothing && filter!(sm -> sm.label == sample_filter, samples)

    for sample in samples
        exposures = get_exposures(db, Int(sample.id))
        for exp_row in exposures
            e_id = Int(exp_row.id)
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

    db      = open_db(p[:experiment_path])
    samples = get_samples(db, 1)
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

function cli_serve(args)
    s = ArgParseSettings(prog = "himalaya serve")
    @add_arg_table! s begin
        "experiment_path"
            required = true
        "--port"
            arg_type = Int
            default  = 8080
        "--host"
            default  = "127.0.0.1"
    end
    p = parse_args(args, s; as_symbols = true)

    db_path = joinpath(p[:experiment_path], "himalaya.db")
    isfile(db_path) || error("no himalaya.db at $db_path — run `himalaya init` first")

    db = open_db(p[:experiment_path])
    println("HimalayaUI serving $(p[:experiment_path]) on http://$(p[:host]):$(p[:port])")
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
    else
        println("Unknown command: $cmd. Available: init, analyze, show, serve")
    end
end
