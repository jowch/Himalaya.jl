using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_trace_routes!()
    @get "/api/exposures/{id}/trace" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT e.filename, e.kind, x.id AS experiment_id, x.analysis_dir
             FROM exposures e JOIN samples s ON s.id = e.sample_id
             JOIN experiments x ON x.id = s.experiment_id
             WHERE e.id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure not found")))
        row = rows[1]
        String(row.kind) == "file" || return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "trace not available for derived exposures")))
        row.filename === missing && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure has no filename")))

        # Resolve via the experiment's configured integration_pattern (same as
        # analyze_exposure!). Hardcoding "{name}.dat" here breaks any experiment
        # whose pattern has trailing tokens like "{name}_tot.dat".
        cfg              = config_from_db(db, Int(row.experiment_id))
        pattern_filename = replace(cfg.integration_pattern, "{name}" => String(row.filename))
        path             = joinpath(String(row.analysis_dir), pattern_filename)
        isfile(path) || return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => ".dat file not found: $path")))

        q, I, σ = load_dat(path)
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:q => q, :I => I, :sigma => σ)))
    end
end
