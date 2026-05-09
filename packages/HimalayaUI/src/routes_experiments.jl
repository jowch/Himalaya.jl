using HTTP, JSON3, DBInterface, Tables, Oxygen

function _experiment_row_to_json(row::NamedTuple)
    d = row_to_json(row)
    cfg_text = get(d, :config, nothing)
    # Defensive: malformed TOML in a single row must not 500 the list endpoint.
    # Fall back to the ASCII default; the UI prettifies it to "Å⁻¹".
    q_units = if cfg_text isa AbstractString && !isempty(cfg_text)
        try
            bl = get(TOML.parse(cfg_text), "beamline", Dict())
            get(bl, "q_units", "A-1")
        catch
            "A-1"
        end
    else
        "A-1"
    end
    d[:q_units] = q_units
    d
end

function register_experiments_routes!()
    @get "/api/experiments" function(req::HTTP.Request)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM experiments ORDER BY id"))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write([_experiment_row_to_json(r) for r in rows]))
    end

    @get "/api/experiments/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(_experiment_row_to_json(rows[1])))
    end

    @patch "/api/experiments/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)

        # Experiment name and path fields are no longer mutable via PATCH.
        # Name is derived from experiment.toml via reingest; path fields
        # (data_dir, analysis_dir, manifest_path) must also go through reingest
        # to stay in sync with the config blob.
        # This route is a defensive surface for future fields only.
        return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment metadata is read-only; rename via experiment.toml + reingest")))
    end

    @post "/api/experiments/{id}/analyze" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT analysis_dir FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))
        analysis_dir = rows[1].analysis_dir

        samples   = get_samples(db, id)
        analyzed  = 0
        skipped   = String[]
        for sm in samples
            for ex in get_exposures(db, Int(sm.id))
                try
                    analyze_exposure!(db, Int(ex.id), String(analysis_dir))
                    analyzed += 1
                catch e
                    push!(skipped, "$(sm.name)/$(ex.filename): $(sprint(showerror, e))")
                end
            end
        end

        log_action!(db, req; action = "analyze",
            entity_type = "experiment", entity_id = id,
            note = "analyzed=$analyzed skipped=$(length(skipped))")

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:analyzed => analyzed, :skipped => skipped)))
    end

    @post "/api/experiments/{id}/reingest" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT path FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))
        exp_path = String(rows[1].path)
        try
            res = reingest!(db, id, exp_path)
            log_action!(db, req; action = "reingest",
                entity_type = "experiment", entity_id = id)
            return HTTP.Response(200,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:status          => String(res.status),
                                 :added_samples   => res.added_samples,
                                 :added_exposures => res.added_exposures,
                                 :manifest_path   => res.manifest_path)))
        catch e
            if e isa ManifestValidationError
                return HTTP.Response(400,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "manifest_invalid",
                                     :violations => [Dict(:kind => string(v.kind),
                                                          :sample_index => v.sample_index,
                                                          :sample_name => v.sample_name,
                                                          :detail => v.detail)
                                                     for v in e.violations])))
            end
            return HTTP.Response(500,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => sprint(showerror, e))))
        end
    end
end
