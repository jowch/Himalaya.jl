using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_experiments_routes!()
    @get "/api/experiments" function(req::HTTP.Request)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM experiments ORDER BY id"))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(rows_to_json(rows)))
    end

    @get "/api/experiments/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(row_to_json(rows[1])))
    end

    @patch "/api/experiments/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)

        fields = Symbol[]
        vals   = Any[]
        for k in (:name, :data_dir, :analysis_dir, :manifest_path)
            if haskey(body, k)
                push!(fields, k)
                push!(vals, body[k])
            end
        end
        if isempty(fields)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "no updatable fields provided")))
        end

        sets = join(["$(string(f)) = ?" for f in fields], ", ")
        DBInterface.execute(db,
            "UPDATE experiments SET $sets WHERE id = ?", vcat(vals, [id]))

        log_action!(db, req; action = "update_experiment",
            entity_type = "experiment", entity_id = id)

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM experiments WHERE id = ?", [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(row_to_json(rows[1])))
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
                    push!(skipped, "$(sm.label)/$(ex.filename): $(sprint(showerror, e))")
                end
            end
        end

        log_action!(db, req; action = "analyze",
            entity_type = "experiment", entity_id = id,
            note = "analyzed=$analyzed skipped=$(length(skipped))")

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:analyzed => analyzed, :skipped => skipped)))
    end
end
