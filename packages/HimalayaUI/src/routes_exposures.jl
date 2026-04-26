using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_exposures_routes!()
    @get "/api/samples/{id}/exposures" function(req::HTTP.Request, id::Int)
        db     = current_db()
        params = HTTP.queryparams(req)
        exclude_rejected = get(params, "exclude_rejected", "false") == "true"

        sql = exclude_rejected ?
            "SELECT * FROM exposures WHERE sample_id = ? AND (status IS NULL OR status != 'rejected') ORDER BY id" :
            "SELECT * FROM exposures WHERE sample_id = ? ORDER BY id"

        exs = Tables.rowtable(DBInterface.execute(db, sql, [id]))
        out = map(exs) do e
            tags = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, key, value, source FROM exposure_tags
                 WHERE exposure_id = ? ORDER BY id", [Int(e.id)]))
            srcs = Tables.rowtable(DBInterface.execute(db,
                "SELECT source_exposure_id, role FROM exposure_sources
                 WHERE averaged_exposure_id = ?", [Int(e.id)]))
            d = row_to_json(e; bool_keys = (:selected,))
            d[:tags]    = rows_to_json(tags)
            d[:sources] = rows_to_json(srcs)
            d
        end
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(out))
    end

    @get "/api/exposures/{id}/image" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT image_path FROM exposures WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure not found")))

        ip = rows[1].image_path
        (ip === nothing || ip isa Missing) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "no image for this exposure")))

        params   = HTTP.queryparams(req)
        is_thumb = get(params, "thumb", "0") == "1"

        img = load_and_lognormalize(String(ip))
        if is_thumb
            img = resize_to_fit(img, 128)
        end
        bytes = encode_png(img)

        HTTP.Response(200,
            ["Content-Type"  => "image/png",
             "Cache-Control" => "max-age=3600"],
            bytes)
    end

    @patch "/api/exposures/{id}/status" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)

        raw_status = get(body, :status, nothing)
        status = raw_status === nothing ? nothing : String(raw_status)

        if status !== nothing && status ∉ ("accepted", "rejected")
            return HTTP.Response(422,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "status must be 'accepted', 'rejected', or null")))
        end

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM exposures WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure not found")))

        DBInterface.execute(db,
            "UPDATE exposures SET status = ? WHERE id = ?", [status, id])

        log_action!(db, req; action = "set_status",
            entity_type = "exposure", entity_id = id,
            note = something(status, "null"))

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => id, :status => status)))
    end

    @patch "/api/exposures/{id}/select" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id FROM exposures WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure not found")))
        sample_id = Int(rows[1].sample_id)

        DBInterface.execute(db,
            "UPDATE exposures SET selected = 0 WHERE sample_id = ?", [sample_id])
        DBInterface.execute(db,
            "UPDATE exposures SET selected = 1 WHERE id = ?", [id])

        log_action!(db, req; action = "select_exposure",
            entity_type = "exposure", entity_id = id)

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => id, :selected => true)))
    end

    @post "/api/exposures/{id}/tags" function(req::HTTP.Request, id::Int)
        db    = current_db()
        body  = json(req)
        key   = String(body.key)
        value = String(body.value)
        res   = DBInterface.execute(db,
            "INSERT INTO exposure_tags (exposure_id, key, value, source)
             VALUES (?, ?, ?, 'manual')", [id, key, value])
        tag_id = Int(DBInterface.lastrowid(res))

        log_action!(db, req; action = "add_tag",
            entity_type = "exposure", entity_id = id, note = "$key=$value")

        HTTP.Response(201, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id=>tag_id, :exposure_id=>id,
                             :key=>key, :value=>value, :source=>"manual")))
    end

    @delete "/api/exposures/{id}/tags/{tag_id}" function(req::HTTP.Request, id::Int, tag_id::Int)
        db = current_db()
        DBInterface.execute(db,
            "DELETE FROM exposure_tags WHERE id = ? AND exposure_id = ?",
            [tag_id, id])
        log_action!(db, req; action = "remove_tag",
            entity_type = "exposure", entity_id = id, note = "tag_id=$tag_id")
        HTTP.Response(204)
    end

    @post "/api/exposures/{id}/analyze" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT e.id, x.analysis_dir
             FROM exposures e JOIN samples s ON s.id = e.sample_id
             JOIN experiments x ON x.id = s.experiment_id
             WHERE e.id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure not found")))

        analyze_exposure!(db, id, String(rows[1].analysis_dir))

        log_action!(db, req; action = "analyze",
            entity_type = "exposure", entity_id = id)

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => id, :analyzed => true)))
    end
end
