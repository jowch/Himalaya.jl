using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_samples_routes!()
    @get "/api/experiments/{id}/samples" function(req::HTTP.Request, id::Int)
        db      = current_db()
        samples = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM samples WHERE experiment_id = ? ORDER BY id", [id]))
        out = map(samples) do sm
            tags = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, key, value, source FROM sample_tags
                 WHERE sample_id = ? ORDER BY id", [Int(sm.id)]))
            d = row_to_json(sm)
            d[:tags] = rows_to_json(tags)
            d
        end
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(out))
    end

    @patch "/api/samples/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        fields, vals = Symbol[], Any[]
        for k in (:name, :notes)
            if haskey(body, k)
                push!(fields, k); push!(vals, body[k])
            end
        end
        if isempty(fields)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "no updatable fields provided")))
        end
        sets = join(["$(string(f)) = ?" for f in fields], ", ")
        DBInterface.execute(db,
            "UPDATE samples SET $sets WHERE id = ?", vcat(vals, [id]))

        log_action!(db, req; action = "update_sample",
            entity_type = "sample", entity_id = id)

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM samples WHERE id = ?", [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(row_to_json(rows[1])))
    end

    @post "/api/samples/{id}/tags" function(req::HTTP.Request, id::Int)
        db    = current_db()
        body  = json(req)
        key   = String(body.key)
        value = String(body.value)
        res   = DBInterface.execute(db,
            "INSERT INTO sample_tags (sample_id, key, value, source)
             VALUES (?, ?, ?, 'manual')",
            [id, key, value])
        tag_id = Int(DBInterface.lastrowid(res))

        log_action!(db, req; action = "add_tag",
            entity_type = "sample", entity_id = id,
            note = "$key=$value")

        HTTP.Response(201, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => tag_id, :sample_id => id,
                             :key => key, :value => value, :source => "manual")))
    end

    @delete "/api/samples/{id}/tags/{tag_id}" function(req::HTTP.Request, id::Int, tag_id::Int)
        db = current_db()
        DBInterface.execute(db,
            "DELETE FROM sample_tags WHERE id = ? AND sample_id = ?",
            [tag_id, id])
        log_action!(db, req; action = "remove_tag",
            entity_type = "sample", entity_id = id,
            note = "tag_id=$tag_id")
        HTTP.Response(204)
    end
end
