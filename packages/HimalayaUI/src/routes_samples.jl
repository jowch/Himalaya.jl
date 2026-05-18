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
        for k in (:display_name, :notes)
            if haskey(body, k)
                v = body[k]
                # Trim leading/trailing whitespace on display_name only.
                v isa AbstractString && k === :display_name && (v = strip(String(v)))
                push!(fields, k); push!(vals, v)
            end
        end
        if isempty(fields)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "no updatable fields provided")))
        end
        return with_idempotency(db, req) do
            sets = join(["$(string(f)) = ?" for f in fields], ", ")
            DBInterface.execute(db,
                "UPDATE samples SET $sets WHERE id = ?", vcat(vals, [id]))

            # Structured payload: the patched fields directly. Frontend
            # applyRemoteToCache spreads this onto the cached sample.
            update_payload = Dict{Symbol, Any}(zip(fields, vals))
            result = apply_event!(InTransaction(), db, req;
                kind = "update_sample",
                entity_type = "sample", entity_id = id,
                payload = update_payload)
            _enqueue_broadcast_from_result!(result, "update_sample", "sample", id)

            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT * FROM samples WHERE id = ?", [id]))
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(row_to_json(rows[1])))
        end
    end

    @post "/api/samples/{id}/tags" function(req::HTTP.Request, id::Int)
        db    = current_db()
        body  = json(req)
        if !haskey(body, :key) || !haskey(body, :value)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "missing required fields: key, value")))
        end
        key   = String(body.key)
        value = String(body.value)
        return with_idempotency(db, req) do
            res   = DBInterface.execute(db,
                "INSERT INTO sample_tags (sample_id, key, value, source)
                 VALUES (?, ?, ?, 'manual')",
                [id, key, value])
            tag_id = Int(DBInterface.lastrowid(res))

            # Look up parent experiment_id so the frontend can invalidate the
            # right samples cache key.
            srows = Tables.rowtable(DBInterface.execute(db,
                "SELECT experiment_id FROM samples WHERE id = ?", [id]))
            exp_id = isempty(srows) ? nothing : Int(srows[1].experiment_id)

            result = apply_event!(InTransaction(), db, req;
                kind = "add_tag",
                entity_type = "sample", entity_id = id,
                payload = Dict(:key => key, :value => value,
                               :tag_id => tag_id, :experiment_id => exp_id))
            _enqueue_broadcast_from_result!(result, "add_tag", "sample", id)

            # Response shape is `SampleTag` exactly (id, key, value, source).
            # The parent FK is implicit from the URL — including it would be
            # cache pollution if the frontend spread the response wholesale.
            HTTP.Response(201, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:id => tag_id,
                                 :key => key, :value => value, :source => "manual")))
        end
    end

    @delete "/api/samples/{id}/tags/{tag_id}" function(req::HTTP.Request, id::Int, tag_id::Int)
        db = current_db()
        return with_idempotency(db, req) do
            # Query parent BEFORE deletion so we can include experiment_id in the
            # event payload for cache invalidation.
            srows = Tables.rowtable(DBInterface.execute(db,
                "SELECT experiment_id FROM samples WHERE id = ?", [id]))
            exp_id = isempty(srows) ? nothing : Int(srows[1].experiment_id)

            DBInterface.execute(db,
                "DELETE FROM sample_tags WHERE id = ? AND sample_id = ?",
                [tag_id, id])
            result = apply_event!(InTransaction(), db, req;
                kind = "remove_tag",
                entity_type = "sample", entity_id = id,
                payload = Dict(:tag_id => tag_id, :experiment_id => exp_id))
            _enqueue_broadcast_from_result!(result, "remove_tag", "sample", id)
            HTTP.Response(204)
        end
    end

    @get "/api/samples/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM samples WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "sample not found")))
        sm   = rows[1]
        tags = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, key, value, source FROM sample_tags
             WHERE sample_id = ? ORDER BY id", [id]))
        d        = row_to_json(sm)
        d[:tags] = rows_to_json(tags)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(d))
    end

    # Batch sample-tag insert (I2.6). One `with_idempotency` transaction, N
    # `sample_tags` rows + N `add_tag` events — the atomic boundary the
    # series scoping step needs (N single `add_tag` calls have none; a
    # mid-batch reload would leave a half-confirmed recipe). Reuses the
    # existing `add_tag` event kind — `source` defaults to 'manual', the
    # scoping step passes 'scoping'.
    @post "/api/samples/tags/batch" function(req::HTTP.Request)
        db   = current_db()
        body = json(req)

        # Validate request shape BEFORE opening the transaction, so a
        # malformed request is a clean 400 with nothing written.
        if !haskey(body, :key)
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "missing required field: key")))
        end
        if !haskey(body, :tags) || isempty(body.tags)
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "missing or empty required field: tags")))
        end
        for t in body.tags
            if !haskey(t, :sample_id) || !haskey(t, :value)
                return HTTP.Response(400,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error =>
                        "each tag requires sample_id and value")))
            end
        end

        key     = String(body.key)
        source  = haskey(body, :source) ? String(body.source) : "manual"
        entries = body.tags

        return with_idempotency(db, req) do
            # `with_idempotency` supplies the single enclosing transaction —
            # every insert + event below commits or rolls back together.
            created = Vector{Dict{Symbol, Any}}()
            for t in entries
                sample_id = Int(t.sample_id)
                value     = String(t.value)
                res = DBInterface.execute(db,
                    "INSERT INTO sample_tags (sample_id, key, value, source)
                     VALUES (?, ?, ?, ?)",
                    [sample_id, key, value, source])
                tag_id = Int(DBInterface.lastrowid(res))

                # Parent experiment_id lets the frontend invalidate the right
                # samples cache key — matches the single-tag route's payload.
                srows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT experiment_id FROM samples WHERE id = ?", [sample_id]))
                exp_id = isempty(srows) ? nothing : Int(srows[1].experiment_id)

                result = apply_event!(InTransaction(), db, req;
                    kind = "add_tag",
                    entity_type = "sample", entity_id = sample_id,
                    payload = Dict(:key => key, :value => value,
                                   :tag_id => tag_id, :experiment_id => exp_id))
                _enqueue_broadcast_from_result!(result, "add_tag",
                                                "sample", sample_id)

                # Each entry is a full SampleTag plus the explicit sample_id —
                # there is no URL param to make the parent FK implicit here.
                push!(created, Dict(:id => tag_id, :sample_id => sample_id,
                                    :key => key, :value => value,
                                    :source => source))
            end
            HTTP.Response(201, ["Content-Type" => "application/json"],
                JSON3.write(created))
        end
    end
end
