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

    # Corpus-wide sample listing — every sample across all experiments, each
    # carrying its `tags` and a `q_units` sourced from the owning experiment's
    # config. Optional `?experiment_id=` narrows to one experiment, so this
    # route can stand in for the experiment-scoped `/api/experiments/{id}/samples`
    # (which remains in place).
    #
    # Three queries total, never N+1: samples, experiment configs, and one
    # batched tag query grouped in memory.
    @get "/api/samples" function(req::HTTP.Request)
        db     = current_db()
        params = HTTP.queryparams(HTTP.URI(req.target))

        # Optional ?experiment_id= filter. A non-integer value is a client
        # error, not something to silently ignore.
        exp_filter = nothing
        if haskey(params, "experiment_id")
            exp_filter = tryparse(Int, params["experiment_id"])
            exp_filter === nothing && return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "experiment_id must be an integer")))
        end

        samples = exp_filter === nothing ?
            Tables.rowtable(DBInterface.execute(db,
                "SELECT * FROM samples ORDER BY id")) :
            Tables.rowtable(DBInterface.execute(db,
                "SELECT * FROM samples WHERE experiment_id = ? ORDER BY id",
                [exp_filter]))

        # experiment_id -> q_units, one TOML parse per experiment (not per sample).
        qunits_by_exp = Dict{Int, String}()
        for er in Tables.rowtable(DBInterface.execute(db,
                "SELECT id, config FROM experiments"))
            qunits_by_exp[Int(er.id)] = _q_units_from_config(er.config)
        end

        # One batched tag query, grouped by sample_id. Skipped entirely when
        # there are no samples — an empty `IN ()` is invalid SQL.
        tags_by_sample = Dict{Int, Vector{Any}}()
        if !isempty(samples)
            ids          = [Int(sm.id) for sm in samples]
            placeholders = join(fill("?", length(ids)), ", ")
            tagrows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, sample_id, key, value, source FROM sample_tags
                 WHERE sample_id IN ($placeholders) ORDER BY id", ids))
            for tr in tagrows
                push!(get!(tags_by_sample, Int(tr.sample_id), Any[]),
                      Dict(:id     => Int(tr.id), :key   => tr.key,
                           :value  => tr.value,   :source => tr.source))
            end
        end

        out = map(samples) do sm
            d          = row_to_json(sm)
            d[:tags]   = get(tags_by_sample, Int(sm.id), Any[])
            # experiment_id is a nullable FK column; tolerate a SQL NULL
            # (e.g. a future ON DELETE SET NULL) rather than throwing in Int().
            d[:q_units] = ismissing(sm.experiment_id) ? "A-1" :
                get(qunits_by_exp, Int(sm.experiment_id), "A-1")
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
end
