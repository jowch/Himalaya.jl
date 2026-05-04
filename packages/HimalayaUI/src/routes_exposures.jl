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
            d[:tags]          = rows_to_json(tags)
            d[:sources]       = rows_to_json(srcs)
            d[:image_version] = image_version_token(e.image_path)
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

        # The frontend appends `?v=<image_version_token>` to the URL, so the
        # URL itself is the cache key. We can mark responses immutable and
        # cache them aggressively — when the underlying TIFF or our
        # processing code changes, the token (and therefore the URL) changes.
        vtoken = image_version_token(ip)
        HTTP.Response(200,
            ["Content-Type"    => "image/png",
             "Cache-Control"   => "private, max-age=31536000, immutable",
             "X-Image-Version" => vtoken],
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

        return with_idempotency(db, req) do
            DBInterface.execute(db,
                "UPDATE exposures SET status = ? WHERE id = ?", [status, id])

            result = apply_event!(InTransaction(), db, req;
                kind = "set_exposure_status",
                entity_type = "exposure", entity_id = id,
                payload = Dict(:status => status))
            _enqueue_broadcast_from_result!(result, "set_exposure_status", "exposure", id)

            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:id => id, :status => status)))
        end
    end

    # `selected` is sample-scoped client state — exactly one exposure per sample is
    # marked selected at a time. Under multiplayer this is intentionally LWW: if Alice
    # and Bob click different exposures within the same sample, the later click wins
    # and SSE broadcasts the resulting state. This route does NOT participate in
    # If-Match conflict resolution (see specs/2026-05-01-multiplayer-instrumentation-design.md).
    @patch "/api/exposures/{id}/select" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id FROM exposures WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure not found")))
        sample_id = Int(rows[1].sample_id)

        return with_idempotency(db, req) do
            DBInterface.execute(db,
                "UPDATE exposures SET selected = 0 WHERE sample_id = ?", [sample_id])
            DBInterface.execute(db,
                "UPDATE exposures SET selected = 1 WHERE id = ?", [id])

            result = apply_event!(InTransaction(), db, req;
                kind = "select_exposure",
                entity_type = "exposure", entity_id = id,
                payload = Dict(:sample_id => sample_id))
            _enqueue_broadcast_from_result!(result, "select_exposure", "exposure", id)

            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:id => id, :selected => true)))
        end
    end

    @post "/api/exposures/{id}/tags" function(req::HTTP.Request, id::Int)
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
                "INSERT INTO exposure_tags (exposure_id, key, value, source)
                 VALUES (?, ?, ?, 'manual')", [id, key, value])
            tag_id = Int(DBInterface.lastrowid(res))

            # Look up parent sample_id so frontend can invalidate the right cache key.
            erows = Tables.rowtable(DBInterface.execute(db,
                "SELECT sample_id FROM exposures WHERE id = ?", [id]))
            sample_id = isempty(erows) ? nothing : Int(erows[1].sample_id)

            result = apply_event!(InTransaction(), db, req;
                kind = "add_tag",
                entity_type = "exposure", entity_id = id,
                payload = Dict(:key => key, :value => value,
                               :tag_id => tag_id, :sample_id => sample_id))
            _enqueue_broadcast_from_result!(result, "add_tag", "exposure", id)

            # Response shape is `ExposureTag` exactly (id, key, value, source);
            # parent FK is implicit from the URL.
            HTTP.Response(201, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:id=>tag_id,
                                 :key=>key, :value=>value, :source=>"manual")))
        end
    end

    @delete "/api/exposures/{id}/tags/{tag_id}" function(req::HTTP.Request, id::Int, tag_id::Int)
        db = current_db()
        return with_idempotency(db, req) do
            # Query parent BEFORE deletion so the event payload carries sample_id.
            erows = Tables.rowtable(DBInterface.execute(db,
                "SELECT sample_id FROM exposures WHERE id = ?", [id]))
            sample_id = isempty(erows) ? nothing : Int(erows[1].sample_id)

            DBInterface.execute(db,
                "DELETE FROM exposure_tags WHERE id = ? AND exposure_id = ?",
                [tag_id, id])
            result = apply_event!(InTransaction(), db, req;
                kind = "remove_tag",
                entity_type = "exposure", entity_id = id,
                payload = Dict(:tag_id => tag_id, :sample_id => sample_id))
            _enqueue_broadcast_from_result!(result, "remove_tag", "exposure", id)
            HTTP.Response(204)
        end
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
        analysis_dir = String(rows[1].analysis_dir)

        # Issue #9: route is now queue-shaped — `with_idempotency` honors the
        # `X-Client-Op-Id` header so a network retry returns the cached
        # response rather than running analyze a second time.
        # `analyze_exposure!` emits its `analyze_run` event with
        # `defer_broadcast=true`; we re-fetch that row below and enqueue a
        # post-commit broadcast carrying the spec'd `post_state` envelope so
        # foreign tabs converge without a refetch round-trip (issue #12).
        return with_idempotency(db, req) do
            # Capture max user_actions id BEFORE analyze runs so we can
            # identify the analyze_run row this call emits (vs. a prior one
            # if analyze_exposure! takes the no-op fast path).
            pre_max_rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT IFNULL(MAX(id), 0) AS m FROM user_actions"))
            pre_max_id = Int(pre_max_rows[1].m)

            analyze_exposure!(db, id, analysis_dir; defer_broadcast = true)

            new_hash = read_inputs_hash(db, id)
            evt = Tables.rowtable(DBInterface.execute(db,
                """SELECT id, user_id, payload
                   FROM user_actions
                   WHERE action = 'analyze_run' AND entity_type = 'exposure'
                     AND entity_id = ? AND id > ?
                   ORDER BY id DESC LIMIT 1""", [id, pre_max_id]))
            if !isempty(evt)
                row = evt[1]
                post_state = Dict{Symbol, Any}(
                    :analysis_inputs_hash => new_hash,
                    :indices              => _serialized_indices_for_broadcast(db, id),
                )
                # The system-emitted row carries NULL client_id/op_id; pass
                # the request's headers in the SSE frame so foreign tabs can
                # still self-echo-filter (own client_id matches).
                _enqueue_post_commit_broadcast!(
                    Int(row.id), "analyze_run", "exposure", id,
                    ismissing(row.user_id) ? nothing : Int(row.user_id),
                    get_client_id(req),
                    get_client_op_id(req),
                    ismissing(row.payload) ? nothing : String(row.payload);
                    post_state = post_state)
            end

            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(
                    :id                   => id,
                    :analyzed             => true,
                    :analysis_inputs_hash => new_hash,
                )))
        end
    end

    @get "/api/exposures/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM exposures WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure not found")))
        ex   = rows[1]
        tags = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, key, value, source FROM exposure_tags
             WHERE exposure_id = ? ORDER BY id", [id]))
        srcs = Tables.rowtable(DBInterface.execute(db,
            "SELECT source_exposure_id, role FROM exposure_sources
             WHERE averaged_exposure_id = ?", [id]))
        d                = row_to_json(ex; bool_keys = (:selected,))
        d[:tags]         = rows_to_json(tags)
        d[:sources]      = rows_to_json(srcs)
        d[:image_version] = image_version_token(ex.image_path)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(d))
    end
end
