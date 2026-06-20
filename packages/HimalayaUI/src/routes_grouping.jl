# packages/HimalayaUI/src/routes_grouping.jl
using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_grouping_routes!()

    # ── Rename a sample ──────────────────────────────────────────────────────
    # PATCH /api/samples/{id}/name
    # Body: { "name": "new label" }
    # Sets samples.name = body.name and name_source = 'user'. Records a
    # sample_renamed event for audit + SSE broadcast.
    # Clones the with_idempotency + apply_event! + _enqueue_broadcast pattern
    # from routes_samples.jl:167–185.
    @patch "/api/samples/{id}/name" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if !haskey(body, :name) || !(body.name isa AbstractString)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "name is required and must be a string")))
        end
        new_name = strip(String(body.name))
        isempty(new_name) && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "name must not be blank")))

        # Resolve the experiment for the SSE payload (frontend invalidates by it).
        exp_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT experiment_id FROM samples WHERE id = ?", [id]))
        isempty(exp_rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "sample not found")))
        exp_id = Int(exp_rows[1].experiment_id)

        return with_idempotency(db, req) do
            DBInterface.execute(db,
                "UPDATE samples SET name = ?, name_source = 'user' WHERE id = ?",
                [new_name, id])
            result = apply_event!(InTransaction(), db, req;
                kind        = "sample_renamed",
                entity_type = "sample",
                entity_id   = id,
                payload     = Dict(:name => new_name, :experiment_id => exp_id))
            _enqueue_broadcast_from_result!(result, "sample_renamed", "sample", id)

            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT * FROM samples WHERE id = ?", [id]))
            isempty(rows) && return HTTP.Response(404,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "sample not found")))
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(row_to_json(rows[1])))
        end
    end

    # ── Move one exposure to a different sample ───────────────────────────────
    # POST /api/exposures/{id}/move
    # Body: { "sample_id": Int }
    # Within-experiment only (spec §9.2 / §5).
    @post "/api/exposures/{id}/move" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if !haskey(body, :sample_id)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "sample_id is required")))
        end
        dest_sample_id = Int(body.sample_id)

        # Validate same-experiment constraint; also capture the source sample id
        # and the experiment id for the SSE payload.
        src_rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT e.sample_id AS from_sample_id,
                      e.experiment_id AS src_exp,
                      s.experiment_id AS dst_exp
               FROM exposures e
               JOIN samples s ON s.id = ?
               WHERE e.id = ?""", [dest_sample_id, id]))
        if isempty(src_rows)
            return HTTP.Response(404,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "exposure or destination sample not found")))
        end
        src_row = first(src_rows)
        if !ismissing(src_row.src_exp) && !ismissing(src_row.dst_exp) &&
                Int(src_row.src_exp) != Int(src_row.dst_exp)
            return HTTP.Response(422,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error =>
                    "cross-experiment moves are not allowed; move would change the resolved analysis_dir")))
        end
        from_sample_id = ismissing(src_row.from_sample_id) ? nothing : Int(src_row.from_sample_id)
        exp_id = ismissing(src_row.dst_exp) ? nothing : Int(src_row.dst_exp)

        return with_idempotency(db, req) do
            result = apply_event!(InTransaction(), db, req;
                kind        = "exposure_moved",
                entity_type = "exposure",
                entity_id   = id,
                payload     = Dict(:sample_id => dest_sample_id,
                                   :from_sample_id => from_sample_id,
                                   :experiment_id => exp_id))
            _enqueue_broadcast_from_result!(result, "exposure_moved", "exposure", id)

            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT * FROM exposures WHERE id = ?", [id]))
            isempty(rows) && return HTTP.Response(404,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "exposure not found")))
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(row_to_json(rows[1])))
        end
    end

    # ── Merge samples: loser → survivor ──────────────────────────────────────
    # POST /api/samples/{id}/merge
    # Body: { "survivor_id": Int, "name"?: String }
    # id = loser. All child rows re-pointed, then loser retired.
    # Whole body is one with_idempotency block (spec §9.3).
    # Undo is session-local (frontend useUndoStack). No server-side undo API.
    # SSE: applyRemoteToCache must refetch loads(id) on merge events (invalidate-only).
    @post "/api/samples/{id}/merge" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if !haskey(body, :survivor_id)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "survivor_id is required")))
        end
        loser_id    = id
        survivor_id = Int(body.survivor_id)
        loser_id == survivor_id && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "loser and survivor must be different")))

        # Validate both samples exist and belong to the same experiment.
        both = Tables.rowtable(DBInterface.execute(db,
            """SELECT id, experiment_id FROM samples WHERE id IN (?, ?)""",
            [loser_id, survivor_id]))
        length(both) != 2 && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "one or both samples not found")))
        if Int(both[1].experiment_id) != Int(both[2].experiment_id)
            return HTTP.Response(422,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "cannot merge samples from different experiments")))
        end
        exp_id = Int(both[1].experiment_id)

        return with_idempotency(db, req) do
            # (i) series_samples dedup — drop loser's membership in any series
            # where the survivor already appears (spec §9.3). Done BEFORE re-point
            # so the surviving UNIQUE(series_id, position) isn't invalidated first.
            DBInterface.execute(db,
                """DELETE FROM series_samples
                   WHERE sample_id = ?
                     AND series_id IN (
                       SELECT series_id FROM series_samples WHERE sample_id = ?
                     )""",
                [loser_id, survivor_id])

            # (ii) sample_tags dedup — drop loser's tag where survivor holds same key;
            # survivor wins on collision (spec §9.3).
            DBInterface.execute(db,
                """DELETE FROM sample_tags
                   WHERE sample_id = ?
                     AND key IN (
                       SELECT key FROM sample_tags WHERE sample_id = ?
                     )""",
                [loser_id, survivor_id])

            # (iii) re-point all child FK tables to survivor.
            # exposures: emit one exposure_moved event per exposure moved.
            # Each per-exposure call is distinct in the idempotency index via its
            # differing entity_id under the (client_op_id, action, entity_id)
            # partial UNIQUE — the outer with_idempotency owns retry-replay.
            loser_exposures = Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM exposures WHERE sample_id = ?", [loser_id]))
            for ex_row in loser_exposures
                ex_id = Int(ex_row.id)
                result = apply_event!(InTransaction(), db, req;
                    kind        = "exposure_moved",
                    entity_type = "exposure",
                    entity_id   = ex_id,
                    payload     = Dict(:sample_id => survivor_id,
                                       :from_sample_id => loser_id,
                                       :experiment_id => exp_id))
                _enqueue_broadcast_from_result!(result, "exposure_moved", "exposure", ex_id)
            end

            # series_samples re-point (survivors of the dedup step).
            DBInterface.execute(db,
                "UPDATE series_samples SET sample_id = ? WHERE sample_id = ?",
                [survivor_id, loser_id])

            # sample_tags re-point (survivors of the dedup step).
            DBInterface.execute(db,
                "UPDATE sample_tags SET sample_id = ? WHERE sample_id = ?",
                [survivor_id, loser_id])

            # sample_messages re-point (no uniqueness constraint — plain UPDATE).
            DBInterface.execute(db,
                "UPDATE sample_messages SET sample_id = ? WHERE sample_id = ?",
                [survivor_id, loser_id])

            # (iv) Retire the loser (sets merged_into_id; does not hard-delete).
            retire_sample!(db, loser_id; merged_into_id=survivor_id)

            # (v) Optional rename of survivor.
            if haskey(body, :name) && body.name isa AbstractString
                new_name = strip(String(body.name))
                if !isempty(new_name)
                    DBInterface.execute(db,
                        "UPDATE samples SET name = ?, name_source = 'user' WHERE id = ?",
                        [new_name, survivor_id])
                    rename_result = apply_event!(InTransaction(), db, req;
                        kind        = "sample_renamed",
                        entity_type = "sample",
                        entity_id   = survivor_id,
                        payload     = Dict(:name => new_name, :experiment_id => exp_id))
                    _enqueue_broadcast_from_result!(rename_result, "sample_renamed", "sample", survivor_id)
                end
            end

            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:loser_id => loser_id, :survivor_id => survivor_id)))
        end
    end

end  # register_grouping_routes!
