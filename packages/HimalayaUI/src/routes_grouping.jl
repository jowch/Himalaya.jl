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

end  # register_grouping_routes!
