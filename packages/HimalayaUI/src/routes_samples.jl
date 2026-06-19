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

        # Per-sample indexing rollup (the R1/#224 contact-sheet seam). Each
        # sample's status comes from its REPRESENTATIVE exposure — highest-id
        # `selected=1`, else highest-id — the same exposure the picker/series
        # index against. Exposures of one sample legitimately disagree (the beam
        # may miss the ordered phase on some frames; a miss reads like form
        # factor), so the representative is a faithful one-line summary rather
        # than a lossy aggregate. We surface that exposure's durable assignment:
        # its dominant (top-scored) confirmed phase when `indexed`, the
        # `form_factor` declaration distinctly, else unindexed. Three bulk
        # queries, no per-sample N+1 (mirrors picker_samples).
        rep_by_sample = Dict{Int, Int}()
        state_by_exp  = Dict{Int, String}()
        phase_by_exp  = Dict{Int, String}()
        if !isempty(samples)
            sids = [Int(sm.id) for sm in samples]
            ph   = join(fill("?", length(sids)), ",")
            grouped_exps = Dict{Int, Vector{NamedTuple}}()
            for e in Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, sample_id, selected FROM exposures
                     WHERE sample_id IN ($ph) ORDER BY sample_id ASC, id ASC", sids))
                push!(get!(grouped_exps, Int(e.sample_id), NamedTuple[]), e)
            end
            for (sid, exps) in grouped_exps
                rep = nothing
                for e in Iterators.reverse(exps)   # highest-id selected wins
                    if e.selected != 0
                        rep = Int(e.id); break
                    end
                end
                rep === nothing && (rep = Int(last(exps).id))   # else highest-id
                rep_by_sample[sid] = rep
            end

            rep_ids = collect(values(rep_by_sample))
            if !isempty(rep_ids)
                rp = join(fill("?", length(rep_ids)), ",")
                for r in Tables.rowtable(DBInterface.execute(db,
                        "SELECT exposure_id, state FROM assignments
                         WHERE exposure_id IN ($rp)", rep_ids))
                    state_by_exp[Int(r.exposure_id)] = String(r.state)
                end
                # Dominant phase: the top-scored assigned index's phase per
                # exposure (mirrors compute_member_snapshot's confirmed_index).
                best_score = Dict{Int, Float64}()
                for r in Tables.rowtable(DBInterface.execute(db,
                        "SELECT m.exposure_id AS eid, i.phase AS phase, i.score AS score
                         FROM assignment_members m JOIN indices i ON i.id = m.index_id
                         WHERE m.exposure_id IN ($rp)", rep_ids))
                    eid = Int(r.eid)
                    sc  = r.score === missing ? -Inf : Float64(r.score)
                    if !haskey(best_score, eid) || sc > best_score[eid]
                        best_score[eid]  = sc
                        phase_by_exp[eid] = String(r.phase)
                    end
                end
            end
        end

        out = map(samples) do sm
            d          = row_to_json(sm)
            d[:tags]   = get(tags_by_sample, Int(sm.id), Any[])
            # experiment_id is a nullable FK column; tolerate a SQL NULL
            # (e.g. a future ON DELETE SET NULL) rather than throwing in Int().
            d[:q_units] = ismissing(sm.experiment_id) ? "A-1" :
                get(qunits_by_exp, Int(sm.experiment_id), "A-1")
            rep = get(rep_by_sample, Int(sm.id), nothing)
            if rep !== nothing
                state = get(state_by_exp, rep, "indexed")
                d[:assignment_state] = state
                # Phase only when actually indexed; form_factor/null carry none.
                d[:phase] = state == "indexed" ? get(phase_by_exp, rep, nothing) : nothing
            end
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
                v = body[k]
                # Trim leading/trailing whitespace on name only.
                v isa AbstractString && k === :name && (v = strip(String(v)))
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
        # `source` is optional and defaults to 'manual'. Validated BEFORE the
        # transaction opens so a malformed request is a clean 400 (mirrors the
        # batch route's validation discipline). The series-scoping path passes
        # an explicit non-'manual' source.
        if haskey(body, :source) && !(body.source isa AbstractString)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "source must be a string")))
        end
        key    = String(body.key)
        value  = String(body.value)
        source = haskey(body, :source) ? String(body.source) : "manual"
        return with_idempotency(db, req) do
            # Single-valued-key invariant (TAG-DEDUP-MODEL): a sample carries at
            # most one tag per key, enforced by the sample_tags_unique_key index.
            # Guard the collision here with a clean 409 — a bare INSERT would hit
            # the index and surface as an unhandled 500. Mirrors the PATCH route's
            # key-collision handling; the frontend TagEditor validates client-side,
            # so this is the server-side backstop.
            coll = Tables.rowtable(DBInterface.execute(db,
                "SELECT 1 FROM sample_tags WHERE sample_id = ? AND key = ?",
                [id, key]))
            isempty(coll) || return HTTP.Response(409,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "sample already has a '$key' tag")))

            res   = DBInterface.execute(db,
                "INSERT INTO sample_tags (sample_id, key, value, source)
                 VALUES (?, ?, ?, ?)",
                [id, key, value, source])
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
                                 :key => key, :value => value, :source => source)))
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

    @patch "/api/samples/{id}/tags/{tag_id}" function(req::HTTP.Request, id::Int, tag_id::Int)
        db   = current_db()
        body = json(req)
        # Partial update: key and/or value, each a string if present.
        has_key = haskey(body, :key); has_val = haskey(body, :value)
        (has_key || has_val) ||
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "provide at least one of: key, value")))
        if (has_key && !(body.key isa AbstractString)) ||
           (has_val && !(body.value isa AbstractString))
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "key and value must be strings")))
        end
        return with_idempotency(db, req) do
            cur = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, key, value, source FROM sample_tags
                 WHERE id = ? AND sample_id = ?", [tag_id, id]))
            isempty(cur) && return HTTP.Response(404, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "tag not found")))
            cur_key = String(cur[1].key)
            new_key = has_key ? String(body.key) : cur_key
            new_val = has_val ? String(body.value) : String(cur[1].value)
            # Single-valued-key rule: only a key *change* can introduce a
            # collision. A value-only edit (or a no-op key write — the client
            # sends the unchanged key alongside a value edit) leaves the key
            # set untouched, so it can't violate the invariant and must not be
            # rejected — including on a sample that already carries a same-key
            # duplicate, which is exactly the state this modal exists to fix.
            if new_key != cur_key
                coll = Tables.rowtable(DBInterface.execute(db,
                    "SELECT 1 FROM sample_tags
                     WHERE sample_id = ? AND key = ? AND id <> ?",
                    [id, new_key, tag_id]))
                isempty(coll) || return HTTP.Response(409, ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "sample already has a '$new_key' tag")))
            end
            srows = Tables.rowtable(DBInterface.execute(db,
                "SELECT experiment_id FROM samples WHERE id = ?", [id]))
            exp_id = (isempty(srows) || ismissing(srows[1].experiment_id)) ?
                     nothing : Int(srows[1].experiment_id)
            DBInterface.execute(db,
                "UPDATE sample_tags SET key = ?, value = ? WHERE id = ? AND sample_id = ?",
                [new_key, new_val, tag_id, id])
            result = apply_event!(InTransaction(), db, req;
                kind = "edit_tag",
                entity_type = "sample", entity_id = id,
                payload = Dict(:tag_id => tag_id, :key => new_key, :value => new_val,
                               :experiment_id => exp_id))
            _enqueue_broadcast_from_result!(result, "edit_tag", "sample", id)
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:id => tag_id, :key => new_key, :value => new_val,
                                 :source => String(cur[1].source))))
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

        # All shape + type validation runs BEFORE the transaction opens, so a
        # malformed request is a clean 400 (never an unguarded-conversion 500)
        # with nothing written — test_route_validation_routing.jl pins this.
        _bad(msg) = HTTP.Response(400, ["Content-Type" => "application/json"],
                                  JSON3.write(Dict(:error => msg)))

        body isa AbstractDict ||
            return _bad("request body must be a JSON object")
        (haskey(body, :key) && body.key isa AbstractString) ||
            return _bad("missing or non-string required field: key")
        (haskey(body, :tags) && body.tags isa AbstractVector && !isempty(body.tags)) ||
            return _bad("missing or empty required field: tags")
        if haskey(body, :source) && !(body.source isa AbstractString)
            return _bad("source must be a string")
        end

        # Each entry must be a {sample_id::Int, value::String} object. A
        # duplicate sample_id is rejected here: every event the batch emits
        # shares one client_op_id, so two entries with the same sample_id
        # would collide on idx_events_unique_op (client_op_id, action,
        # entity_id) — the second add_tag INSERT is swallowed as an
        # idempotent retry, leaving a tag row with no durable event.
        seen = Set{Int}()
        for t in body.tags
            (t isa AbstractDict && haskey(t, :sample_id) && haskey(t, :value)) ||
                return _bad("each tag requires sample_id and value")
            sid = t.sample_id
            (sid isa Integer || (sid isa AbstractFloat && isinteger(sid))) ||
                return _bad("each tag sample_id must be an integer")
            t.value isa AbstractString ||
                return _bad("each tag value must be a string")
            sid_int = Int(sid)
            sid_int in seen && return _bad("duplicate sample_id in batch: $sid_int")
            push!(seen, sid_int)
        end

        key     = String(body.key)
        source  = haskey(body, :source) ? String(body.source) : "manual"
        entries = body.tags

        return with_idempotency(db, req) do
            # `with_idempotency` supplies the single enclosing transaction —
            # every insert/update + event below commits or rolls back together.
            created = Vector{Dict{Symbol, Any}}()
            for t in entries
                sample_id = Int(t.sample_id)
                value     = String(t.value)

                # Parent experiment_id lets the frontend invalidate the right
                # samples cache key — matches the single-tag route's payload.
                # Hoisted to the top: all three branches (insert/update/no-op)
                # need it for the event payload.
                srows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT experiment_id FROM samples WHERE id = ?", [sample_id]))
                exp_id = (isempty(srows) || ismissing(srows[1].experiment_id)) ?
                         nothing : Int(srows[1].experiment_id)

                # Upsert on (sample_id, key): update the oldest matching row if
                # the key already exists, insert only when absent. Ordering by id
                # ensures a stable winner when duplicates are already present
                # (LO-TAGDUP recovery path) and the SELECT covers source so we
                # can preserve the existing tag's provenance on update.
                existing = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, value, source FROM sample_tags
                     WHERE sample_id = ? AND key = ? ORDER BY id LIMIT 1",
                    [sample_id, key]))

                if isempty(existing)
                    # --- INSERT branch: key absent → add new tag ---
                    res = DBInterface.execute(db,
                        "INSERT INTO sample_tags (sample_id, key, value, source)
                         VALUES (?, ?, ?, ?)",
                        [sample_id, key, value, source])
                    tag_id       = Int(DBInterface.lastrowid(res))
                    actual_source = source   # new row carries the batch source
                    result = apply_event!(InTransaction(), db, req;
                        kind = "add_tag",
                        entity_type = "sample", entity_id = sample_id,
                        payload = Dict(:key => key, :value => value,
                                       :tag_id => tag_id, :experiment_id => exp_id))
                    _enqueue_broadcast_from_result!(result, "add_tag",
                                                    "sample", sample_id)
                elseif String(existing[1].value) != value
                    # --- UPDATE branch: key present but value changed ---
                    # Do NOT touch `source`: scoping must not overwrite a manual
                    # tag's provenance (the whole point of the upsert semantics).
                    tag_id        = Int(existing[1].id)
                    actual_source = String(existing[1].source)
                    DBInterface.execute(db,
                        "UPDATE sample_tags SET value = ? WHERE id = ?",
                        [value, tag_id])
                    result = apply_event!(InTransaction(), db, req;
                        kind = "edit_tag",
                        entity_type = "sample", entity_id = sample_id,
                        payload = Dict(:tag_id => tag_id, :key => key,
                                       :value => value, :experiment_id => exp_id))
                    _enqueue_broadcast_from_result!(result, "edit_tag",
                                                    "sample", sample_id)
                else
                    # --- NO-OP branch: key present, value unchanged ---
                    # No write, no event — the batch response still carries one
                    # entry per input row (contract preserved below).
                    tag_id        = Int(existing[1].id)
                    actual_source = String(existing[1].source)
                end

                # Each entry is a full SampleTag plus the explicit sample_id —
                # there is no URL param to make the parent FK implicit here.
                # One shared push for all three branches: uses the resolved
                # tag_id, the new value (or unchanged value on no-op), and the
                # ACTUAL stored source (batch source on insert, existing tag's
                # source on update/no-op).
                push!(created, Dict(:id => tag_id, :sample_id => sample_id,
                                    :key => key, :value => value,
                                    :source => actual_source))
            end
            HTTP.Response(201, ["Content-Type" => "application/json"],
                JSON3.write(created))
        end
    end
end
