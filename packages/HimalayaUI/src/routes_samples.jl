using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_samples_routes!()
    @get "/api/samples" function(req::HTTP.Request)
        db     = current_db()
        params = HTTP.queryparams(HTTP.URI(req.target))

        exp_filter = nothing
        if haskey(params, "experiment_id")
            exp_filter = tryparse(Int, params["experiment_id"])
            exp_filter === nothing && return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "experiment_id must be an integer")))
        end

        samples = exp_filter === nothing ?
            sqlrows(DBInterface.execute(db,
                "SELECT * FROM samples ORDER BY id")) :
            sqlrows(DBInterface.execute(db,
                "SELECT * FROM samples WHERE experiment_id = ? ORDER BY id",
                [exp_filter]))

        qunits_by_setup = Dict{Int, String}()
        for sr in sqlrows(DBInterface.execute(db,
                "SELECT id, q_units FROM setups"))
            qunits_by_setup[Int(sr.id)] = ismissing(sr.q_units) ? "A-1" : String(sr.q_units)
        end

        rep_by_sample = Dict{Int, Int}()
        sids = [Int(sm.id) for sm in samples]
        if !isempty(samples)
            ph   = join(fill("?", length(sids)), ",")
            grouped_exps = Dict{Int, Vector{NamedTuple}}()
            for e in sqlrows(DBInterface.execute(db,
                    "SELECT id, sample_id, status FROM exposures
                     WHERE sample_id IN ($ph) ORDER BY sample_id ASC, id ASC", sids))
                push!(get!(grouped_exps, Int(e.sample_id), NamedTuple[]), e)
            end
            for (sid, exps) in grouped_exps
                rep = representative_exposure_id(exps)
                rep === nothing || (rep_by_sample[sid] = rep)
            end
        end
        call_states  = sample_call_states(db, sids)
        call_members = sample_call_members(db, sids)

        out = map(samples) do sm
            d          = row_to_json(sm)
            d[:q_units] = ismissing(sm.setup_id) ? "A-1" :
                get(qunits_by_setup, Int(sm.setup_id), "A-1")
            d[:representative_exposure_id] = get(rep_by_sample, Int(sm.id), nothing)
            state = get(call_states, Int(sm.id), "indexed")
            summary = call_summary(state, get(call_members, Int(sm.id), []))
            d[:assignment_state] = state
            d[:phase] = summary.phase
            d[:call_phases] = summary.call_phases
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
                v isa AbstractString && k === :name && (v = strip(String(v)))
                push!(fields, k); push!(vals, v)
            end
        end
        if :name in fields
            push!(fields, :name_source); push!(vals, "user")
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

            update_payload = Dict{Symbol, Any}(zip(fields, vals))
            result = apply_event!(InTransaction(), db, req;
                kind = "update_sample",
                entity_type = "sample", entity_id = id,
                payload = update_payload)
            _enqueue_broadcast_from_result!(result, "update_sample", "sample", id)

            rows = sqlrows(DBInterface.execute(db,
                "SELECT * FROM samples WHERE id = ?", [id]))
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(row_to_json(rows[1])))
        end
    end
end
