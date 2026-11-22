using HTTP, JSON3, DBInterface, Tables, Oxygen, SQLite
using Himalaya

function _assignment_body(db::SQLite.DB, exposure_id::Int)
    srow = Tables.rowtable(DBInterface.execute(db,
        "SELECT sample_id FROM exposures WHERE id = ?", [exposure_id]))
    sid = isempty(srow) ? missing : srow[1].sample_id
    state = ismissing(sid) ? "indexed" : sample_call_state(db, Int(sid))
    members = Tables.rowtable(DBInterface.execute(db,
        "SELECT index_id FROM assignment_members WHERE exposure_id = ? ORDER BY index_id",
        [exposure_id]))
    Dict(:exposure_id => exposure_id,
         :state       => state,
         :members     => [Int(m.index_id) for m in members])
end

_assignment_post_state(b::AbstractDict) =
    Dict(:state => b[:state], :members => b[:members])

function predicted_q_for_phase(phase_name::AbstractString, basis::Float64)::Vector{Float64}
    bare = last(split(phase_name, '.'))
    P = try
        getfield(Himalaya, Symbol(bare))
    catch
        return Float64[]
    end
    (P isa Type && P <: Himalaya.Phase) || return Float64[]
    ratios = Himalaya.phaseratios(P; normalize=true)
    [basis * r for r in ratios]
end

function _bonnet_for_index(db::SQLite.DB, exposure_id::Integer,
                           phase_name::AbstractString, lattice_d, index_id::Integer)
    (lattice_d === nothing || ismissing(lattice_d)) && return nothing
    P = resolve_phase(phase_name)
    P === nothing && return nothing
    a = Float64(lattice_d)
    a > 0 || return nothing

    anchors = Tables.rowtable(DBInterface.execute(db,
        """SELECT i.phase, i.lattice_d
           FROM assignment_members m JOIN indices i ON i.id = m.index_id
           WHERE m.exposure_id = ? AND i.id != ?
             AND i.lattice_d IS NOT NULL
             AND i.phase IN ('Pn3m', 'Im3m', 'Ia3d')""",
        [Int(exposure_id), Int(index_id)]))

    for anc in anchors
        Pa = resolve_phase(String(anc.phase))
        Pa === nothing && continue
        Pa === P && continue
        pred = Himalaya.bonnet_lattice(Pa, Float64(anc.lattice_d), P)
        pred === nothing && continue
        return Dict(:predicted_a => pred,
                    :consistent  => abs(a - pred) <= 0.02 * pred)
    end
    nothing
end

function register_analysis_routes!()
    @get "/api/exposures/{id}/indices" function(req::HTTP.Request, id::Int)
        db = current_db()
        indices = Tables.rowtable(DBInterface.execute(db,
            """SELECT id, exposure_id, phase, basis, score, r_squared, lattice_d,
                  status, kind, inputs_hash
           FROM indices WHERE exposure_id = ? ORDER BY score DESC""", [id]))
        out = map(indices) do ix
            peak_rows = Tables.rowtable(DBInterface.execute(db,
                """SELECT ip.peak_id, ip.ratio_position, ip.residual,
                          COALESCE(ap.q, pc.q) AS q_observed
                   FROM index_peaks ip
                   LEFT JOIN auto_peaks ap     ON ap.id = ip.peak_id AND ip.peak_kind = 'auto'
                   LEFT JOIN peak_curations pc ON pc.id = ip.peak_id AND ip.peak_kind = 'curation'
                   WHERE ip.index_id = ? ORDER BY ip.ratio_position""",
                [Int(ix.id)]))
            predicted = predicted_q_for_phase(String(ix.phase), Float64(ix.basis))
            d = row_to_json(ix)
            d[:peaks]        = rows_to_json(peak_rows)
            d[:predicted_q]  = predicted
            d[:bonnet]       = _bonnet_for_index(db, id, String(ix.phase), ix.lattice_d, Int(ix.id))
            d
        end
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(out))
    end

    @get "/api/exposures/{id}/assignment" function(req::HTTP.Request, id::Int)
        db = current_db()
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(_assignment_body(db, id)))
    end

    @get "/api/samples/{id}/call" function(req::HTTP.Request, id::Int)
        db = current_db()
        _sample_exists(db, id) || return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "sample not found")))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(_sample_call_body(db, id)))
    end

    @post "/api/samples/{id}/call/state" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if !haskey(body, :state)
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "missing field: state")))
        end
        state = String(body.state)
        if !(state in ("indexed", "form_factor", "null"))
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "invalid state: $state")))
        end
        _sample_exists(db, id) || return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "sample not found")))
        return with_idempotency(db, req) do
            result = apply_event!(InTransaction(), db, req;
                kind        = "sample_call_set_state",
                entity_type = "sample",
                entity_id   = id,
                payload     = Dict(:state => state))
            mates = align_call_state!(db, req, bridge_mates(db, id), state; broadcast = false)
            b = _sample_call_body(db, id)
            _enqueue_broadcast_from_result!(result, "sample_call_set_state", "sample", id;
                post_state = Dict(:call => _sample_call_post_state(b)))
            broadcast_call_states!(db, mates)
            b[:event_id]    = result.event_id
            b[:view_row_id] = result.view_row_id
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(b))
        end
    end

    @post "/api/exposures/{id}/assignment/members" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if !haskey(body, :index_id)
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "missing field: index_id")))
        end
        local index_id::Int
        try
            index_id = Int(body.index_id)
        catch
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "index_id must be an integer")))
        end
        return with_idempotency(db, req) do
            isempty(Tables.rowtable(DBInterface.execute(db,
                "SELECT 1 FROM indices WHERE id = ? AND exposure_id = ?", [index_id, id]))) &&
                return HTTP.Response(409, ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "this index changed after re-analysis, so pick it again")))
            result = apply_event!(InTransaction(), db, req;
                kind        = "assignment_add",
                entity_type = "exposure",
                entity_id   = id,
                payload     = Dict(:index_id => index_id))
            mates = _align_mates_to_indexed!(db, req, id)
            b = _assignment_body(db, id)
            _enqueue_broadcast_from_result!(result, "assignment_add", "exposure", id;
                post_state = Dict(:assignment => _assignment_post_state(b)))
            broadcast_call_states!(db, mates)
            b[:event_id]    = result.event_id
            b[:view_row_id] = result.view_row_id
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(b))
        end
    end

    @delete "/api/exposures/{id}/assignment/members/{index_id}" function(req::HTTP.Request, id::Int, index_id::Int)
        db = current_db()
        return with_idempotency(db, req) do
            result = apply_event!(InTransaction(), db, req;
                kind        = "assignment_remove",
                entity_type = "exposure",
                entity_id   = id,
                payload     = Dict(:index_id => index_id))
            b = _assignment_body(db, id)
            _enqueue_broadcast_from_result!(result, "assignment_remove", "exposure", id;
                post_state = Dict(:assignment => _assignment_post_state(b)))
            b[:event_id]    = result.event_id
            b[:view_row_id] = result.view_row_id
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(b))
        end
    end

    @post "/api/exposures/{id}/speculative" function(req::HTTP.Request, id::Int)
        db = current_db()
        body = json(req)
        for field in (:phase, :anchor_peak_id, :anchor_ratio)
            if !haskey(body, field)
                return HTTP.Response(400,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "missing field: $(field)")))
            end
        end
        local phase_name::String
        local anchor_peak_id::Int
        local anchor_ratio::Int
        try
            phase_name     = String(body.phase)
            anchor_peak_id = Int(body.anchor_peak_id)
            anchor_ratio   = Int(body.anchor_ratio)
        catch
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "phase must be string; anchor_peak_id and anchor_ratio must be integers")))
        end
        return with_idempotency(db, req) do
            additional      = haskey(body, :additional) ? body.additional : []

            P = resolve_phase(phase_name)
            P === nothing && return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "unknown phase: $phase_name")))

            ratio_to_peak = Dict{Int, Int}(anchor_ratio => anchor_peak_id)
            for entry in additional
                rp  = Int(entry.ratio_position)
                pid = Int(entry.peak_id)
                haskey(ratio_to_peak, rp) && return HTTP.Response(400,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "duplicate ratio_position $rp")))
                ratio_to_peak[rp] = pid
            end

            local nid::Int
            try
                nid = insert_speculative_index!(db, id, P, ratio_to_peak)
            catch e
                return HTTP.Response(400,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => sprint(showerror, e))))
            end
            payload = Dict(:index_id => nid)
            result = apply_event!(InTransaction(), db, req;
                kind        = "speculative_created",
                entity_type = "exposure",
                entity_id   = id,
                payload     = payload)

            _enqueue_post_commit_broadcast!(
                Int(result.event_id), "speculative_created", "exposure", Int(id),
                result.user_id, result.client_id, result.client_op_id,
                result.payload_json)

            rows = Tables.rowtable(DBInterface.execute(db,
                """SELECT id, exposure_id, phase, basis, score, r_squared, lattice_d,
                  status, kind, inputs_hash
           FROM indices WHERE id = ?""", [nid]))
            ix = rows[1]
            peak_rows = Tables.rowtable(DBInterface.execute(db,
                """SELECT ip.peak_id, ip.ratio_position, ip.residual,
                          COALESCE(ap.q, pc.q) AS q_observed
                   FROM index_peaks ip
                   LEFT JOIN auto_peaks ap     ON ap.id = ip.peak_id AND ip.peak_kind = 'auto'
                   LEFT JOIN peak_curations pc ON pc.id = ip.peak_id AND ip.peak_kind = 'curation'
                   WHERE ip.index_id = ? ORDER BY ip.ratio_position""", [nid]))
            predicted = predicted_q_for_phase(String(ix.phase), Float64(ix.basis))
            d = row_to_json(ix)
            d[:peaks]       = rows_to_json(peak_rows)
            d[:predicted_q] = predicted
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(d))
        end
    end

    @post "/api/exposures/{id}/custom-index" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        for field in (:phase, :basis)
            if !haskey(body, field)
                return HTTP.Response(400, ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "missing field: $(field)")))
            end
        end
        local phase_name::String
        local basis::Float64
        local ratios::Union{Vector{Float64}, Nothing}
        try
            phase_name = String(body.phase)
            basis      = Float64(body.basis)
            raw = haskey(body, :ratios) ? body.ratios : nothing
            raw === nothing || raw isa AbstractVector ||
                error("ratios must be an array")
            ratios     = raw === nothing ? nothing :
                         Float64[Float64(r) for r in raw]
        catch
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "phase must be string; basis a number; ratios an array of numbers")))
        end
        if ratios !== nothing && (isempty(ratios) || any(r -> !(r > 0), ratios))
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "ratios must be a non-empty array of positive numbers")))
        end
        P = resolve_phase(phase_name)
        P === nothing && return HTTP.Response(400, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "unknown phase: $phase_name")))
        basis > 0 || return HTTP.Response(400, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "basis must be positive")))

        return with_idempotency(db, req) do
            nid = insert_custom_index!(db, id, P, basis; drawn_ratios = ratios)

            sc = apply_event!(InTransaction(), db, req;
                kind        = "speculative_created",
                entity_type = "exposure",
                entity_id   = id,
                payload     = Dict(:index_id => nid))
            _enqueue_post_commit_broadcast!(
                Int(sc.event_id), "speculative_created", "exposure", Int(id),
                sc.user_id, sc.client_id, sc.client_op_id, sc.payload_json)

            aa = apply_event!(InTransaction(), db, req;
                kind        = "assignment_add",
                entity_type = "exposure",
                entity_id   = id,
                payload     = Dict(:index_id => nid))
            mates = _align_mates_to_indexed!(db, req, id)
            ab = _assignment_body(db, id)
            _enqueue_broadcast_from_result!(aa, "assignment_add", "exposure", id;
                post_state = Dict(:assignment =>
                    Dict(:state => ab[:state], :members => ab[:members])))
            broadcast_call_states!(db, mates)

            rows = Tables.rowtable(DBInterface.execute(db,
                """SELECT id, exposure_id, phase, basis, score, r_squared, lattice_d,
                      status, kind, inputs_hash
               FROM indices WHERE id = ?""", [nid]))
            ix = rows[1]
            peak_rows = Tables.rowtable(DBInterface.execute(db,
                """SELECT ip.peak_id, ip.ratio_position, ip.residual,
                          COALESCE(ap.q, pc.q) AS q_observed
                   FROM index_peaks ip
                   LEFT JOIN auto_peaks ap     ON ap.id = ip.peak_id AND ip.peak_kind = 'auto'
                   LEFT JOIN peak_curations pc ON pc.id = ip.peak_id AND ip.peak_kind = 'curation'
                   WHERE ip.index_id = ? ORDER BY ip.ratio_position""", [nid]))
            predicted = predicted_q_for_phase(String(ix.phase), Float64(ix.basis))
            d = row_to_json(ix)
            d[:peaks]       = rows_to_json(peak_rows)
            d[:predicted_q] = predicted
            d[:bonnet]      = _bonnet_for_index(db, id, String(ix.phase), ix.lattice_d, nid)
            d[:event_id]    = sc.event_id
            d[:view_row_id] = sc.view_row_id
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(d))
        end
    end

    @delete "/api/indices/{id}" function(req::HTTP.Request, id::Int)
        db = current_db()
        return with_idempotency(db, req) do
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, kind, exposure_id FROM indices WHERE id = ?", [id]))
            isempty(rows) && return HTTP.Response(404,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "index not found")))
            kind_val = String(rows[1].kind)
            kind_val == "speculative" || return HTTP.Response(403,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "only speculative indices can be deleted; use group exclusion for auto indices")))

            exposure_id = Int(rows[1].exposure_id)
            was_member = !isempty(Tables.rowtable(DBInterface.execute(db,
                "SELECT 1 FROM assignment_members WHERE exposure_id = ? AND index_id = ?",
                [exposure_id, id])))

            DBInterface.execute(db,
                "DELETE FROM index_group_members WHERE index_id = ?", [id])
            DBInterface.execute(db,
                "DELETE FROM index_peaks WHERE index_id = ?", [id])
            DBInterface.execute(db,
                "DELETE FROM indices WHERE id = ?", [id])

            payload = Dict(:index_id => id)
            result = apply_event!(InTransaction(), db, req;
                kind        = "speculative_deleted",
                entity_type = "exposure",
                entity_id   = exposure_id,
                payload     = payload)

            _enqueue_post_commit_broadcast!(
                Int(result.event_id), "speculative_deleted", "exposure",
                Int(exposure_id),
                result.user_id, result.client_id, result.client_op_id,
                result.payload_json)

            if was_member
                ar = apply_event!(InTransaction(), db, req;
                    kind        = "assignment_remove",
                    entity_type = "exposure",
                    entity_id   = exposure_id,
                    payload     = Dict(:index_id => id))
                ab = _assignment_body(db, exposure_id)
                _enqueue_broadcast_from_result!(ar, "assignment_remove", "exposure", exposure_id;
                    post_state = Dict(:assignment => _assignment_post_state(ab)))
            end

            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:deleted => id)))
        end
    end
end
