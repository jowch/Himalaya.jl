const _SETUP_PATCH_FIELDS = ("name", "included")

function _setup_patch(db, row, body)
    unknown = sort!([String(k) for k in keys(body) if !(String(k) in _SETUP_PATCH_FIELDS)])
    isempty(unknown) || return "not a setup field: $(join(unknown, ", ")); only name and included can be changed"
    sets = String[]; params = Any[]
    changed = Dict{Symbol,Any}(:experiment_id => Int(row.experiment_id))
    if haskey(body, :name) || haskey(body, "name")
        v = _body_get(body, :name)
        v === nothing || v isa AbstractString || return "name must be a string or null"
        display = v === nothing || isempty(strip(v)) || strip(v) == row.folder ? nothing : String(strip(v))
        push!(sets, "display_name = ?"); push!(params, something(display, missing))
        changed[:display_name] = display
        changed[:name] = something(display, String(row.folder))
    end
    if haskey(body, :included) || haskey(body, "included")
        v = _body_get(body, :included)
        v isa Bool || return "included must be true or false"
        if !v && Int(row.included) == 1
            others = sqlrows(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM setups WHERE experiment_id = ? AND included = 1 AND id != ?",
                [Int(row.experiment_id), Int(row.id)]))
            Int(first(others).c) == 0 && return "keep at least one setup shown"
        end
        push!(sets, "included = ?"); push!(params, Int(v))
        changed[:included] = v
    end
    isempty(sets) && return "no patchable fields; supply name or included"
    (sets, params, changed)
end

function register_setups_routes!()
    @patch "/api/setups/{id}" function(req::HTTP.Request, id::Int)
        db = current_db()
        row = setup_row(db, id)
        row === nothing && return _json_error(404, Dict(:error => "setup not found"))
        patch = _setup_patch(db, row, json(req))
        patch isa String && return _json_error(400, Dict(:error => patch))
        sets, params, changed = patch
        scan = get(changed, :included, false) === true && Int(row.included) == 0
        resp = with_idempotency(db, req) do
            DBInterface.execute(db, "UPDATE setups SET $(join(sets, ", ")) WHERE id = ?", [params..., id])
            result = apply_event!(InTransaction(), db, req;
                kind = "setup_updated", entity_type = "setup", entity_id = id, payload = changed)
            _enqueue_broadcast_from_result!(result, "setup_updated", "setup", id)
            HTTP.Response(200, ["Content-Type" => "application/json"],
                to_json(_setup_row_to_json(setup_row(db, id))))
        end
        if scan && resp.status == 200
            eid = Int(row.experiment_id)
            _set_ingest_status!(db, eid, "scanning")
            Threads.@spawn run_ingest!(db, eid; phase = "rescan")
        end
        resp
    end
end
