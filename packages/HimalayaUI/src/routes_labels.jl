using HTTP, JSON3, DBInterface, Tables, Oxygen

const MAX_LABEL_LENGTH = 64

function labels_for_experiment(db, exp_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT l.id, l.sample_id, l.label, l.state, l.source
             FROM sample_labels l JOIN samples s ON s.id = l.sample_id
            WHERE s.experiment_id = ? AND s.merged_into_id IS NULL AND l.state != 'removed'
            ORDER BY l.sample_id, l.id""", [exp_id]))
    [Dict(:id => Int(r.id), :sample_id => Int(r.sample_id), :label => String(r.label),
          :state => String(r.state), :source => String(r.source)) for r in rows]
end

function normalize_label(raw)
    raw isa AbstractString || return nothing
    label = strip(raw)
    (isempty(label) || length(label) > MAX_LABEL_LENGTH || any(c -> c == '\n' || c == '\r', label)) && return nothing
    String(label)
end

function _label_state(db, sample_id::Int, label::String)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT state FROM sample_labels WHERE sample_id = ? AND label = ? COLLATE NOCASE",
        [sample_id, label]))
    isempty(rows) ? nothing : String(rows[1].state)
end

function _emit_label_event!(db, req, kind::String, sample_id::Int, label::String, exp_id::Int)
    result = apply_event!(InTransaction(), db, req;
        kind = kind, entity_type = "sample", entity_id = sample_id,
        payload = Dict(:label => label, :experiment_id => exp_id))
    _enqueue_broadcast_from_result!(result, kind, "sample", sample_id)
end

function add_label!(db, req, sample_id::Int, label::String, exp_id::Int)
    _label_state(db, sample_id, label) == "confirmed" && return
    _emit_label_event!(db, req, "label_added", sample_id, label, exp_id)
end

function remove_label!(db, req, sample_id::Int, label::String, exp_id::Int)
    state = _label_state(db, sample_id, label)
    (state === nothing || state == "removed") && return
    _emit_label_event!(db, req, "label_removed", sample_id, label, exp_id)
end

function _live_sample_experiment(db, sample_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT experiment_id FROM samples WHERE id = ? AND merged_into_id IS NULL", [sample_id]))
    isempty(rows) || ismissing(rows[1].experiment_id) ? nothing : Int(rows[1].experiment_id)
end

function register_label_routes!()
    @get "/api/experiments/{id}/labels" function(req::HTTP.Request, id::Int)
        _json_ok(labels_for_experiment(current_db(), id))
    end

    @post "/api/samples/{id}/labels" function(req::HTTP.Request, id::Int)
        db = current_db()
        body = try
            json(req)
        catch
            return _json_error(400, "body must be JSON")
        end
        label = normalize_label(body isa AbstractDict ? get(body, "label", nothing) : nothing)
        label === nothing && return _json_error(400, "label must be 1 to $MAX_LABEL_LENGTH characters on one line")
        exp_id = _live_sample_experiment(db, id)
        exp_id === nothing && return _json_error(404, "sample not found")

        return with_idempotency(db, req) do
            add_label!(db, req, id, label, exp_id)
            _json_ok(labels_for_experiment(db, exp_id))
        end
    end

    @delete "/api/samples/{id}/labels" function(req::HTTP.Request, id::Int)
        db = current_db()
        params = HTTP.queryparams(HTTP.URI(req.target))
        label = normalize_label(get(params, "label", nothing))
        label === nothing && return _json_error(400, "label is required")
        exp_id = _live_sample_experiment(db, id)
        exp_id === nothing && return _json_error(404, "sample not found")

        return with_idempotency(db, req) do
            remove_label!(db, req, id, label, exp_id)
            _json_ok(labels_for_experiment(db, exp_id))
        end
    end
end
