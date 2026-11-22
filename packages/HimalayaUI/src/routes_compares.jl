using HTTP, JSON3, DBInterface, Tables, Oxygen

function compares_for_experiment(db, exp_id::Int)
    compares = Tables.rowtable(DBInterface.execute(db,
        """SELECT c.sample_id, c.ref_sample_id, u.username
             FROM sample_compares c
             JOIN samples s ON s.id = c.sample_id
             JOIN samples r ON r.id = c.ref_sample_id
             LEFT JOIN users u ON u.id = c.added_by
            WHERE s.experiment_id = ? AND s.merged_into_id IS NULL AND r.merged_into_id IS NULL
            ORDER BY c.sample_id, c.id""", [exp_id]))
    subs = Tables.rowtable(DBInterface.execute(db,
        """SELECT t.sample_id, t.ref_sample_id
             FROM sample_subtractions t
             JOIN samples s ON s.id = t.sample_id
             JOIN samples r ON r.id = t.ref_sample_id
            WHERE s.experiment_id = ? AND s.merged_into_id IS NULL AND r.merged_into_id IS NULL
            ORDER BY t.sample_id""", [exp_id]))
    Dict(
        :compares => [Dict(:sample_id => Int(r.sample_id), :ref_sample_id => Int(r.ref_sample_id),
                           :added_by => ismissing(r.username) ? nothing : String(r.username)) for r in compares],
        :subtractions => [Dict(:sample_id => Int(r.sample_id), :ref_sample_id => Int(r.ref_sample_id)) for r in subs],
    )
end

_int_list(raw) = raw === nothing ? Int[] :
    (raw isa AbstractVector && all(x -> x isa Integer, raw)) ? Int.(collect(raw)) : nothing

_current_compares(db, sample_id::Int) = Set(Int(r.ref_sample_id) for r in Tables.rowtable(DBInterface.execute(db,
    "SELECT ref_sample_id FROM sample_compares WHERE sample_id = ?", [sample_id])))

function _current_subtraction(db, sample_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT ref_sample_id FROM sample_subtractions WHERE sample_id = ?", [sample_id]))
    isempty(rows) ? nothing : Int(rows[1].ref_sample_id)
end

const MAX_COMPARES = 5

function compare_changes(current::Set{Int}, subtraction::Union{Int, Nothing},
                         add::Vector{Int}, remove::Vector{Int}, sub_given::Bool, sub)
    adding = unique(filter(r -> r ∉ current && r ∉ remove, add))
    removing = unique(filter(r -> r ∈ current, remove))
    payload = Dict{Symbol, Any}()
    isempty(adding) || (payload[:add] = adding)
    isempty(removing) || (payload[:remove] = removing)
    if sub_given
        sub != subtraction && (payload[:subtraction] = sub)
    elseif subtraction !== nothing && subtraction ∈ removing
        payload[:subtraction] = nothing
    end
    payload
end

function register_compare_routes!()
    @get "/api/experiments/{id}/compares" function(req::HTTP.Request, id::Int)
        _json_ok(compares_for_experiment(current_db(), id))
    end

    @patch "/api/samples/{id}/compares" function(req::HTTP.Request, id::Int)
        db = current_db()
        body = try
            json(req)
        catch
            return _json_error(400, "body must be JSON")
        end
        body isa AbstractDict || return _json_error(400, "body must be a JSON object")
        add = _int_list(get(body, "add", nothing))
        remove = _int_list(get(body, "remove", nothing))
        (add === nothing || remove === nothing) && return _json_error(400, "add and remove must be lists of sample ids")
        sub_given = haskey(body, "subtraction")
        sub = sub_given ? body["subtraction"] : nothing
        (sub === nothing || sub isa Integer) || return _json_error(400, "subtraction must be a sample id or null")
        sub = sub === nothing ? nothing : Int(sub)
        refs = unique(vcat(add, sub === nothing ? Int[] : [sub]))
        id ∈ refs && return _json_error(400, "a sample cannot be compared with itself")
        ids = unique(vcat(id, refs))
        marks = join(fill("?", length(ids)), ", ")
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, experiment_id FROM samples WHERE id IN ($marks) AND merged_into_id IS NULL", ids))
        length(rows) == length(ids) || return _json_error(404, "sample not found")
        length(unique(Int(r.experiment_id) for r in rows)) == 1 ||
            return _json_error(422, "samples are in different experiments")
        exp_id = Int(rows[1].experiment_id)

        return with_idempotency(db, req) do
            current = _current_compares(db, id)
            payload = compare_changes(current, _current_subtraction(db, id),
                                      add, remove, sub_given, sub)
            length(current) + length(get(payload, :add, Int[])) - length(get(payload, :remove, Int[])) > MAX_COMPARES &&
                return _json_error(422, "a sample has at most $(MAX_COMPARES) comparisons")
            if !isempty(payload)
                payload[:experiment_id] = exp_id
                result = apply_event!(InTransaction(), db, req;
                    kind = "compares_changed", entity_type = "sample", entity_id = id, payload = payload)
                _enqueue_broadcast_from_result!(result, "compares_changed", "sample", id)
            end
            _json_ok(compares_for_experiment(db, exp_id))
        end
    end
end
