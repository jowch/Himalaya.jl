using HTTP, JSON3, DBInterface, Tables, Oxygen

function compares_for_experiment(db, exp_id::Int)
    compares = sqlrows(DBInterface.execute(db,
        """SELECT c.sample_id, c.ref_sample_id, u.username
             FROM sample_compares c
             JOIN samples s ON s.id = c.sample_id
             LEFT JOIN users u ON u.id = c.added_by
            WHERE s.experiment_id = ?
            ORDER BY c.sample_id, c.id""", [exp_id]))
    subs = sqlrows(DBInterface.execute(db,
        """SELECT t.sample_id, t.ref_sample_id
             FROM sample_subtractions t
             JOIN samples s ON s.id = t.sample_id
            WHERE s.experiment_id = ?
            ORDER BY t.sample_id""", [exp_id]))
    Dict(
        :compares => [Dict(:sample_id => Int(r.sample_id), :ref_sample_id => Int(r.ref_sample_id),
                           :added_by => ismissing(r.username) ? nothing : String(r.username)) for r in compares],
        :subtractions => [Dict(:sample_id => Int(r.sample_id), :ref_sample_id => Int(r.ref_sample_id))
                          for r in subs if !ismissing(r.ref_sample_id)],
        :backgrounds => backgrounds_for_experiment(db, exp_id, subs),
    )
end

function backgrounds_for_experiment(db, exp_id::Int, subs)
    stored = Dict(Int(r.sample_id) => (ismissing(r.ref_sample_id) ? nothing : Int(r.ref_sample_id)) for r in subs)
    [begin
         default = run.is_buffer || isempty(run.run_buffers) ? nothing : first(run.run_buffers)
         Dict(:sample_id => id,
              :ref_sample_id => haskey(stored, id) ? stored[id] : default,
              :source => haskey(stored, id) ? "user" : "default",
              :default_ref_sample_id => default,
              :is_buffer => run.is_buffer,
              :run_buffer_ids => run.run_buffers)
     end
     for (id, run) in sort!(collect(autosampler_runs(db, exp_id)); by = first)]
end

_int_list(raw) = raw === nothing ? Int[] :
    (raw isa AbstractVector && all(x -> x isa Integer, raw)) ? Int.(collect(raw)) : nothing

_current_compares(db, sample_id::Int) = Set(Int(r.ref_sample_id) for r in sqlrows(DBInterface.execute(db,
    "SELECT ref_sample_id FROM sample_compares WHERE sample_id = ?", [sample_id])))

function _current_subtraction(db, sample_id::Int)
    rows = sqlrows(DBInterface.execute(db,
        "SELECT ref_sample_id FROM sample_subtractions WHERE sample_id = ?", [sample_id]))
    isempty(rows) ? nothing : ismissing(rows[1].ref_sample_id) ? :none : Int(rows[1].ref_sample_id)
end

const MAX_COMPARES = 5

function compare_changes(current::Set{Int}, subtraction::Union{Int, Symbol, Nothing},
                         add::Vector{Int}, remove::Vector{Int}, sub_given::Bool, sub;
                         managed::Bool = false)
    adding = unique(filter(r -> r ∉ current && r ∉ remove, add))
    removing = unique(filter(r -> r ∈ current, remove))
    payload = Dict{Symbol, Any}()
    isempty(adding) || (payload[:add] = adding)
    isempty(removing) || (payload[:remove] = removing)
    if managed
        if sub_given
            wanted = sub == "default" ? nothing : sub === nothing ? :none : sub
            wanted == subtraction || (payload[:subtraction] = wanted === :none ? "none" : wanted)
        end
        return payload
    end
    subtraction = subtraction === :none ? nothing : subtraction
    if sub_given
        wanted = sub == "default" ? nothing : sub
        wanted != subtraction && (payload[:subtraction] = wanted)
    elseif subtraction !== nothing && (subtraction ∈ removing ||
                                       (subtraction ∉ current && !(isempty(adding) && isempty(removing))))
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
        (sub === nothing || sub isa Integer || sub == "default") ||
            return _json_error(400, "subtraction must be a sample id, null or \"default\"")
        sub = sub isa Integer ? Int(sub) : sub
        refs = unique(vcat(add, sub isa Int ? [sub] : Int[]))
        id ∈ refs && return _json_error(400, "a sample cannot be compared with itself")
        ids = unique(vcat(id, refs))
        marks = join(fill("?", length(ids)), ", ")
        rows = sqlrows(DBInterface.execute(db,
            "SELECT id, experiment_id FROM samples WHERE id IN ($marks)", ids))
        length(rows) == length(ids) || return _json_error(404, "sample not found")
        length(unique(Int(r.experiment_id) for r in rows)) == 1 ||
            return _json_error(422, "samples are in different experiments")
        exp_id = Int(rows[1].experiment_id)

        return with_idempotency(db, req) do
            current = _current_compares(db, id)
            payload = compare_changes(current, _current_subtraction(db, id),
                                      add, remove, sub_given, sub;
                                      managed = id in autosampler_samples(db, [id]))
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
