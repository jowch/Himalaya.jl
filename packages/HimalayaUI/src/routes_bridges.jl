using HTTP, JSON3, DBInterface, Tables, Oxygen

_json_ok(body) = HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(body))

function bridge_groups_for_experiment(db, exp_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT b.sample_id, b.group_id
             FROM sample_bridges b JOIN samples s ON s.id = b.sample_id
            WHERE s.experiment_id = ? AND s.merged_into_id IS NULL
            ORDER BY b.group_id, b.sample_id""", [exp_id]))
    sizes = Dict{Int, Int}()
    for r in rows
        sizes[Int(r.group_id)] = get(sizes, Int(r.group_id), 0) + 1
    end
    [Dict(:sample_id => Int(r.sample_id), :group_id => Int(r.group_id))
     for r in rows if sizes[Int(r.group_id)] >= 2]
end

function _bridge_group_of(db, sample_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT group_id FROM sample_bridges WHERE sample_id = ?", [sample_id]))
    isempty(rows) ? nothing : Int(rows[1].group_id)
end

_bridge_members(db, group_id::Int) = [Int(r.sample_id) for r in Tables.rowtable(DBInterface.execute(db,
    "SELECT sample_id FROM sample_bridges WHERE group_id = ? ORDER BY sample_id", [group_id]))]

_next_bridge_group_id(db) = Int(first(Tables.rowtable(DBInterface.execute(db,
    "SELECT COALESCE(MAX(group_id), 0) + 1 AS g FROM sample_bridges"))).g)

function _live_sample_count(db, sample_ids::Vector{Int})
    isempty(sample_ids) && return 0
    marks = join(fill("?", length(sample_ids)), ", ")
    Int(first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM samples WHERE id IN ($marks) AND merged_into_id IS NULL",
        sample_ids))).n)
end

function bridge_moves(ga::Union{Int, Nothing}, gb::Union{Int, Nothing}, a::Int, b::Int,
                      members_of::Function, fresh::Function)
    ga === nothing && gb === nothing && (g = fresh(); return [(a, g), (b, g)])
    ga === nothing && return [(a, gb)]
    gb === nothing && return [(b, ga)]
    ga == gb && return Tuple{Int, Int}[]
    [(s, ga) for s in members_of(gb)]
end

function _emit_bridge_event!(db, req, kind::String, sample_id::Int, payload)
    result = apply_event!(InTransaction(), db, req;
        kind = kind, entity_type = "sample", entity_id = sample_id, payload = payload)
    _enqueue_broadcast_from_result!(result, kind, "sample", sample_id)
end

function bridge_samples!(db, req, a::Int, b::Int, exp_id::Int)
    moves = bridge_moves(_bridge_group_of(db, a), _bridge_group_of(db, b), a, b,
                         g -> _bridge_members(db, g), () -> _next_bridge_group_id(db))
    for (sid, g) in moves
        _emit_bridge_event!(db, req, "sample_bridged", sid,
            Dict(:group_id => g, :experiment_id => exp_id))
    end
end

function register_bridge_routes!()
    @get "/api/experiments/{id}/bridges" function(req::HTTP.Request, id::Int)
        _json_ok(bridge_groups_for_experiment(current_db(), id))
    end

    @post "/api/samples/{id}/bridge" function(req::HTTP.Request, id::Int)
        db = current_db()
        body = try
            json(req)
        catch
            return _json_error(400, "body must be JSON")
        end
        raw = body isa AbstractDict ? get(body, "to_sample_id", nothing) : nothing
        to = raw isa Integer ? Int(raw) : nothing
        to === nothing && return _json_error(400, "to_sample_id is required")
        to == id && return _json_error(400, "a sample cannot be bridged to itself")
        both = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, experiment_id FROM samples WHERE id IN (?, ?) AND merged_into_id IS NULL", [id, to]))
        length(both) == 2 || return _json_error(404, "sample not found")
        Int(both[1].experiment_id) == Int(both[2].experiment_id) ||
            return _json_error(422, "samples are in different experiments")
        exp_id = Int(both[1].experiment_id)

        return with_idempotency(db, req) do
            bridge_samples!(db, req, id, to, exp_id)
            _json_ok(bridge_groups_for_experiment(db, exp_id))
        end
    end

    @delete "/api/samples/{id}/bridge" function(req::HTTP.Request, id::Int)
        db = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT experiment_id FROM samples WHERE id = ?", [id]))
        isempty(rows) && return _json_error(404, "sample not found")
        exp_id = Int(rows[1].experiment_id)

        return with_idempotency(db, req) do
            g = _bridge_group_of(db, id)
            if g !== nothing
                rest = filter(!=(id), _bridge_members(db, g))
                leaving = _live_sample_count(db, rest) <= 1 ? vcat(id, rest) : [id]
                for sid in leaving
                    _emit_bridge_event!(db, req, "sample_unbridged", sid,
                        Dict(:group_id => g, :experiment_id => exp_id))
                end
            end
            _json_ok(bridge_groups_for_experiment(db, exp_id))
        end
    end
end
