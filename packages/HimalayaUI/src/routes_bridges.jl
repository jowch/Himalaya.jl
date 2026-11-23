using HTTP, JSON3, DBInterface, Tables, Oxygen

_json_ok(body) = HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(body))

function bridge_groups_for_experiment(db, exp_id::Int)
    rows = sqlrows(DBInterface.execute(db,
        """SELECT b.sample_id, b.group_id
             FROM sample_bridges b JOIN samples s ON s.id = b.sample_id
            WHERE s.experiment_id = ?
            ORDER BY b.group_id, b.sample_id""", [exp_id]))
    sizes = Dict{Int, Int}()
    for r in rows
        sizes[Int(r.group_id)] = get(sizes, Int(r.group_id), 0) + 1
    end
    [Dict(:sample_id => Int(r.sample_id), :group_id => Int(r.group_id))
     for r in rows if sizes[Int(r.group_id)] >= 2]
end

function _bridge_group_of(db, sample_id::Int)
    rows = sqlrows(DBInterface.execute(db,
        "SELECT group_id FROM sample_bridges WHERE sample_id = ?", [sample_id]))
    isempty(rows) ? nothing : Int(rows[1].group_id)
end

_bridge_members(db, group_id::Int) = [Int(r.sample_id) for r in sqlrows(DBInterface.execute(db,
    "SELECT sample_id FROM sample_bridges WHERE group_id = ? ORDER BY sample_id", [group_id]))]

_next_bridge_group_id(db) = Int(first(sqlrows(DBInterface.execute(db,
    "SELECT COALESCE(MAX(group_id), 0) + 1 AS g FROM sample_bridges"))).g)

function _live_sample_count(db, sample_ids::Vector{Int})
    isempty(sample_ids) && return 0
    marks = join(fill("?", length(sample_ids)), ", ")
    Int(first(sqlrows(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM samples WHERE id IN ($marks)",
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

function bridge_mates(db, sample_id::Int)
    g = _bridge_group_of(db, sample_id)
    g === nothing && return Int[]
    [Int(r.id) for r in sqlrows(DBInterface.execute(db,
        """SELECT s.id FROM sample_bridges b JOIN samples s ON s.id = b.sample_id
            WHERE b.group_id = ? AND s.id <> ?
            ORDER BY s.id""", [g, sample_id]))]
end

function _confirmed_shared_labels(db, sample_ids::Vector{Int})
    isempty(sample_ids) && return Dict{String, String}()
    marks = join(fill("?", length(sample_ids)), ", ")
    out = Dict{String, String}()
    for r in sqlrows(DBInterface.execute(db,
            "SELECT label FROM sample_labels WHERE sample_id IN ($marks) AND state = 'confirmed' ORDER BY id",
            sample_ids))
        label = String(r.label)
        shared_label(label) && get!(out, lowercase(label), label)
    end
    out
end

function _share_labels!(db, req, sample_id::Int, exp_id::Int)
    group = vcat(sample_id, bridge_mates(db, sample_id))
    labels = collect(values(_confirmed_shared_labels(db, group)))
    for sid in group, label in labels
        _add_label_one!(db, req, sid, label, exp_id)
    end
end

function bridged_label_drift(db)
    sets(sid) = Set(keys(_confirmed_shared_labels(db, [sid])))
    sort([g for (g, ms) in _live_bridge_groups(db) if length(unique(sets.(ms))) > 1])
end

function _live_bridge_groups(db)
    groups = Dict{Int, Vector{Int}}()
    for r in sqlrows(DBInterface.execute(db,
            "SELECT group_id, sample_id FROM sample_bridges ORDER BY group_id, sample_id"))
        push!(get!(groups, Int(r.group_id), Int[]), Int(r.sample_id))
    end
    filter(kv -> length(last(kv)) >= 2, groups)
end

function migrate_share_bridged_labels!(db::SQLite.DB)
    _migrated(db, MIGRATION_BRIDGED_LABELS) && return nothing
    without_review_pins() do
    SQLite.transaction(db) do
        req = _system_request()
        for (_, members) in sort(collect(_live_bridge_groups(db)); by = first)
            exp_id = Int(first(sqlrows(DBInterface.execute(db,
                "SELECT experiment_id FROM samples WHERE id = ?", [first(members)]))).experiment_id)
            labels = sort(collect(values(_confirmed_shared_labels(db, members))))
            for sid in members, label in labels
                _label_state(db, sid, label) == "confirmed" && continue
                apply_event!(InTransaction(), db, req; kind = "label_added", entity_type = "sample",
                    entity_id = sid, payload = Dict(:label => label, :experiment_id => exp_id))
            end
        end
        _record_migration!(db, MIGRATION_BRIDGED_LABELS)
    end end
    nothing
end

function bridge_samples!(db, req, a::Int, b::Int, exp_id::Int; keep = nothing)
    plan = bridge_call_plan(db, a, b, keep)
    plan.ok || return plan
    moves = bridge_moves(_bridge_group_of(db, a), _bridge_group_of(db, b), a, b,
                         g -> _bridge_members(db, g), () -> _next_bridge_group_id(db))
    for (sid, g) in moves
        _emit_bridge_event!(db, req, "sample_bridged", sid,
            Dict(:group_id => g, :experiment_id => exp_id))
    end
    apply_bridge_call_plan!(db, req, plan)
    _share_labels!(db, req, a, exp_id)
    plan
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
        both = sqlrows(DBInterface.execute(db,
            "SELECT id, experiment_id FROM samples WHERE id IN (?, ?)", [id, to]))
        length(both) == 2 || return _json_error(404, "sample not found")
        Int(both[1].experiment_id) == Int(both[2].experiment_id) ||
            return _json_error(422, "samples are in different experiments")
        exp_id = Int(both[1].experiment_id)

        keep = body isa AbstractDict ? get(body, "keep_call", nothing) : nothing
        keep isa Integer && (keep = Int(keep))
        return with_idempotency(db, req) do
            plan = bridge_samples!(db, req, id, to, exp_id; keep)
            plan.ok || return _json_error(plan.status, plan.body[:error]; filter(kv -> first(kv) != :error, plan.body)...)
            _json_ok(bridge_groups_for_experiment(db, exp_id))
        end
    end

    @delete "/api/samples/{id}/bridge" function(req::HTTP.Request, id::Int)
        db = current_db()
        rows = sqlrows(DBInterface.execute(db,
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
