using DBInterface, Tables, SQLite

function _group_call(db, sample_ids::Vector{Int})
    states = sample_call_states(db, sample_ids)
    own = _own_call_members(db, sample_ids)
    members = reduce(vcat, (get(own, s, []) for s in sample_ids); init = [])
    st = [get(states, s, "indexed") for s in sample_ids]
    state = !isempty(members) ? "indexed" :
        "form_factor" in st ? "form_factor" :
        "null" in st ? "null" : "indexed"
    (state = state, members = members, own = own,
     phases = sort(unique([m.phase for m in members])),
     called = export_call_state(state, members) != "uncalled")
end

_side(db, sid::Int) = vcat(sid, bridge_mates(db, sid))

_same_call(x, y) = x.state == y.state && x.phases == y.phases

function _call_conflict_body(db, a::Int, b::Int, ca, cb)
    name(sid) = String(first(Tables.rowtable(DBInterface.execute(db,
        "SELECT name FROM samples WHERE id = ?", [sid]))).name)
    Dict(:error => "these samples have different calls; say which to keep",
         :calls => [Dict(:sample_id => sid, :name => name(sid), :state => c.state, :call_phases => c.phases)
                    for (sid, c) in ((a, ca), (b, cb))])
end

function bridge_call_plan(db, a::Int, b::Int, keep)
    sa, sb = _side(db, a), _side(db, b)
    a in sb && return (ok = true, sides = nothing)
    ca, cb = _group_call(db, sa), _group_call(db, sb)
    conflict = ca.called && cb.called && !_same_call(ca, cb)
    if keep === nothing
        conflict && return (ok = false, status = 409, body = _call_conflict_body(db, a, b, ca, cb))
        winner = ca.called ? ca : cb.called ? cb : nothing
        return (ok = true, sides = (sa, sb), target = winner === nothing ? nothing : winner.state, losers = Int[])
    end
    if keep == "both"
        conflict && (ca.state != "indexed" || cb.state != "indexed") &&
            return (ok = false, status = 422, body = Dict(:error => "both calls can be kept only when both are phase calls"))
        winner = ca.called ? ca : cb.called ? cb : nothing
        return (ok = true, sides = (sa, sb), target = winner === nothing ? nothing : winner.state, losers = Int[])
    end
    keep == a || keep == b ||
        return (ok = false, status = 400, body = Dict(:error => "keep_call must be \"both\" or one of the two sample ids"))
    (wc, losers) = keep == a ? (ca, sb) : (cb, sa)
    (ok = true, sides = (sa, sb), target = wc.called ? wc.state : nothing,
     losers = conflict && wc.state == "indexed" ? losers : Int[])
end

function _emit_call_state!(db, req, sid::Int, state::AbstractString; broadcast::Bool = true)
    result = apply_event!(InTransaction(), db, req;
        kind = "sample_call_set_state", entity_type = "sample", entity_id = sid,
        payload = Dict(:state => String(state)))
    broadcast && _enqueue_broadcast_from_result!(result, "sample_call_set_state", "sample", sid;
        post_state = Dict(:call => _sample_call_post_state(_sample_call_body(db, sid))))
    result
end

function _emit_member_removed!(db, req, exposure_id::Int, index_id::Int)
    result = apply_event!(InTransaction(), db, req;
        kind = "assignment_remove", entity_type = "exposure", entity_id = exposure_id,
        payload = Dict(:index_id => index_id))
    _enqueue_broadcast_from_result!(result, "assignment_remove", "exposure", exposure_id;
        post_state = Dict(:assignment => _assignment_post_state(_assignment_body(db, exposure_id))))
end

function apply_bridge_call_plan!(db, req, plan)
    plan.sides === nothing && return nothing
    for (_, ms) in _own_call_members(db, plan.losers), m in ms
        _emit_member_removed!(db, req, m.exposure_id, m.index_id)
    end
    plan.target === nothing && return nothing
    align_call_state!(db, req, vcat(plan.sides...), plan.target)
end

function align_call_state!(db, req, sample_ids::Vector{Int}, state::AbstractString; broadcast::Bool = true)
    states = sample_call_states(db, sample_ids)
    for sid in sample_ids
        get(states, sid, "indexed") == state || _emit_call_state!(db, req, sid, state; broadcast)
    end
    nothing
end

function bridged_call_drift(db)
    out = Int[]
    for (g, ms) in sort(collect(_live_bridge_groups(db)); by = first)
        states = sample_call_states(db, ms)
        length(unique(get(states, s, "indexed") for s in ms)) > 1 && push!(out, g)
    end
    out
end

function migrate_share_bridged_calls!(db::SQLite.DB)
    _migrated(db, MIGRATION_BRIDGED_CALLS) && return nothing
    SQLite.transaction(db) do
        req = _system_request()
        for (_, ms) in sort(collect(_live_bridge_groups(db)); by = first)
            c = _group_call(db, ms)
            c.called && align_call_state!(db, req, ms, c.state; broadcast = false)
        end
        _record_migration!(db, MIGRATION_BRIDGED_CALLS)
    end
    nothing
end

function _align_mates_to_indexed!(db, req, exposure_id::Int)
    row = Tables.rowtable(DBInterface.execute(db,
        "SELECT sample_id FROM exposures WHERE id = ?", [exposure_id]))
    (isempty(row) || row[1].sample_id === missing) && return nothing
    align_call_state!(db, req, bridge_mates(db, Int(row[1].sample_id)), "indexed")
end
