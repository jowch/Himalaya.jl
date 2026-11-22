using DBInterface, Tables, SQLite

function sample_call_state(db::SQLite.DB, sample_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT state FROM sample_calls WHERE sample_id = ?", [sample_id]))
    isempty(rows) ? "indexed" : String(rows[1].state)
end

_sample_exists(db::SQLite.DB, sample_id::Int) = !isempty(Tables.rowtable(DBInterface.execute(db,
    "SELECT 1 FROM samples WHERE id = ?", [sample_id])))

function _call_groups(db::SQLite.DB, sample_ids::AbstractVector{<:Integer})
    out = Dict(Int(sid) => [Int(sid)] for sid in sample_ids)
    isempty(sample_ids) && return out
    ph = join(fill("?", length(sample_ids)), ",")
    for r in Tables.rowtable(DBInterface.execute(db, """
            SELECT b.sample_id AS sid, m.sample_id AS mate
              FROM sample_bridges b
              JOIN sample_bridges m ON m.group_id = b.group_id AND m.sample_id <> b.sample_id
             WHERE b.sample_id IN ($ph)
             ORDER BY m.sample_id""", collect(Int, sample_ids)))
        push!(out[Int(r.sid)], Int(r.mate))
    end
    out
end

function sample_call_members(db::SQLite.DB, sample_ids::AbstractVector{<:Integer})
    groups = _call_groups(db, sample_ids)
    own = _own_call_members(db, unique(reduce(vcat, values(groups); init = Int[])))
    out = Dict{Int, valtype(own)}()
    for (sid, group) in groups
        ms = reduce(vcat, (get(own, g, eltype(valtype(own))[]) for g in group); init = eltype(valtype(own))[])
        isempty(ms) || (out[sid] = sort(ms; by = m -> m.index_id))
    end
    out
end

function _own_call_members(db::SQLite.DB, sample_ids::AbstractVector{<:Integer})
    out = Dict{Int, Vector{NamedTuple{(:index_id, :exposure_id, :phase, :score), Tuple{Int, Int, String, Float64}}}}()
    isempty(sample_ids) && return out
    ph = join(fill("?", length(sample_ids)), ",")
    for r in Tables.rowtable(DBInterface.execute(db, """
            SELECT e.sample_id AS sid, m.index_id AS iid, m.exposure_id AS eid,
                   i.phase AS phase, i.score AS score
              FROM assignment_members m
              JOIN exposures e ON e.id = m.exposure_id
              JOIN indices i ON i.id = m.index_id
             WHERE e.sample_id IN ($ph)
             ORDER BY m.index_id""", collect(Int, sample_ids)))
        push!(get!(out, Int(r.sid), eltype(valtype(out))[]),
            (index_id = Int(r.iid), exposure_id = Int(r.eid), phase = String(r.phase),
             score = r.score === missing ? -Inf : Float64(r.score)))
    end
    out
end

function sample_call_states(db::SQLite.DB, sample_ids::AbstractVector{<:Integer})
    out = Dict{Int, String}()
    isempty(sample_ids) && return out
    ph = join(fill("?", length(sample_ids)), ",")
    for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id, state FROM sample_calls WHERE sample_id IN ($ph)",
            collect(Int, sample_ids)))
        out[Int(r.sample_id)] = String(r.state)
    end
    out
end

function call_summary(state::AbstractString, members)
    (state == "indexed" && !isempty(members)) || return (phase = nothing, call_phases = String[])
    (phase = argmax(m -> m.score, members).phase,
     call_phases = sort(unique([m.phase for m in members])))
end

export_call_state(state::AbstractString, members) =
    state == "null" || (state == "indexed" && isempty(members)) ? "uncalled" : String(state)

function _sample_call_body(db::SQLite.DB, sample_id::Int)
    state = sample_call_state(db, sample_id)
    ms = get(sample_call_members(db, [sample_id]), sample_id, [])
    sm = call_summary(state, ms)
    Dict(:sample_id   => sample_id,
         :state       => state,
         :phase       => sm.phase,
         :call_phases => sm.call_phases,
         :members     => [Dict(:index_id => m.index_id, :exposure_id => m.exposure_id,
                               :phase => m.phase) for m in ms])
end

_sample_call_post_state(b::AbstractDict) =
    Dict(:state => b[:state], :members => b[:members])

function _upsert_sample_call!(db::SQLite.DB, sample_id::Int, state::AbstractString)
    DBInterface.execute(db,
        """INSERT INTO sample_calls (sample_id, state) VALUES (?, ?)
           ON CONFLICT(sample_id) DO UPDATE SET state = excluded.state""",
        [sample_id, String(state)])
    nothing
end

function _set_sample_call_state!(db::SQLite.DB, sample_id::Int, state::AbstractString;
                                 only_exposure::Union{Int, Nothing} = nothing)
    _upsert_sample_call!(db, sample_id, state)
    state == "indexed" && return nothing
    if only_exposure === nothing
        DBInterface.execute(db,
            """DELETE FROM assignment_members
                WHERE exposure_id IN (SELECT id FROM exposures WHERE sample_id = ?)""",
            [sample_id])
    else
        DBInterface.execute(db,
            "DELETE FROM assignment_members WHERE exposure_id = ?", [only_exposure])
    end
    nothing
end

function _set_exposure_sample_state!(db::SQLite.DB, exposure_id::Int, state::AbstractString)
    DBInterface.execute(db,
        """INSERT INTO sample_calls (sample_id, state)
           SELECT sample_id, ? FROM exposures
            WHERE id = ? AND sample_id IS NOT NULL
           ON CONFLICT(sample_id) DO UPDATE SET state = excluded.state""",
        [String(state), exposure_id])
    nothing
end

function _index_sample_if_called!(db::SQLite.DB, exposure_id::Int)
    has = Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM assignment_members WHERE exposure_id = ? LIMIT 1", [exposure_id]))
    isempty(has) || _set_exposure_sample_state!(db, exposure_id, "indexed")
    nothing
end

_call_state_effect(kind::AbstractString, payload) =
    kind == "assignment_add" ? "indexed" :
    kind == "assignment_set_state" ? String(payload.state) : nothing

_call_rows(db::SQLite.DB, where_sql::AbstractString, params) =
    collect(Tables.rowtable(DBInterface.execute(db,
        """SELECT id, action, entity_type, entity_id, payload FROM user_actions
            WHERE payload IS NOT NULL AND ($where_sql)""", params)))

function _rebuild_exposure_with_call!(db::SQLite.DB, exposure_id::Int)
    rows = _call_rows(db, "entity_type = 'exposure' AND entity_id = ?", [exposure_id])
    srow = Tables.rowtable(DBInterface.execute(db,
        "SELECT sample_id FROM exposures WHERE id = ?", [exposure_id]))
    sid = isempty(srow) || srow[1].sample_id === missing ? nothing : Int(srow[1].sample_id)
    if sid !== nothing
        append!(rows, _call_rows(db,
            "entity_type = 'sample' AND entity_id = ? AND action = 'sample_call_set_state'", [sid]))
        append!(rows, _call_rows(db,
            """entity_type = 'exposure' AND entity_id <> ?
               AND entity_id IN (SELECT id FROM exposures WHERE sample_id = ?)
               AND action IN ('assignment_add', 'assignment_set_state')""", [exposure_id, sid]))
    end
    sort!(rows; by = r -> r.id)
    for r in rows
        kind = String(r.action)
        payload = JSON3.read(String(r.payload))
        if r.entity_type == "sample"
            _set_sample_call_state!(db, sid, String(payload.state); only_exposure = exposure_id)
        elseif Int(r.entity_id) != exposure_id
            _upsert_sample_call!(db, sid, _call_state_effect(kind, payload))
        else
            update_view_for_event!(db, kind, exposure_id, payload, Int(r.id))
        end
    end
    nothing
end

function _rebuild_sample_with_call!(db::SQLite.DB, sample_id::Int)
    rows = _call_rows(db, "entity_type = 'sample' AND entity_id = ?", [sample_id])
    append!(rows, _call_rows(db,
        """entity_type = 'exposure'
           AND entity_id IN (SELECT id FROM exposures WHERE sample_id = ?)
           AND action IN ('assignment_add', 'assignment_remove', 'assignment_set_state')""",
        [sample_id]))
    sort!(rows; by = r -> r.id)
    for r in rows
        update_view_for_event!(db, String(r.action), Int(r.entity_id),
                               JSON3.read(String(r.payload)), Int(r.id))
    end
    nothing
end
