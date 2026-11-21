using DBInterface, Tables, SQLite

"""
    sample_call_state(db, sample_id) -> String

The sample's phase-call state (#350): 'indexed' unless the user marked the
sample as form factor or null.
"""
function sample_call_state(db::SQLite.DB, sample_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT state FROM sample_calls WHERE sample_id = ?", [sample_id]))
    isempty(rows) ? "indexed" : String(rows[1].state)
end

_sample_exists(db::SQLite.DB, sample_id::Int) = !isempty(Tables.rowtable(DBInterface.execute(db,
    "SELECT 1 FROM samples WHERE id = ?", [sample_id])))

"""
    sample_call_members(db, sample_ids) -> Dict{Int, Vector}

Every member of each sample's call, from any of its frames, kept or dropped,
as `(index_id, exposure_id, phase, score)`.
"""
function sample_call_members(db::SQLite.DB, sample_ids::AbstractVector{<:Integer})
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

"""
    call_summary(state, members) -> (; phase, call_phases)

`phase` is the top-scored member's phase, the one a status chip shows.
`call_phases` lists every phase in the call, sorted. Both are empty unless
the state is 'indexed'.
"""
function call_summary(state::AbstractString, members)
    (state == "indexed" && !isempty(members)) || return (phase = nothing, call_phases = String[])
    (phase = argmax(m -> m.score, members).phase,
     call_phases = sort(unique([m.phase for m in members])))
end

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
