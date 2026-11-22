"""
    phase_call(db, sample_id) -> NamedTuple

The sample's phase call, as HimalayaUI's Call list shows it:
`(; sample_id, state, phase, call_phases, members)`.

`state ∈ {"indexed", "form_factor", "null"}` (`"indexed"` when the sample has
no call row). `members` are the confirmed indices across all of the sample's
exposures, as `(; index_id, exposure_id, phase, score)` ordered by index id.
`phase` is the best-scoring member's phase and `call_phases` the sorted
distinct member phases; both are empty (`nothing`, `String[]`) unless the state
is `"indexed"` and there is at least one member.

Mirrors `HimalayaUI._sample_call_body`; the contract test guards drift.
"""
function phase_call(db::SQLite.DB, sample_id::Integer)
    sid = Int(sample_id)
    srows = Tables.rowtable(DBInterface.execute(db,
        "SELECT state FROM sample_calls WHERE sample_id = ?", [sid]))
    state = isempty(srows) ? "indexed" : String(srows[1].state)
    members = [(index_id = Int(r.index_id), exposure_id = Int(r.exposure_id),
                phase = String(r.phase), score = r.score === missing ? -Inf : Float64(r.score))
               for r in Tables.rowtable(DBInterface.execute(db, """
                   SELECT m.index_id, m.exposure_id, i.phase, i.score
                     FROM assignment_members m
                     JOIN exposures e ON e.id = m.exposure_id
                     JOIN indices i   ON i.id = m.index_id
                    WHERE e.sample_id = ?
                    ORDER BY m.index_id""", [sid]))]
    called = state == "indexed" && !isempty(members)
    phase = called ? argmax(m -> m.score, members).phase : nothing
    call_phases = called ? sort(unique([m.phase for m in members])) : String[]
    return (; sample_id = sid, state, phase, call_phases, members)
end

"""
    sample_labels(db, sample_id) -> Vector{<:NamedTuple}

The sample's labels that are not removed, as `(; id, label, state, source)` in
the order they were added. `state ∈ {"confirmed", "suggested"}`: a suggested
label (for example an automatic Buffer guess) is not a user's decision, so
filter on `state == "confirmed"` when that matters.
"""
sample_labels(db::SQLite.DB, sample_id::Integer) = Tables.rowtable(DBInterface.execute(db,
    "SELECT id, label, state, source FROM sample_labels WHERE sample_id = ? AND state != 'removed' ORDER BY id",
    [Int(sample_id)]))

"""
    subtraction_background(db, sample_id) -> Union{Int, Nothing}

The id of the sample chosen as this sample's subtraction background (usually a
buffer), or `nothing` when none is set.
"""
function subtraction_background(db::SQLite.DB, sample_id::Integer)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT ref_sample_id FROM sample_subtractions WHERE sample_id = ?", [Int(sample_id)]))
    isempty(rows) ? nothing : Int(rows[1].ref_sample_id)
end

"""
    bridged_samples(db, sample_id) -> Vector{Int}

The other samples a user has bridged with this one (marked as the same
specimen, for example after a mistaken split), ordered by id. Empty when the
sample has no bridge.
"""
bridged_samples(db::SQLite.DB, sample_id::Integer) = [Int(r.sample_id) for r in
    Tables.rowtable(DBInterface.execute(db, """
        SELECT b.sample_id FROM sample_bridges b
          JOIN sample_bridges me ON me.group_id = b.group_id
         WHERE me.sample_id = ? AND b.sample_id != me.sample_id
         ORDER BY b.sample_id""", [Int(sample_id)]))]
