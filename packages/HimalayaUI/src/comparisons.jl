using SHA, JSON3, SQLite, DBInterface, Tables
using Dates

"""
    compute_content_hash(db, comparison_id) -> String

SHA-256 over a canonical serialization of the comparison's identity:
title, description, and members ordered by `display_order` (each member's
exposure_id, render parameters, peak_display, snapshot).

Used by the dispatcher (`comparison_created` / `comparison_submitted`) to
populate `comparisons.content_hash`. Mirrors the client-side
`contentHash.ts` so chat citations like `@comparison:42@<hash8>` resolve
identically across server and browser.

NULL columns serialize as JSON `null`; member ordering is deterministic
(`display_order ASC, id ASC`). Snapshot and peak_display are read as their
already-canonicalized JSON strings — re-parsing them through `JSON3.read`
+ re-emit would re-key them and create drift between server and client.

Returns the lowercase hex digest with a `sha256:` prefix to match the
spec's hash format.
"""
function compute_content_hash(db::SQLite.DB, comparison_id::Integer)::String
    cmp_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT title, description FROM comparisons WHERE id = ?", [Int(comparison_id)]))
    isempty(cmp_rows) && error("compute_content_hash: comparison $comparison_id not found")
    cmp = cmp_rows[1]
    title = String(cmp.title)
    description = ismissing(cmp.description) ? nothing : String(cmp.description)

    member_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT id, exposure_id, display_order, band_height, y_offset,
                  normalization, color_override, label_override,
                  q_window_min, q_window_max, peak_display, snapshot
           FROM comparison_members
           WHERE comparison_id = ?
           ORDER BY display_order ASC, id ASC""", [Int(comparison_id)]))

    members = Vector{Any}(undef, length(member_rows))
    for (i, m) in enumerate(member_rows)
        # Already-canonical JSON columns are re-parsed so the encoded form
        # is structural (a nested object) rather than a quoted string —
        # otherwise an embedded quote in the JSON literal would shift the
        # hash on round-trips.
        snap_obj = ismissing(m.snapshot) ? nothing : JSON3.read(String(m.snapshot))
        peak_obj = ismissing(m.peak_display) ? nothing : JSON3.read(String(m.peak_display))
        members[i] = Dict{Symbol,Any}(
            :exposure_id    => ismissing(m.exposure_id)   ? nothing : Int(m.exposure_id),
            :display_order  => Int(m.display_order),
            :band_height    => Float64(m.band_height),
            :y_offset       => Float64(m.y_offset),
            :normalization  => String(m.normalization),
            :color_override => ismissing(m.color_override) ? nothing : String(m.color_override),
            :label_override => ismissing(m.label_override) ? nothing : String(m.label_override),
            :q_window_min   => ismissing(m.q_window_min)   ? nothing : Float64(m.q_window_min),
            :q_window_max   => ismissing(m.q_window_max)   ? nothing : Float64(m.q_window_max),
            :peak_display   => peak_obj,
            :snapshot       => snap_obj,
        )
    end

    payload = Dict{Symbol,Any}(
        :title       => title,
        :description => description,
        :members     => members,
    )
    canonical = JSON3.write(payload)
    "sha256:" * bytes2hex(SHA.sha256(canonical))
end

"""
    member_ids_for_comparison(db, comparison_id) -> Set{Int}

Set of current `comparison_members.id` values for a comparison. Used by the
`comparison_submitted` dispatcher diff to identify rows to DELETE (in DB,
absent from payload).
"""
function member_ids_for_comparison(db::SQLite.DB, comparison_id::Integer)::Set{Int}
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM comparison_members WHERE comparison_id = ?",
        [Int(comparison_id)]))
    Set{Int}(Int(r.id) for r in rows)
end

"""
    user_id_for_event(db, event_id) -> Union{Int, Nothing}

Resolve `user_actions.user_id` for an event row id. Returns `nothing` for
NULL `user_id` (system-emitted events) or unknown event ids.
"""
function user_id_for_event(db::SQLite.DB, event_id::Integer)::Union{Int, Nothing}
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT user_id FROM user_actions WHERE id = ?", [Int(event_id)]))
    isempty(rows) && return nothing
    ismissing(rows[1].user_id) ? nothing : Int(rows[1].user_id)
end

"""
    comparison_now_iso() -> String

Server-local ISO-8601 UTC timestamp used by the comparison dispatcher to
populate `created_at` / `updated_at`. Matches the format used by
`broadcast_event!` so timeline scans line up across tables.
"""
function comparison_now_iso()::String
    Dates.format(Dates.now(Dates.UTC), Dates.dateformat"yyyy-mm-ddTHH:MM:SS.sssZ")
end
