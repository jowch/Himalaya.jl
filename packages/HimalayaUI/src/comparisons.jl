using SHA, JSON3, SQLite, DBInterface, Tables
using Dates

"""
    CONFIRMED_INDEX_R2_GATE

Shared R² hard-gate for `confirmed_index` snapshots. Mirrors the threshold
used by the frontend `PhasePanel` so an index that clears the UI hide-low-R²
filter is the same one that lands in the comparison snapshot. Bumping this
should be a deliberate edit in lockstep with `frontend/.../PhasePanel.tsx`.
"""
const CONFIRMED_INDEX_R2_GATE = 0.98

"""
    canonical_json(x) -> String

Deterministic JSON serialization for content-hash inputs. **Object keys
are sorted alphabetically** at every level so the byte stream is
unambiguous across implementations — Julia `Dict{Symbol,Any}` iteration
order is hash-based and not portable to JS.

Numbers are emitted in their JSON3 default form (which matches the JS
`JSON.stringify` output for IEEE-754 doubles). `nothing` becomes `null`.
Strings, booleans, integers all delegate to `JSON3.write`.

Mirrors `canonical_json` in `frontend/src/lib/comparison/contentHash.ts`
— the cross-language fixture in `test/contentHash.test.ts` /
`test_comparisons.jl` pins them to identical bytes for a fixed input.
"""
function canonical_json(x::AbstractDict)::String
    keys_sorted = sort(collect(keys(x)); by = k -> string(k))
    parts = String[]
    for k in keys_sorted
        push!(parts, JSON3.write(string(k)) * ":" * canonical_json(x[k]))
    end
    "{" * join(parts, ",") * "}"
end

# JSON3.Object behaves like an ordered Dict but iteration is insertion-
# ordered (driven by the source bytes). Re-emit through the AbstractDict
# branch by converting keys → strings and recursing on each value.
function canonical_json(x::JSON3.Object)::String
    keys_sorted = sort(collect(keys(x)); by = k -> string(k))
    parts = String[]
    for k in keys_sorted
        push!(parts, JSON3.write(string(k)) * ":" * canonical_json(x[k]))
    end
    "{" * join(parts, ",") * "}"
end

canonical_json(x::AbstractVector)::String =
    "[" * join((canonical_json(v) for v in x), ",") * "]"
canonical_json(x::JSON3.Array)::String =
    "[" * join((canonical_json(v) for v in x), ",") * "]"
canonical_json(::Nothing)::String = "null"
canonical_json(::Missing)::String = "null"
canonical_json(x::Bool)::String = x ? "true" : "false"

# Numbers: emit in a JS-compatible form. JS `JSON.stringify(1.0)` → "1"
# (no trailing `.0`), but Julia JSON3.write(1.0) → "1.0". Cross-language
# parity requires we collapse integer-valued floats to their integer
# representation.
function canonical_json(x::AbstractFloat)::String
    isnan(x) || isinf(x) ? "null" : (
        # Integer-valued? Emit without ".0" so it matches JS JSON.stringify.
        (isfinite(x) && trunc(x) == x && abs(x) < 1e16) ?
            string(Int(x)) :
            JSON3.write(x))
end

# Fallback: integers, strings — JSON3 already produces canonical output
# that matches JSON.stringify byte-for-byte.
canonical_json(x)::String = JSON3.write(x)

"""
    compute_content_hash(db, comparison_id) -> String

SHA-256 over a canonical serialization of the comparison's identity:
title, description, and members ordered by `display_order` (each member's
exposure_id, render parameters, peak_display, snapshot, **and
`display_order` itself** — included as a tuple field so the hash is
unambiguous when display_order ties exist; see `canonical_json`).

Used by the dispatcher (`comparison_created` / `comparison_submitted`) to
populate `comparisons.content_hash`. Mirrors the client-side
`contentHash.ts` so chat citations like `@comparison:42@<hash8>` resolve
identically across server and browser.

NULL columns serialize as JSON `null`; member ordering is deterministic
(`display_order ASC, id ASC`). Snapshot and peak_display are read as their
already-canonicalized JSON strings — re-parsing them through `JSON3.read`
+ re-emit would normally re-key them, but we route through
`canonical_json` which alphabetically sorts at every nesting level so the
hash is invariant across re-parse cycles.

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
        # is structural (a nested object) rather than a quoted string.
        # `canonical_json` re-sorts keys, so re-parse drift is impossible.
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
    "sha256:" * bytes2hex(SHA.sha256(canonical_json(payload)))
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

"""
    current_content_hash(db, comparison_id) -> Union{String, Nothing}

Read the stored `content_hash` for an existing comparison. Returns `nothing`
if the comparison doesn't exist (the 409-conflict path uses this to detect
a freshly-created or already-deleted comparison vs a hash-mismatched one).
"""
function current_content_hash(db::SQLite.DB, comparison_id::Integer)::Union{String, Nothing}
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT content_hash FROM comparisons WHERE id = ?", [Int(comparison_id)]))
    isempty(rows) && return nothing
    ismissing(rows[1].content_hash) ? nothing : String(rows[1].content_hash)
end

"""
    is_author(db, comparison_id, user_id) -> Bool

Author check for the route gates. Returns `false` if the comparison doesn't
exist, if `user_id` is `nothing` (unauthenticated request), or if
`comparisons.created_by` is NULL (orphan author — see spec §Authorship) or
differs from `user_id`.
"""
function is_author(db::SQLite.DB, comparison_id::Integer,
                   user_id::Union{Integer, Nothing})::Bool
    user_id === nothing && return false
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT created_by FROM comparisons WHERE id = ?", [Int(comparison_id)]))
    isempty(rows) && return false
    ismissing(rows[1].created_by) && return false
    Int(rows[1].created_by) == Int(user_id)
end

"""
    compute_member_snapshot(db, exposure_id) -> Dict{Symbol, Any}

Server-side derivation of the per-member snapshot shape per spec
§"Derived analysis state and staleness". Returns:

```
Dict(
  :effective_peaks => [
    Dict(:id, :q, :intensity, :sharpness, :source) for each effective peak
  ],
  :confirmed_index => Dict(:id, :phase, :lattice_d, :r_squared, :ngc, :peak_ids) | nothing,
  :analysis_inputs_hash => String | nothing,
)
```

`source` is `"auto"` or `"manual"` (matches `GET /api/exposures/:id/peaks`);
manual peaks carry `intensity = nothing`.

`confirmed_index` is the highest-scored member of the active custom
`index_groups` row whose `r_squared >= CONFIRMED_INDEX_R2_GATE`. Returns
`nothing` if there's no custom group for the exposure, no member meets the
gate, or no exposure exists.

This helper is the source of truth for the dispatcher's
`comparison_created` fallback (when the client omits a snapshot for a new
member) AND the `is_member_stale` fresh side. Reads from the DB — never the
trace file directly — so it requires `auto_peaks` and `peak_curations` rows
to be in place. Add curations get `intensity = nothing` (matches
`get_peaks_for_exposure`).

NOTE: The peak-set union (auto − exclude ∪ add) is a parallel SQL
implementation of `effective_peaks(db, exposure_id, q, I)` from
`pipeline.jl` — same q-tolerance (`MAX(1e-6, |q|*0.001)`), same union
semantics. We keep the parallel implementation here rather than calling
`effective_peaks` directly so this function can run without trace file I/O
(load_dat) — important on the snapshot-create hot path. The contract is
pinned by the regression test "compute_member_snapshot agrees with
effective_peaks (auto + exclude + add)" in test_comparisons.jl; if you
change the q-tolerance or union semantics in either implementation, both
must move together or that test will fail.
"""
function compute_member_snapshot(db::SQLite.DB, exposure_id::Integer)::Dict{Symbol, Any}
    eid = Int(exposure_id)

    # Mirror routes_peaks.jl::GET /api/exposures/:id/peaks: union of auto
    # peaks (joined against exclude curations) with `add` curations. Auto
    # peaks suppressed when an exclude curation matches their q. Sorted by q.
    peak_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT a.id, a.q, a.intensity, a.sharpness, 'auto' AS source
           FROM auto_peaks a
           LEFT JOIN peak_curations c
               ON c.exposure_id = a.exposure_id
              AND c.kind = 'exclude'
              AND ABS(c.q - a.q) <= MAX(1e-6, ABS(a.q) * 0.001)
           WHERE a.exposure_id = ? AND c.q IS NULL
           UNION ALL
           SELECT id, q, NULL AS intensity, NULL AS sharpness, 'manual' AS source
           FROM peak_curations
           WHERE exposure_id = ? AND kind = 'add'
           ORDER BY q""", [eid, eid]))

    effective_peaks = Vector{Any}(undef, length(peak_rows))
    for (i, r) in enumerate(peak_rows)
        effective_peaks[i] = Dict{Symbol, Any}(
            :id        => Int(r.id),
            :q         => Float64(r.q),
            :intensity => ismissing(r.intensity) ? nothing : Float64(r.intensity),
            # NULL → 0.0 (NOT nothing). The client-side
            # `computeMemberSnapshot` already coerces `null → 0` (snapshot.ts:60)
            # so the JS `MemberSnapshotPeak.sharpness: number` type contract
            # holds. If the server emits `nothing`, a GET-fetched snapshot
            # rehashed client-side diverges from the locally-computed hash
            # and `content_hash` parity breaks for the chat-citation
            # `@comparison:42@<hash8>` resolver.
            :sharpness => ismissing(r.sharpness) ? 0.0 : Float64(r.sharpness),
            :source    => String(r.source),
        )
    end

    # confirmed_index: highest-scored member of the active custom group
    # that clears the R² gate. Reads from index_groups + index_group_members
    # joined against indices. `inputs_hash` tracking belongs to the staleness
    # signal, not this snapshot — we just record what was confirmed.
    confirmed_index = nothing
    confirmed_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT i.id, i.phase, i.basis, i.lattice_d, i.r_squared, i.score
           FROM index_groups g
           JOIN index_group_members m ON m.group_id = g.id
           JOIN indices i ON i.id = m.index_id
           WHERE g.exposure_id = ? AND g.kind = 'custom' AND g.active = 1
             AND i.r_squared IS NOT NULL AND i.r_squared >= ?
           ORDER BY i.score DESC NULLS LAST, i.id ASC
           LIMIT 1""", [eid, CONFIRMED_INDEX_R2_GATE]))
    if !isempty(confirmed_rows)
        ix = confirmed_rows[1]
        ix_id = Int(ix.id)
        peak_id_rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT peak_id FROM index_peaks
               WHERE index_id = ? ORDER BY ratio_position""", [ix_id]))
        peak_ids = [Int(r.peak_id) for r in peak_id_rows]
        phase_str = ismissing(ix.phase) ? "" : String(ix.phase)
        lattice_d = ismissing(ix.lattice_d) ? nothing : Float64(ix.lattice_d)
        r_squared = ismissing(ix.r_squared) ? nothing : Float64(ix.r_squared)
        ngc = _ngc_for_phase(phase_str, lattice_d)
        confirmed_index = Dict{Symbol, Any}(
            :id        => ix_id,
            :phase     => phase_str,
            :lattice_d => lattice_d,
            :r_squared => r_squared,
            :ngc       => ngc,
            :peak_ids  => peak_ids,
        )
    end

    Dict{Symbol, Any}(
        :effective_peaks      => effective_peaks,
        :confirmed_index      => confirmed_index,
        :analysis_inputs_hash => read_inputs_hash(db, eid),
    )
end

"""
    is_member_stale(db, member_row) -> Bool

Compare the snapshot's `analysis_inputs_hash` against the live exposure's
`analysis_inputs_hash`. Returns `false` for orphan members
(`exposure_id IS NULL`) — there's no live exposure to drift from. Returns
`true` whenever either side is present but doesn't match the other; both
NULLs are considered fresh.
"""
function is_member_stale(db::SQLite.DB, member_row)::Bool
    ismissing(member_row.exposure_id) && return false
    eid = Int(member_row.exposure_id)
    snapshot_hash = nothing
    if !ismissing(member_row.snapshot)
        try
            obj = JSON3.read(String(member_row.snapshot))
            if haskey(obj, :analysis_inputs_hash) && obj.analysis_inputs_hash !== nothing
                snapshot_hash = String(obj.analysis_inputs_hash)
            end
        catch
            # Malformed snapshot JSON: treat as stale so the user sees the
            # warning rather than silently rendering with bogus data.
            return true
        end
    end
    current_hash = read_inputs_hash(db, eid)
    snapshot_hash != current_hash
end

"""
    member_ids_for_comparison(db, comparison_id) -> Set{Int}

(Already defined above; restated here for the helper inventory.)
"""

"""
    fetch_comparison_with_members(db, comparison_id) -> Union{Dict, Nothing}

Full nested response shape used by `GET /api/comparisons/:id` and the 409
`current_state` body. Returns `nothing` if the comparison doesn't exist.

The response shape:

```
Dict(
  :id, :title, :description, :content_hash,
  :created_by, :created_at, :updated_at,
  :forked_from_id, :forked_at_hash,
  :forked_from_title => parent's title or nothing,
  :members => [Dict(:id, :comparison_id, :exposure_id, :display_order,
                    :band_height, :y_offset, :normalization,
                    :color_override, :label_override,
                    :q_window_min, :q_window_max,
                    :peak_display, :snapshot, :is_stale,
                    :created_by, :created_at), …],
)
```

`peak_display` and `snapshot` are returned as parsed JSON (not strings).
`is_stale` is computed per member.
"""
function fetch_comparison_with_members(db::SQLite.DB, comparison_id::Integer)
    cid = Int(comparison_id)
    cmp_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT id, title, description, content_hash, created_by,
                  created_at, updated_at, forked_from_id, forked_at_hash
           FROM comparisons WHERE id = ?""", [cid]))
    isempty(cmp_rows) && return nothing
    cmp = cmp_rows[1]

    # Look up the parent's title (if any) so the lineage badge can render
    # "Forked from <title>" without a second round-trip.
    forked_from_title = nothing
    if !ismissing(cmp.forked_from_id)
        parent_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT title FROM comparisons WHERE id = ?", [Int(cmp.forked_from_id)]))
        if !isempty(parent_rows) && !ismissing(parent_rows[1].title)
            forked_from_title = String(parent_rows[1].title)
        end
    end

    member_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT id, comparison_id, exposure_id, display_order,
                  band_height, y_offset, normalization,
                  color_override, label_override,
                  q_window_min, q_window_max,
                  peak_display, snapshot, created_by, created_at
           FROM comparison_members
           WHERE comparison_id = ?
           ORDER BY display_order ASC, id ASC""", [cid]))

    members = Vector{Any}(undef, length(member_rows))
    for (i, m) in enumerate(member_rows)
        snap_obj = ismissing(m.snapshot) ? nothing : JSON3.read(String(m.snapshot))
        peak_obj = ismissing(m.peak_display) ? nothing : JSON3.read(String(m.peak_display))
        members[i] = Dict{Symbol, Any}(
            :id             => Int(m.id),
            :comparison_id  => Int(m.comparison_id),
            :exposure_id    => ismissing(m.exposure_id)    ? nothing : Int(m.exposure_id),
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
            :is_stale       => is_member_stale(db, m),
            :created_by     => ismissing(m.created_by)     ? nothing : Int(m.created_by),
            :created_at     => ismissing(m.created_at)     ? nothing : String(m.created_at),
        )
    end

    Dict{Symbol, Any}(
        :id              => Int(cmp.id),
        :title           => ismissing(cmp.title) ? "" : String(cmp.title),
        :description     => ismissing(cmp.description) ? nothing : String(cmp.description),
        :content_hash    => ismissing(cmp.content_hash) ? "" : String(cmp.content_hash),
        :created_by      => ismissing(cmp.created_by) ? nothing : Int(cmp.created_by),
        :created_at      => ismissing(cmp.created_at) ? nothing : String(cmp.created_at),
        :updated_at      => ismissing(cmp.updated_at) ? nothing : String(cmp.updated_at),
        :forked_from_id  => ismissing(cmp.forked_from_id) ? nothing : Int(cmp.forked_from_id),
        :forked_at_hash  => ismissing(cmp.forked_at_hash) ? nothing : String(cmp.forked_at_hash),
        :forked_from_title => forked_from_title,
        :members         => members,
    )
end

"""
    recently_used_exposures(db, user_id; limit=20) -> Vector{Int}

Picker "Recently used" derivation per spec §"Recently used" derivation.
Queries `comparison_members` joined against the user's pick history
(`created_by` = `user_id`), grouped by `exposure_id`, sorted by latest pick
time. Excludes orphans (`exposure_id IS NULL`).
"""
function recently_used_exposures(db::SQLite.DB, user_id::Integer; limit::Int = 20)::Vector{Int}
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT exposure_id, MAX(created_at) AS last_picked
           FROM comparison_members
           WHERE created_by = ? AND exposure_id IS NOT NULL
           GROUP BY exposure_id
           ORDER BY last_picked DESC
           LIMIT ?""", [Int(user_id), Int(limit)]))
    [Int(r.exposure_id) for r in rows]
end

"""
    picker_samples(db, experiment_id) -> Vector{Dict{Symbol, Any}}

Picker primary list per spec §"PR1 — sample-first picker → Backend".

For each sample in `experiment_id`, returns:
  :sample              => sample row (Symbol-keyed Dict, with :tags vector)
  :indexing_exposure_id => Int or nothing — `selected = 1` else MAX(id) else nothing
  :all_exposures       => Vector of {:id, :filename, :selected::Bool}, ORDER BY id ASC

Three bulk queries (no JOIN'd Cartesian flatten, no per-sample N+1):
(1) `WHERE experiment_id = ?` for samples,
(2) `WHERE sample_id IN (...)` for exposures,
(3) `WHERE sample_id IN (...)` for sample_tags.
Empty experiment ⇒ [].
"""
function picker_samples(db::SQLite.DB, experiment_id::Integer)::Vector{Dict{Symbol, Any}}
    # Explicit column list (PR #96 review): keep the picker JSON shape
    # deliberate so a future column added to `samples` doesn't auto-leak
    # into the picker payload.
    samples = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, experiment_id, name, label, notes
         FROM samples WHERE experiment_id = ? ORDER BY id",
        [Int(experiment_id)]))
    isempty(samples) && return Dict{Symbol, Any}[]

    sample_ids   = [Int(s.id) for s in samples]
    placeholders = join(fill("?", length(sample_ids)), ",")
    exposures    = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, sample_id, filename, selected FROM exposures
         WHERE sample_id IN ($placeholders) ORDER BY sample_id ASC, id ASC",
        sample_ids))
    tag_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, sample_id, key, value, source FROM sample_tags
         WHERE sample_id IN ($placeholders) ORDER BY sample_id ASC, id ASC",
        sample_ids))

    # Group exposures + tags by sample_id in Julia (avoids a Cartesian JOIN dedup).
    grouped_exps = Dict{Int, Vector{NamedTuple}}()
    for e in exposures
        push!(get!(grouped_exps, Int(e.sample_id), NamedTuple[]), e)
    end
    grouped_tags = Dict{Int, Vector{NamedTuple}}()
    for t in tag_rows
        push!(get!(grouped_tags, Int(t.sample_id), NamedTuple[]), t)
    end

    out = Dict{Symbol, Any}[]
    for sm in samples
        sid  = Int(sm.id)
        exps = get(grouped_exps, sid, NamedTuple[])

        # Resolve indexing exposure: highest-id selected wins (defensive
        # against legacy multi-selected data); else highest-id overall;
        # else nothing for zero-exposure samples. `exps` is sorted id ASC,
        # so iterating in reverse yields highest-id first.
        idx_id = nothing
        for e in Iterators.reverse(exps)
            if e.selected != 0
                idx_id = Int(e.id); break
            end
        end
        if idx_id === nothing && !isempty(exps)
            idx_id = Int(last(exps).id)   # last == max by id ASC ordering
        end

        # Sample → row_to_json + bulk-grouped tags. Drop sample_id from each
        # tag dict (it's redundant once grouped under the sample).
        tags_for_sm = get(grouped_tags, sid, NamedTuple[])
        tag_dicts = map(tags_for_sm) do t
            d = row_to_json(t)
            delete!(d, :sample_id)
            d
        end
        sample_dict = row_to_json(sm)
        sample_dict[:tags] = tag_dicts

        all_exp = [row_to_json(e; bool_keys = (:selected,)) for e in exps]
        push!(out, Dict{Symbol, Any}(
            :sample               => sample_dict,
            :indexing_exposure_id => idx_id,
            :all_exposures        => all_exp,
        ))
    end
    out
end

"""
    comparisons_for_experiment(db, experiment_id) -> Vector{Dict}

Per-experiment listing per spec §REST API. Returns the comparisons that
have at least one member whose exposure → sample → experiment chain
points at `experiment_id`. Sorted by latest `user_actions.created_at` for
that comparison (so newly-edited rows float to the top).
"""
function comparisons_for_experiment(db::SQLite.DB, experiment_id::Integer)::Vector{Dict{Symbol, Any}}
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT DISTINCT c.id, c.title, c.description, c.content_hash,
                  c.created_by, c.created_at, c.updated_at,
                  c.forked_from_id, c.forked_at_hash,
                  COALESCE((SELECT MAX(ua.timestamp) FROM user_actions ua
                            WHERE ua.entity_type = 'comparison'
                              AND ua.entity_id = c.id),
                           c.updated_at) AS last_event_at
           FROM comparisons c
           JOIN comparison_members cm ON cm.comparison_id = c.id
           JOIN exposures e ON e.id = cm.exposure_id
           JOIN samples s ON s.id = e.sample_id
           WHERE s.experiment_id = ?
           ORDER BY last_event_at DESC, c.id DESC""", [Int(experiment_id)]))
    _comparison_listing_rows(rows)
end

"""
    comparisons_listing(db) -> Vector{Dict}

Global listing per `GET /api/comparisons`. Sorted by latest event time.
Includes orphan comparisons (members with NULL exposure_id) and
zero-member comparisons (which can exist transiently between create and
member insert; defensive).
"""
function comparisons_listing(db::SQLite.DB)::Vector{Dict{Symbol, Any}}
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT c.id, c.title, c.description, c.content_hash,
                  c.created_by, c.created_at, c.updated_at,
                  c.forked_from_id, c.forked_at_hash,
                  COALESCE((SELECT MAX(ua.timestamp) FROM user_actions ua
                            WHERE ua.entity_type = 'comparison'
                              AND ua.entity_id = c.id),
                           c.updated_at) AS last_event_at
           FROM comparisons c
           ORDER BY last_event_at DESC, c.id DESC"""))
    _comparison_listing_rows(rows)
end

# Helper: shared listing-row projection. Returns the lightweight per-row
# shape used by both `comparisons_listing` and `comparisons_for_experiment`
# (no member nesting — listing rows are summaries; clients fetch
# `GET /api/comparisons/:id` for details).
function _comparison_listing_rows(rows)::Vector{Dict{Symbol, Any}}
    out = Vector{Dict{Symbol, Any}}(undef, length(rows))
    for (i, r) in enumerate(rows)
        out[i] = Dict{Symbol, Any}(
            :id              => Int(r.id),
            :title           => ismissing(r.title) ? "" : String(r.title),
            :description     => ismissing(r.description) ? nothing : String(r.description),
            :content_hash    => ismissing(r.content_hash) ? "" : String(r.content_hash),
            :created_by      => ismissing(r.created_by) ? nothing : Int(r.created_by),
            :created_at      => ismissing(r.created_at) ? nothing : String(r.created_at),
            :updated_at      => ismissing(r.updated_at) ? nothing : String(r.updated_at),
            :forked_from_id  => ismissing(r.forked_from_id) ? nothing : Int(r.forked_from_id),
            :forked_at_hash  => ismissing(r.forked_at_hash) ? nothing : String(r.forked_at_hash),
        )
    end
    out
end

"""
    forks_of_comparison(db, comparison_id) -> Vector{Dict}

List of comparisons whose `forked_from_id` points at this id. Returns the
same row shape as `comparisons_listing`.
"""
function forks_of_comparison(db::SQLite.DB, comparison_id::Integer)::Vector{Dict{Symbol, Any}}
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT c.id, c.title, c.description, c.content_hash,
                  c.created_by, c.created_at, c.updated_at,
                  c.forked_from_id, c.forked_at_hash,
                  COALESCE((SELECT MAX(ua.timestamp) FROM user_actions ua
                            WHERE ua.entity_type = 'comparison'
                              AND ua.entity_id = c.id),
                           c.updated_at) AS last_event_at
           FROM comparisons c
           WHERE c.forked_from_id = ?
           ORDER BY last_event_at DESC, c.id DESC""", [Int(comparison_id)]))
    _comparison_listing_rows(rows)
end
