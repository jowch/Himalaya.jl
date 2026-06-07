using SQLite, DBInterface, Tables, JSON3, SHA

# ─────────────────────────────────────────────────────────────────────────────
# Series business logic (I2.2, #165) — adapted from `comparisons.jl`.
#
# `series.jl` and `comparisons.jl` are both `include`d into the `HimalayaUI`
# module, so this file reuses the in-module generic helpers `canonical_json`,
# `user_id_for_event`, `comparison_now_iso`, `compute_member_snapshot`, and
# `is_member_stale` directly rather than duplicating them. Phase 3 (#175 / I3.6)
# relocates those when it deletes `comparisons.jl` — out of scope here.
#
# Two departures from the comparison originals: there is no `is_author` gate,
# and the `last_event_at` mixed-timestamp sort bug (#76) is fixed with a SQL
# `datetime()` wrapper in `series_listing` / `forks_of_series`.
# ─────────────────────────────────────────────────────────────────────────────

"""
    series_listing(db) -> Vector{Dict}

Corpus-wide listing for `GET /api/series`. Includes zero-member and
orphan-member series defensively. The `last_event_at` projection wraps the
coalesced timestamp in SQLite `datetime()` so the space-separated
`user_actions.timestamp` and the `T`-separated `series.updated_at` normalise
to one comparable form — fixing sort bug #76 rather than copying it from
`comparisons.jl:669`. The projected value is the normalised form, so it is a
valid client sort key, not merely a display hint.
"""
function series_listing(db::SQLite.DB)::Vector{Dict{Symbol, Any}}
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT s.id, s.title, s.description, s.content_hash,
                  s.created_by, s.created_at, s.updated_at,
                  s.forked_from_id, s.forked_at_hash,
                  s.view_grouping_mode, s.view_show_peak_ticks, s.view_show_peak_labels,
                  s.ordering_variable,
                  datetime(COALESCE((SELECT MAX(ua.timestamp) FROM user_actions ua
                                     WHERE ua.entity_type = 'series'
                                       AND ua.entity_id = s.id), s.updated_at))
                      AS last_event_at,
                  u.username AS author_username,
                  (SELECT COUNT(*) FROM series_members sm
                   WHERE sm.series_id = s.id) AS member_count,
                  (SELECT GROUP_CONCAT(json_extract(sm.snapshot, '\$.confirmed_index.phase')
                                       || '#' || sm.display_order, '|')
                   FROM series_members sm
                   WHERE sm.series_id = s.id
                     AND json_extract(sm.snapshot, '\$.confirmed_index.phase') IS NOT NULL)
                      AS member_phases_concat,
                  EXISTS (
                    SELECT 1 FROM series_members sm
                    JOIN exposures e ON e.id = sm.exposure_id
                    WHERE sm.series_id = s.id
                      AND sm.exposure_id IS NOT NULL
                      AND json_extract(sm.snapshot, '\$.analysis_inputs_hash')
                          IS NOT e.analysis_inputs_hash
                  ) AS has_stale_members,
                  (SELECT COUNT(DISTINCT sa.experiment_id) > 1
                   FROM series_members sm
                   JOIN exposures ex ON ex.id = sm.exposure_id
                   JOIN samples sa   ON sa.id = ex.sample_id
                   WHERE sm.series_id = s.id) AS spans_experiments
           FROM series s
           LEFT JOIN users u ON u.id = s.created_by
           ORDER BY last_event_at DESC, s.id DESC"""))
    _series_listing_rows(rows)
end

# Lightweight per-row listing projection — the shape `series_listing` and
# `forks_of_series` return. Adapted from `_comparison_listing_rows`. Reuses the
# in-module `_topk_phases` / `_count_distinct_phases` phase-token helpers.
function _series_listing_rows(rows)::Vector{Dict{Symbol, Any}}
    out = Vector{Dict{Symbol, Any}}(undef, length(rows))
    for (i, r) in enumerate(rows)
        phases_str = ismissing(r.member_phases_concat) ? "" : String(r.member_phases_concat)
        member_phases = _topk_phases(phases_str, 3)
        out[i] = Dict{Symbol, Any}(
            :id                    => Int(r.id),
            :title                 => ismissing(r.title) ? "" : String(r.title),
            :description           => ismissing(r.description) ? nothing : String(r.description),
            :content_hash          => ismissing(r.content_hash) ? "" : String(r.content_hash),
            :created_by            => ismissing(r.created_by) ? nothing : Int(r.created_by),
            :created_at            => ismissing(r.created_at) ? nothing : String(r.created_at),
            :updated_at            => ismissing(r.updated_at) ? nothing : String(r.updated_at),
            :forked_from_id        => ismissing(r.forked_from_id) ? nothing : Int(r.forked_from_id),
            :forked_at_hash        => ismissing(r.forked_at_hash) ? nothing : String(r.forked_at_hash),
            :view_grouping_mode    => ismissing(r.view_grouping_mode) ? nothing : String(r.view_grouping_mode),
            :view_show_peak_ticks  => ismissing(r.view_show_peak_ticks) ? nothing : Bool(r.view_show_peak_ticks),
            :view_show_peak_labels => ismissing(r.view_show_peak_labels) ? nothing : Bool(r.view_show_peak_labels),
            :last_event_at         => ismissing(r.last_event_at) ? nothing : String(r.last_event_at),
            :author_username       => ismissing(r.author_username) ? nothing : String(r.author_username),
            :member_count          => Int(r.member_count),
            :member_phases         => member_phases,
            :member_phase_count    => _count_distinct_phases(phases_str),
            :has_stale_members     => Bool(r.has_stale_members),
            :ordering_variable     => ismissing(r.ordering_variable) ? nothing : String(r.ordering_variable),
            # Cross-experiment = members resolve to >1 distinct samples.experiment_id.
            # Valid because q is absolute (Å⁻¹); see redesign-notes architecture decision 1.
            :spans_experiments     => !ismissing(r.spans_experiments) && Bool(r.spans_experiments),
        )
    end
    out
end

"""
    forks_of_series(db, series_id) -> Vector{Dict}

Series whose `forked_from_id` points at this id. Same row shape as
`series_listing`; same `datetime()` `last_event_at` fix (#76).
"""
function forks_of_series(db::SQLite.DB, series_id::Integer)::Vector{Dict{Symbol, Any}}
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT s.id, s.title, s.description, s.content_hash,
                  s.created_by, s.created_at, s.updated_at,
                  s.forked_from_id, s.forked_at_hash,
                  s.view_grouping_mode, s.view_show_peak_ticks, s.view_show_peak_labels,
                  s.ordering_variable,
                  datetime(COALESCE((SELECT MAX(ua.timestamp) FROM user_actions ua
                                     WHERE ua.entity_type = 'series'
                                       AND ua.entity_id = s.id), s.updated_at))
                      AS last_event_at,
                  u.username AS author_username,
                  (SELECT COUNT(*) FROM series_members sm
                   WHERE sm.series_id = s.id) AS member_count,
                  (SELECT GROUP_CONCAT(json_extract(sm.snapshot, '\$.confirmed_index.phase')
                                       || '#' || sm.display_order, '|')
                   FROM series_members sm
                   WHERE sm.series_id = s.id
                     AND json_extract(sm.snapshot, '\$.confirmed_index.phase') IS NOT NULL)
                      AS member_phases_concat,
                  EXISTS (
                    SELECT 1 FROM series_members sm
                    JOIN exposures e ON e.id = sm.exposure_id
                    WHERE sm.series_id = s.id
                      AND sm.exposure_id IS NOT NULL
                      AND json_extract(sm.snapshot, '\$.analysis_inputs_hash')
                          IS NOT e.analysis_inputs_hash
                  ) AS has_stale_members,
                  (SELECT COUNT(DISTINCT sa.experiment_id) > 1
                   FROM series_members sm
                   JOIN exposures ex ON ex.id = sm.exposure_id
                   JOIN samples sa   ON sa.id = ex.sample_id
                   WHERE sm.series_id = s.id) AS spans_experiments
           FROM series s
           LEFT JOIN users u ON u.id = s.created_by
           WHERE s.forked_from_id = ?
           ORDER BY last_event_at DESC, s.id DESC""", [Int(series_id)]))
    _series_listing_rows(rows)
end

"""
    fetch_series_with_plate(db, series_id) -> Union{Dict, Nothing}

Full nested response shape for `GET /api/series/:id` and the `series_plate_committed`
`post_state` envelope. Returns `nothing` if the series does not exist.

The series row carries the recipe columns (`ordering_variable`, `order_rule`,
`state`); `:members` is the plate (`series_members`, frozen per-exposure
snapshots) ordered by `display_order`; `:samples` is the recipe
(`series_samples`) ordered by `position`. `peak_display` / `snapshot` are
returned as parsed JSON; `is_stale` is computed per member.
"""
function fetch_series_with_plate(db::SQLite.DB, series_id::Integer)
    sid = Int(series_id)
    s_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT id, title, description, content_hash, created_by,
                  created_at, updated_at, forked_from_id, forked_at_hash,
                  view_grouping_mode, view_show_peak_ticks, view_show_peak_labels,
                  ordering_variable, order_rule, state
           FROM series WHERE id = ?""", [sid]))
    isempty(s_rows) && return nothing
    s = s_rows[1]

    # Parent title for the lineage badge (one extra round-trip avoided).
    forked_from_title = nothing
    if !ismissing(s.forked_from_id)
        parent_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT title FROM series WHERE id = ?", [Int(s.forked_from_id)]))
        if !isempty(parent_rows) && !ismissing(parent_rows[1].title)
            forked_from_title = String(parent_rows[1].title)
        end
    end

    member_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT id, series_id, exposure_id, display_order,
                  band_height, y_offset, normalization,
                  color_override, label_override,
                  q_window_min, q_window_max,
                  peak_display, snapshot, created_by, created_at
           FROM series_members
           WHERE series_id = ?
           ORDER BY display_order ASC, id ASC""", [sid]))
    members = Vector{Any}(undef, length(member_rows))
    for (i, m) in enumerate(member_rows)
        snap_obj = ismissing(m.snapshot) ? nothing : JSON3.read(String(m.snapshot))
        peak_obj = ismissing(m.peak_display) ? nothing : JSON3.read(String(m.peak_display))
        members[i] = Dict{Symbol, Any}(
            :id             => Int(m.id),
            :series_id      => Int(m.series_id),
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

    sample_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT id, series_id, sample_id, position, pinned, excluded
           FROM series_samples
           WHERE series_id = ?
           ORDER BY position ASC, id ASC""", [sid]))
    samples = [Dict{Symbol, Any}(
        :id        => Int(r.id),
        :series_id => Int(r.series_id),
        :sample_id => Int(r.sample_id),
        :position  => Int(r.position),
        :pinned    => Bool(r.pinned),
        :excluded  => Bool(r.excluded),
    ) for r in sample_rows]

    Dict{Symbol, Any}(
        :id                    => Int(s.id),
        :title                 => ismissing(s.title) ? "" : String(s.title),
        :description           => ismissing(s.description) ? nothing : String(s.description),
        :content_hash          => ismissing(s.content_hash) ? "" : String(s.content_hash),
        :created_by            => ismissing(s.created_by) ? nothing : Int(s.created_by),
        :created_at            => ismissing(s.created_at) ? nothing : String(s.created_at),
        :updated_at            => ismissing(s.updated_at) ? nothing : String(s.updated_at),
        :forked_from_id        => ismissing(s.forked_from_id) ? nothing : Int(s.forked_from_id),
        :forked_at_hash        => ismissing(s.forked_at_hash) ? nothing : String(s.forked_at_hash),
        :forked_from_title     => forked_from_title,
        :view_grouping_mode    => ismissing(s.view_grouping_mode) ? nothing : String(s.view_grouping_mode),
        :view_show_peak_ticks  => ismissing(s.view_show_peak_ticks) ? nothing : Bool(s.view_show_peak_ticks),
        :view_show_peak_labels => ismissing(s.view_show_peak_labels) ? nothing : Bool(s.view_show_peak_labels),
        :ordering_variable     => ismissing(s.ordering_variable) ? nothing : String(s.ordering_variable),
        :order_rule            => String(s.order_rule),
        :state                 => String(s.state),
        :members               => members,
        :samples               => samples,
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Mutating-route helpers (I2.2, Task 5)
# ─────────────────────────────────────────────────────────────────────────────

"""
    series_exists(db, series_id) -> Bool

The 404 existence probe for the mutating routes. A dedicated probe is required
because a draft series carries `content_hash IS NULL` by design, so a NULL hash
cannot distinguish "missing" from "draft" — see `current_series_content_hash`.
"""
function series_exists(db::SQLite.DB, series_id::Integer)::Bool
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 AS one FROM series WHERE id = ?", [Int(series_id)]))
    !isempty(rows)
end

"""
    current_series_content_hash(db, series_id) -> Union{String, Nothing}

The stored `content_hash`. Returns `nothing` for a missing series AND for an
uncommitted draft (drafts have NULL `content_hash`). Used only for the `commit`
409 optimistic-concurrency check — never as the existence probe (that is
`series_exists`).
"""
function current_series_content_hash(db::SQLite.DB, series_id::Integer)::Union{String, Nothing}
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT content_hash FROM series WHERE id = ?", [Int(series_id)]))
    isempty(rows) && return nothing
    ismissing(rows[1].content_hash) ? nothing : String(rows[1].content_hash)
end

"""
    compute_series_content_hash(db, series_id) -> String

`sha256:`-prefixed hash of the series **plate** — title, description, and the
`series_members` rows. The recipe (`series_samples`) is deliberately excluded
(master plan §5.1): `content_hash` reflects the committed plate only, so
`series_recipe_updated` never touches it.
"""
function compute_series_content_hash(db::SQLite.DB, series_id::Integer)::String
    s_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT title, description FROM series WHERE id = ?", [Int(series_id)]))
    isempty(s_rows) && error("compute_series_content_hash: series $series_id not found")
    s = s_rows[1]
    # series.title is nullable (a degenerate draft placeholder can carry NULL);
    # "" is the sentinel. compute_content_hash for comparisons has no guard
    # here because comparison titles are always set.
    title = ismissing(s.title) ? "" : String(s.title)
    description = ismissing(s.description) ? nothing : String(s.description)

    member_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT id, exposure_id, display_order, band_height, y_offset,
                  normalization, color_override, label_override,
                  q_window_min, q_window_max, peak_display, snapshot
           FROM series_members
           WHERE series_id = ?
           ORDER BY display_order ASC, id ASC""", [Int(series_id)]))
    members = Vector{Any}(undef, length(member_rows))
    for (i, m) in enumerate(member_rows)
        # Already-canonical JSON columns are re-parsed so the encoded form
        # is structural (a nested object) rather than a quoted string.
        # `canonical_json` re-sorts keys, so re-parse drift is impossible.
        snap_obj = ismissing(m.snapshot) ? nothing : JSON3.read(String(m.snapshot))
        peak_obj = ismissing(m.peak_display) ? nothing : JSON3.read(String(m.peak_display))
        members[i] = Dict{Symbol,Any}(
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
        )
    end
    payload = Dict{Symbol,Any}(
        :title       => title,
        :description => description,
        :members     => members,
    )
    "sha256:" * bytes2hex(SHA.sha256(canonical_json(payload)))
end

# ─────────────────────────────────────────────────────────────────────────────
# Mutating-route helpers (I2.2, Task 6)
# ─────────────────────────────────────────────────────────────────────────────

"""
    _series_sample_payload(m_in, default_position) -> Dict{Symbol, Any}

Normalize one recipe entry from a request body into the `series_samples`
payload shape. `sample_id` is required (a recipe row with no target is
unrenderable). `position` defaults to `default_position` — the caller passes
the entry's 0-based index, so a `samples` array with omitted positions yields
sequential positions rather than every entry colliding on `position = 0`
(`UNIQUE(series_id, position)`). `pinned` / `excluded` default to `false`.
"""
function _series_sample_payload(m_in, default_position::Integer)
    Dict{Symbol, Any}(
        :sample_id => Int(m_in[:sample_id]),
        :position  => haskey(m_in, :position) ? Int(m_in[:position]) : Int(default_position),
        :pinned    => haskey(m_in, :pinned)   ? Bool(m_in[:pinned])   : false,
        :excluded  => haskey(m_in, :excluded) ? Bool(m_in[:excluded]) : false,
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Batch trace loader (Phase-4, 4a) — unblocks the Series-folio CardFigure.
# ─────────────────────────────────────────────────────────────────────────────

"""
    series_member_traces(db, series_id) -> Dict{Int,Any}

Resolve and load every member trace of `series_id`, keyed by exposure_id, in
display order. Members with no exposure, a derived (non-"file") exposure, no
filename, or a missing `.dat` on disk are SKIPPED (omitted from the map) — the
batch route degrades gracefully; the folio's `toWaterfallRows` renders an empty
row for any absent exposure. `config_from_db` is memoized per experiment.

Note the deliberate skip-vs-throw asymmetry: a member whose `.dat` is *absent*
is skipped, but a member whose `.dat` is *present-but-corrupt* (unparseable /
<3 columns) intentionally surfaces as a 500 via `load_dat` — corruption is NOT
swallowed (it should fail loudly, not silently vanish from the map).
"""
function series_member_traces(db::SQLite.DB, series_id::Integer)
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT e.id AS exposure_id, e.filename, e.kind,
                  x.id AS experiment_id, x.analysis_dir
           FROM series_members sm
           JOIN exposures e   ON e.id = sm.exposure_id
           JOIN samples s     ON s.id = e.sample_id
           JOIN experiments x ON x.id = s.experiment_id
           WHERE sm.series_id = ? AND sm.exposure_id IS NOT NULL
           ORDER BY sm.display_order ASC, sm.id ASC""", [series_id]))

    out = Dict{Int,Any}()
    cfg_cache = Dict{Int,Any}()
    for row in rows
        row.kind === missing && continue
        String(row.kind) == "file" || continue
        row.filename === missing && continue
        eid = Int(row.experiment_id)
        cfg = get!(cfg_cache, eid) do
            config_from_db(db, eid)
        end
        pattern = replace(cfg.integration_pattern, "{name}" => String(row.filename))
        path    = joinpath(String(row.analysis_dir), pattern)
        isfile(path) || continue
        q, I, σ = load_dat(path)
        out[Int(row.exposure_id)] = Dict(:q => q, :I => I, :sigma => σ)
    end
    return out
end

"""
    _series_member_payload(db, m_in) -> Dict{Symbol, Any}

Normalize one plate-member entry from a request body into the payload shape the
`series_plate_committed` dispatcher expects. Fills `snapshot` from
`compute_member_snapshot(db, …)` when the client omitted it AND an
`exposure_id` is present; an orphan member (no `exposure_id`) with no snapshot
gets a minimal one so the `NOT NULL` `json_valid` CHECK on
`series_members.snapshot` is satisfied. Adapted verbatim from
`_comparison_member_payload` (the `series_members` and `comparison_members`
column shapes are identical).
"""
function _series_member_payload(db::SQLite.DB, m_in)
    eid = haskey(m_in, :exposure_id) ? m_in[:exposure_id] : nothing
    snap_in = haskey(m_in, :snapshot) ? m_in[:snapshot] : nothing
    snap = if snap_in !== nothing
        snap_in
    elseif eid !== nothing
        compute_member_snapshot(db, Int(eid))
    else
        # Orphan member (no exposure_id) — minimal stub so the NOT NULL
        # json_valid CHECK on series_members.snapshot is satisfied.
        Dict{Symbol, Any}(
            :effective_peaks      => Any[],
            :confirmed_index      => nothing,
            :analysis_inputs_hash => nothing,
        )
    end
    Dict{Symbol, Any}(
        :id             => haskey(m_in, :id) ? m_in[:id] : nothing,
        :exposure_id    => eid,
        :display_order  => haskey(m_in, :display_order) ? Int(m_in[:display_order]) : 0,
        :band_height    => haskey(m_in, :band_height)   ? Float64(m_in[:band_height]) : 1.0,
        :y_offset       => haskey(m_in, :y_offset)      ? Float64(m_in[:y_offset])    : 0.0,
        :normalization  => haskey(m_in, :normalization) ? String(m_in[:normalization]) : "none",
        :color_override => haskey(m_in, :color_override) ? m_in[:color_override] : nothing,
        :label_override => haskey(m_in, :label_override) ? m_in[:label_override] : nothing,
        :q_window_min   => haskey(m_in, :q_window_min)   ? m_in[:q_window_min]   : nothing,
        :q_window_max   => haskey(m_in, :q_window_max)   ? m_in[:q_window_max]   : nothing,
        :peak_display   => haskey(m_in, :peak_display)   ? m_in[:peak_display]   : nothing,
        :snapshot       => snap,
    )
end
