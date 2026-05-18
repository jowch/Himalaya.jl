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
                  ) AS has_stale_members
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
        )
    end
    out
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
