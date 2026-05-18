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
