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

Corpus-wide listing for `GET /api/series`. (Filled in by Task 2.)
"""
function series_listing(db::SQLite.DB)::Vector{Dict{Symbol, Any}}
    return Dict{Symbol, Any}[]
end
