# Series REST Routes + Business Logic Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `series.jl` (business logic) and `routes_series.jl` (REST routes) so the series data model has a full, frozen route surface — shipping invisibly, with empty `series` tables and no UI consumer.

**Architecture:** `routes_series.jl` mirrors `routes_comparisons.jl`; `series.jl` adapts logic from `comparisons.jl`. Both are `include`d into the `HimalayaUI` module, so they reuse in-module helpers (`canonical_json`, `user_id_for_event`, `comparison_now_iso`, `compute_member_snapshot`, `is_member_stale`, `_json_error`, `_view_fields_error`, `_topk_phases`, `_count_distinct_phases`) directly. Two departures from the Compare-page originals: the `is_author` gate is dropped, and the `last_event_at` mixed-timestamp sort bug (#76) is fixed with a SQL `datetime()` wrapper. The series mutating routes emit events whose dispatcher branches do not exist until #166–#168; `update_view_for_event!` silently no-ops unknown kinds, so those routes are written in full but produce no view rows until then.

**Tech Stack:** Julia, Oxygen.jl (routing), SQLite.jl, JSON3.jl, stdlib `Test`.

**Spec:** `docs/superpowers/specs/2026-05-18-series-routes-design.md` — read it before starting.

---

## Reference files (read before starting)

- `packages/HimalayaUI/src/comparisons.jl` — the business-logic original being adapted.
- `packages/HimalayaUI/src/routes_comparisons.jl` — the route original being mirrored.
- `packages/HimalayaUI/src/db.jl:723` — `migrate_series!`, the series schema.
- `packages/HimalayaUI/src/events.jl:305` — `update_view_for_event!` (no-ops unknown kinds at line 426).
- `packages/HimalayaUI/test/test_routes_comparisons.jl` — the route-test pattern (`with_test_server`).

## File structure

- **Create** `packages/HimalayaUI/src/series.jl` — business logic: `series_listing`, `_series_listing_rows`, `forks_of_series`, `fetch_series_with_plate`, `compute_series_content_hash`, `current_series_content_hash`, `series_exists`, `_series_member_payload`, `_series_sample_payload`.
- **Create** `packages/HimalayaUI/src/routes_series.jl` — `register_series_routes!()` with all 12 routes.
- **Create** `packages/HimalayaUI/test/test_routes_series.jl` — Julia route + logic tests.
- **Modify** `packages/HimalayaUI/src/HimalayaUI.jl` — add two `include` lines.
- **Modify** `packages/HimalayaUI/src/server.jl` — add `register_series_routes!()` to `register_routes!`.
- **Modify** `packages/HimalayaUI/test/runtests.jl` — add `include("test_routes_series.jl")`.

## Conventions for this plan

- All commands run from `packages/HimalayaUI/` unless stated otherwise.
- **Per-task test run** (fast — just the new file plus its server helper):
  ```bash
  julia --project=. -e 'using HimalayaUI, Test, HTTP, JSON3, SQLite, DBInterface, Tables; include("test/test_http.jl"); include("test/test_routes_series.jl")'
  ```
  `test/test_http.jl` defines `with_test_server`; including it first is required for the server-backed tests. Tasks 2–5 test pure functions and do not need the server, but the command above works for every task.
- Internal (non-exported) functions are reached as `HimalayaUI.<name>` in tests.
- Commit after every task. End every commit message with:
  ```
  Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
  ```

---

## Task 1: Scaffold the files, wiring, and `GET /api/series`

**Files:**
- Create: `packages/HimalayaUI/src/series.jl`
- Create: `packages/HimalayaUI/src/routes_series.jl`
- Create: `packages/HimalayaUI/test/test_routes_series.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl`
- Modify: `packages/HimalayaUI/src/server.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/test_routes_series.jl`:

```julia
using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# ─────────────────────────────────────────────────────────────────────────────
# I2.2 (#165) — series REST routes + business logic.
#
# The mutating routes emit events whose dispatcher branches do not land until
# #166–#168; `update_view_for_event!` no-ops unknown kinds, so those routes
# write a `user_actions` row but no view rows. Behavioural round-trips for the
# mutating routes belong to #166–#168 and #170; this file covers GET
# round-trips, the `last_event_at` sort fix, the messages round-trip (works
# today), and HTTP-level smoke tests for the mutating routes.
# ─────────────────────────────────────────────────────────────────────────────

# A DB with one experiment / sample / exposure — a valid FK target for
# `series_members.exposure_id` and `series_samples.sample_id`.
function _series_test_db(tmp::String)
    db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    DBInterface.execute(db, """INSERT INTO experiments
        (id, name, path, data_dir, analysis_dir)
        VALUES (10, 'exp', '/x', '/x/d', '/x/a')""")
    DBInterface.execute(db,
        "INSERT INTO samples (id, experiment_id, name) VALUES (100, 10, 'sA')")
    DBInterface.execute(db,
        "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (1000, 100, 'JC001', 1)")
    db
end

@testset "Series routes" begin

    @testset "GET /api/series — empty corpus" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                resp = HTTP.get("$base/api/series", ["X-Username" => "alice"])
                @test resp.status == 200
                @test JSON3.read(resp.body) == []
            end
            close(db)
        end
    end

end
```

- [ ] **Step 2: Run the test to verify it fails**

Run (from `packages/HimalayaUI/`):
```bash
julia --project=. -e 'using HimalayaUI, Test, HTTP, JSON3, SQLite, DBInterface, Tables; include("test/test_http.jl"); include("test/test_routes_series.jl")'
```
Expected: FAIL — `GET /api/series` returns 404 (route not registered) so `resp.status == 200` fails, or an `HTTP.StatusError` is thrown.

- [ ] **Step 3: Create `series.jl` with a minimal `series_listing`**

Create `packages/HimalayaUI/src/series.jl`:

```julia
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
```

- [ ] **Step 4: Create `routes_series.jl` with `register_series_routes!()` and `GET /api/series`**

Create `packages/HimalayaUI/src/routes_series.jl`:

```julia
using HTTP, JSON3, DBInterface, Tables, Oxygen, SQLite

# ─────────────────────────────────────────────────────────────────────────────
# I2.2 (#165) — REST routes for the series model. Mirrors
# `routes_comparisons.jl` (corpus-wide — there is no `/api/experiments/{eid}/series`).
#
# Every mutating route wraps in `with_idempotency(db, req)` and uses
# `apply_event!(InTransaction(), …)`. Body-shape validation runs BEFORE
# `with_idempotency` so a malformed payload returns an uncached 400.
#
# No route carries an `is_author` / 403 gate (architecture decision 3).
# Existence (404) and optimistic-concurrency (409) checks remain.
#
# `_json_error` and `_view_fields_error` are reused from `routes_comparisons.jl`
# (same module). Phase 3 (#175 / I3.6) relocates them when it deletes that file.
# ─────────────────────────────────────────────────────────────────────────────

function register_series_routes!()
    # ── Listing ─────────────────────────────────────────────────────────────

    @get "/api/series" function(req::HTTP.Request)
        db = current_db()
        rows = series_listing(db)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(rows))
    end
end
```

- [ ] **Step 5: Wire the includes and registration**

In `packages/HimalayaUI/src/HimalayaUI.jl`, add `include("series.jl")` immediately after `include("comparisons.jl")` (line 15):

```julia
include("comparisons.jl")
include("series.jl")
include("events.jl")
```

In the same file, add `include("routes_series.jl")` immediately after `include("routes_comparisons.jl")` (line 27):

```julia
include("routes_comparisons.jl")
include("routes_series.jl")
include("routes_picker.jl")
```

In `packages/HimalayaUI/src/server.jl`, add `register_series_routes!()` immediately after `register_comparisons_routes!()` (line 129):

```julia
    register_comparisons_routes!()
    register_series_routes!()
    register_picker_routes!()
```

In `packages/HimalayaUI/test/runtests.jl`, add `include("test_routes_series.jl")` immediately after `include("test_routes_comparisons.jl")` (line 42):

```julia
    include("test_routes_comparisons.jl")
    include("test_routes_series.jl")
    include("test_picker_routes.jl")
```

- [ ] **Step 6: Run the test to verify it passes**

Run (from `packages/HimalayaUI/`):
```bash
julia --project=. -e 'using HimalayaUI, Test, HTTP, JSON3, SQLite, DBInterface, Tables; include("test/test_http.jl"); include("test/test_routes_series.jl")'
```
Expected: PASS — `Series routes` testset, 2 tests (`resp.status == 200`, body is `[]`).

- [ ] **Step 7: Commit**

```bash
cd /home/jonathanchen/projects/Himalaya.jl/.claude/worktrees/series-routes-i2.2
git add packages/HimalayaUI/src/series.jl packages/HimalayaUI/src/routes_series.jl \
        packages/HimalayaUI/test/test_routes_series.jl \
        packages/HimalayaUI/src/HimalayaUI.jl packages/HimalayaUI/src/server.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "$(printf 'feat: scaffold series.jl + routes_series.jl with GET /api/series\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 2: `series_listing` + `_series_listing_rows` + the `last_event_at` sort fix (#76)

**Files:**
- Modify: `packages/HimalayaUI/src/series.jl`
- Modify: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Append inside the `@testset "Series routes"` block in `test/test_routes_series.jl`, before the closing `end`:

```julia
    @testset "series_listing — last_event_at sort is recency-correct (#76)" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            # Series A (id 1): updated_at is a T-separated ISO string with a Z
            # suffix; no events. Series B (id 2): an older updated_at, but a
            # `user_actions` event whose space-separated timestamp is MORE
            # recent than A's updated_at. True recency order: B before A.
            DBInterface.execute(db, """INSERT INTO series (id, title, updated_at, state)
                VALUES (1, 'A', '2026-05-01T00:00:00.000Z', 'committed')""")
            DBInterface.execute(db, """INSERT INTO series (id, title, updated_at, state)
                VALUES (2, 'B', '2026-04-01T00:00:00.000Z', 'committed')""")
            DBInterface.execute(db, """INSERT INTO user_actions
                (action, entity_type, entity_id, timestamp)
                VALUES ('series_recipe_updated', 'series', 2, '2026-05-10 12:00:00')""")

            listing = HimalayaUI.series_listing(db)
            @test length(listing) == 2
            # Bug #76: lexical sort would put A first ('T' > ' '). The
            # datetime() wrapper sorts by true recency — B first.
            @test listing[1][:id] == 2
            @test listing[2][:id] == 1
            # The projected last_event_at is normalised (no 'T'/'Z'), so it is
            # itself a valid client sort key — the bug is closed end-to-end.
            @test !occursin('T', listing[1][:last_event_at])
            @test !occursin('Z', listing[1][:last_event_at])
            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the per-task command. Expected: FAIL — `series_listing` still returns `[]`, so `length(listing) == 2` fails.

- [ ] **Step 3: Implement `series_listing` and `_series_listing_rows`**

In `series.jl`, replace the placeholder `series_listing` body with the real query, and add `_series_listing_rows` after it:

```julia
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
```

- [ ] **Step 4: Run the test to verify it passes**

Run the per-task command. Expected: PASS — the new testset's 5 tests pass; the Task 1 testset still passes.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/series.jl packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(printf 'feat: series_listing with the datetime() last_event_at sort fix (#76)\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 3: `fetch_series_with_plate` + `GET /api/series/{id}`

**Files:**
- Modify: `packages/HimalayaUI/src/series.jl`
- Modify: `packages/HimalayaUI/src/routes_series.jl`
- Modify: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Append inside the `@testset "Series routes"` block, before its closing `end`:

```julia
    @testset "GET /api/series/{id}" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                # Missing series → 404.
                resp404 = HTTP.get("$base/api/series/999", ["X-Username" => "alice"];
                                   status_exception = false)
                @test resp404.status == 404

                # Seed a series with one recipe row and one plate member.
                DBInterface.execute(db, """INSERT INTO series
                    (id, title, state) VALUES (5, 'S5', 'draft')""")
                DBInterface.execute(db, """INSERT INTO series_samples
                    (series_id, sample_id, position, pinned, excluded)
                    VALUES (5, 100, 0, 1, 0)""")
                DBInterface.execute(db, """INSERT INTO series_members
                    (series_id, exposure_id, display_order, snapshot, created_at)
                    VALUES (5, 1000, 0, '{"effective_peaks":[],"confirmed_index":null,"analysis_inputs_hash":null}', '2026-05-01T00:00:00.000Z')""")

                resp = HTTP.get("$base/api/series/5", ["X-Username" => "alice"])
                @test resp.status == 200
                got = JSON3.read(resp.body, Dict{Symbol, Any})
                @test got[:id] == 5
                @test got[:title] == "S5"
                @test got[:state] == "draft"
                @test length(got[:members]) == 1
                @test got[:members][1][:exposure_id] == 1000
                @test length(got[:samples]) == 1
                @test got[:samples][1][:sample_id] == 100
                @test got[:samples][1][:pinned] == true
            end
            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the per-task command. Expected: FAIL — `GET /api/series/{id}` is not registered (404 for the seeded series too, or `series not found` only by accident); `fetch_series_with_plate` is undefined.

- [ ] **Step 3: Implement `fetch_series_with_plate`**

Append to `series.jl`:

```julia
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
        :order_rule            => ismissing(s.order_rule) ? "manual" : String(s.order_rule),
        :state                 => ismissing(s.state) ? "committed" : String(s.state),
        :members               => members,
        :samples               => samples,
    )
end
```

- [ ] **Step 4: Add the `GET /api/series/{id}` route**

In `routes_series.jl`, inside `register_series_routes!()`, after the `GET "/api/series"` block:

```julia
    # ── Detail ──────────────────────────────────────────────────────────────

    @get "/api/series/{id}" function(req::HTTP.Request, id::Int)
        db = current_db()
        out = fetch_series_with_plate(db, id)
        out === nothing && return _json_error(404, "series not found")
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
    end
```

- [ ] **Step 5: Run the test to verify it passes**

Run the per-task command. Expected: PASS — the `GET /api/series/{id}` testset passes; earlier testsets still pass.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/series.jl packages/HimalayaUI/src/routes_series.jl \
        packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(printf 'feat: fetch_series_with_plate + GET /api/series/{id}\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 4: `forks_of_series` + `GET /api/series/{id}/forks`

**Files:**
- Modify: `packages/HimalayaUI/src/series.jl`
- Modify: `packages/HimalayaUI/src/routes_series.jl`
- Modify: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Append inside the `@testset "Series routes"` block:

```julia
    @testset "GET /api/series/{id}/forks" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                # No forks → empty array.
                resp0 = HTTP.get("$base/api/series/1/forks", ["X-Username" => "alice"])
                @test resp0.status == 200
                @test JSON3.read(resp0.body) == []

                # Parent (id 7) + a fork (id 8, forked_from_id = 7).
                DBInterface.execute(db, """INSERT INTO series (id, title, state)
                    VALUES (7, 'parent', 'committed')""")
                DBInterface.execute(db, """INSERT INTO series
                    (id, title, state, forked_from_id)
                    VALUES (8, 'child', 'committed', 7)""")

                resp = HTTP.get("$base/api/series/7/forks", ["X-Username" => "alice"])
                @test resp.status == 200
                forks = JSON3.read(resp.body, Vector{Dict{Symbol, Any}})
                @test length(forks) == 1
                @test forks[1][:id] == 8
                @test forks[1][:forked_from_id] == 7
            end
            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the per-task command. Expected: FAIL — `/api/series/{id}/forks` is not registered; `forks_of_series` is undefined.

- [ ] **Step 3: Implement `forks_of_series`**

Append to `series.jl`:

```julia
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
           WHERE s.forked_from_id = ?
           ORDER BY last_event_at DESC, s.id DESC""", [Int(series_id)]))
    _series_listing_rows(rows)
end
```

- [ ] **Step 4: Add the `GET /api/series/{id}/forks` route**

In `routes_series.jl`, after the `GET "/api/series/{id}"` block:

```julia
    @get "/api/series/{id}/forks" function(req::HTTP.Request, id::Int)
        db = current_db()
        rows = forks_of_series(db, id)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(rows))
    end
```

- [ ] **Step 5: Run the test to verify it passes**

Run the per-task command. Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/series.jl packages/HimalayaUI/src/routes_series.jl \
        packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(printf 'feat: forks_of_series + GET /api/series/{id}/forks\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 5: `compute_series_content_hash` + `current_series_content_hash` + `series_exists`

These three helpers feed the mutating routes (Tasks 6–9). They are pure functions tested directly.

**Files:**
- Modify: `packages/HimalayaUI/src/series.jl`
- Modify: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Append inside the `@testset "Series routes"` block:

```julia
    @testset "content-hash + existence helpers" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            # series_exists / current_series_content_hash on a missing id.
            @test HimalayaUI.series_exists(db, 404) == false
            @test HimalayaUI.current_series_content_hash(db, 404) === nothing

            # A draft series: it exists, but its content_hash is NULL — the
            # reason series_exists must be a separate probe (spec §4).
            DBInterface.execute(db, """INSERT INTO series (id, title, state)
                VALUES (3, 'draft-s', 'draft')""")
            @test HimalayaUI.series_exists(db, 3) == true
            @test HimalayaUI.current_series_content_hash(db, 3) === nothing

            # compute_series_content_hash is plate-only: adding a recipe row
            # (series_samples) must NOT change the hash; adding a plate member
            # (series_members) must.
            h_empty = HimalayaUI.compute_series_content_hash(db, 3)
            DBInterface.execute(db, """INSERT INTO series_samples
                (series_id, sample_id, position) VALUES (3, 100, 0)""")
            @test HimalayaUI.compute_series_content_hash(db, 3) == h_empty
            DBInterface.execute(db, """INSERT INTO series_members
                (series_id, exposure_id, display_order, snapshot, created_at)
                VALUES (3, 1000, 0, '{"effective_peaks":[],"confirmed_index":null,"analysis_inputs_hash":null}', '2026-05-01T00:00:00.000Z')""")
            @test HimalayaUI.compute_series_content_hash(db, 3) != h_empty

            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the per-task command. Expected: FAIL — `series_exists`, `current_series_content_hash`, `compute_series_content_hash` are undefined.

- [ ] **Step 3: Implement the three helpers**

Append to `series.jl`:

```julia
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
```

- [ ] **Step 4: Run the test to verify it passes**

Run the per-task command. Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/series.jl packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(printf 'feat: series content-hash + existence helpers\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 6: `POST /api/series` — create a draft

The mutating routes write a `user_actions` row but no view rows until #166–#168 (`update_view_for_event!` no-ops the unknown `series_created` kind). The test is therefore an HTTP-level smoke test: 4xx on a bad body, 201 on a good one, and a `user_actions` row written.

**Files:**
- Modify: `packages/HimalayaUI/src/series.jl`
- Modify: `packages/HimalayaUI/src/routes_series.jl`
- Modify: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Append inside the `@testset "Series routes"` block:

```julia
    @testset "POST /api/series — create draft (smoke)" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                # Missing title → 400 (uncached validation error).
                resp400 = HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:samples => []));
                    status_exception = false)
                @test resp400.status == 400

                # samples present but not an array → 400.
                resp400b = HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:title => "t", :samples => "nope"));
                    status_exception = false)
                @test resp400b.status == 400

                # Well-formed create → 201; a series row + a user_actions row.
                resp = HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(
                        :title => "My series",
                        :order_rule => "ascending",
                        :samples => [Dict(:sample_id => 100, :position => 0,
                                          :pinned => false, :excluded => false)])))
                @test resp.status == 201
                created = JSON3.read(resp.body, Dict{Symbol, Any})
                new_id = created[:id]
                @test new_id isa Integer
                # Degenerate until #166: the dispatcher no-ops, so the body is
                # the placeholder projection — empty members, empty samples.
                @test created[:members] == []
                @test created[:samples] == []
                # The durable event row IS written.
                ev = Tables.rowtable(DBInterface.execute(db,
                    "SELECT action FROM user_actions WHERE entity_type='series' AND entity_id=?",
                    [new_id]))
                @test length(ev) == 1
                @test ev[1].action == "series_created"
            end
            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the per-task command. Expected: FAIL — `POST /api/series` is not registered.

- [ ] **Step 3: Implement `_series_sample_payload`**

Append to `series.jl`:

```julia
"""
    _series_sample_payload(m_in) -> Dict{Symbol, Any}

Normalize one recipe entry from a request body into the `series_samples`
payload shape. `sample_id` is required (a recipe row with no target is
unrenderable); `position` / `pinned` / `excluded` default to `0` / `false` /
`false`.
"""
function _series_sample_payload(m_in)
    Dict{Symbol, Any}(
        :sample_id => Int(m_in[:sample_id]),
        :position  => haskey(m_in, :position) ? Int(m_in[:position]) : 0,
        :pinned    => haskey(m_in, :pinned)   ? Bool(m_in[:pinned])   : false,
        :excluded  => haskey(m_in, :excluded) ? Bool(m_in[:excluded]) : false,
    )
end
```

- [ ] **Step 4: Add the `POST /api/series` route**

In `routes_series.jl`, after the `GET "/api/series"` block (before the detail routes), add:

```julia
    # ── Create (draft) ──────────────────────────────────────────────────────

    @post "/api/series" function(req::HTTP.Request)
        db   = current_db()
        body = json(req)

        # Validate the body shape BEFORE with_idempotency so a malformed
        # payload returns an uncached 400 (validation toast), not a cached 500.
        if !haskey(body, :title)
            return _json_error(400, "missing required field: title")
        end
        if haskey(body, :samples) && body.samples !== nothing &&
                !(body.samples isa AbstractVector)
            return _json_error(400, "samples must be an array")
        end
        verr = _view_fields_error(body)
        verr === nothing || return verr

        title = String(body.title)
        description = haskey(body, :description) && body.description !== nothing ?
                      String(body.description) : nothing
        forked_from_id = haskey(body, :forked_from_id) && body.forked_from_id !== nothing ?
                         Int(body.forked_from_id) : nothing
        forked_at_hash = haskey(body, :forked_at_hash) && body.forked_at_hash !== nothing ?
                         String(body.forked_at_hash) : nothing
        ordering_variable = haskey(body, :ordering_variable) && body.ordering_variable !== nothing ?
                            String(body.ordering_variable) : nothing
        order_rule = haskey(body, :order_rule) && body.order_rule !== nothing ?
                     String(body.order_rule) : nothing
        view_grouping_mode = haskey(body, :view_grouping_mode) && body.view_grouping_mode !== nothing ?
            String(body.view_grouping_mode) : nothing
        view_show_peak_ticks  = haskey(body, :view_show_peak_ticks)  ? body.view_show_peak_ticks  : nothing
        view_show_peak_labels = haskey(body, :view_show_peak_labels) ? body.view_show_peak_labels : nothing
        samples_in = haskey(body, :samples) && body.samples !== nothing ? body.samples : ()

        return with_idempotency(db, req) do
            # Mint the AUTOINCREMENT id with a NULL-only placeholder row. The
            # `series_created` dispatcher (#166) upserts at this id and sets
            # `state='draft'`; until it lands the placeholder keeps the schema
            # default `state='committed'` — accepted (the row is degenerate).
            res = DBInterface.execute(db, "INSERT INTO series DEFAULT VALUES")
            new_id = Int(DBInterface.lastrowid(res))

            samples_payload = [_series_sample_payload(m) for m in samples_in]
            payload = Dict{Symbol, Any}(
                :title                 => title,
                :description           => description,
                :forked_from_id        => forked_from_id,
                :forked_at_hash        => forked_at_hash,
                :ordering_variable     => ordering_variable,
                :order_rule            => order_rule,
                :view_grouping_mode    => view_grouping_mode,
                :view_show_peak_ticks  => view_show_peak_ticks,
                :view_show_peak_labels => view_show_peak_labels,
                :samples               => samples_payload,
            )
            result = apply_event!(InTransaction(), db, req;
                kind        = "series_created",
                entity_type = "series",
                entity_id   = new_id,
                payload     = payload)

            out = fetch_series_with_plate(db, new_id)
            out === nothing && error(
                "post-write fetch_series_with_plate returned nothing for id=$(new_id)")
            # series_created carries NO post_state envelope (master-plan §5.2);
            # foreign tabs reconcile via replay-as-rerun (#166).
            _enqueue_broadcast_from_result!(result, "series_created",
                                            "series", new_id)
            HTTP.Response(201, ["Content-Type" => "application/json"], JSON3.write(out))
        end
    end
```

- [ ] **Step 5: Run the test to verify it passes**

Run the per-task command. Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/series.jl packages/HimalayaUI/src/routes_series.jl \
        packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(printf 'feat: POST /api/series — create a draft series\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 7: `PATCH /api/series/{id}` — recipe edit

**Files:**
- Modify: `packages/HimalayaUI/src/routes_series.jl`
- Modify: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Append inside the `@testset "Series routes"` block:

```julia
    @testset "PATCH /api/series/{id} — recipe edit (smoke)" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                # Missing series → 404.
                resp404 = HTTP.patch("$base/api/series/999",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:samples => []));
                    status_exception = false)
                @test resp404.status == 404

                # samples not an array → 400.
                DBInterface.execute(db, """INSERT INTO series (id, title, state)
                    VALUES (12, 'edit-me', 'draft')""")
                resp400 = HTTP.patch("$base/api/series/12",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:samples => "nope"));
                    status_exception = false)
                @test resp400.status == 400

                # Well-formed recipe edit → 200; a user_actions row written.
                resp = HTTP.patch("$base/api/series/12",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(
                        :order_rule => "descending",
                        :samples => [Dict(:sample_id => 100, :position => 0)])))
                @test resp.status == 200
                ev = Tables.rowtable(DBInterface.execute(db,
                    "SELECT action FROM user_actions WHERE entity_type='series' AND entity_id=12"))
                @test length(ev) == 1
                @test ev[1].action == "series_recipe_updated"
            end
            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the per-task command. Expected: FAIL — `PATCH /api/series/{id}` is not registered.

- [ ] **Step 3: Add the `PATCH /api/series/{id}` route**

In `routes_series.jl`, after the `GET "/api/series/{id}/forks"` block:

```julia
    # ── Recipe edit ─────────────────────────────────────────────────────────

    @patch "/api/series/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if haskey(body, :samples) && body.samples !== nothing &&
                !(body.samples isa AbstractVector)
            return _json_error(400, "samples must be an array")
        end
        ordering_variable = haskey(body, :ordering_variable) && body.ordering_variable !== nothing ?
                            String(body.ordering_variable) : nothing
        order_rule = haskey(body, :order_rule) && body.order_rule !== nothing ?
                     String(body.order_rule) : nothing
        samples_in = haskey(body, :samples) && body.samples !== nothing ? body.samples : ()

        return with_idempotency(db, req) do
            # No author gate (architecture decision 3). `series_exists` — not
            # `current_series_content_hash` — is the existence probe: a draft
            # has a NULL content_hash by design.
            if !series_exists(db, id)
                return _json_error(404, "series not found")
            end
            samples_payload = [_series_sample_payload(m) for m in samples_in]
            payload = Dict{Symbol, Any}(
                :ordering_variable => ordering_variable,
                :order_rule        => order_rule,
                :samples           => samples_payload,
            )
            result = apply_event!(InTransaction(), db, req;
                kind        = "series_recipe_updated",
                entity_type = "series",
                entity_id   = id,
                payload     = payload)

            out = fetch_series_with_plate(db, id)
            out === nothing && error(
                "post-write fetch_series_with_plate returned nothing for id=$(id)")
            # No post_state envelope (master-plan §5.2).
            _enqueue_broadcast_from_result!(result, "series_recipe_updated",
                                            "series", id)
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
        end
    end
```

- [ ] **Step 4: Run the test to verify it passes**

Run the per-task command. Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_series.jl packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(printf 'feat: PATCH /api/series/{id} — recipe edit\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 8: `POST /api/series/{id}/commit` — commit the plate

Carries a `post_state` envelope (the one series event that does — master-plan §5.2). Keeps the `expected_content_hash`→409 optimistic-concurrency check. No author gate.

**Files:**
- Modify: `packages/HimalayaUI/src/series.jl`
- Modify: `packages/HimalayaUI/src/routes_series.jl`
- Modify: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Append inside the `@testset "Series routes"` block:

```julia
    @testset "POST /api/series/{id}/commit (smoke + no gate + 409)" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                _member = Dict(:exposure_id => 1000, :display_order => 0,
                               :snapshot => Dict(:effective_peaks => Any[],
                                                 :confirmed_index => nothing,
                                                 :analysis_inputs_hash => nothing))

                # Missing series → 404.
                resp404 = HTTP.post("$base/api/series/999/commit",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:members => [_member]));
                    status_exception = false)
                @test resp404.status == 404

                # Seed a committed series authored by alice (id 1), with a
                # stored content_hash so the 409 path can be exercised.
                DBInterface.execute(db, "INSERT INTO users (id, username) VALUES (1, 'alice')")
                DBInterface.execute(db, """INSERT INTO series
                    (id, title, state, created_by, content_hash)
                    VALUES (20, 'committed-s', 'committed', 1, 'sha256:deadbeef')""")

                # No is_author gate: bob (not the author) commits → NOT 403.
                resp = HTTP.post("$base/api/series/20/commit",
                    ["X-Username" => "bob", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:members => [_member]));
                    status_exception = false)
                @test resp.status != 403
                @test resp.status == 200
                ev = Tables.rowtable(DBInterface.execute(db,
                    "SELECT action FROM user_actions WHERE entity_type='series' AND entity_id=20"))
                @test any(r -> r.action == "series_plate_committed", ev)

                # Conflict: a wrong expected_content_hash → 409.
                resp409 = HTTP.post("$base/api/series/20/commit",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:members => [_member],
                                     :expected_content_hash => "sha256:WRONG"));
                    status_exception = false)
                @test resp409.status == 409
                conflict = JSON3.read(resp409.body, Dict{Symbol, Any})
                @test conflict[:error] == "conflict"
            end
            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the per-task command. Expected: FAIL — `POST /api/series/{id}/commit` is not registered.

- [ ] **Step 3: Implement `_series_member_payload`**

Append to `series.jl`:

```julia
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
```

- [ ] **Step 4: Add the `POST /api/series/{id}/commit` route**

In `routes_series.jl`, after the `PATCH "/api/series/{id}"` block:

```julia
    # ── Commit the plate (the old "submit") ─────────────────────────────────

    @post "/api/series/{id}/commit" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if !haskey(body, :members) || !(body.members isa AbstractVector)
            return _json_error(400, "members must be an array")
        end
        verr = _view_fields_error(body)
        verr === nothing || return verr
        expected_hash = haskey(body, :expected_content_hash) &&
                        body.expected_content_hash !== nothing ?
                        String(body.expected_content_hash) : nothing

        return with_idempotency(db, req) do
            # Existence (404) before the conflict check (409) — HTTP semantics.
            # No author gate (architecture decision 3).
            if !series_exists(db, id)
                return _json_error(404, "series not found")
            end
            # Optimistic-concurrency check (NOT the author gate): the stored
            # hash must match the client's expected_content_hash, else 409.
            current_hash = current_series_content_hash(db, id)
            if expected_hash !== nothing && current_hash !== expected_hash
                current_state = fetch_series_with_plate(db, id)
                return HTTP.Response(409, ["Content-Type" => "application/json"],
                    JSON3.write(Dict(
                        :error         => "conflict",
                        :current_hash  => current_hash,
                        :current_state => current_state,
                    )))
            end

            members_payload = [_series_member_payload(db, m) for m in body.members]
            payload = Dict{Symbol, Any}(:members => members_payload)
            result = apply_event!(InTransaction(), db, req;
                kind        = "series_plate_committed",
                entity_type = "series",
                entity_id   = id,
                payload     = payload)

            out = fetch_series_with_plate(db, id)
            out === nothing && error(
                "post-write fetch_series_with_plate returned nothing for id=$(id)")
            # series_plate_committed is the one series event carrying a
            # post_state envelope (master-plan §5.2).
            _enqueue_broadcast_from_result!(result, "series_plate_committed",
                                            "series", id; post_state = out)
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
        end
    end
```

- [ ] **Step 5: Run the test to verify it passes**

Run the per-task command. Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/series.jl packages/HimalayaUI/src/routes_series.jl \
        packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(printf 'feat: POST /api/series/{id}/commit — commit the plate\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 9: `DELETE /api/series/{id}`

**Files:**
- Modify: `packages/HimalayaUI/src/routes_series.jl`
- Modify: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Append inside the `@testset "Series routes"` block:

```julia
    @testset "DELETE /api/series/{id} (smoke + no gate)" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                # Missing series → 404.
                resp404 = HTTP.delete("$base/api/series/999", ["X-Username" => "alice"];
                                      status_exception = false)
                @test resp404.status == 404

                # Series authored by alice; bob (not the author) deletes it —
                # no is_author gate, so NOT 403.
                DBInterface.execute(db, "INSERT INTO users (id, username) VALUES (1, 'alice')")
                DBInterface.execute(db, """INSERT INTO series
                    (id, title, state, created_by) VALUES (30, 'doomed', 'committed', 1)""")
                resp = HTTP.delete("$base/api/series/30", ["X-Username" => "bob"];
                                   status_exception = false)
                @test resp.status != 403
                @test resp.status == 200
                deleted = JSON3.read(resp.body, Dict{Symbol, Any})
                @test deleted[:id] == 30
                @test deleted[:deleted] == true
                ev = Tables.rowtable(DBInterface.execute(db,
                    "SELECT action FROM user_actions WHERE entity_type='series' AND entity_id=30"))
                @test any(r -> r.action == "series_deleted", ev)
            end
            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the per-task command. Expected: FAIL — `DELETE /api/series/{id}` is not registered.

- [ ] **Step 3: Add the `DELETE /api/series/{id}` route**

In `routes_series.jl`, after the `POST "/api/series/{id}/commit"` block:

```julia
    # ── Delete ──────────────────────────────────────────────────────────────

    @delete "/api/series/{id}" function(req::HTTP.Request, id::Int)
        db = current_db()
        return with_idempotency(db, req) do
            # No author gate (architecture decision 3).
            if !series_exists(db, id)
                return _json_error(404, "series not found")
            end
            result = apply_event!(InTransaction(), db, req;
                kind        = "series_deleted",
                entity_type = "series",
                entity_id   = id,
                payload     = Dict{Symbol, Any}(:id => id))
            # No post_state envelope (master-plan §5.2).
            _enqueue_broadcast_from_result!(result, "series_deleted", "series", id)
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:id => id, :deleted => true,
                                 :event_id => result.event_id)))
        end
    end
```

- [ ] **Step 4: Run the test to verify it passes**

Run the per-task command. Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_series.jl packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(printf 'feat: DELETE /api/series/{id}\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 10: Messages — `GET` / `POST /api/series/{id}/messages`

Fully functional in I2.2: `post_message` is an existing dispatcher kind, and the route writes the `series_messages` row directly.

**Files:**
- Modify: `packages/HimalayaUI/src/routes_series.jl`
- Modify: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Append inside the `@testset "Series routes"` block:

```julia
    @testset "series messages — full round-trip" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                DBInterface.execute(db, """INSERT INTO series (id, title, state)
                    VALUES (40, 'chatty', 'committed')""")

                # Empty thread.
                resp0 = HTTP.get("$base/api/series/40/messages", ["X-Username" => "alice"])
                @test resp0.status == 200
                @test JSON3.read(resp0.body) == []

                # POST without X-Username → 401.
                resp401 = HTTP.post("$base/api/series/40/messages",
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:body => "hi"));
                    status_exception = false)
                @test resp401.status == 401

                # POST with an empty body → 400.
                resp400 = HTTP.post("$base/api/series/40/messages",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:body => "   "));
                    status_exception = false)
                @test resp400.status == 400

                # POST a real message → 201, then GET sees it.
                resp = HTTP.post("$base/api/series/40/messages",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:body => "first post")))
                @test resp.status == 201
                msg = JSON3.read(resp.body, Dict{Symbol, Any})
                @test msg[:body] == "first post"
                @test msg[:series_id] == 40

                resp2 = HTTP.get("$base/api/series/40/messages", ["X-Username" => "alice"])
                thread = JSON3.read(resp2.body, Vector{Dict{Symbol, Any}})
                @test length(thread) == 1
                @test thread[1][:body] == "first post"
                @test thread[1][:author] == "alice"
            end
            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the per-task command. Expected: FAIL — the message routes are not registered.

- [ ] **Step 3: Add the message routes**

In `routes_series.jl`, after the `DELETE "/api/series/{id}"` block:

```julia
    # ── Chat thread ─────────────────────────────────────────────────────────

    @get "/api/series/{id}/messages" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db, """
            SELECT m.id, m.series_id, m.author_id,
                   u.username AS author,
                   m.body, m.created_at
            FROM series_messages m
            LEFT JOIN users u ON u.id = m.author_id
            WHERE m.series_id = ?
            ORDER BY m.created_at ASC, m.id ASC
        """, [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(rows_to_json(rows)))
    end

    @post "/api/series/{id}/messages" function(req::HTTP.Request, id::Int)
        db       = current_db()
        username = get_username(req)
        if username === nothing
            return _json_error(401, "X-Username header required")
        end
        body = json(req)
        text = haskey(body, :body) ? strip(String(body.body)) : ""
        if isempty(text)
            return _json_error(400, "message body required")
        end

        return with_idempotency(db, req) do
            # The route writes the message row directly; the dispatcher is a
            # no-op for `post_message`. `entity_type='series_message'`
            # differentiates it from the comparison / sample message paths.
            author_id = get_or_create_user!(db, username)
            res = DBInterface.execute(db,
                "INSERT INTO series_messages (series_id, author_id, body) VALUES (?, ?, ?)",
                [id, author_id, text])
            msg_id = Int(DBInterface.lastrowid(res))

            row = Tables.rowtable(DBInterface.execute(db, """
                SELECT m.id, m.series_id, m.author_id,
                       u.username AS author,
                       m.body, m.created_at
                FROM series_messages m
                LEFT JOIN users u ON u.id = m.author_id
                WHERE m.id = ?
            """, [msg_id]))[1]
            msg_json = row_to_json(row)

            result = apply_event!(InTransaction(), db, req;
                kind        = "post_message",
                entity_type = "series_message",
                entity_id   = msg_id,
                payload     = msg_json)
            _enqueue_broadcast_from_result!(result, "post_message",
                                            "series_message", msg_id)
            HTTP.Response(201, ["Content-Type" => "application/json"],
                JSON3.write(msg_json))
        end
    end
```

- [ ] **Step 4: Run the test to verify it passes**

Run the per-task command. Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_series.jl packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(printf 'feat: GET/POST /api/series/{id}/messages\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 11: Pins — `GET /api/users/me/series-pins`, `POST` / `DELETE /api/series/{id}/pin`

`GET /api/users/me/series-pins` is a plain read (works today); the pin/unpin routes emit `series_pinned` / `series_unpinned` (`entity_type='user'`), degenerate until #168.

**Files:**
- Modify: `packages/HimalayaUI/src/routes_series.jl`
- Modify: `packages/HimalayaUI/test/test_routes_series.jl`

- [ ] **Step 1: Write the failing test**

Append inside the `@testset "Series routes"` block:

```julia
    @testset "series pins (smoke)" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                # series-pins without X-Username → 401.
                resp401 = HTTP.get("$base/api/users/me/series-pins";
                                   status_exception = false)
                @test resp401.status == 401

                # Fresh user → empty pin list.
                resp0 = HTTP.get("$base/api/users/me/series-pins", ["X-Username" => "alice"])
                @test resp0.status == 200
                @test JSON3.read(resp0.body) == []

                # Pin a missing series → 404.
                resp404 = HTTP.post("$base/api/series/999/pin", ["X-Username" => "alice"];
                                    status_exception = false)
                @test resp404.status == 404

                # Pin an existing series → 200, a user_actions row written.
                DBInterface.execute(db, """INSERT INTO series (id, title, state)
                    VALUES (50, 'pin-me', 'committed')""")
                respPin = HTTP.post("$base/api/series/50/pin", ["X-Username" => "alice"])
                @test respPin.status == 200
                pinned = JSON3.read(respPin.body, Dict{Symbol, Any})
                @test pinned[:series_id] == 50
                @test pinned[:pinned] == true
                evP = Tables.rowtable(DBInterface.execute(db,
                    "SELECT action FROM user_actions WHERE entity_type='user'"))
                @test any(r -> r.action == "series_pinned", evP)

                # Unpin → 200.
                respUnpin = HTTP.delete("$base/api/series/50/pin", ["X-Username" => "alice"])
                @test respUnpin.status == 200
                unpinned = JSON3.read(respUnpin.body, Dict{Symbol, Any})
                @test unpinned[:pinned] == false
            end
            close(db)
        end
    end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the per-task command. Expected: FAIL — the pin routes are not registered.

- [ ] **Step 3: Add the pin routes**

In `routes_series.jl`, after the message routes — inserted before the single
closing `end` of `register_series_routes!()`, exactly like every prior task:

```julia
    # ── Pins ────────────────────────────────────────────────────────────────
    #
    # Per-user pinned series. Pin/unpin events are stored with
    # `entity_type='user'`, `entity_id=user_id` (the `comparison_pinned`
    # precedent) — the affected series rides in the payload as `series_id`.
    # The dispatcher branches land in #168; until then these routes write a
    # `user_actions` row but no `series_pins` row.

    @get "/api/users/me/series-pins" function(req::HTTP.Request)
        db = current_db()
        user_id = get_user_id_for_request(db, req)
        if user_id === nothing
            return _json_error(401, "X-Username header required")
        end
        rows = Tables.rowtable(DBInterface.execute(db, """
            SELECT series_id FROM series_pins
            WHERE user_id = ?
            ORDER BY pinned_at DESC, series_id DESC
        """, [user_id]))
        ids = Int[Int(r.series_id) for r in rows]
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(ids))
    end

    @post "/api/series/{id}/pin" function(req::HTTP.Request, id::Int)
        db = current_db()
        user_id = get_user_id_for_request(db, req)
        if user_id === nothing
            return _json_error(401, "X-Username header required")
        end
        return with_idempotency(db, req) do
            if !series_exists(db, id)
                return _json_error(404, "series not found")
            end
            result = apply_event!(InTransaction(), db, req;
                kind        = "series_pinned",
                entity_type = "user",
                entity_id   = user_id,
                payload     = Dict{Symbol, Any}(:series_id => id))
            _enqueue_broadcast_from_result!(result, "series_pinned", "user", user_id)
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:series_id => id, :pinned => true)))
        end
    end

    @delete "/api/series/{id}/pin" function(req::HTTP.Request, id::Int)
        db = current_db()
        user_id = get_user_id_for_request(db, req)
        if user_id === nothing
            return _json_error(401, "X-Username header required")
        end
        return with_idempotency(db, req) do
            # Idempotent at the SQL layer once #168 lands — unpinning a
            # never-pinned series is a no-op DELETE. The event is still
            # recorded so the cross-tab broadcast fires either way.
            result = apply_event!(InTransaction(), db, req;
                kind        = "series_unpinned",
                entity_type = "user",
                entity_id   = user_id,
                payload     = Dict{Symbol, Any}(:series_id => id))
            _enqueue_broadcast_from_result!(result, "series_unpinned", "user", user_id)
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:series_id => id, :pinned => false)))
        end
    end
```

Note: the three `@get`/`@post`/`@delete` blocks above are inserted *inside*
`register_series_routes!()`, before its existing closing `end`. Do not add a
new `end` — the function keeps the single closing `end` it has had since Task 1.

- [ ] **Step 4: Run the test to verify it passes**

Run the per-task command. Expected: PASS — all 11 series-route testsets pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_series.jl packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(printf 'feat: series pin routes (GET series-pins, POST/DELETE pin)\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 12: Full-suite verification

No code change — confirm the new file integrates and nothing regressed.

- [ ] **Step 1: Run the full backend suite**

From `packages/HimalayaUI/`:
```bash
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-series-test.out 2>&1
```
This is slow (5–10 min). When it finishes, inspect the result:
```bash
grep -E 'Test Summary|FAIL|Error|Series routes' /tmp/jl-series-test.out
```
Expected: the `Series routes` testset reports all tests passing; the overall summary shows no failures. (Per `memory/himalaya_suite_flakes.md`, the fast-skip P99 latency test and the `GET /api/health` port-race test flake intermittently and unrelatedly — if only those two fail, that is not a regression. Any `series`-related failure IS a regression.)

- [ ] **Step 2: Confirm the frontend is unaffected**

No frontend file was touched. Confirm:
```bash
cd /home/jonathanchen/projects/Himalaya.jl/.claude/worktrees/series-routes-i2.2
git diff --name-only 41070a4..HEAD -- packages/HimalayaUI/frontend
```
Expected: no output (zero frontend files changed).

- [ ] **Step 3: Confirm the acceptance criteria**

Re-read `docs/superpowers/specs/2026-05-18-series-routes-design.md` §7 and confirm each box against the suite output:
- All series routes round-trip on empty `series` tables — GET routes + messages fully; mutating routes at the HTTP level.
- The `last_event_at` sort test passes (Task 2).
- No `is_author` gate — the non-author commit (Task 8) and delete (Task 9) tests pass without 403.
- Julia route tests pass.
- The frontend is unaffected (Step 2).

- [ ] **Step 4: Final commit (only if Step 1–3 surfaced a fix)**

If Steps 1–3 were all clean, there is nothing to commit — Task 11 already committed the last code change. If a fix was needed, commit it:
```bash
git add -A
git commit -m "$(printf 'fix: <describe the regression fix>\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Self-review notes (for the implementer)

- **Spec coverage:** every spec §5 route maps to a task (5.1→T1, 5.2→T6, 5.3→T3, 5.4→T4, 5.5→T7, 5.6→T8, 5.7→T9, 5.8→T10, 5.9→T11); every spec §4 helper maps to a task (`series_listing`/`_series_listing_rows`→T2, `forks_of_series`→T4, `fetch_series_with_plate`→T3, `compute_series_content_hash`/`current_series_content_hash`/`series_exists`→T5, `_series_member_payload`→T8, `_series_sample_payload`→T6); spec §6 tests are distributed across the tasks and confirmed in T12.
- **`apply_event!` signature:** all mutating routes call `apply_event!(InTransaction(), db, req; kind, entity_type, entity_id, payload)` and read `result.event_id` — identical to `routes_comparisons.jl`.
- **No new module-level helpers were invented:** `current_db`, `json`, `with_idempotency`, `apply_event!`, `InTransaction`, `_enqueue_broadcast_from_result!`, `_json_error`, `_view_fields_error`, `get_username`, `get_user_id_for_request`, `get_or_create_user!`, `rows_to_json`, `row_to_json`, `canonical_json`, `compute_member_snapshot`, `is_member_stale`, `_topk_phases`, `_count_distinct_phases` all already exist in the `HimalayaUI` module.
