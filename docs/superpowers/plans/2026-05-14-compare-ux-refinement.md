# Compare UX refinement — Phase 1 (heavy-lifting) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship the structural refinements to the Compare page from
`docs/superpowers/specs/2026-05-14-compare-ux-refinement-design.md` — a
single-mode page with an inline-edit gesture surface, a calmer header,
revamped member rows, and a sidebar that actually tells you what's in each
comparison.

**Spec:** [`docs/superpowers/specs/2026-05-14-compare-ux-refinement-design.md`](../specs/2026-05-14-compare-ux-refinement-design.md).

**Deliberate scope split (read before reviewing).** This plan covers the
*heavy-lifting* subset of the spec. The remaining items are deferred to a
follow-up Phase 2 plan, by mutual agreement between the spec author and
the implementer (recorded here so spec/plan reviewers don't flag them as
omissions):

| Status | Scope |
|---|---|
| **In this plan (Phase 1)** | Phase A (backend projection + migration), Phase B (route flattening), Phase C (title strip + status surface + toolbar + Save pill + merge into `Compare.tsx`), Phase E (member row redesign), Phase F (sidebar redesign). §G visual-language vocabulary is applied **inline** as each component is touched (not a separate phase). §H deletes happen naturally inside B/C. |
| **Deferred to Phase 2** | Phase D (right-rail Chat ↔ + Add traces tabs — interim picker chip ships in Phase 1). The external add-trace drag flow from spec §7.3 (reorder ships; cross-surface exposure drops wait). Color picker swatch row in the expanded member row. Cmd+Z undo stack. Description inline-edit (today's input box still works). Sidebar `⋯` overflow per spec §4/§8.2 (pin star stays in Phase 1). "Reset to author's view" right-click affordance per spec §6.4. Inline-modal replacement for `window.prompt` in C-14 fork-title gesture. Actor attribution on the `viewing-stale` banner (spec §6.2 updated to drop `@user` in Phase 1). Queue-interleave coverage that drives `replayCoordinator.handleRemoteEvent` directly for the comparison whole-row overwrite (A-9 Step 3d ships sequential-write coverage only). |
| **Resolved at impl time** | The two spec §14 open questions (picker scope on `/compare/all`, orphan-member handling for `has_stale_members`). |

**Architecture.** Five phases land in order; A→B→C→E→F. Each phase is
independently shippable with green tests. Phase A is additive only
(backend migration + projection); Phase B drops the `/edit` URL; Phase C
folds `ComparePage` + `ComparePageEdit` into a single `Compare.tsx` with
the inline-edit gesture surface; Phase E reworks the member gutter rows;
Phase F redesigns sidebar rows against Phase A's projection.

**Tech stack.** Backend: Julia 1.10, SQLite.jl, Oxygen.jl 1.10.x. Frontend:
React 18, TypeScript 5 (strict, `exactOptionalPropertyTypes: true`),
TanStack Query, Zustand, Tailwind CSS v4, Observable Plot, Vitest, React
Testing Library, Playwright.

---

## File structure inventory

These files are touched by this plan. Each phase task lists the
sub-changes it owns.

**Create**
- `packages/HimalayaUI/frontend/src/pages/Compare.tsx` — single-mode page, replaces `ComparePage.tsx` + `ComparePageEdit.tsx`.
- `packages/HimalayaUI/frontend/src/hooks/useCompareMode.ts` — tagged-union `viewing | viewing-stale | editing-mine | editing-as-fork-of | creating-blank | creating-from-fork`.
- `packages/HimalayaUI/frontend/src/hooks/useCompareDraftDirty.ts` — boolean signal sourced from active draft.
- `packages/HimalayaUI/frontend/src/components/InlineEditableText.tsx` — text → input on click; Enter commits, Esc cancels.
- `packages/HimalayaUI/frontend/src/components/CompareTitleStrip.tsx` — title + meta line (author / age / forked-from).
- `packages/HimalayaUI/frontend/src/components/CompareStatusSurface.tsx` — banner stack (NeedsReview, ServerUpdated, Saved).
- `packages/HimalayaUI/frontend/src/components/SavePill.tsx` — context-aware Save changes / Save as fork pill.
- `packages/HimalayaUI/frontend/src/components/CompareToolbar.tsx` — Group ▾ / Annotations ▾ / ⋯ more / Export + SavePill mount.
- `packages/HimalayaUI/frontend/src/components/RowActionZone.tsx` — right-edge action zone primitive (`⋯` + `⋮⋮` cue).
- `packages/HimalayaUI/frontend/src/lib/comparison/dragThreshold.ts` — `pointerdown → 4px → drag` guard helper.
- `packages/HimalayaUI/frontend/src/lib/comparison/relativeTime.ts` — `"edited 2h ago"` formatter helper.
- Test files paired to each new component.

**Modify**
- `packages/HimalayaUI/src/db.jl` — `migrate_compare_view_choices!`.
- `packages/HimalayaUI/src/comparisons.jl` — extend `_comparison_listing_rows` projector; add `view_*` columns to `fetch_comparison_with_members`.
- `packages/HimalayaUI/src/routes_comparisons.jl` — accept + echo `view_*` fields on submit.
- `packages/HimalayaUI/src/events.jl` — `_update_view_for_comparison_submitted!` writes `view_*` fields.
- `packages/HimalayaUI/frontend/src/api.ts` — extend `ComparisonSummary`, `Comparison`, `SaveComparisonBody`, `ComparisonMemberInput`.
- `packages/HimalayaUI/frontend/src/state.ts` — view-choice draft setters; nothing removed.
- `packages/HimalayaUI/frontend/src/lib/comparison/draft.ts` — add `viewGroupingMode` / `viewShowPeakTicks` / `viewShowPeakLabels` to `ActiveDraft`.
- `packages/HimalayaUI/frontend/src/lib/comparison/draftFactories.ts` — wire view fields into draft factory functions.
- `packages/HimalayaUI/frontend/src/lib/queue/mutators/saveComparison.ts` — include view fields in the request body.
- `packages/HimalayaUI/frontend/src/lib/comparison/routes.ts` — drop `edit?: boolean` option; update `comparePath` callers.
- `packages/HimalayaUI/frontend/src/components/AppShell.tsx` — route `/compare/...*` to `Compare.tsx`; redirect `*/edit` to bare.
- `packages/HimalayaUI/frontend/src/components/ComparisonSidebar.tsx` — consume new projection, replace `"3 traces"` placeholder, new row layout.
- `packages/HimalayaUI/frontend/src/components/MemberMetaRow.tsx` — collapse/expand, right-edge actions, click-vs-drag threshold.
- `packages/HimalayaUI/frontend/src/components/MemberMetaGutter.tsx` — visible inter-row gap, drop indicator, plot mirror, threshold-guarded reorder.
- `packages/HimalayaUI/frontend/src/components/BandResizeDivider.tsx` — replaced by gap-based resize; component file simplifies or is folded into `MemberMetaGutter.tsx`.

**Delete (at end of Phase C)**
- `packages/HimalayaUI/frontend/src/pages/ComparePage.tsx` — folded into `Compare.tsx`.
- `packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx` — folded into `Compare.tsx`.
- `packages/HimalayaUI/frontend/src/components/LineageBadge.tsx` — folded into `CompareTitleStrip.tsx`.
- `packages/HimalayaUI/frontend/src/components/NeedsReviewBadge.tsx` — folded into `CompareStatusSurface.tsx`.
- `packages/HimalayaUI/frontend/src/components/ForksPopover.tsx` — content moved into the toolbar's ⋯-more dropdown.

---

## Pre-flight

### Task PF-1: Worktree sanity check

**Files:** (none — environmental)

- [ ] **Step 1: Confirm worktree.**

```bash
pwd
git rev-parse --abbrev-ref HEAD
git status -sb
```

Expected: working directory ends in `.claude/worktrees/compare-ux-refinement-spec`; branch starts with `worktree-`; tree clean apart from this plan.

- [ ] **Step 2: First-time setup (only if instantiate hasn't been run).**

```bash
(cd packages/HimalayaUI/frontend && npm install)
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.instantiate()'
```

- [ ] **Step 3: Establish green baseline.**

```bash
cd packages/HimalayaUI/frontend && npm run build
```

Expected: TypeScript compiles, Vite build succeeds.

```bash
cd packages/HimalayaUI/frontend && npm test -- --run > /tmp/vitest-baseline.out 2>&1
grep -E "Tests:|Failed|passed" /tmp/vitest-baseline.out | tail -5
```

Expected: all Vitest tests pass.

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-baseline.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-baseline.out | tail -10
```

Expected: all Julia tests pass (5–10 min).

If anything is red on baseline, STOP and surface to the user before touching code.

---

## Phase A — Backend listing projection + view-choice migration

Phase goal: extend the comparisons listing payload with `author_username`,
`member_count`, `member_phases`, `has_stale_members`, `last_event_at`; add
the three view-choice columns to `comparisons` and round-trip them through
GET / POST / SSE. Frontend types updated to match.

### Task A-1: Add `migrate_compare_view_choices!` to `db.jl`

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl`
- Test: `packages/HimalayaUI/test/test_db.jl`

- [ ] **Step 1: Write the failing migration test.**

Append to `packages/HimalayaUI/test/test_db.jl`:

```julia
@testset "migrate_compare_view_choices!" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "himalaya.db"))
        # All three view-choice columns must exist on a freshly opened DB.
        cols = Set(r.name for r in
            Tables.rowtable(DBInterface.execute(db,
                "PRAGMA table_info(comparisons)")))
        @test "view_grouping_mode" in cols
        @test "view_show_peak_ticks" in cols
        @test "view_show_peak_labels" in cols

        # Each is nullable (default NULL).
        DBInterface.execute(db,
            "INSERT INTO comparisons (id, title) VALUES (1, 'x')")
        row = Tables.rowtable(DBInterface.execute(db,
            "SELECT view_grouping_mode, view_show_peak_ticks, view_show_peak_labels
             FROM comparisons WHERE id = 1"))[1]
        @test ismissing(row.view_grouping_mode)
        @test ismissing(row.view_show_peak_ticks)
        @test ismissing(row.view_show_peak_labels)
        close(db)
    end
end
```

- [ ] **Step 2: Run the test to verify it fails.**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI, Test; include("packages/HimalayaUI/test/test_db.jl")' 2>&1 | grep -E "view_grouping_mode|view_show_peak|migrate_compare_view"
```

Expected: failure — columns don't exist yet.

- [ ] **Step 3: Add the migration function in `db.jl`.**

Locate `function migrate_compare_relax_nullability!` (the migration immediately preceding ours in run order). After its definition closes, add:

```julia
"""
    migrate_compare_view_choices!(db)

Add `view_grouping_mode`, `view_show_peak_ticks`, `view_show_peak_labels`
columns to `comparisons` so the author's view choices round-trip across
viewers (spec §6.4). All NULL on existing rows; the frontend falls back to
per-tab Zustand defaults when NULL.

Idempotent: each `ALTER TABLE ... ADD COLUMN` is wrapped in a try/catch
that treats "duplicate column name" as success — the same pattern used
by other additive migrations.
"""
function migrate_compare_view_choices!(db::SQLite.DB)
    for stmt in (
        "ALTER TABLE comparisons ADD COLUMN view_grouping_mode TEXT",
        "ALTER TABLE comparisons ADD COLUMN view_show_peak_ticks INTEGER",
        "ALTER TABLE comparisons ADD COLUMN view_show_peak_labels INTEGER",
    )
        try
            DBInterface.execute(db, stmt)
        catch err
            msg = sprint(showerror, err)
            occursin("duplicate column name", lowercase(msg)) || rethrow()
        end
    end
end
```

- [ ] **Step 4: Wire it into `create_schema!`.**

Locate the call site of `migrate_compare_relax_nullability!(db)` in
`create_schema!`. Add immediately after it:

```julia
    # 2026-05-14 — Compare UX refinement (spec §6.4): persist view choices
    # on the comparison so they round-trip across viewers.
    migrate_compare_view_choices!(db)
```

- [ ] **Step 5: Run the test to verify it passes.**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI, Test; include("packages/HimalayaUI/test/test_db.jl")' 2>&1 | grep -E "Test Summary|did not pass|fail" | head -5
```

Expected: PASS.

- [ ] **Step 6: Commit.**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_db.jl
git commit -m "migrate: add comparisons.view_* columns (Compare UX A-1)"
```

### Task A-2: Update `migrate_compare_view_choices!` idempotency on existing DB

**Files:**
- Test: `packages/HimalayaUI/test/test_db.jl`

- [ ] **Step 1: Write the legacy-DB test.**

Append to `test_db.jl`:

```julia
@testset "migrate_compare_view_choices! is idempotent" begin
    mktempdir() do tmp
        # Open once → migration runs.
        db = open_db(joinpath(tmp, "himalaya.db"))
        close(db)
        # Open again → migration runs again; must not throw.
        db = open_db(joinpath(tmp, "himalaya.db"))
        cols = Set(r.name for r in
            Tables.rowtable(DBInterface.execute(db,
                "PRAGMA table_info(comparisons)")))
        @test "view_grouping_mode" in cols
        close(db)
    end
end
```

- [ ] **Step 2: Run.**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI, Test; include("packages/HimalayaUI/test/test_db.jl")' 2>&1 | tail -5
```

Expected: PASS (the try/catch on duplicate column name is what makes this pass).

- [ ] **Step 3: Commit.**

```bash
git add packages/HimalayaUI/test/test_db.jl
git commit -m "test: migrate_compare_view_choices! idempotency (Compare UX A-2)"
```

### Task A-3: Extend `_comparison_listing_rows` projection

**Files:**
- Modify: `packages/HimalayaUI/src/comparisons.jl`
- Test: `packages/HimalayaUI/test/test_comparisons.jl`

Goal: each Dict returned by `_comparison_listing_rows` gains five new keys:
`author_username`, `member_count`, `member_phases`, `has_stale_members`,
`last_event_at`. The `last_event_at` value is already computed by the
listing queries' COALESCE subquery — `_comparison_listing_rows` just needs
to project it. The other four are new and require modifying the
`comparisons_for_experiment` and `comparisons_listing` SELECT statements.

- [ ] **Step 1: Write the failing test.**

Append to `packages/HimalayaUI/test/test_comparisons.jl`:

```julia
@testset "listing projection — Compare UX A-3" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "himalaya.db"))

        # Minimal fixture: user, experiment, sample, exposure, peak,
        # comparison with one member whose snapshot pins a Pn3m index.
        DBInterface.execute(db, "INSERT INTO users (id, username) VALUES (1, 'alice')")
        DBInterface.execute(db, """INSERT INTO experiments (id, name, root_dir, analysis_dir)
                                   VALUES (10, 'exp', '/x', '/x/a')""")
        DBInterface.execute(db, """INSERT INTO samples (id, experiment_id, name)
                                   VALUES (100, 10, 'sA')""")
        DBInterface.execute(db, """INSERT INTO exposures
                                   (id, sample_id, filename, exposure_type, selected)
                                   VALUES (1000, 100, 'JC001', 'simple', 1)""")
        DBInterface.execute(db, """INSERT INTO comparisons
                                   (id, title, content_hash, created_by, updated_at)
                                   VALUES (1, 'Cubic vs Hex', 'h', 1, '2026-05-14T10:00:00Z')""")
        snap = JSON3.write(Dict(
            :confirmed_index => Dict(:id => 1, :phase => "Pn3m"),
            :inputs_hash => "ih1",
        ))
        DBInterface.execute(db, """INSERT INTO comparison_members
                                   (comparison_id, exposure_id, display_order,
                                    band_height, y_offset, normalization, snapshot)
                                   VALUES (1, 1000, 0, 1.0, 0.0, 'max', ?)""",
                            [snap])

        rows = HimalayaUI.comparisons_listing(db)
        @test length(rows) == 1
        r = rows[1]
        @test r[:author_username]    == "alice"
        @test r[:member_count]       == 1
        @test r[:member_phases]      == ["Pn3m"]
        @test r[:has_stale_members]  == false
        @test r[:last_event_at] isa Union{String, Nothing}
        close(db)
    end
end
```

- [ ] **Step 2: Run to verify failure.**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI, Test; include("packages/HimalayaUI/test/test_comparisons.jl")' 2>&1 | grep -E "did not pass|Test Failed" | head -5
```

Expected: failures on the five new key assertions.

- [ ] **Step 3: Extend the listing SELECTs.**

In `packages/HimalayaUI/src/comparisons.jl`, replace each of
`comparisons_for_experiment` and `comparisons_listing` SELECT bodies. The
existing SELECT already projects `last_event_at` as a COALESCE subquery.
Add four more projections.

For `comparisons_listing`, the SELECT becomes:

```julia
function comparisons_listing(db::SQLite.DB)::Vector{Dict{Symbol, Any}}
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT c.id, c.title, c.description, c.content_hash,
                  c.created_by, c.created_at, c.updated_at,
                  c.forked_from_id, c.forked_at_hash,
                  c.view_grouping_mode, c.view_show_peak_ticks, c.view_show_peak_labels,
                  COALESCE((SELECT MAX(ua.timestamp) FROM user_actions ua
                            WHERE ua.entity_type = 'comparison'
                              AND ua.entity_id = c.id), c.updated_at) AS last_event_at,
                  u.username AS author_username,
                  (SELECT COUNT(*) FROM comparison_members cm
                   WHERE cm.comparison_id = c.id) AS member_count,
                  (SELECT GROUP_CONCAT(json_extract(cm.snapshot, '$.confirmed_index.phase')
                                       || '#' || cm.display_order, '|')
                   FROM comparison_members cm
                   WHERE cm.comparison_id = c.id
                     AND json_extract(cm.snapshot, '$.confirmed_index.phase') IS NOT NULL) AS member_phases_concat,
                  EXISTS (
                    SELECT 1 FROM comparison_members cm
                    JOIN exposures e ON e.id = cm.exposure_id
                    WHERE cm.comparison_id = c.id
                      AND cm.exposure_id IS NOT NULL
                      AND json_extract(cm.snapshot, '$.inputs_hash')
                          != e.analysis_inputs_hash
                  ) AS has_stale_members
           FROM comparisons c
           LEFT JOIN users u ON u.id = c.created_by
           ORDER BY last_event_at DESC, c.id DESC"""))
    _comparison_listing_rows(rows)
end
```

Apply the same projection extension to `comparisons_for_experiment`. Full
replacement:

```julia
function comparisons_for_experiment(db::SQLite.DB, experiment_id::Integer)::Vector{Dict{Symbol, Any}}
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT DISTINCT c.id, c.title, c.description, c.content_hash,
                  c.created_by, c.created_at, c.updated_at,
                  c.forked_from_id, c.forked_at_hash,
                  c.view_grouping_mode, c.view_show_peak_ticks, c.view_show_peak_labels,
                  COALESCE((SELECT MAX(ua.timestamp) FROM user_actions ua
                            WHERE ua.entity_type = 'comparison'
                              AND ua.entity_id = c.id), c.updated_at) AS last_event_at,
                  u.username AS author_username,
                  (SELECT COUNT(*) FROM comparison_members cm2
                   WHERE cm2.comparison_id = c.id) AS member_count,
                  (SELECT GROUP_CONCAT(json_extract(cm2.snapshot, '$.confirmed_index.phase')
                                       || '#' || cm2.display_order, '|')
                   FROM comparison_members cm2
                   WHERE cm2.comparison_id = c.id
                     AND json_extract(cm2.snapshot, '$.confirmed_index.phase') IS NOT NULL) AS member_phases_concat,
                  EXISTS (
                    SELECT 1 FROM comparison_members cm2
                    JOIN exposures e2 ON e2.id = cm2.exposure_id
                    WHERE cm2.comparison_id = c.id
                      AND cm2.exposure_id IS NOT NULL
                      AND json_extract(cm2.snapshot, '$.inputs_hash')
                          != e2.analysis_inputs_hash
                  ) AS has_stale_members
           FROM comparisons c
           LEFT JOIN users u ON u.id = c.created_by
           JOIN comparison_members cm ON cm.comparison_id = c.id
           JOIN exposures e ON e.id = cm.exposure_id
           JOIN samples s ON s.id = e.sample_id
           WHERE s.experiment_id = ?
           ORDER BY last_event_at DESC, c.id DESC""", [Int(experiment_id)]))
    _comparison_listing_rows(rows)
end
```

**Notes on this SQL:**
- The outer `JOIN comparison_members cm` is the experiment-scope filter and
  produces multiple rows per comparison; `SELECT DISTINCT` collapses them.
  The aggregate subqueries use a distinct alias `cm2` so they don't
  collide with the outer `cm`.
- The `LEFT JOIN users u` is single-row (`c.created_by` is a FK), so it
  doesn't multiply rows and `DISTINCT` is still cheap.
- The phase aggregate now carries `phase || '#' || display_order` so the
  Julia side has both the phase string AND the source `display_order`
  for the §8.1 tiebreaker (frequency desc, then `display_order` asc).
  Apply the same change to `comparisons_listing` above.

- [ ] **Step 3b: Extend `forks_of_comparison` SELECT.**

`forks_of_comparison` (in `comparisons.jl`) also feeds rows into
`_comparison_listing_rows`. Once Step 4 extends the projector to read
the five new aggregates (`author_username`, `member_count`,
`member_phases_concat`, `has_stale_members`, `view_grouping_mode` /
`view_show_peak_ticks` / `view_show_peak_labels`), `forks_of_comparison`
must SELECT them too — otherwise the projector hits
`MethodError`/`UndefVarError` on missing `NamedTuple` fields. Apply the
same projection extension shown in Step 3 to `forks_of_comparison`'s
SELECT.

- [ ] **Step 3c: Pin the `forks_of_comparison` contract.**

Add a testset to `test_comparisons.jl` that fixtures a parent
comparison plus one fork, calls `HimalayaUI.forks_of_comparison(db,
parent_id)`, and asserts the new keys land on the fork row:

```julia
@testset "forks_of_comparison projects new aggregates (Compare UX A-3)" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "himalaya.db"))
        # …fixture: user, experiment, sample, exposure, parent comparison + fork
        # …with one member each (Pn3m and Hex respectively)
        forks = HimalayaUI.forks_of_comparison(db, parent_id)
        @test length(forks) == 1
        f = forks[1]
        @test f[:author_username]   == "alice"
        @test f[:member_count]      == 1
        @test f[:member_phases]     == ["Hex"]
        @test f[:has_stale_members] == false
        @test haskey(f, :view_grouping_mode)
        @test haskey(f, :last_event_at)
        close(db)
    end
end
```

Spec §10 calls this out explicitly: "Plus a `forks_of_comparison`
projection test — that helper goes through `_comparison_listing_rows`
too, so the new fields surface on fork listings as a silent side-effect.
Pin the contract."

- [ ] **Step 4: Update `_comparison_listing_rows` to surface the new fields.**

Find the helper. Its existing body iterates rows and emits a Dict for
each. Where it returns the Dict, add the five new keys. The
`member_phases_concat` field is a `|`-joined string; split into a unique-by-frequency
top-3 list with a `+N more` tail emitted client-side (we just send the
deduped list here; client decides truncation).

Pseudocode for the iteration:

```julia
function _comparison_listing_rows(rows)::Vector{Dict{Symbol, Any}}
    out = Vector{Dict{Symbol, Any}}()
    for r in rows
        # Phase list: split + dedup preserving order of first occurrence,
        # then keep top 3 by frequency.
        phases_str = ismissing(r.member_phases_concat) ? "" : String(r.member_phases_concat)
        member_phases = _topk_phases(phases_str, 3)

        push!(out, Dict{Symbol, Any}(
            :id                    => Int(r.id),
            :title                 => ismissing(r.title) ? nothing : r.title,
            :description           => ismissing(r.description) ? nothing : r.description,
            :content_hash          => ismissing(r.content_hash) ? nothing : r.content_hash,
            :created_by            => ismissing(r.created_by) ? nothing : Int(r.created_by),
            :created_at            => ismissing(r.created_at) ? nothing : r.created_at,
            :updated_at            => ismissing(r.updated_at) ? nothing : r.updated_at,
            :forked_from_id        => ismissing(r.forked_from_id) ? nothing : Int(r.forked_from_id),
            :forked_at_hash        => ismissing(r.forked_at_hash) ? nothing : r.forked_at_hash,
            :view_grouping_mode    => ismissing(r.view_grouping_mode) ? nothing : r.view_grouping_mode,
            :view_show_peak_ticks  => ismissing(r.view_show_peak_ticks) ? nothing : Bool(r.view_show_peak_ticks),
            :view_show_peak_labels => ismissing(r.view_show_peak_labels) ? nothing : Bool(r.view_show_peak_labels),
            :last_event_at         => ismissing(r.last_event_at) ? nothing : r.last_event_at,
            :author_username       => ismissing(r.author_username) ? nothing : r.author_username,
            :member_count          => Int(r.member_count),
            :member_phases         => member_phases,
            :has_stale_members     => Bool(r.has_stale_members),
        ))
    end
    out
end

"""
    _topk_phases(concat, k) -> Vector{String}

Parses a `|`-joined list of `"<phase>#<display_order>"` tokens
(NULL-filtered upstream) and returns the top-K phases by frequency.
Tiebreak by smallest `display_order` first-seen (deterministic per spec
§8.1). Empty / missing input → empty vector.
"""
function _topk_phases(concat::AbstractString, k::Integer)::Vector{String}
    isempty(concat) && return String[]
    parts = filter(!isempty, split(concat, '|'))
    counts        = Dict{String, Int}()
    first_seen_do = Dict{String, Int}()
    for p in parts
        # Each token looks like "Pn3m#3" — split on the final '#'.
        s_token = String(p)
        idx = findlast(==('#'), s_token)
        idx === nothing && continue  # malformed; skip
        phase = s_token[1:idx-1]
        do_num = tryparse(Int, s_token[idx+1:end])
        do_num === nothing && continue
        counts[phase] = get(counts, phase, 0) + 1
        if !haskey(first_seen_do, phase) || do_num < first_seen_do[phase]
            first_seen_do[phase] = do_num
        end
    end
    # Sort by (descending frequency, ascending first-seen display_order)
    phases = collect(keys(counts))
    sort!(phases, by = p -> (-counts[p], first_seen_do[p]))
    return phases[1:min(k, length(phases))]
end
```

- [ ] **Step 5: Run the tests to verify pass.**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI, Test; include("packages/HimalayaUI/test/test_comparisons.jl")' 2>&1 | grep -E "Test Summary|did not pass" | head -3
```

Expected: PASS.

- [ ] **Step 5b: Extend `test_route_response_shapes.jl`.**

The route-shape contract test at `packages/HimalayaUI/test/test_route_response_shapes.jl`
uses `assert_keys` against an EXACT key list — `assert_keys` fails if
unexpected keys appear. After A-3 the listing route returns eight new
keys (`view_grouping_mode`, `view_show_peak_ticks`, `view_show_peak_labels`,
`last_event_at`, `author_username`, `member_count`, `member_phases`,
`has_stale_members`); the existing `assert_keys` call near line 698
will fail until the expected key list is extended.

Locate the testset that asserts the GET `/api/comparisons` listing
shape (the one that constructs an `expected` Vector{Symbol} and calls
`assert_keys(row, expected)`). Extend the `expected` list with the
eight new keys:

```julia
expected = [
    :id, :title, :description, :content_hash, :created_by, :created_at,
    :updated_at, :forked_from_id, :forked_at_hash,
    :view_grouping_mode, :view_show_peak_ticks, :view_show_peak_labels,
    :last_event_at, :author_username, :member_count, :member_phases,
    :has_stale_members,
]
```

Six-layer contract testing (per `docs/contract-testing.md`): without
this, the listing-helper test (Step 1) passes while the route-shape
test fails — the same projection landed at two layers and only one
sees the keys.

- [ ] **Step 5c: Re-run both contract tests.**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI, Test; include("packages/HimalayaUI/test/test_comparisons.jl"); include("packages/HimalayaUI/test/test_route_response_shapes.jl")' 2>&1 | grep -E "Test Summary|did not pass" | head -5
```

Expected: PASS.

- [ ] **Step 6: Commit.**

```bash
git add packages/HimalayaUI/src/comparisons.jl packages/HimalayaUI/test/test_comparisons.jl packages/HimalayaUI/test/test_route_response_shapes.jl
git commit -m "feat(comparisons): listing projects author/count/phases/stale (Compare UX A-3)"
```

### Task A-4: Extend `fetch_comparison_with_members` to return view-choice fields

**Files:**
- Modify: `packages/HimalayaUI/src/comparisons.jl`
- Test: `packages/HimalayaUI/test/test_comparisons.jl`

- [ ] **Step 1: Write the failing test.**

Append to `test_comparisons.jl`:

```julia
@testset "fetch_comparison_with_members projects view_* — Compare UX A-4" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "himalaya.db"))
        DBInterface.execute(db, """INSERT INTO comparisons
                                   (id, title, content_hash, view_grouping_mode,
                                    view_show_peak_ticks, view_show_peak_labels)
                                   VALUES (1, 'x', 'h', 'byPhase', 1, 0)""")
        result = HimalayaUI.fetch_comparison_with_members(db, 1)
        @test result !== nothing
        @test result[:view_grouping_mode]    == "byPhase"
        @test result[:view_show_peak_ticks]  == true
        @test result[:view_show_peak_labels] == false
        close(db)
    end
end
```

- [ ] **Step 2: Run to verify fail.**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI, Test; include("packages/HimalayaUI/test/test_comparisons.jl")' 2>&1 | grep -E "did not pass" | head -3
```

Expected: failures on the three `view_*` assertions.

- [ ] **Step 3: Update `fetch_comparison_with_members`.**

In `comparisons.jl`, edit the SELECT inside `fetch_comparison_with_members`:

```julia
    cmp_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT id, title, description, content_hash, created_by,
                  created_at, updated_at, forked_from_id, forked_at_hash,
                  view_grouping_mode, view_show_peak_ticks, view_show_peak_labels
           FROM comparisons WHERE id = ?""", [cid]))
```

And in the Dict that the function returns, add the three projections:

```julia
        :view_grouping_mode    => ismissing(cmp.view_grouping_mode) ? nothing : cmp.view_grouping_mode,
        :view_show_peak_ticks  => ismissing(cmp.view_show_peak_ticks) ? nothing : Bool(cmp.view_show_peak_ticks),
        :view_show_peak_labels => ismissing(cmp.view_show_peak_labels) ? nothing : Bool(cmp.view_show_peak_labels),
```

(immediately before the `:members => members,` entry — preserves the
existing key order convention).

- [ ] **Step 4: Run.**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI, Test; include("packages/HimalayaUI/test/test_comparisons.jl")' 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 4b: Extend `test_route_response_shapes.jl` per-id contract — all THREE sites.**

`test_route_response_shapes.jl` has THREE `assert_keys` sites against
the closed Comparison shape. A-4 adds three view_* keys to
`fetch_comparison_with_members`'s return; each site needs the same
extension or one will regress silently:

1. **POST `/api/comparisons → full Comparison shape`** (around line 526).
2. **GET `/api/comparisons/:id`** (sibling testset, same expected list).
3. **POST `/api/comparisons/:id/submit → 409 conflict body's nested
   `current_state`** (around line 610-615). The 409 body carries
   `current_state = fetch_comparison_with_members(...)`, so the
   conflict path is a fourth caller of the same projection and
   `assert_keys` checks for *exact* key membership.

Locate each testset's `expected = [...]` and extend with the three
view_* keys (insert before `:members` to preserve the existing
key-order convention):

```julia
expected = [
    :id, :title, :description, :content_hash, :created_by, :created_at,
    :updated_at, :forked_from_id, :forked_at_hash, :forked_from_title,
    :view_grouping_mode, :view_show_peak_ticks, :view_show_peak_labels,
    :members,
]
```

For the 409 conflict-body site, the assertion is nested
(`assert_keys(body[:current_state], expected)`) — same `expected` list
applies, but make sure the `body[:current_state]` indirection is preserved.

Six-layer rule: the listing-projection layer (A-3) and the per-id
projection layer (A-4 — three call sites) BOTH have `assert_keys`
contracts; ALL FOUR sites need the extension or one regresses
silently and the round-1 reviewer's concern resurfaces.

- [ ] **Step 4c: Re-run contract tests.**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI, Test; include("packages/HimalayaUI/test/test_comparisons.jl"); include("packages/HimalayaUI/test/test_route_response_shapes.jl")' 2>&1 | grep -E "Test Summary|did not pass" | head -5
```

Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/src/comparisons.jl packages/HimalayaUI/test/test_comparisons.jl packages/HimalayaUI/test/test_route_response_shapes.jl
git commit -m "feat(comparisons): fetch_comparison_with_members projects view_* (Compare UX A-4)"
```

### Task A-5: Wire view-choice fields through submit + SSE round-trip

**Files:**
- Modify: `packages/HimalayaUI/src/routes_comparisons.jl`
- Modify: `packages/HimalayaUI/src/events.jl`
- Create: `packages/HimalayaUI/test/test_routes_comparisons.jl` (new file)
- Modify: `packages/HimalayaUI/test/runtests.jl` (register the new test file)

- [ ] **Step 0: Scaffold `test_routes_comparisons.jl` + client helpers.**

The existing test harness in `test_http.jl` exposes
`with_test_server(f, db::SQLite.DB)` where `f(port, base)` receives a
listening port + base URL — see callers across `test_routes_*.jl`.
The A-5/A-6 testsets in this plan reference `with_test_server() do client`,
`post_json(client, path, body)`, `get_json(client, path)`, and
`_minimal_member_payload()` — none of those helpers exist yet.

**Two equally-valid landing options:**

a. **Inline the canonical pattern** (recommended — matches every other
   `test_routes_*.jl` in the suite). Rewrite the testsets below to:

   ```julia
   @testset "submit echoes view_grouping_mode (Compare UX A-5)" begin
       mktempdir() do tmp
           db = open_db(joinpath(tmp, "himalaya.db"))
           # …fixture: experiment, sample, exposure, optionally a parent comparison…
           with_test_server(db) do port, base
               # Create
               resp = HTTP.post("$base/api/comparisons",
                   ["X-Username" => "alice", "Content-Type" => "application/json"],
                   JSON3.write(Dict(:title => "t", :members => [_minimal_member_payload()])))
               created = JSON3.read(resp.body, Dict{Symbol, Any})
               # Submit with view fields…
               resp2 = HTTP.post("$base/api/comparisons/$(created[:id])/submit", …)
               submitted = JSON3.read(resp2.body, Dict{Symbol, Any})
               @test submitted[:view_grouping_mode] == "byPhase"
           end
           close(db)
       end
   end
   ```

   This avoids inventing a new helper layer and uses HTTP.jl directly,
   matching the existing pattern.

b. **Define thin helpers at the top of `test_routes_comparisons.jl`:**

   ```julia
   using HTTP, JSON3, SQLite, Test
   using HimalayaUI

   function post_json(base::String, path::String, body::AbstractDict; username::String = "alice")
       resp = HTTP.post(string(base, path),
           ["X-Username" => username, "Content-Type" => "application/json"],
           JSON3.write(body); status_exception = false)
       (status = resp.status, body = JSON3.read(resp.body, Dict{Symbol, Any}))
   end

   function get_json(base::String, path::String; username::String = "alice")
       resp = HTTP.get(string(base, path), ["X-Username" => username]; status_exception = false)
       (status = resp.status, body = JSON3.read(resp.body, Dict{Symbol, Any}))
   end

   _minimal_member_payload() = Dict(:exposure_id => …, :display_order => 0,
       :band_height => 1.0, :y_offset => 0.0, :normalization => "max", :snapshot => "{}")
   ```

   Then call them inside `with_test_server(db) do port, base`:
   `post_json(base, "/api/comparisons", body)`.

**Pick option (a) — direct HTTP.post.** It avoids creating a fourth
test-helper API in the codebase and keeps the failure mode visible
(no helper between the test and the wire). Rewrite Steps 1, 6, and the
A-6 testset (and the A-5 Step 5c regression sets in
`test_route_response_shapes.jl`) to use this pattern. The
`_minimal_member_payload` helper is small enough to inline at the top
of `test_routes_comparisons.jl` as a `const` or `function`, so it can
be reused across the file.

**Register the file.** Add
`include("test_routes_comparisons.jl")` to `packages/HimalayaUI/test/runtests.jl`
in the order matching other route tests (alphabetical).

- [ ] **Step 1: Write the failing route test.**

Append a testset to `test_routes_comparisons.jl` (created in Step 0)
that posts a submit payload including the three view-choice fields
and asserts they appear on the returned JSON. **Use the canonical
pattern from Step 0** — `with_test_server(db) do port, base; HTTP.post(…)`.
The pseudocode below shows the high-level shape; flesh out fixtures
inline:

```julia
@testset "submit echoes view_grouping_mode (Compare UX A-5)" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "himalaya.db"))
        # …fixture: user, experiment, sample, exposure, _minimal_member_payload()
        #   helper at top of file returns a Dict for one member tied to the exposure…

        with_test_server(db) do port, base
            # Create a bare comparison + minimal member.
            resp = HTTP.post("$base/api/comparisons",
                ["X-Username" => "alice", "Content-Type" => "application/json"],
                JSON3.write(Dict(:title => "t", :members => [_minimal_member_payload()])))
            created = JSON3.read(resp.body, Dict{Symbol, Any})
            cid = created[:id]
            hash = created[:content_hash]

            # Submit with view fields.
            resp2 = HTTP.post("$base/api/comparisons/$(cid)/submit",
                ["X-Username" => "alice", "Content-Type" => "application/json"],
                JSON3.write(Dict(
                    :title => "t",
                    :members => [_minimal_member_payload()],
                    :expected_content_hash => hash,
                    :view_grouping_mode    => "byPhase",
                    :view_show_peak_ticks  => true,
                    :view_show_peak_labels => false,
                )))
            submitted = JSON3.read(resp2.body, Dict{Symbol, Any})
            @test submitted[:view_grouping_mode]    == "byPhase"
            @test submitted[:view_show_peak_ticks]  == true
            @test submitted[:view_show_peak_labels] == false

            # GET round-trip.
            resp3 = HTTP.get("$base/api/comparisons/$(cid)", ["X-Username" => "alice"])
            fetched = JSON3.read(resp3.body, Dict{Symbol, Any})
            @test fetched[:view_grouping_mode] == "byPhase"
        end
        close(db)
    end
end
```

`_minimal_member_payload()` is defined at the top of
`test_routes_comparisons.jl` per Step 0 (it's a small helper, not a
shared API surface).

- [ ] **Step 2: Run to verify fail.**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI, Test; include("packages/HimalayaUI/test/test_routes_comparisons.jl")' 2>&1 | grep -E "did not pass" | head -3
```

Expected: route does not yet accept or echo the fields.

- [ ] **Step 3: Update `routes_comparisons.jl` payload dict-build.**

In both POST `/api/comparisons` and POST `/api/comparisons/:id/submit`,
the handler builds a `Dict{Symbol, Any}` from the JSON body before
passing it to `apply_event!(InTransaction(), ..., payload = payload)`
inside a `with_idempotency` block. Today's code pre-extracts and
`String()`-coerces `title`, `description`, `forked_from_id`,
`forked_at_hash` into locals (see `routes_comparisons.jl:107-113`) and
then builds the payload Dict from those locals at `:225-229` (submit
handler; create handler is structurally identical). View fields are
silently dropped before they ever reach the dispatcher.

**Preserve the existing local-extraction style — add three more locals,
don't replace the pattern.** Mixing the existing
pre-coerced locals (`title`, `forked_from_id`) with inline
`get(body, ...)` calls would muddy the route-validation-vs-dispatcher-write
separation. Add `view_grouping_mode` next to the existing
String-coerced locals, plus the two booleans next to it:

```julia
# Add after the existing local extractions for title / description /
# forked_from_id / forked_at_hash (around routes_comparisons.jl:107-113).
view_grouping_mode = haskey(body, :view_grouping_mode) && body.view_grouping_mode !== nothing ?
    String(body.view_grouping_mode) : nothing
view_show_peak_ticks  = haskey(body, :view_show_peak_ticks)  ? body.view_show_peak_ticks  : nothing
view_show_peak_labels = haskey(body, :view_show_peak_labels) ? body.view_show_peak_labels : nothing
```

Then extend the existing payload Dict-build to thread them through:

```julia
payload = Dict{Symbol, Any}(
    :title => title,
    :description => description,
    :members => members,
    :view_grouping_mode    => view_grouping_mode,
    :view_show_peak_ticks  => view_show_peak_ticks,
    :view_show_peak_labels => view_show_peak_labels,
)
```

The `haskey(body, …) && body.X !== nothing` form is required so an
omitted field stays as `nothing` (dispatcher's `COALESCE`/bare-`?`
distinguishes "omitted" from "explicit null"); a present-but-null
field maps to `nothing` too because JSON3 NULL ⇒ Julia `nothing`.
Note: the Step 4 dispatcher already `String()`-coerces; doing the
coercion in BOTH the route and the dispatcher is defensive but cheap.
Apply the SAME extension to BOTH handlers (create + submit).

- [ ] **Step 4: Update `events.jl::_update_view_for_comparison_submitted!`.**

This is the submit-path handler. Locate the function. Its current shape:

1. Parses the payload member list and computes the existing/new/delete
   diff sets.
2. Loops over the diffs running INSERT / UPDATE / DELETE on
   `comparison_members`.
3. **AFTER the member-diff loop**, calls
   `new_hash = compute_content_hash(db, Int(entity_id))` to derive the
   post-diff hash.
4. Currently runs **two separate `UPDATE comparisons SET` statements** —
   one for title/description and one for `content_hash` + `updated_at`.

**Replace both UPDATEs (the ones at the END of the function, AFTER the
hash computation) with a single combined UPDATE that also writes the
view-choice columns:**

```julia
    # After the member-diff loop AND after `new_hash = compute_content_hash(...)`.
    # NOTE: `String(...)` coercion matches existing behaviour at events.jl:587-590.
    # JSON3 string subtypes must be normalised before write or the column ends up
    # holding the raw JSON3 type, which breaks hash-on-rename assertions.
    title_val = haskey(payload, :title) && payload.title !== nothing ? String(payload.title) : nothing
    desc_val  = haskey(payload, :description) && payload.description !== nothing ? String(payload.description) : nothing
    vgm  = haskey(payload, :view_grouping_mode) && payload.view_grouping_mode !== nothing ? String(payload.view_grouping_mode) : nothing
    vspt = haskey(payload, :view_show_peak_ticks)  ? payload.view_show_peak_ticks  : nothing
    vspl = haskey(payload, :view_show_peak_labels) ? payload.view_show_peak_labels : nothing

    DBInterface.execute(db,
        """UPDATE comparisons SET
              title = COALESCE(?, title),
              description = COALESCE(?, description),
              view_grouping_mode = ?,
              view_show_peak_ticks = ?,
              view_show_peak_labels = ?,
              content_hash = ?, updated_at = ?
           WHERE id = ?""",
        [title_val, desc_val, vgm, vspt, vspl, new_hash, ts, entity_id])
```

`title`/`description` use `COALESCE` (omit-keeps-current). View fields
use bare `?` — clients can null them out (e.g., "reset to default").

**Hash computation note.** `compute_content_hash` already hashes only
title / description / member identity (per spec §6.4). It does NOT
include the view-choice columns. No change to `hash.jl`.

- [ ] **Step 5: Update `events.jl::_update_view_for_comparison_created!` — BOTH branches.**

This handler has TWO branches (see `events.jl:459-487`):

- **INSERT branch (replay path)** — fires when the comparison row was
  deleted and we're replaying a `comparison_created` event into a fresh
  schema. Lines 459-465 in current code.
- **UPDATE branch (live path)** — fires when the route handler already
  minted the id via `INSERT INTO comparisons DEFAULT VALUES` at
  `routes_comparisons.jl:124` and the dispatcher is filling in title /
  description / hash. Lines 476-487 in current code.

Today's live POST `/api/comparisons` ALWAYS hits the UPDATE branch. If
we extend only the INSERT branch, the create path will silently drop
view choices in production while replay tests pass — a class-of-bug we
already burned on with the M2.2 peak/curation enrichment.

**Compute the three view-choice values once at the top of the function**
(same `String()` coercion as Step 4):

```julia
    vgm  = haskey(payload, :view_grouping_mode) && payload.view_grouping_mode !== nothing ? String(payload.view_grouping_mode) : nothing
    vspt = haskey(payload, :view_show_peak_ticks)  ? payload.view_show_peak_ticks  : nothing
    vspl = haskey(payload, :view_show_peak_labels) ? payload.view_show_peak_labels : nothing
```

**Extend the INSERT branch** (replay path) to include the three new
columns and bindings:

```julia
    DBInterface.execute(db,
        """INSERT INTO comparisons
              (id, title, description, content_hash, created_by,
               created_at, updated_at, forked_from_id, forked_at_hash,
               view_grouping_mode, view_show_peak_ticks, view_show_peak_labels)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        [entity_id, title_val, desc_val, new_hash, created_by,
         ts, ts, forked_from_id_val, forked_at_hash_val,
         vgm, vspt, vspl])
```

**Extend the UPDATE branch** (live path) — add the three view_* columns
to the SET clause:

```julia
    DBInterface.execute(db,
        """UPDATE comparisons SET
              title = ?, description = ?, content_hash = ?,
              forked_from_id = ?, forked_at_hash = ?,
              view_grouping_mode = ?,
              view_show_peak_ticks = ?,
              view_show_peak_labels = ?,
              updated_at = ?
           WHERE id = ?""",
        [title_val, desc_val, new_hash,
         forked_from_id_val, forked_at_hash_val,
         vgm, vspt, vspl,
         ts, entity_id])
```

(Bindings arranged to match the column order; cross-check against the
existing UPDATE in the code so neither branch flips a binding.)

This is the second half of A-5; without extending BOTH branches, a
`/compare/new` save that sets view choices loses them silently in
production even if replay-path tests are green.

- [ ] **Step 5b: Thread `post_state` through both broadcast call sites in `routes_comparisons.jl`.**

The submit/create handlers today call
`_enqueue_broadcast_from_result!(result, "comparison_submitted", "comparison", id)`
WITHOUT a `post_state=` kwarg (see `routes_comparisons.jl:147` for
`comparison_created` and `:234` for `comparison_submitted`). The
broadcast helper in `events.jl:235-243` defaults `post_state = nothing`,
so the emitted SSE frame omits the `post_state` envelope entirely (see
the gate at `events.jl:734`: `post_state === nothing ||
(fields[:post_state] = post_state)`).

This is the same class-of-bug as M2.2 peak/curation enrichment: a route
emits a successful response while the SSE frame carries no payload, and
foreign clients cannot reconcile their caches without an extra refetch.
A-6's test in the next task asserts `post_state[:view_grouping_mode]
== "distinct"` and A-9's cache-shape contract requires a `post_state`
to splice from — both fail unless this step lands.

**Change both call sites** to compute the projection and pass it through:

```julia
# routes_comparisons.jl — POST /api/comparisons (around line 147)
post = fetch_comparison_with_members(db, id)
# `fetch_comparison_with_members` returns Union{Dict, Nothing}; post-write
# inside the same tx it MUST be a Dict, otherwise an earlier step
# silently failed. Surface the invariant violation rather than emitting
# a broken SSE frame.
post === nothing && error("post-write fetch_comparison_with_members returned nothing for id=$(id)")
_enqueue_broadcast_from_result!(result, "comparison_created", "comparison", id;
                                post_state = post)

# routes_comparisons.jl — POST /api/comparisons/:id/submit (around line 234)
post = fetch_comparison_with_members(db, id)
post === nothing && error("post-write fetch_comparison_with_members returned nothing for id=$(id)")
_enqueue_broadcast_from_result!(result, "comparison_submitted", "comparison", id;
                                post_state = post)
```

**Where to call `fetch_comparison_with_members`.** Inside the
`with_idempotency` block (i.e., before `_enqueue_broadcast_from_result!`
returns), so the projection is read against the same transaction as the
writes. Reading post-commit would race other writers under
`parallel = true`. The `_enqueue_broadcast_from_result!` helper queues
the broadcast for the post-commit flush — it doesn't emit immediately —
so this stays compatible with InTransaction non-broadcast semantics.

**Rollback semantics of the `nothing`-guard.** `error(...)` inside
`with_idempotency` aborts the savepoint and bubbles out of the route
handler as a 500. That IS the desired failure mode here — we'd rather
roll back the write than emit a broken SSE frame to every connected
client. The trade-off is that any post-commit `nothing` return (which
should be impossible under a single tx) discards the write entirely.
If implementation surfaces a real reason this fires, surface the id
in the error message (already done) so the failure is debuggable from
the route log.

**Why this used to be a contingency.** The earlier draft of A-6 listed
this as "if the test fails, also update the broadcast site." Promoting
to an explicit A-5 step removes the contingent path: A-6 tests the
broadcast contract; A-9 splices from the broadcast payload; both depend
on `post_state` being present, and that depends on this step.

- [ ] **Step 5c: Pin HTTP response body to `fetch_comparison_with_members` — both POST handlers.**

A-9 splices the SSE `post_state` into `queryKeys.comparison(cid)`.
`saveComparisonMutator.onSuccess` writes the HTTP response into the
same cache key. If the HTTP body shape diverges from the SSE
`post_state` shape (different field set, different `forked_from_title`
resolution, etc.), an HTTP-wins race will clobber an SSE-correct row.

**Make both POST handlers serialize `fetch_comparison_with_members(db,
id)` as the HTTP body.** Today each handler builds a small inline Dict
and passes it to `JSON3.write(...)` inside `HTTP.Response(...)`.
Replace that inline Dict with the projection so HTTP and SSE write
identical shapes:

```julia
# routes_comparisons.jl — POST /api/comparisons (around line 130-150)
with_idempotency(db, req) do
    # …existing INSERT INTO comparisons DEFAULT VALUES; capture rowid
    #   into local `id` via `DBInterface.lastrowid(res)` BEFORE the
    #   apply_event! call (the dispatcher's `entity_id = id` kwarg
    #   requires it). Then call apply_event!(InTransaction(), ..., entity_id = id)
    #   to get the `result` NamedTuple…
    post = fetch_comparison_with_members(db, id)
    post === nothing && error("post-write fetch_comparison_with_members returned nothing for id=$(id)")
    # `_enqueue_broadcast_from_result!` MUST receive the `result`
    # NamedTuple — it reads .event_id / .user_id / .client_id /
    # .client_op_id / .payload_json off its first arg (events.jl:235-243).
    # The projection rides as the post_state kwarg.
    _enqueue_broadcast_from_result!(result, "comparison_created", "comparison", id;
                                    post_state = post)
    # `with_idempotency`'s do-block MUST return an HTTP.Response (see
    # idempotency.jl:83; idempotency.jl:128 reads .status and .body
    # off the return value). Serialise the projection as the body.
    return HTTP.Response(201, ["Content-Type" => "application/json"], JSON3.write(post))
end
```

(Apply the same change to the submit handler: replace its inline
response Dict with `post = fetch_comparison_with_members(db, id)`,
guard against `nothing`, broadcast via `result`, return
`HTTP.Response(200, …, JSON3.write(post))`.)

**API gotchas (load-bearing — getting these wrong throws at runtime):**
- The first arg of `_enqueue_broadcast_from_result!` is the **`result`
  NamedTuple from `apply_event!`**, NOT the post_state projection.
  Passing `post` directly would trigger `getproperty` errors on
  missing `.event_id` etc.
- The do-block return is **`HTTP.Response`**, NOT a NamedTuple. Every
  existing caller in `routes_comparisons.jl` follows this pattern;
  match it.

**Why this is correct.** Both write paths now use the same projection;
HTTP-wins and SSE-wins races converge on identical cache rows.
`forked_from_title` and other server-computed fields survive either
path. Replay-as-rerun is unaffected (still rolls back on the same
cache key).

**Regression test.** Append to `test_route_response_shapes.jl`.

**Scaffolding note.** The pseudocode below uses `with_test_server() do client`
+ `post_json` / `get_json` as SHORTHAND. The actual implementation uses
the canonical pattern from A-5 Step 0: `mktempdir() do tmp; db = open_db(...);
with_test_server(db) do port, base; HTTP.post("$base/api/...",
["X-Username" => "alice", "Content-Type" => "application/json"],
JSON3.write(body)); end; close(db); end`. The two A-5 regression sets
land in `test_route_response_shapes.jl` (NOT `test_routes_comparisons.jl`)
because `assert_keys` is a file-local helper in that file; co-locating
keeps the helper in scope. Apply the same scaffolding rewrite to A-5
Step 6 and the A-6 SSE test.

```julia
@testset "submit response body uses fetch_comparison_with_members shape (Compare UX A-5)" begin
    # …mktempdir + open_db + fixture (user/experiment/sample/exposure)…
    with_test_server(db) do port, base
        # Create + submit with view_* via HTTP.post; deserialise via JSON3.read.
        resp = HTTP.post("$base/api/comparisons", …, JSON3.write(Dict(:title => "t",
            :members => [_minimal_member_payload()])))
        created = JSON3.read(resp.body, Dict{Symbol, Any})
        hash = created[:content_hash]

        resp2 = HTTP.post("$base/api/comparisons/$(created[:id])/submit", …,
            JSON3.write(Dict(:title => "t", :members => [_minimal_member_payload()],
                :expected_content_hash => hash, :view_grouping_mode => "byPhase")))
        submitted = JSON3.read(resp2.body, Dict{Symbol, Any})

        # Body has the same key set as the per-id GET (per A-4 Step 4b).
        expected = [
            :id, :title, :description, :content_hash, :created_by, :created_at,
            :updated_at, :forked_from_id, :forked_at_hash, :forked_from_title,
            :view_grouping_mode, :view_show_peak_ticks, :view_show_peak_labels,
            :members,
        ]
        assert_keys(submitted, expected)
        @test submitted[:view_grouping_mode] == "byPhase"
    end
end

@testset "create response body uses fetch_comparison_with_members shape (Compare UX A-5)" begin
    # …same mktempdir + with_test_server(db) do port, base… scaffolding…
    resp = HTTP.post("$base/api/comparisons", …,
        JSON3.write(Dict(:title => "fresh", :members => [_minimal_member_payload()],
            :view_grouping_mode => "distinct")))
    created = JSON3.read(resp.body, Dict{Symbol, Any})
    expected = [
        :id, :title, :description, :content_hash, :created_by, :created_at,
        :updated_at, :forked_from_id, :forked_at_hash, :forked_from_title,
        :view_grouping_mode, :view_show_peak_ticks, :view_show_peak_labels,
        :members,
    ]
    assert_keys(created, expected)
    @test created[:view_grouping_mode] == "distinct"
end
```

Without this step the queue reviewer's round-2 concern stands: own-op
HTTP `onSuccess` could overwrite an SSE-correct cache row with an
HTTP-incomplete row whenever the route's response Dict isn't a strict
superset of `fetch_comparison_with_members`.

- [ ] **Step 6: Extend the test to cover the create path too.**

In the same testset (or a sibling), POST `/api/comparisons` with view
fields set and assert they ride through to GET. The create path was a
"Should-fix" callout from the plan reviewer — without an explicit test,
the Step 5 INSERT change can silently regress.

```julia
@testset "create path echoes view fields (Compare UX A-5b)" begin
    # Scaffolding: same mktempdir + open_db + with_test_server(db) do port, base
    # pattern as A-5 Step 0. Pseudocode below uses shorthand.
    resp = HTTP.post("$base/api/comparisons", …,
        JSON3.write(Dict(
            :title => "fresh", :members => [_minimal_member_payload()],
            :view_grouping_mode    => "byPhase",
            :view_show_peak_ticks  => false,
            :view_show_peak_labels => true,
        )))
    created = JSON3.read(resp.body, Dict{Symbol, Any})
    @test created[:view_grouping_mode]    == "byPhase"
    @test created[:view_show_peak_ticks]  == false
    @test created[:view_show_peak_labels] == true
    resp2 = HTTP.get("$base/api/comparisons/$(created[:id])", ["X-Username" => "alice"])
    fetched = JSON3.read(resp2.body, Dict{Symbol, Any})
    @test fetched[:view_grouping_mode] == "byPhase"
end
```

- [ ] **Step 7: Run to verify pass.**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI, Test; include("packages/HimalayaUI/test/test_routes_comparisons.jl")' 2>&1 | grep -E "Test Summary|did not pass" | tail -5
```

Expected: PASS.

- [ ] **Step 8: Commit.**

```bash
git add packages/HimalayaUI/src/routes_comparisons.jl packages/HimalayaUI/src/events.jl packages/HimalayaUI/test/test_routes_comparisons.jl
git commit -m "feat(comparisons): submit + create accept + echo view_* fields (Compare UX A-5)"
```

### Task A-6: Pin contract — SSE frame includes view-choice fields

**Files:**
- Test: `packages/HimalayaUI/test/test_idempotency_replay_invariant.jl`

- [ ] **Step 1: Add an SSE-capture testset.**

Use the in-process subscriber pattern from CLAUDE.md "In-process SSE
subscriber for testing". `_capture_sse_during` returns a NamedTuple with
a `.pending` Channel of raw frame strings; the test must `take!` and
JSON-parse to inspect the post_state.

```julia
using JSON3  # if not already imported at the top of the file

@testset "SSE frame includes view_* fields on comparison_submitted (Compare UX A-6)" begin
    # Scaffolding: same canonical pattern as A-5 Step 0 — mktempdir +
    # open_db + with_test_server(db) do port, base. Pseudocode below
    # abbreviates fixture + setup.
    resp = HTTP.post("$base/api/comparisons", ["X-Username" => "alice", "Content-Type" => "application/json"],
        JSON3.write(Dict(:title => "t", :members => [_minimal_member_payload()])))
    created = JSON3.read(resp.body, Dict{Symbol, Any})
    cid = created[:id]
    hash = created[:content_hash]

    frames = _capture_sse_during("comparison_submitted") do
        HTTP.post("$base/api/comparisons/$(cid)/submit",
            ["X-Username" => "alice", "Content-Type" => "application/json"],
            JSON3.write(Dict(
                :title => "t",
                :members => [_minimal_member_payload()],
                :expected_content_hash => hash,
                :view_grouping_mode    => "distinct",
                :view_show_peak_ticks  => false,
                :view_show_peak_labels => true,
            )))
    end

    # `_capture_sse_during` returns a drained `Vector{String}` of
    # captured frames (and closes its internal channel before
    # returning — see test_idempotency_replay_invariant.jl:41-65).
    # Each frame is shaped "event: curation\ndata: {json}\n\n".
    @test length(frames) == 1                              # exactly one broadcast per op
    data_line = match(r"^data:\s*(.+)$"m, frames[1]).captures[1]
    parsed = JSON3.read(data_line, Dict{Symbol, Any})
    @test parsed[:kind] == "comparison_submitted"
    post_state = parsed[:post_state]
    @test post_state[:view_grouping_mode]    == "distinct"
    @test post_state[:view_show_peak_ticks]  == false
    @test post_state[:view_show_peak_labels] == true
end
```

- [ ] **Step 2: Run.**

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI, Test; include("packages/HimalayaUI/test/test_idempotency_replay_invariant.jl")' 2>&1 | tail -5
```

Expected: PASS. The `post_state` envelope was threaded through both
broadcast call sites in A-5 Step 5b — this test is the contract that
pins it stays that way.

**Note on the exactly-one-broadcast pin.** The `@test length(frames)
== 1` assertion in Step 1 above IS the second-frame guard — drives
directly off `_capture_sse_during`'s drained-Vector return. A naive
`!isready(sub.pending)` immediately after a `take!` would be racey
(the post-commit queue pushes onto the channel from a separate task;
a second in-flight push could land milliseconds later, returning a
false PASS), but the helper already handles the settle delay
internally and closes the channel before returning, so the Vector
length is the deterministic signal.

Without this assertion, a regression that adds a redundant
`_enqueue_broadcast_from_result!` call (or a future helper that
broadcasts inside InTransaction) silently produces a duplicate frame
that the test as written would not catch — own-op deduplication on
the client masks the duplication until it bites a foreign tab.

- [ ] **Step 3: Commit.**

```bash
git add packages/HimalayaUI/test/test_idempotency_replay_invariant.jl
git commit -m "test(sse): view_* fields land in comparison_submitted post_state (Compare UX A-6)"
```

### Task A-7: Extend the TypeScript types in `api.ts`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/api.ts`
- Test: (compile-time only; `npm run build`)

- [ ] **Step 1: Extend `ComparisonSummary`.**

Edit `packages/HimalayaUI/frontend/src/api.ts`. Replace the existing
interface with:

```ts
export interface ComparisonSummary {
  id: number;
  title: string;
  description: string | null;
  content_hash: string;
  created_by: number | null;
  created_at: string | null;
  updated_at: string | null;
  forked_from_id: number | null;
  forked_at_hash: string | null;
  view_grouping_mode: string | null;
  view_show_peak_ticks: boolean | null;
  view_show_peak_labels: boolean | null;
  last_event_at: string | null;
  author_username: string | null;
  member_count: number;
  member_phases: string[];
  has_stale_members: boolean;
}
```

- [ ] **Step 2: Extend `Comparison`.**

Add the three `view_*` fields and `last_event_at` to the `Comparison`
interface (it doesn't need `author_username` / `member_count` /
`member_phases` / `has_stale_members` because the full GET already
includes `members[]`, but `view_*` and `last_event_at` should be present
for consistency):

```ts
export interface Comparison {
  id: number;
  title: string;
  description: string | null;
  content_hash: string;
  created_by: number | null;
  created_at: string | null;
  updated_at: string | null;
  forked_from_id: number | null;
  forked_at_hash: string | null;
  forked_from_title: string | null;
  view_grouping_mode: string | null;
  view_show_peak_ticks: boolean | null;
  view_show_peak_labels: boolean | null;
  last_event_at: string | null;
  members: ComparisonMember[];
}
```

- [ ] **Step 3: Extend `SaveComparisonBody` to accept view fields.**

Locate the existing `SaveComparisonBody` type. Add:

```ts
export interface SaveComparisonBody {
  title: string;
  description?: string;
  members: ComparisonMemberInput[];
  expected_content_hash?: string;
  forked_from_id?: number | null;
  forked_at_hash?: string | null;
  view_grouping_mode?: "bySample" | "byPhase" | "distinct" | null;
  view_show_peak_ticks?: boolean | null;
  view_show_peak_labels?: boolean | null;
}
```

- [ ] **Step 4: Verify build.**

```bash
cd packages/HimalayaUI/frontend && npm run build 2>&1 | tail -20
```

Expected: TypeScript compiles. Any existing call site that destructures
`ComparisonSummary` and treats it as a closed shape may need a small
update — fix as the type-checker flags.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/api.ts
git commit -m "types(api): extend Comparison(Summary) with view_*, projection fields (Compare UX A-7)"
```

### Task A-8: Cache-shape contract test for the new listing fields

**Files:**
- Test: `packages/HimalayaUI/frontend/test/queue/cache-shape.test.ts`

- [ ] **Step 1: Add a cache-shape assertion.**

```ts
import { describe, it, expect } from "vitest";
import type { ComparisonSummary } from "../../src/api";

describe("comparison listing cache shape — Compare UX A-8", () => {
  it("includes the new projection fields", () => {
    const row: ComparisonSummary = {
      id: 1,
      title: "x",
      description: null,
      content_hash: "h",
      created_by: 1,
      created_at: null,
      updated_at: null,
      forked_from_id: null,
      forked_at_hash: null,
      view_grouping_mode: null,
      view_show_peak_ticks: null,
      view_show_peak_labels: null,
      last_event_at: "2026-05-14T10:00:00Z",
      author_username: "alice",
      member_count: 3,
      member_phases: ["Pn3m", "Hex"],
      has_stale_members: false,
    };
    // Compile-time anchor; surface a value-time check too:
    expect(row.member_count).toBe(3);
    expect(row.member_phases).toEqual(["Pn3m", "Hex"]);
    expect(row.author_username).toBe("alice");
    // Pin each new view_* field at the value layer — protects against a
    // future refactor that renames (e.g. viewGroupingMode) while leaving
    // the old name aliased; compile-time-only checks slip through.
    expect(row.view_grouping_mode).toBeNull();
    expect(row.view_show_peak_ticks).toBeNull();
    expect(row.view_show_peak_labels).toBeNull();
    expect(row.last_event_at).toBe("2026-05-14T10:00:00Z");
    expect(row.has_stale_members).toBe(false);
  });

  it("accepts populated view_* values too", () => {
    const row: ComparisonSummary = {
      id: 2, title: "y", description: null, content_hash: "h2",
      created_by: 1, created_at: null, updated_at: null,
      forked_from_id: null, forked_at_hash: null,
      view_grouping_mode: "byPhase",
      view_show_peak_ticks: true,
      view_show_peak_labels: false,
      last_event_at: null, author_username: null,
      member_count: 0, member_phases: [], has_stale_members: false,
    };
    expect(row.view_grouping_mode).toBe("byPhase");
    expect(row.view_show_peak_ticks).toBe(true);
    expect(row.view_show_peak_labels).toBe(false);
  });
});
```

- [ ] **Step 2: Run.**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/queue/cache-shape.test.ts
```

Expected: PASS.

- [ ] **Step 3: Commit.**

```bash
git add packages/HimalayaUI/frontend/test/queue/cache-shape.test.ts
git commit -m "test(cache): pin ComparisonSummary projection shape (Compare UX A-8)"
```

### Task A-9: applyRemoteToCache splices `post_state` into the comparison cache

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts`
- Test: `packages/HimalayaUI/frontend/test/queue/applyRemoteToCache.compare.test.ts` (new)

**Why this is not a passthrough whitelist tweak.** The current
`comparison_submitted` branch (see `applyRemoteToCache.ts:207-221`)
deliberately INVALIDATES `queryKeys.comparison(id)` +
`queryKeys.comparisonMembers(id)` + `["comparisons"]` rather than
splicing from the payload — the existing comment explains that the
response shape (with `is_stale`, `forked_from_title`, …) is
server-computed and replicating it client-side would drift.

A-5 Step 5b now ships `post_state = fetch_comparison_with_members(db,
id)` on every `comparison_submitted` / `comparison_created` SSE frame.
That `post_state` IS the same projection the GET endpoint returns
(including the server-computed `forked_from_title` plus the new view_*
fields), so splicing is now safe AND required — A-9's contract test
asserts `qc.getQueryData(queryKeys.comparison(id)).view_grouping_mode
== "byPhase"`, which can only land if the splice writes through.

- [ ] **Step 1: New test file.**

```ts
import { describe, it, expect } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import { queryKeys } from "../../src/queries";
import type { Comparison } from "../../src/api";

describe("applyRemoteToCache: comparison_submitted view_* — Compare UX A-9", () => {
  it("writes view_* fields into the cache", () => {
    const qc = new QueryClient();
    const cid = 7;
    const baseline: Comparison = {
      id: cid, title: "x", description: null, content_hash: "h0",
      created_by: 1, created_at: null, updated_at: null,
      forked_from_id: null, forked_at_hash: null, forked_from_title: null,
      view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
      last_event_at: null, members: [],
    };
    qc.setQueryData(queryKeys.comparison(cid), baseline);

    applyRemoteToCache(qc, {
      id: 99, kind: "comparison_submitted",
      entity_type: "comparison", entity_id: cid,
      actor: "bob", client_id: "tab2", client_op_id: "op2",
      ts: "2026-05-14T10:00:00Z",
      payload: null,
      post_state: {
        ...baseline,
        title: "y",
        content_hash: "h1",
        view_grouping_mode: "byPhase",
        view_show_peak_ticks: true,
        view_show_peak_labels: false,
      },
    });

    const updated = qc.getQueryData<Comparison>(queryKeys.comparison(cid));
    expect(updated?.view_grouping_mode).toBe("byPhase");
    expect(updated?.view_show_peak_ticks).toBe(true);
    expect(updated?.view_show_peak_labels).toBe(false);
  });
});
```

- [ ] **Step 2: Run.** Expected: FAIL — the current branch invalidates,
so the cache read after `applyRemoteToCache` returns the baseline (not
the post_state with `view_grouping_mode: "byPhase"`).

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/queue/applyRemoteToCache.compare.test.ts
```

- [ ] **Step 3: Switch the `comparison_submitted` / `comparison_created` branch from invalidate to splice.**

Open `packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts`.
Find the existing branch (around line 207-221). Replace the
invalidate-only logic with:

```ts
case "comparison_submitted":
case "comparison_created": {
  const cid = remote.entity_id;
  if (remote.post_state != null) {
    // post_state IS fetch_comparison_with_members(db, cid) — the same
    // projection /api/comparisons/:id returns, including server-computed
    // forked_from_title and the new view_* fields (per A-4 + A-5 step 5b).
    // Splicing is safe; we can drop the cache-invalidation round-trip.
    const post = remote.post_state as Comparison;  // narrow once
    qc.setQueryData(queryKeys.comparison(cid), post);
    if (Array.isArray(post.members)) {
      qc.setQueryData(queryKeys.comparisonMembers(cid), post.members);
    }
  } else {
    // Fallback: pre-A-5 frame without post_state. Keep the safety net so
    // a partial deploy can't strand the cache.
    qc.invalidateQueries({ queryKey: queryKeys.comparison(cid) });
    qc.invalidateQueries({ queryKey: queryKeys.comparisonMembers(cid) });
  }
  // The listing cache always invalidates — too many denormalised fields
  // (member_count, member_phases, last_event_at) for a manual splice.
  qc.invalidateQueries({ queryKey: ["comparisons"] });
  break;
}
```

- [ ] **Step 3b: Add a `forked_from_title` survival assertion to the test.**

Extend the testset to verify the server-computed field is preserved
across the splice — this is the contract the previous invalidate-only
design existed to protect:

```ts
it("preserves forked_from_title through the splice", () => {
  const qc = new QueryClient();
  const cid = 7;
  const baseline: Comparison = {
    id: cid, title: "x", description: null, content_hash: "h0",
    created_by: 1, created_at: null, updated_at: null,
    forked_from_id: 3, forked_at_hash: "h_parent",
    forked_from_title: "Parent",  // server-computed
    view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
    last_event_at: null, members: [],
  };
  qc.setQueryData(queryKeys.comparison(cid), baseline);
  applyRemoteToCache(qc, {
    id: 99, kind: "comparison_submitted",
    entity_type: "comparison", entity_id: cid,
    actor: "bob", client_id: "tab2", client_op_id: "op2",
    ts: "2026-05-14T10:00:00Z", payload: null,
    post_state: { ...baseline, title: "renamed", content_hash: "h1" },
  });
  const updated = qc.getQueryData<Comparison>(queryKeys.comparison(cid));
  expect(updated?.forked_from_title).toBe("Parent");
  expect(updated?.title).toBe("renamed");
});
```

- [ ] **Step 3c: Run.** Expected: PASS.

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/queue/applyRemoteToCache.compare.test.ts
```

- [ ] **Step 3d: Pin whole-row-overwrite ordering for the comparison cache key.**

Comparison is the first whole-row overwrite added to
`applyRemoteToCache` — peaks/indices use list-merges. The test below
exercises the bare-sequential-write contract (optimistic write,
then foreign splice, then final-state convergence). It does NOT
exercise `replayCoordinator.handleRemoteEvent` directly — that
genuine replay-as-rerun coverage lives in `test/queue/replayCoordinator.test.ts`
and stays the source of truth for the rollback-then-re-run-onMutate
invariant. This test adds the comparison-specific row-shape pin.

A follow-up Should-fix worth promoting to a separate step (not done
here): drive a parallel test through `useQueueMutation.mutate(...)`
with a real mid-flight queue so `replayCoordinator.handleRemoteEvent`
fires and the optimistic flicker can be measured. Tracked in the
plan's Out of scope list as "queue interleave coverage for
comparison whole-row overwrite."

```ts
describe("applyRemoteToCache comparison whole-row overwrite — Compare UX A-9", () => {
  it("foreign-wins final state with an optimistic own-op write in the cache", () => {
    const qc = new QueryClient();
    const cid = 7;
    const baseline: Comparison = {
      id: cid, title: "x", description: null, content_hash: "h0",
      created_by: 1, created_at: null, updated_at: null,
      forked_from_id: null, forked_at_hash: null, forked_from_title: null,
      view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
      last_event_at: null, members: [],
    };
    qc.setQueryData(queryKeys.comparison(cid), baseline);

    // 1) Own-op onMutate writes optimistic title "x-mine" into cache.
    qc.setQueryData(queryKeys.comparison(cid), { ...baseline, title: "x-mine" });

    // 2) Foreign frame arrives.
    applyRemoteToCache(qc, {
      id: 99, kind: "comparison_submitted",
      entity_type: "comparison", entity_id: cid,
      actor: "bob", client_id: "tab2", client_op_id: "op-foreign",
      ts: "2026-05-14T10:00:00Z", payload: null,
      post_state: { ...baseline, title: "x-foreign", content_hash: "h1" },
    });

    // 3) Final state — foreign-wins. NOTE: own-op write at step (1) is
    //    NOT driven through useQueueMutation here; this test only pins
    //    the sequential-write contract for the new whole-row cache key.
    //    Genuine rollback/re-run-onMutate coverage against this cache
    //    shape lives in replayCoordinator.test.ts.
    const after = qc.getQueryData<Comparison>(queryKeys.comparison(cid));
    expect(after?.title).toBe("x-foreign");
    expect(after?.content_hash).toBe("h1");
  });
});
```

If `replayCoordinator`'s rollback path then runs against this cache
and the own-op's `restore` closure tries to revert to "x-mine", the
queue's "re-run onMutate forward in insertion order" pass rebuilds
the optimistic state from the new baseline. Final-state convergence
is preserved; this test pins the intermediate write order.

- [ ] **Step 4: Commit.**

```bash
git add packages/HimalayaUI/frontend/test/queue/applyRemoteToCache.compare.test.ts packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts
git commit -m "feat(queue): splice post_state into comparison cache (Compare UX A-9)"
```

### Task A-10: Save mutator includes view fields in request body

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/saveComparison.ts`
- Test: `packages/HimalayaUI/frontend/test/lib/queue/mutators/saveComparison.test.ts` (new)

- [ ] **Step 1: Failing test.**

```ts
import { describe, it, expect, vi } from "vitest";
import { saveComparisonMutator } from "../../../../src/lib/queue/mutators/saveComparison";

describe("saveComparison passes view_* into request body — Compare UX A-10", () => {
  it("forwards view_grouping_mode and friends", async () => {
    const spy = vi.spyOn(globalThis, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ id: 1, members: [], content_hash: "h" }),
                   { status: 200, headers: { "content-type": "application/json" } }));
    await saveComparisonMutator.request({
      id: 1, title: "t", members: [],
      expected_content_hash: "h0",
      view_grouping_mode: "byPhase",
      view_show_peak_ticks: true,
      view_show_peak_labels: false,
      username: "alice", clientId: "tab", clientOpId: "op",
    });
    const init = spy.mock.calls[0]?.[1];
    const body = JSON.parse(String(init?.body ?? "{}"));
    expect(body.view_grouping_mode).toBe("byPhase");
    expect(body.view_show_peak_ticks).toBe(true);
    expect(body.view_show_peak_labels).toBe(false);
    spy.mockRestore();
  });
});
```

- [ ] **Step 2: Run.**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/lib/queue/mutators/saveComparison.test.ts
```

Expected: fail — the mutator's `buildBody` doesn't yet include view fields.

- [ ] **Step 3: Update `saveComparison.ts`.**

Extend the `SaveComparisonInput` type and `buildBody`:

```ts
export interface SaveComparisonInput {
  id?: number;
  title: string;
  description?: string | null;
  members: ComparisonMemberInput[];
  expected_content_hash?: string;
  forked_from_id?: number | null;
  forked_at_hash?: string | null;
  view_grouping_mode?: "bySample" | "byPhase" | "distinct" | null;
  view_show_peak_ticks?: boolean | null;
  view_show_peak_labels?: boolean | null;
}

function buildBody(p: SaveComparisonInput): SaveComparisonBody {
  const out: SaveComparisonBody = { title: p.title, members: p.members };
  if (p.description !== undefined) out.description = p.description;
  if (p.expected_content_hash !== undefined) out.expected_content_hash = p.expected_content_hash;
  if (p.forked_from_id !== undefined) out.forked_from_id = p.forked_from_id;
  if (p.forked_at_hash !== undefined) out.forked_at_hash = p.forked_at_hash;
  if (p.view_grouping_mode    !== undefined) out.view_grouping_mode    = p.view_grouping_mode;
  if (p.view_show_peak_ticks  !== undefined) out.view_show_peak_ticks  = p.view_show_peak_ticks;
  if (p.view_show_peak_labels !== undefined) out.view_show_peak_labels = p.view_show_peak_labels;
  return out;
}
```

- [ ] **Step 4: Run.**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/lib/queue/mutators/saveComparison.test.ts
```

Expected: PASS.

- [ ] **Step 4b: Layer-6 (`onSuccess`) test — HTTP response writes view_* into the cache.**

Six-layer contract testing (per `docs/contract-testing.md`): A-9 pins
the SSE-driven path; this step pins the HTTP-response-wins-race path.
Without it, a future refactor that drops view_* from the API response
type would silently omit them from `qc.setQueryData(queryKeys.comparison(id), response)`
until the SSE frame arrived.

Append to `saveComparison.test.ts`:

```ts
import { QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../../../../src/queries";

describe("saveComparison onSuccess reconciliation — Compare UX A-10", () => {
  it("writes view_* fields AND members into both cache keys from the HTTP response", () => {
    const qc = new QueryClient();
    const fakeMembers = [{ exposure_id: 100, display_order: 0 }] as any;
    const fakeResponse = {
      id: 1, title: "t", description: null, content_hash: "h1",
      created_by: 1, created_at: null, updated_at: null,
      forked_from_id: null, forked_at_hash: null, forked_from_title: null,
      view_grouping_mode: "byPhase",
      view_show_peak_ticks: true,
      view_show_peak_labels: false,
      last_event_at: null, members: fakeMembers,
    };
    saveComparisonMutator.onSuccess(
      { id: 1, title: "t", members: [],
        view_grouping_mode: "byPhase",
        view_show_peak_ticks: true,
        view_show_peak_labels: false,
        username: "alice", clientId: "tab", clientOpId: "op" },
      fakeResponse,
      qc,
    );
    // Layer-6 contract: BOTH cache keys must reflect the response.
    // saveComparison.ts:84 already writes both; this assertion pins
    // it so a future refactor that drops one write is caught here.
    const cachedComparison = qc.getQueryData(queryKeys.comparison(1));
    expect((cachedComparison as typeof fakeResponse).view_grouping_mode).toBe("byPhase");
    expect((cachedComparison as typeof fakeResponse).view_show_peak_ticks).toBe(true);
    expect((cachedComparison as typeof fakeResponse).view_show_peak_labels).toBe(false);
    const cachedMembers = qc.getQueryData(queryKeys.comparisonMembers(1));
    expect(cachedMembers).toEqual(fakeMembers);
  });
});
```

If `saveComparisonMutator.onSuccess` currently strips fields from the
response when writing into the cache, fix it as part of this step.

- [ ] **Step 4c: Mutator-registry resolution test.**

Replay-as-rerun depends on `mutatorRegistry.ts` resolving an OpPayload
back to the right mutator. Build an OpPayload that includes view_*
fields and confirm `resolveMutator(payload)` returns
`saveComparisonMutator`:

```ts
import { resolveMutator } from "../../../src/lib/queue/mutatorRegistry";

describe("mutatorRegistry resolves saveComparison with view_* — Compare UX A-10", () => {
  it("returns saveComparisonMutator for a save payload carrying view_*", () => {
    const payload = {
      op_kind: "comparison_submitted" as const,  // or whatever the registry keys on
      input: {
        id: 1, title: "t", members: [],
        view_grouping_mode: "byPhase",
        view_show_peak_ticks: true,
        view_show_peak_labels: false,
        username: "alice", clientId: "tab", clientOpId: "op",
      },
    };
    expect(resolveMutator(payload)).toBe(saveComparisonMutator);
  });
});
```

Confirms the registry's discriminator (likely keyed on op_kind, not
input shape) still routes correctly after the input type widens. If
the registry path doesn't match this shape, adapt to the actual API —
the goal is to pin "view_* widening does not break replay routing."

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/saveComparison.ts packages/HimalayaUI/frontend/test/lib/queue/mutators/saveComparison.test.ts
git commit -m "feat(queue): saveComparison rides view_* fields (Compare UX A-10)"
```

### Task A-11: Phase A regression sweep

**Files:** (none — runs the suites)

- [ ] **Step 1: Run Julia suite.**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-A.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-A.out | tail -20
```

Expected: clean.

- [ ] **Step 2: Run Vitest suite.**

```bash
cd packages/HimalayaUI/frontend && npm test -- --run > /tmp/vitest-A.out 2>&1
grep -E "Tests:|Failed|passed" /tmp/vitest-A.out | tail -5
```

Expected: clean.

- [ ] **Step 3: Build.**

```bash
cd packages/HimalayaUI/frontend && npm run build 2>&1 | tail -10
```

Expected: clean.

If any of A-1…A-10 has regressed, fix the regression in a follow-up commit
before starting Phase B.

---

## Phase B — Route flattening

Phase goal: drop the `/edit` URL segment. `/compare/:id` is the only route
shape going forward. The merge into a single component happens in Phase C;
in Phase B we just redirect and tidy `comparePath()`.

### Task B-1: Add the redirect Route

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/AppShell.tsx`
- Test: `packages/HimalayaUI/frontend/test/compareRouteFlatten.test.tsx` (new)

- [ ] **Step 1: Failing route test.**

```tsx
import { describe, it, expect } from "vitest";
import { render } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import { AppShell } from "../src/components/AppShell";

function LocationProbe(): JSX.Element {
  const loc = useLocation();
  return <div data-testid="probe">{loc.pathname}</div>;
}

describe("Compare /edit redirect — Compare UX B-1", () => {
  it("redirects /experiments/:eid/compare/:id/edit → /experiments/:eid/compare/:id", async () => {
    const { findByTestId } = render(
      <MemoryRouter initialEntries={["/experiments/10/compare/5/edit"]}>
        <Routes>
          <Route path="*" element={
            <>
              <LocationProbe />
              <AppShell />
            </>
          }/>
        </Routes>
      </MemoryRouter>,
    );
    const probe = await findByTestId("probe");
    expect(probe.textContent).toBe("/experiments/10/compare/5");
  });
});
```

- [ ] **Step 2: Run.**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/compareRouteFlatten.test.tsx
```

Expected: fail (no redirect yet).

- [ ] **Step 3: Implement the redirect in AppShell.**

In `AppShell.tsx`'s `<Routes>`, add a redirect Route BEFORE the existing
`/experiments/:eid/compare/:id/edit` and `/compare/all/:id/edit` routes
(or instead of them — they both then 301 to the bare URL):

```tsx
import { Navigate } from "react-router-dom";

// Add inside the Routes block, immediately after the Compare routes:
<Route
  path="/experiments/:eid/compare/:id/edit"
  element={<EditToBareRedirect />}
/>
<Route
  path="/compare/all/:id/edit"
  element={<EditToBareRedirect />}
/>
```

Above the component, define the helper:

```tsx
import { useLocation } from "react-router-dom";

function EditToBareRedirect(): JSX.Element {
  // Use react-router's useLocation, NOT window.location — under MemoryRouter
  // (which the test uses) JSDOM's window.location.pathname stays at "/" and
  // the redirect would land in the wrong place. useLocation reads from the
  // router's internal store and is correct in both browser and test envs.
  const loc = useLocation();
  const here = loc.pathname.replace(/\/edit\/?$/, "");
  return <Navigate to={here} replace />;
}
```

The existing routes that mount `ComparePageEdit` for `/edit` and `/new` URLs
remain — `ComparePageEdit` still mounts for `/new` (a draft surface that
isn't `/edit`); only `/edit` is redirected. `/new` will be folded into
`Compare.tsx` in Phase C.

- [ ] **Step 4: Run.**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/compareRouteFlatten.test.tsx
```

Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/AppShell.tsx packages/HimalayaUI/frontend/test/compareRouteFlatten.test.tsx
git commit -m "feat(routes): redirect Compare /edit → bare URL (Compare UX B-1)"
```

### Task B-2: Drop `edit?: boolean` from `comparePath`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/comparison/routes.ts`
- Modify: callers of `comparePath({ ..., edit: true })`
- Test: `packages/HimalayaUI/frontend/test/comparePath.test.ts` (extend or create)

- [ ] **Step 1: Locate callers.**

```bash
grep -rn "edit: true\|edit:true\|\.edit\b" packages/HimalayaUI/frontend/src/ \
  --include="*.ts" --include="*.tsx" | grep -v node_modules | grep -i compare
```

Expected: at least the `EditOrForkButton::onEdit` callsite inside
`ComparePage.tsx`. There may be more. List them.

- [ ] **Step 2: Failing test.**

In `test/comparePath.test.ts` (create if needed):

```ts
import { describe, it, expect } from "vitest";
import { comparePath } from "../src/lib/comparison/routes";

describe("comparePath drops edit option — Compare UX B-2", () => {
  it("returns a bare comparison URL (no /edit)", () => {
    expect(comparePath({ scope: "experiment", eid: 10, id: 5 }))
      .toBe("/experiments/10/compare/5");
    expect(comparePath({ scope: "all", id: 5 }))
      .toBe("/compare/all/5");
  });

  it("supports isNew (no /edit)", () => {
    expect(comparePath({ scope: "experiment", eid: 10, isNew: true }))
      .toBe("/experiments/10/compare/new");
  });
});
```

- [ ] **Step 3: Run.**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/comparePath.test.ts
```

Expected: passes today (the bare paths already work). The point of this
test is to lock the shape before we drop the option.

- [ ] **Step 4: Drop `edit` from `ComparePathOpts` and `comparePath`.**

In `routes.ts`, remove the `edit?: boolean | undefined` field from
`ComparePathOpts` and the branch that appends `/edit` from the body of
`comparePath`. (If your TS compiler proves the call site doesn't pass
`edit`, you can remove the parameter entirely.)

- [ ] **Step 5: Update callers.**

Replace any `comparePath({ ..., edit: true })` invocation with a navigation
to the bare path. `ComparePage.tsx::EditOrForkButton::onEdit` is the main
hit:

```tsx
const onEdit = (): void => {
  loadDraft(comparison, qc);
  navigate(comparePath({ scope, eid: experimentId, id: comparison.id }));
};
```

This is technically a behavioral change: previously clicking "Edit" took
the user to `/edit`. Now it stays at the bare URL but seeds a draft into
Zustand. The redirect from B-1 keeps any older URL working; this change
is the modern path.

- [ ] **Step 6: Run.**

```bash
cd packages/HimalayaUI/frontend && npm run build 2>&1 | tail -5
cd packages/HimalayaUI/frontend && npm test -- --run > /tmp/vitest-B2.out 2>&1
grep -E "Tests:|Failed|passed" /tmp/vitest-B2.out | tail -5
```

Expected: clean.

- [ ] **Step 7: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/lib/comparison/routes.ts packages/HimalayaUI/frontend/src/pages/ComparePage.tsx packages/HimalayaUI/frontend/test/comparePath.test.ts
git commit -m "refactor(routes): drop comparePath({edit}) option (Compare UX B-2)"
```

---

## Phase C — Title strip + status surface + toolbar + Save pill + merge into `Compare.tsx`

This is the largest phase. We build new presentational components
bottom-up, then merge `ComparePage` + `ComparePageEdit` into `Compare.tsx`
that uses them.

### Task C-1: `dragThreshold.ts` helper + tests

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/comparison/dragThreshold.ts`
- Test: `packages/HimalayaUI/frontend/test/lib/comparison/dragThreshold.test.ts`

- [ ] **Step 1: Failing test.**

```ts
import { describe, it, expect } from "vitest";
import { makeDragThresholdState, type DragOutcome } from "../../../src/lib/comparison/dragThreshold";

describe("dragThreshold — Compare UX C-1", () => {
  it("treats < 4px move as click", () => {
    const s = makeDragThresholdState({ thresholdPx: 4 });
    s.onPointerDown(10, 10);
    s.onPointerMove(12, 12);
    const outcome: DragOutcome = s.onPointerUp(12, 12);
    expect(outcome).toBe("click");
  });

  it("treats >= 4px move as drag", () => {
    const s = makeDragThresholdState({ thresholdPx: 4 });
    s.onPointerDown(10, 10);
    s.onPointerMove(15, 10);
    expect(s.onPointerUp(15, 10)).toBe("drag");
  });

  it("uses Manhattan distance for the threshold check (axis-agnostic)", () => {
    const s = makeDragThresholdState({ thresholdPx: 4 });
    s.onPointerDown(0, 0);
    // 3px in x AND 3px in y → euclidean ~4.24, manhattan 6 → drag either way
    s.onPointerMove(3, 3);
    expect(s.onPointerUp(3, 3)).toBe("drag");
  });
});
```

- [ ] **Step 2: Run.**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/lib/comparison/dragThreshold.test.ts
```

Expected: fail (module doesn't exist).

- [ ] **Step 3: Implement `dragThreshold.ts`.**

```ts
/**
 * Drag-vs-click disambiguation helper (spec §7.3 / Compare UX refinement).
 *
 * Usage:
 *   const state = makeDragThresholdState({ thresholdPx: 4 });
 *   onPointerDown:  state.onPointerDown(x, y)
 *   onPointerMove:  if (state.onPointerMove(x, y) === "drag-start") startDrag();
 *   onPointerUp:    const outcome = state.onPointerUp(x, y)
 *                    // outcome = "click" or "drag"
 *
 * Returns "click" when the pointer-up lands within `thresholdPx` Manhattan
 * distance of the down point AND no intervening move crossed the threshold.
 * Returns "drag" otherwise.
 */
export type DragOutcome = "click" | "drag";

export interface DragThresholdState {
  onPointerDown(x: number, y: number): void;
  onPointerMove(x: number, y: number): "drag-start" | "below-threshold";
  onPointerUp(x: number, y: number): DragOutcome;
  isDragging(): boolean;
}

export function makeDragThresholdState(
  opts: { thresholdPx: number },
): DragThresholdState {
  let downX: number | null = null;
  let downY: number | null = null;
  let dragging = false;

  return {
    onPointerDown(x, y) { downX = x; downY = y; dragging = false; },
    onPointerMove(x, y) {
      if (downX === null || downY === null) return "below-threshold";
      const dx = Math.abs(x - downX);
      const dy = Math.abs(y - downY);
      if (!dragging && (dx + dy) >= opts.thresholdPx) {
        dragging = true;
        return "drag-start";
      }
      return "below-threshold";
    },
    onPointerUp(x, y) {
      const wasDragging = dragging;
      if (!wasDragging && downX !== null && downY !== null) {
        const dx = Math.abs(x - downX);
        const dy = Math.abs(y - downY);
        if ((dx + dy) >= opts.thresholdPx) {
          dragging = true;
        }
      }
      const outcome: DragOutcome = dragging ? "drag" : "click";
      downX = null; downY = null; dragging = false;
      return outcome;
    },
    isDragging() { return dragging; },
  };
}
```

- [ ] **Step 4: Run.**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/lib/comparison/dragThreshold.test.ts
```

Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/lib/comparison/dragThreshold.ts packages/HimalayaUI/frontend/test/lib/comparison/dragThreshold.test.ts
git commit -m "feat(lib): drag threshold helper (Compare UX C-1)"
```

### Task C-2: `relativeTime.ts` helper + tests

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/comparison/relativeTime.ts`
- Test: `packages/HimalayaUI/frontend/test/lib/comparison/relativeTime.test.ts`

- [ ] **Step 1: Failing test.**

```ts
import { describe, it, expect } from "vitest";
import { relativeTime } from "../../../src/lib/comparison/relativeTime";

describe("relativeTime — Compare UX C-2", () => {
  const NOW = new Date("2026-05-14T12:00:00Z").getTime();
  const at = (s: string) => relativeTime(s, NOW);

  it("renders 'just now' for <60s", () => {
    expect(at("2026-05-14T11:59:30Z")).toBe("just now");
  });
  it("renders minutes", () => {
    expect(at("2026-05-14T11:55:00Z")).toBe("5m ago");
  });
  it("renders hours", () => {
    expect(at("2026-05-14T10:00:00Z")).toBe("2h ago");
  });
  it("renders days up to 30d", () => {
    expect(at("2026-05-10T12:00:00Z")).toBe("4d ago");
  });
  it("falls back to a date for >30d", () => {
    expect(at("2026-01-01T12:00:00Z")).toMatch(/\d{4}-\d{2}-\d{2}/);
  });
  it("returns null for null input", () => {
    expect(relativeTime(null, NOW)).toBeNull();
  });
});
```

- [ ] **Step 2: Run.**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/lib/comparison/relativeTime.test.ts
```

Expected: fail.

- [ ] **Step 3: Implement.**

```ts
/**
 * `relativeTime("2026-05-14T10:00:00Z")` → "2h ago".
 *
 * Granularity matches typical app UX: "just now" < 1min,
 * "Nm ago" < 1h, "Nh ago" < 24h, "Nd ago" < 30d, ISO date thereafter.
 */
export function relativeTime(iso: string | null, nowMs: number = Date.now()): string | null {
  if (iso === null) return null;
  const t = Date.parse(iso);
  if (Number.isNaN(t)) return null;
  const diffMs = nowMs - t;
  const s = Math.floor(diffMs / 1000);
  if (s < 60) return "just now";
  const m = Math.floor(s / 60);
  if (m < 60) return `${m}m ago`;
  const h = Math.floor(m / 60);
  if (h < 24) return `${h}h ago`;
  const d = Math.floor(h / 24);
  if (d <= 30) return `${d}d ago`;
  return new Date(t).toISOString().slice(0, 10);
}
```

- [ ] **Step 4: Run.**

Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/lib/comparison/relativeTime.ts packages/HimalayaUI/frontend/test/lib/comparison/relativeTime.test.ts
git commit -m "feat(lib): relativeTime helper (Compare UX C-2)"
```

### Task C-3: `InlineEditableText` component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/InlineEditableText.tsx`
- Test: `packages/HimalayaUI/frontend/test/InlineEditableText.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { InlineEditableText } from "../src/components/InlineEditableText";

describe("InlineEditableText — Compare UX C-3", () => {
  it("renders text by default, switches to input on click", () => {
    render(<InlineEditableText value="hello" onCommit={() => {}} />);
    const text = screen.getByText("hello");
    fireEvent.click(text);
    expect(screen.getByRole("textbox")).toHaveValue("hello");
  });

  it("commits on Enter", () => {
    const onCommit = vi.fn();
    render(<InlineEditableText value="hello" onCommit={onCommit} />);
    fireEvent.click(screen.getByText("hello"));
    const input = screen.getByRole("textbox");
    fireEvent.change(input, { target: { value: "world" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(onCommit).toHaveBeenCalledWith("world");
  });

  it("cancels on Esc", () => {
    const onCommit = vi.fn();
    render(<InlineEditableText value="hello" onCommit={onCommit} />);
    fireEvent.click(screen.getByText("hello"));
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "world" } });
    fireEvent.keyDown(screen.getByRole("textbox"), { key: "Escape" });
    expect(onCommit).not.toHaveBeenCalled();
    expect(screen.getByText("hello")).toBeInTheDocument();
  });

  it("commits on blur", () => {
    const onCommit = vi.fn();
    render(<InlineEditableText value="hello" onCommit={onCommit} />);
    fireEvent.click(screen.getByText("hello"));
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "world" } });
    fireEvent.blur(screen.getByRole("textbox"));
    expect(onCommit).toHaveBeenCalledWith("world");
  });

  it("carries data-interactable='edit'", () => {
    render(<InlineEditableText value="hello" onCommit={() => {}} />);
    expect(screen.getByText("hello").closest("[data-interactable]"))
      .toHaveAttribute("data-interactable", "edit");
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

```tsx
import { useEffect, useRef, useState } from "react";

interface Props {
  value: string;
  onCommit: (newValue: string) => void;
  /** Optional placeholder when value is empty */
  placeholder?: string;
  /** className for the rest-state element */
  className?: string;
  /** className for the input */
  inputClassName?: string;
  /** Test selector */
  testId?: string;
  /**
   * When true, renders as plain text — no click handler, no hover affordance,
   * no `data-interactable`. Used by review-mode `CompareTitleStrip` to render
   * the title as read-only without silently dropping clicks. Default `false`.
   */
  readOnly?: boolean;
}

/**
 * Click text → input. Enter commits, Esc cancels (restores prior value),
 * blur commits. Carries `data-interactable="edit"` for the
 * visual-language vocabulary (spec §4).
 *
 * When `readOnly` is true, behaves as a static `<span>` with no interactive
 * affordance. Caller can flip this based on `useCompareMode().kind ===
 * "viewing"` (or `"viewing-stale"`).
 */
export function InlineEditableText({
  value, onCommit, placeholder, className, inputClassName, testId,
  readOnly = false,
}: Props): JSX.Element {
  const [editing, setEditing] = useState(false);
  const [draft, setDraft] = useState(value);
  const inputRef = useRef<HTMLInputElement>(null);

  useEffect(() => { if (!editing) setDraft(value); }, [value, editing]);
  useEffect(() => {
    if (editing && inputRef.current) {
      inputRef.current.focus();
      inputRef.current.select();
    }
  }, [editing]);

  if (readOnly) {
    return (
      <span data-testid={testId} className={className}>
        {value === "" ? <span className="text-fg-dim">{placeholder ?? ""}</span> : value}
      </span>
    );
  }

  if (editing) {
    return (
      <span data-interactable="edit">
        <input
          ref={inputRef}
          type="text"
          data-testid={testId}
          value={draft}
          placeholder={placeholder}
          onChange={(e) => setDraft(e.target.value)}
          onKeyDown={(e) => {
            if (e.key === "Enter") { onCommit(draft); setEditing(false); }
            else if (e.key === "Escape") { setDraft(value); setEditing(false); }
          }}
          onBlur={() => { onCommit(draft); setEditing(false); }}
          className={inputClassName ?? "bg-transparent border-b border-border outline-none"}
        />
      </span>
    );
  }

  return (
    <span
      data-interactable="edit"
      data-testid={testId}
      onClick={() => setEditing(true)}
      role="button"
      tabIndex={0}
      onKeyDown={(e) => { if (e.key === "Enter") setEditing(true); }}
      className={`${className ?? ""} cursor-text hover:underline decoration-dotted underline-offset-4`}
    >
      {value === "" ? <span className="text-fg-dim">{placeholder ?? ""}</span> : value}
    </span>
  );
}
```

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/InlineEditableText.tsx packages/HimalayaUI/frontend/test/InlineEditableText.test.tsx
git commit -m "feat(ui): InlineEditableText component (Compare UX C-3)"
```

### Task C-4: Extend `ActiveDraft` with view-choice fields

**Build-order note.** This task was originally C-5; it ran AFTER the
`useCompareMode` hook (originally C-4). The plan reviewer caught that the
hook's test fixture references `viewGroupingMode` / `viewShowPeakTicks` /
`viewShowPeakLabels` on `ActiveDraft`, which TypeScript-strict mode
rejects until those fields exist. Swapping the order — draft extension
FIRST, then the hook — eliminates the compile-time miss. Task numbering
is updated accordingly.

**State-source resolution (load-bearing).** Today's `state.ts` owns
`groupingMode: GroupingMode` as Zustand client state with a named
action `setGroupingMode`. After C-4 the same conceptual value lives on
`ActiveDraft.viewGroupingMode` (client-state proxy for the saved
field) AND on the server record (`Comparison.view_grouping_mode` via
TanStack Query). Three sources of truth ⇒ the `<GroupingModeToggle />`
mutation lands on the wrong store half the time.

**Decision (apply before extending draft.ts):** delete the Zustand
`groupingMode` field and its `setGroupingMode` action. Introduce a
single canonical helper `effectiveGroupingMode(draft, comparison)`
that returns the active value for any render context:

```ts
// lib/comparison/effectiveGroupingMode.ts (new)
import type { ActiveDraft } from "./draft";
import type { Comparison } from "../../api";

export type GroupingMode = "bySample" | "byPhase" | "distinct";

export function effectiveGroupingMode(
  draft: ActiveDraft | null,
  comparison: Comparison | undefined,
): GroupingMode {
  // Draft owns the value during edit (incl. viewer escape hatch where
  // toggling without a draft creates one).
  if (draft?.viewGroupingMode !== undefined) return draft.viewGroupingMode;
  // Otherwise the server record is the source of truth.
  return (comparison?.view_grouping_mode as GroupingMode | null) ?? "bySample";
}
```

This matches the CLAUDE.md state-split rule: server state lives in
TanStack Query; client UI state that mirrors server values lives on
the draft, not in a separate Zustand field. The helper centralises
the precedence so consumers never spell out the fallback chain.

**Audit (preliminary work before Step 1).** Run a full grep to
enumerate every consumer of the Zustand `groupingMode`:

```bash
grep -rn "groupingMode\|setGroupingMode" packages/HimalayaUI/frontend/src \
  --include="*.ts" --include="*.tsx"
```

The current count is **7 files** with at least one usage each:

- `src/state.ts` — the field + action (delete both).
- `src/components/GroupingModeToggle.tsx` — read + write. Rewire:
  read via `effectiveGroupingMode(draft, comparison)`; write via the
  draft's `viewGroupingMode` setter (and auto-create a draft if none
  exists, per spec §6.4 viewer escape hatch).
- `src/pages/ComparePage.tsx` — 3 references including a prop pass
  through to `MultiTracePlot`. Replace with
  `effectiveGroupingMode(draft, comparison)` computed once per render
  near the top of the page body; thread the resulting value down.
- `src/pages/ComparePageEdit.tsx` — 2 references (incl. prop pass at
  the plot mount). Same rewire as above.
- `src/components/MultiTracePlot.tsx` — receives `groupingMode` as a
  PROP (not Zustand). Stays prop-driven; just confirm the prop is now
  sourced from the helper at the call site.
- `src/components/MemberTraceLayer.tsx` — same as MultiTracePlot:
  prop-driven, no rewire, just verify the upstream uses the helper.
- `src/lib/figure-export/adapters/multiTraceAdapter.ts` — the figure
  export path also takes `groupingMode`. Confirm its caller computes
  via the helper. The figure-export contract test should still pass
  unchanged (it asserts the rendered output, not the source of the
  value).

This pre-flight wiring lives in a new task (Step 0 below) so the
Zustand delete in Step 3 doesn't blow up seven call sites at once.

- [ ] **Step 0: Introduce `effectiveGroupingMode` helper + thread through both pages.**

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/comparison/effectiveGroupingMode.ts`
- Modify: `packages/HimalayaUI/frontend/src/pages/ComparePage.tsx` (3 refs)
- Modify: `packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx` (2 refs)
- Modify: `packages/HimalayaUI/frontend/src/components/GroupingModeToggle.tsx` (read + write)
- Modify: `packages/HimalayaUI/frontend/src/lib/figure-export/adapters/multiTraceAdapter.ts` (caller)
- Verify-only (no edit expected — prop-driven, just confirm the caller's source):
  `packages/HimalayaUI/frontend/src/components/MultiTracePlot.tsx`,
  `packages/HimalayaUI/frontend/src/components/MemberTraceLayer.tsx`.
- Test: `packages/HimalayaUI/frontend/test/lib/comparison/effectiveGroupingMode.test.ts`

The helper test:

```ts
import { describe, it, expect } from "vitest";
import { effectiveGroupingMode } from "../../../src/lib/comparison/effectiveGroupingMode";

describe("effectiveGroupingMode — Compare UX C-4 Step 0", () => {
  it("returns draft.viewGroupingMode when set", () => {
    const draft = { /* …minimal ActiveDraft… */ viewGroupingMode: "byPhase" } as any;
    const comparison = { view_grouping_mode: "distinct" } as any;
    expect(effectiveGroupingMode(draft, comparison)).toBe("byPhase");
  });

  it("falls back to comparison.view_grouping_mode when draft has undefined", () => {
    const draft = { /* …minimal ActiveDraft… */ viewGroupingMode: undefined } as any;
    const comparison = { view_grouping_mode: "distinct" } as any;
    expect(effectiveGroupingMode(draft, comparison)).toBe("distinct");
  });

  it("returns 'bySample' when both are null/undefined (default)", () => {
    expect(effectiveGroupingMode(null, undefined)).toBe("bySample");
    expect(effectiveGroupingMode(null, { view_grouping_mode: null } as any)).toBe("bySample");
  });
});
```

After this step lands, every consumer reads via the helper but the
Zustand field is still in place (no consumer reads it directly any
more — but `state.ts` still declares it). The next task (Step 1 below,
the original C-4 Step 1) then runs the failing-test loop for `ActiveDraft.viewGroupingMode`,
and Step 3+ deletes the Zustand `groupingMode`/`setGroupingMode`
declarations. Splitting into two commits keeps each commit
build-green.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/comparison/draft.ts`
- Modify: `packages/HimalayaUI/frontend/src/lib/comparison/draftFactories.ts`
- Test: `packages/HimalayaUI/frontend/test/lib/comparison/draft.test.ts` (new or extend)

- [ ] **Step 1: Failing test.**

```ts
import { describe, it, expect } from "vitest";
import { emptyDraft } from "../../../src/lib/comparison/draft";
import { loadDraftFromComparison } from "../../../src/lib/comparison/draftFactories";
import type { Comparison } from "../../../src/api";

describe("ActiveDraft includes view choices — Compare UX C-4", () => {
  it("emptyDraft seeds view choices as undefined", () => {
    const d = emptyDraft();
    expect(d.viewGroupingMode).toBeUndefined();
    expect(d.viewShowPeakTicks).toBeUndefined();
    expect(d.viewShowPeakLabels).toBeUndefined();
  });

  it("loadDraftFromComparison projects existing view choices", () => {
    const c: Comparison = {
      id: 1, title: "x", description: null, content_hash: "h",
      created_by: null, created_at: null, updated_at: null,
      forked_from_id: null, forked_at_hash: null, forked_from_title: null,
      view_grouping_mode: "byPhase", view_show_peak_ticks: false, view_show_peak_labels: true,
      last_event_at: null, members: [],
    };
    const d = loadDraftFromComparison(c);
    expect(d.viewGroupingMode).toBe("byPhase");
    expect(d.viewShowPeakTicks).toBe(false);
    expect(d.viewShowPeakLabels).toBe(true);
  });
});
```

(Verify the actual signature of `loadDraftFromComparison`; pattern
follows `draftFactories.ts`.)

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Update `draft.ts`.**

```ts
export interface ActiveDraft {
  id: number | undefined;
  baseHash: string | undefined;
  title: string;
  description: string;
  members: DraftMember[];
  forkedFromId: number | undefined;
  forkedAtHash: string | undefined;
  /** "bySample" | "byPhase" | "distinct" | undefined (means: inherit/server-default) */
  viewGroupingMode: "bySample" | "byPhase" | "distinct" | undefined;
  viewShowPeakTicks: boolean | undefined;
  viewShowPeakLabels: boolean | undefined;
}
```

Update `emptyDraft`:

```ts
export function emptyDraft(): ActiveDraft {
  return {
    id: undefined, baseHash: undefined,
    title: "", description: "", members: [],
    forkedFromId: undefined, forkedAtHash: undefined,
    viewGroupingMode: undefined,
    viewShowPeakTicks: undefined,
    viewShowPeakLabels: undefined,
  };
}
```

- [ ] **Step 4: Update `draftFactories.ts`.**

In `loadDraftFromComparison` (and any sibling factories that build an
ActiveDraft from a Comparison), copy the three view fields with `??` to
coalesce server `null` to `undefined`:

```ts
viewGroupingMode:    c.view_grouping_mode    as ActiveDraft["viewGroupingMode"] ?? undefined,
viewShowPeakTicks:   c.view_show_peak_ticks  ?? undefined,
viewShowPeakLabels:  c.view_show_peak_labels ?? undefined,
```

- [ ] **Step 5: Update `saveComparison` callsites to forward view fields.**

Wherever `useSaveComparison.mutate({...})` is built today (mainly
`ComparePageEdit.tsx::handleSave` and the `ConflictModal` submit path):

```ts
const payload: SaveComparisonInput = {
  id: draft.id, title: draft.title, description: draft.description,
  members: builtMembers,
  expected_content_hash: draft.baseHash,
  forked_from_id: draft.forkedFromId,
  forked_at_hash: draft.forkedAtHash,
  view_grouping_mode:    draft.viewGroupingMode    ?? null,
  view_show_peak_ticks:  draft.viewShowPeakTicks   ?? null,
  view_show_peak_labels: draft.viewShowPeakLabels  ?? null,
};
```

**Do NOT include `username`, `clientId`, or `clientOpId` in the
caller-built payload.** `useQueueMutation.mutationFn` mints
`client_op_id` per `mutate()` call and injects username + clientId
(see `useQueueMutation.ts:126`). Capturing `newClientOpId()` at payload
construction time shares one idempotency key across all retries of all
mutations — the canonical mutation-queue gotcha from CLAUDE.md. The
existing `ComparePageEdit.tsx::handleSave` call site already omits
these fields; preserve that.

- [ ] **Step 6: Run.**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/lib/comparison/draft.test.ts
cd packages/HimalayaUI/frontend && npm run build 2>&1 | tail -10
```

Expected: PASS + build clean.

- [ ] **Step 7: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/lib/comparison/draft.ts packages/HimalayaUI/frontend/src/lib/comparison/draftFactories.ts packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx packages/HimalayaUI/frontend/src/components/ConflictModal.tsx packages/HimalayaUI/frontend/test/lib/comparison/draft.test.ts
git commit -m "feat(draft): ActiveDraft carries view choices (Compare UX C-4)"
```

### Task C-5: `useCompareMode` hook

**Files:**
- Create: `packages/HimalayaUI/frontend/src/hooks/useCompareMode.ts`
- Test: `packages/HimalayaUI/frontend/test/useCompareMode.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
import { describe, it, expect } from "vitest";
import { renderHook } from "@testing-library/react";
import { useCompareMode } from "../src/hooks/useCompareMode";
import { useAppState } from "../src/state";
import type { Comparison } from "../src/api";

const baseComparison: Comparison = {
  id: 1, title: "x", description: null, content_hash: "h",
  created_by: 5, created_at: null, updated_at: null,
  forked_from_id: null, forked_at_hash: null, forked_from_title: null,
  view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
  last_event_at: null, members: [],
};

describe("useCompareMode — Compare UX C-5", () => {
  it("returns 'viewing' when no draft + no comparison", () => {
    useAppState.setState({ activeDraft: null });
    const { result } = renderHook(() => useCompareMode({
      comparison: undefined, currentUserId: 5,
    }));
    expect(result.current.kind).toBe("viewing");
  });

  it("returns 'viewing' when comparison loaded and no draft", () => {
    useAppState.setState({ activeDraft: null });
    const { result } = renderHook(() => useCompareMode({
      comparison: baseComparison, currentUserId: 5,
    }));
    expect(result.current.kind).toBe("viewing");
  });

  it("returns 'editing-mine' when draft id matches and user is author", () => {
    useAppState.setState({
      activeDraft: { id: 1, baseHash: "h", title: "x", description: "",
        members: [], forkedFromId: undefined, forkedAtHash: undefined,
        viewGroupingMode: undefined, viewShowPeakTicks: undefined, viewShowPeakLabels: undefined },
    });
    const { result } = renderHook(() => useCompareMode({
      comparison: baseComparison, currentUserId: 5,
    }));
    expect(result.current.kind).toBe("editing-mine");
  });

  it("returns 'editing-as-fork-of' when draft id matches but user is NOT author", () => {
    useAppState.setState({
      activeDraft: { id: 1, baseHash: "h", title: "x", description: "",
        members: [], forkedFromId: undefined, forkedAtHash: undefined,
        viewGroupingMode: undefined, viewShowPeakTicks: undefined, viewShowPeakLabels: undefined },
    });
    const { result } = renderHook(() => useCompareMode({
      comparison: baseComparison, currentUserId: 99,
    }));
    expect(result.current.kind).toBe("editing-as-fork-of");
  });

  it("returns 'creating-blank' when draft id is undefined and no fork lineage", () => {
    useAppState.setState({
      activeDraft: { id: undefined, baseHash: undefined, title: "n", description: "",
        members: [], forkedFromId: undefined, forkedAtHash: undefined,
        viewGroupingMode: undefined, viewShowPeakTicks: undefined, viewShowPeakLabels: undefined },
    });
    const { result } = renderHook(() => useCompareMode({
      comparison: undefined, currentUserId: 5,
    }));
    expect(result.current.kind).toBe("creating-blank");
  });

  it("returns 'creating-from-fork' when draft id is undefined and forkedFromId is set", () => {
    useAppState.setState({
      activeDraft: { id: undefined, baseHash: undefined, title: "Copy of x", description: "",
        members: [], forkedFromId: 42, forkedAtHash: "h_parent",
        viewGroupingMode: undefined, viewShowPeakTicks: undefined, viewShowPeakLabels: undefined },
    });
    const { result } = renderHook(() => useCompareMode({
      comparison: undefined, currentUserId: 5,
    }));
    expect(result.current.kind).toBe("creating-from-fork");
    if (result.current.kind === "creating-from-fork") {
      expect(result.current.parentId).toBe(42);
    }
  });

  it("returns 'viewing-stale' when no draft and staleAgainstHash is set", () => {
    useAppState.setState({ activeDraft: null });
    const { result } = renderHook(() => useCompareMode({
      comparison: baseComparison, currentUserId: 5,
      staleAgainstHash: "h_prev",
    }));
    expect(result.current.kind).toBe("viewing-stale");
  });
});
```

Currently `currentUserId` is exposed via a `useCurrentUserId()` hook, not
a direct Zustand field — so the hook takes it as a parameter to stay
pure and testable. `Compare.tsx` passes `useCurrentUserId()` in.

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

```ts
import { useMemo } from "react";
import { useAppState } from "../state";
import type { Comparison } from "../api";

export type CompareMode =
  | { kind: "viewing" }
  | { kind: "viewing-stale"; previousHash: string }
  | { kind: "editing-mine"; draftId: number }
  | { kind: "editing-as-fork-of"; parentId: number }
  | { kind: "creating-blank" }
  | { kind: "creating-from-fork"; parentId: number };

export function useCompareMode(
  opts: {
    comparison: Comparison | undefined;
    currentUserId: number | undefined;
    /** Set when a foreign SSE event drifted content_hash since last read. */
    staleAgainstHash?: string | undefined;
  },
): CompareMode {
  const draft = useAppState((s) => s.activeDraft);

  return useMemo<CompareMode>(() => {
    if (draft === null) {
      if (opts.staleAgainstHash !== undefined) {
        return { kind: "viewing-stale", previousHash: opts.staleAgainstHash };
      }
      return { kind: "viewing" };
    }
    if (draft.id === undefined) {
      if (draft.forkedFromId !== undefined) {
        return { kind: "creating-from-fork", parentId: draft.forkedFromId };
      }
      return { kind: "creating-blank" };
    }
    // Draft is tied to an existing comparison id.
    const author = opts.comparison?.created_by ?? null;
    const isAuthor = author !== null
      && opts.currentUserId !== undefined
      && opts.currentUserId === author;
    if (isAuthor) return { kind: "editing-mine", draftId: draft.id };
    return { kind: "editing-as-fork-of", parentId: draft.id };
  }, [draft, opts.currentUserId, opts.comparison?.created_by, opts.staleAgainstHash]);
}
```

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/hooks/useCompareMode.ts packages/HimalayaUI/frontend/test/useCompareMode.test.tsx
git commit -m "feat(hooks): useCompareMode tagged-union (Compare UX C-5)"
```

### Task C-6: `useCompareDraftDirty` hook

**Files:**
- Create: `packages/HimalayaUI/frontend/src/hooks/useCompareDraftDirty.ts`
- Test: `packages/HimalayaUI/frontend/test/useCompareDraftDirty.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
import { describe, it, expect } from "vitest";
import { renderHook } from "@testing-library/react";
import { useCompareDraftDirty } from "../src/hooks/useCompareDraftDirty";
import { useAppState } from "../src/state";

describe("useCompareDraftDirty — Compare UX C-6", () => {
  it("returns false when there is no draft", () => {
    useAppState.setState({ activeDraft: null });
    const { result } = renderHook(() => useCompareDraftDirty());
    expect(result.current).toBe(false);
  });
  it("returns true when there is any draft", () => {
    useAppState.setState({
      activeDraft: { id: undefined, baseHash: undefined,
        title: "x", description: "", members: [],
        forkedFromId: undefined, forkedAtHash: undefined,
        viewGroupingMode: undefined, viewShowPeakTicks: undefined, viewShowPeakLabels: undefined },
    });
    const { result } = renderHook(() => useCompareDraftDirty());
    expect(result.current).toBe(true);
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

```ts
import { useAppState } from "../state";

/**
 * Returns `true` whenever there is any active draft. The simplification —
 * "presence of a draft = dirty" — is intentional. Drafts only exist as
 * the result of an actual edit gesture (per spec §3, baseHash capture
 * happens on first edit), so this is sufficient.
 */
export function useCompareDraftDirty(): boolean {
  return useAppState((s) => s.activeDraft !== null);
}
```

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/hooks/useCompareDraftDirty.ts packages/HimalayaUI/frontend/test/useCompareDraftDirty.test.tsx
git commit -m "feat(hooks): useCompareDraftDirty (Compare UX C-6)"
```

### Task C-7: `SavePill` component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/SavePill.tsx`
- Test: `packages/HimalayaUI/frontend/test/SavePill.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { SavePill } from "../src/components/SavePill";

describe("SavePill — Compare UX C-7", () => {
  it("hides when not dirty", () => {
    const { queryByTestId } = render(
      <SavePill dirty={false} mode={{ kind: "viewing" }} onSave={() => {}} isSaving={false}/>);
    expect(queryByTestId("save-pill")).toBeNull();
  });

  it("shows 'Save changes' for editing-mine", () => {
    render(<SavePill dirty={true} mode={{ kind: "editing-mine", draftId: 1 }} onSave={() => {}} isSaving={false}/>);
    expect(screen.getByTestId("save-pill")).toHaveTextContent("Save changes");
  });

  it("shows 'Save as fork…' for editing-as-fork-of", () => {
    render(<SavePill dirty={true} mode={{ kind: "editing-as-fork-of", parentId: 1 }} onSave={() => {}} isSaving={false}/>);
    expect(screen.getByTestId("save-pill")).toHaveTextContent(/Save as fork/);
  });

  it("shows 'Save' for creating-blank", () => {
    render(<SavePill dirty={true} mode={{ kind: "creating-blank" }} onSave={() => {}} isSaving={false}/>);
    expect(screen.getByTestId("save-pill")).toHaveTextContent("Save");
  });

  it("shows 'Save fork' for creating-from-fork (post-morph)", () => {
    render(<SavePill dirty={true} mode={{ kind: "creating-from-fork", parentId: 1 }} onSave={() => {}} isSaving={false}/>);
    expect(screen.getByTestId("save-pill")).toHaveTextContent(/Save fork/);
  });

  it("calls onSave when clicked", () => {
    const onSave = vi.fn();
    render(<SavePill dirty={true} mode={{ kind: "editing-mine", draftId: 1 }} onSave={onSave} isSaving={false}/>);
    fireEvent.click(screen.getByTestId("save-pill"));
    expect(onSave).toHaveBeenCalled();
  });

  it("disables and shows 'Saving…' while isSaving", () => {
    render(<SavePill dirty={true} mode={{ kind: "editing-mine", draftId: 1 }} onSave={() => {}} isSaving={true}/>);
    expect(screen.getByTestId("save-pill")).toBeDisabled();
    expect(screen.getByTestId("save-pill")).toHaveTextContent("Saving…");
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

```tsx
import type { CompareMode } from "../hooks/useCompareMode";

interface Props {
  dirty: boolean;
  mode: CompareMode;
  onSave: () => void;
  isSaving: boolean;
}

const COPY: Record<CompareMode["kind"], string> = {
  "viewing":             "",  // hidden anyway
  "viewing-stale":       "",  // hidden anyway
  "editing-mine":        "Save changes",
  "editing-as-fork-of":  "Save as fork…",
  "creating-blank":      "Save",
  "creating-from-fork":  "Save fork",  // post-morph; user already confirmed the title
};

export function SavePill({ dirty, mode, onSave, isSaving }: Props): JSX.Element | null {
  if (!dirty) return null;
  return (
    <button
      type="button"
      data-testid="save-pill"
      data-interactable="button"
      data-mode={mode.kind}
      disabled={isSaving}
      onClick={onSave}
      title="Save (Cmd+Enter)"
      className="px-3 py-1 rounded bg-accent text-bg disabled:opacity-50 text-sm font-medium"
    >
      {isSaving ? "Saving…" : COPY[mode.kind]}
    </button>
  );
}
```

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/SavePill.tsx packages/HimalayaUI/frontend/test/SavePill.test.tsx
git commit -m "feat(ui): SavePill with mode-aware copy (Compare UX C-7)"
```

### Task C-8: `CompareTitleStrip` component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/CompareTitleStrip.tsx`
- Test: `packages/HimalayaUI/frontend/test/CompareTitleStrip.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { CompareTitleStrip } from "../src/components/CompareTitleStrip";

const wrap = (ui: JSX.Element) => <MemoryRouter>{ui}</MemoryRouter>;

describe("CompareTitleStrip — Compare UX C-8", () => {
  it("renders title and meta", () => {
    render(wrap(<CompareTitleStrip
      title="Cubic vs Hex"
      description={null}
      memberCount={4}
      authorUsername="alice"
      isCurrentUserAuthor={false}
      lastEventAt={null}
      forkedFromTitle={null}
      forkedFromHref={null}
      onTitleChange={() => {}}
      onDescChange={() => {}}
      now={Date.parse("2026-05-14T12:00:00Z")}
    />));
    expect(screen.getByText("Cubic vs Hex")).toBeInTheDocument();
    expect(screen.getByText(/by alice/)).toBeInTheDocument();
    expect(screen.getByText(/4 traces/)).toBeInTheDocument();
  });

  it("renders 'by you' when current user is the author", () => {
    render(wrap(<CompareTitleStrip
      title="t" description={null}
      memberCount={1}
      authorUsername="alice"
      isCurrentUserAuthor={true}
      lastEventAt={null} forkedFromTitle={null} forkedFromHref={null}
      onTitleChange={() => {}} onDescChange={() => {}}
      now={Date.now()}
    />));
    expect(screen.getByText(/by you/)).toBeInTheDocument();
  });

  it("renders forked-from link when set", () => {
    render(wrap(<CompareTitleStrip
      title="t" description={null}
      memberCount={1} authorUsername={null} isCurrentUserAuthor={false}
      lastEventAt={null}
      forkedFromTitle="Parent"
      forkedFromHref="/compare/all/9"
      onTitleChange={() => {}} onDescChange={() => {}}
      now={Date.now()}
    />));
    const link = screen.getByText("Parent");
    expect(link.closest("a")).toHaveAttribute("href", "/compare/all/9");
  });

  it("commits a title edit", () => {
    const onTitleChange = vi.fn();
    render(wrap(<CompareTitleStrip
      title="old" description={null}
      memberCount={1} authorUsername={null} isCurrentUserAuthor={false}
      lastEventAt={null} forkedFromTitle={null} forkedFromHref={null}
      onTitleChange={onTitleChange}
      onDescChange={() => {}}
      now={Date.now()}
    />));
    fireEvent.click(screen.getByText("old"));
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "new" } });
    fireEvent.keyDown(screen.getByRole("textbox"), { key: "Enter" });
    expect(onTitleChange).toHaveBeenCalledWith("new");
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

```tsx
import { Link } from "react-router-dom";
import { InlineEditableText } from "./InlineEditableText";
import { relativeTime } from "../lib/comparison/relativeTime";

interface Props {
  title: string;
  description: string | null;
  memberCount: number;
  authorUsername: string | null;
  isCurrentUserAuthor: boolean;
  lastEventAt: string | null;
  forkedFromTitle: string | null;
  forkedFromHref: string | null;
  onTitleChange: (s: string) => void;
  onDescChange: (s: string) => void;
  /**
   * When true, title and description render as plain text (no click-to-edit
   * affordance). Used in review mode so an idle click doesn't silently
   * discard the user's typing. Default `false`.
   */
  readOnly?: boolean;
  now?: number; // injected for tests
}

export function CompareTitleStrip(p: Props): JSX.Element {
  const byline = p.isCurrentUserAuthor
    ? "by you"
    : p.authorUsername !== null
      ? `by ${p.authorUsername}`
      : "by —";
  const rel = relativeTime(p.lastEventAt, p.now);
  return (
    <div data-testid="compare-title-strip" className="flex flex-col gap-1">
      <h2 className="text-lg font-semibold text-fg">
        <InlineEditableText
          value={p.title}
          onCommit={p.onTitleChange}
          placeholder="Untitled comparison"
          testId="compare-title"
          readOnly={p.readOnly ?? false}
        />
      </h2>
      <div className="text-sm text-fg-dim flex items-center gap-1 flex-wrap">
        <span>{byline}</span>
        <span aria-hidden="true">·</span>
        <span>{rel === null ? "—" : `edited ${rel}`}</span>
        <span aria-hidden="true">·</span>
        <span>{p.memberCount} traces</span>
        {p.forkedFromHref !== null && p.forkedFromTitle !== null && (
          <>
            <span aria-hidden="true">·</span>
            <span>
              forked from{" "}
              <Link to={p.forkedFromHref} className="text-accent hover:underline">
                {p.forkedFromTitle}
              </Link>
            </span>
          </>
        )}
      </div>
      {p.description !== null && (
        <div className="text-sm text-fg-dim">
          <InlineEditableText
            value={p.description ?? ""}
            onCommit={p.onDescChange}
            placeholder="Add a description…"
            testId="compare-description"
            readOnly={p.readOnly ?? false}
          />
        </div>
      )}
    </div>
  );
}
```

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/CompareTitleStrip.tsx packages/HimalayaUI/frontend/test/CompareTitleStrip.test.tsx
git commit -m "feat(ui): CompareTitleStrip (Compare UX C-8)"
```

### Task C-9: `CompareStatusSurface` component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/CompareStatusSurface.tsx`
- Test: `packages/HimalayaUI/frontend/test/CompareStatusSurface.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { CompareStatusSurface } from "../src/components/CompareStatusSurface";

describe("CompareStatusSurface — Compare UX C-9", () => {
  it("renders nothing when there are no banners", () => {
    const { container } = render(
      <CompareStatusSurface
        needsReview={null} serverUpdate={null} savedAt={null}
      />,
    );
    expect(container.querySelector("[data-testid='compare-status-surface']")).toBeNull();
  });

  it("renders a needs-review banner with action", () => {
    const onReanalyze = vi.fn();
    render(<CompareStatusSurface
      needsReview={{ onReanalyze }}
      serverUpdate={null}
      savedAt={null}
    />);
    expect(screen.getByText(/Needs review/)).toBeInTheDocument();
    fireEvent.click(screen.getByTestId("compare-status-resnapshot"));
    expect(onReanalyze).toHaveBeenCalled();
  });

  it("renders a server-update banner with acknowledge action", () => {
    const onAcknowledge = vi.fn();
    render(<CompareStatusSurface
      needsReview={null}
      serverUpdate={{ previousHash: "h_prev", onAcknowledge }}
      savedAt={null}
    />);
    expect(screen.getByText(/updated since you last viewed/i)).toBeInTheDocument();
    fireEvent.click(screen.getByTestId("compare-status-acknowledge"));
    expect(onAcknowledge).toHaveBeenCalled();
  });

  it("renders 'Saved' after a recent save", () => {
    render(<CompareStatusSurface needsReview={null} serverUpdate={null} savedAt={Date.now()}/>);
    expect(screen.getByText("Saved")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

```tsx
interface NeedsReview { onReanalyze: () => void; }
// `previousHash` is the hash the user last acknowledged; the banner
// surfaces while it differs from the current `content_hash`.
// `onAcknowledge` rebases the acknowledged hash to the current value
// (see `useStaleAgainstHash` in C-15 Step 4b).
interface ServerUpdate { previousHash: string; onAcknowledge: () => void; }

interface Props {
  needsReview: NeedsReview | null;
  serverUpdate: ServerUpdate | null;
  savedAt: number | null;
}

export function CompareStatusSurface(p: Props): JSX.Element | null {
  if (!p.needsReview && !p.serverUpdate && p.savedAt === null) return null;
  return (
    <div data-testid="compare-status-surface" className="flex flex-col gap-1">
      {p.needsReview && (
        <div className="px-3 py-2 rounded border border-warning bg-warning/10 text-sm text-fg flex items-center gap-2">
          <span aria-hidden="true">⚠</span>
          <span>
            Needs review — analysis changed since last submit.
          </span>
          <button
            type="button"
            data-testid="compare-status-resnapshot"
            data-interactable="button"
            onClick={p.needsReview.onReanalyze}
            className="ml-auto px-2 py-0.5 rounded border border-border hover:bg-bg-hover text-xs"
          >
            Re-snapshot
          </button>
        </div>
      )}
      {p.serverUpdate && (
        <div
          data-testid="compare-status-server-update"
          className="px-3 py-2 rounded border border-border bg-bg-subtle text-sm text-fg flex items-center gap-2"
        >
          <span>Server-side updated since you last viewed — save may conflict.</span>
          <button
            type="button"
            data-testid="compare-status-acknowledge"
            data-interactable="button"
            onClick={p.serverUpdate.onAcknowledge}
            className="ml-auto px-2 py-0.5 rounded border border-border hover:bg-bg-hover text-xs"
          >
            Acknowledge
          </button>
        </div>
      )}
      {showSaved && (
        <div className="px-3 py-2 rounded border border-success bg-success/10 text-sm text-fg">
          Saved
        </div>
      )}
    </div>
  );
}
```

**`showSaved` lifecycle.** `Date.now() - p.savedAt < 2000` in the JSX
makes the visibility a function of wall-clock time but the component
has no re-render trigger after the 2s window — the pill would stick
forever once shown. Compute via `useEffect`:

```tsx
const [showSaved, setShowSaved] = useState(false);
useEffect(() => {
  if (p.savedAt === null) {
    setShowSaved(false);
    return;
  }
  const elapsed = Date.now() - p.savedAt;
  if (elapsed >= 2000) {
    setShowSaved(false);
    return;
  }
  setShowSaved(true);
  const handle = setTimeout(() => setShowSaved(false), 2000 - elapsed);
  return () => clearTimeout(handle);
}, [p.savedAt]);
```

Also extend the test to assert the pill disappears after the window
(use `vi.useFakeTimers()` + `vi.advanceTimersByTime(2000)`):

```ts
it("hides the Saved pill after 2s", async () => {
  vi.useFakeTimers();
  const savedAt = Date.now();
  render(<CompareStatusSurface needsReview={null} serverUpdate={null} savedAt={savedAt} />);
  expect(screen.getByText("Saved")).toBeInTheDocument();
  vi.advanceTimersByTime(2100);
  expect(screen.queryByText("Saved")).toBeNull();
  vi.useRealTimers();
});
```

Without the effect-driven lifecycle, the existing "renders 'Saved'
after a recent save" assertion passes deceptively — the pill never
goes away because React doesn't re-render on wall-clock changes.

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/CompareStatusSurface.tsx packages/HimalayaUI/frontend/test/CompareStatusSurface.test.tsx
git commit -m "feat(ui): CompareStatusSurface banner stack (Compare UX C-9)"
```

### Task C-10: `CompareToolbar` component (Group / Annot / ⋯ more / Export + SavePill)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/CompareToolbar.tsx`
- Test: `packages/HimalayaUI/frontend/test/CompareToolbar.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { CompareToolbar } from "../src/components/CompareToolbar";

const wrap = (ui: JSX.Element) => <MemoryRouter>{ui}</MemoryRouter>;

describe("CompareToolbar — Compare UX C-10", () => {
  it("mounts the children controls", () => {
    render(wrap(<CompareToolbar
      groupingControl={<button data-testid="g">G</button>}
      annotationControl={<button data-testid="a">A</button>}
      forksCount={2}
      onCopyLink={() => {}}
      onDelete={() => {}}
      onDiscardChanges={null}
      onFork={() => {}}
      exportControl={<button data-testid="x">X</button>}
      saveControl={<button data-testid="s">S</button>}
    />));
    expect(screen.getByTestId("g")).toBeInTheDocument();
    expect(screen.getByTestId("a")).toBeInTheDocument();
    expect(screen.getByTestId("x")).toBeInTheDocument();
    expect(screen.getByTestId("s")).toBeInTheDocument();
  });

  it("opens the more menu and triggers actions", () => {
    const onCopyLink = vi.fn();
    render(wrap(<CompareToolbar
      groupingControl={null} annotationControl={null}
      forksCount={0} onCopyLink={onCopyLink} onDelete={() => {}}
      onDiscardChanges={null} onFork={() => {}}
      exportControl={null} saveControl={null}
    />));
    fireEvent.click(screen.getByTestId("compare-toolbar-more"));
    fireEvent.click(screen.getByText("Copy link"));
    expect(onCopyLink).toHaveBeenCalled();
  });

  it("includes 'Discard changes' only when dirty", () => {
    const { queryByText, rerender } = render(wrap(<CompareToolbar
      groupingControl={null} annotationControl={null}
      forksCount={0} onCopyLink={() => {}} onDelete={() => {}}
      onDiscardChanges={null} onFork={() => {}}
      exportControl={null} saveControl={null}
    />));
    fireEvent.click(screen.getByTestId("compare-toolbar-more"));
    expect(queryByText("Discard changes")).toBeNull();

    rerender(wrap(<CompareToolbar
      groupingControl={null} annotationControl={null}
      forksCount={0} onCopyLink={() => {}} onDelete={() => {}}
      onDiscardChanges={() => {}} onFork={() => {}}
      exportControl={null} saveControl={null}
    />));
    fireEvent.click(screen.getByTestId("compare-toolbar-more"));
    expect(screen.getByText("Discard changes")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

```tsx
import { useEffect, useRef, useState } from "react";
import { useFocusTrap } from "../hooks/useFocusTrap";

interface Props {
  groupingControl: React.ReactNode;
  annotationControl: React.ReactNode;
  forksCount: number;
  onCopyLink: () => void;
  onDelete: () => void;
  onDiscardChanges: (() => void) | null;
  onFork: () => void;
  exportControl: React.ReactNode;
  saveControl: React.ReactNode;
}

export function CompareToolbar(p: Props): JSX.Element {
  const [open, setOpen] = useState(false);
  const menuRef = useRef<HTMLDivElement>(null);
  // Trap Tab/Shift+Tab inside the open menu — canonical pattern from
  // NavModal / OnboardingFlow. Restores focus to the trigger on close.
  useFocusTrap(menuRef, open);

  useEffect(() => {
    if (!open) return;
    const onDoc = (e: MouseEvent) => {
      if (!menuRef.current?.contains(e.target as Node)) setOpen(false);
    };
    const onEsc = (e: KeyboardEvent) => { if (e.key === "Escape") setOpen(false); };
    document.addEventListener("mousedown", onDoc);
    document.addEventListener("keydown", onEsc);
    return () => {
      document.removeEventListener("mousedown", onDoc);
      document.removeEventListener("keydown", onEsc);
    };
  }, [open]);

  return (
    <div data-testid="compare-toolbar" className="flex items-center gap-2 flex-wrap">
      {p.groupingControl}
      {p.annotationControl}
      <div className="relative" ref={menuRef}>
        <button
          type="button"
          data-testid="compare-toolbar-more"
          data-interactable="button"
          aria-expanded={open}
          aria-haspopup="menu"
          onClick={() => setOpen((v) => !v)}
          className="px-2 py-0.5 rounded text-xs text-fg-dim hover:text-fg hover:bg-bg-hover border border-transparent hover:border-border"
        >
          ⋯ more
        </button>
        {open && (
          <div role="menu" className="absolute z-50 mt-1 right-0 min-w-[200px] card border border-border bg-bg-elevated shadow-lg p-1">
            <button
              type="button"
              data-testid="compare-toolbar-forks"
              data-interactable="button"
              onClick={() => { setOpen(false); }}
              className="w-full text-left px-2 py-1 hover:bg-bg-hover text-sm"
            >
              Forks ({p.forksCount})
            </button>
            <button type="button" onClick={() => { p.onCopyLink(); setOpen(false); }} className="w-full text-left px-2 py-1 hover:bg-bg-hover text-sm">Copy link</button>
            <button type="button" onClick={() => { p.onFork(); setOpen(false); }} className="w-full text-left px-2 py-1 hover:bg-bg-hover text-sm">Fork</button>
            <button type="button" onClick={() => { p.onDelete(); setOpen(false); }} className="w-full text-left px-2 py-1 hover:bg-bg-hover text-sm text-danger">Delete</button>
            {p.onDiscardChanges !== null && (
              <button type="button" onClick={() => { p.onDiscardChanges!(); setOpen(false); }} className="w-full text-left px-2 py-1 hover:bg-bg-hover text-sm">Discard changes</button>
            )}
          </div>
        )}
      </div>
      <span className="ml-auto inline-flex items-center gap-2">
        {p.exportControl}
        {p.saveControl}
      </span>
    </div>
  );
}
```

(Note: this currently does NOT mount a Forks dropdown — that's wired in
Task C-14 where the toolbar receives a `ForksList` prop, replacing the
"Forks (N)" button that closes the menu without doing anything. Left as a
stub here for visual completeness; the click handler is a no-op.)

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/CompareToolbar.tsx packages/HimalayaUI/frontend/test/CompareToolbar.test.tsx
git commit -m "feat(ui): CompareToolbar shell with ⋯-more menu (Compare UX C-10)"
```

### Task C-11: Create the `Compare.tsx` skeleton

**Files:**
- Create: `packages/HimalayaUI/frontend/src/pages/Compare.tsx`
- Modify: `packages/HimalayaUI/frontend/src/components/AppShell.tsx`
- Test: `packages/HimalayaUI/frontend/test/Compare.routing.test.tsx`

This task introduces the merged page but, **for safety, keeps it
behaviourally identical** to the today's split pages by delegating to
`ComparePageEdit` when a draft is active and `ComparePage` when not.
Subsequent tasks (C-12…C-15) progressively replace the delegated bodies
with the new title strip / status surface / toolbar / save pill, finally
letting us delete `ComparePage.tsx` and `ComparePageEdit.tsx`.

- [ ] **Step 1: Failing routing test.**

```tsx
import { describe, it, expect } from "vitest";
import { render } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { AppShell } from "../src/components/AppShell";

describe("Compare unified routing — Compare UX C-11", () => {
  it("mounts Compare.tsx for /compare/all/:id", async () => {
    const { findByTestId } = render(
      <MemoryRouter initialEntries={["/compare/all/5"]}>
        <AppShell />
      </MemoryRouter>,
    );
    // Compare.tsx exports a testid that ComparePage/ComparePageEdit do not.
    const el = await findByTestId("compare-page");
    expect(el).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement `Compare.tsx`.**

```tsx
import { useParams, useLocation } from "react-router-dom";
import { useAppState } from "../state";
import { ComparePage } from "./ComparePage";
import { ComparePageEdit } from "./ComparePageEdit";

/**
 * Compare.tsx — unified single-mode shell for both /compare/:id and
 * /compare/new. Subsequent Compare UX tasks (C-12 onward) progressively
 * replace the delegated bodies with the new title strip / toolbar / save
 * pill, and ultimately delete ComparePage.tsx + ComparePageEdit.tsx.
 *
 * Until then this is a passthrough that switches based on draft presence
 * (matches the today's URL-driven split exactly, just at a different
 * decision boundary). The mode flip the user is suffering today (URL
 * /edit vs bare) is gone immediately — adding a trace no longer requires
 * a route change.
 */
export function Compare(): JSX.Element {
  const params = useParams<{ eid?: string; id?: string }>();
  const location = useLocation();
  const isNew = location.pathname.endsWith("/new");
  const hasDraft = useAppState((s) => s.activeDraft !== null);

  // The decision logic is: if URL says /new OR there's an active draft,
  // mount the edit body; otherwise mount the review body.
  const Body = (isNew || hasDraft) ? ComparePageEdit : ComparePage;

  // `className="contents"` makes this wrapper participate as a
  // display:contents container — it preserves the testid handle without
  // taking a grid track. WorkspaceGrid mounts children directly into the
  // 3-column subgrid via `display: grid`; a regular <div> wrapper would
  // collapse all three tracks into the first cell.
  return (
    <div data-testid="compare-page" className="contents">
      <Body />
    </div>
  );
}
```

**Layout invariant.** The wrapper MUST be `display: contents` (Tailwind
`contents`) — confirm via the Compare.routing test by adding a sibling
assertion that the rendered grid still has three children at its
WorkspaceGrid root, not one collapsed cell.

- [ ] **Step 4: Update AppShell.tsx routes.**

In `AppShell.tsx`, swap the imports + the `<Routes>` block. The four
existing routes (`/experiments/:eid/compare`, `/experiments/:eid/compare/:id`,
`/experiments/:eid/compare/new`, `/compare/all/...`) all route to
`<Compare />`. The two `/edit` redirect routes from Task B-1 stay.
`ComparePage` and `ComparePageEdit` still exist as imports inside
`Compare.tsx`.

```tsx
import { Compare } from "../pages/Compare";

// Remove: import { ComparePage } from "../pages/ComparePage";
// Remove: import { ComparePageEdit } from "../pages/ComparePageEdit";

// In the Routes block, replace the four ComparePage/ComparePageEdit
// elements with <Compare /> instead. Keep the redirect routes.
```

- [ ] **Step 5: Run.** Expected: PASS.

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/Compare.routing.test.tsx
cd packages/HimalayaUI/frontend && npm run build 2>&1 | tail -5
```

- [ ] **Step 6: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/pages/Compare.tsx packages/HimalayaUI/frontend/src/components/AppShell.tsx packages/HimalayaUI/frontend/test/Compare.routing.test.tsx
git commit -m "feat(pages): Compare.tsx unified shell (delegating) (Compare UX C-11)"
```

### Task C-12: Replace ComparePage review-mode header with the new components

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/ComparePage.tsx`
- Test: `packages/HimalayaUI/frontend/test/Compare.header.review.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { MemoryRouter } from "react-router-dom";
// ... [imports for ComparePage and any needed mocks]

describe("Compare header (review mode) — Compare UX C-12", () => {
  it("renders the new title strip + toolbar (not the legacy badges)", async () => {
    // Render the page with a loaded comparison in cache via QueryClient.setQueryData.
    // Assertions:
    expect(await screen.findByTestId("compare-title-strip")).toBeInTheDocument();
    expect(await screen.findByTestId("compare-toolbar")).toBeInTheDocument();
    // The legacy badges are gone:
    expect(screen.queryByTestId("comparison-lineage")).toBeNull();
    expect(screen.queryByTestId("comparison-forks-trigger")).toBeNull();
  });
});
```

(NOTE: the full fixture setup mirrors existing Compare component tests
under `packages/HimalayaUI/frontend/test/ComparePage.*.test.tsx`. Copy
that scaffolding.)

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Replace the review-mode header in `ComparePage.tsx`.**

Find the JSX block in `ComparePage.tsx` that renders the header
(`<div data-testid="compare-review-header"...>`). Replace it with mounts
for `CompareTitleStrip`, `CompareStatusSurface`, and `CompareToolbar`,
passing the needed handlers and state. The `EditOrForkButton` callsite
disappears (Save pill from C-7 will cover the gesture in C-14). The
`LineageBadge`, `NeedsReviewBadge`, `ForksPopover` callsites disappear.

The full replacement (large) — replace the existing header block with:

```tsx
<CompareTitleStrip
  title={comparison.title}
  description={comparison.description}
  memberCount={members.length}
  authorUsername={authorUsername}  // fetched via useUser(comparison.created_by)
  isCurrentUserAuthor={isAuthor}
  lastEventAt={comparison.last_event_at}
  forkedFromTitle={comparison.forked_from_title}
  forkedFromHref={forkedFromHref /* compute from forked_from_id */}
  onTitleChange={() => {}}   // unreachable when readOnly is true
  onDescChange={() => {}}
  readOnly                    // review mode — InlineEditableText renders as plain text
/>
<CompareStatusSurface
  needsReview={isStale ? { onReanalyze: handleReanalyze } : null}
  serverUpdate={null}
  savedAt={null}
/>
<CompareToolbar
  groupingControl={<GroupingModeToggle />}
  annotationControl={<AnnotationToggles />}
  forksCount={forks?.length ?? 0}
  onCopyLink={handleCopyLink}
  onDelete={handleDelete}
  onDiscardChanges={null}  // review mode → never dirty → no discard
  onFork={handleFork}
  exportControl={<FigureExportControls .../>}
  saveControl={null}  // review mode → no save
/>
```

Delete the now-unused imports: `LineageBadge`, `NeedsReviewBadge`,
`ForksPopover`, `EditOrForkButton` (inline component). The `EditOrForkButton`
function body deletes too.

- [ ] **Step 4: Run.**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/Compare.header.review.test.tsx
cd packages/HimalayaUI/frontend && npm run build 2>&1 | tail -5
```

Expected: PASS + build clean. (Some previous tests for LineageBadge / NeedsReviewBadge / ForksPopover may need updates — fix them as part of this task: the affected E2E tests can be migrated to the new selectors `data-testid="compare-title-strip"`, `data-testid="compare-toolbar-more"`, etc.)

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/pages/ComparePage.tsx packages/HimalayaUI/frontend/test/Compare.header.review.test.tsx
git commit -m "refactor(compare): swap review header for TitleStrip/Status/Toolbar (Compare UX C-12)"
```

### Task C-13: Replace ComparePageEdit header with new components + SavePill

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx`
- Test: `packages/HimalayaUI/frontend/test/Compare.header.edit.test.tsx`

- [ ] **Step 1: Failing test (Save pill replaces triplet).**

```tsx
describe("Compare edit header — Compare UX C-13", () => {
  it("renders Save pill in place of Save/Cancel/Discard triplet", async () => {
    // ... fixture setup that mounts ComparePageEdit with an active draft
    expect(await screen.findByTestId("save-pill")).toBeInTheDocument();
    expect(screen.queryByTestId("comparison-cancel")).toBeNull();
    expect(screen.queryByTestId("comparison-discard")).toBeNull();
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Replace the edit-mode header in `ComparePageEdit.tsx`.**

Locate the `editCenter` block (around line 290+). Replace its top section
(the title input + Save/Cancel/Discard buttons) with `CompareTitleStrip`
bound to `setDraftTitle` / `setDraftDescription`, `CompareStatusSurface`,
and `CompareToolbar` with the new `SavePill` mounted in the `saveControl`
slot. Pass `onDiscardChanges={handleDiscard}` so "Discard changes" surfaces
inside the ⋯-more menu.

The `<input data-testid="compare-edit-title">` element + the
`data-testid="comparison-save"` / `comparison-cancel` /
`comparison-discard` triplet are all replaced. The Cmd+Enter / Esc key
bindings on the SavePill keep working — wire them in the page's
`useGlobalShortcuts` registration block.

The Discard pop-confirm: today's `handleDiscard` should be guarded by a
small `confirm()` since it's now a less-prominent gesture (inside a
menu).

- [ ] **Step 4: Run.** Expected: PASS.

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/Compare.header.edit.test.tsx
cd packages/HimalayaUI/frontend && npm run build 2>&1 | tail -5
```

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx packages/HimalayaUI/frontend/test/Compare.header.edit.test.tsx
git commit -m "refactor(compare): edit header → Title/Status/Toolbar/SavePill (Compare UX C-13)"
```

### Task C-14: Wire the editing-as-fork-of Save pill flow

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx` (and/or wherever `handleSave` lives)
- Modify: `packages/HimalayaUI/frontend/src/lib/comparison/draftFactories.ts`
- Modify: `packages/HimalayaUI/frontend/src/state.ts` (add `setDraftForkOf` action)
- Test: `packages/HimalayaUI/frontend/test/Compare.savePill.fork.test.tsx`

Non-author edits land here. The SavePill copy reads "Save as fork…"; on
click, prompt for a fork title (default `Copy of <parent title>`), morph
the draft into a fork (clear id, clear baseHash, set forkedFromId +
forkedAtHash), and submit via the create path. The existing
`startFork(comparison, qc)` factory in `draftFactories.ts` already does
most of this; we extract the rename step.

**Mode state machine note** (per spec §3 update). The pre-prompt state is
`editing-as-fork-of` (draft.id is set, points at the parent). After
`setDraftForkOf` morphs the draft, the next render observes `draft.id ===
undefined` AND `draft.forkedFromId !== undefined`, which is the
`creating-from-fork` variant of the union (already implemented by C-5).
The submit then routes through `saveComparisonMutator` with no
`expected_content_hash` (create path), so the fork commits without a
conflict possibility.

The `useCompareMode` hook from C-5 already routes both
`creating-from-fork` and `creating-blank` correctly; C-14 adds no hook
changes — only the `setDraftForkOf` action and the `handleSave` morph
branch below.

- [ ] **Step 1: Failing test.**

```tsx
describe("Compare Save-as-fork flow — Compare UX C-14", () => {
  it("submits via create path with fork lineage when user is non-author", async () => {
    // Set draft.id = 7, currentUserId != comparison.created_by, click SavePill
    // -> prompt UI appears -> confirm "My fork" -> assert request body has
    //    forked_from_id = 7, forked_at_hash = "<baseHash>", id absent
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

In `handleSave`, branch on `mode.kind`. When `editing-as-fork-of`:

```ts
async function handleSave() {
  if (mode.kind === "editing-as-fork-of" && draft !== null) {
    const baseTitle = comparison?.title ?? "comparison";
    // window.prompt note: Playwright auto-dismisses dialogs by default.
    // Any e2e covering this flow MUST register a dialog handler before
    // triggering the gesture:
    //   page.on("dialog", d => d.accept("My fork"));
    // Alternative (recommended if e2e coverage grows): replace with an
    // inline modal that uses useFocusTrap — Phase 2 follow-up.
    const proposed = window.prompt("Title for your fork:", `Copy of ${baseTitle}`);
    if (proposed === null) return;  // cancelled
    // Morph the draft in-place:
    setDraftForkOf({
      newTitle: proposed.trim() === "" ? `Copy of ${baseTitle}` : proposed,
      sourceId: mode.parentId,
      sourceHash: draft.baseHash ?? "",
    });
    // setDraftForkOf is a new Zustand action that clears `id` + `baseHash`,
    // sets `forkedFromId` and `forkedAtHash`, and sets the new title.
    //
    // CRITICAL: re-read the morphed draft from the store BEFORE building
    // the save payload below. If we reuse a payload object built earlier
    // (e.g., from a closure that captured `draft` pre-morph), the second
    // submission inherits the original draft's id and the create-path
    // route never fires. Always build the payload fresh after the morph.
    const fresh = useAppState.getState().activeDraft;
    if (fresh === null) return;
    save.mutate(buildSavePayload(fresh));
    return;
  }
  // ... non-fork path: build payload from `draft` directly and mutate.
}
```

**Why this matters for idempotency.** `useQueueMutation.mutationFn`
mints `client_op_id` PER `mutate()` CALL (per CLAUDE.md), so the two
back-to-back submissions inherent in the morph-then-save flow naturally
get distinct keys — but only if each `mutate()` invocation gets its own
SaveComparisonInput object. If a closure captures a single payload and
both code paths reuse it, the two ops collide on one
`idempotent_responses` row and the second submission silently no-ops.
Build the payload immediately before each `mutate()`.

Add `setDraftForkOf` to `state.ts`:

```ts
setDraftForkOf: (p: { newTitle: string; sourceId: number; sourceHash: string }) => {
  set((s) => {
    if (s.activeDraft === null) return s;
    return { activeDraft: {
      ...s.activeDraft,
      id: undefined,
      baseHash: undefined,
      forkedFromId: p.sourceId,
      forkedAtHash: p.sourceHash,
      title: p.newTitle,
    } };
  });
},
```

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx packages/HimalayaUI/frontend/src/state.ts packages/HimalayaUI/frontend/test/Compare.savePill.fork.test.tsx
git commit -m "feat(compare): SavePill 'Save as fork' flow (Compare UX C-14)"
```

### Task C-15: Inline merge — fold `ComparePage.tsx` + `ComparePageEdit.tsx` bodies into `Compare.tsx`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/Compare.tsx`
- Delete: `packages/HimalayaUI/frontend/src/pages/ComparePage.tsx`
- Delete: `packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx`
- Test: existing tests cover behaviour.

At this point both pages share the same header components (C-12, C-13)
and Save semantics (C-14). The remaining differences are mostly bookkeeping
(handlers, refs, slot wiring). This task is the merge.

- [ ] **Step 1: Move the bodies.**

Copy the contents of `ComparePage.tsx` and `ComparePageEdit.tsx` into
`Compare.tsx`. Reconcile state: use `useCompareMode` (C-5) +
`useCompareDraftDirty` (C-6) for branching; share the same plot mount,
gutter mount, and right-rail mount (chat OR picker chip — picker chip is
interim per the Phase 2 deferral).

The picker-mount in "edit-mode-was-here" reads:

```tsx
{mode.kind !== "viewing" && (
  <div className="...">
    {/* Interim picker chip — Phase 2 will move this to a right-rail tab */}
    <ComparisonPickerPanel ... />
  </div>
)}
```

- [ ] **Step 2: Update `Compare.tsx`'s import block.**

Remove imports of `ComparePage` and `ComparePageEdit`. Add direct imports
of the components those pages used.

- [ ] **Step 3: Delete the two page files.**

```bash
git rm packages/HimalayaUI/frontend/src/pages/ComparePage.tsx packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx
```

- [ ] **Step 4: Wire the §7.1 empty state.**

Per spec §7.1, when `members.length === 0` the center pane renders a
three-element empty block: "No traces yet." headline, a `+ Add traces`
button that auto-flips the right rail to the picker, and "Or drag
exposures from the picker." copy. The whole region also acts as a drop
target for external exposure drags.

Add a failing test FIRST, then implement:

```tsx
// test/Compare.emptyState.test.tsx
describe("Compare empty state — Compare UX C-15", () => {
  it("renders the empty block when a draft has zero members", () => {
    // Seed Zustand activeDraft with empty members
    // render(<Compare />)
    expect(screen.getByTestId("compare-empty-state")).toBeInTheDocument();
    expect(screen.getByText(/No traces yet/i)).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /Add traces/i })).toBeInTheDocument();
    expect(screen.getByText(/Or drag exposures/i)).toBeInTheDocument();
  });

  it("clicking + Add traces opens the interim picker chip", () => {
    // Phase 1 does NOT yet have a right-rail Chat/Add tab system
    // (that's Phase 2 per spec §5). The Phase 1 surface is the interim
    // ComparisonPickerPanel chip wired into Compare.tsx — toggle it
    // via the existing Zustand action `setPickerOpen(true)` (or
    // whatever the interim chip uses).
    // 1) render(<Compare />) with an empty draft.
    // 2) Click the "+ Add traces" button.
    // 3) Assert the ComparisonPickerPanel is now visible / open
    //    (data-testid="comparison-picker-panel") OR the Zustand
    //    pickerOpen field flipped to true. Don't assert against a
    //    "right-rail tab" state that doesn't exist in Phase 1.
  });

  it("renders as a drop target for external exposure drags", () => {
    // dragOver the empty region with a synthetic exposure DataTransfer payload
    // assert the empty region carries data-drop-target="true"
    // (or that the underlying drop handler is wired — verify via mock)
  });
});
```

Then add the empty block to `Compare.tsx`:

```tsx
{members.length === 0 && (
  <div
    data-testid="compare-empty-state"
    onDragOver={handleExposureDragOver}
    onDrop={handleExposureDrop}
    className="..."
  >
    <h2>No traces yet.</h2>
    <button type="button" data-interactable="button" onClick={handleAddTraces}>
      + Add traces
    </button>
    <p>Or drag exposures from the picker.</p>
  </div>
)}
```

**`handleAddTraces` Phase 1 wiring.** The right-rail Chat ↔ + Add
tabs are Phase 2 (spec §5). For Phase 1, route the button to the
interim `ComparisonPickerPanel` chip that already mounts inside
`Compare.tsx` (see C-15 Step 1's picker-mount block — `mode.kind !==
"viewing"` opens the panel).

**Mode-aware draft creation (load-bearing).** The handler must branch
on `mode.kind` — calling `startNewDraft()` unconditionally would mint
a BLANK draft (`id: undefined`, `members: []`, no `forkedFromId`) and
the next save would create an unrelated new comparison instead of
editing or forking the one being viewed. The draft must be seeded
from the current comparison so its identity survives.

```ts
const handleAddTraces = () => {
  if (mode.kind === "viewing" || mode.kind === "viewing-stale") {
    // Seed the draft from the loaded comparison so the user is
    // editing-or-forking THIS comparison, not creating a new one.
    // `loadDraftFromComparison` (draftFactories.ts) builds an
    // ActiveDraft with `id: comparison.id`, members copied, etc.
    // The author-vs-non-author distinction is resolved later when
    // SavePill copy renders via useCompareMode → editing-mine vs
    // editing-as-fork-of.
    if (comparison !== undefined) {
      useAppState.getState().setActiveDraft(loadDraftFromComparison(comparison));
    }
  } else if (mode.kind === "creating-blank" || mode.kind === "creating-from-fork") {
    // Already on a draft surface — the picker chip should already be
    // visible; nothing to do.
    return;
  } else {
    // editing-mine / editing-as-fork-of: same as creating — already a
    // draft surface.
    return;
  }
};
```

The handler reads `mode` and `comparison` from the same scope the rest
of `Compare.tsx` does. If `comparison` is undefined (still loading),
the button is hidden anyway because the empty state only renders when
`members.length === 0` AND there IS a comparison (or a draft) to attach
to.

(The drop-target wiring uses the same handler the gutter uses; if
external-add is fully Phase 2, the empty-state drop handler can route
to the same "interim picker" path that other drag sites use.)

- [ ] **Step 4b: Wire SSE → `viewing-stale` detection.**

`useCompareMode` accepts a `staleAgainstHash?: string` prop that
returns `{ kind: "viewing-stale", previousHash }`. The hook does NOT
listen to SSE itself; `Compare.tsx` must compute the prop.

Add a small hook that watches the comparison cache row's `content_hash`
and tracks the last hash this client *acknowledged* (i.e., last seen
while no draft was active). When the cache `content_hash` drifts beyond
the acknowledged value WHILE no draft is active, the prop fires:

```tsx
function useStaleAgainstHash(comparison: Comparison | undefined): {
  previousHash: string | undefined;
  acknowledge: () => void;
} {
  const draft = useAppState((s) => s.activeDraft);
  // useState (not useRef) — first-sighting + rebase events MUST trigger
  // a re-render so the banner appears/clears deterministically.
  const [acked, setAcked] = useState<string | undefined>(undefined);

  useEffect(() => {
    if (comparison === undefined) return;
    // First sighting (no acked yet) — anchor to the current hash whether
    // we're viewing or editing. Spec §3 baseHash-on-first-edit is owned
    // by `draft.baseHash`, not this hook.
    //
    // KNOWN EDGE CASE (Phase 1 acceptable): if React batches the initial
    // cache load AND a fast SSE drift into a single render, the first
    // observed `comparison.content_hash` IS the post-drift value, and
    // the drift is silently swallowed. Recovering the pre-drift hash
    // would require subscribing to the React-Query cache observer
    // directly (rather than reading the rendered prop) — architecturally
    // larger; deferred to Phase 2 alongside actor attribution.
    if (acked === undefined) {
      setAcked(comparison.content_hash);
      return;
    }
    // No more auto-advancement here. Two events advance `acked`:
    //   1. `acknowledge()` (user clicked the "Reload"/"Acknowledge" CTA).
    //   2. Draft-clear after own-op completion (separate effect below).
    // Foreign drift while viewing leaves `acked` pinned so the banner
    // stays visible until the user dismisses or acks.
  }, [comparison?.content_hash, acked]);

  // Rebase on own-op completion: when the user's own draft completes
  // (Zustand activeDraft transitions non-null → null), the new hash IS
  // the truth they just authored. Without this rebase, every save
  // would flash `viewing-stale` for one render.
  //
  // `prevDraftRef` initialised to `draft` (NOT `null`) so the initial
  // mount doesn't fire a spurious rebase, and so React 18 StrictMode's
  // double-invoke of this effect's cleanup+re-run is idempotent: the
  // ref re-receives the same `draft` value on the second pass.
  const prevDraftRef = useRef<typeof draft>(draft);
  useEffect(() => {
    if (prevDraftRef.current !== null && draft === null && comparison !== undefined) {
      setAcked(comparison.content_hash);
    }
    prevDraftRef.current = draft;
  }, [draft, comparison?.content_hash]);

  const acknowledge = useCallback(() => {
    if (comparison !== undefined) setAcked(comparison.content_hash);
  }, [comparison?.content_hash]);

  if (draft !== null) return { previousHash: undefined, acknowledge };  // editing — own draft owns the hash
  if (comparison === undefined) return { previousHash: undefined, acknowledge };
  if (acked === undefined) return { previousHash: undefined, acknowledge };
  return {
    previousHash: acked === comparison.content_hash ? undefined : acked,
    acknowledge,
  };
}
```

**Three behavioural pins** that future refactors must preserve:
1. **First-sighting anchor.** `acked` initialises to the first observed
   `content_hash` — anchors detection so the banner doesn't flash on
   page load.
2. **User acknowledgement clears the banner.** `acknowledge()` rebases
   `acked` to the current hash; banner disappears on next render. The
   `CompareStatusSurface` "Reload" / "Acknowledge" button wires to this.
3. **Own-op completion does NOT flash the banner.** When the user
   saves their own draft, `activeDraft` transitions to `null` and the
   cache `content_hash` drifts simultaneously. The own-op-rebase
   effect re-anchors before the next render so `viewing-stale` is not
   emitted for a single frame.

Pass the result into `useCompareMode({ ..., staleAgainstHash })` and
into `CompareStatusSurface`'s `serverUpdate` prop. **Confirm C-9 ships
the matching `ServerUpdate` interface** (already updated in C-9 Step 3
to `{ previousHash, onAcknowledge }` — the two changes land together
and share the same wire shape):

```ts
interface ServerUpdate {
  previousHash: string;
  onAcknowledge: () => void;
}
```

And mount:

```tsx
const { previousHash, acknowledge } = useStaleAgainstHash(comparison);
const mode = useCompareMode({ comparison, currentUserId, staleAgainstHash: previousHash });
// ...
<CompareStatusSurface
  needsReview={...}
  serverUpdate={previousHash !== undefined
    ? { previousHash, onAcknowledge: acknowledge }
    : null}
  savedAt={...}
/>
```

The `acknowledge` callback is stable (`useCallback` + the same
`comparison.content_hash` dep) so passing it down doesn't re-mount
`CompareStatusSurface` each render.

Add a failing test that mounts `Compare.tsx`, sets initial cache, then
updates the cache via `qc.setQueryData` (simulating SSE arrival), and
asserts the banner appears + `useCompareMode` returns
`viewing-stale`:

```tsx
// test/Compare.viewingStale.test.tsx
describe("Compare viewing-stale wiring — Compare UX C-15", () => {
  it("flips to viewing-stale when content_hash drifts while no draft is active", async () => {
    // 1) Seed qc with comparison content_hash = "h1"
    // 2) render(<Compare />), assert CompareStatusSurface has no banner
    // 3) qc.setQueryData(queryKeys.comparison(id), {...prev, content_hash: "h2"})
    // 4) Await re-render, assert findByTestId("compare-status-server-update")
    //    is present (the inner banner — coarser "compare-status-surface"
    //    would also succeed when only NeedsReview or Saved renders).
    // 5) Assert useCompareMode emitted kind === "viewing-stale"
  });

  it("banner clears when the user acknowledges (Reload / Acknowledge)", async () => {
    // 1) Seed qc with comparison content_hash = "h1", no draft.
    // 2) Drift via qc.setQueryData → content_hash = "h2". Assert banner present.
    // 3) Click the acknowledge CTA (data-testid="compare-status-acknowledge").
    // 4) Assert banner is gone AND useCompareMode emits kind === "viewing".
    //    The acked hash is rebased to "h2"; subsequent renders against "h2"
    //    no longer return previousHash.
  });

  it("does NOT flash viewing-stale during own-op save completion", async () => {
    // 1) Seed qc with content_hash = "h1", set activeDraft (own draft).
    // 2) Simulate own-op save success: qc.setQueryData → content_hash = "h2"
    //    AND useAppState.setState({ activeDraft: null }) (same React batch).
    // 3) Assert across renders that useCompareMode NEVER emits "viewing-stale".
    //    The own-op-rebase effect anchors acked to "h2" before the next paint.
  });
});
```

Without this step the `viewing-stale` variant of the union is
unreachable in practice — the hook accepts the prop but nothing
computes it. Spec §3 explicitly names this contract. The three test
cases pin the three behavioural rules from the hook's "Three
behavioural pins" callout above.

- [ ] **Step 4c: Add the `baseHash` capture-moment contract test.**

Per spec §3, `draft.baseHash` is captured at *first edit*, NOT at
edit-mode entry. Foreign updates arriving BEFORE the first edit refresh
the cache (so the draft is built from the latest comparison); foreign
updates AFTER the first edit do NOT mutate `baseHash` (so the user's
conflict-detection anchor is stable).

Create `test/Compare.baseHashCaptureMoment.test.tsx`:

```tsx
describe("baseHash capture moment — Compare UX C-15 / spec §3", () => {
  it("captures baseHash on first edit, not at mode entry", async () => {
    // 1) Seed qc with comparison content_hash = "h1", no draft.
    // 2) User clicks the title → setDraftTitle("renamed") — first edit.
    // 3) Assert useAppState.getState().activeDraft?.baseHash === "h1".
  });

  it("foreign update BEFORE first edit refreshes the draft baseline", () => {
    // 1) Seed qc with content_hash = "h1", no draft.
    // 2) Foreign SSE updates qc to content_hash = "h2" (no draft yet, so
    //    activeDraft stays null and the next 'first edit' will use h2).
    // 3) User edits title → setDraftTitle("renamed").
    // 4) Assert activeDraft.baseHash === "h2".
  });

  it("foreign update AFTER first edit does NOT mutate baseHash", () => {
    // 1) Seed qc with content_hash = "h1", no draft.
    // 2) User edits title → activeDraft.baseHash === "h1".
    // 3) Foreign SSE updates qc to content_hash = "h2".
    // 4) Assert activeDraft.baseHash is STILL "h1" (the user's anchor).
    //    Assert mode kind transitions to viewing-stale or conflict-pending
    //    per the §3 state machine.
  });
});
```

This is the load-bearing test for spec §3's edit-mode entry rules.
Spec L96-L97 names this file explicitly: `Compare.baseHashCaptureMoment.test.tsx`.

- [ ] **Step 5: Run.**

```bash
cd packages/HimalayaUI/frontend && npm run build 2>&1 | tail -15
cd packages/HimalayaUI/frontend && npm test -- --run > /tmp/vitest-C15.out 2>&1
grep -E "Tests:|Failed|passed" /tmp/vitest-C15.out | tail -5
```

Expected: clean build, all tests pass. Any leftover test referencing
`ComparePage`/`ComparePageEdit` by component name needs an update to
import `Compare` instead.

- [ ] **Step 6: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/pages/Compare.tsx packages/HimalayaUI/frontend/test/Compare.emptyState.test.tsx packages/HimalayaUI/frontend/test/Compare.viewingStale.test.tsx packages/HimalayaUI/frontend/test/Compare.baseHashCaptureMoment.test.tsx
git commit -m "refactor(compare): merge ComparePage(+Edit) + empty/stale/baseHash (Compare UX C-15)"
```

### Task C-16: Delete folded components

**Files:**
- Delete: `packages/HimalayaUI/frontend/src/components/LineageBadge.tsx` (content lives in `CompareTitleStrip`)
- Delete: `packages/HimalayaUI/frontend/src/components/NeedsReviewBadge.tsx` (lives in `CompareStatusSurface`)
- Delete: `packages/HimalayaUI/frontend/src/components/ForksPopover.tsx` (folded into toolbar — Task C-17)

- [ ] **Step 1: Confirm no stale imports.**

```bash
grep -rn "LineageBadge\|NeedsReviewBadge\|ForksPopover" packages/HimalayaUI/frontend/src \
  --include="*.ts" --include="*.tsx" | grep -v node_modules
```

Expected: no results.

If any references remain in test files, update them to assert on the new
testids (`compare-title-strip`, `compare-status-resnapshot`, etc.) or
delete the tests if they're component-internal.

- [ ] **Step 2: Delete the three files.**

```bash
git rm packages/HimalayaUI/frontend/src/components/LineageBadge.tsx packages/HimalayaUI/frontend/src/components/NeedsReviewBadge.tsx packages/HimalayaUI/frontend/src/components/ForksPopover.tsx
```

- [ ] **Step 3: Run.**

```bash
cd packages/HimalayaUI/frontend && npm run build 2>&1 | tail -10
cd packages/HimalayaUI/frontend && npm test -- --run > /tmp/vitest-C16.out 2>&1
grep -E "Tests:|Failed|passed" /tmp/vitest-C16.out | tail -3
```

Expected: clean.

- [ ] **Step 4: Commit.**

```bash
git commit -m "cleanup(compare): delete folded badge/popover components (Compare UX C-16)"
```

### Task C-17: Forks dropdown inside toolbar ⋯ more

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/CompareToolbar.tsx`
- Test: `packages/HimalayaUI/frontend/test/CompareToolbar.test.tsx`

- [ ] **Step 1: Failing test.**

Add to `CompareToolbar.test.tsx`:

```tsx
it("renders forks list inside ⋯ more dropdown when provided", () => {
  render(wrap(<CompareToolbar
    groupingControl={null} annotationControl={null}
    forksList={[
      { id: 9, title: "fork-9", href: "/compare/all/9" },
      { id: 10, title: "fork-10", href: "/compare/all/10" },
    ]}
    onCopyLink={() => {}} onDelete={() => {}}
    onDiscardChanges={null} onFork={() => {}}
    exportControl={null} saveControl={null}
  />));
  fireEvent.click(screen.getByTestId("compare-toolbar-more"));
  expect(screen.getByText("fork-9")).toBeInTheDocument();
  expect(screen.getByText("fork-10")).toBeInTheDocument();
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Replace `forksCount` prop with `forksList`.**

In `CompareToolbar`, change `forksCount: number` to
`forksList: Array<{ id: number; title: string; href: string }>`. Replace
the "Forks (N)" button with an expand-on-hover sub-section showing each
fork as a link. Update `Compare.tsx` to derive `forksList` from
`useComparisonForks(id)`.

**Test migration callout.** Task C-10's tests construct `CompareToolbar`
with the old `forksCount: number` prop. As part of this step, also
update those tests to use `forksList: Array<...>` (a count test becomes
a list-length test, etc.). The commit below stages
`packages/HimalayaUI/frontend/test/CompareToolbar.test.tsx` alongside
the component change so the suite stays green between C-10 and C-17.

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/CompareToolbar.tsx packages/HimalayaUI/frontend/src/pages/Compare.tsx packages/HimalayaUI/frontend/test/CompareToolbar.test.tsx
git commit -m "feat(compare): forks list inside toolbar ⋯-more menu (Compare UX C-17)"
```

### Task C-18: Phase C regression sweep

**Files:** (none — runs the suites)

- [ ] **Step 1.** `npm run build` clean.
- [ ] **Step 2.** `npm test -- --run` clean.
- [ ] **Step 3.** Julia suite still clean (no backend changes since Phase A, but worth re-running before Phase E).

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-C.out 2>&1
grep -E "did not pass|fail" /tmp/jl-C.out | head -10
```

Expected: clean.

If anything regressed, fix before continuing.

---

## Phase E — Member row redesign

Phase goal: rebuild `MemberMetaRow` + `MemberMetaGutter` per spec §7.2-7.4.
Collapsed/expanded rows, right-edge action zone, grab-anywhere drag with
4px threshold, visible inter-row resize gap.

### Task E-1: `RowActionZone` primitive

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/RowActionZone.tsx`
- Test: `packages/HimalayaUI/frontend/test/RowActionZone.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { RowActionZone } from "../src/components/RowActionZone";

describe("RowActionZone — Compare UX E-1", () => {
  it("renders ⋯ overflow and ⋮⋮ drag cue", () => {
    render(<RowActionZone onOverflow={() => {}}/>);
    expect(screen.getByTestId("row-action-overflow")).toBeInTheDocument();
    expect(screen.getByTestId("row-action-drag-cue")).toBeInTheDocument();
  });
  it("dispatches overflow click", () => {
    const onOverflow = vi.fn();
    render(<RowActionZone onOverflow={onOverflow}/>);
    fireEvent.click(screen.getByTestId("row-action-overflow"));
    expect(onOverflow).toHaveBeenCalled();
  });
  it("⋮⋮ is signage (no click handler runs)", () => {
    const onOverflow = vi.fn();
    render(<RowActionZone onOverflow={onOverflow}/>);
    fireEvent.click(screen.getByTestId("row-action-drag-cue"));
    expect(onOverflow).not.toHaveBeenCalled();
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

```tsx
interface Props {
  onOverflow: () => void;
}

export function RowActionZone({ onOverflow }: Props): JSX.Element {
  return (
    <div className="inline-flex items-center gap-1 text-fg-dim">
      <button
        type="button"
        data-testid="row-action-overflow"
        data-interactable="button"
        onClick={(e) => { e.stopPropagation(); onOverflow(); }}
        className="px-1 hover:text-fg hover:bg-bg-hover rounded"
        title="Row actions"
      >
        ⋯
      </button>
      <span
        data-testid="row-action-drag-cue"
        aria-hidden="true"
        className="px-1 select-none text-fg-dim cursor-grab"
      >
        ⋮⋮
      </span>
    </div>
  );
}
```

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/RowActionZone.tsx packages/HimalayaUI/frontend/test/RowActionZone.test.tsx
git commit -m "feat(ui): RowActionZone primitive (Compare UX E-1)"
```

### Task E-2: Collapse/expand state on `MemberMetaRow`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/MemberMetaRow.tsx`
- Modify: `packages/HimalayaUI/frontend/src/components/MemberMetaGutter.tsx`
- Test: `packages/HimalayaUI/frontend/test/MemberMetaRow.collapse.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
describe("MemberMetaRow collapse/expand — Compare UX E-2", () => {
  it("renders collapsed view by default", () => {
    // ... fixture-mount one MemberMetaRow
    expect(screen.getByTestId("member-meta-row")).toHaveAttribute("data-expanded", "false");
  });

  it("expands when clicked on the row body", () => {
    fireEvent.click(screen.getByTestId("member-meta-row-body"));
    expect(screen.getByTestId("member-meta-row")).toHaveAttribute("data-expanded", "true");
  });

  it("collapses when clicked again", () => {
    fireEvent.click(screen.getByTestId("member-meta-row-body"));
    fireEvent.click(screen.getByTestId("member-meta-row-body"));
    expect(screen.getByTestId("member-meta-row")).toHaveAttribute("data-expanded", "false");
  });

  it("only one row expanded at a time across siblings", () => {
    // Mount two rows in the same gutter; expand row A; click row B; assert A collapsed.
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

`MemberMetaRow` gains `expanded: boolean` and `onToggleExpand: () => void`
props (controlled). The single-row-expanded invariant lives on the
parent (`MemberMetaGutter`) which holds a single
`expandedMemberId: number | null` state and threads `expanded={expandedMemberId === member.id}` into each row + `onToggleExpand` that
sets / clears that state.

Pseudocode for `MemberMetaRow`:

```tsx
return (
  <div
    data-testid="member-meta-row"
    data-expanded={expanded}
    data-member-id={member.id}
    data-interactable="expand"
  >
    <div
      data-testid="member-meta-row-body"
      onClick={onToggleExpand}
      className="..."
    >
      <span aria-hidden="true">{expanded ? "▾" : "▸"}</span>
      {/* label + meta */}
    </div>
    <RowActionZone onOverflow={() => setOverflowOpen((v) => !v)} />
    {expanded && (
      <div>
        {/* existing per-member control widgets (label/color/norm/q-win/peaks) */}
      </div>
    )}
  </div>
);
```

The expanded body keeps the existing per-member control widgets — we are
NOT redesigning those individually in this plan (color swatch upgrade is
Phase 2 deferral). They simply move from "always visible" to "visible
when expanded".

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/MemberMetaRow.tsx packages/HimalayaUI/frontend/src/components/MemberMetaGutter.tsx packages/HimalayaUI/frontend/test/MemberMetaRow.collapse.test.tsx
git commit -m "feat(member-row): collapse/expand with single-expanded invariant (Compare UX E-2)"
```

### Task E-3: Grab-anywhere drag with threshold

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/MemberMetaRow.tsx`
- Test: `packages/HimalayaUI/frontend/test/MemberMetaRow.drag.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, fireEvent } from "@testing-library/react";
import { MemberMetaRow } from "../src/components/MemberMetaRow";

// JSDOM note: `PointerEvent` is not a real constructor in JSDOM, but
// fireEvent.pointerDown/Move/Up dispatches an event whose `clientX/Y`
// are populated from the second arg. The component's handler reads
// `e.clientX` directly, so this is sufficient.

describe("MemberMetaRow drag-vs-click threshold — Compare UX E-3", () => {
  it("treats a pointer-up within 4px as toggle-expand", () => {
    const onToggleExpand = vi.fn();
    const onDragStart = vi.fn();
    const { getByTestId } = render(
      <MemberMetaRow
        /* …member fixture… */
        onToggleExpand={onToggleExpand}
        onDragStart={onDragStart}
      />,
    );
    const row = getByTestId("member-meta-row");
    fireEvent.pointerDown(row, { clientX: 10, clientY: 10, pointerId: 1 });
    fireEvent.pointerMove(row, { clientX: 12, clientY: 10, pointerId: 1 });
    fireEvent.pointerUp(row,   { clientX: 12, clientY: 10, pointerId: 1 });
    expect(onToggleExpand).toHaveBeenCalledTimes(1);
    expect(onDragStart).not.toHaveBeenCalled();
  });

  it("treats a pointer-move beyond 4px as drag (NOT toggle-expand)", () => {
    const onToggleExpand = vi.fn();
    const onDragStart = vi.fn();
    const { getByTestId } = render(
      <MemberMetaRow
        /* …member fixture (id = 7)… */
        onToggleExpand={onToggleExpand}
        onDragStart={onDragStart}
      />,
    );
    const row = getByTestId("member-meta-row");
    fireEvent.pointerDown(row, { clientX: 10, clientY: 10, pointerId: 1 });
    fireEvent.pointerMove(row, { clientX: 20, clientY: 10, pointerId: 1 });  // |Δx|=10 > 4
    fireEvent.pointerUp(row,   { clientX: 20, clientY: 10, pointerId: 1 });
    expect(onDragStart).toHaveBeenCalledWith(7);
    expect(onToggleExpand).not.toHaveBeenCalled();
  });

  it("uses Manhattan distance — diagonal movement past threshold counts as drag", () => {
    const onToggleExpand = vi.fn();
    const onDragStart = vi.fn();
    const { getByTestId } = render(
      <MemberMetaRow /* …member… */
                     onToggleExpand={onToggleExpand} onDragStart={onDragStart} />,
    );
    const row = getByTestId("member-meta-row");
    fireEvent.pointerDown(row, { clientX: 10, clientY: 10, pointerId: 1 });
    fireEvent.pointerMove(row, { clientX: 13, clientY: 13, pointerId: 1 });  // |Δx|+|Δy|=6 > 4
    fireEvent.pointerUp(row,   { clientX: 13, clientY: 13, pointerId: 1 });
    expect(onDragStart).toHaveBeenCalledTimes(1);
    expect(onToggleExpand).not.toHaveBeenCalled();
  });
});
```

**Why these assertions matter.** A naive test that fires only
`fireEvent.click(row)` would pass `toggleExpand was called` trivially —
it can't distinguish "drag suppressed expand" from "no events fired."
Explicit `pointerDown → pointerMove → pointerUp` with `clientX/Y` is
the only way to exercise the threshold gate. We also assert `onDragStart`
is called with the member id, not just "called" — confirms the
handler reads from the row's bound member rather than a stale closure.

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement using `dragThreshold.ts` (Task C-1).**

In `MemberMetaRow`, replace the simple onClick toggle with a
pointer-event handler that uses `makeDragThresholdState`. On
`pointer-up → "click"`: call `onToggleExpand`. On
`pointermove → "drag-start"`: invoke the dataTransfer setup the gutter
currently uses for reorder. (Or fire an `onDragStart(memberId)` prop into
the gutter.)

Pseudocode:

```tsx
const dragState = useRef(makeDragThresholdState({ thresholdPx: 4 }));

const onPointerDown = (e: React.PointerEvent) => {
  dragState.current.onPointerDown(e.clientX, e.clientY);
  // ... no drag-start yet
};
const onPointerMove = (e: React.PointerEvent) => {
  const step = dragState.current.onPointerMove(e.clientX, e.clientY);
  if (step === "drag-start") onDragStart(member.id);
};
const onPointerUp = (e: React.PointerEvent) => {
  const outcome = dragState.current.onPointerUp(e.clientX, e.clientY);
  if (outcome === "click") onToggleExpand();
};
```

NOTE: the existing HTML5 drag-and-drop in `MemberMetaGutter` uses
`dragstart` events. The interaction between manual pointer threshold and
the browser's native `dragstart` needs care: one common approach is to
call `e.preventDefault()` on `pointerdown`, manage drag manually via a
"floating" portal, and emit drop events at pointer-up positions over the
gutter. An alternative — keep the current HTML5 wiring — is to make the
row's `draggable` attribute toggle on / off via the threshold (set
`draggable=true` only AFTER threshold crossing). Implementer should
choose the simpler path that doesn't regress E2E.

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/MemberMetaRow.tsx packages/HimalayaUI/frontend/test/MemberMetaRow.drag.test.tsx
git commit -m "feat(member-row): grab-anywhere + 4px threshold (Compare UX E-3)"
```

### Task E-4: Visible inter-row resize gap

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/MemberMetaGutter.tsx`
- Modify: `packages/HimalayaUI/frontend/src/components/BandResizeDivider.tsx` (simplify or fold into gutter)
- Test: `packages/HimalayaUI/frontend/test/MemberMetaGutter.resize.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, fireEvent } from "@testing-library/react";
import { MemberMetaGutter } from "../src/components/MemberMetaGutter";

// JSDOM has no layout engine — `getComputedStyle().height` returns ""
// for any element without an inline style, and `getBoundingClientRect`
// returns zero rects. Assert on `data-state` (component-controlled) and
// inline `style.height` (Tailwind class strings are forbidden per
// CLAUDE.md and would change with styling refactors anyway).

describe("Inter-row resize gap — Compare UX E-4", () => {
  it("renders the gap between members in rest state", () => {
    const { getByTestId } = render(<MemberMetaGutter /* …two-member fixture… */ />);
    const gap = getByTestId("member-resize-gap");
    expect(gap).toHaveAttribute("data-state", "rest");
    // Inline style is set by the component (4px rest, 12px hover) so we
    // can assert it without depending on a CSS layout pass.
    expect(gap.style.height).toBe("4px");
  });

  it("flips to hover state on pointer-enter and resets on pointer-leave", () => {
    const { getByTestId } = render(<MemberMetaGutter /* …two-member fixture… */ />);
    const gap = getByTestId("member-resize-gap");
    fireEvent.pointerEnter(gap);
    expect(gap).toHaveAttribute("data-state", "hover");
    expect(gap.style.height).toBe("12px");
    fireEvent.pointerLeave(gap);
    expect(gap).toHaveAttribute("data-state", "rest");
    expect(gap.style.height).toBe("4px");
  });

  it("dispatches band-height update on drag (Δy → onResize)", () => {
    const onResize = vi.fn();
    const { getByTestId } = render(
      <MemberMetaGutter /* …two-member fixture… */ onResize={onResize} />,
    );
    const gap = getByTestId("member-resize-gap");
    fireEvent.pointerDown(gap, { clientX: 0, clientY: 100, pointerId: 1 });
    fireEvent.pointerMove(gap, { clientX: 0, clientY: 110, pointerId: 1 });
    fireEvent.pointerUp(gap,   { clientX: 0, clientY: 110, pointerId: 1 });
    expect(onResize).toHaveBeenCalled();
    // Last call delivers the cumulative dy = 10
    const last = onResize.mock.calls.at(-1);
    expect(last?.[0]).toMatchObject({ dy: 10 });
  });

  it("tints the band on hover via data-active-band on the plot overlay", () => {
    const { getByTestId, getByTestIdQueryFn } = render(
      <>
        <MemberMetaGutter /* …two-member fixture, gap above member 7… */ />
        <div data-testid="member-trace" data-member-id="7" />
      </>,
    );
    const gap = getByTestId("member-resize-gap");
    fireEvent.pointerEnter(gap);
    // Component publishes the active member-id; the plot overlay subscribes
    // and sets data-active-band. Use a single source of truth via a context
    // mock or assert directly on the published value.
    const overlay = getByTestId("member-trace");
    expect(overlay).toHaveAttribute("data-active-band", "7");
  });
});
```

**Why these assertions matter.** Tailwind class strings are forbidden
(CLAUDE.md: "never assert on Tailwind class strings"); `getComputedStyle`
is unusable under JSDOM. The `data-state` + inline-style pattern lets
the test pin behaviour without a layout pass. The "tints the band"
case verifies the cross-component invariant from spec §7.4 that
otherwise has no test coverage.

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

Replace the `BandResizeDivider`'s 1px hairline with a `<div>` that
renders 4px (rest) / 12px (hover) tall, between member rows. The cursor
is `row-resize` while hovering. Wire pointer events to dispatch the
existing band-height adjustment logic (currently in
`MemberMetaGutter.tsx`).

Also: while a gap is in active hover OR being dragged, set a
`data-active-band={memberIdAbove}` attribute on the plot's per-band
overlays so the bands tint accent on hover (visual coupling). The plot
overlays already exist with `data-testid="member-trace"` and
`data-member-id`; the gutter just needs to publish the active member-id
via a context or callback.

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/MemberMetaGutter.tsx packages/HimalayaUI/frontend/src/components/BandResizeDivider.tsx packages/HimalayaUI/frontend/test/MemberMetaGutter.resize.test.tsx
git commit -m "feat(gutter): visible inter-row resize gap with band tint (Compare UX E-4)"
```

### Task E-5: Drop indicator + plot mirror on reorder drag

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/MemberMetaGutter.tsx`
- Test: `packages/HimalayaUI/frontend/test/MemberMetaGutter.reorder.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
describe("Drop indicator + plot mirror — Compare UX E-5", () => {
  it("renders a 2px drop indicator at the hover position during reorder drag", () => {
    // dragOver between row A and row B → assert <div data-testid="drop-indicator">
    // present at the gap.
  });
  it("plot bands mirror the indicator", () => {
    // assert the plot's per-band overlay at the same y has data-drop-target="true"
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

In `MemberMetaGutter.tsx`, on `dragOver` events compute the nearest
inter-row position and render `<div data-testid="drop-indicator">` at
that y. Publish the same y to a context that the plot overlays subscribe
to so they tint accent on the matching band.

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/MemberMetaGutter.tsx packages/HimalayaUI/frontend/test/MemberMetaGutter.reorder.test.tsx
git commit -m "feat(gutter): drop indicator + plot mirror on reorder (Compare UX E-5)"
```

### Task E-6: Phase E regression sweep

**Files:** (none — runs the suites)

- [ ] **Step 1.** `npm run build`.
- [ ] **Step 2.** `npm test -- --run` clean.
- [ ] **Step 3.** Manual smoke (dev loop) — open `/compare/:id`, expand a row, reorder, resize a band, save. Capture any UI rough edges as follow-up issues.

---

## Phase F — Sidebar redesign

Phase goal: render the new row layout (title / phase summary / author + age,
right-edge action zone) consuming the Phase A projection. New empty states.

### Task F-1: Replace the placeholder string with phase summary

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ComparisonSidebar.tsx`
- Test: `packages/HimalayaUI/frontend/test/ComparisonSidebar.projection.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
describe("Sidebar row projection — Compare UX F-1", () => {
  it("renders phases · traces, author byline, relative age", async () => {
    // Mock useComparisons to return a row with new fields.
    const row: ComparisonSummary = {
      id: 1, title: "Cubic vs Hex",
      description: null, content_hash: "h",
      created_by: 5, created_at: null,
      updated_at: "2026-05-14T11:00:00Z",
      forked_from_id: null, forked_at_hash: null,
      view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
      last_event_at: "2026-05-14T11:00:00Z",
      author_username: "alice",
      member_count: 4,
      member_phases: ["Pn3m", "Hex", "Lam"],
      has_stale_members: false,
    };
    // render with mocked query data + currentUserId !== 5
    expect(screen.getByText(/Pn3m · Hex · Lam · 4 traces/)).toBeInTheDocument();
    expect(screen.getByText(/by alice/)).toBeInTheDocument();
    expect(screen.getByText(/edited 1h ago/)).toBeInTheDocument();
  });

  it("renders 'by you' when current user is the author", async () => {
    // currentUserId === created_by → assert "by you"
  });

  it("renders ⚠ when has_stale_members is true", async () => {
    // assert <span data-testid="sidebar-stale-warn"> present
  });

  it("renders '+N more' for more than 3 phases", async () => {
    // member_phases length 5 → "phase1 · phase2 · phase3 · +2 more · 5 traces"
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

Replace the row body inside `ComparisonSidebar.tsx`:

```tsx
function formatPhaseSummary(phases: string[], total: number): string {
  if (phases.length === 0) return `${total} traces`;
  if (phases.length <= 3) return `${phases.join(" · ")} · ${total} traces`;
  return `${phases.slice(0, 3).join(" · ")} · +${phases.length - 3} more · ${total} traces`;
}

// inside the .map((c) => ...) row:
const isMine = currentUserId !== undefined && c.created_by === currentUserId;
const byline = isMine
  ? "by you"
  : c.author_username !== null ? `by ${c.author_username}` : "by —";
const rel = relativeTime(c.last_event_at, Date.now());

return (
  <li
    key={c.id}
    data-active={c.id === activeComparisonId ? "true" : undefined}
    data-interactable="button"
    onClick={() => onPickRow(c.id)}
    className="..."
  >
    <div className="font-medium truncate">{c.title}</div>
    <div className="text-xs text-fg-dim truncate">
      {formatPhaseSummary(c.member_phases, c.member_count)}
    </div>
    <div className="text-xs text-fg-dim truncate">
      {byline} · {rel === null ? "—" : `edited ${rel}`}
    </div>
    {c.has_stale_members && (
      <span data-testid="sidebar-stale-warn" className="absolute top-1 right-8 text-warning">⚠</span>
    )}
    {/* Pin star + ⋯ overflow — existing pin button stays; add ⋯ overflow on the right edge */}
  </li>
);
```

- [ ] **Step 4: Migrate the client-side sort key from `updated_at` to `last_event_at`.**

Spec §8.4 (L506-508): "Pinned-first … then `last_event_at` desc
(**replaces `updated_at` desc** — `last_event_at` covers edits and chat
activity)."

A-3's projection now ships `last_event_at` on every listing row, but
the sidebar's client-side sort is still keyed on `updated_at`. Update
`ComparisonSidebar.tsx` (or wherever `useComparisons()` results are
ordered before render) to:

1. Group by pinned status (pinned-first), then
2. Within each group sort by `last_event_at` desc (string-compare ISO
   timestamps; SQLite emits them in the canonical `YYYY-MM-DDTHH:MM:SSZ`
   format so lexicographic = chronological).

Add a failing test FIRST:

```tsx
describe("Sidebar sort key — Compare UX F-1 / spec §8.4", () => {
  it("sorts pinned-first, then by last_event_at desc (not updated_at)", () => {
    const rows: ComparisonSummary[] = [
      // … three rows: pinned A (last_event_at:T2, updated_at:T0),
      //   unpinned B (last_event_at:T3, updated_at:T1),
      //   unpinned C (last_event_at:T1, updated_at:T3) …
    ];
    // Mock useComparisons → rows; pin A only.
    // Render <ComparisonSidebar />.
    const items = screen.getAllByRole("listitem");
    // Expected order: A (pinned), B (last_event_at T3), C (last_event_at T1).
    // If sort were still updated_at-based: A, C (T3), B (T1) — wrong.
    expect(items[0]).toHaveTextContent("A");
    expect(items[1]).toHaveTextContent("B");
    expect(items[2]).toHaveTextContent("C");
  });
});
```

Then implement the sort change. The pinned check uses the existing
`useComparisonPins` hook. The whole rule is:

```tsx
const pins = useComparisonPins();
const sorted = useMemo(() => {
  return [...rows].sort((a, b) => {
    const ap = pins.has(a.id) ? 0 : 1;
    const bp = pins.has(b.id) ? 0 : 1;
    if (ap !== bp) return ap - bp;
    // Both pinned or both unpinned — sort by last_event_at desc.
    // Nulls land LAST: a row with no events is conceptually older than
    // any timestamped row, so it sorts after.
    // (Naive `bt.localeCompare(at)` with `?? ""` puts nulls FIRST under
    // descending sort because `""` lexicographically precedes any ISO
    // timestamp. Branch explicitly to invert that.)
    const at = a.last_event_at;
    const bt = b.last_event_at;
    if (at === null && bt === null) return 0;
    if (at === null) return 1;   // a is null → after b
    if (bt === null) return -1;  // b is null → after a
    return bt.localeCompare(at);
  });
}, [rows, pins]);
```

Without this step, the projection ships `last_event_at` but the UI
still orders on `updated_at` — chat activity moves a row to nowhere.

- [ ] **Step 5: Run.** Expected: PASS.

- [ ] **Step 6: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/ComparisonSidebar.tsx packages/HimalayaUI/frontend/test/ComparisonSidebar.projection.test.tsx
git commit -m "feat(sidebar): phase summary + author + age + stale-warn + last_event_at sort (Compare UX F-1)"
```

### Task F-2: Draft indicator on sidebar rows

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ComparisonSidebar.tsx`
- Test: `packages/HimalayaUI/frontend/test/ComparisonSidebar.draft.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
describe("Sidebar draft indicator — Compare UX F-2", () => {
  it("renders • + (draft) suffix on the row matching active draft id", () => {
    // useAppState.setState({ activeDraft: { id: 5, ... } })
    // assert text contains "(draft)" on the row for id=5
    // assert <span data-testid="sidebar-draft-dot"> present
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

Read `useAppState((s) => s.activeDraft)`. For each row, if
`activeDraft?.id === c.id`, append `(draft)` to the title text and render
the `•` dot prefix. Change byline to "by you · just now".

```tsx
const draftId = useAppState((s) => s.activeDraft?.id);
// ...
const isDraft = draftId === c.id;
const title = isDraft ? `${c.title} (draft)` : c.title;
```

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/ComparisonSidebar.tsx packages/HimalayaUI/frontend/test/ComparisonSidebar.draft.test.tsx
git commit -m "feat(sidebar): draft indicator (Compare UX F-2)"
```

### Task F-3: Empty state copy

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ComparisonSidebar.tsx`
- Test: `packages/HimalayaUI/frontend/test/ComparisonSidebar.empty.test.tsx`

- [ ] **Step 1: Failing test.**

```tsx
describe("Sidebar empty states — Compare UX F-3", () => {
  it("shows the no-comparisons-yet state + + New button", () => {
    // useComparisons → []
    expect(screen.getByText(/No comparisons in this experiment yet/)).toBeInTheDocument();
    expect(screen.getByTestId("sidebar-empty-new")).toBeInTheDocument();
  });
  it("shows search-empty state + clear button", () => {
    // useComparisons → [row]
    // user types in search → no matches
    expect(screen.getByText(/No matches/)).toBeInTheDocument();
    expect(screen.getByTestId("sidebar-empty-clear")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run.** Expected: fail.

- [ ] **Step 3: Implement.**

Two branches in the rendering logic:

```tsx
if (rows.length === 0) {
  return (
    <div className="px-4 py-8 text-center text-fg-muted">
      <p className="mb-3 text-sm">No comparisons in this experiment yet.</p>
      <button data-testid="sidebar-empty-new" onClick={onNew} className="...">
        + New comparison
      </button>
    </div>
  );
}
if (filtered.length === 0 && search.trim() !== "") {
  return (
    <div className="px-4 py-8 text-center text-fg-muted">
      <p className="mb-3 text-sm">No matches for "{search}".</p>
      <button data-testid="sidebar-empty-clear" onClick={() => setSearch("")} className="...">
        Clear search
      </button>
    </div>
  );
}
```

- [ ] **Step 4: Run.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/components/ComparisonSidebar.tsx packages/HimalayaUI/frontend/test/ComparisonSidebar.empty.test.tsx
git commit -m "feat(sidebar): empty state copy (Compare UX F-3)"
```

### Task F-4: Phase F regression sweep + e2e

**Files:** (none — runs suites + new spec)

- [ ] **Step 1.** `npm run build` clean.
- [ ] **Step 2.** `npm test -- --run` clean.
- [ ] **Step 3.** Add a Playwright mock spec.

```bash
touch packages/HimalayaUI/frontend/e2e/compare-sidebar-projection.spec.ts
```

Content:

```ts
import { test, expect } from "@playwright/test";

test("sidebar shows phase summary and author + age", async ({ page }) => {
  // Mock /api/experiments/.../comparisons to return a row with the new shape.
  // ... see e2e/AGENTS.md for patterns
  // assert phase summary + byline + age render.
});
```

```bash
cd packages/HimalayaUI/frontend && npm run e2e -- --grep "sidebar shows phase summary"
```

Expected: PASS.

- [ ] **Step 4. Commit.**

```bash
git add packages/HimalayaUI/frontend/e2e/compare-sidebar-projection.spec.ts
git commit -m "test(e2e): sidebar projection spec (Compare UX F-4)"
```

---

## Final wrap-up

### Task FIN-1: Run full suite

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-final.out 2>&1
grep -E "did not pass|fail" /tmp/jl-final.out | head -20

cd packages/HimalayaUI/frontend && npm test -- --run > /tmp/vitest-final.out 2>&1
grep -E "Tests:|Failed|passed" /tmp/vitest-final.out | tail -3

cd packages/HimalayaUI/frontend && npm run build 2>&1 | tail -10

cd packages/HimalayaUI/frontend && npm run e2e > /tmp/e2e-final.out 2>&1
grep -E "passed|failed" /tmp/e2e-final.out | tail -5
```

All four must be clean.

### Task FIN-2: Spec sync note

Open `docs/superpowers/specs/2026-05-14-compare-ux-refinement-design.md`.
At the very top of the doc, add a one-line note:

```markdown
**Implementation status:** Phase 1 (heavy-lifting) shipped per
`docs/superpowers/plans/2026-05-14-compare-ux-refinement.md`. Phase 2
items deferred — see the plan's "Deliberate scope split" preamble for
the list.
```

Commit:

```bash
git add docs/superpowers/specs/2026-05-14-compare-ux-refinement-design.md
git commit -m "docs(spec): note Phase 1 ship + Phase 2 deferral"
```

### Task FIN-3: Open the PR

Branch is already `worktree-compare-ux-refinement-spec`. Push and open
the PR:

```bash
git push -u origin worktree-compare-ux-refinement-spec
gh pr create --title "Compare UX refinement Phase 1: single-mode page + sidebar projection" --body "$(cat <<'EOF'
## Summary
- Drop /compare/:id/edit route; everything happens at /compare/:id with an inline-edit gesture surface.
- Single Save pill replaces the Save/Cancel/Discard triplet; copy adapts to author (Save changes) vs non-author (Save as fork…).
- Three-tier header (title strip / status surface / toolbar) replaces today's overpacked review header. EditOrForkButton, LineageBadge, NeedsReviewBadge, ForksPopover folded into the new surfaces (files deleted).
- Member rows collapse/expand with a right-edge action zone (⋯ + ⋮⋮ cue); grab-anywhere drag with 4px threshold; visible inter-row resize gap with band-tint on hover.
- Sidebar rows show phase summary + author + relative age + ⚠ for stale + draft indicator, consuming a new listing projection (author_username, member_count, member_phases, has_stale_members, last_event_at). Sort key migrates from `updated_at` to `last_event_at` (chat activity moves a row).
- View choices (grouping mode, peak ticks/labels) persist on the comparison via three new nullable columns and round-trip through GET / submit / SSE.
- Backend `comparison_submitted` / `comparison_created` SSE frames now carry `post_state = fetch_comparison_with_members(db, id)`; frontend `applyRemoteToCache` switches from invalidate to splice, preserving server-computed fields like `forked_from_title` while threading view_* through to the cache without a refetch round-trip.

## Out of scope (Phase 2)
- Right-rail Chat ↔ + Add traces tabs (spec §5).
- External add-trace drag flow + drop-anywhere (spec §7.3, §7.5).
- Color swatch row in expanded member rows.
- Cmd+Z undo (spec §3).
- Description inline-edit ergonomics (C-8 emits the editable surface;
  full keyboard / autosave behaviour deferred).
- "Reset to author's view" right-click affordance on the toolbar
  (spec §6.4 viewer escape hatch — surfaces in Phase 2 alongside the
  color swatch).
- Sidebar row `⋯` overflow menu (spec §4 / §8.2). F-1 leaves the
  existing pin button in place; the `⋯` joins in Phase 2 when the
  Forks-from-row gesture lands.
- Picker scope on `/compare/all` (spec §14 open question).
- Inline-modal replacement for `window.prompt` in C-14 (functional via
  `window.prompt` in Phase 1; an inline modal with `useFocusTrap` is
  the long-term path).
- Actor attribution on the `viewing-stale` banner (spec §6.2 updated
  to drop `@user` in Phase 1; `useStaleAgainstHash` derives the
  banner from cache-hash drift, not from SSE event payloads. Adding
  actor would require threading the SSE `actor` field through a
  separate signal).
- Queue-interleave coverage driving `replayCoordinator.handleRemoteEvent`
  directly for the comparison whole-row overwrite (A-9 Step 3d ships
  sequential-write coverage only; genuine replay-as-rerun coverage
  lives in `replayCoordinator.test.ts`).

See `docs/superpowers/plans/2026-05-14-compare-ux-refinement.md`
"Deliberate scope split" preamble.

## Test plan
- [ ] Julia suite green (`packages/HimalayaUI` tests)
- [ ] Vitest green (frontend)
- [ ] Playwright mocked e2e green
- [ ] Manual: open existing comparison → edit title → drag-reorder → save changes
- [ ] Manual: open someone else's comparison → edit → save → fork prompt → fork lands with correct lineage

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## Spec ↔ Plan coverage check (self-review)

| Spec section | Covered by | Notes |
|---|---|---|
| §1 Motivation | (descriptive; no task needed) | |
| §2 Goals | A–F | All five goals satisfied: one-mode page (B+C), single Save pill (C-7+C-13+C-14), interaction vocab (C-8/C-9/C-10/E-1/E-2/E-3/F-1, applied inline), informative sidebar (F-1), view-choice round-trip (A-1/A-4/A-5/A-7/A-10/C-4). |
| §2 Non-goals | (descriptive; no task needed) | Preserved by Phase 1 not touching `lib/queue/`, `comparison_messages`, mutation queue framework, or `comparison_pins`. |
| §3 Routing & mode model | B-1, B-2, C-5, C-6, C-11, C-14, C-15 (Step 4b SSE→viewing-stale wiring w/ acknowledge + own-op-rebase, Step 4c baseHashCaptureMoment test) | `/edit` redirect, comparePath option drop, useCompareMode (all six tagged-union variants), useCompareDraftDirty, page merge, setDraftForkOf mode transition. SavePill flow (C-7) instantiates the mode-aware copy. C-15 Step 4b's `useStaleAgainstHash` uses `useState` (not `useRef`) and rebases on both user-acknowledge AND own-op completion to prevent banner flicker on save. Step 4c pins `baseHash` capture-on-first-edit per spec §3 L96-97. |
| §4 Visual language | C-3, E-1, plus `data-interactable` attributes added inline in every new component | Right-edge action zone primitive in E-1. Sidebar `⋯` overflow → Phase 2. |
| §5 Right rail Chat / Add tabs | **Phase 2** (deferred). Interim chip ships with C-15. C-15 Step 4's `+ Add traces` button wires to the interim `ComparisonPickerPanel` chip, not a Phase-2 tab system. | |
| §6 Header tiers | C-8 (title strip; `readOnly` prop wired in C-12 so review mode doesn't drop edits silently), C-9 (status surface w/ `ServerUpdate` interface = `{previousHash, onAcknowledge}` + Saved-pill useEffect lifecycle), C-10/C-17 (toolbar + ⋯-more dropdown with `useFocusTrap`), C-12 (review-side replace), C-13 (edit-side replace). | View choices §6.4 → A-1/A-4/A-5 (incl. Step 5b post_state broadcast threading w/ `nothing`-handling + Step 5c HTTP body = `fetch_comparison_with_members`)/A-7/A-10/C-4 (incl. Step 0 `effectiveGroupingMode` helper + Zustand `groupingMode` deletion across 6 files). "Reset to author's view" right-click → Phase 2. |
| §7 Plot center | §7.1 empty state — C-15 Step 4 (failing-test scaffolding for empty block + `+ Add traces` wired to interim picker chip + drop target). §7.2 collapse/expand — E-2. §7.3 drag — E-3 (reorder only, explicit JSDOM-pointer-event test scaffolding). §7.4 resize — E-4 (`data-state` + inline style assertions, no `getComputedStyle`). §7.5 external drop — **Phase 2** (deferred). | |
| §8 Sidebar | F-1 (projection consumed + Step 4 sort-key migration `updated_at` → `last_event_at` per spec §8.4, nulls-land-last branch corrected), F-2 (draft indicator), F-3 (empty states) | Pin star stays as today via `useComparisonPins`; `⋯` overflow → Phase 2. |
| §9 Migration | A-1, A-2, B-1, B-2 | Idempotency tested in A-2. |
| §10 Testing | Tests authored per phase | Six-layer coverage: A-3/A-4 (route emit + `test_route_response_shapes.jl` `assert_keys` extension Steps 5b/4b — covers ALL THREE per-id sites incl. 409 conflict body), A-6 (SSE + Step 2b drain-and-count second-frame assertion), A-9 (applyRemoteToCache splice + `forked_from_title` survival + Step 3d replay-ordering interleave), A-8 (cache shape with value-time view_* assertions), A-10 (onMutate request body + Step 4b onSuccess reconciliation + Step 4c registry resolution), A-7 (compile-time response shape via TS types), A-5 Step 5c regression test (submit + create response body shape = `fetch_comparison_with_members`). Plus `forks_of_comparison` projection contract test (A-3 Steps 3b/3c). |
| §11 Build sequence | Phase order A → B → C → E → F | One phase per shippable unit. |
| §12 Risks | C-1 (threshold), C-5 (state machine), A-8 (cache shape), A-9 (splice vs invalidate, replay-ordering), C-4 Step 0 (`groupingMode` rewire blast radius). | Inline-edit on Esc resolved in C-3. |
| §13 Files touched | File-structure inventory at top | |
| §14 Open questions | Picker scope on /compare/all — Phase 2. Orphan `has_stale_members` — A-3 SQL handles via `cm.exposure_id IS NOT NULL` predicate. Esc-on-page (does NOT discard) — handled implicitly by C-3 (Esc inside `InlineEditableText` blurs without discard); page-level handler intentionally absent. | |

---

## Quick reference — gotchas pulled from CLAUDE.md / AGENTS.md hierarchy

- **Worktree.** All work in `.claude/worktrees/compare-ux-refinement-spec`. Don't `cd` to the original repo.
- **Test slowness.** Julia suite is 5–10 min. Capture to file once, grep the file. Don't re-run with different greps.
- **`Tables.rowtable` NULL semantics.** `ismissing(row.field)` not `=== nothing`. See A-3 / A-4.
- **`exactOptionalPropertyTypes: true`.** Optional fields declared as `T | undefined` (not `T?`) — see ActiveDraft in C-4.
- **Zustand named actions.** Don't `useAppState.setState(...)` directly outside tests; add named actions instead. See C-14's `setDraftForkOf`.
- **TanStack Query keys.** Don't fabricate keys; reuse `queryKeys.*` from `queries.ts`.
- **`apply_event!(InTransaction(), …)` inside `with_idempotency`.** Never the public `apply_event!(db, req)` — that opens a nested savepoint and broadcasts immediately. See A-5.
- **`data-interactable` attribute** on every new interactive surface for E2E selectors (per spec §4). Don't assert on Tailwind classes.
- **Singleton Oxygen routes.** `@get "/path/{id}" function(req::HTTP.Request, id::Int) ... end`.
- **`AGENTS.md` hierarchy.** Read the AGENTS.md nearest the file you're touching. The CLAUDE.md change on 2026-05-14 extracted module-specific rules into `src/AGENTS.md`, `packages/HimalayaUI/src/AGENTS.md`, `packages/HimalayaUI/frontend/src/AGENTS.md`, etc.
