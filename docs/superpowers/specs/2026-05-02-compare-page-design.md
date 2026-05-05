# Compare Page — Design Spec

**Date:** 2026-05-02 (revised 2026-05-05 for mutation-queue + per-tab `client_id`; revised 2026-05-05 for draft/submit model)
**Status:** Draft (revising as we go)

## Context

HimalayaUI has two functioning pages today: **Index** (curate one exposure's peaks and phase) and **Inspect** (curate which exposures of a sample to keep). The third page in the sidebar is **Compare**, currently a stub. The use case is the figure scientists actually publish: pick several traces — from one sample or across samples or even across experiments — overlay them on a shared `q` axis with adjustable y-offsets, and discuss what the differences mean.

Two distinguishing properties from the existing pages:

- **Many traces, not one.** Index and Inspect are both single-exposure-focused. Compare's primary unit is a *set* of exposures with per-member render parameters (offset, normalization, color, label, q-window).
- **Persisted artefact, not transient state.** A comparison is a thing a user creates, names, returns to, discusses with collaborators, and cites in chat.

Three architectural pieces shipped after the original web-app design spec was written, and a v1 Compare must compose with all three from day one:

- **Event log + dispatcher** ([docs/event-log.md](../../event-log.md)) — `update_view_for_event!` is the sole writer to view tables.
- **Mutation queue + idempotency** ([docs/mutation-queue.md](../../mutation-queue.md)) — every mutating route wraps in `with_idempotency`; events inside use the in-tx variant.
- **Six-layer contract testing** ([docs/contract-testing.md](../../contract-testing.md)) — adding any new `OpKind` requires the full test matrix.

## Mental model: comparisons are like GitHub PRs

A comparison is a **user-curated artefact**, not a collaborative workspace. The closest analogy is a GitHub pull request:

- One author drafts the comparison locally, iterating on members / offsets / normalization until they're happy.
- They submit it. Submission is the publish boundary — before submission, no one else can see the draft; after submission, it's a discussable artefact.
- The submitted comparison opens with a discussion thread. Other users comment on the figure, the author may iterate further (a new "edit pass" → re-submit), and chat references resolve against the submitted state.
- Edit mode and review mode have different interface surfaces. Edit mode is rich with member controls; review mode is plot + chat, no controls.

Implications:

- **Atomic submission is the only write path.** Within edit mode, member changes are pure client state. Submission is one mutation: one `client_op_id`, one idempotency boundary, one SSE frame, one `content_hash` change.
- **Live collaboration on member parameters is intentionally out of scope.** Two users can each fork their own variant of a comparison and discuss them in chat. They cannot co-edit offsets in real time. This is acceptable because the workflow is figure-crafting (one person, minutes-to-hours), not collaborative whiteboarding.
- **Per-field event log goes away.** `user_actions` records `comparison_submitted` with a full snapshot payload, not "Alice changed offset of trace 3." Acceptable trade because individual offset tweaks are not the kind of evidence we mine for parameter justification (unlike peak edits).
- **Conflict resolution is whole-snapshot, and only same-author cross-tab.** Author-only edit (see "Authorship and forking") means cross-user conflicts can't happen by construction. The remaining case: an author with the same comparison open in two edit-mode tabs submits from both. The second submit detects the stale baseline via `expected_content_hash` and shows a coarse "your other tab already changed this — overwrite / discard / fork" choice.

## Goals

- Persisted, named, navigable comparisons with discussion threads.
- Multi-trace overlay rendered against a shared `q` axis, with per-member offset / normalization / color / label / q-window.
- A clean **edit mode** (rich member controls, draft-only changes) vs. **review mode** (plot + chat, read-only) split.
- Atomic submission with idempotent retries and conflict detection via `expected_content_hash`.
- Mention-target parity with peak / index / exposure / sample (`@comparison:42`).
- Cross-experiment membership: a comparison can include exposures from any experiment in the central DB.

## Non-goals

- Live collaborative editing of member parameters (see "Mental model" above).
- Per-field audit history of edit-mode tweaks. Only submission boundaries are durable.
- Cross-experiment **peak/index** comparison (we are comparing rendered traces, not curation state).
- Statistical analysis of compared traces (cross-correlation, alignment scoring, residuals). Compare is a figure builder.
- Folder/tag organization. Flat list + search + pin until evidence demands more.
- Server-side draft persistence (drafts live in client `sessionStorage` — survives reload, scoped to one tab; cross-device drafts are deferred).
- Merge-back from forks (PR-style merge). Forks are independent artifacts in v1; if the original author wants changes from a fork, they manually copy.
- Export (PDF / SVG / PNG). Browser print is acceptable for v1.

## Pain points the design must avoid

- **Context drift.** A chat message saying "see the offset on JC001-120" must mean the same thing when re-read a week later. The submission boundary anchors this — chat references resolve against a `content_hash` that only changes on submit.
- **Membership churn.** Deleting an exposure that is a member of a published comparison must not silently rewrite the figure or break chat history.
- **Selection ergonomics at N=8.** Picking eight traces via a tree drill-down is painful. Picking via `@`-mention is fluent for two but awkward for eight. The picker modal (see "Selection UX") is what addresses this.
- **Lost drafts.** A user spending 20 minutes refining a comparison cannot lose it to a tab refresh. Hence `sessionStorage` mirroring of the draft state.

## Schema

Two view tables, both written exclusively by `update_view_for_event!` per the dispatcher contract.

```sql
CREATE TABLE comparisons (
  id              INTEGER PRIMARY KEY AUTOINCREMENT,      -- mention-target rule
  title           TEXT NOT NULL,
  description     TEXT,
  content_hash    TEXT NOT NULL,                          -- updated atomically on submit
  created_by      INTEGER REFERENCES users(id) ON DELETE SET NULL,
  created_at      TEXT NOT NULL,
  updated_at      TEXT NOT NULL,                          -- = last submission time
  forked_from_id  INTEGER REFERENCES comparisons(id) ON DELETE SET NULL,  -- lineage pointer
  forked_at_hash  TEXT                                    -- parent's content_hash at fork time
);

CREATE INDEX idx_comparisons_forked_from ON comparisons(forked_from_id);

CREATE TABLE comparison_members (
  id              INTEGER PRIMARY KEY AUTOINCREMENT,
  comparison_id   INTEGER NOT NULL REFERENCES comparisons(id) ON DELETE CASCADE,
  exposure_id     INTEGER REFERENCES exposures(id) ON DELETE SET NULL,  -- placeholder on delete
  display_order   INTEGER NOT NULL,
  band_height     REAL    NOT NULL DEFAULT 1.0,            -- vertical band ratio; see "Plot rendering"
  y_offset        REAL    NOT NULL DEFAULT 0,              -- baseline nudge within the band
  normalization   TEXT    NOT NULL DEFAULT 'none',         -- 'none' | 'max' | 'area' | 'qwindow'
  color_override  TEXT,                                    -- nullable; phase palette by default
  label_override  TEXT,                                    -- nullable; falls back to exposure name
  q_window_min    REAL,
  q_window_max    REAL,
  peak_display    TEXT,                                    -- JSON; per-peak display state, see "Plot rendering"
  created_by      INTEGER REFERENCES users(id) ON DELETE SET NULL,  -- drives picker "recently used"
  created_at      TEXT    NOT NULL                         -- set by dispatcher on insert
);

CREATE INDEX idx_comparison_members_by_comparison ON comparison_members(comparison_id, display_order);
CREATE INDEX idx_events_by_comparison ON user_actions(entity_id) WHERE entity_type = 'comparison';

CREATE TABLE comparison_messages (
  id            INTEGER PRIMARY KEY,
  comparison_id INTEGER NOT NULL REFERENCES comparisons(id) ON DELETE CASCADE,
  author_id     INTEGER REFERENCES users(id) ON DELETE SET NULL,
  body          TEXT NOT NULL,
  created_at    DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_comparison_messages_comparison ON comparison_messages(comparison_id, created_at);
```

Notes:

- `comparisons.id` is `AUTOINCREMENT` because comparisons are mention targets, same rule as `experiments`, `samples`, `exposures`, `peaks`, `indices` in CLAUDE.md.
- `comparison_members.exposure_id` is `ON DELETE SET NULL` (not `CASCADE`) so an exposure deletion leaves a visible "missing exposure" placeholder in the figure rather than silently mutating it. Chat references stay intact.
- `created_by` is `ON DELETE SET NULL` per the user-FK rule.
- **`comparisons` has no `experiment_id` column.** A comparison is a pure artifact — its scope is implied by its members, not a stored field. Sidebar listing and per-experiment views are derived from membership at query time (see "Sidebar listing" below). Empty comparisons are not persisted: a comparison only enters the DB when its first submission carries at least one member.
- **`forked_from_id` is `ON DELETE SET NULL`** — when a parent is deleted, forks survive as independent artifacts but lose their lineage pointer. Preserves the fork's chat history and content; doesn't cascade-delete derivative work. The reverse direction is also fine: deleting a fork has no effect on its parent.
- **`forked_at_hash` is immutable** — captures the parent's `content_hash` at fork time. Enables a future "parent has changed since you forked" indicator without requiring it in v1.
- `idx_events_by_comparison` is the parallel of `idx_events_by_exposure` — required for `rebuild_views_from_log!(db, comparison_id)` to perform.
- **`comparison_messages` is parallel to the existing `sample_messages`** — same shape, different scope. Chat history cascades on comparison delete (chat is part of the artifact's discussion, not independently durable). Generalizing both into a single polymorphic `messages` table would be cleaner long-term but is out of scope for this v1; doing it as part of Compare's plan would be feature creep.

**Cross-experiment membership.** Members may belong to any experiment in the central DB. A comparison "appears in" experiment X's sidebar if any of its members reference an exposure under X — the comparison can appear in multiple experiments' sidebars simultaneously, which is the correct cross-cutting view.

## Event kinds

Three view-producing kinds. All flow through `apply_event!(InTransaction(), …)` inside a `with_idempotency` body.

| Kind                   | Dispatcher action                                                                                  |
|------------------------|----------------------------------------------------------------------------------------------------|
| `comparison_created`   | `INSERT INTO comparisons(...)` + `INSERT INTO comparison_members(...)` for each member in payload. Computes initial `content_hash`. `entity_type='comparison'`, `entity_id` is the new id. |
| `comparison_submitted` | Diff payload snapshot vs current `comparison_members` for the comparison; `DELETE` rows no longer present, `INSERT` new ones, `UPDATE` changed ones. Recomputes `content_hash` + bumps `updated_at`. May also update `title` / `description`. |
| `comparison_deleted`   | `DELETE FROM comparisons WHERE id=entity_id`. Cascades to `comparison_members`.                    |

Plus one log-only kind:

| Kind                   | Purpose                                                                                            |
|------------------------|----------------------------------------------------------------------------------------------------|
| `add_message` (extended) | Existing kind gains `entity_type='comparison'` semantics. Comparison chat threads use a new `comparison_messages` table parallel to the existing `sample_messages` (which stays sample-scoped). The `add_message` event handler routes to whichever table based on the event's `entity_type`. |

That's it. No `comparison_member_added`, `_removed`, `_updated`, `_reordered`, or `comparison_renamed` — title, description, and full membership snapshot all ride inside `comparison_submitted` (and inside the initial `comparison_created` for new comparisons).

### Submission diff

`comparison_submitted` carries the full client-side draft snapshot. The dispatcher computes the diff against the current view-table state inside the transaction:

```julia
function update_view_for_comparison_submitted!(db, entity_id, payload)
    new_members = Dict(m.id => m for m in payload.members if m.id !== nothing)
    new_anonymous = [m for m in payload.members if m.id === nothing]
    current_ids = Set(member_ids_for_comparison(db, entity_id))

    # Delete members not in payload
    for id in setdiff(current_ids, keys(new_members))
        DBInterface.execute(db, "DELETE FROM comparison_members WHERE id=?", (id,))
    end
    # Update existing members
    for (id, m) in new_members
        if id in current_ids
            DBInterface.execute(db, "UPDATE comparison_members SET ... WHERE id=?", (..., id))
        end
    end
    # Insert new members (no id in payload)
    for m in new_anonymous
        DBInterface.execute(db, "INSERT INTO comparison_members(...) VALUES (...)", ...)
    end
    # Update comparison metadata + content_hash
    new_hash = compute_content_hash(db, entity_id)
    DBInterface.execute(db, "UPDATE comparisons SET title=?, description=?, content_hash=?, updated_at=? WHERE id=?",
                        (payload.title, payload.description, new_hash, now_iso(), entity_id))
    return nothing
end
```

The route around this body checks `expected_content_hash` against the pre-update value and returns 409 on mismatch.

### Conflict detection

`POST /api/comparisons/:id/submit` accepts an `expected_content_hash` field in the body. The route, inside the `with_idempotency` body but before the dispatcher call:

```julia
current = current_content_hash(db, comparison_id)
if !isnothing(payload.expected_content_hash) && current != payload.expected_content_hash
    return HTTP.Response(409, JSON3.write((;
        error = "conflict",
        current_hash = current,
        current_state = fetch_comparison_with_members(db, comparison_id),
    )))
end
```

`expected_content_hash = nothing` means "I'm not checking" (used by the very first submission of a freshly-created comparison, where there is no prior state). Subsequent submissions always include it.

Client receives 409 → shows the conflict modal (see "Conflict modal" under Frontend below). No automatic merge; the user picks.

## Mutation queue integration

One queue mutator, not seven. The skeleton:

```ts
// lib/queue/mutators/saveComparison.ts
const saveComparison: Mutator<SaveComparisonInput, ComparisonScope, ComparisonResponse> = {
  kind: "comparison_save",
  onMutate(payload, qc) {
    // No optimistic write. Submission shows a spinner, not an
    // optimistic membership swap — the user has been seeing the
    // draft locally throughout edit mode.
    return { restore: () => {} };
  },
  request(payload, signal) {
    // POST /api/comparisons/:id/submit  (or POST /api/comparisons for create — id absent means create)
  },
  onSuccess(payload, response, qc) {
    qc.setQueryData(queryKeys.comparison(response.id), response);
    qc.setQueryData(queryKeys.comparisonMembers(response.id), response.members);
    // Invalidate every per-experiment comparisons listing — comparisons are
    // membership-derived, so the response.members can touch any experiment.
    qc.invalidateQueries({ queryKey: ["comparisons"], exact: false });
  },
};
```

Plus a thin `deleteComparison` mutator. That's the queue surface for v1.

Things that *don't* apply because of the draft/submit model:

- **No negative placeholder ids.** Members aren't sent until submit; the cache reflects either the saved state or the local draft (Zustand), never an in-flight optimistic row.
- **No `useComparisonHasPendingOps` granularity.** A comparison either has a pending `comparison_save` or it doesn't; the existing `useMutationState` hook is enough.
- **No per-field `applyRemoteToCache` branches.** A foreign `comparison_submitted` event invalidates `comparison(id)` + `comparisonMembers(id)` and refetches. Member-level field merges aren't a thing.

What still applies:

- **`with_idempotency` wraps the submit route.** Retries return the cached response without re-executing.
- **`apply_event!(InTransaction(), …)` inside.** SSE frame is post-commit-broadcast.
- **Six-layer contract tests** for `comparison_save` and `comparison_delete` only — that's 2 kinds × 6 = 12 cases, not 7 × 6 = 42.

## Hash memoization

`comparisons.content_hash` is SHA-256 over a canonical serialization of:

- `(title, description)`
- members ordered by `display_order`, each as `(exposure_id, band_height, y_offset, normalization, color_override, label_override, q_window_min, q_window_max, peak_display)`

Computed by the dispatcher at the end of any view-mutating event.

**Citation behavior.** Because `content_hash` only changes at submission boundaries, chat references are simple:

- `@comparison:42` resolves to "the current state."
- `@comparison:42@<hash8>` resolves to a frozen snapshot. If `@hash8` differs from the current `content_hash`, the chip shows a small "(changed)" indicator and a hover tooltip with a diff summary.
- The chat compose path **eagerly** writes `@comparison:42@<currenthash8>` at insertion time. Eager resolution is unambiguous in the draft/submit model because content hash only changes at submission boundaries — there are no intermediate hashes.

## Authorship and forking

A comparison has one author (`comparisons.created_by`). Authorship is the gating concept for editing; non-authors cannot edit, but they can fork.

### Author-only edit

The submit and delete routes check `comparisons.created_by == user_id_for(req)` and return `403 Forbidden` if the requesting user isn't the author. Frontend hides the Edit button (and the Delete affordance) for non-author viewers, replacing them with a Fork button.

This is enforced under the existing `X-Username` trust model — the same auditing-but-not-authentication regime that governs `created_by` everywhere else in the codebase. When real authentication lands, the route gate becomes strictly enforced rather than label-checked; no spec change needed.

**Orphaned comparisons.** A comparison whose `created_by` is `NULL` (because the author's user row was deleted) cannot be edited by anyone. It becomes a read-only artifact, fork-only.

### Forking

Anyone (including the author) can fork any comparison. A fork is structurally a new comparison created via the existing `POST /api/comparisons` route, with the payload carrying:

```json
{
  "title": "Copy of <original title>",
  "description": "<inherited from parent>",
  "members": [...],                       // deep copy of parent's members at fork time
  "forked_from_id": <parent_id>,
  "forked_at_hash": "<parent's current content_hash>"
}
```

The dispatcher writes a `comparison_created` event with these fields — no new event kind, no new contract test row. The new comparison's `created_by` is the forker, `id` is fresh, and `content_hash` is recomputed (initially identical to the parent's if nothing was changed before submission, distinct after any edit).

After fork: the user lands directly in **edit mode** of the new comparison. Forking implies "I want to change something" — landing in review mode of an unchanged copy would be a useless intermediate.

**Lineage display.** A forked comparison shows a badge in the review-mode header:

```
Forked from "Pn3m → Im3m transition" by Alice  [view parent →]
```

The parent's title hyperlinks to its review page. If the parent is deleted (`forked_from_id IS NULL`), the badge becomes:

```
Forked from a deleted comparison
```

A reverse-direction "Forks (N) →" link near the badge opens a popover listing the comparison's own forks. Cheap (the `idx_comparisons_forked_from` index makes the query fast).

**No merge.** Fork is one-directional. Bob's edits to a fork of Alice's comparison stay on Bob's fork; there's no PR-style merge-back in v1. If Alice wants to incorporate Bob's changes, she manually copies the relevant edits.

**Cross-cutting workflow.** Conflict modal's "Copy to new" affordance becomes "Fork" — it's the same action, now generalized as the v1 escape valve for any "I want to make changes but the original isn't mine" or "I want to branch my own work" scenario.

## Multiplayer

What stays multiplayer:

- **Chat thread.** `add_message` events broadcast normally. Two users (or any number) can discuss a comparison live regardless of who authored it. Comments are gated only by user identity, not by authorship.
- **Submission visibility.** When the author submits, all other viewers' review-mode pages invalidate and refetch via SSE → the new figure appears.
- **Per-tab `client_id` (PR #32).** Two tabs of the same user have distinct `client_id`s; the originating tab's edit refreshes the other.

What does not become multiplayer:

- **Edit-mode draft state.** Strictly local (Zustand + sessionStorage). Because edit is author-only, only one user is ever in edit mode at a time per comparison.
- **Co-editing.** Two users cannot edit the same comparison concurrently (author-only edit makes this structurally impossible). Cross-user collaboration happens via fork.

**Conflict cases that remain:**

- **Same-author cross-tab.** Alice has tab A in edit mode and tab B in edit mode (or saves from one, then submits from the other). The 409 path with `expected_content_hash` handles this: the second submit sees a stale baseline and the conflict modal opens.

R5b's deferred field-level optimistic concurrency is *not* needed — whole-snapshot conflict is much easier to present, and author-only edit narrows the cases that hit it.

## REST API

```
GET    /api/experiments/:eid/comparisons                  -- list comparisons that touch this experiment (membership-derived)
GET    /api/comparisons                                   -- list all comparisons (cross-experiment)
POST   /api/comparisons                                   -- create + initial submit; body must include >=1 member
                                                          -- optionally carries forked_from_id + forked_at_hash for forks
GET    /api/comparisons/:id                               -- fetch + members + lineage info
GET    /api/comparisons/:id/forks                         -- list forks of this comparison (cheap via index)
POST   /api/comparisons/:id/submit                        -- author-only (403 otherwise); idempotent; 409 on hash mismatch
DELETE /api/comparisons/:id                               -- author-only (403 otherwise); idempotent
GET    /api/comparisons/:id/messages                      -- chat thread (any user)
POST   /api/comparisons/:id/messages                      -- any user; add_message with entity_type='comparison'; idempotent
```

No member-level routes — all member changes ride inside `POST /submit`. The create route is *not* nested under an experiment because comparisons aren't owned by experiments; it's a top-level resource.

`GET /api/experiments/:eid/comparisons` is a derived listing:

```sql
SELECT DISTINCT c.*
FROM comparisons c
JOIN comparison_members cm ON cm.comparison_id = c.id
JOIN exposures e ON e.id = cm.exposure_id
JOIN samples s ON s.id = e.sample_id
WHERE s.experiment_id = :eid
ORDER BY ...;
```

Sorted by `MAX(created_at) FROM user_actions WHERE entity_type='comparison' AND entity_id = c.id` (latest event time).

## Frontend

### Mode boundary

A comparison page is in one of two modes:

- **Review mode (default).** Plot + chat + read-only header. No member controls. Mention previews, hover details, plot toolbar (log/linear, peak overlays). An "Edit" button enters edit mode.
- **Edit mode.** Plot + member panel + add-member affordance + title/description inputs. No chat — comments would compose against state nobody else can see. "Save" / "Cancel" / "Discard draft" buttons in a sticky footer. The plot updates live from local state.

The boundary is enforced at the route level: `/compare/:id` is review mode; `/compare/:id/edit` is edit mode (URL state survives reload, plus `sessionStorage` carries the draft body). Cancel returns to `/compare/:id`. Save submits, then returns to review mode on success.

For the **initial create** flow, the user lands directly in edit mode for a draft that exists only client-side (Zustand + sessionStorage). The first Save with at least one member = `comparison_created` and the comparison enters the DB. Cancelling never-submitted creation discards. An empty draft never produces a server row.

### Sidebar listing

The sidebar list of comparisons is **membership-derived**, not stored. When viewing experiment X, the sidebar shows comparisons whose members touch any exposure in X — meaning a comparison spanning experiments A and B appears in both sidebars simultaneously. This matches the cross-cutting nature of comparisons: a figure that involves experiment X's data is genuinely relevant whether the user is in X or in another experiment that shares some of the same traces.

A header toggle in the sidebar switches between scoped and global views:

- **This experiment** (default): comparisons that touch the current experiment.
- **All experiments**: every comparison in the central DB.

Sort by latest event (`MAX(created_at)` from `user_actions` for the comparison's entity stream). Search/pin within whichever scope is active.

### Page layout

Compare mirrors the Index page's three-card composition for visual continuity. Sidebar (collapsible comparison list, scope-toggle in header) + main column (plot panel + chat). The chat card has the same shape, height behavior, and component (`ChatCard`) as on Index — users should not have to relearn where it is or how it works.

```
┌── Sidebar ──┬── Plot Panel (wide card) ─────────────────┐
│ comparisons │  ┌─ plot ─┐ ┌── per-trace metadata ────┐  │
│ list        │  │  thin  │ │  trace 1: phase, a, R²,K │  │
│ + New       │  │  tall  │ │  trace 2: …               │  │
│ search      │  │  3:10  │ │  trace 3: …               │  │
│ pin         │  │ aspect │ │  …                        │  │
│             │  └────────┘ └───────────────────────────┘  │
│             ├── Chat Card ───────────────────────────────┤
│             │  same shape as Index/Inspect              │
└─────────────┴────────────────────────────────────────────┘
```

The plot panel is a single card. Inside it, two regions:

1. **Plot column.** Thin tall plot, fixed aspect ratio `0.3` (3:10 W:H) for v1 — hardcoded based on the workflow not benefiting from on-the-fly tuning. The `aspect` prop is plumbed in the component for future surfacing if feedback shows it's needed. With 8 traces and a 600px-wide panel, the plot column is ~180px wide × ~600px tall, giving each trace a vertical band of ~75px.
2. **Per-trace metadata gutter.** Wider than the plot. Each trace's metadata row is absolute-positioned so its vertical center aligns with the trace's y-baseline in the plot — one source of truth for the y-stacking math, both halves stay in lockstep.

The metadata row is a single line per trace (matches the thin-plot aesthetic): trace label, phase chip, lattice parameter `a`, R², and `K` for cubic phases. Hover or click expands into a detail card overlay. Number-of-peaks is intentionally *not* in this row — by the time a comparison exists, the indexing decision is already made; surfacing peak count here would re-litigate something settled. Peak count remains accessible via the hover/expand path for completeness.

**Reflow.** On narrow viewports the sidebar collapses behind a button; the plot panel stays intact and the metadata gutter wraps below the plot column past a breakpoint where the gutter would otherwise compress below readable. Chat moves below the plot panel, as on other pages. Same reflow story as Index — no new mental model.

**Mode differences in the same layout footprint.**

- *Review mode:* metadata gutter is read-only; chat is mounted; annotation toggles in header.
- *Edit mode:* metadata gutter rows expose drag handles on the left and per-row controls (offset, normalization, color, label, q-window) inline or in an expand-on-click sub-panel; chat slot is replaced by the member-add affordance and submit/cancel buttons; annotation toggles are hidden (everything is shown so the user can edit).

### Plot rendering

One shared Observable `<Plot>` instance with a single q-scale and a single zoom domain — so brush + double-click reset zoom uniformly across all traces. Marks are organized as a per-member render function: each member contributes its line + peak ticks + labels into the band of the y-axis allocated to it. Zustand draft updates trigger a re-mark of just that member's contribution.

Conceptually:

```ts
function MultiTraceMarks({ members, draft, panelHeight }) {
  const totalRatio = sum(members.map(m => m.band_height));
  let cumulative = 0;
  return members.map((m) => {
    const yBand = [
      panelHeight * (cumulative / totalRatio),
      panelHeight * ((cumulative + m.band_height) / totalRatio),
    ];
    cumulative += m.band_height;
    return (
      <MemberTraceLayer
        key={m.id}
        member={m}
        yBand={yBand}
        peakDisplay={draft.peakDisplayFor(m.id)}
        onPeakClick={editMode ? cyclePeakDisplay : undefined}
      />
    );
  });
}
```

The y-stacking math (which member occupies which band) is shared with `MemberMetaRow` in the gutter so reorders stay aligned automatically.

#### Band heights

Each member's vertical band is sized by `band_height` (default `1.0`). Heights are *ratios*, not absolute units — the rendered pixel height of band `i` is `panel_height × (band_height_i / Σ band_heights)`. Default all-`1.0` produces uniform stacking; bumping one member's `band_height` to `2.0` allocates it twice the room (useful for crowded peaks or very dynamic traces). The figure always fills the panel — there is no scrollable overflow regardless of how the heights are rebalanced.

Each member's baseline is computed from preceding bands:

```
y_baseline_i = panel_height × (Σ_{k < i} band_height_k / Σ band_height) + y_offset_i
```

`y_offset` is a fine-tuning nudge within the band (rarely used; defaults to `0`).

#### Working band vs total band

The total band (sized by `band_height`) is the *envelope* allocated to a member. Within it, peaks and annotations render into a smaller **working band** — a centered sub-region defined as a fraction of the total (`working_band_fraction`, hardcoded at `0.7` for v1; not yet a schema field, defer until tunability is asked for).

```
total band:    [────────────────────────────]
                 ↑15%               ↑15%
                  headroom          headroom
working band:    [────────────────]    ← peaks + annotations live here
```

**Why this matters in SAXS.** A direct-beam dropoff at low `q` can be orders of magnitude above the diffraction peaks, even after qwindow-style exclusion. Without headroom, the dropoff would either crush the peaks (full-band normalization with no exclusion) or visually overlap the working region of the next trace up (working-band-fills-total). The headroom band absorbs the tail without invading the next trace's working band — adjacent peaks never collide regardless of how dramatic the dropoff is.

**Tail clipping rule.** A trace's rendered curve is allowed to enter the full total band (working + headroom). Anything that would extend *beyond* the total band is clipped at the band envelope. Adjacent traces' working bands are guaranteed disjoint; only their headroom regions can host tail bleed, and even those don't intersect because the bands tile y without overlap.

**Working band as the normalization target.** The default normalization maps each trace's reference value to the top of the *working* band, not the total band. The trace's lower-q tail (if dramatically larger than the reference) extends naturally into the headroom and clips at the band edge if needed. This decouples "where peaks render" from "how dynamic the trace is."

The "reference value" is computed adaptively per member:

```
if member has detected peaks within q_window:
    reference = max(peak_intensity for peaks within q_window)    # peak-fit
else:
    reference = max(signal within q_window)                       # signal-fit
```

The two cases handle distinct workflows:

- *Peak-bearing traces* (typical sample exposures): peak-fit ensures the brightest peak hits the top of the working band, so peak-position comparison reads cleanly across the stack.
- *Peakless traces* (buffer, control, exploratory exposures with no indexable features): signal-fit ensures the trace's overall envelope fills the working band so the user can still see its shape — without this fallback, a buffer trace with weak features and no detected peaks would have nothing meaningful to normalize against and might collapse to a flat line.

The fallback is automatic; users don't pick between the two manually. Per-member `normalization` enum (`max` / `qwindow` / `none` / `area`) modulates the *strategy* but the peak-vs-signal selection inside the strategy is implicit.

#### Note on absolute intensity (SAXS context)

Absolute intensities are scientifically meaningful in SAXS more broadly — scattering strength carries information about concentration, contrast, and structure factor magnitude. Defaulting to peaks-fit-working-band normalization is a deliberate concession to **this tool's purpose**: Compare exists for peak-position-and-shape analysis across traces, where the absolute scale is secondary to alignment.

Users who need absolute-scale comparison have two escape valves:

- *Per-member `normalization = 'none'`* — render the trace at its raw scale within its band, clipped at the band envelope. Useful when the absolute intensity ratio between two traces is the question being asked.
- *Per-member `band_height` adjustments* — manually allocate more vertical room to a strong scatterer so its raw shape doesn't crash into the band ceiling.

This is not a limitation of the model; it's a calibration of the default to the typical workflow. If users frequently need to compare absolute intensities, that's evidence to surface a "preserve absolute scale" header toggle, but we shouldn't pre-build it.

#### Resize affordance — drag handles between metadata rows

In edit mode, a thin grabbable horizontal divider sits *between* each adjacent pair of metadata rows in the gutter. Cursor changes to row-resize on hover. Dragging the divider:

- **Push semantics in normalized space.** Drag down by `Δ` (in normalized ratio units): trace above grows by `Δ`, trace below shrinks by `Δ`. The total `Σ band_heights` stays constant, so the figure's pixel layout is stable — no reflow of unrelated traces.
- **Clamp at a minimum.** Each `band_height` is clamped at a small positive minimum (e.g. `0.1`) so a band can't be dragged out of existence. To remove a trace entirely, the user removes the member.
- **Snap to neighbors.** When the dragged divider passes within a small dead zone of "match the band on the other side," it snaps — so producing two equal-height bands stays easy.

A "**Reset heights**" button in the header zeros out customizations: every `band_height` returns to `1.0`. Useful for recovering uniform stacking after experimentation.

Dragging only resizes; reorder happens via the separate row-grip handle on the metadata row itself (see "Member reorder"). The two affordances live at different positions on the row — divider lives *between* rows, grip lives *on* the row — so the user gestures don't collide.

#### Existing concerns

Log/linear x toggle, auto-fit y-floor at `p01·0.5` from `PlotCard::computeFit` applied per-member to each trace's normalized signal scaled into its band. Q-window scoping (`q_window_min/max`) clips the rendered curve and is *separate* from the global zoom — per-member q-window is part of the saved comparison state, global zoom is per-tab UI.

In **edit mode**, the plot reads from local draft state. In **review mode**, from `queryKeys.comparisonMembers(id)`.

### Peak click semantics (edit mode)

Peaks render as small triangle markers on each trace. Clicking cycles the peak's display state for that member:

```
shown (default) → labeled → hidden → shown
```

`alt+click` jumps directly to hidden for users with many peaks. State persisted in `comparison_members.peak_display`:

```json
{
  "hidden": [<peak_id>, ...],
  "labeled": [<peak_id>, ...]
}
```

**Peak identity is stored as peak id, not q-value.** This is purely an internal tracking concern — peak ids never appear as visual labels on the plot. Identity stays stable across re-fetches; if the underlying exposure is reanalyzed and a peak id disappears (the auto-vs-manual id-instability gotcha from Plan 7 docs), the rendering gracefully falls back to "shown" for that peak. No q-value snapping at the comparison layer; we trust upstream peak curation in the Index page and focus Compare's logic on rendering.

Peak ids surface visually only via:

- Hover tooltip on the triangle (q-value, optional id, prominence)
- A debug "show peak ids" toggle (developer-facing; not in v1 user UI)

### Label placement

Labeled peaks render their text **vertically above** the triangle marker. Leader lines connect text to triangle when collision-avoidance offsets the label horizontally — same "puppet string" idea as scatter-label dodging. Implementation: collect all labeled peaks in q order per trace, run a 1D dodge layout pass, emit text + leader-line marks.

Auto-rotation of label text (vertical to fit narrow plot column) is on; horizontal labels would frequently exceed the thin plot's width.

### Zoom and pan

Same control surface as the existing `TraceViewer`: brush to zoom q-window, double-click to reset, mouse-wheel to pan/zoom. Reuse the same hooks where possible. Note the semantic difference vs Index: on Compare the brush is *only* zoom — peak editing happens via clicks on triangles, never via a brush selection.

### Annotation toggles (review mode)

Header control with three checkboxes; all default on. Toggling produces the clean view used for figure capture or close inspection of the curves themselves:

- **Peak ticks** — triangles + labels.
- **Phase / Miller index ticks** — predicted-ratio overlays from `phaseratios`.
- **Per-trace labels** — the small label rendered on or near each trace baseline.

State is per-tab Zustand only; not persisted on the comparison. The comparison stores *what to show*; the toggles control whether to render the show-list right now. Hidden in edit mode (everything renders so the user can edit it).

### Member reorder

Drag handle on each metadata row, *not* on the plot. Plot interactions are reserved for peak clicks. Same affordance as Inspect's thumbnail drag. On drop, the draft's member ordering updates; both plot and gutter re-render aligned (shared y-stacking math).

### Trace coloring

Default: **one color per sample**, deterministic from `sample_id`. Computed at render time as `palette[hash(sample_id) % palette.length]` against a small categorical palette (~10–12 distinguishable hues, accessibility-checked). No schema work — colors are derived, not stored.

Implications:

- **Reordering doesn't shift colors.** Color is a property of the sample, not the slot.
- **Same-sample multiple traces share the color.** Acceptable for v1 — replicate-shading isn't supported elsewhere in Himalaya yet, so introducing it here would be ahead of the rest of the codebase. If two traces from the same sample become hard to distinguish, the user breaks the tie via per-member `color_override`.
- **No cross-experiment equivalence.** Two samples with the same name in different experiments have different `samples.id`s and get different colors. There's no canonical-sample concept in the codebase to lean on; pretending we could detect equivalence here would be inventing judgment the system doesn't have.
- **Orphan placeholders (deleted exposure → `exposure_id IS NULL`)** render in a neutral gray. The per-member `color_override` still wins if set.

#### Grouping mode toggle

Header control in both modes (review and edit) selects how default colors are computed:

```
Color: [By sample ▾]  → By sample (default) / By phase / Distinct per trace
```

- *By sample* — `palette[hash(sample_id) % len]`. v1 default.
- *By phase* — confirmed-phase color from [phases.ts](../../../packages/HimalayaUI/frontend/src/phases.ts); unindexed traces fall back to neutral gray.
- *Distinct per trace* — palette walked by `display_order`. Useful when every trace needs to be visually unique regardless of grouping.

State is **per-tab Zustand**, not persisted on the comparison — same as the annotation toggles. The comparison stores `color_override` per member (which always wins); the grouping mode controls only the default render for members without an override. Different viewers may want different views ("show me phases" vs "show me samples") of the same comparison without re-saving it.

#### Per-member override

The metadata row's expanded state in edit mode includes a color swatch grid (palette choices). Click to set; "Reset color" clears the override and reverts to the active grouping mode's computed default. Custom hex input deferred — palette-only for v1 keeps figures visually consistent within a comparison.

### Selection UX

Adding members is the most-used edit-mode action. The primary path is a **searchable picker modal**; the secondary path is a **warm-add accelerator** from the Inspect / Sample pages. `@`-mention is reserved for chat references and is *not* used to add members (keeping mention semantically pure: it's for talking *about* curated things, not for curating).

#### Picker modal (primary)

Triggered by the "Add traces" button in edit mode.

```
┌── Add traces to comparison ─────────────────────────────────┐
│  Search: [_______________________]                          │
│  Filters: [Experiment: this ▾] [Tag: ▾]                    │
│                                                             │
│  Recently used                                              │
│  ┌──────────────────────────────────────────────────────┐  │
│  │ ☐ JC-buffer-01    PBS pH 7.4, 25°C                  │  │
│  │ ☐ JC-lipid-DOPC   DOPC bicelles                      │  │
│  └──────────────────────────────────────────────────────┘  │
│                                                             │
│  All exposures                                              │
│  ┌──────────────────────────────────────────────────────┐  │
│  │ ☐ JC001-120   Sample A1                              │  │
│  │     50% DOPC + 50% chol, 10 mM CaCl₂, hydrated 24h  │  │
│  │ ☐ JC001-121   Sample A1                              │  │
│  │     50% DOPC + 50% chol, 10 mM CaCl₂, hydrated 48h  │  │
│  │ ☑ already added (locked)                             │  │
│  │ JC001-130   Sample A2                                │  │
│  │     50% POPC, control                                │  │
│  │ …                                                     │  │
│  └──────────────────────────────────────────────────────┘  │
│                          [Cancel] [Add 3 selected]          │
└────────────────────────────────────────────────────────────┘
```

**Row content.** Exposure name, sample name (`samples.label` if present else `samples.name`), sample notes (`samples.notes`, truncated to ~2 lines / ~120 chars; full text on hover). Pick-by-sample-identity is the actual mental model — scientists recognize traces by what the sample was, not by trace properties. Phase chips deliberately omitted from v1; trace-property metadata can be added later if it turns out to be missed in practice. (UI copy may say "description" rather than "notes" to be friendlier; the underlying column is `samples.notes`.)

**Filters.** Two only:

- *Experiment.* Multi-select; defaults to the experiment context the user was in when they opened the picker (the comparison page's current `:eid` route param, if any; otherwise unscoped). Cross-experiment picks require explicitly broadening this chip.
- *Tag.* Multi-select of `(key, value)` pairs from `sample_tags` (the actual schema is key-value, not flat tag names). Pairs are displayed as `key:value` chips in the dropdown. Multi-select uses OR semantics across selected pairs (any sample matching any selected pair passes the filter). Includes a synthetic "untagged" option (samples with zero rows in `sample_tags`). Hides itself entirely if no tags exist anywhere in scope.

No phase / status / R² / peak-count filters in v1. They re-litigate decisions already made; user is picking by sample knowledge, not by trace evaluation.

**Sort under "All exposures".** Alphabetical by exposure name. Predictable scrolling at scale; works without any temporal mental model.

**"Recently used" section.** Per-user, server-side history. The top of the modal surfaces the user's most recently picked exposures (across all comparisons, across all experiments) — explicit support for the long-tail workflow of repeated controls (buffer, lipid control, water, etc.) that get added to many comparisons.

**Already-added members.** Shown in the list as `☑ already added (locked)` — visible (so the user knows the trace is in the comparison) but not actionable. Prevents double-add.

**Pick order = append order.** Selecting `[120, 130, 045]` in that click order appends them to `display_order` in that order. The user reorders via drag in the metadata gutter afterward.

**Search.** Substring match across `exposures.name`, `samples.label`, `samples.name`, `samples.notes`, and `sample_tags.value`. Case-insensitive, no fuzzy matching for v1 — debounced SQL `LIKE '%query%'` is fine at lab scale.

**Empty state.** "No exposures match. Try clearing filters or broadening to all experiments." Direct nudge toward the most likely missing scope.

**Keyboard.** Arrow keys move selection cursor; space toggles the focused row's checkbox; Enter adds all selected; Esc cancels. Search input has focus on open; tab moves to the filter chips and then into the list.

**Component reuse.** The list rows share shape with how exposures appear elsewhere. Pull a `<ExposureListRow>` component used by both the Inspect page's exposure list and this picker so the visual identity is consistent and one set of changes propagates.

#### "Recently used" derivation

No separate `picker_recent` table. The `comparison_members.created_by` and `comparison_members.created_at` columns (already in the main schema above) are populated by the dispatcher on every member `INSERT` — both during `comparison_created` and during `comparison_submitted`'s diff-insert path. The "recently used by user X" query is then:

```sql
SELECT cm.exposure_id, MAX(cm.created_at) AS last_picked
FROM comparison_members cm
WHERE cm.created_by = :user_id
  AND cm.exposure_id IS NOT NULL
GROUP BY cm.exposure_id
ORDER BY last_picked DESC
LIMIT 20;
```

Crosses experiments naturally (no `experiment_id` filter). Picks the user removed but didn't re-add are absent — "recently used" shouldn't include things the user actively rejected.

#### Warm-add accelerator (secondary)

On the Inspect page, an "Add to comparison" affordance per exposure (or as a multi-select toolbar action) opens a small dropdown:

```
Add to:
  ○ Recent comparison: "Pn3m → Im3m transition"   (current draft)
  ○ Pick a comparison...   → opens the picker modal
  ○ + New comparison        → creates draft + adds the exposure(s)
```

This catches the "I was already in Inspect, I want to throw this trace into a comparison" workflow without forcing a context switch through the modal. Adding to "current draft" is silent — the comparison page reflects it next time the user navigates there; if they're already on Compare in another tab, SSE refresh handles it.

#### Future: starred picks

A "Starred" section above "Recently used" for deliberately-curated common picks. Same shape as recently-used, but explicitly chosen by the user. Not v1; flagged here so the picker layout reserves visual space for it.

### State model

- **Server state (TanStack Query):** comparisons list, comparison detail (saved state), members (saved), messages. Query keys: `queryKeys.comparisons(experimentId)`, `queryKeys.comparison(id)`, `queryKeys.comparisonMembers(id)`, `queryKeys.comparisonMessages(id)`.
- **Client state (Zustand):**
  - Edit-mode draft per comparison id (or a single "active draft" slot, since only one comparison can be in edit mode at a time per tab).
  - Sidebar expansion, drag-in-progress state.
  - Mode (read from URL, mirrored for cross-component reads).
  - All via named actions per CLAUDE.md frontend gotchas.
- **Persistence:** the draft slot is mirrored to `sessionStorage` with a schema version. Reload resumes editing in place. Closing the tab loses the draft (acceptable for v1; server-side drafts are a v2 concern).

### Mention integration

`@comparison:N` works the same as other mention types. `MentionPicker` gains a `Comparisons` group, scoped to the active experiment by default with a `@comparison:` prefix to widen to all experiments. Hover preview shows title + member count. Compose-time eager hash resolution (see "Hash memoization" above).

Mentions disabled or routed differently in edit-mode chat? — N/A, edit mode has no chat.

### Skeleton loading

Per the boneyard rules in CLAUDE.md, every load-gated card on Compare gets a `<Skeleton>` wrapper with `loading={query.isLoading}`, a `fixture`, a `fallback` mirroring HintText, and `flex-1 min-h-0 flex flex-col` on the wrapper.

### Infrastructure banner

Already global (`InfrastructureBanner` mounted at App.tsx). No Compare-specific work — `comparison_save` participates automatically since it goes through `useQueueMutation`.

### Conflict modal

When `POST /submit` returns 409, the client opens a modal. With author-only edit, the conflict path is now strictly *same-author cross-tab* (you can't collide with yourself across users since you're the only one with edit access):

> **You submitted a different version of this comparison from another tab.**
>
> *Side-by-side diff: tab A's version vs your current draft.*
>
> [Discard my changes] [Overwrite with mine] [Fork to a new comparison]

"Overwrite with mine" re-submits with `expected_content_hash` set to the server's current hash. "Fork to a new" creates a new comparison from the local draft with `forked_from_id` set to the original — the same fork mechanism described in "Authorship and forking" above.

## Open questions

1. ~~**Selection UX.**~~ Resolved: searchable picker modal as primary path (filters: experiment + tag; row content: name + sample name + sample description; per-user recently-used surfaced) + warm-add accelerator from Inspect/Sample pages. `@`-mention reserved for chat references. See "Selection UX" above.
2. ~~**Offset control.**~~ Resolved: per-member `band_height` (ratio, default `1.0`) sized via drag handles between metadata rows in the gutter. `y_offset` retained as fine-tune nudge within the band. See "Plot rendering — Band heights" and "Resize affordance" above.
3. ~~**Color/grouping default.**~~ Resolved: by sample (deterministic palette indexed by `sample_id`), with a per-tab grouping-mode toggle (by sample / by phase / distinct) and per-member `color_override`. No replicate-shading; same-sample multi-trace shares the color. See "Trace coloring" above.
4. ~~**Normalization defaults.**~~ Resolved: peaks-fit-working-band by default, with adaptive peak-fit-vs-signal-fit fallback (peak-bearing traces normalize on peak max in qwindow; peakless traces fall back to signal max in qwindow). Working band = inner ~70% of total band; headroom absorbs direct-beam dropoff without invading neighboring working bands. SAXS context note documents the tradeoff explicitly. See "Plot rendering — Working band vs total band" and "Note on absolute intensity" above.
5. ~~**Comparison home experiment when members span experiments.**~~ Resolved: comparisons have no `experiment_id` field. Sidebar listings are membership-derived — a comparison appears in every experiment's sidebar that has at least one member touching that experiment. Empty comparisons are never persisted; the server only sees a comparison after its first submission with ≥1 member. Cross-experiment "All comparisons" view available via `GET /api/comparisons`. See "Schema notes" and "REST API" above.
6. ~~**Versioning / forking.**~~ Resolved: author-only edit; anyone can fork via `POST /api/comparisons` with `forked_from_id` + `forked_at_hash`. Fork uses the existing `comparison_created` event kind (no new event surface). Lineage badge in review-mode header; reverse-direction "Forks (N)" popover. No merge-back in v1. See "Authorship and forking" above.
7. ~~**Discard-draft prompt on stale resume.**~~ Resolved: no separate prompt. With author-only edit, the only way for a draft to become stale is the same user submitting from a parallel tab — already handled by the 409 conflict modal. `sessionStorage` (not `localStorage`) means drafts don't survive tab close, so cross-session staleness doesn't exist. Just proceed; 409 handles the rare case.
8. ~~**Aspect-ratio control surfacing.**~~ Resolved: hardcode `0.3` (3:10 W:H) for v1, no header control. Author's experience suggests on-the-fly width tuning is rarely needed; without a workflow that demands it, the control would be UI clutter. Revisit after collecting user feedback — if multiple users hit datasets where the default reads poorly, surface a header control. The `aspect` prop stays plumbed in the component so the upgrade path is just adding the control, not refactoring.

## Implementation sequence (sketch)

1. Schema (incl. `peak_display`) + dispatcher branches for `comparison_created` / `comparison_submitted` / `comparison_deleted` + `rebuild_views_from_log!` extension + dispatcher round-trip test.
2. Backend routes wrapped in `with_idempotency`; conflict-detection path; `test_route_response_shapes.jl` rows; idempotency-replay invariant rows.
3. `saveComparison` + `deleteComparison` mutators; six-layer contract tests for each (12 cases total).
4. Frontend list page (with scope toggle: this experiment / all) + create flow (edit mode for client-only draft until first submit).
5. Picker modal (search + experiment/tag filters + recently-used + already-added locks); warm-add affordances from Inspect/Sample.
6. Multi-trace plot component: shared Plot instance with per-member mark layers, y-stacking math driven by `band_height` ratios, log/linear x toggle, q-window brush zoom, double-click reset.
7. Metadata gutter rows (label / phase / `a` / R² / `K`); drag-handle reorder *on* row, resize divider *between* rows; review-mode read-only vs edit-mode controls; "Reset heights" header button.
8. Edit-mode peak click cycle (shown → labeled → hidden → shown); label placement above triangle with leader-line dodge; `peak_display` persistence.
9. Review-mode annotation toggles (peak ticks / phase ticks / labels).
10. Review-mode chat thread + `@comparison` mention resolution + content_hash citation.
11. Author-only edit gates (route 403 + frontend hide-Edit-show-Fork); fork creation flow; lineage badge + forks popover.
12. Conflict modal (overwrite / discard / fork).
13. Sidebar polish, search, pin, skeleton coverage, draft sessionStorage persistence.

## See also

- [event-log.md](../../event-log.md) — dispatcher contract, hash memoization, SSE broadcast semantics.
- [mutation-queue.md](../../mutation-queue.md) — queue architecture, `client_op_id` lifecycle, deferred-promise pattern, `with_idempotency`.
- [contract-testing.md](../../contract-testing.md) — six-layer testing rule.
- [himalayaui-design.md](../../himalayaui-design.md) — design philosophy.
- [chat-mention-design.md](2026-04-29-chat-mention-design.md) — `@`-mention system that comparisons extend.
- [2026-05-01-multiplayer-instrumentation-design.md](2026-05-01-multiplayer-instrumentation-design.md) — Plan 7 design that this builds on.
