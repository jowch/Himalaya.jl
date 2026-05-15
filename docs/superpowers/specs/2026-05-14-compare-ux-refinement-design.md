# Compare UX refinement — design

**Status:** design (2026-05-14)
**Supersedes (UX layer only):** `docs/superpowers/specs/2026-05-02-compare-page-design.md`
— this refinement keeps the original data model, mutation queue contracts,
snapshot-at-submit semantics, and conflict resolution flow intact. It restructures
the page surface, the editing model, and the listing payload.

---

## 1. Motivation

The Compare page works but is hard to live in. Concrete pains today:

1. **Heavy mode flip.** `/compare/:id` (review) and `/compare/:id/edit` are
   separate React routes with different right-rail contents (ChatCard ↔ hint
   card) and different header controls. Adding a trace forces a route change;
   sharing a URL requires picking the right one.
2. **Sidebar rows are uninformative.** `ComparisonSidebar.tsx` hard-codes
   `"3 traces"` as a literal string. No phase summary, no member count, no
   author, no relative age. Picking the right comparison from a list is
   guess-and-click.
3. **Review header is overpacked.** Seven controls (`GroupingModeToggle`,
   `AnnotationToggles`, `NeedsReviewBadge`, `EditOrForkButton`, `LineageBadge`,
   `ForksPopover`, `FigureExportControls`) wrap into one row that mixes
   toggles, status badges, and rare navigation popovers, with no visual
   hierarchy.
4. **Two entry points to "add a trace."** Edit-mode's picker panel and
   Inspect's `WarmAddMenu` are different surfaces solving the same job.
5. **Author intent doesn't round-trip.** `groupingMode`, `showPeakTicks`,
   and `showPeakLabels` live in per-tab Zustand. A comparison crafted "By
   phase" opens "By sample" for the next viewer.
6. **No inline rename / no autosave.** Title is a giant text input in the
   plot center; Save / Cancel / Discard compete as three same-weight buttons.
7. **Member gutter affordances are non-obvious.** Drag handle is a small
   grip; band-resize dividers are 1px hairlines visible only on hover;
   per-member controls cram into a 280px column.

This spec refines the page around one organizing idea: **one calm page that
becomes editable in place, with consistent interactivity vocabulary.**

## 2. Goals and non-goals

**Goals**
- Eliminate the review ↔ edit route flip. One page, always in-place editable.
- One commit gesture (a "Save changes" pill) that appears only when there's
  something to save and adapts copy for author vs. non-author (fork).
- Make affordances discoverable through a consistent hover-reveal visual
  language across rows, the gutter, the sidebar, and the title bar.
- Make the sidebar rows informative enough to navigate without clicking
  blindly.
- Persist the author's view choices (grouping mode, peak ticks/labels) on
  the comparison so they round-trip to other viewers.

**Non-goals**
- No change to snapshot-at-submit semantics, the `with_idempotency` flow,
  the mutation queue framework, or the conflict-resolution UI.
- No change to the picker's internal data model — it's the same sample-first
  list, moved into a right-rail tab.
- No change to fork semantics. Non-author edits still produce a fork via
  `forked_from_id` / `forked_at_hash`; only the gesture surface changes.
- No new event kinds. No change to `comparison_messages` chat surface
  beyond moving it into a tab.
- No change to `comparison_pins`, the square-aspect plot constraint
  (issue #81), or the mutation queue framework files (`lib/queue/`).
- Not introducing autosave. The explicit Save pill remains the commit
  gesture; this keeps snapshot-at-submit cleanly aligned with one user
  intent.

## 3. Routing and mode model

**Routes after this change**
- `/experiments/:eid/compare` — sidebar with no comparison selected.
- `/experiments/:eid/compare/new` — fresh draft, no saved comparison.
- `/experiments/:eid/compare/:id` — one comparison, may have a draft.
- `/compare/all` and `/compare/all/:id` — global scope, same shapes.

The `/edit` segment is dropped. A `<Route path="*/edit" />` handler in the
router maps any legacy URL to its non-`/edit` counterpart via `<Navigate>`,
preserving `sessionStorage`-resident drafts on the way through.

**Mode model — there is no global "edit mode"**
- The page is always editable in place: title, description, member order,
  band height, normalization, per-member controls.
- A "Save changes" pill at the top-right of the plot center appears only
  when `dirty === true`. Until then the page is calm — no Save / Cancel /
  Discard triplet.
- **`baseHash` capture moment.** The first edit gesture on the page (title
  edit, member add, drag-reorder, etc.) captures the *current cached*
  `content_hash` into `draft.baseHash`. Before any edit, there is no draft.
  Foreign SSE updates that arrive *before* first edit refresh the cached
  hash (so the would-be capture tracks the latest server truth);
  foreign updates *after* first edit do NOT mutate `draft.baseHash`, and
  the 409 conflict modal is the recovery path. This subtly differs from
  the prior spec's "baseHash captured at edit-mode entry" rule — it has
  to, because there is no edit-mode entry anymore. A contract test
  (`Compare.baseHashCaptureMoment.test.tsx`) pins this.
- **Esc** blurs the active editor (does not discard).
- **Cmd+Z** undoes the last draft change locally (Zustand stores a tiny
  undo stack per draft; depth 20). Inside the inline title `<input>`, the
  browser's native undo wins (we don't intercept Cmd+Z while a text
  input owns focus).
- **Discard changes** lives inside a small popover on the Save pill, not as
  a permanent button.
- **Cross-tab / SSE coexistence.** While a draft is active, foreign
  `comparison_submitted` events for the same id surface an inline
  *"Server-side updated since you last viewed — save may conflict"*
  affordance in the status surface (§6). The viewer can acknowledge
  the drift (rebases the local "last seen" hash and clears the
  banner) or proceed to save normally, which surfaces the existing
  409 conflict modal on collision (see "Two fork triggers" below).
  Note: the actor name is intentionally NOT surfaced in Phase 1 —
  `useStaleAgainstHash` derives the banner from cache-hash drift, not
  from SSE event payloads; adding actor attribution requires threading
  the SSE `actor` field through a separate signal and is deferred to
  Phase 2.

**Author vs. non-author — two distinct fork triggers**

The Save pill expresses *intent*; the 409 conflict modal handles *recovery*.
These are intentionally separate.

| Trigger | When | UX |
|---|---|---|
| **Non-author intent fork** | Pill is clicked and the current user is not `comparison.created_by`. Pill copy is `Save as fork…`. | Inline prompt for the fork title (default `"Copy of <parent title>"`); confirm → morph draft to a create-with-lineage → POST `/api/comparisons`. Today's `EditOrForkButton` + `startFork` factory cover this flow today; we just move the gesture surface from a separate button to the pill. |
| **Conflict-recovery fork** | The author's submit hits a 409 because foreign edits landed after `draft.baseHash` froze. Existing `ConflictModal` fires (unchanged). | The modal offers Overwrite / Fork / Discard. The Fork branch reuses the same `startFork` factory, so the outcome is identical to the intent path — just reached after a server collision. |

Non-authors never see the `ConflictModal` on save, because they always
route through the intent-fork path before submit — there's no
`expected_content_hash` riding on the create call, so the server can't
409 them.

**View choices and `content_hash`.** `view_grouping_mode`,
`view_show_peak_ticks`, and `view_show_peak_labels` are **not** part of
`content_hash` — they are presentation preferences, not identity. They
ride in the submit payload (and persist server-side), but a viewer
toggling Annotations does NOT race against the author's save via 409.
The fields are owned by `SaveComparisonBody` (in `api.ts`) and forwarded
by `saveComparison`'s `buildBody`. `compute_content_hash` continues to
hash only title, description, and member identity / display fields.

**`useCompareMode` tagged union — full state space**

```ts
type CompareMode =
  // No draft. Read-only surface for author or non-author alike.
  | { kind: "viewing" }
  // No draft, but the cached comparison's content_hash has drifted from
  // the previously-read one due to a foreign SSE update. Same surface as
  // "viewing" but a status banner advertises the change. Becomes
  // "viewing" silently if the user clicks "reload" / re-fetches.
  | { kind: "viewing-stale"; previousHash: string }
  // Author has a draft against an existing comparison id.
  | { kind: "editing-mine"; draftId: number }
  // Non-author has a draft against an existing comparison id. On
  // SavePill click, transitions to "creating-from-fork" via
  // setDraftForkOf (which clears draft.id and sets forkedFromId).
  | { kind: "editing-as-fork-of"; parentId: number }
  // /compare/new with no fork lineage. Bare draft.
  | { kind: "creating-blank" }
  // /compare/new with a parent — either the post-morph state after
  // setDraftForkOf, or WarmAddMenu's "+ New comparison from this fork".
  | { kind: "creating-from-fork"; parentId: number };
```

**Transitions worth naming explicitly:**
- `editing-as-fork-of` → `creating-from-fork`: triggered by the SavePill
  click prompt. `setDraftForkOf({ newTitle, parentId, parentHash })`
  clears `draft.id` and `draft.baseHash`, sets `forkedFromId` /
  `forkedAtHash`, and updates the title. The next render sees `draft.id
  === undefined` and falls into `creating-from-fork`. The submit
  routes via the create path (no `expected_content_hash`), so the
  fork commits without a conflict possibility.
- `viewing` → `viewing-stale`: triggered by `applyRemoteToCache`
  handling a `comparison_submitted` for the open comparison when no
  draft exists. The hook detects `content_hash` drift.
- `viewing-deleted`: not in v1 scope. If the comparison is deleted by
  another user, the next GET will 404 and the page mounts the
  `StaleUrlPage` fallback.

This union is the single source for header copy, pill copy, and
permissions.

## 4. Visual language for interactivity

**Rule of thumb:** rest state is calm (no chrome). Affordances surface on
hover. The whole Compare surface adopts one consistent rest / hover / active
language.

| Interaction type | Rest | Hover | Active / Pressed |
|---|---|---|---|
| Static interactable (button, chip, link, tab) | `text-fg-dim`, no border | `text-fg`, `bg-bg-hover`, `border-border` | `bg-bg-subtle`, `text-fg`, accent ring or outline |
| Draggable row | `text-fg`, no chrome | `bg-bg-hover`, faint 1px ring around row, cursor: `grab` | cursor: `grabbing`, row lifted (small shadow), source row at 60% opacity |
| Inline-editable text (title, label, meta) | text only | underline-dotted, pencil ✎ at right edge, cursor: `text` | input chrome appears, Esc cancels / Enter commits |
| Click-to-expand row | text + glyph only | chevron ▸ fades in at left | chevron rotates ▾, row expands |
| Drop target (during drag) | invisible | n/a (pointer is captured) | 2px accent ring on target gap or band, faint accent bg tint |

All tokens already exist in `styles.css` `@theme` (`--color-bg-hover`,
`--color-bg-subtle`, `--color-border`, `--color-fg-dim`, etc.) — no new
design tokens are required.

**Right-edge action convention.** Action affordances (drag cue, overflow
menu, pin star) live on the trailing (right) edge of any row in a gutter
(member gutter and comparison sidebar both). Content (label, meta) flows
left-to-center; the leading edge stays free for hover-revealed chevrons.

**Selectors.** Each interactive surface carries `data-interactable={kind}`
on its outermost element (`button`, `drag`, `edit`, `expand`, `droptarget`)
so E2E and unit tests can assert on the type without coupling to Tailwind
class strings (per CLAUDE.md frontend gotcha).

## 5. Right rail — Chat / + Add traces tabs

The right rail in `WorkspaceGrid`'s right slot becomes a tabbed surface
with two tabs:

- **Chat** (the existing `ChatCard`, comparison-scoped, visible at all
  times — no longer hidden when editing).
- **+ Add traces** (the existing `ComparisonPickerBody`, sample-first
  search + tag filter + recently picked + checkbox per row).

**Default tab**
- No draft active → Chat.
- Draft active (user has started editing, or this is `/compare/new`) → +
  Add traces.

**Auto-switch rule.** Clicking any `+` empty-state button on the page
auto-flips the rail to + Add traces. Drag-and-drop targets (§7) accept
exposure drops whether or not + Add traces is the active tab; landing a
drop auto-flips the tab as a discovery aid.

**Picker tab content additions**
- A short "Drop exposures here" affordance at the top of the panel,
  advertising drag-and-drop from Inspect's thumbnail strip and from the
  sidebar.
- Footer line "Already added: N" clickable to scroll to the already-added
  section.

**WarmAddMenu unification.** Inspect's `WarmAddMenu` collapses to one
button (`+ Add to comparison`) with a thin dropdown:
- `+ New comparison` — discard existing draft, start a fresh one with the
  active exposure pre-added; navigate to `/experiments/:eid/compare/new`.
- `Add to "<active draft title>"` — visible only when a draft exists.
- Both paths land on `/compare/:id` with the right rail forced to + Add
  traces.

**State.** Right-rail tab is per-tab Zustand (`compareRightTab: "chat" |
"add"`), not URL-persisted. Persisted to `sessionStorage` (not the
localStorage-backed `himalaya-ui:state` Zustand persist channel) so it
survives reload but does not leak between tabs. Same persistence channel
as the existing draft (`himalaya-ui:compare-draft`).

## 6. Header, status surface, title strip — three-tier layout

The plot-center header collapses into three vertically stacked tiers:

```
┌─ Title strip ────────────────────────────────────────────────────────┐
│  Cubic vs Hex (replicate set)                ✎                       │
│  by Alex · edited 2h ago · 4 traces · forked from "Cubic phase ref"  │
└──────────────────────────────────────────────────────────────────────┘
┌─ Status surface (only when something to say) ────────────────────────┐
│  ⚠ Needs review — analysis changed since last submit. [Re-snapshot]  │
└──────────────────────────────────────────────────────────────────────┘
┌─ Toolbar ────────────────────────────────────────────────────────────┐
│  [Group: By sample ▾]  [Annotations ▾]   ⋯ more   ⤓ Export  [Save]   │
└──────────────────────────────────────────────────────────────────────┘
```

### 6.1 Title strip

- Title is rendered as text. Click flips it to an inline `<input>` styled
  to look like the heading. Esc cancels; Enter or blur commits to draft.
- A subtle pencil (✎) hover-glyph signals editability — no persistent
  bordered input box.
- The description sits below as muted text with the same single-click-to-edit
  pattern.
- **Meta row** under the title (single dim line, separated by `·`):
  `by <author> · edited <relative> · <N> traces`. If a fork, append
  `· forked from <parent title>` linked to the parent (this absorbs
  `LineageBadge.tsx`, which can be deleted).
- `created_by === currentUserId` renders as `by you` instead of the username.

### 6.2 Status surface

- One stacked column under the title strip; shows transient state. Collapses
  to 0px when empty.
- Banner types: `Needs review` (warning amber, matching `StaleIndicesBanner`'s
  visual vocabulary), `Server-side updated since you last viewed — save may
  conflict` (info; actor attribution deferred to Phase 2 — see §3
  cross-tab note), `Saved` (success, auto-dismisses after 2s).
- The Re-snapshot action lives inside the Needs Review banner. Uses the
  same call site as today's `NeedsReviewBadge.tsx::re-snapshot`.

### 6.3 Toolbar

- **Left:** `[Group ▾]` and `[Annotations ▾]` — view controls, rendered as
  dropdown buttons (compact, less visual noise than today's pill rows).
- **Center:** `⋯ more` opens a small dropdown containing rare actions:
  `Forks (N)` (today's `ForksPopover` content; the file can be deleted or
  shrunk to a list renderer), `Copy link`, `Delete`, and `Discard changes`
  (only when dirty).
- **Right:** `⤓ Export` (today's `FigureExportControls`) and the **Save
  pill** (only when dirty).

### 6.4 Persisted view choices

Three nullable columns are added to `comparisons`:

| Column | Type | Notes |
|---|---|---|
| `view_grouping_mode` | TEXT NULL | `"bySample"`, `"byPhase"`, `"distinct"`, or NULL |
| `view_show_peak_ticks` | INTEGER NULL | 0/1 or NULL |
| `view_show_peak_labels` | INTEGER NULL | 0/1 or NULL |

The frontend reads these as defaults; a user override on the page is a
draft change (counts toward dirty), and commits on save. NULL on existing
comparisons; frontend falls back to per-tab Zustand defaults when NULL.

**Non-author view-choice behaviour.** A non-author toggling Group ▾ /
Annotations ▾ on someone else's comparison creates a draft (per §3, any
edit gesture captures `baseHash` and seeds a draft). On SavePill click,
copy is `Save as fork…` — view-choice changes are part of the fork
payload, not an in-place mutation of the parent. This is intentional:
viewers cannot mutate the author's persisted view choices without
explicit fork intent. **Goal 5 is therefore preserved for authors and
not unintentionally broken for viewers** — viewers either save-as-fork
(preserving their view choices on the new comparison) or use the
right-click "Reset to author's view" affordance to discard the local
draft and snap back to the persisted defaults.

No new event kind. The view-choice fields ride in the existing
`comparison_submitted` payload and the `comparison_created` payload
(both create and submit paths INSERT/UPDATE the three columns).

**`content_hash` does not include view-choice fields** (see §3 "View
choices and `content_hash`"). This keeps cross-tab Annotation toggles
from racing against the author's save via 409.

## 7. Plot center — empty state, member rows, drag-and-drop

### 7.1 Empty state

```
┌────────────────────────────────────────────────────┐
│                                                    │
│         No traces yet.                             │
│                                                    │
│         [+ Add traces]                             │
│                                                    │
│         Or drag exposures from a sample list       │
│         into this area to add them.                │
│                                                    │
└────────────────────────────────────────────────────┘
```

- The `+ Add traces` button auto-flips the right rail to + Add traces.
- The empty-state region is itself a drop target for exposure drags.
- No `MemberMetaGutter` while empty — full center width.

### 7.2 Member rows — collapsed and expanded

```
COLLAPSED (default)
┌────────────────────────────────────────────────────────┐
│    Sample-A · E257 (overridden label)        ⋯  ⋮⋮     │  ← right-edge actions
│    Pn3m · R² 0.99 · q [0.05–0.35]                      │
└────────────────────────────────────────────────────────┘

EXPANDED (click row body to toggle)
┌────────────────────────────────────────────────────────┐
│ ▾  Sample-A · E257                           ⋯  ⋮⋮     │
│    Pn3m · R² 0.99                                      │
│   ─────────────────────────────────────────────────    │
│   Label   [Sample-A · E257                       ]     │  inline-editable
│   Color   ● ● ● ● ● ● ● ●  ↩ reset                     │  swatch row
│   Norm    [Max ▾]                                      │
│   q-win   [0.05] – [0.35]    [Reset]                   │
│   Peaks   ☑ ticks   ☑ labels   ☐ hide noise            │
└────────────────────────────────────────────────────────┘
```

- **Collapsed view** shows: label (inline-editable), one-line meta
  (phase, R², q-window), and the right-edge action zone (`⋯` overflow,
  `⋮⋮` drag cue).
- **Click handling.** A click on the row body (anywhere outside the
  right-edge actions and the inline label) toggles expansion. A click on
  the label text enters inline-edit instead of toggling. A click on the
  right-edge actions (`⋯`, `⋮⋮`) does what those affordances do. **One row
  expanded at a time by default** to keep vertical alignment with the plot
  bands predictable.
- The expanded form-controls region is not draggable (typing/picking
  shouldn't surprise-drag). Only the top strip (label + meta) of the
  expanded row remains a drag target.
- A chevron `▸/▾` fades in on the left edge on hover; rotates when
  expanded.

### 7.3 Drag, reorder, and the threshold guard

The Compare page has **two distinct drag flows** that must coexist:

| Flow | Source | Drop target | Effect | Code path |
|---|---|---|---|---|
| **Reorder** | An existing member row (collapsed body or expanded top strip) | An inter-row gap in the same gutter | `setMemberOrder([...])` on the draft | Reuses today's HTML5 drag-and-drop wiring in `MemberMetaGutter.tsx`; `dataTransfer` carries the member id |
| **External add** | An exposure in any list outside Compare (Inspect thumbnail strip, picker row, sidebar sample expansion) | Plot center, gutter empty space, or an inter-row gap | Adds the exposure to the draft as a new member, at the y where the drop landed (or appended) | New `dataTransfer` payload format `{ kind: "exposure", id }`; the gutter and plot accept both payload kinds and dispatch by `kind` |

Behaviour shared by both flows:

- **Grab anywhere.** The collapsed row body is the drag target for reorder;
  the `⋮⋮` glyph on the right edge is signage, not a hit target.
- **Drag-vs-click disambiguation.** `pointerdown → pointermove > 4px →
  drag begins`. `pointerup < 4px` runs the click handler (toggle expand,
  or commit inline-edit if the down was on label text).
- **Drag-vs-text-selection in the label.** The same threshold handles it:
  a click on the label enters inline-edit. Text-selection inside a closed-row
  label by mouse-drag is intentionally not supported — click to edit, select
  inside the input. Matches Notion/Linear row patterns.
- **Drop indicator.** A 2px accent line in the gutter at the drop position;
  the plot's band stripes briefly mirror the indicator at the same y.
  Same visual treatment for both flows.

### 7.4 Band resize

- The hairline `BandResizeDivider` is replaced by a visible **gap between
  rows** (4px → 12px on hover) that carries a `row-resize` cursor.
- Hovering the gap highlights both adjacent bands in the plot (the visual
  coupling makes "I'm resizing the boundary between A and B" unambiguous).
- The existing `computeYBands` math is unchanged; the inter-band padding
  constant moves from implicit 0 to a small explicit value.

### 7.5 Drop-anywhere add

- Exposures in any list (Inspect thumbnail strip, picker list, sidebar
  sample expansion) are draggable.
- The plot center and the gutter accept exposure drops while + Add traces
  is open or not. Dropping onto the center auto-opens + Add traces as a
  discovery aid; the drop succeeds regardless.
- A drop is equivalent to clicking the picker checkbox: add to draft,
  positioned at the y where the drop landed or appended.

## 8. Sidebar redesign and listing projection

### 8.1 New listing payload fields

Added to the JSON returned by `GET /api/comparisons` and
`GET /api/experiments/:id/comparisons`. Both go through `comparisons_listing` /
`comparisons_for_experiment` in `comparisons.jl` and the shared
`_comparison_listing_rows` projector.

| Field | Source | Notes |
|---|---|---|
| `author_username: string \| null` | `LEFT JOIN users ON users.id = comparisons.created_by` | Null when author deleted (FK SET NULL) |
| `member_count: int` | `(SELECT COUNT(*) FROM comparison_members WHERE comparison_id = c.id)` | Bounded by N |
| `member_phases: string[]` | `json_extract(comparison_members.snapshot, '$.confirmed_index.phase')`, aggregated | Up to 3 unique phases, ordered by frequency desc, ties broken by **first-seen `display_order` asc** (deterministic). Members without a confirmed index contribute nothing. Rest summarized as `+N more` client-side. |
| `has_stale_members: bool` | `EXISTS(... snapshot_inputs_hash != e.analysis_inputs_hash)` | Same logic that drives the per-member `is_stale` flag |
| `last_event_at: string \| null` | Already computed as a server sort key — now projected | Frontend uses it for the "edited X ago" display |

Performance for an experiment with K comparisons is K row reads plus K
member-aggregate subqueries — well inside the 50ms budget on WAL-mode SQLite
for typical K (≤ few hundred). If this regresses, promote to a denormalized
`comparisons.member_phase_summary TEXT` column refreshed at save (no caller
API change).

### 8.2 Sidebar row layout

```
┌─ ComparisonSidebar ─────────────────────────┐
│  [This experiment] [All]            + New   │
│  ┌──────────────────────────────────────┐   │
│  │ Search comparisons                   │   │
│  └──────────────────────────────────────┘   │
│  ─────────────────────────────────────────  │
│ ┃ Cubic vs Hex                    ★  ⋯     │  ← active (accent left-border)
│ ┃ Pn3m · Hex · Lam · 4 traces              │
│ ┃ by alex · edited 2h ago                  │
│  ─────────────────────────────────────────  │
│   Replicate set                  ⚠  ☆  ⋯   │  ← stale; ☆ hover-only
│   Pn3m · Pn3m · 2 traces                    │
│   by jc · edited 1d ago                     │
│  ─────────────────────────────────────────  │
│   Cubic study (draft)            •  ☆  ⋯   │  ← unsaved draft indicator
│   Pn3m · +2 more · 5 traces                 │
│   by you · just now                         │
└─────────────────────────────────────────────┘
```

- **Line 1:** title (single-line truncate, unchanged).
- **Line 2:** `phase₁ · phase₂ · phase₃ [· +N more] · K traces` — replaces
  the hard-coded `"3 traces"`. Members without a confirmed index contribute
  nothing to the dedup list.
- **Line 3:** `by <username> · edited <relative>`. If
  `created_by === currentUserId`, render `by you`. Relative-time formatting
  reuses any existing helper or `Intl.RelativeTimeFormat`.
- **Right-edge action zone:** `⚠` (when `has_stale_members`), `★/☆` pin
  (always rendered; unpinned star is dim and hover-only), `⋯` overflow
  (Delete, Copy link, Fork). Per §4 vocabulary, only the pinned star and
  the `⚠` chip stay always-visible.
- **Active marker:** accent left-border (1px → 3px), reuses
  `data-active="true"`.
- **Draft indicator:** `•` dot prefix and `(draft)` suffix on the title
  when the comparison has an open draft in this tab's sessionStorage.
  Age line reads `by you · just now`.

### 8.3 Empty states

- **No comparisons in this experiment:** `"No comparisons in this
  experiment yet."` + a prominent `+ New comparison` button.
- **No search match:** `"No matches for '<query>'. Clear search."`

### 8.4 Sort and grouping

- Pinned-first (preserves pin order from `useComparisonPins`), then
  `last_event_at` desc (replaces `updated_at` desc — `last_event_at` covers
  edits and chat activity).
- No per-author or per-tag grouping in v1.

## 9. Migration

**Frontend route flattening**
- `<Route path="*/edit" />` handler redirects via `<Navigate>` driven by
  `useLocation()` (NOT `window.location`, which lags under `MemoryRouter`
  in JSDOM tests).
  Tested in `frontend/test/compareRouteFlatten.test.tsx`.
- `comparePath()` in `lib/comparison/routes.ts` loses the `edit?: boolean`
  option (note: the option is named `edit`, not `isEdit`, in the current
  signature). Callers to update:
  - `pages/ComparePage.tsx::EditOrForkButton::onEdit` (becomes obsolete in
    Phase C — the page itself becomes editable in place).
  - `pages/ComparePageEdit.tsx` callsites that build cancel-to-review URLs.
  - `components/ConflictModal.tsx` "edit" navigations.
  - `components/WarmAddMenu.tsx` warm-add post-navigations.
  - `components/ComparisonSidebar.tsx` row-click and `+ New` navigations
    (already non-edit; verify no leftover `edit:true` flags).
- Drafts in `sessionStorage` (`himalaya-ui:compare-draft`,
  `COMPARE_DRAFT_VERSION = 1`) survive untouched — keyed by id, not URL.

**Backend schema**
- One additive migration `migrate_compare_view_choices!` in `db.jl` adds
  the three nullable view-choice columns on `comparisons`. All existing
  rows are NULL; frontend falls back to per-tab Zustand defaults.
- No `comparison_members` change. Phase summary uses `json_extract` on
  the existing snapshot column.
- No data backfill.

**No new event kinds.** The view-choice fields ride in the existing
`comparison_submitted` payload.

## 10. Testing strategy — six layers per `docs/contract-testing.md`

| Layer | Test file | What it pins |
|---|---|---|
| 1. Route emit | `test/test_routes_comparisons.jl` | Listing payload carries `author_username`, `member_count`, `member_phases`, `has_stale_members`, `last_event_at`. **Both** `POST /api/comparisons` (create) AND `POST /api/comparisons/:id/submit` accept + echo `view_grouping_mode` / `view_show_peak_ticks` / `view_show_peak_labels` — the create path is NOT optional. |
| 2. SSE payload | `test/test_idempotency_replay_invariant.jl` | Both `comparison_submitted` AND `comparison_created` SSE frames include view-choice fields when set. |
| 3. `applyRemoteToCache` | `frontend/test/queue/applyRemoteToCache.compare.test.ts` | Foreign `comparison_submitted` updates listing rows' new fields. |
| 4. Cache shape | `frontend/test/queue/cache-shape.test.ts` | `queryKeys.comparison(id)` and `queryKeys.comparisonsListing(scope)` shapes include the new fields. |
| 5. onMutate | `frontend/test/lib/queue/mutators/saveComparison.test.ts` | View-choice fields ride into `SaveComparisonInput` → `SaveComparisonBody` for both create and submit. |
| 6. Response shape | `test/test_route_response_shapes.jl` | Canonical JSON snapshot for the listing endpoint AND for `GET /api/comparisons/:id`. |

**Plus** a `forks_of_comparison` projection test — that helper goes
through `_comparison_listing_rows` too, so the new fields surface on
fork listings as a silent side-effect. Pin the contract.

**UX-specific tests**

- `frontend/test/Compare.singleMode.test.tsx` — pill appears on first edit,
  disappears on successful save, retains on 409.
- `frontend/test/Compare.dragThreshold.test.tsx` — 3px move = click; 5px move
  = drag. Edge cases around inline-edit on label.
- `frontend/test/ComparisonSidebar.projection.test.tsx` — phase summary
  truncation (`+N more`), author byline, draft indicator, hover-revealed
  pin / overflow.
- `frontend/test/Compare.rightRail.test.tsx` — Chat ↔ + Add traces tab
  switch; auto-switch on first dirty change; warm-add lands on + Add traces.
- `frontend/test/Compare.titleInline.test.tsx` — click title → input;
  Esc cancels; Enter commits.
- `frontend/test/Compare.viewChoicePersistence.test.tsx` — toggle grouping
  mode → dirty → save → reload → new viewer sees the same mode.
- `frontend/e2e/compare-single-mode.spec.ts` (mocked) — end-to-end:
  open comparison → edit title → add trace via drag → reorder → Save →
  URL never changes.
- `frontend/e2e/live/compare-fork-pill.spec.ts` (live) — non-author edits →
  pill says "Save as fork…" → fork created with correct lineage.

## 11. Build sequence

Each phase is independently shippable with green tests before the next
begins (per CLAUDE.md TDD convention).

| Phase | Scope | Reversible? |
|---|---|---|
| **A. Backend projection + migration** | New listing fields; view-choice columns; contract tests. | Yes (additive only) |
| **B. Route flattening** | Drop `/edit`, 301 redirect, merge `ComparePage` + `ComparePageEdit` into a single `Compare.tsx`. Behavior unchanged otherwise. | Reverts cleanly |
| **C. Title strip + toolbar + status surface** | §6. Save pill replaces Save/Cancel/Discard triplet. `EditOrForkButton`, `LineageBadge`, `ForksPopover` folded into the new header/toolbar. | Yes |
| **D. Right-rail tabs + warm-add unification** | §5. ChatCard + ComparisonPickerBody under tabs. `WarmAddMenu` simplified. | Yes |
| **E. Member row redesign** | §7. Collapse / expand, right-edge action zone, grab-anywhere drag with threshold, visible resize gap, drop-anywhere. `MemberMetaRow.tsx` / `MemberMetaGutter.tsx` rewritten. | Yes |
| **F. Sidebar redesign** | §8. New row layout consuming Phase A's projection; draft indicator; new empty states. | Yes |
| **G. Visual-language pass** | §4 tokens consistently across both gutters and the sidebar; `data-interactable` attributes added. | Yes |
| **H. Cleanup + deletes** | Remove dead code: `ComparePageEdit.tsx`, the `EditOrForkButton` inline function, possibly `LineageBadge.tsx` and `WarmAddMenu.tsx`. | Final step |

## 12. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Drag-vs-click threshold misfires (single click registers as drag) | Dedicated unit test (`Compare.dragThreshold.test.tsx`); 4px threshold from common practice; live spec covers the gesture in a real browser. |
| Single-page state machine handles view / draft / fork-draft / new draft | `useCompareMode()` returns a tagged union; single source for header copy, pill copy, and permissions. |
| Foreign SSE updates land while user is mid-edit | Existing `baseHash` freeze + 409 conflict modal already cover this. New addition: inline "Server-side updated since you last viewed" notice in the status surface (actor attribution deferred to Phase 2); modal still fires on Save. |
| Phase-summary projection latency for hot listings | Promote to denormalized `comparisons.member_phase_summary TEXT` column refreshed at save. No caller API change. |
| Visual-language pass leaks Tailwind classes into tests | All new interactive surfaces carry `data-interactable={kind}` for E2E selectors. |

## 13. Files touched (preview, non-exhaustive)

**Create**
- `packages/HimalayaUI/frontend/src/pages/Compare.tsx` — merged single-mode page.
- `packages/HimalayaUI/frontend/src/hooks/useCompareMode.ts` — tagged-union mode helper.
- `packages/HimalayaUI/frontend/src/hooks/useCompareDraftDirty.ts` — dirty-state signal.
- New tests per §10.

**Modify**
- `packages/HimalayaUI/src/comparisons.jl` — listing projection.
- `packages/HimalayaUI/src/routes_comparisons.jl` — accept view-choice fields on submit.
- `packages/HimalayaUI/src/db.jl` — `migrate_compare_view_choices!`.
- `packages/HimalayaUI/frontend/src/api.ts` — extend `ComparisonSummary` + `Comparison`.
- `packages/HimalayaUI/frontend/src/state.ts` — `compareRightTab` + draft helpers.
- `packages/HimalayaUI/frontend/src/components/ComparisonSidebar.tsx` — new row layout.
- `packages/HimalayaUI/frontend/src/components/MemberMetaRow.tsx` — collapse / expand, right-edge action zone, grab-anywhere drag.
- `packages/HimalayaUI/frontend/src/components/MemberMetaGutter.tsx` — visible resize gap, drop indicator, plot mirror.
- `packages/HimalayaUI/frontend/src/components/ComparisonPickerPanel.tsx` — drop-target advertisement; integration as a right-rail tab.
- `packages/HimalayaUI/frontend/src/components/WarmAddMenu.tsx` — simplified dropdown.
- `packages/HimalayaUI/frontend/src/lib/comparison/routes.ts` — drop `isEdit`.
- `packages/HimalayaUI/frontend/src/lib/queue/mutators/saveComparison.ts` — extend payload with view-choice fields.

**Delete or fold**
- `packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx` — folded into `Compare.tsx`.
- `packages/HimalayaUI/frontend/src/components/LineageBadge.tsx` — folded into title meta row.
- `packages/HimalayaUI/frontend/src/components/ForksPopover.tsx` — content folded into the `⋯ more` toolbar dropdown; file may stay as a list renderer.

## 14. Open questions for implementation time

- **Undo depth.** Cmd+Z undoes the last draft change up to depth 20 — is
  there a reason to make this deeper? Defer until usage.
- **Color picker swatch set.** The expanded member row's color swatch row
  needs a defined palette. Suggested: reuse `COMPARE_PALETTE` from
  `lib/comparison/coloring.ts` and add `↩ reset` to clear `color_override`.
- **Inline-edit on Esc with unsaved draft.** Today, Esc anywhere on
  `ComparePageEdit` is captured by `useFocusTrap` inside modals and
  otherwise unmanaged. Confirm Esc-on-page does NOT discard draft (only
  blurs the active editor) at implementation time.
- **Touch device behavior.** Pointer threshold is the same on touch
  (`pointerdown` is touch-agnostic). Long-press as alternate drag-start
  may help on tablet but is not in v1 scope.
- **Picker scope on `/compare/all`.** Global comparisons can mix exposures
  across experiments; the existing picker `ComparisonPickerBody` is
  experiment-scoped via `usePickerSamples(experimentId)`. Resolving "what
  does + Add traces show in global scope" needs one of: (a) a multi-select
  experiment switcher at the top of the picker; (b) the picker scopes to
  the experiment of the first-added member; (c) defer cross-experiment
  comparisons from + Add traces in v1 and require warm-add from Inspect
  for any cross-experiment addition. Recommend (b) for v1 — simplest, and
  almost-always-correct for the common single-experiment case.
- **`has_stale_members` for orphan members.** A member with
  `exposure_id IS NULL` (post-`ON DELETE SET NULL`) has no live exposure
  to hash-compare against. Following the existing per-member `is_member_stale`
  function in `comparisons.jl`, treat orphans as non-stale for this flag
  (they're flagged separately at render time as a deleted-exposure
  placeholder). Confirm at implementation time.
