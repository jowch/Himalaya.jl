# HimalayaUI — Design Philosophy

This document captures **how we think about HimalayaUI**, separately from the
specific shape it has today. The first half is the durable stuff — principles
we keep returning to. The second half describes the choices made for the
current iteration; expect those to drift as we learn more.

---

## Part 1 — How we think about it

### 1.1 Scientific instrument, not dashboard

HimalayaUI is a tool a scientist uses for *minutes-to-hours* at a stretch
while looking carefully at one trace at a time. That's a different posture
from the *seconds-at-a-glance* read of a status dashboard.

Implication: the trace itself is the loudest thing on the page; everything
else recedes until it is asked for. Persistent chrome that competes with
the plot — repeated metric strips, always-on sidebars, dense filters — is
the failure mode we are most allergic to.

When in doubt: **does this widget help the user *see* the data, or is it
asking them to track something else?**

### 1.2 Two workflow steps, visually separate

There are exactly two things the user does on the Focus page (`/sample/:id`):

1. **Pick peaks.** Foreground, active edit. The user is making a claim
   about what's a real diffraction peak vs noise.
2. **Assign indices.** Interpretation. Given the peaks, what phase
   indexing best explains them?

These are different cognitive moves and they should look different.
Peaks read as **bright/neon** (active edit, demanding attention). Indices
read as **muted/earthy** (interpretation context, present but quiet).
Hue alone tells the eye which thing it is looking at.

### 1.3 Fade to neutral, not to a paler version of the color

When attention shifts to one element (e.g. hovering a candidate index),
everything else must *get out of the way*. Multiplying alpha leaves "faded
orange" still distinctly orange — competing for the eye's color signal at
the same hue, just dimmer.

We fade to neutral gray (`var(--color-ink-soft)`) instead. This removes the
color signal entirely, so the hovered phase becomes the only chromatic
element on the canvas. The faded ticks are still readable; they just
stop *speaking color*.

The exception is curation state (e.g. excluded auto peaks): those keep
their identity color at low opacity, because the visual marker carries
meaning ("you excluded this"), not just context.

### 1.4 Identity at boundaries: meaning, not PK

State that lives across reanalysis (the user's "active set" of indices,
their excluded auto peaks) must be tracked by **what it means** — phase +
basis, q-value of the excluded peak — not by database auto-increment IDs.

In the common case the implementation now preserves PKs across reanalysis
(see §2.6), so the principle is mostly invisible. But it asserts itself
the moment peaks genuinely shift: the user's curation is then re-resolved
against the new candidates by *meaning*. If no near match survives the
peak edit, the pick is **honestly dropped** — that change in peaks
invalidated this indexing, and it's better to surface the absence than
to silently carry stale state.

### 1.5 Server state vs client state

TanStack Query owns server state (experiments, samples, exposures, peaks,
indices, groups). Zustand owns client state (active sample,
hovered index, theme, persisted scope, modal step). They never overlap.

This split is load-bearing. Mutations invalidate scoped query keys
(`peaks(id)`, `indices(id)`, `groups(id)`); cache freshness comes from
the server, not from client diff-tracking.

When the server changes shape, the client's optimistic mirrors don't
go stale because there are none. When the user navigates, the client
state updates instantly without a network round-trip.

### 1.6 Keyboard as a first-class surface

A scientist in flow uses both hands. The interface yields to that with:

- `/`, `⌘K` — open the navigation modal
- `,` / `.` — step backward / forward through samples
- `T` — toggle theme
- `Esc` — close the active modal

These are not "power user" extras layered on at the end. The modal-based
nav and the title-as-button design exist *so that* the keyboard path is
parallel to the click path. If a feature only works by mouse, that's a
warning that we missed something.

### 1.7 Modal as a deliberate moment

Switching samples is a context shift, not a continuous interaction.
A persistent sidebar of samples invites users to scan instead of decide.
A modal forces a moment: "I am switching context now," then dismisses
itself, leaving the workspace clean.

This is also why the modal cascades (experiment → sample) instead of
flattening: the cascade matches the user's mental model of *which
question they're answering at each step*, and Backspace/chip-removal
makes "going back" cheap.

### 1.8 Title as primary action

The page title is not chrome. It is the affordance that tells you
**what scope you are looking at** *and* lets you change it. Clicking
the title opens the nav modal at the right step. There is no separate
"change sample" button.

This collapses two roles into one element and removes a class of
"how do I switch samples?" confusion.

### 1.9 Defer rather than design vaguely

When a feature is real-but-later (exposure triage, tag editing, peer
review workflow, hkl labels), we **delete the UI and keep the schema
stable**. Half-built features rot, persistent placeholder widgets
fill the page with apologies, and frozen design committees produce
bad designs.

The deletes leave clear marks in git history; the schema columns sit
dormant until a future iteration revisits them. Adding the UI back is
a tracked task, not a re-litigation.

### 1.10 Texture creates depth, not noise

The grained obsidian background is doing work: it gives the floating
cards something to *float on*. A solid black would make the cards feel
glued to the screen; a gradient would compete with the plot.

Subtle texture (fractal-noise SVG at low opacity, blended) reads as
*surface* without reading as *content*. The cards then have somewhere
to be, instead of being everything.

### 1.11 Honest dropping over silent staleness

When the user's curation no longer matches reality (e.g. they edited
peaks and a previously-active indexing now doesn't fit), we drop the
stale state. We do not carry it forward as if it were still meaningful,
and we do not show a broken intermediate state with members pointing
at deleted rows.

The user prefers seeing "the thing you picked is gone" to seeing
"something invisible is wrong."

### 1.12 Iteration over completeness

Ship the simpler thing. Live with it. Notice what actually breaks. Fix
that. The current iteration is *not* the spec — it's a checkpoint.
Several decisions in §2 below started out the opposite of where they
landed (Miller inset top-right → bottom-left → outside the plot, since
removed entirely; circles above triangles → vlines from top; faded color
→ faded gray). Each reversal cost an hour and saved a permanent worse state.

---

## Part 2 — How we chose to do it (this iteration)

These are the concrete choices in the current build. They are not
*principles*; they are working answers. Expect changes here.

### 2.1 Layout

- **Focus workspace** (`/sample/:id`): a trace hero (the loud plot) plus a
  detector panel, a phase-call rail, and a notes margin/drawer, on a
  max-width workspace. The plot is the loudest thing on the page; the rail
  and panels recede until asked for.
- **Vertical structure**: app utility row → page body. The shell is
  URL-routed (`CorpusShell` / `CorpusTopbar`); the page title moved into the
  workspace's top strip. The earlier Index/Compare tab rocker was removed in
  the greenfield cutover (merge `dcac451`, PR #281).
- **Card-header utility** (`.card-header`, height 56px, 1rem padding,
  `flex items-center`) shared between the plot card's title strip and
  the indices card's "Index choices" header so their top edges line up.
- **Title strip breadcrumb:** the experiment name renders in `text-ink-soft`
  even when set; the sample name uses `text-ink`. This is intentional —
  the experiment is context (the container) and the sample is the leaf the
  user is actively working on. De-emphasising the parent draws the eye to
  the sample without hiding the experiment path.
- **Series stage** (`/series`): the former Compare page was retired and
  folded into Series, which renders multi-trace overlays through a
  folio / scoping / builder set of pages. Deep `/compare*` links redirect to
  `/series`.

### 2.2 Typography

- **Plus Jakarta Sans** everywhere. We tried mono for kbd hints,
  timestamps, and stat labels — it added a second visual rhythm without
  earning its complexity. Sans-only reads more cohesive.
- The `--font-mono` CSS variable is defined in the theme and used by the
  `.text-data` / `.text-data-strong` roles (`styles.css`) for genuinely
  monospace data surfaces (e.g. q-values, flag labels).

### 2.3 Color

- OKLCH throughout. Perceptually uniform, so chroma and lightness move
  together predictably across hues.
- **Peaks are bright/saturated**. Auto = `--color-accent` (ice blue, ~220°).
  Manual = `--color-peak-manual` (neon magenta, ~340°), well clear of every
  index hue.
- **Indices are vivid but calm**, chroma 0.10–0.13, lightness 0.76–0.80,
  inspired by the v3 design reference's sage/amber/violet anchors. The
  palette in `phases.ts` covers Pn3m (amber), Im3m (sage), Ia3d (violet),
  Fm3m (coral), Fd3m (rose-purple), Hexagonal (seafoam teal), Lamellar
  (periwinkle), Square (chartreuse). Hues stay clear of ±20° of the peak
  hues (220° accent, 340° manual) and the high-chroma warning zone (~75°).
- **Faded annotations render in `--color-ink-soft`** (neutral gray) at
  reduced opacity, *not* at the phase color with reduced opacity.
- `:root { color-scheme: dark }` (and `light` on the matching theme
  override) re-themes native form controls and scrollbars without
  per-pseudo-element styling.

### 2.4 Trace plot

- A bespoke d3 trace-plot engine at the core (`print/plot/TracePlot.tsx`,
  built on `d3-scale`/`d3-shape`; axis formatting is hand-rolled, not
  `d3-format`). React owns the SVG
  tree directly — trace path, overlay marks (cursor, ticks, peak triangles)
  and axes are all declarative React elements, so there is no imperative
  `replaceChildren` to fight on each render. Observable Plot was fully
  retired (2026-06-13); `@observablehq/plot` is no longer a dependency.
- **X axis label**: `q (Å⁻¹)` (`print/plot/Axis.tsx`) with a hand-rolled tick
  formatter (`lib/plot/formatAxis.ts`) that uses plain decimal or scientific
  notation — never SI-prefix abbreviations, which would render e.g. 0.040 as
  the nonsensical "40 m" (milli-).
- **Wheel scroll** zooms around the cursor; **double-click** resets to
  full range. The visible q-range is shared with the numeric inputs in the
  plot card's title strip (the controlled `xDomain` / `setXDomain` state on
  `TracePlate`; there is no discrete `QRange` component).
- **Click empty plot space** adds a manual peak at the exact clicked q.
  **Click within ~10 pixels of an existing peak triangle** removes it.
  No q-snap — the user zooms in for precision.
- **Cursor crosshair** (solid vertical rule + follow-dot in `--color-ink-soft`)
  only appears inside the plot interior, gated by the plot's margin
  constants. Neutral gray so it doesn't compete with the phase-coloured
  ticks.
- **Peak markers** (down-triangles) sit ~30% above the trace line, not
  on it, so they don't get visually swallowed by the curve. Predicted-q
  ticks that match a peak terminate ~7px above the triangle so the two
  markers don't visually fuse.
- **Predicted-q ticks** render in two places:
  - A **track row** above the plot (`TRACK_Y_TOP`–`TRACK_Y_BOTTOM` band
    inside the bumped `MARGIN_TOP`) carries the persistent phase-colour
    swatches at 55% opacity. This is where colour lives by default.
  - **Plot vlines** inside the data area are neutral gray
    (`--color-ink-soft`) at 35% opacity by default — they show *where* an
    index would land without competing with the trace data.
  - On hover, the hovered index's ticks go solid full-opacity in *both*
    places (track and plot), and any other indices fade to gray 30% in
    both places. The track row doubles as a legend-by-position.

**Why bespoke d3 rather than Observable Plot.** The decisive reason was
*inverted projection ownership*: Plot owned the q→pixel mapping, so every
interactive layer had to reverse-engineer the projection back out of Plot's
rendered DOM (`plotEl.scale('x')`), forcing a ~420-line imperative bridge
(`Plot.plot()` → `replaceChildren`). The owned engine exports the projection
(`print/plot/projection.ts` — `makeProjection`/`makeAxis`, backed by
`d3-scale`), so every layer reads the same scale declaratively.

**The comb / residual panel is a separate SVG subsystem, not part of
`<TracePlot>`.** It lives in `print/comb/` (`CombScaffold`, `CombChart`,
`ResidualChart`; `CombsPanel` in `print/components/`) and rides the shared
q-link spine (predicted-q, `phaseColor`) but is never feature-gated into
`TracePlot` — an explicit boundary. Within it:

- **`CombChart` and `ResidualChart` are siblings, not one component with a
  `view` prop.** They share only the log-q x-projection and the rowed
  scaffold; their glyph geometry (teeth/carets/rings vs zero-line/band/scatter)
  is entirely different, so a unified `view`-prop component would violate
  one-responsibility.
- **Both are fixed to log-q**, with no per-panel log/linear toggle and *no*
  coupling to the hero trace's log/linear state. The comb is a ratio diagram
  (teeth at √N multiples of q₁): even in log, crowded in linear.
- **The hovered-q "hot" highlight uses a wider stem + larger ringed cap in
  the row's own phase colour — never the terracotta accent.** This diverges
  from the old mockup CSS (`.tooth.hot { stroke: var(--accent) }`): the
  accent is reserved for the active assignment, not hover.
- **`ResidualChart` uses a fixed symmetric y-domain (±`RESID_DOMAIN` = ±3%)**,
  not auto-scaled per row, with a tolerance band at ±`RESID_BAND` = ±2.2%.
  Display tiers: within band → filled dot; outside band but within domain →
  hollow dot; beyond domain → clamped to the edge.
- **`CombLegend` legends the trace-plot peak-glyph vocabulary** (auto
  triangle / manual diamond / predicted-absent caret / excluded), *not* the
  comb's intrinsic tooth/caret/ring vocabulary — even though `CombsPanel`
  places it below `CombChart`. The two vocabularies deliberately differ.

**The series waterfall composes N single-trace `TracePlot`s** (one per
member), not one shared coordinate space. This per-member partitioning is the
deliberate replacement for the legacy `computeYBands`/`computeReference`/
`applyNormalization` machinery in `lib/comparison/` (those files still exist
but `src/print/` imports none of them). The vertical **`PhaseStrip`
companion** (`orientation="vertical"`) and the waterfall stack share the same
low→high `display_order`; pixel alignment between them is the *plate's*
responsibility — `WaterfallChart` does not embed the strip.

### 2.5 Animation

Three primitives, all 120–140ms ease:

- `.anim-pal-in` — backdrop / overlay fade (opacity 0 → 1).
- `.anim-pal-scale` — popover entry: `translateY(-6px) scale(0.98)` → none.
- `.anim-overlay … overlay-fade-in` — one-shot fade-in for newly mounted
  SVG overlay nodes (peak triangles, ticks). We use a CSS animation rather
  than transitions because the overlay rebuilds on every render — there is
  no element-to-element transition to fire.

The CSS transitions on `*` (color/background/border/opacity, 120ms ease-out)
soften ambient theme/hover changes without us having to opt-in per element.

### 2.6 Backend: active-set preservation

We used to do a snapshot-delete-recreate-reattach dance every time a peak
edit triggered reanalysis. Plan 7 replaced it with two layered mechanisms
that keep the §1.4 principle intact while paying for themselves in
multiplayer-readiness and cold-start latency:

1. **Curation lives in its own table.** `peak_curations` (kind `add` or
   `exclude`) is separate from the machine output in `auto_peaks`.
   `effective_peaks(db, exposure_id, q, I)` synthesises the working set
   as `auto − excludes ∪ adds`. Reanalysis no longer destroys curation
   rows, so there is nothing to snapshot for excluded peaks.
2. **Auto peaks are diff-updated, not recreated.** `diff_update_auto_peaks!`
   matches old and new findpeaks output by q within tolerance and issues
   `UPDATE` for survivors and targeted `INSERT`/`DELETE` for the rest.
   Auto peak IDs are preserved across reruns, which means
   `index_peaks(peak_id)` references stay valid for the common case and
   never need re-resolution.
3. **Indices fall back to phase+basis matching only when peaks shift
   beyond tolerance.** When `diff_update_auto_peaks!` reports that some
   referenced peaks vanished or moved past `MEMBER_REATTACH_RELTOL = 0.05`
   (5% of basis), `_persist_analysis_inner!` re-resolves the affected
   speculative indices against the new candidate set by phase + closest
   basis. No match within tolerance → the pick is dropped (honest
   staleness, §1.11).
4. **Hash memoization skips work entirely when inputs haven't changed.**
   `analyze_exposure!` SHA-256s the trace bytes (`trace_hash`) and the
   effective peak set (`analysis_inputs_hash`); when both match the
   stored hash and the corresponding rows already exist, findpeaks and
   indexpeaks skip. The active set survives untouched because no writes
   fire. See [`docs/event-log.md`](event-log.md) for the full hash
   contract.

Cache reconciliation flows through the mutation queue (Plan 8): own-op
confirmations and foreign SSE events are funnelled through
`handleRemoteEvent` / `applyRemoteToCache.ts` (`lib/queue/`), which updates
the `peaks` / `indices` / `groups` caches for the affected exposure so the
Active set updates immediately after a peak edit triggers auto-reanalysis.
SSE multiplayer fan-out (§2.10) takes the same path. See
[`docs/mutation-queue.md`](mutation-queue.md).

### 2.7 Persistence

Zustand `persist` middleware → `localStorage`. Persists `username`
(+ name), `activeExperimentId`, `activeSampleId`, `tutorialSeen` (the
dual-nav `activePage` and the `theme` toggle were both retired — routing
is URL-based and there is one fixed palette). Refresh lands the user back
on the same scope; switching
machines starts them over (intentional — server doesn't know about
"current view").

### 2.8 Onboarding

`OnboardingFlow.tsx` triggers when `username` is missing:

1. Pick existing user, or `+ New user`.
2. **Tutorial slides** (4 short steps) only on the new-user path. Returning
   users skip. `tutorialSeen` persists.

### 2.9 Backend shape

- SQLite-per-experiment, schema in `db.jl`. Oxygen.jl REST routes,
  one file per resource (`routes_experiments.jl`, `routes_peaks.jl`, …).
- **Peaks split (Plan 7).** The legacy `peaks` table is gone. Machine
  output lives in `auto_peaks`; user curation lives in `peak_curations`
  (kind `add` or `exclude`). `index_peaks` carries a `peak_kind`
  discriminator so a single index can reference both. The `/api/peaks`
  route synthesises the legacy union via `UNION ALL` for the frontend.
- **Structured event log (Plan 7).** `user_actions` is now the
  source-of-truth log: `payload` (JSON), `undoes_event_id`, indexed
  per-exposure. View tables (`peak_curations`, `index_group_members`)
  are written exclusively by the dispatcher in `events.jl`. See
  [`docs/event-log.md`](event-log.md) for the dispatcher contract,
  hash memoization, and SSE multiplayer.
- **Chat (`sample_messages`).** FK `author_id → users.id` with
  `ON DELETE SET NULL`. Persists as a parked data plane — the chat
  presentation (ChatCard / @-mention subsystem) was retired 2026-05-29 and
  there is no UI on top of it.
- Exposure, tag, and notes endpoints are intact even though the UI
  doesn't surface them right now (see §1.9).

### 2.10 Multiplayer (Plan 7 R5a)

- `GET /api/events` is a Server-Sent Events stream; every client opens
  one EventSource on mount.
- `apply_event!` fires `broadcast_event!` *after* its transaction
  commits, so subscribers never see a rolled-back event. Process death
  between commit and broadcast loses the frame but not the event
  (durable in `user_actions`); clients reconcile on EventSource
  reconnect via TanStack Query refetch.
- The frontend SSE handler lives in `print/App.tsx` (the `PrintApp` root):
  it opens one EventSource and routes every `curation` frame to
  `handleRemoteEvent` (`lib/queue/replayCoordinator.ts`), the mutation-queue
  replay-as-rerun path that reconciles the `peaks` / `indices` / `groups`
  caches for the affected exposure. See [`docs/mutation-queue.md`](mutation-queue.md).
  **Self-echo filter:** the guard keys on a **per-tab `client_id`**
  (`lib/clientId.ts`) — a frame whose `client_id` matches the local tab's id
  is treated as own-op, so the tab's own optimistic UI isn't clobbered by a
  replay of its own write. Two tabs of the same user are distinct subscribers
  (the guard is per-tab, not per-username).
- **Conflict resolution was cancelled (2026-06-03).** The `If-Match` +
  409-retry path (R5b) and the Series conflict-resolution modal were both
  dropped — there is no conflict UI. Last-write-wins is permanent; the
  conflict story is instead being addressed by edit-tracking → undo/redo →
  versioning (designed in Layer 4). See
  [`docs/event-log.md`](event-log.md) §"Conflict resolution".

### 2.11 Figure export

- **Content-WYSIWYG, style-transformed, via a tiny split-button.** The
  `ExportButton` (`print/components/`) is a split control: the main face
  copies a clipboard PNG; the chevron opens Download-as-PNG / Download-as-SVG.
  There is no preview, no options dialog, no sheet or modal — the deliberate
  non-design (the deferred `ExportSheet`/`ConflictModal` were collapsed away).
- **The presentation styling is the adapter's job, not the button's.** The
  call site hands `useFigureExport` a style-agnostic `renderSvg: () => string`
  thunk that returns the clean scientific-presentation SVG (white background,
  GraphPad-style figure); the hook and button never know the styling. So the
  on-screen plot can stay dark/interactive while the exported figure is print
  styling. (Hook re-entry is guarded by an `inFlight` ref — see the frontend
  `AGENTS.md`.)

---

## Part 3 — Open questions for next iterations

Things we know we don't yet have a good answer for:

- **Exposure triage.** A sample often has 5–20 exposures and the UI
  currently auto-picks the first by id. A Lightroom-style filmstrip
  ("good / bad / maybe", keyboard-driven) is a natural fit, but we
  haven't designed it.
- **Reviewer workflow.** Multi-user science: is this analysis "approved"?
  By whom? The substrate is now in place — `user_actions` is a
  structured event log with payloads and `undoes_event_id` (Plan 7 R4),
  so a per-exposure audit view or "promote to approved" workflow can
  be built without further schema work. The UI for promotion / review
  is the missing piece.
- **Per-peak hkl labels on the plot.** Useful, but visual budget is
  tight. We need to know which use-cases actually demand them before
  spending the ink.
- **Color-blind accessibility for phase hues.** Eight hues at similar
  chroma are not safe for all forms of color vision. A dash-pattern
  channel (ticks) and shape channel (Miller dots) could make the phase
  distinction redundant with hue. Not done.
- **Chat threads / mentions / reactions.** The chat UI was retired
  2026-05-29 (presentation deleted, including the @-mention subsystem); only
  the parked `sample_messages` data plane remains. The original experiment —
  per-sample conversation — never gathered enough evidence that scientists
  want it. The open question is whether to revive it at all before layering
  threads / mentions / reactions on top.

See [`docs/future-feature-ideas.md`](future-feature-ideas.md) for the
running list of deferred work, including these and analysis-engine ideas
that aren't UI-shaped (extended lattice types, sub-pixel peak positions,
background subtraction in the pipeline).
