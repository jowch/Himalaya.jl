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

There are exactly two things the user does on the Index page:

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

We fade to neutral gray (`var(--color-fg-dim)`) instead. This removes the
color signal entirely, so the hovered phase becomes the only chromatic
element on the canvas. The faded ticks are still readable; they just
stop *speaking color*.

The exception is curation state (e.g. excluded auto peaks): those keep
their identity color at low opacity, because the visual marker carries
meaning ("you excluded this"), not just context.

### 1.4 Identity at boundaries: meaning, not PK

State that lives across reanalysis (the user's "active set" of indices)
must be tracked by **what it means** — phase + basis — not by database
auto-increment IDs. When `persist_analysis!` deletes and re-inserts every
candidate row, the PKs are gone, but the *meaning* of the user's pick is
still there: "I chose the Pn3m candidate near basis 0.5."

We translate at the boundary. The custom-group's `(phase, basis)` tuples
are snapshotted before delete, then re-attached to the closest fresh
candidate of the same phase within a tolerance. If no near match
survives the peak edit, the pick is **honestly dropped** — that change
in peaks invalidated this indexing, and it's better to surface the
absence than to silently carry stale state.

### 1.5 Server state vs client state

TanStack Query owns server state (experiments, samples, exposures, peaks,
indices, groups, messages). Zustand owns client state (active sample,
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
landed (Miller inset top-right → bottom-left → outside the plot;
circles above triangles → vlines from top; faded color → faded gray).
Each reversal cost an hour and saved a permanent worse state.

---

## Part 2 — How we chose to do it (this iteration)

These are the concrete choices in the current build. They are not
*principles*; they are working answers. Expect changes here.

### 2.1 Layout

- **Three cards on the Index page**, ratio 22 / 56 / 22 (chat / plot / indices)
  on a max-width 1600px workspace. Below 1100px the chat reflows below the
  other two.
- **Vertical structure**: app utility row (44px) → tab rocker row → page
  body. The tab rocker (Index / Compare) sits where a per-page title used
  to live; the actual page title moved into the plot card's top strip.
- **Card-header utility** (`.card-header`, height 56px, 1rem padding,
  `flex items-center`) shared between the plot card's title strip and
  the indices card's "Index choices" header so their top edges line up.
- **Compare page** is a placeholder; the tab exists, the content does not.

### 2.2 Typography

- **Plus Jakarta Sans** everywhere. We tried mono for kbd hints,
  timestamps, and stat labels — it added a second visual rhythm without
  earning its complexity. Sans-only reads more cohesive.
- The `--font-mono` CSS variable is still defined in the theme so we can
  re-introduce it for genuinely-monospace surfaces (raw data dumps, SMILES
  strings) if those ever land. Currently unreferenced.

### 2.3 Color

- OKLCH throughout. Perceptually uniform, so chroma and lightness move
  together predictably across hues.
- **Peaks are bright/saturated**. Auto = `--color-accent` (ice blue, ~220°).
  Manual = `--color-peak-manual` (neon magenta, ~340°), well clear of every
  index hue.
- **Indices are muted/earthy**, chroma 0.07–0.10, lightness 0.70–0.78. The
  palette in `phases.ts` covers Pn3m (terracotta), Im3m (sage), Ia3d
  (periwinkle), Fm3m (mauve), Fd3m (dusty rose), Hexagonal (ochre),
  Lamellar (dusty teal), Square (amber-tan). No phase is within ±25° of
  the peak hues.
- **Faded annotations render in `--color-fg-dim`** (neutral gray) at
  reduced opacity, *not* at the phase color with reduced opacity.
- `:root { color-scheme: dark }` (and `light` on the matching theme
  override) re-themes native form controls and scrollbars without
  per-pseudo-element styling.

### 2.4 Trace plot

- Observable Plot at the core. We wrap it in a `host > plotContainer +
  overlayRef` nested DOM so React can manage the overlay SVG (cursor,
  ticks, peak triangles) without `replaceChildren` wiping it on each
  render.
- **X axis label**: `q (Å⁻¹)` with a custom `tickFormat` that uses plain
  decimal or scientific notation. Plot's default SI formatter rendered
  e.g. 0.040 as "40 m" (milli-) which made no sense in this context.
- **Wheel scroll** zooms around the cursor; **double-click** resets to
  full range. The visible q-range is shared with numeric inputs in the
  plot card's title strip (the `QRange` controls).
- **Click empty plot space** adds a manual peak at the exact clicked q.
  **Click within ~10 pixels of an existing peak triangle** removes it.
  No q-snap — the user zooms in for precision.
- **Cursor crosshair** (dashed vertical rule + follow-dot tracking the
  trace) only appears inside the plot interior, gated by the plot's
  margin constants.
- **Peak markers** (down-triangles) sit ~30% above the trace line, not
  on it, so they don't get visually swallowed by the curve.
- **Predicted-q ticks** for the active group render as short top-of-frame
  marks. Hovering an alternative index dims the active ticks to neutral
  gray and adds the hovered index's ticks at full chromatic strength.

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

`pipeline.jl::persist_analysis!` does the snapshot-delete-recreate-reattach
dance described in §1.4:

1. Snapshot custom-group members as `(group_id, phase, basis)` tuples.
2. Delete prior auto peaks, indices, auto-group rows, and stale custom-
   group memberships.
3. Re-insert candidates with fresh PKs.
4. For each snapshotted member, find the closest fresh candidate of the
   same phase. If the basis delta is within `MEMBER_REATTACH_RELTOL = 0.05`
   (5% of the snapshotted basis), re-attach. Otherwise drop.
5. If any custom group survives with at least one re-attached member,
   the freshly-created auto group is demoted to `active = 0`. The
   curated set stays the active set.

`queries.ts::invalidateExposure` invalidates `peaks`, `indices`, **and
`groups`** so the right-rail Active set updates immediately after any peak
edit triggers auto-reanalysis.

### 2.7 Persistence

Zustand `persist` middleware → `localStorage`. Persists `username`,
`activeExperimentId`, `activeSampleId`, `activePage`, `theme`,
`tutorialSeen`. Refresh lands the user back on the same scope; switching
machines starts them over (intentional — server doesn't know about
"current view").

### 2.8 Onboarding

`OnboardingFlow.tsx` triggers when `username` is missing:

1. Pick existing user, or `+ New user`.
2. **Tutorial slides** (4 short steps) only on the new-user path. Returning
   users skip. `tutorialSeen` persists.

### 2.9 Backend shape (unchanged from previous iterations)

- SQLite-per-experiment, schema in `db.jl`. Oxygen.jl REST routes,
  one file per resource (`routes_experiments.jl`, `routes_peaks.jl`, …).
- New this iteration: `sample_messages` table (FK `author_id → users.id`,
  `ON DELETE SET NULL`) + `routes_messages.jl`. Backs the per-sample
  ChatCard.
- Exposure, tag, and notes endpoints are intact even though the UI
  doesn't surface them right now (see §1.9).

---

## Part 3 — Open questions for next iterations

Things we know we don't yet have a good answer for:

- **Compare page content.** The tab is a placeholder. What's the actual
  primary task on Compare? Stacked / waterfall plots? Side-by-side
  per-sample cards? We don't know yet, and putting in something generic
  would lock in the wrong shape.
- **Exposure triage.** A sample often has 5–20 exposures and the UI
  currently auto-picks the first by id. A Lightroom-style filmstrip
  ("good / bad / maybe", keyboard-driven) is a natural fit, but we
  haven't designed it.
- **Reviewer workflow.** Multi-user science: is this analysis "approved"?
  By whom? The `users` and `user_actions` tables exist; the UI for
  promotion / review does not.
- **Per-peak hkl labels on the plot.** Useful, but visual budget is
  tight. We need to know which use-cases actually demand them before
  spending the ink.
- **Color-blind accessibility for phase hues.** Eight hues at similar
  chroma are not safe for all forms of color vision. A dash-pattern
  channel (ticks) and shape channel (Miller dots) could make the phase
  distinction redundant with hue. Not done.
- **Chat threads / mentions / reactions.** ChatCard is intentionally a
  flat list right now. We added the message thread to test whether
  per-sample conversation is even a thing scientists do. Wait for
  evidence before layering features on.

See [`docs/future-feature-ideas.md`](future-feature-ideas.md) for the
running list of deferred work, including these and analysis-engine ideas
that aren't UI-shaped (extended lattice types, sub-pixel peak positions,
background subtraction in the pipeline).
