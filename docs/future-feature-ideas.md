# Future Feature Ideas

Ideas that are intentionally out of scope for current development but worth
preserving for later planning. Each entry says *what* and *why later* — not
a full design.

For the design philosophy and the current iteration's choices, see
[`himalayaui-design.md`](himalayaui-design.md).

---

## Analysis engine

### Extended lattice types

Monoclinic and tetragonal lattice indexing. These require fitting a lattice
parameter *vector* (2+ free parameters) rather than a single basis, which
changes how `indexpeaks` and `score` work internally. Will require extending
the `Phase` type hierarchy and the indexing engine. Design the extension points
when the need is concrete.

### Sub-pixel peak positions

Parabolic interpolation around detected peak maxima for more precise q values.
Currently peaks are returned at grid positions.

### Background subtraction in pipeline

Automated SNIP or asymmetric least-squares background subtraction as a
pre-processing step, particularly useful for traces with steeply falling
backgrounds. Currently handled manually (background-subtracted exposures are
ingested as derived exposures in the DB).

### Auto-best-exposure selection

Heuristics for automatically selecting the best exposure for a given sample:
- Detect flares (anomalous low-q intensity spikes)
- Score exposures by signal quality (peak prominence, background level)
- Flag poor-signal exposures automatically

### `Fm3m` in `indexpeaks` dispatch

The all-phases loop in `src/index.jl` omits `Fm3m`. The phase is defined and
`minpeaks`/`phaseratios` exist, but `indexpeaks` can never return an `Fm3m`
index. Known pre-existing gap.

---

## Index page — annotation depth

### Peak classification on ticks

Color-code predicted-q ticks by whether they are explained / unexplained /
excluded. Currently all ticks for the active group render in the phase color;
adding a "missing peak" or "excluded peak" channel would let the user see
fit gaps at a glance.

### `hkl` labels above active-group ticks

Render the Miller indices `(hkl)` next to each predicted-q tick on the trace.
Useful for publication figures and for users learning the indexing. Visual
budget is tight, so probably hover-revealed rather than always-on.

### Manual index entry modal

Currently an indexing must come from `indexpeaks`. A scientist who knows
the answer should be able to enter `Pn3m, basis 0.523` directly and have
it scored against the current peaks without going through candidate
generation.

### Color-blind accessibility on phase hues

Eight earthy phases at similar chroma are not distinguishable under all
forms of color vision. Add a redundant channel — dash patterns on ticks,
shape on Miller dots — so phase identity doesn't depend on hue alone.

### Multi-exposure overlay on the trace

Overlay multiple exposures (or derived exposures) on the same plot with
configurable I-offsets. Lets the user see flares and cross-exposure
consistency at a glance.

---

## Index page — data triage

### Exposure-triage page

Lightroom-style filmstrip of all exposures in a sample. Keyboard-driven
good / bad / maybe rating. The current Index page auto-picks the first
exposure by id; this page would let the user curate which exposure the
analysis runs on.

### Tag editing UI

Backend (`sample_tags`) and routes (`routes_tags.jl`) are intact, but the
UI was dropped from the three-card redesign. Re-introduce when we know
what tag-driven workflows actually look like — probably alongside the
summary table or the cross-experiment comparison page.

### Stale-indices banner

Detect when peak edits have happened since the last reanalysis and surface
a banner offering to re-run. The auto-reanalysis chain in `queries.ts`
covers most cases now, but a manual fallback is useful for batch edits.

---

## Compare page (currently a placeholder)

### Stacked / waterfall comparison

Multi-sample trace overlay with configurable I-offset between traces.
Publication-quality SVG export. Useful for visualizing phase transitions
across a sample series.

### Summary table

Full-screen tabular view of all samples in an experiment: confirmed phases,
lattice parameters, R², score. Sortable, filterable by tag. Export to
CSV/JSON.

### Cross-experiment sample comparison

Query `sample_tags` and analysis results across multiple experiment
databases to compare samples with similar compositions (e.g., same
lipid/peptide system across beamtimes). The SQLite-per-experiment schema
is migration-friendly: a thin aggregation layer can open multiple DBs and
JOIN across them via SQLite's `ATTACH` mechanism.

---

## Multi-user / review

### Reviewer workflow

`reviewed` status transitions on indexings, with a reviewer attribution.
The `users` and `user_actions` tables exist; the UI for promotion / review
does not.

### Per-user audit view

`user_actions` rows already record who did what. A view that shows a
user's recent actions (peak edits, group changes, message posts) would
support both review and onboarding ("show me what alice has been doing").

### Chat: mentions, threads, reactions

Current ChatCard is intentionally a flat per-sample list. Wait for usage
data before adding `@user`, threading, or emoji reactions — these are
expensive to design well and easy to design badly.

### Avatars

Currently a 2-letter initials chip in the utility cluster. If multi-user
becomes more central, real avatars (uploaded or Gravatar-style) could
distinguish authors in chat and audit views.

---

## Onboarding & docs

### Tutorial content beyond the four intro slides

Current tutorial covers title-button, three-card layout, sample stepping,
and active-set editing. Add slides for keyboard shortcuts, the chat card,
and theme toggle as those features grow more load-bearing.

### Empty / error state mock-ups

Default Hint/Loading text is fine for development. Real empty states ("no
exposures yet — drop a `.dat` here") and error states ("analysis failed;
click to see traceback") need design.

---

## Operations / configuration

### Beamline-config editor

Each experiment has parameters (wavelength, sample-detector distance,
calibrant) that currently live outside the DB. Expose these in the
experiment metadata and add an editor UI.

### Manifest editor

Manifest CSV is currently edited externally and re-imported. An in-app
editor (with sample-name validation against the data dir) would close
the loop.

### Export UI

Backend has `routes_export.jl` (JSON / CSV). No UI surface yet — add a
download button in the Compare or summary-table page when those exist.

### Derived-exposure construction

The schema supports derived exposures (e.g. background-subtracted) but
there's no UI to create them. Probably belongs on the exposure-triage
page when that lands.

---

## Frontend infrastructure

### Code splitting

Vite warns the bundle is >500 kB. Probably split: NavModal +
OnboardingFlow + ComparePage as lazy chunks; pull `@observablehq/plot`
out of the main bundle if possible.

### `ResizeObserver` on TraceViewer

Currently the plot's container width is read at render time. Window
resizes don't re-run the effect unless data changes. Add a `ResizeObserver`
when this becomes a noticeable problem.

### Light theme rebalancing

The OKLCH palette is dark-tuned. The light theme inherits the same chroma
values but at lower lightness, which is too saturated for paper-white
surfaces. Re-tune chroma per theme when light becomes more than a
curiosity.

---

## Style / design

### Print / export layout

A dedicated, chrome-free layout for screen capture or PDF export of a
single plot card with annotations baked in. Useful for putting a figure
into a paper or talk.

### Per-sample cover image

A thumbnail of the trace plot, used in the nav modal's sample list and
(eventually) the Compare page. Generated server-side or rendered to
canvas client-side and cached.
