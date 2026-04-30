# Future Feature Ideas

Ideas that are intentionally out of scope for current development but worth
preserving for later planning. Each entry says *what* and *why later* — not
a full design.

For the design philosophy and the current iteration's choices, see
[`himalayaui-design.md`](himalayaui-design.md).

---

## Analysis engine

### Extended lattice types

Monoclinic (`a, b, c, β`) and tetragonal (`a, c`) lattice indexing. These
require fitting a lattice parameter *vector* rather than a single basis,
which changes how `indexpeaks`, `Index`, `score`, and `fit` work internally.
Design the extension points when the need is concrete; notes below are the
starting point for that work.

**Why the current code can't do this as-is.** The whole engine rests on one
load-bearing assumption: a phase is a 1D ratio series scaled by a single
scalar `basis`. Cubic, hexagonal, square, and lamellar all have one free
lattice parameter, so peak ratios are constants stored in
[src/phase.jl:62-75](src/phase.jl:62). Tetragonal ratios depend on `c/a`;
monoclinic ratios depend on three ratios and an angle. The vectorized
matcher at [src/index.jl:134-148](src/index.jl:134) — `observed_ratios = X * (1 ./ B)` —
collapses search to 1D *only* because the basis is a scalar. With two
parameters, no scalar division turns a peak into a comparable ratio.

**Recommended approach: dispatch-driven, not parallel implementations.**
Add a parameter-count axis to the `Phase` hierarchy:

```julia
abstract type OneParamPhase  <: Phase end
abstract type TwoParamPhase  <: Phase end
abstract type FourParamPhase <: Phase end

abstract type Cubic     <: OneParamPhase end
abstract type Lamellar  <: OneParamPhase end
abstract type Hexagonal <: OneParamPhase end
abstract type Square    <: OneParamPhase end
abstract type Tetragonal <: TwoParamPhase end
abstract type Monoclinic <: FourParamPhase end
```

Every existing `where {P<:Phase}` method still matches — nothing breaks.
The new abstract layer is a dispatch axis disguised as a hierarchy.

Then split the matcher by trait:

```julia
function indexpeaks(::Type{P}, peaks, sharpness, domain, tol, requiremin) where {P<:OneParamPhase}
    # exactly today's body — fast vectorized X * (1./B) path
end

function indexpeaks(::Type{P}, peaks, sharpness, domain, tol, requiremin) where {P<:TwoParamPhase}
    # subset-fitting body: pick peak pairs, assign trial (hkl)₁,(hkl)₂,
    # solve linearly for (a, c) on q² = A(h²+k²) + B·l², score predictions
end
```

The top-level enumerator at [src/index.jl:98-100](src/index.jl:98) just adds
the new phases to the loop — each phase routes itself to the correct method
via dispatch. Cubic/hex/lamellar pay zero added cost.

The same pattern applies to `predictpeaks`, `fit`, `score`'s coverage weights,
`remove_subsets`, and `==` — one method per parameter-count class instead of
runtime branches. The hexagonal `λ = 2/√3` special case at
[src/index.jl:238](src/index.jl:238) is a foreshadowing of this exact pattern.

**The `Index` struct change is unavoidable.**
[src/index.jl:7-11](src/index.jl:7) declares `basis::Real` concretely.
Generalize once: `basis::Real` → `params::Tuple` (or `NTuple{N,Real}` keyed
on phase). Provide `basis(idx) = idx.params[1]` as a back-compat accessor —
the codebase already uses `basis(index)` everywhere it reads the value
([src/index.jl:23](src/index.jl:23)), so most call sites are insulated.

**DB schema bump.** `indices.basis FLOAT` → either `indices.params TEXT`
(JSON tuple) or split into typed columns (`lattice_a`, `lattice_b`,
`lattice_c`, `lattice_beta`). One-time migration: scalar → `[scalar]`.
Schema lives in [packages/HimalayaUI/src/db.jl](packages/HimalayaUI/src/db.jl).

**Prerequisite: predicate-based selection rules.** Both the one-param and
multi-param matchers need to enumerate allowed `(hkl)` per phase. Today's
cubic phases get this for free because `phaseratios` is the pre-enumerated
result. Tetragonal/monoclinic need explicit `is_allowed(::Type{P}, h, k, l)`
predicates. This is the deferred refactor described in **Predicate-based
phase ratios** below — it becomes a hard prerequisite when this work is
picked up. Doing the predicate refactor first is the natural opening move.

**What dispatch does not solve.**

- *Algorithmic cost.* The `TwoParamPhase` matcher does peak-pair enumeration
  × HKL-pair assignment. Combinatorial work; the language can't help.
- *Score semantics.* `score`'s coverage weights peaks by ordinal rank
  `1/r` in the canonical ratio list. With parameter-dependent ordering,
  "rank" stops being canonical — define `peak_rank_weights(::Type{P}, params)`
  per phase. Cubic returns today's weights; tetragonal computes them from
  the q-ordering induced by `(a, c)`. The dispatch hook is clean; the
  modeling decision still has to be made.

**End-state.** `indexpeaks(peaks, sharpness)` keeps the same signature, the
same outer loop, the same `remove_subsets` post-processing. Cubic/hex/
lamellar retain the vectorized `X * (1./B)` matcher byte-for-byte. New
phases plug in by adding (a) an HKL-allowed predicate, (b) a
`predict_qs(P, params)` method, (c) participating in the multi-param
`indexpeaks` method. The `Index` struct generalizes once. The DB grows a
column or a JSON blob. Score and fit get one extra method per
parameter-count class.

### Predicate-based phase ratios with on-demand extension

Today `phaseratios(P)` returns a hardcoded `Vector{Float64}` truncated by hand
at the practical detection ceiling for each phase
([src/phase.jl:62-75](src/phase.jl:62)). Proposed shape: a new method
`phaseratios(P, n)` that returns the first `n` allowed ratios — slicing the
hardcoded prefix when `n` is small, generating additional ratios from
selection-rule predicates when `n` exceeds the prefix.

Not done today because: real SAXS data rarely shows more peaks than the
truncation; the hand-curated arrays are directly auditable against the
*International Tables for Crystallography*; and there's no observed pressure.
The static lists are also a useful provenance anchor — every value has a
citation, no value came from "the predicate, probably."

Revisit when: a 6th+ phase is added (predicate scaling beats hand-tabulation
at that point); data routinely shows peaks past the current truncation; or
someone gets a reflection condition wrong because the predicate isn't visible
in the source.

Deeper variant worth considering at that point: encode the selection rules
themselves as per-phase predicates (e.g. a `CubicSelectionRules` struct
capturing body-centering, axial conditions, `0kl` rule, `hhl` rule), make
those the single source of truth, and delete the hardcoded arrays once a test
confirms the generator reproduces them exactly. Resist doing this preemptively.

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
