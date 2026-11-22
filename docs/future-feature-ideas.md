# Future Feature Ideas

Ideas that are intentionally out of scope for current development but worth
preserving for later planning. Each entry says *what* and *why later* — not
a full design.

For the product principles and design system, see [`PRODUCT.md`](../PRODUCT.md)
and [`DESIGN.md`](../DESIGN.md).

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
[src/phase.jl](../src/phase.jl). Tetragonal ratios depend on `c/a`;
monoclinic ratios depend on three ratios and an angle. The vectorized
matcher at [src/index.jl](../src/index.jl) — `observed_ratios = X * (1 ./ B)` —
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

The top-level enumerator at [src/index.jl](../src/index.jl) just adds
the new phases to the loop — each phase routes itself to the correct method
via dispatch. Cubic/hex/lamellar pay zero added cost.

The same pattern applies to `predictpeaks`, `fit`, `score`'s coverage weights,
`remove_subsets`, and `==` — one method per parameter-count class instead of
runtime branches. The hexagonal `λ = 2/√3` special case at
[src/index.jl](../src/index.jl) is a foreshadowing of this exact pattern.

**The `Index` struct change is unavoidable.**
[src/index.jl](../src/index.jl) declares `basis::Real` concretely.
Generalize once: `basis::Real` → `params::Tuple` (or `NTuple{N,Real}` keyed
on phase). Provide `basis(idx) = idx.params[1]` as a back-compat accessor —
the codebase already uses `basis(index)` everywhere it reads the value
([src/index.jl](../src/index.jl)), so most call sites are insulated.

**DB schema bump.** `indices.basis FLOAT` → either `indices.params TEXT`
(JSON tuple) or split into typed columns (`lattice_a`, `lattice_b`,
`lattice_c`, `lattice_beta`). One-time migration: scalar → `[scalar]`.
Schema lives in [packages/HimalayaUI/src/db.jl](../packages/HimalayaUI/src/db.jl).

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
([src/phase.jl](../src/phase.jl)). Proposed shape: a new method
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
backgrounds. Today the workspace subtracts a chosen buffer sample for display
only, scaled by beam monitors; nothing is subtracted in the pipeline.

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

## Workspace — annotation depth

### Peak classification on ticks

Color-code predicted-q ticks by whether they are explained / unexplained /
excluded. Currently all ticks for the active group render in the phase color;
adding a "missing peak" or "excluded peak" channel would let the user see
fit gaps at a glance.

### `hkl` labels above active-group ticks

Render the Miller indices `(hkl)` next to each predicted-q tick on the trace.
Useful for publication figures and for users learning the indexing. Visual
budget is tight, so probably hover-revealed rather than always-on.

### Color-blind accessibility on phase hues

Eight earthy phases at similar chroma are not distinguishable under all
forms of color vision. Add a redundant channel — dash patterns on ticks,
shape on Miller dots — so phase identity doesn't depend on hue alone.

---

## Peak move on the trace plot

The workspace plot lets a click add a manual peak, switch an auto peak off or
remove a manual one. There is no way to move a peak: q is immutable.

**Why later — the hard backend blocker.** Peaks are stored as two tables
([db.jl](../packages/HimalayaUI/src/db.jl)): `auto_peaks` (analysis
output, **regenerated on every reanalysis**) and `peak_curations` (the manual
delta layer, `kind IN ('exclude','add')`). Critically, **curations are
keyed by q, not by a stable peak identity**: an exclude curation suppresses
"whatever auto peak sits at q=X" (`events.jl:328`); the
effective list is `auto_peaks LEFT JOIN peak_curations` with adds unioned in.
So:
- There is **no representation for "a moved auto peak"** — an auto peak's q is
  owned by analysis. Move collapses to exclude-at-old-q + add-at-new-q: lossy
  (provenance flips to a manual `add` with `intensity = nothing`, identity and
  history lost), and manual-only (can't truly move an auto peak).
- Because exclusion is **q-matched**, a reanalysis that nudges an auto peak's
  q by epsilon can orphan its curation.

Move-done-right requires the curation model rethought first: a stable
per-peak identity that survives reanalysis, or a `shift`/`move` curation kind
that reanalysis honors. That backend design is the real prerequisite, and is
why this is one redesign (plot card + curation model), not a `PATCH /peaks/:id
{q}` bolt-on. Pick it up when the curation/reanalysis model is on the table.

---

## Cross-sample comparison

(Series and the Compare page were deleted 2026-11-21. Anything like them gets a
fresh design; the old code is only in git history.)

### Summary table

Full-screen tabular view of all samples in an experiment: confirmed phases,
lattice parameters, R², score. Sortable, filterable by label. Export to
CSV/JSON.

### Cross-experiment sample comparison

Query sample labels and analysis results across experiments to compare
samples with similar compositions (e.g., same lipid/peptide system across
beamtimes). Every experiment already lives in the one central DB, so this
is a query and a UI, not an aggregation layer.

---

## Multi-user / review

### Versioning and a durable history

Multiplayer is last-write-wins with no conflict UI (see `docs/event-log.md`
§"Conflict resolution"). Per-tab undo/redo shipped in the workspace (PR #403):
undo reads the server first and restores only what the edit changed. Still
open: versioning (automatic snapshots or named saves) and a durable,
cross-session history built from `user_actions`.

### Reviewer workflow

`reviewed` status transitions on indexings, with a reviewer attribution.
The `users` and `user_actions` tables exist; the UI for promotion / review
does not.

### Per-user audit view

`user_actions` rows already record who did what. A view that shows a
user's recent actions (peak edits, calls, drops, labels) would
support both review and onboarding ("show me what alice has been doing").

---

## Onboarding & docs

### Tutorial: the ⌘K palette

The three intro slides ("Add a beamtime", "One sample at a time", "Call it
and move on") cover the workspace keys and comparisons. They don't mention
the ⌘K palette, which now holds Rescan, settings, bridges and plot views.

### Empty / error state mock-ups

Default Hint/Loading text is fine for development. Real empty states ("no
frames yet — drop a `.dat` here") and error states ("analysis failed;
click to see traceback") need design.

---

## Operations / configuration

### Beamline-config editor

Each setup stores its geometry (energy, flight path, beam centre, pixel
size) on the `setups` row, set when the setup is added. There is no UI to
edit it afterwards, and no calibrant record.

### Export UI

Backend has `routes_export.jl` (JSON / CSV), which `bin/verify` uses. No UI
surface yet; a summary table (above) is the natural home for a download.

### Multi-frame overlay on the trace

Overlay several frames of one sample on the same plot with configurable
I-offsets, to see flares and frame-to-frame consistency at a glance. Compare
overlays other samples' representative frames only.

### Derived-exposure construction

The schema supports derived exposures (e.g. background-subtracted) but
there's no UI to create them. It fits with combining frames in the workspace
(per-frame include toggles, time-weighted sum or median).

### Reshoot grouping for new scans

Auto-grouping derives samples from `(load, slot)` only, so the **same
specimen re-shot in a later load becomes two samples** — `load` is part of
sample identity. In the manifest-era dev data this is frequent (~20% of
samples span multiple loads). This is true for migrated data too: the
production migration applies the same `(load, slot)` partition a fresh scan
would, so a reshot specimen splits into per-load samples there as well (it
reuses the old row for the earliest cell to keep the human name, and mints
a new sample for the later cell). Re-uniting a reshoot's cells into one
specimen is the deferred capability, and it is the *same* `merged_into_id`
action for migrated and new data. The hard part is **detecting** which
cells to suggest merging: PRP files carry no specimen name/barcode, and
filenames don't reliably share a specimen token across a reshoot's frame
ranges — the only in-data signal is same-`horizontal_position`-across-loads,
which is weak (a holder gets reloaded with different specimens at the same
slots). So this can't be auto-derived reliably; it needs careful design.

Realistic directions, cheapest first: (1) manual merge with good UX
(`merged_into_id` exists; the grouping-review merge route was deleted 2026-11-21, and today a bridge links the two samples instead); (2)
*suggested* merges from position-coincidence, human-confirmed, never auto-applied (a wrong
auto-merge fuses distinct specimens — worse than a missed one); (3) capture
specimen identity at acquisition time (an optional label/barcode), which
recovers the capability deterministically but re-adds a manifest-ish input
the rework deliberately removed. Decision deferred until we have
post-rework usage data: the 20% rate is from manifest-era acquisition, and
whether the frictionless workflow reshoots as often is unknown. If it does,
that's the evidence for (3) — no heuristic recovers a signal the
acquisition never wrote down.

---

## Frontend infrastructure

### Code splitting

Vite warns the bundle is >500 kB. Probably split: the palette +
OnboardingFlow as lazy chunks.

### Semantic spacing / density scale

There's no `@theme` spacing-token system; layouts use literal Tailwind class
strings. A semantic spacing scale would let density be tuned coherently
app-wide rather than per-surface.

---

## Style / design

### Print / export layout

A dedicated, chrome-free layout for screen capture or PDF export of a
single plot card with annotations baked in. Useful for putting a figure
into a paper or talk.
