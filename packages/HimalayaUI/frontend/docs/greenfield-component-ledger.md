# Greenfield "The Print" — Component Ledger

**Created:** 2026-06-02 · **Branch:** `worktree-greenfield-ui-rebuild` (unmerged) · **Status:** living doc

## What this is

A complete, layer-by-layer inventory of every component the 7 redesign mockups
(`docs/redesign-mockups/*.html`, repo root) demand, with build status. It exists so the
**bottom-up** build draws from a checklist instead of re-discovering scope batch-by-batch.

It was built from a due-diligence sweep of all 7 mockups (2026-06-02) cross-referenced
against the built layers in `src/print/`. The earlier roadmap was a bullet list and missed
several components that don't decompose cleanly into "renderer + primitives" — those are
flagged ⚠️ below.

## How to use it

- **Before starting a batch:** pick gaps (⬜) from a layer, derive each spec from the cited
  mockup region, build bottom-up (a component's dependencies must be ✅ before you build it).
- **When a component lands:** flip ⬜ → ✅ here in the same commit, and note its file path.
- **When scope changes:** record the decision in the *Out-of-scope & deferred registry* so it
  isn't re-litigated.
- This is **not** a build-order lock. Order is decided per batch; the dependency arrows only
  say what must exist *first*, not *when*.

## Status legend

| Mark | Meaning |
|---|---|
| ✅ | Built (file exists, has stories + tests) |
| ⬜ | Gap — not built |
| ⛔ | Out of scope — decided, do not build |
| ⏸ | Deferred — intentionally parked, may revisit |
| ⚠️ | Was missing from the original roadmap bullet-list (late-discovery risk) |

---

## Layer 0 — Primitives (`src/print/ui/`, design-guard-exempt)

**42 built — the foundation is complete.** Appearance lives here; consumers pass placement only.

`Badge` · `BonnetBadge` · `Button` · `Card` · `CheckCircle` · `Chip` · `Dot` · `EmptyState` ·
`FacetChip` · `FilterChip` · `GripHandle` · `HintText` · `IconButton` · `Input` · `KbKey` ·
`KbLegend` · `Kicker` · `Menu` · `MetaList` · `ModalShell` · `NoticePill` · `PeakGlyph` · `PhaseChip` ·
`PhaseLabel` · `PhaseStrip` · `ProgressBar` · `RejectOverlay` · `ScoreBar` · `SearchInput` ·
`SegmentedControl` · `SignalBars` · `Slider` · `StageTabs` · `Swatch` · `TagEditor` · `TagList` ·
`TagPill` · `Toast` · `ToggleSwitch` · `Tooltip` · `TopBar` · `Wordmark` — all ✅

> Batch 5 (2026-06-02) added **`NoticePill`** as a refactor-on-contact primitive (folio
> `SeriesCard` kick-row pill): `tone: "new" | "draft"` — `new` is the accent-tinted "+N new
> match" recipe signal (`color-mix(in oklab, var(--color-accent) 14%, transparent)`), `draft`
> is the dashed-faint marker. It authors the tint/dashed appearance the placement-only card
> can't, so it lives in `ui/`. `data-testid="notice-pill"`, `data-tone`.

> Batch 3 (2026-06-02) added two refactor-on-contact primitives so phase color stays out of the
> placement-only composite layer (mirrors `ScoreBar`/`PhaseChip`): **`Swatch`** (phase-colored 9px
> color sample, `shape: "square"|"circle"`, `data-swatch`/`data-phase`/`data-shape`) and
> **`PhaseLabel`** (phase-colored text wrapper, placement `className` sizes it). `Dot` was left
> untouched (its `tone` enum stays status-pure; a circular color sample is a `Swatch`).

**Primitive-level watch list** (decide at refactor-on-contact, may need a new primitive or a variant):

| Need | Mockup source | Likely resolution | Status |
|---|---|---|---|
| `Field` (label + value + chevron dropdown affordance) | series rail "Ordering variable" / scoping `OrderField` | new small primitive, or compose `Menu` trigger | ⬜ decide on contact |
| Number-input variant of `Input` (lattice param spinner) | focus custom-index modal | ✅ **resolved (Batch 8)**: `Input` gained `mono?` + already forwards `type="number"`/`min`/`max`/`step` via `...rest`; `SegmentedControl` gained `stretch`; `Slider` gained `ariaLabel` (bare accessible slider) | ✅ |
| `UnsetPlaceholder` ("○ Not indexed") | sample-table status cell | compose `Dot` + text — no new primitive | n/a |

---

## Layer 1 — Renderers (SVG cores: own pixels + interactions)

| Component | What | Composes / engine | Mockup | Status | File |
|---|---|---|---|---|---|
| `TracePlot` | 1D integration trace hero (peaks, axis, log/linear) | d3 projection + marks | focus-plot | ✅ | `plot/TracePlot.tsx` |
| `Axis` / `PlotFrame` | shared axis + frame chrome | — | focus-plot | ✅ | `plot/` |
| `DetectorImage` | 2D Debye–Scherrer rings | SVG + LUT | focus-workspace, sample-table | ✅ | `detector/DetectorImage.tsx` |
| `CombChart` / `ResidualChart` | reflection comb + indexing-space residual | SVG | focus-plot (combs panel) | ✅ | `comb/` |
| `WaterfallChart` | stacked-offset multi-trace hero (+ peak-q cursor) | composes N `TracePlot` | series-plot, series-builder | ✅ | `waterfall/WaterfallChart.tsx` |
| ⚠️ **`Sparkline`** | tiny 76×28 trace thumbnail | small `TracePlot` config or micro-renderer | series-scoping (`SampleRow`), series-builder (`TraceList`) | ✅ | `plot/Sparkline.tsx` |
| `CardFigure` (mini-waterfall) | frozen mini-waterfall inside a series card | frozen per-row `TracePlot` stack (reuses the engine) | series-folio (`SeriesCard`) | ✅ | `waterfall/CardFigure.tsx` |
| `CustomPreview` viz | static predicted-teeth-vs-observed-peaks render in custom-index modal | hand-rolled 1-row SVG + `makeAxis` (NOT `CombChart`) | focus-plot (custom-index modal) | ✅ | `comb/CustomPreview.tsx` |
| `CleanFigure` | export idiom (plain GraphPad render, no Print branding) | standalone SVG (Arial, plain hex); reuses only `makeAxis`/`axisTicks` math | series-plot (export sheet) | ✅ | `export/CleanFigure.tsx` |
| ⛔ `HeatmapChart` | row-per-sample binned heatmap (Waterfall⇄Heatmap toggle) | — | series-builder/plot | ⛔ **out of scope (2026-06-02)** | — |
| ⏸ `TrackingLine` / cross-trace migration | polyline tracking one reflection across the stack + ghost rings | overlay on `WaterfallChart` | series-builder ("Track reflections") | ⏸ **deferred** (declined in waterfall honing — buildable client-side, but readability call) | — |

---

## Layer 2 — Tier-1 composites (`src/print/components/`, primitives only)

| Component | What | Composes | Mockup | Status | File |
|---|---|---|---|---|---|
| `Stepper` | topbar sample nav | `IconButton`×2 + label | focus | ✅ | `Stepper.tsx` |
| `PlateHeader` | plate title block (kicker+title+sub + tools slot) | `Kicker` + text | focus, series | ✅ | `PlateHeader.tsx` |
| `ToolBar` | plot controls (scale toggle + buttons) | `SegmentedControl` + `Button` | focus, series | ✅ | `ToolBar.tsx` |
| `CandidateRow` | phase candidate (name+why+score+bonnet+add) | `ScoreBar`+`BonnetBadge`+`Button`+`PhaseChip` | focus, scoping | ✅ | `CandidateRow.tsx` |
| `CombLegend` | comb vocabulary key | `PeakGlyph` + text | focus | ✅ | `CombLegend.tsx` |
| `RepresentativeBox` | set-representative box | `Button` + caption | loupe | ✅ | `RepresentativeBox.tsx` |
| `Verdict` | kept/dropped state row | `Dot` + `Button` | loupe | ✅ | `Verdict.tsx` |
| ★ **`Thumbnail`** | mini `DetectorImage` + frame-no + rep-dot + reject-X | `DetectorImage`+`Dot`+`RejectOverlay`+text | focus, sample-table, loupe (**3 sites**) | ✅ | `components/Thumbnail.tsx` |
| ★ **`ThumbnailGallery`** | filmstrip wrapper (exposure switcher / loupe strip) | `Thumbnail`×N | focus DetectorPanel, sample-table row, loupe | ✅ | `components/ThumbnailGallery.tsx` |
| `PhaseBlock` | assigned-phase line (name+score+bar+remove) | `PhaseLabel`+`ScoreBar`+`IconButton` | focus assignment cart | ✅ | `PhaseBlock.tsx` |
| `FolioHeader` | folio page header (kicker+title+sub+count) | `Kicker`+text | series-folio | ✅ | `FolioHeader.tsx` |
| `AutoGroup` | confident-expert summary box (★ + text + actions) | star svg+text+link-buttons | scoping, builder | ✅ | `AutoGroup.tsx` |
| `MemberRow` | trace-list row (grip+swatch+name+phase) — **builder `.trow` only** | `GripHandle`+`Swatch`+`PhaseChip` | series-builder | ✅ | `MemberRow.tsx` |
| `ReadingRow` | phases-present line (dot+name+span+lattice) | `Swatch`+`PhaseLabel`+text | series-plot (reading panel) | ✅ | `ReadingRow.tsx` |
| `SeriesMemberRow` | series-**plot** `.member` row (gradient/dashed/solid swatch + per-phase colored names + variable + lattice·q₁ data line) — the split-out plot counterpart to builder `MemberRow` (Batch 9) | `Swatch`(coexist/empty)+`PhaseLabel`+text | series-plot (member list) | ✅ | `SeriesMemberRow.tsx` |
| `KeptCell` / `TagsCell` / `StatusCell` / `SpecCell` | sample-table cell composites | `TagList`/`PhaseChip`/`CheckCircle`/`Dot`/text | sample-table | ✅ | `components/{SpecCell,KeptCell,StatusCell}.tsx` (TagsCell = TagList primitive, no new file) |
| `RailSection` | rail section (faint `Kicker` label + count + children + note) | `Kicker`(faint)+text | focus, series-builder | ✅ | `components/RailSection.tsx` (Batch 8; reused by `BuilderRail` later) |
| `ModalHead` | generic modal header (accent `Kicker` + serif title + × close) | `Kicker`(accent)+`IconButton` | focus modal, export | ✅ | `components/ModalHead.tsx` (Batch 8; reused by `ExportSheet`) |
| `ModalFieldRow` | labeled modal control row (`text-label` 100px + control slot + mono `labelSuffix`) | text+slot | focus modal | ✅ | `components/ModalFieldRow.tsx` (Batch 8) |
| `ModalFooter` | generic modal footer (note slot + actions slot) | slots | focus modal, export | ✅ | `components/ModalFooter.tsx` (Batch 8; reused by `ExportSheet`) |
| `FitMetadata` | "Lands on N of M peaks · {a}={v} Å[· snapped]" line | text | focus modal | ✅ | `components/FitMetadata.tsx` (Batch 8) |
| `LatticeParamControl` | bare `Slider` + mono number `Input` + unit | `Slider`+`Input`(mono) | focus modal | ✅ | `components/LatticeParamControl.tsx` (Batch 8) |

---

## Layer 3 — Tier-2 panels (renderer + tier-1 composites)

| Component | What | Composes | Mockup | Status |
|---|---|---|---|---|
| `TracePlate` | Focus hero panel (lifted plate: header + toolbar + trace) — **consumer-owned scale toggle + armed +Peak + scroll-zoom** — `components/TracePlate.tsx` (Batch 7, +7.1 `xDomain` scroll-zoom pass-through). Legend lives in `CombsPanel`, not here (per mockup `.combs-legend`). | `Card`(elevated)+`PlateHeader`+`ToolBar`+`TracePlot` | focus-plot | ✅ |
| `DetectorPanel` | detector frame (`bg-frame-edge` aspect-square) + **phase-ring overlay** + header exposure strip + hint — `components/DetectorPanel.tsx` (Batch 7, +7.1 `DetectorRings` overlay + `hoveredQ`/`onHoverQ` triple-link). Presentational: `src`/`rings`/`calibration`/exposure-switcher passed in. | `Card`+`PanelHeader`+`DetectorImage`(full)+`DetectorRings`+`ThumbnailGallery`(xs) | focus-plot | ✅ |
| `CombsPanel` | comb/indexing-space view toggle + legend footer — `components/CombsPanel.tsx` (Batch 7). Presentational: `view`/`hoveredQ` props. Uses peak-glyph `CombLegend` (comb-tooth legend vocab deferred). | `Card`+`PanelHeader`+`CombChart`/`ResidualChart`+`SegmentedControl`+`CombLegend` | focus-plot | ✅ |
| `PanelHeader` | shared `.panel-head` (faint `Kicker` label + optional right tools slot) — `components/PanelHeader.tsx` (Batch 7). Reused by `DetectorPanel`/`CombsPanel`. | `Kicker`(faint) | focus-plot | ✅ |
| `AssignmentCart` | phase-call plate (coexistence tag when >1 · full-bleed dividers · empty state) — `components/AssignmentCart.tsx` (Batch 6; **Batch 8** added the `.as-foot` "+ custom index…" trigger via optional `onCustomIndex`, shown in both empty + filled states). | `PhaseBlock`×N | focus-workspace | ✅ |
| `AssignmentRail` | the full Focus assignment right-sidebar (`<aside class="rail">` → two `RailSection`s over slot-passed cart + candidate-list) — `components/AssignmentRail.tsx` (Batch 8). Presentational: `assignment`/`candidates` are slots, counts/notes are props. | `RailSection`×2 + `AssignmentCart` + `CandidateList` | focus-workspace | ✅ |
| `CandidateList` | ranked candidate list (thin `flex-col gap-2` wrapper) — `components/CandidateRow.tsx` (co-built with `CandidateRow`, Batch 1) | `CandidateRow`×N | focus | ✅ |
| `SheetTable` | contact-sheet plate (elevated `Card` + aligned 5-col header + slotted rows + empty state) — `components/SheetTable.tsx` (Batch 11). Children-slotting (Gallery pattern); header reuses the exported `SAMPLE_TABLE_COLS` so columns align with every row. | `Card`+`Kicker`×5 + `SampleTableRow`×N (slotted) | sample-table | ✅ | `SheetTable.tsx` |
| ★ **`SampleTableRow`** | one sample row (screened+spec+exposures+kept+tags+status) — `components/SampleTableRow.tsx` | `CheckCircle`+`ThumbnailGallery`+cell composites | sample-table | ✅ |
| `BigFrame` | loupe large detector frame (+ dropped tag, big reject-X) — `components/BigFrame.tsx` (Batch 10). Square frame, `bg-accent`/`text-plate` dropped pill, `text-frame-tag` caption, dims image + shows `RejectOverlay` when `rejected`. | `DetectorImage`+`RejectOverlay`+text | loupe | ✅ | `BigFrame.tsx` |
| `LoupeSidePanel` | loupe aside (metadata + verdict + rep box + tags + keys) — `components/LoupeSidePanel.tsx` (Batch 10). Presentational; exports `LOUPE_KEYS` default; section headers = `Kicker tone="faint"`. | `MetaList`+`Verdict`+`RepresentativeBox`+`TagList`+`KbLegend` | loupe | ✅ | `LoupeSidePanel.tsx` |
| `SeriesCard` | folio card (mini-waterfall + meta + phase strip + pill) | `CardFigure`+`PhaseStrip`+`NoticePill`+text | series-folio | ✅ `components/SeriesCard.tsx` |
| `Gallery` | masonry wall of series cards | `SeriesCard`×N (CSS multi-column) | series-folio | ✅ `components/Gallery.tsx` |
| `ScopePlate` | scoping worksheet (editable title + order field + sample rows + preview) | `PlateHeader`+scoping-`SampleRow`×N+`PhaseStrip` | series-scoping | ⬜ |
| scoping `SampleRow` | ordered member row (grip+sparkline+name+parsed value+flag) | `GripHandle`+`Sparkline`+`FlagButton`+text | series-scoping | ⬜ |
| `MemberList` | series-**plot** reading member list (data-driven + controlled hover) — `components/MemberList.tsx` (Batch 9). Composes `SeriesMemberRow` (NOT builder `MemberRow`); builder reorderable list deferred to `BuilderRail`. | `SeriesMemberRow`×N | series-plot | ✅ | `MemberList.tsx` |
| `SeriesPlate` | series-plot figure plate (header + scale toggle + waterfall + foot legend) — `components/SeriesPlate.tsx` (Batch 9). Mirrors `TracePlate`. Heatmap rep-toggle **omitted** (HeatmapChart out-of-scope); vertical phase-strip companion **deferred** (needs L1 renderer aligned to WaterfallChart row geometry). | `PlateHeader`+`ToolBar`+`WaterfallChart`+`Swatch` | series-plot | ✅ | `SeriesPlate.tsx` |
| `ReadingPanel` | derived phases-present summary (`Card` + rows + coex/form-factor notes) — `components/ReadingPanel.tsx` (Batch 9) | `Card`+`ReadingRow`×N | series-plot | ✅ | `ReadingPanel.tsx` |
| `BuilderRail` | builder sidebar (ordering/representation/display/traces sections) | `RailSection`×N + `Field`/`Slider`/`MemberList` | series-builder | ⬜ |

---

## Modals & overlays (hidden-by-default — highest late-discovery risk)

`ModalShell` primitive ✅ provides the scrim + shell. The **content** of each modal below is a gap.

### ★ `CustomIndexModal` — focus-plot (latest focus mockup), fully decomposed

Speculative-phase entry sheet. Triggered from the assignment cart's "custom index…" footer.

**✅ BUILT in Batch 8** (2026-06-02) — `components/CustomIndexModal.tsx`, fully presentational (symmetry/param/preview/fit + handlers are all props; the `SYMS` table + snap math live in the future Focus page).

| Sub-component | What | Composes | Status |
|---|---|---|---|
| `CustomIndexModal` (shell) | scrim + sheet container | `ModalShell` | ✅ `components/CustomIndexModal.tsx` |
| `ModalHead` | kicker + title + close button | `Kicker`(accent)+`IconButton` | ✅ `components/ModalHead.tsx` |
| `SymmetrySelector` | pick crystal symmetry (Pn3m/Im3m/Ia3d/Lamellar/Hexagonal) | `SegmentedControl`(`stretch`) | ✅ **folded into `SegmentedControl`** — no separate file (it carries no logic the control doesn't) |
| `ModalFieldRow` | labeled control row (the `.cs-row` label+control) | text+slot | ✅ `components/ModalFieldRow.tsx` |
| `LatticeParamControl` | dual control: slider + number field + unit | `Slider`+`Input`(mono,number)+text | ✅ `components/LatticeParamControl.tsx` |
| `CustomPreview` viz | live predicted-teeth vs observed-peaks render | hand-rolled SVG (Batch 4) | ✅ `comb/CustomPreview.tsx` |
| `FitMetadata` | "landed on N peaks · a = … · snapped" line | text | ✅ `components/FitMetadata.tsx` |
| `ModalFooter` | note + Cancel / Add-to-Assignment actions | `Button`×2 | ✅ `components/ModalFooter.tsx` |

### Other overlays

| Component | What | Composes | Mockup | Status |
|---|---|---|---|---|
| ⚠️ **`ExportSheet`** | full-screen export dialog wrapping the clean render | `ModalShell`+`ModalHead`+`CleanFigure`+`ModalFooter` | series-plot | ⬜ |
| **`CullBar`** | floating batch-reject action bar (bottom-center; visible on multi-select) — `components/CullBar.tsx` (Batch 11). Presentational (`count`/`show`/handlers); dark `bg-ink` pill, `accent` Drop + `ghostInverse` Restore/Clear; `aria-hidden`+`tabIndex=-1` when hidden. | count + `Button`×3 (`accent`/`ghostInverse`) | sample-table | ✅ | `CullBar.tsx` |
| ⚠️ **`RailBack`** | floating tab to restore collapsed rail (full-bleed mode) | `IconButton`+label | series-builder | ⬜ |
| ⚠️ **`Dock`** | floating offset-slider pod (full-bleed mode) | `Slider`+label | series-builder | ⬜ |
| `ConflictModal` | member-conflict resolution (exists in legacy Series) | `ModalShell`+content | series-builder | ⬜ verify vs mockups |

---

## Layer 4 — Page shells (`src/print/pages/`, tier-3)

All ⬜ — deferred until composite + renderer layers are proven in Storybook (bottom-up).

| Page | Mockup(s) | Panels it assembles |
|---|---|---|
| Focus workspace | `focus-workspace.html`, `2026-05-29-focus-plot.html` | `TracePlate` · `DetectorPanel` · `CombsPanel` · `AssignmentCart` · `CandidateList` · `CustomIndexModal` |
| Samples / corpus | `sample-table.html` | `SheetTable` · `CullBar` · (loupe: `BigFrame` · `LoupeSidePanel` · `ThumbnailGallery`) |
| Series folio | `series-folio.html` | `FolioHeader` · `Gallery`(`SeriesCard`) · `SearchInput` · `FilterChip`s · `SegmentedControl` |
| Series scoping | `series-scoping.html` | `ScopePlate` · scoping `SampleRow`s · `AutoGroup` · `CandidateRow`s |
| Series builder | `series-builder.html`, `2026-05-29-series-plot.html` | `SeriesPlate`(`WaterfallChart`) · `BuilderRail` · `MemberList` · `ReadingPanel` · `ExportSheet` · `RailBack`/`Dock` |

---

## Out-of-scope & deferred registry (decisions, so they aren't re-litigated)

| Item | Decision | Date | Rationale |
|---|---|---|---|
| `HeatmapChart` (Waterfall⇄Heatmap toggle) | ⛔ **out of scope** | 2026-06-02 | Not needed for the rebuild; waterfall is the canonical Series representation. |
| `TrackingLine` / cross-trace migration track + ghost rings | ⏸ **deferred** | 2026-06-02 | Buildable client-side from `confirmed_index.phase + lattice_d` (legacy `buildAnchorMap` proves it), but q-proximity matching is fragile and the connector is hard to read. The peak-q cursor delivers ~90% of the reading value. Revisit only if users ask. |
| `MemberRow` series-plot `.member` variant (gradient coexistence swatch + inline colored phase names + lattice data) | ✅ **split done (Batch 9)** | 2026-06-02 | RESOLVED: built as the new `SeriesMemberRow` leaf (not a `MemberRow` variant — structurally distinct: no grip/chip; gradient/dashed/solid swatch + per-phase colored names + var + lattice·q₁ data line). `Swatch` gained the promised `coexistWith` gradient + `empty` (form-factor dashed) + `size` modes (the `Swatch.tsx` stub). Builder `MemberRow` (Batch 3) stays the `.trow`. |
| Batch 10 — Loupe slice (`BigFrame` + `LoupeSidePanel`) | ✅ **loupe surface built** | 2026-06-02 | The loupe view's two composites, derived from `sample-table.html` `.loupe-shell`. **Pure composition — ZERO new primitives / no `ui/` refactor-on-contact:** every appearance detail traced to an existing token or primitive (the `.frame-tag` caption uses the `--color-frame-tag` token added back in R0c specifically for this; the solid "Dropped" pill uses `bg-accent`+`text-plate` token classes, NOT `NoticePill` which is accent-*tinted*; the big ✕ reuses `RejectOverlay`). `BigFrame` = square `bg-frame-edge` frame, `DetectorImage size="full"` dimmed `opacity-40` when rejected, `data-rejected` flag. `LoupeSidePanel` is **presentational** (Batch 7 contract): no `useState`, all state lifted; `shortcuts` defaults to exported `LOUPE_KEYS`; the page-sim assembly (`LoupeAssembly.stories.tsx`) owns selected-frame/dropped-set/rep state. **Layer-4 deferral (consistent w/ Batch 7/9):** the `.loupe-plate` shell (bordered plate + `loupe-back` + serif `loupe-head` + 2-col body) + `loupe-strip` are page assembly — a combined `LoupeAssembly` story simulates the page (2-col grid + `ThumbnailGallery` filmstrip) for fidelity without building the page component. frontend-reviewer clean on both (no issues). Gate green (lint:design + tsc + build-storybook); BigFrame (kept+dropped), LoupeSidePanel, and the assembly all visually verified vs the mockup. |
| Batch 9 — Series reading slice (`SeriesMemberRow`/`MemberList`/`ReadingPanel`/`SeriesPlate`) | ✅ **series-plot reading surface built** | 2026-06-02 | The series-plot rail panels + figure plate, derived from `2026-05-29-series-plot.html`. **`SeriesPlate` mirrors `TracePlate`** (`Card elevated`+`PlateHeader`+`ToolBar`+`WaterfallChart`+foot). **Two deliberate omissions (controls don't lie / no unbuilt deps):** the waterfall⇄heatmap rep-toggle is **omitted** (only the log/linear scale `SegmentedControl` ships — `HeatmapChart` is out-of-scope per the row above); the vertical `#phasestrip` companion is **deferred** — it needs an L1 renderer aligned to `WaterfallChart`'s internal (measured, bandHeight-weighted) row geometry, which the chart doesn't expose, so it's a separate slice (either expose row y-positions or render the strip inside the chart). **Panels presentational** (Batch 7 contract): `MemberList` hover (`hoveredKey`/`onHoverMember`) + `SeriesPlate` scale/hover are page-owned; stories simulate the page. `MemberList` is the series-plot reading list (composes `SeriesMemberRow`); the builder reorderable list is still `BuilderRail`'s. frontend-reviewer clean (no issues). Full gate green (518 print tests / 71 files, lint:design + tsc clean, build-storybook OK); all four surfaces visually verified vs the mockup. |
| Type-scale snaps: `FolioHeader` title 31px→`text-display` (27); `PhaseBlock` name 22px→`text-headline` (19) | ✅ **snap to named role** | 2026-06-02 | The mockups use off-scale serif sizes (31/22px) with no named role; `text-[31px]` is guard-banned. Composites snap to the **nearest named role** rather than add a one-off scale step. If fidelity review objects, the fix is a new `--text-*` step in `styles.css` (the authoring layer), never an arbitrary in a composite. Folio count 26px = `text-headline-lg` is exact. |
| `NoticePill` is a NEW primitive, not a `Chip`/`Badge` variant (Batch 5) | ✅ **new ui/ primitive** | 2026-06-02 | The folio card-kick pill (`+N new match` accent-tint / dashed `Draft`) is a distinct 10px tinted/dashed look — `Chip` is a full-border pill, `Badge` a flat mono count, neither fits. The accent tint (`color-mix(…var(--color-accent) 14%…)`) + dashed border are appearance the placement-only `SeriesCard` cannot author, so they live in `ui/` (guard-exempt). No guard-allowlist edit needed — `ui/` was already exempt. `SeriesCard` stays placement-only by composing it. |
| `SeriesCard` figure region sits on PLATE, not paper-sunk (Batch 5 review fix) | ✅ **plate surface + paper-sunk hairline** | 2026-06-02 | First pass painted `bg-paper-sunk` across the whole figure pad; the mockup `.card-fig` keeps the pad on the card's plate surface with only a 1px paper-sunk bottom line — and `CardFigure`'s inner plots use `paperColor="var(--color-plate)"`, so a sunk pad created a plate-on-sunk seam. Fixed to `border-b border-paper-sunk` (the `--color-paper-sunk` Tailwind v4 `@theme` token generates `border-paper-sunk`, guard-clean). `cursor-pointer` also made conditional on `onClick`. |
| `src/print/export/` added to the design-guard `isExcluded` allowlist (Batch 4) | ✅ **new exempt renderer dir** | 2026-06-02 | `CleanFigure` is the export idiom — it DELIBERATELY sheds the Print look (Arial + literal hex `#c8841f`/`#4a4ba8`/`#111`, no `var(--color-*)` tokens), so it must author literal appearance like the other renderer dirs (`plot`/`comb`/`detector`). Edit is a single trailing-slash-bounded prefix in `scripts/check-design.mjs:isExcluded` — proven to still flag `print/components/` and not over-match `print/exporters-*`. Mirrors the legacy `lib/figure-export/**` export-palette allowlist. `CardFigure`/`CustomPreview` did NOT need exemption (literal-free, `phaseColor()` + tokens only). |
| `AssignmentCart` owns count→tag + positional-divider + empty (Batch 6) | ✅ **container owns cross-child state** | 2026-06-02 | The mockup `.phasecall` conditional `.pc-tag` ("Coexistence · N phases", only when >1 phase) and the `.pc-block + .pc-block` divider are state a lone `PhaseBlock` can't know — they belong to the container. `AssignmentCart` derives count from `Children.count`, renders the tag at count>1, and uses the Gallery-proven `Children.map` + positional `border-t border-hair` on a full-width `px-4` wrapper for **full-bleed** dividers (NOT `divide-y` on a padded wrapper — that insets the hairline). Wrappers carry `data-testid="cart-block" data-divider` so the load-bearing divider is asserted without class-string tests. `CandidateList` was already built (Batch 1, in `CandidateRow.tsx`) so this slice was just `AssignmentCart` + an optional content-only `series` prop on `PhaseBlock` (`.pcb-series`, guard-clean: `font-mono`/`text-caption`/`text-ink-*` only). |
| Focus panels are PRESENTATIONAL; state lifted to the consumer (Batch 7) | ✅ **panels render props, the page owns state** | 2026-06-02 | `TracePlate` (scale + `addPeakArmed`), `DetectorPanel` (current `src` + exposure-switcher element), and `CombsPanel` (`view` + `hoveredQ`) hold NO internal `useState` — every interaction datum is a prop with a consumer-owned handler. This is "container owns cross-child state" pushed one level up: the cross-PANEL state (the q-hover link that lights peak↔ring↔tooth across all three) can only live in the page, so panels must not capture it. Stories simulate the page with a small `useState` demo wrapper. The chrome is the `Card` primitive: `elevated` → the lifted hero plate (`TracePlate`), flat default → the two lower panels. |
| Focus lower panels share `PanelHeader`; `CombLegend` tooth-vocab deferred (Batch 7) | ✅ **shared `.panel-head` seam + deferred legend** | 2026-06-02 | `DetectorPanel`/`CombsPanel` share `PanelHeader` (faint `Kicker` label + right tools slot) — the `.panel-head`. The `.panel-h` label = `Kicker tone="faint"` exactly (uppercase 700 tracked ink-faint). `TracePlate` does NOT use it (the hero uses the big serif `PlateHeader`). The mockup `.combs-legend` is a comb-TOOTH vocabulary (predicted&observed / predicted-absent / leftover); the built `CombLegend` is the peak-GLYPH vocab (auto/manual/predicted-absent/excluded). `CombsPanel` uses the existing `CombLegend` for now — the tooth-specific legend is a DEFERRED `CombLegend` refinement (do not block a slice on it). `Thumbnail` gained an `xs` (≈30px) size for the dense Focus exposure strip via its existing numeric inline-style `SIZE_PX` idiom (guard-clean; `Thumbnail` lives in `components/`, NOT `ui/`). `exactOptionalPropertyTypes: true` forces conditional-spread of optional props forwarded to renderers whose matching prop is not `\| undefined` (`TracePlate.interaction`, `DetectorPanel.lutVariant`, `CombsPanel.hoveredQ`/`onHoverQ`). |
| Batch 7.1 — Focus-plates interaction-completeness pass | ✅ **wired the missing panel interactions** | 2026-06-02 | An interaction audit of the slice found three panel-level gaps (all fixed): **(a) TracePlate scroll-zoom** — `TracePlot` emitted `onXDomain` but the plate forwarded no controlled `xDomain`, so the zoom never rendered (TracePlot is fully controlled, no internal zoom state); added the `xDomain` pass-through so the page round-trips wheel→`onXDomain`→`xDomain`. **(b) DetectorPanel phase rings** — the panel rendered only the grayscale `DetectorImage`; the whole point (rings in phase colour) was missing. Now overlays `DetectorRings` via `buildRingPlacements(rings, calibration ?? null)` (null cal → centered fallback; `.det-box` made `relative` so the `absolute inset-0` overlay lands; `imageAspect` stays 1 on the square frame), and exposes `hoveredQ`/`onHoverQ` for the cross-panel "triple lights up" link. **(c) Gallery overflow** — the `PanelHeader` tools slot was width-unbounded so a long exposure strip widened the header; added `min-w-0` to the tools wrapper (+ `flex-shrink-0` on the label) so an `overflow-x-auto` `ThumbnailGallery` scrolls in place. **TRACE q-LINK DEFERRED TO PAGE ASSEMBLY (user decision):** `TracePlot` emits no outward hovered-q (only internal `hoverId`), so hovering a *trace peak* can't yet drive the ring / tooth. The panels already expose `hoveredQ`/`onHoverQ`; the `TracePlot` `onHoverQ` emit + incoming-`hoveredQ` highlight is a Layer-1 engine addition to wire when the Focus *page* is built (the page owns the shared q). **Plus four bottom-up rendering-correctness fixes (commits 2147470/6a38a91/4bd055e):** (d) **DetectorPanel real beam center** — added a `beamCenter` override (measured center + presentational fallback radii — the real intermediate state) + `imageAspect` (frame box ratio so the real portrait image fills with no letterbox and rings register); the story now uses sample-37's real geometry (beam ≈ {0.430, 0.198}) like the `DetectorRings` storyboard, so rings overlay the actual diffraction arcs. (e) **DetectorPanel frame bound** — `maxWidth 320` (mockup det-size) / `minWidth 160` inline-style (geometry, guard-clean) so the frame — esp. the empty placeholder — can't blow up to fill an unbounded container. (f) **`TracePlot` trace+annotation clip** (Layer-1 engine) — trace line, ±σ band, peaks, and labels are clipped to the plot rect (unique `useId()` clipPath) so q outside the visible window (post-zoom) stops at the axes instead of overdrawing spines/ticks; the clip top is left open to the SVG edge (`y = -margins.top`) so peak labels keep their headroom — only left/right/bottom (the axis edges) bound the annotations. New `TracePlate` `Zoomed` story demonstrates it. |
| Batch 11 — Contact-sheet slice (`SheetTable` + `CullBar`) | ✅ **corpus contact-sheet surface built** | 2026-06-03 | The contact-sheet view's two composites, derived from `sample-table.html`. **One refactor-on-contact (the affirmative case):** `CullBar` is a dark `bg-ink` floating pill whose Restore/Clear must read as light-muted-on-dark — no existing `Button` variant fits the dark surface, so a `ghostInverse` variant was added in the exempt `ui/` layer (`text-paper/70 hover:text-paper` — color-channel alpha so `transition-colors` eases it and the focus ring stays full-strength; NOT element `opacity`). Reject reuses `variant="accent"`. **`SheetTable` is children-slotting** (the `Gallery` pattern — the page maps samples→`SampleTableRow`, the plate only lays them out + adds the aligned header): the header grid applies the **exported `SAMPLE_TABLE_COLS`** (NOT a re-derived `grid-cols-[…]`) so header↔row columns provably align; `Card elevated` + `overflow-hidden` clips the header border + rows to the rounded corners; labels are `Kicker tone="faint"`. **`CullBar` presentational** (Batch 7 contract): `count`/`show`/handlers are props; when hidden it's `opacity-0 pointer-events-none` AND `aria-hidden` + `tabIndex=-1` on the buttons (so the destructive Drop/Clear can't be keyboard/AT-reached on an invisible bar). **Radius deviation:** the mockup's 10px pill → `rounded-md` (5px), the system's one radius step, as `Toast` already does. **Layer-4 deferral (consistent w/ Batch 7/9/10):** the `.sheet-shell` page shell, the `.head` (Contact-sheet kicker + sub + the `6/9 screened` `ProgressBar`), `.kb-legend`, the view-seg sheet⇄loupe toggle, the Beamtime facet, and the cross-row selection-Set + keyboard wiring are page assembly — a `ContactSheetAssembly` story simulates the page (head + `ProgressBar` + `SheetTable` mapping `SampleTableRow`s + `CullBar` driven by a real `useState<Set>` selection) for fidelity without building the page component. frontend-reviewer clean on `CullBar` (one IMPORTANT fixed — hidden-bar tab order) and `SheetTable` (no issues). Gate green (lint:design + tsc + build-storybook); SheetTable, CullBar, and the assembly all visually verified vs the mockup. |
| Batch 8 — Focus assignment rail + custom-index modal (both PRESENTATIONAL) | ✅ **state lifted to the page** | 2026-06-02 | The two clusters the user named ("whole assignment side bar" + "custom index modal"), built as one coupled slice — the rail's `.as-foot` footer is what opens the modal. **`AssignmentRail`** wraps two `RailSection`s over slot-passed `AssignmentCart`+`CandidateList` (counts/notes as props). **`CustomIndexModal`** holds NO `useState` — symmetry/paramValue/previewSeries/fit + every handler are props; the `SYMS` table (per-symmetry lattice defaults/min/max), snap-to-peak math, and predicted-teeth computation stay in the future Focus page. **`SymmetrySelector` folded into `SegmentedControl`** (+ new `stretch` prop) — no separate file; it carries no logic the control doesn't. **The `.as-foot` "+ custom index…" trigger lives inside `AssignmentCart`** (deferred from Batch 6) via optional `onCustomIndex`, shown in both empty + filled states. **Three additive refactor-on-contact `ui/` extensions:** `Input.mono` (mono value font), `SegmentedControl.stretch` (full-width equal segments), `Slider.ariaLabel` (accessible bare slider, NO visible label row — the modal row label lives in `ModalFieldRow`). **`ModalHead`/`ModalFieldRow`/`ModalFooter` are GENERIC** — reused by the later `ExportSheet`. **`.text-label` is the EXACT named role for the `.cs-flabel` field label** (11.5px/500/uppercase/0.06em/ink-soft); serif title snaps 20→19px (`text-headline`). `CustomPreview` (Batch 4) wrapped in the `.cs-previewwrap` paper-sunk plate. frontend-reviewer clean (no BLOCKER/IMPORTANT; one NIT fixed — always-true min/max conditional spreads dropped in `LatticeParamControl`). Full gate green (612 print tests / 81 files, lint:design + tsc clean, build-storybook OK); both surfaces visually verified vs `focus-workspace.html` + `.custom-sheet`. With Batch 6's cart + Batch 7's plates, the **Focus page's entire panel + modal layer is now built**. |

---

## Coverage summary

- **Layer 0 primitives:** 42/42 ✅ (foundation complete; +`Swatch` +`PhaseLabel`, Batch 3; +`NoticePill`, Batch 5)
- **Layer 1 renderers:** 9 ✅ — **layer COMPLETE** (Batch 4 closed `CardFigure`/`CustomPreview`/`CleanFigure`) · 1 ⛔ · 1 ⏸
- **Layer 2 tier-1 composites:** 18 ✅ — layer complete + 1 (Batch 3 closed `PhaseBlock`/`FolioHeader`/`AutoGroup`/`MemberRow`/`ReadingRow`; **Batch 9 added `SeriesMemberRow`** — the split-out series-plot member leaf)
- **Layer 3 tier-2 panels:** ~22 ✅ (Batch 5: `SeriesCard`/`Gallery`; Batch 6: `AssignmentCart` + `CandidateList`; Batch 7: `TracePlate`/`DetectorPanel`/`CombsPanel` + `PanelHeader`; Batch 8: `AssignmentRail` + `RailSection` + the 6 modal leaves; Batch 9: `MemberList` + `ReadingPanel` + `SeriesPlate`; Batch 10: `BigFrame` + `LoupeSidePanel`; **Batch 11: `SheetTable`**) · ~2 ⬜ (`ScopePlate`+scoping `SampleRow`, `BuilderRail`)
- **Modals/overlays:** 1 cluster ✅ (`CustomIndexModal` + 6 leaves, Batch 8; `ModalShell` chrome ✅) · **`CullBar` ✅ (Batch 11)** · ~4 gap clusters ⬜ (`ExportSheet`, `RailBack`, `Dock`, `ConflictModal`)
- **Layer 4 pages:** 0/5 ⬜

Layers 0, 1, and 2 are now all complete. The build frontier is the **Layer 3 tier-2 panels**, which are now fully unblocked (every renderer + tier-1 composite they compose exists):

- ~~`AssignmentCart`←`PhaseBlock` · `CandidateList`←`CandidateRow`~~ ✅ **DONE (Batch 6)** · ~~`ReadingPanel`←`ReadingRow` · `MemberList`←`SeriesMemberRow`~~ ✅ **DONE (Batch 9)**
- ~~`Gallery`←`SeriesCard`←(`CardFigure`+`PhaseStrip`)~~ ✅ **DONE (Batch 5)** · ~~`SheetTable`←`SampleTableRow`~~ ✅ **DONE (Batch 11)**
- ~~`TracePlate` · `DetectorPanel` · `CombsPanel` (+ shared `PanelHeader`)~~ ✅ **DONE (Batch 7)** · ~~`SeriesPlate`←(`PlateHeader`+`ToolBar`+`WaterfallChart`)~~ ✅ **DONE (Batch 9)**
- ~~loupe panels (`BigFrame`, `LoupeSidePanel`)~~ ✅ **DONE (Batch 10)** · scoping (`ScopePlate`+`SampleRow`), `BuilderRail`
- modals/overlays: `CustomIndexModal`←`CustomPreview` · `ExportSheet`←`CleanFigure` · ~~`CullBar`~~ ✅ **DONE (Batch 11)** · plus `RailBack`/`Dock`/`ConflictModal`

**Batch 7 (2026-06-02) closed the Focus-plates slice** — the Focus hero trio (`TracePlate`/`DetectorPanel`/`CombsPanel`) + the shared `PanelHeader`, plus a `Thumbnail` `xs` size (refactor-on-contact). All three panels are presentational (scale/exposure/view/hovered-q lifted to the consumer). With Batch 6's assignment rail, **the entire Focus workspace's panel layer is now built — the Focus *page* (Layer 4) is assemblable.** Remaining next vertical slices, each fully unblocked:

1. ~~**Series reading** — `MemberList`←`SeriesMemberRow` + `ReadingPanel`←`ReadingRow` + `SeriesPlate`←`WaterfallChart`~~ ✅ **DONE (Batch 9)** — the series-plot reading surface (rail panels + figure plate). Remaining series-**builder** surface (`BuilderRail`, builder `MemberList` reorder, `ExportSheet`/`RailBack`/`Dock`) is its own slice.
2. ~~**Contact sheet** — `SheetTable`←`SampleTableRow` + `CullBar`~~ ✅ **DONE (Batch 11)** — the corpus contact-sheet grid + floating cull bar. The Layer-4 contact-sheet *page* (sheet-shell + head/progress + kb-legend + sheet⇄loupe toggle + selection/keyboard wiring) is assemblable. Remaining L3 frontier: **series scoping** (`ScopePlate`+scoping `SampleRow`, needs new `FlagButton`) and **series builder** (`BuilderRail`, needs new `Field` + reorderable member list).
3. ~~**Loupe** — `BigFrame`+`LoupeSidePanel`~~ ✅ **DONE (Batch 10)** — the loupe view's panels (big detector frame + aside). The Layer-4 loupe *page* (plate shell + filmstrip + back/head) is assemblable. Remaining: **scoping** (`ScopePlate`+scoping `SampleRow`).
