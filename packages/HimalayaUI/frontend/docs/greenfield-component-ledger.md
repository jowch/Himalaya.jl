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
| Number-input variant of `Input` (lattice param spinner) | focus custom-index modal | extend `Input` (`type="number"` + unit suffix) | ⬜ decide on contact |
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
| `KeptCell` / `TagsCell` / `StatusCell` / `SpecCell` | sample-table cell composites | `TagList`/`PhaseChip`/`CheckCircle`/`Dot`/text | sample-table | ✅ | `components/{SpecCell,KeptCell,StatusCell}.tsx` (TagsCell = TagList primitive, no new file) |

---

## Layer 3 — Tier-2 panels (renderer + tier-1 composites)

| Component | What | Composes | Mockup | Status |
|---|---|---|---|---|
| `TracePlate` | Focus hero panel (lifted plate: header + toolbar + trace) — **consumer-owned scale toggle + armed +Peak + scroll-zoom** — `components/TracePlate.tsx` (Batch 7, +7.1 `xDomain` scroll-zoom pass-through). Legend lives in `CombsPanel`, not here (per mockup `.combs-legend`). | `Card`(elevated)+`PlateHeader`+`ToolBar`+`TracePlot` | focus-plot | ✅ |
| `DetectorPanel` | detector frame (`bg-frame-edge` aspect-square) + **phase-ring overlay** + header exposure strip + hint — `components/DetectorPanel.tsx` (Batch 7, +7.1 `DetectorRings` overlay + `hoveredQ`/`onHoverQ` triple-link). Presentational: `src`/`rings`/`calibration`/exposure-switcher passed in. | `Card`+`PanelHeader`+`DetectorImage`(full)+`DetectorRings`+`ThumbnailGallery`(xs) | focus-plot | ✅ |
| `CombsPanel` | comb/indexing-space view toggle + legend footer — `components/CombsPanel.tsx` (Batch 7). Presentational: `view`/`hoveredQ` props. Uses peak-glyph `CombLegend` (comb-tooth legend vocab deferred). | `Card`+`PanelHeader`+`CombChart`/`ResidualChart`+`SegmentedControl`+`CombLegend` | focus-plot | ✅ |
| `PanelHeader` | shared `.panel-head` (faint `Kicker` label + optional right tools slot) — `components/PanelHeader.tsx` (Batch 7). Reused by `DetectorPanel`/`CombsPanel`. | `Kicker`(faint) | focus-plot | ✅ |
| `AssignmentCart` | phase-call plate (coexistence tag when >1 · full-bleed dividers · empty state) — `components/AssignmentCart.tsx` (Batch 6). Custom-index footer deferred to the `CustomIndexModal` slice. | `PhaseBlock`×N | focus-workspace | ✅ |
| `CandidateList` | ranked candidate list (thin `flex-col gap-2` wrapper) — `components/CandidateRow.tsx` (co-built with `CandidateRow`, Batch 1) | `CandidateRow`×N | focus | ✅ |
| `SheetTable` | contact-sheet grid (header + rows) | `TableHeader`+`SampleTableRow`×N | sample-table | ⬜ |
| ★ **`SampleTableRow`** | one sample row (screened+spec+exposures+kept+tags+status) — `components/SampleTableRow.tsx` | `CheckCircle`+`ThumbnailGallery`+cell composites | sample-table | ✅ |
| `BigFrame` | loupe large detector frame (+ dropped tag, big reject-X) | `DetectorImage`+text+`RejectOverlay` | loupe | ⬜ |
| `LoupeSidePanel` | loupe aside (metadata + verdict + rep box + tags + keys) | `MetaList`+`Verdict`+`RepresentativeBox`+`TagList`+`KbLegend` | loupe | ⬜ |
| `SeriesCard` | folio card (mini-waterfall + meta + phase strip + pill) | `CardFigure`+`PhaseStrip`+`NoticePill`+text | series-folio | ✅ `components/SeriesCard.tsx` |
| `Gallery` | masonry wall of series cards | `SeriesCard`×N (CSS multi-column) | series-folio | ✅ `components/Gallery.tsx` |
| `ScopePlate` | scoping worksheet (editable title + order field + sample rows + preview) | `PlateHeader`+scoping-`SampleRow`×N+`PhaseStrip` | series-scoping | ⬜ |
| scoping `SampleRow` | ordered member row (grip+sparkline+name+parsed value+flag) | `GripHandle`+`Sparkline`+`FlagButton`+text | series-scoping | ⬜ |
| `MemberList` | builder/plot member list (swatch+value+phase, reorderable) | `MemberRow`/`TRow`×N | series-builder, series-plot | ⬜ |
| `SeriesPlate` | series builder/plot figure plate (header + waterfall + caption) | `PlateHeader`+`ToolBar`+`WaterfallChart` | series-builder, series-plot | ⬜ |
| `ReadingPanel` | derived phases-present summary | `ReadingRow`×N | series-plot | ⬜ |
| `BuilderRail` | builder sidebar (ordering/representation/display/traces sections) | `RailSection`×N + `Field`/`Slider`/`MemberList` | series-builder | ⬜ |

---

## Modals & overlays (hidden-by-default — highest late-discovery risk)

`ModalShell` primitive ✅ provides the scrim + shell. The **content** of each modal below is a gap.

### ★ `CustomIndexModal` — focus-plot (latest focus mockup), fully decomposed

Speculative-phase entry sheet. Triggered from the assignment cart's "custom index…" footer.

| Sub-component | What | Composes | Status |
|---|---|---|---|
| `CustomIndexModal` (shell) | scrim + sheet container | `ModalShell` | ⬜ |
| `ModalHead` | kicker + title + close button | `Kicker`+`IconButton` | ⬜ |
| `SymmetrySelector` | pick crystal symmetry (Pn3m/Im3m/Ia3d/Lamellar/Hexagonal…) | `SegmentedControl` | ⬜ |
| `LatticeParamControl` | dual control: slider + number field + unit | `Slider`+`Input`(number)+text | ⬜ |
| `CustomPreview` viz | live predicted-teeth vs observed-peaks render | `CombChart` reuse (Layer 1 gap) | ⬜ |
| `FitMetadata` | "landed on N peaks · a = … · snapped" line | text | ⬜ |
| `ModalFooter` | note + Cancel / Add-to-Assignment actions | `Button`×2 | ⬜ |

### Other overlays

| Component | What | Composes | Mockup | Status |
|---|---|---|---|---|
| ⚠️ **`ExportSheet`** | full-screen export dialog wrapping the clean render | `ModalShell`+`ModalHead`+`CleanFigure`+`ModalFooter` | series-plot | ⬜ |
| ⚠️ **`CullBar`** | floating batch-reject action bar (slides up on multi-select) | count + `Button`×N | sample-table | ⬜ |
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
| `MemberRow` series-plot `.member` variant (gradient coexistence swatch + inline colored phase names + lattice data) | ⏸ **split out** | 2026-06-02 | `MemberRow` (Batch 3) is the series-**builder** `.trow` only. The series-plot `.member` is a different row (gradient swatch for coexistence, per-phase inline names) — belongs with the later MemberList/SeriesPlate work. `Swatch` will gain a gradient/coexistence mode then (stub comment already in `Swatch.tsx`). |
| Type-scale snaps: `FolioHeader` title 31px→`text-display` (27); `PhaseBlock` name 22px→`text-headline` (19) | ✅ **snap to named role** | 2026-06-02 | The mockups use off-scale serif sizes (31/22px) with no named role; `text-[31px]` is guard-banned. Composites snap to the **nearest named role** rather than add a one-off scale step. If fidelity review objects, the fix is a new `--text-*` step in `styles.css` (the authoring layer), never an arbitrary in a composite. Folio count 26px = `text-headline-lg` is exact. |
| `NoticePill` is a NEW primitive, not a `Chip`/`Badge` variant (Batch 5) | ✅ **new ui/ primitive** | 2026-06-02 | The folio card-kick pill (`+N new match` accent-tint / dashed `Draft`) is a distinct 10px tinted/dashed look — `Chip` is a full-border pill, `Badge` a flat mono count, neither fits. The accent tint (`color-mix(…var(--color-accent) 14%…)`) + dashed border are appearance the placement-only `SeriesCard` cannot author, so they live in `ui/` (guard-exempt). No guard-allowlist edit needed — `ui/` was already exempt. `SeriesCard` stays placement-only by composing it. |
| `SeriesCard` figure region sits on PLATE, not paper-sunk (Batch 5 review fix) | ✅ **plate surface + paper-sunk hairline** | 2026-06-02 | First pass painted `bg-paper-sunk` across the whole figure pad; the mockup `.card-fig` keeps the pad on the card's plate surface with only a 1px paper-sunk bottom line — and `CardFigure`'s inner plots use `paperColor="var(--color-plate)"`, so a sunk pad created a plate-on-sunk seam. Fixed to `border-b border-paper-sunk` (the `--color-paper-sunk` Tailwind v4 `@theme` token generates `border-paper-sunk`, guard-clean). `cursor-pointer` also made conditional on `onClick`. |
| `src/print/export/` added to the design-guard `isExcluded` allowlist (Batch 4) | ✅ **new exempt renderer dir** | 2026-06-02 | `CleanFigure` is the export idiom — it DELIBERATELY sheds the Print look (Arial + literal hex `#c8841f`/`#4a4ba8`/`#111`, no `var(--color-*)` tokens), so it must author literal appearance like the other renderer dirs (`plot`/`comb`/`detector`). Edit is a single trailing-slash-bounded prefix in `scripts/check-design.mjs:isExcluded` — proven to still flag `print/components/` and not over-match `print/exporters-*`. Mirrors the legacy `lib/figure-export/**` export-palette allowlist. `CardFigure`/`CustomPreview` did NOT need exemption (literal-free, `phaseColor()` + tokens only). |
| `AssignmentCart` owns count→tag + positional-divider + empty (Batch 6) | ✅ **container owns cross-child state** | 2026-06-02 | The mockup `.phasecall` conditional `.pc-tag` ("Coexistence · N phases", only when >1 phase) and the `.pc-block + .pc-block` divider are state a lone `PhaseBlock` can't know — they belong to the container. `AssignmentCart` derives count from `Children.count`, renders the tag at count>1, and uses the Gallery-proven `Children.map` + positional `border-t border-hair` on a full-width `px-4` wrapper for **full-bleed** dividers (NOT `divide-y` on a padded wrapper — that insets the hairline). Wrappers carry `data-testid="cart-block" data-divider` so the load-bearing divider is asserted without class-string tests. `CandidateList` was already built (Batch 1, in `CandidateRow.tsx`) so this slice was just `AssignmentCart` + an optional content-only `series` prop on `PhaseBlock` (`.pcb-series`, guard-clean: `font-mono`/`text-caption`/`text-ink-*` only). |
| Focus panels are PRESENTATIONAL; state lifted to the consumer (Batch 7) | ✅ **panels render props, the page owns state** | 2026-06-02 | `TracePlate` (scale + `addPeakArmed`), `DetectorPanel` (current `src` + exposure-switcher element), and `CombsPanel` (`view` + `hoveredQ`) hold NO internal `useState` — every interaction datum is a prop with a consumer-owned handler. This is "container owns cross-child state" pushed one level up: the cross-PANEL state (the q-hover link that lights peak↔ring↔tooth across all three) can only live in the page, so panels must not capture it. Stories simulate the page with a small `useState` demo wrapper. The chrome is the `Card` primitive: `elevated` → the lifted hero plate (`TracePlate`), flat default → the two lower panels. |
| Focus lower panels share `PanelHeader`; `CombLegend` tooth-vocab deferred (Batch 7) | ✅ **shared `.panel-head` seam + deferred legend** | 2026-06-02 | `DetectorPanel`/`CombsPanel` share `PanelHeader` (faint `Kicker` label + right tools slot) — the `.panel-head`. The `.panel-h` label = `Kicker tone="faint"` exactly (uppercase 700 tracked ink-faint). `TracePlate` does NOT use it (the hero uses the big serif `PlateHeader`). The mockup `.combs-legend` is a comb-TOOTH vocabulary (predicted&observed / predicted-absent / leftover); the built `CombLegend` is the peak-GLYPH vocab (auto/manual/predicted-absent/excluded). `CombsPanel` uses the existing `CombLegend` for now — the tooth-specific legend is a DEFERRED `CombLegend` refinement (do not block a slice on it). `Thumbnail` gained an `xs` (≈30px) size for the dense Focus exposure strip via its existing numeric inline-style `SIZE_PX` idiom (guard-clean; `Thumbnail` lives in `components/`, NOT `ui/`). `exactOptionalPropertyTypes: true` forces conditional-spread of optional props forwarded to renderers whose matching prop is not `\| undefined` (`TracePlate.interaction`, `DetectorPanel.lutVariant`, `CombsPanel.hoveredQ`/`onHoverQ`). |
| Batch 7.1 — Focus-plates interaction-completeness pass | ✅ **wired the missing panel interactions** | 2026-06-02 | An interaction audit of the slice found three panel-level gaps (all fixed): **(a) TracePlate scroll-zoom** — `TracePlot` emitted `onXDomain` but the plate forwarded no controlled `xDomain`, so the zoom never rendered (TracePlot is fully controlled, no internal zoom state); added the `xDomain` pass-through so the page round-trips wheel→`onXDomain`→`xDomain`. **(b) DetectorPanel phase rings** — the panel rendered only the grayscale `DetectorImage`; the whole point (rings in phase colour) was missing. Now overlays `DetectorRings` via `buildRingPlacements(rings, calibration ?? null)` (null cal → centered fallback; `.det-box` made `relative` so the `absolute inset-0` overlay lands; `imageAspect` stays 1 on the square frame), and exposes `hoveredQ`/`onHoverQ` for the cross-panel "triple lights up" link. **(c) Gallery overflow** — the `PanelHeader` tools slot was width-unbounded so a long exposure strip widened the header; added `min-w-0` to the tools wrapper (+ `flex-shrink-0` on the label) so an `overflow-x-auto` `ThumbnailGallery` scrolls in place. **REMAINING ENGINE GAP (flagged, not yet built): `TracePlot` emits no OUTWARD hovered-q** — it has only internal `hoverId`. So hovering a *trace peak* can't drive the detector ring / comb tooth (the other two directions work). Closing the full triple-link needs a Layer-1 `TracePlot` addition: emit `onHoverQ` (peak hover → q) AND highlight on incoming `hoveredQ`. Deferred pending a scope call — it touches the trace-plot engine, not the panel layer. |

---

## Coverage summary

- **Layer 0 primitives:** 42/42 ✅ (foundation complete; +`Swatch` +`PhaseLabel`, Batch 3; +`NoticePill`, Batch 5)
- **Layer 1 renderers:** 9 ✅ — **layer COMPLETE** (Batch 4 closed `CardFigure`/`CustomPreview`/`CleanFigure`) · 1 ⛔ · 1 ⏸
- **Layer 2 tier-1 composites:** 17 ✅ — **layer COMPLETE** (Batch 3 closed `PhaseBlock`/`FolioHeader`/`AutoGroup`/`MemberRow`/`ReadingRow`)
- **Layer 3 tier-2 panels:** 8 ✅ (Batch 5: `SeriesCard`/`Gallery`; Batch 6: `AssignmentCart` + `CandidateList` confirmed-already-built; Batch 7: `TracePlate`/`DetectorPanel`/`CombsPanel` + shared `PanelHeader`) · ~10 ⬜
- **Modals/overlays:** 0 ✅ (`ModalShell` chrome ✅) · ~6 gap clusters ⬜
- **Layer 4 pages:** 0/5 ⬜

Layers 0, 1, and 2 are now all complete. The build frontier is the **Layer 3 tier-2 panels**, which are now fully unblocked (every renderer + tier-1 composite they compose exists):

- ~~`AssignmentCart`←`PhaseBlock` · `CandidateList`←`CandidateRow`~~ ✅ **DONE (Batch 6)** · `ReadingPanel`←`ReadingRow` · `MemberList`←`MemberRow`
- ~~`Gallery`←`SeriesCard`←(`CardFigure`+`PhaseStrip`)~~ ✅ **DONE (Batch 5)** · `SheetTable`←`SampleTableRow`
- ~~`TracePlate` · `DetectorPanel` · `CombsPanel` (+ shared `PanelHeader`)~~ ✅ **DONE (Batch 7)** · `SeriesPlate`←(`PlateHeader`+`ToolBar`+`WaterfallChart`)
- loupe panels (`BigFrame`, `LoupeSidePanel`), scoping (`ScopePlate`+`SampleRow`), `BuilderRail`
- modals/overlays: `CustomIndexModal`←`CustomPreview` · `ExportSheet`←`CleanFigure` · plus `CullBar`/`RailBack`/`Dock`/`ConflictModal`

**Batch 7 (2026-06-02) closed the Focus-plates slice** — the Focus hero trio (`TracePlate`/`DetectorPanel`/`CombsPanel`) + the shared `PanelHeader`, plus a `Thumbnail` `xs` size (refactor-on-contact). All three panels are presentational (scale/exposure/view/hovered-q lifted to the consumer). With Batch 6's assignment rail, **the entire Focus workspace's panel layer is now built — the Focus *page* (Layer 4) is assemblable.** Remaining next vertical slices, each fully unblocked:

1. **Series reading** — `MemberList`←`MemberRow` + `ReadingPanel`←`ReadingRow` + `SeriesPlate`←`WaterfallChart` (the series-builder reading surface).
2. **Contact sheet** — `SheetTable`←`SampleTableRow` + `CullBar` (the corpus grid).
3. **Loupe / scoping** — `BigFrame`+`LoupeSidePanel`, or `ScopePlate`+`SampleRow`.
