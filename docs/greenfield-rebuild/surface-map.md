# Surface map (mockup → new page → components)

Source mockups: `docs/redesign-mockups/`
Target tree: `src/print/` (pages under `src/print/pages/`, components under `src/print/components/`)

---

| Mockup | New page (`src/print/pages/`) | Key components (`src/print/components/`) | Renderer(s) |
|---|---|---|---|
| `sample-table.html` (contact-sheet view) | `ContactSheetPage` | `TopBar`, `StageTabs`, `FacetChip`, `DataTable`, `SampleRow`, `ScreenedMark`, `ThumbnailStrip`, `Thumbnail`, `RejectOverlay`, `TagList`, `TagPill`, `PhaseChip`, `ProgressBar`, `KbLegend`, `KbKey`, `FloatingCullBar` | — |
| `sample-table.html` (loupe section — same file, own route) | `LoupePage` | `TopBar`, `StageTabs`, `SegmentedControl` (contact-sheet/loupe toggle), `Card` (loupe-plate), `ThumbnailStrip`, `Thumbnail`, `RejectOverlay`, `MetaList`, `SignalBars`, `VerdictCard`, `RepresentativeBox`, `TagList`, `TagPill`, `KbLegend`, `KbKey`, `Kicker`, `IconButton` | `DetectorImage` (big-frame) |
| `focus-workspace.html` | `FocusWorkspacePage` | `TopBar`, `StageTabs`, `Stepper`, `Card` (trace plate, detector panel, reflections panel), `Kicker`, `SegmentedControl` (log/linear), `Button` (Auto-fit, + Peak), `ExposureSwitcher`, `PhasecallBlock`, `CandidateRow`, `ScoreBar`, `NotesMargin`, `NoteEntry`, `HintText`, `Dot`, `PhaseChip`, `KbKey` | `TraceViewer` (trace SVG hero), `DetectorImage`, `ReflectionsTable` |
| `2026-05-29-focus-plot.html` | `FocusWorkspacePage` (replaces/augments the above — same route `/sample/:id`) | `TopBar`, `StageTabs`, `Stepper`, `Card` (trace plate, detector panel, combs panel), `Kicker`, `SegmentedControl` (log/linear, comb/resid), `Button` (Auto-fit, + Peak), `ExposureSwitcher`, `PeakGlyph`, `CandidateRow`, `BonnetBadge`, `PhasecallBlock` (assignment cart variant), `ScoreBar`, `CombPanel`, `CombLegend`, `CustomIndexSheet`, `ModalShell`, `HintText` | `TraceViewer` (new comb-aware), `DetectorImage`, `CombsRenderer` |
| `series-folio.html` | `SeriesFolioPage` | `TopBar`, `StageTabs`, `Button` (+ New series), `Kicker`, `SearchInput`, `FilterChip`, `SegmentedControl` (sort), `FolioCard`, `MiniWaterfall`, `PhaseStrip`, `TagPill`, `EmptyState` | — |
| `series-scoping.html` | `SeriesScopingPage` | `TopBar`, `StageTabs`, `Button` (ghost Discard), `Kicker`, `Card` (scope-plate), `AutogroupBanner`, `OrderingField`, `ScopingSampleRow`, `GripHandle`, `Sparkline`, `PhaseStrip`, `SegmentedControl` (implied variable picker), `HintText`, `Dot`, `KbKey` | — |
| `series-builder.html` + `2026-05-29-series-plot.html` | `SeriesBuilderPage` | `TopBar`, `StageTabs`, `Button` (Duplicate, Export), `Card` (plate), `Kicker`, `TagPill` (fig-tag metadata pills), `MiniWaterfall`/`WaterfallPlot`, `SegmentedControl` (waterfall/heatmap, log/linear), `AutogroupBanner`, `OrderingField`, `Slider`, `ToggleSwitch`, `TraceListRow`, `GripHandle`, `PhaseChip`, `Button` (+ Add sample, Copy as PNG), `HintText`, `RegionBracket`, `FloatingDock`, `RailbackTab`, `PhaseStrip`, `PhaseReadingPanel`, `MemberList`, `ExportButton` | `MultiTracePlot` (waterfall + heatmap renderers), `PhaseStrip` (companion strip) |

---

## Page notes

### `ContactSheetPage` (`/`)
The main table/grid surface. Two sub-views (contact-sheet and loupe) share the route and toggle via a `SegmentedControl` in the topbar. The loupe is currently shown in the same file but should be a separate routed sub-view (or distinct route `/sample/:id/loupe`).

### `LoupePage` (`/sample/:id/loupe` or inline toggle on ContactSheetPage)
The full-bleed single-exposure inspection surface. Detector image dominates the left column; metadata, verdict, rep-picker, and tags are in the right aside. Navigation is purely keyboard (← →, X, R, Esc).

### `FocusWorkspacePage` (`/sample/:id`)
The indexing workspace. The 2026-05-29-focus-plot mockup is the redesigned version of the focus-workspace (same page/route). It replaces the reflections table with the CombPanel (comb teeth on a shared q ruler), replaces the "phase call" card with an "assignment cart" (empty start, populated by adding candidates), and adds the BonnetBadge flag on candidates. The old `ReflectionsTable` is replaced by `CombsRenderer`; `TraceViewer` gains the new `PeakGlyph` (triangle-down / diamond / caret) encoding.

### `SeriesFolioPage` (`/series`)
The masonry wall of saved series. Each card is a `FolioCard` containing a `MiniWaterfall` (the frozen figure). The wall height is variable per card (masonry column layout). Draft cards use dashed border + lower opacity.

### `SeriesScopingPage` (`/series/new` or `/series/new?from=<sampleIds>`)
Entered from multi-select on the contact sheet. Centered worksheet plate. Reviews Himalaya's auto-grouping proposal, confirms the ordering variable, resolves any flagged parsed values, and gates the "Confirm & build" action until all values are confirmed.

### `SeriesBuilderPage` (`/series/:id`)
The series figure and editing rail. The 2026-05-29-series-plot mockup is the redesigned version of the series-builder (same route). Key differences: the on-screen view is a thinking instrument, export is decoupled (clean scientific render, not a screen-grab), and the phase-strip companion runs vertically alongside the waterfall SVG. Rail collapses to full-bleed with a floating dock for the offset slider. Both waterfall and heatmap representations share the same `MultiTracePlot` renderer via a `repr` prop.

---

## Shared layout primitives (not page-specific)

| Component | Used by |
|---|---|
| `TopBar` | All pages |
| `StageTabs` | All pages |
| `Kicker` | ContactSheetPage, FocusWorkspacePage, SeriesFolioPage, SeriesScopingPage, SeriesBuilderPage, LoupePage |
| `SegmentedControl` | ContactSheetPage (sheet/loupe), FocusWorkspacePage (log/lin, comb/resid), SeriesFolioPage (sort), SeriesBuilderPage (waterfall/heatmap, log/lin) |
| `Card` | FocusWorkspacePage (panels), SeriesScopingPage (scope-plate), LoupePage (loupe-plate) |
| `Button` | All pages |
| `HintText` | FocusWorkspacePage, SeriesScopingPage, SeriesBuilderPage |
| `PhaseChip` | ContactSheetPage, FocusWorkspacePage, SeriesBuilderPage, LoupePage |
| `PhaseStrip` | SeriesFolioPage, SeriesScopingPage, SeriesBuilderPage |
| `ScoreBar` | FocusWorkspacePage |
| `Dot` | All pages (topbar tab active dot) |
| `GripHandle` | SeriesScopingPage, SeriesBuilderPage |
| `PeakGlyph` | FocusWorkspacePage (focus-plot version), SeriesBuilderPage (track anchors) |
| `EmptyState` | SeriesFolioPage |
