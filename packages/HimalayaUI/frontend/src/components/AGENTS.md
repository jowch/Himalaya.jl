# frontend/src/components — UI components

React components for the three-card Index workspace, Inspect page, Compare page, modals, and shared primitives.

## Where things live

- **Shell + chrome**: `AppShell`, `AppHeader`, `TabRocker`, `TitleButton`, `UtilityCluster`, `WorkspaceGrid`, `BandResizeDivider`, `AnnotationToggles`.
- **Index workspace cards**: `ChatCard`, `PlotCard` + `TraceViewer`, `IndicesCard`, `PhasePanel`, `StaleIndicesBanner`, `SpeculativeBuilder`, `MillerPlot`, `Pn3mIcon`.
- **Inspect**: `DetectorImage`, `DetectorImageCard`, `ThumbnailGallery`, `SampleMetadataCard`.
- **Mentions**: `MentionChip`, `MentionCompose`, `MentionPicker`.
- **Compare**: `MultiTracePlot`, `MemberTraceLayer`, `MemberMetaRow`, `MemberMetaGutter`, `ComparisonSidebar`, `ComparisonPickerBody`, `ComparisonPickerPanel`, `SamplePickerRow`, `GroupingModeToggle`, `WarmAddMenu`, `CompareTitleStrip`, `CompareStatusSurface`, `CompareToolbar`, `SavePill`, `InlineEditableText`, `RowActionZone`.
- **Onboarding / modals / errors**: `OnboardingFlow`, `NavModal`, `ConflictModal`, `StaleUrlPage`, `ResolvingFallback`, `InfrastructureBanner`, `FigureExportControls`.
- **Primitives**: `components/ui/`.

## Observable Plot inside React

The plot element has a runtime `.scale(name).invert(px)` method that isn't in DOM types; cast with `(el as unknown as { scale: ... })`. Used by `TraceViewer` to translate click pixel coords to q values.

## `QNumInput` (exported from `PlotCard.tsx`)

Focus-gated controlled input: external `value` prop changes are synced to draft state only when the input is **not** focused, preventing wheel-zoom events from interrupting mid-edit. Follow this pattern for any numeric input that can be updated by external events.

## `StaleIndicesBanner` (mounted in `PhasePanel`)

Renders when *any* index's `inputs_hash` differs from its exposure's current `analysis_inputs_hash` (hash-derived, not a `status` enum — that was removed in Plan 7 R3). The Re-analyze button posts to `/api/exposures/:id/analyze`, which recomputes hashes; matching hashes hide the banner. New routes that change the effective peak set surface the banner automatically because hashes drift; no extra UI wiring needed.

**Plan 8 update**: also gated on `useExposureHasPendingPeakOps` (returns null while a peak op is in flight) to mask cross-entity refetch races during queue mutations. Debounce reduced from 2000ms → 150ms because synchronous reanalyze in the curation handler closes the stale window deterministically. See `docs/event-log.md` §2 + §3a.

## `DetectorImage` / `ImageBitmap.close()`

- **`ImageBitmap.close()` neuters width/height to 0** per the Web spec. Capture dims *before* closing: `const { width, height } = bitmap; bitmap.close();`. Regression test in `test/DetectorImage.test.tsx` uses getter-based mocks that simulate the neutering — keep it green when touching `DetectorImage.tsx`.
- **Auto-rotates to landscape** via a `ResizeObserver` on the wrapper div: when `containerAspect > imageAspect * 1.25`, rotate the canvas 90° and JS-set `maxWidth`/`maxHeight` to swap the layout box (CSS-only doesn't cut it because `transform: rotate` doesn't change a canvas's bounding box). Re-evaluates inside `renderImage` after each new image so swapping exposures with different aspects re-checks. JSDOM lacks `ResizeObserver` — the stub in `test/setup.ts` keeps unit tests honest.

## `TraceViewer` auto-fit (floor-only)

`PlotCard::computeFit` sets `yDomain = [max(p01·0.5, fullMax/1e5), fullMax·1.2]`:

- **Bottom**: 1st percentile of *positive* in-window intensities scaled down — suppresses dead-pixel zeros while keeping the low-signal tail visible. Clamped at `fullMax/1e5` so a single near-zero pixel can't blow the y-range past five log decades.
- **Top**: full trace max (so peaks-vs-beam relative scale stays visible without resetting).
- When peaks exist, x is also tightened to `[firstPeak·0.7, lastPeak·1.3]`.
- Auto-fires on `activeExposureId` change. Double-click → `onReset` clears both axes.

## Focus trapping in modals

`hooks/useFocusTrap.ts` exports `useFocusTrap(containerRef, active)`. Call inside any modal/overlay that should keep Tab focus within its bounds. It intercepts Tab/Shift+Tab to cycle among focusable children and restores the previously-focused element on cleanup. `NavModal` and `OnboardingFlow` already use it.

## Anti-patterns

- Don't reach for `transform: rotate` and expect bounding boxes to follow — `DetectorImage`'s ResizeObserver + JS-set max-w/h pattern is load-bearing.
- Don't introspect tooltip text or Tailwind classes in tests — use `data-*` attributes.
