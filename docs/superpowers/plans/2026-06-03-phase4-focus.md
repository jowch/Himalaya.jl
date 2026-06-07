# Phase-4 Plan 2 — `PageFrame` + Focus workspace cutover

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish the `PageFrame` layout primitive (and refactor the Loupe onto it), then replace the legacy Focus workspace (`/sample/:sampleId`) with a greenfield `src/print/pages/FocusPage.tsx` assembled from the already-built print composites + the d3 `TracePlot` engine.

**Architecture:** The greenfield Focus page mounts body-only inside the carried `CorpusShell` `<Outlet>`, wrapped in `<PageFrame width="focus">`. It owns the cross-child state the presentational panels expect (active exposure, scale, xDomain, armed-add, the shared `hoveredQ` q-link, preview/speculative/custom-index modal state) and wires the CARRIED query hooks/mutators. Three legacy components with no greenfield equivalent (`StaleIndicesBanner`, `SpeculativeBuilder`, `FocusNotesMargin`) are ported to presentational `src/print/` versions; the d3 `TracePlot` gains the outward `onHoverQ` emit + incoming `hoveredQ` input that the q-link triple needs.

**Tech Stack:** React 18, TypeScript (`exactOptionalPropertyTypes: true`), TanStack Query, react-router-dom, Zustand, Vitest + RTL, Playwright (mocked e2e), boneyard-js skeletons, Tailwind closed-look design-guard (`npm run lint:design`), d3 (the greenfield plot engine).

**Provenance discipline (from the strategy spec):** every file is NEW (`src/print/**`), OLD (`src/pages/**`, `src/components/**` — deleted per-route), or CARRIED (`queries.ts`/`api.ts`/`state.ts`/`lib/**`/hooks — imported as-is). `src/print/pages/FocusPage.tsx` imports **only NEW + CARRIED, never OLD**. Verify with a grep before the final gate.

**Standing constraints (do not violate):**
- Commit ONLY specifically-named files. NEVER `git add -A`/`git add .`.
- NEVER stage `src/bones/registry.ts`, `src/bones/contact-sheet.bones.json`, or anything under `docs/superpowers/plans/`.
- Every commit's exact last line: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
- Typecheck: `npx tsc --noEmit -p tsconfig.build.json`. Work from `packages/HimalayaUI/frontend`.
- Design guard: `src/print/**` (incl. `pages/`) is NOT exempt — placement/layout/named-token/arbitrary-sizing classes only; no arbitrary `text-[…]`/`rounded-[…]`/raw-colour/side-stripe literals. `src/print/ui/**` and `src/print/plot/**` ARE exempt (appearance-authoring layers).
- Dispatch implementers SEQUENTIALLY (avoid `.git/index.lock` collisions).

---

## Reference: verified APIs (mapped from live source — confirm line numbers at execution; they drift)

### CARRIED hooks (`src/queries.ts`) + model (`src/api.ts`)
```ts
// queries.ts — all exposureId-scoped; queries disabled when id === undefined
useTrace(exposureId?)        -> UseQueryResult<Trace>          // Trace = { q:number[], I:number[], sigma:number[] }
usePeaks(exposureId?)        -> UseQueryResult<Peak[]>
useIndices(exposureId?)      -> UseQueryResult<IndexEntry[]>
useAssignment(exposureId?)   -> UseQueryResult<Assignment>     // durable 3-state cart
useAddPeak(exposureId)       -> .mutate(q: number)             // optimistic placeholder id < 0
useRemovePeak(exposureId)    -> .mutate(peakId: number)
useSetPeakExcluded(exposureId) -> .mutate({ peakId, excluded })
useAddAssignmentPhase(exposureId)    -> .mutate(indexId: number)
useRemoveAssignmentPhase(exposureId) -> .mutate(indexId: number)
useSetAssignmentState(exposureId)    -> .mutate(state: AssignmentState)
useCommitCustomIndex(exposureId)     -> .mutate(phase: string, basis: number)
useReanalyzeExposure(exposureId)     -> .mutate({})            // triggers hash refresh
useSpeculativeSnap(exposureId, phase?, anchorPeakId?, anchorRatio?) -> UseQueryResult<SpeculativeSnap[]>
useCreateSpeculative(exposureId)     -> useQueueMutation result
useDeleteIndex(exposureId)           -> .mutate(indexId: number)
useExposureHasPendingPeakOps(exposureId) -> boolean            // gates stale banner + speculative snap
// corpus/notes (cache-coherence: notes from useSamples(experimentId), NOT useCorpusSamples)
useCorpusSamples(); useExperiment(experimentId); useExposures(sampleId?); useSamples(experimentId);
useUpdateSample(experimentId, sampleId) -> .mutate({ ...partial Sample })
// derive: lib/assignment.ts
deriveActiveIndices(assignment, indices) -> IndexEntry[]       // active set from durable cart membership

// api.ts shapes
interface Peak { id; exposure_id; q; intensity:number|null; prominence:number|null; sharpness:number|null; source:"auto"|"manual"; excluded:boolean; }
interface IndexEntry { id; exposure_id; phase:string; basis:number; score:number|null; r_squared:number|null; lattice_d:number|null; ngc:number|null; status:"candidate"; kind:"auto"|"speculative"; inputs_hash:string|null; peaks:IndexPeakRef[]; predicted_q:number[]; bonnet?:{predicted_a:number;consistent:boolean}|null; }
type AssignmentState = "indexed" | "form_factor" | "null";
interface Assignment { exposure_id:number; state:AssignmentState; members:number[]; }  // index ids ascending
interface Trace { q:number[]; I:number[]; sigma:number[]; }
```

### Zustand (`src/state.ts`) — fields the legacy Focus surface uses
`activeSampleId`, `activeExposureId`, `setActiveExposure`, `activeExperimentId`, `hoveredIndexId`, `hoveredPeakId`, `previewIndexId`, `setPreviewIndex`, `notesDrawerOpen`, `closeNotesDrawer`, `openNavModal`. (The greenfield page keeps the route param as the source of truth for the sample; `hoveredQ` becomes **local page state**, not Zustand — see Task 3 rationale.)

### Greenfield composites already built (`src/print/components/`, presentational)
```ts
// TracePlate.tsx
interface TracePlateProps { kicker?; title; subtitle?; trace:TraceModel; scale:"log"|"lin"; onScaleChange:(n)=>void; onAutoFit?; addPeakArmed?:boolean; onToggleAddPeak?; interaction?:TracePlotInteraction|false; xDomain?:[number,number]|null; plotHeight?:number; className?; }
// AssignmentRail.tsx
interface AssignmentRailProps { assignment:ReactNode; assignmentCount?:ReactNode; assignmentNote?:ReactNode; candidates:ReactNode; candidatesNote?:ReactNode; className?; }
// AssignmentCart.tsx
interface AssignmentCartProps { children?:ReactNode; empty?:ReactNode; onCustomIndex?:()=>void; className?; }
// PhaseBlock.tsx
interface PhaseBlockProps { phase:string; score:number; meta:ReactNode; onRemove?:()=>void; series?:ReactNode; className?; }
// CandidateRow.tsx
interface CandidateRowProps { phase:string; score:number|null; why:ReactNode; bonnet?:boolean; selected?:boolean; onToggle?:()=>void; className?; }
interface CandidateListProps { children:ReactNode; className?; }
// CombsPanel.tsx
type CombView = "comb"|"resid";
interface CombsPanelProps { assigned:CombSeries[]; leftover?:number[]; view:CombView; onViewChange:(n)=>void; hoveredQ?:number; onHoverQ?:(q?)=>void; label?:ReactNode; className?; }
// DetectorPanel.tsx
interface DetectorPanelProps { src:string|null; label?; tools?:ReactNode; lutVariant?; rings?:RingInput[]; calibration?; beamCenter?; imageAspect?; hoveredQ?:number; onHoverQ?:(q?)=>void; hint?; className?; }
// CustomIndexModal.tsx (greenfield — fully presentational)
interface CustomIndexModalProps { open:boolean; onClose; onCancel; onAdd; symmetries:readonly string[]; symmetry:string; onSymmetryChange; paramName:string; paramValue:string; paramMin:number; paramMax:number; paramStep?:number; onParamChange; unit?; previewSeries:CombSeries; observed:number[]; onSelectObserved?; fit:{landed:number;total:number;snapped?:boolean}; className?; }
// RailSection.tsx { label; count?; note?; children?; className? }
// CombLegend.tsx { items?; className? }   PlateHeader, ToolBar, PanelHeader, Stepper exist.
```

### Greenfield engine (`src/print/plot/TracePlot.tsx`)
```ts
interface TracePlotProps { trace:TraceModel; xDomain?:[number,number]|null; xType?:"log"|"linear"; yType?; axes?; xLabel?; yLabel?; interaction?:TracePlotInteraction|false; overlay?:(ctx:PlotContext)=>ReactNode; height?; width?; paperColor?; className?; "data-testid"?; show?:{peaks?;labels?;band?}; highlightPeakIds?:ReadonlySet<number>; yHeadroom?; }
interface TracePlotInteraction { onXDomain:(d:[number,number]|null)=>void; onAddPeak?:(q)=>void; onClickPeak?:(peakId,altKey)=>void; onReset?:()=>void; hitTolerancePx?; }
interface TraceModel { trace:Trace; peaks:PlotPeak[]; phase:string|null; }
interface PlotPeak { id; q; intensity?; source:"auto"|"manual"; excluded?; predictedAbsent?; hot?; color?; label?; }
// GAP (Task 3): no outward onHoverQ emit; no incoming hoveredQ input. highlightPeakIds dims peaks NOT in the set (the losingPeakIds inverse).
```

### Legacy port sources (read these; do NOT import them into print/)
- `src/components/StaleIndicesBanner.tsx` — props `{exposureId?, debounceMs?=150}`; stale = count of `IndexEntry.inputs_hash !== Exposure.analysis_inputs_hash`; gated by `useExposureHasPendingPeakOps`; fires `useReanalyzeExposure.mutate({})`.
- `src/components/SpeculativeBuilder.tsx` — props `{exposureId, onClose}`; internal state `(phase, anchorPeakId, anchorRatio, included)`; uses `usePeaks`/`useSpeculativeSnap`/`useCreateSpeculative`; renders in `ModalShell`.
- `src/components/FocusNotesMargin.tsx` — props `{sample, onSaveNotes}`; q-ref regex `q\s*≈\s*\d+(?:\.\d+)?`; focus-gated textarea + draft state.
- `src/components/PlotCard.tsx` — the trace-hero data logic to port to page adapters: xDomain/yDomain state, `losingPeakIds` useMemo (~lines 250–263), trace→TraceModel mapping.
- `src/components/PhasePanel.tsx` — rail body logic (phase call card + candidate selection + custom-index commit + speculative disclosure).
- `src/components/FocusDetectorPanel.tsx` — phase-coloured ring-spec construction for `DetectorPanel.rings`.
- `src/components/CombPanel.tsx` — comb/resid series construction for `CombsPanel.assigned`/`leftover`.
- `src/components/FocusWorkspaceLayout.tsx` — the 3-col grid + drawer composition to mirror.

---

## File map
- **Create** `src/print/components/PageFrame.tsx` + `test/print-components/PageFrame.test.tsx` — layout primitive + width map.
- **Modify** `src/print/pages/LoupePage.tsx` — refactor onto `<PageFrame width="loupe">`.
- **Modify** `src/print/plot/TracePlot.tsx` (+ its test) — add `onHoverQ` emit + `hoveredQ` input.
- **Create** `src/print/components/{StaleBanner,SpeculativeDialog,NotesMargin}.tsx` (+ tests) — the three greenfield ports.
- **Create** `src/print/pages/focusAdapters.ts` + `test/print-pages/focusAdapters.test.ts` — pure page adapters (TraceModel, losingPeakIds complement, ring specs, comb series, custom-index SYMS+snap).
- **Create** `src/print/pages/FocusPage.tsx` + `test/print-pages/FocusPage.test.tsx` — the assembled page.
- **Modify** `src/components/AppRoutes.tsx` — repoint `/sample/:sampleId`.
- **Delete** legacy Focus-only components (grep-guarded): `FocusWorkspacePage`, `FocusWorkspaceLayout`, `FocusPlotHeader`, `FocusDetectorPanel`, `FocusNotesMargin`, `CombPanel`, `IndicesCard` (+ their tests). **Keep** shared: `PlotCard`/`TraceViewer`/`PlotSurface`/`PhasePanel`/`StaleIndicesBanner`/`SpeculativeBuilder`/`CustomIndexModal`(legacy)/`DetectorImage`/`DetectorRingOverlay` **until** they are provably unimported (some are Index-shared — grep first; defer any still-referenced deletion).
- **Capture/commit** `src/bones/focus.bones.json` (only this file; never `registry.ts`).

---

## Task 1: `PageFrame` layout primitive + width map

**Files:** Create `src/print/components/PageFrame.tsx`, `test/print-components/PageFrame.test.tsx`

- [ ] **Step 1: Write the failing test** — `test/print-components/PageFrame.test.tsx`
```tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { PageFrame, PAGE_WIDTHS } from "../../src/print/components/PageFrame";

describe("PageFrame", () => {
  it("renders children inside a centered frame", () => {
    render(<PageFrame width="loupe"><p>body</p></PageFrame>);
    const frame = screen.getByTestId("page-frame");
    expect(frame).toHaveTextContent("body");
    expect(frame.className).toContain("mx-auto");
  });
  it("applies the named max-width class for each surface", () => {
    for (const [name, cls] of Object.entries(PAGE_WIDTHS)) {
      const { unmount } = render(<PageFrame width={name as keyof typeof PAGE_WIDTHS}>x</PageFrame>);
      expect(screen.getByTestId("page-frame").className).toContain(cls);
      unmount();
    }
  });
  it("forwards a placement className", () => {
    render(<PageFrame width="sheet" className="py-7">x</PageFrame>);
    expect(screen.getByTestId("page-frame").className).toContain("py-7");
  });
});
```

- [ ] **Step 2: Run → fail** — `npm test -- print-components/PageFrame` → cannot resolve module.

- [ ] **Step 3: Implement** — `src/print/components/PageFrame.tsx`
```tsx
import type { ReactNode } from "react";

/** Single source of truth for per-surface content widths (mockup values).
 *  Changing a surface's width is a one-line edit here. The topbar app-frame
 *  reads `PAGE_WIDTHS.sheet` so the chrome stays coherent. */
export const PAGE_WIDTHS = {
  loupe: "max-w-[1100px]",
  sheet: "max-w-[1240px]",
  focus: "max-w-[1160px]",
  folio: "max-w-[1380px]",
  scoping: "max-w-[760px]",
  builder: "max-w-[1180px]",
} as const;

export type PageWidth = keyof typeof PAGE_WIDTHS;

function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

/** Centers + caps a page body at its surface width. Greenfield pages wrap their
 *  body in this instead of a hand-rolled `mx-auto max-w-[…]`. PLACEMENT-ONLY
 *  className (padding/margin). */
export function PageFrame({
  width,
  className,
  children,
}: {
  width: PageWidth;
  className?: string;
  children: ReactNode;
}): JSX.Element {
  return (
    <div data-testid="page-frame" className={cx("mx-auto w-full", PAGE_WIDTHS[width], className)}>
      {children}
    </div>
  );
}
```

- [ ] **Step 4: Run → pass** — `npm test -- print-components/PageFrame`.
- [ ] **Step 5: Gate** — `npx tsc --noEmit -p tsconfig.build.json`; `npm run lint:design` (arbitrary `max-w-[…]` is sizing = allowed) → exit 0.
- [ ] **Step 6: Commit**
```bash
git add src/print/components/PageFrame.tsx test/print-components/PageFrame.test.tsx
git commit -m "$(printf 'Phase-4: PageFrame layout primitive + per-surface width map\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 2: Refactor the Loupe onto `PageFrame` + share the width with the topbar

**Files:** Modify `src/print/pages/LoupePage.tsx`, `src/components/CorpusTopbar.tsx`

- [ ] **Step 1:** In `LoupePage.tsx`, replace the two body wrappers `<div data-testid="loupe-page" className="mx-auto max-w-[1100px] px-8 py-7">` (the not-found branch and the main branch) with `<PageFrame width="loupe" className="px-8 py-7"><div data-testid="loupe-page">…</div></PageFrame>` — keep `data-testid="loupe-page"` on an inner div so existing tests/e2e still find it. Add `import { PageFrame } from "../components/PageFrame";`. Update the boneyard `LOUPE_FIXTURE` wrapper the same way if it carried the literal.

- [ ] **Step 2:** In `CorpusTopbar.tsx`, replace the literal `max-w-[1240px]` on the `corpus-topbar-inner` div with the shared source: `import { PAGE_WIDTHS } from "../../print/components/PageFrame";` and use `` className={`mx-auto flex h-full w-full ${PAGE_WIDTHS.sheet} items-center gap-4 px-6`} ``. (Cross-tree import: `src/components` → `src/print` is allowed; the boundary rule forbids the reverse — print importing legacy.)

- [ ] **Step 3: Run** — `npm test -- print-pages/LoupePage CorpusTopbar` → green (loupe-page testid still present; topbar unchanged in behaviour).
- [ ] **Step 4: Gate** — `npx tsc --noEmit -p tsconfig.build.json`; `npm run lint:design`; `npm run e2e -- corpus-culling` → all green.
- [ ] **Step 5: Commit**
```bash
git add src/print/pages/LoupePage.tsx src/components/CorpusTopbar.tsx
git commit -m "$(printf 'Phase-4: adopt PageFrame in the Loupe + source the topbar width from it\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 3: TracePlot q-link — outward `onHoverQ` emit + incoming `hoveredQ` input

**Files:** Modify `src/print/plot/TracePlot.tsx` + its test (`test/print-plot/TracePlot.test.tsx` — confirm exact path).

**Rationale:** the q-link triple (trace ⇄ comb ⇄ detector) shares one hovered-q. `CombsPanel`/`DetectorPanel` already expose `hoveredQ`/`onHoverQ`; `TracePlot` has only *internal* hover. The Focus page owns the shared `hoveredQ` state; TracePlot must both emit (on peak hover) and accept (to light a peak hovered elsewhere). `src/print/plot/**` is design-guard-exempt.

- [ ] **Step 1: Write failing tests** — assert: (a) hovering within `hitTolerancePx` of a peak fires `onHoverQ(peak.q)`; pointer-leave fires `onHoverQ(undefined)`. (b) Passing `hoveredQ` near a peak's q marks that peak `hot` (same visual as internal hover) — assert via the peak glyph's hot data-attribute/testid the engine already uses for internal hover (read `PlotPeaks`/`TracePlot` for the exact hot marker, mirror it). Use the existing TracePlot test harness/fixtures.

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement** — add to `TracePlotProps`: `onHoverQ?: (q: number | undefined) => void;` and `hoveredQ?: number;`. In the pointer-move handler that sets internal `hoverId`, when a peak is within tolerance call `onHoverQ?.(peak.q)`; in pointer-leave call `onHoverQ?.(undefined)`. Derive the effective hot peak as `internalHoverId ?? peakNearest(hoveredQ)` so an external `hoveredQ` lights the matching glyph (tolerance in q-space: nearest peak with `|peak.q - hoveredQ|` minimal within the same px tolerance projected to q). Keep internal hover authoritative when both present. Do not change existing props/behaviour.

- [ ] **Step 4: Run → pass.**
- [ ] **Step 5: Gate** — `npx tsc …`; `npm run lint:design`; run the full plot suite `npm test -- print-plot` (ensure no regression to scroll-zoom/add/click).
- [ ] **Step 6: Commit**
```bash
git add src/print/plot/TracePlot.tsx test/print-plot/TracePlot.test.tsx
git commit -m "$(printf 'Phase-4: TracePlot q-link — emit onHoverQ + accept hoveredQ for the focus triple\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 4: Port `StaleBanner` (greenfield, presentational)

**Files:** Create `src/print/components/StaleBanner.tsx` + `test/print-components/StaleBanner.test.tsx`. Port source: `src/components/StaleIndicesBanner.tsx` (read it).

**Design:** presentational — the page lifts the staleness computation + gating. Props:
```ts
interface StaleBannerProps {
  staleCount: number;        // 0 → render nothing
  pending?: boolean;         // useExposureHasPendingPeakOps — disables the button
  onReanalyze: () => void;
  className?: string;
}
```
- [ ] **Step 1:** Test — `staleCount===0` renders nothing; `staleCount>0` renders the count + a "Re-analyze" button (assert by role/text, not class); clicking fires `onReanalyze`; `pending` disables the button. Mirror the legacy copy.
- [ ] **Step 2:** Run → fail.
- [ ] **Step 3:** Implement from `src/print/ui` primitives (`NoticePill`/`Button`/`HintText` — read the legacy banner's copy + the greenfield ui/ barrel for the right primitive). Placement-only; appearance in `ui/`.
- [ ] **Step 4:** Run → pass. **Step 5:** Gate (tsc + lint:design). **Step 6:** Commit `Phase-4 Focus: port StaleBanner (presentational)`.

---

## Task 5: Port `SpeculativeDialog` (greenfield, presentational)

**Files:** Create `src/print/components/SpeculativeDialog.tsx` + test. Port source: `src/components/SpeculativeBuilder.tsx`.

**Design:** presentational modal — state lifts to the page (phase, anchorPeakId, anchorRatio, included). Props (derive exact shape from the legacy component's rendered controls):
```ts
interface SpeculativeDialogProps {
  open: boolean;
  phases: readonly string[]; phase: string; onPhaseChange: (p: string) => void;
  peaks: { id: number; q: number }[]; anchorPeakId: number | undefined; onAnchorChange: (id: number) => void;
  anchorRatio: number; onAnchorRatioChange: (r: number) => void;
  snap: SpeculativeSnap[];                 // preview rows from useSpeculativeSnap
  pending?: boolean;
  onCreate: () => void; onClose: () => void;
  className?: string;
}
```
- [ ] **Step 1:** Test the presentational contract (open gates render; changing phase/anchor fires callbacks; Create fires `onCreate`; `pending` disables Create). Use `ModalShell` from `src/print/ui`.
- [ ] **Step 2–6:** fail → implement (compose `ModalShell`/`Field`/`SegmentedControl`/`Button` from `src/print/ui`) → pass → gate → commit `Phase-4 Focus: port SpeculativeDialog (presentational)`.

---

## Task 6: Port `NotesMargin` (greenfield)

**Files:** Create `src/print/components/NotesMargin.tsx` + test. Port source: `src/components/FocusNotesMargin.tsx`.

**Design:** draft + focus-gate state is purely local UI → keep it internal. Props:
```ts
interface NotesMarginProps {
  notes: string | null;
  onSaveNotes: (notes: string) => void;   // fire on blur
  onHoverQ?: (q: number | undefined) => void;   // q-ref hover → q-link (optional)
  className?: string;
}
```
- [ ] **Step 1:** Test — renders existing notes; q-ref tokens (regex `q\s*≈\s*\d+(?:\.\d+)?`) render as highlighted spans; editing then blurring fires `onSaveNotes(draft)`; the textarea is focus-gated (mirror legacy). Assert by role/text.
- [ ] **Step 2–6:** fail → implement (compose `src/print/ui` text primitives; keep the q-ref highlight + focus-gate logic from the legacy source) → pass → gate → commit `Phase-4 Focus: port NotesMargin`.

---

## Task 7: Focus page adapters (pure)

**Files:** Create `src/print/pages/focusAdapters.ts` + `test/print-pages/focusAdapters.test.ts`.

Pure functions (no React). Derive exact logic from the cited legacy sources:
```ts
toTraceModel(trace, peaks, activeIndices, phase): TraceModel          // Peak[] → PlotPeak[] (source/excluded/predictedAbsent/color); phase from assignment
complementPeakIds(allPeakIds: number[], losing: Set<number>): Set<number>   // highlightPeakIds = all − losing  (the dim inverse; PlotCard.tsx ~250-263)
losingPeakIds(hoveredCandidate, activeIndices, peaks): Set<number>    // peaks an add-of-hoveredCandidate would orphan
toDetectorRings(activeIndices, calibration): RingInput[]              // phase-coloured Debye–Scherrer rings (FocusDetectorPanel)
toCombSeries(activeIndices, peaks): { assigned: CombSeries[]; leftover: number[] }   // CombPanel construction
// custom-index: the SYMS table + snap math lifted out of the legacy CustomIndexModal
CUSTOM_SYMS: readonly { name: string; paramName: string; unit: string; min: number; max: number; step?: number }[];
customIndexPreview(symmetry, paramValue, observed): { previewSeries: CombSeries; fit: { landed; total; snapped? } };
```
- [ ] **Step 1:** Unit-test each with small fixtures (especially `complementPeakIds`, `losingPeakIds`, and `customIndexPreview` fit counting — these are load-bearing). **Step 2–6:** fail → implement (port the math from `PlotCard.tsx`/`FocusDetectorPanel.tsx`/`CombPanel.tsx`/legacy `CustomIndexModal.tsx`) → pass → gate → commit `Phase-4 Focus: pure page adapters`.

> SAXS-physics note: `losingPeakIds`/`toCombSeries`/ring construction touch indexing/phase math. After this task, request a `saxs-physics-reviewer` pass before the assembly task.

---

## Task 8: Assemble `FocusPage`

**Files:** Create `src/print/pages/FocusPage.tsx` + `test/print-pages/FocusPage.test.tsx`.

Mirror the legacy 3-column grid (`FocusWorkspaceLayout.tsx`): `xl:grid-cols-[1fr_348px_250px]` — work area (TracePlate over [DetectorPanel | CombsPanel]) · AssignmentRail · NotesMargin (xl+; drawer below xl via `ModalShell variant="drawer"`). Wrap the body in `<PageFrame width="focus">`.

**The page owns** (lift state here): `scale`, `addArmed`, `xDomain`, the shared `hoveredQ` (local `useState`, fed to TracePlate/CombsPanel/DetectorPanel/NotesMargin), `combView`, `previewIndexId`, speculative-dialog state, custom-index-modal state. Sample id = route param (`useParams`); active exposure from Zustand `activeExposureId` (carry `useSyncActiveSampleFromRoute` + the auto-pick). Notes from `useSamples(experimentId)` (cache-coherence — NOT `useCorpusSamples`).

**Wire:** `useTrace`/`usePeaks`/`useIndices`/`useAssignment` + `deriveActiveIndices`; mutators `useAddPeak`/`useRemovePeak`/`useSetPeakExcluded`/`useAddAssignmentPhase`/`useRemoveAssignmentPhase`/`useSetAssignmentState`/`useCommitCustomIndex`/`useReanalyzeExposure`/`useCreateSpeculative`/`useDeleteIndex`. Stale banner: compute `staleCount` from indices vs `exposure.analysis_inputs_hash`, pass to `StaleBanner` with `useExposureHasPendingPeakOps`. `losingPeakIds` → `complementPeakIds` → `TracePlate`/`TracePlot.highlightPeakIds`.

- [ ] **Step 1: Write failing tests** (mock `../../src/queries` + `../../src/state` + `boneyard-js/react` like the Loupe test): page renders TracePlate + AssignmentRail + DetectorPanel + CombsPanel; selecting a candidate calls `useAddAssignmentPhase.mutate(indexId)`; arming + clicking the plot adds a peak; the q-link round-trips (hovering a comb row sets `hoveredQ`, lighting the trace peak); `+ custom index…` opens the modal; not-found + no-exposure states. Assert by testid/role/behaviour — never class strings.
- [ ] **Step 2: Run → fail.**
- [ ] **Step 3: Implement** `FocusPage.tsx` (NEW + CARRIED imports only). `Skeleton name="focus"` with a greenfield fixture + `fallback`.
- [ ] **Step 4: Run → pass.**
- [ ] **Step 5: Gate** — tsc; lint:design; `npm test -- print-pages/FocusPage`.
- [ ] **Step 6: Commit** `Phase-4 Focus: assemble greenfield FocusPage`.

---

## Task 9: Repoint the route + delete the OLD

**Files:** Modify `src/components/AppRoutes.tsx`; delete legacy Focus-only files.

- [ ] **Step 1:** Repoint — `grep -n "FocusWorkspacePage" src/components/AppRoutes.tsx`; change the import to `import { FocusPage } from "../print/pages/FocusPage";` and the route element to `<FocusPage />` (keep the `/sample/:sampleId` path). tsc green.
- [ ] **Step 2: Grep-guard before each delete** — `grep -rn "FocusWorkspacePage\|FocusWorkspaceLayout\|FocusPlotHeader\|FocusDetectorPanel\|FocusNotesMargin\|CombPanel\|IndicesCard" src test e2e`. Delete only files with **no remaining importers outside themselves/their tests**. `PhasePanel`/`PlotCard`/`TraceViewer`/`PlotSurface`/`StaleIndicesBanner`/`SpeculativeBuilder`/legacy `CustomIndexModal`/`DetectorImage`/`DetectorRingOverlay` are **shared with the Index/other surfaces** — confirm unimported before deleting; if still referenced, DEFER (note it) — do not break other surfaces.
- [ ] **Step 3:** `git rm` the confirmed-orphaned Focus-only files + their tests.
- [ ] **Step 4: Verify whole** — `npx tsc …` (no dangling imports); `npm test -- focus Focus` (only greenfield focus tests run + pass).
- [ ] **Step 5: Commit** `Phase-4 Focus: route /sample/:sampleId -> greenfield + delete legacy Focus`.

---

## Task 10: Boneyard recapture

**Files:** Capture/commit `src/bones/focus.bones.json` ONLY (never `registry.ts`).

- [ ] **Step 1:** `npm run dev`, open `/sample/<a real sample id>` against the live/mocked backend, trigger the cold-loading state so the plugin writes `src/bones/focus.bones.json`. **Step 2:** Verify the geometry matches the 3-col layout. **Step 3:** `git add src/bones/focus.bones.json` (confirm `registry.ts` is NOT staged — if wiring requires a `registry.ts` change, STOP and surface to the human) → commit `Phase-4 Focus: capture greenfield focus skeleton`.

---

## Task 11: Final gate + visual fidelity

- [ ] **Step 1: Provenance grep** — `grep -rn 'from "\.\./\.\./components\|from "\.\./\.\./pages\|/src/components/\|/src/pages/' src/print/pages/FocusPage.tsx src/print/pages/focusAdapters.ts` → no matches (NEW + CARRIED only).
- [ ] **Step 2: Full local gate** — `npm run lint:design`; `npx tsc --noEmit -p tsconfig.build.json`; `npm test -- print-pages print-components print-plot`; `npm run e2e` (mocked, incl. any focus spec); `npm run build` → all green.
- [ ] **Step 3: Manual visual fidelity** — `npm run dev`, open `/sample/:id`, compare against `docs/redesign-mockups/focus-workspace.html` + `docs/redesign-mockups/2026-05-29-focus-plot.html`: 3-col layout, trace hero with armed-add + auto-fit + log/lin, the q-link triple lighting across trace/comb/detector, assignment cart + candidate rows + custom-index modal, stale banner, notes margin/drawer. Fix placement-only issues in the page, appearance in `ui/` primitives.
- [ ] **Step 4: saxs-physics-reviewer pass** on the adapters (Task 7) + any comb/ring/index math. Then done — `/sample/:sampleId` runs greenfield under the carried shell; PageFrame is established for Samples (plan 3).

---

## Self-review checklist (run before declaring done)
- **Spec coverage:** PageFrame (T1) + Loupe refactor + topbar width source (T2) · q-link engine (T3) · 3 ports (T4–6) · adapters incl. losingPeakIds/custom-index (T7) · assembly (T8) · repoint+delete (T9) · boneyard (T10) · gate+fidelity+physics review (T11). ✅
- **Type consistency:** hook/mutator names + `.mutate` shapes match the verified `queries.ts`; composite prop names match the verified `src/print/components/` interfaces; `PlotPeak`/`TraceModel`/`IndexEntry`/`Assignment` field names match `api.ts`/`TracePlot.tsx`. ✅
- **Provenance:** FocusPage imports NEW + CARRIED only; T11.1 enforces by grep. Shared-component deletes are grep-guarded + deferred if still referenced. ✅
- **Verify-at-execution:** legacy port line numbers + the exact TracePlot hot-marker testid drift — implementers confirm against live source (the reference section says so). ✅
