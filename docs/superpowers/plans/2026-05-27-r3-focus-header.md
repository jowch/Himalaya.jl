# R3 — Focus header (serif sample title + provenance) Implementation Plan

> **For agentic workers:** TDD, one commit per step group.

**Goal:** Give the focus-workspace trace plate a proper header — terracotta "Integration" kicker, Newsreader serif sample name, mono `smp · beamtime · representative exposure` subline — seeded from the route's sample, replacing the leftover experiment-picker placeholder ("pick an experiment / click to change").

**Architecture:** Add an opt-in `headerSlot?: JSX.Element` prop to `PlotCard`. When supplied, `TitleStrip` renders that node on the left in place of the experiment-picker button; the right-side q-controls/export cluster is unchanged. A new `FocusPlotHeader` component renders the kicker/serif/sub stack from sample + experiment + exposure data. `FocusWorkspaceLayout` (which already resolves `corpusSample` and `experimentId`) builds the `FocusPlotHeader` and passes it as `headerSlot`. Default (prop-less) PlotCard keeps the picker — opt-in, so no other usage changes.

**Tech Stack:** React, TanStack Query, Vitest + RTL, Tailwind ("The Print" tokens: `text-print-accent` kicker, `text-display` serif, `text-meta`/mono sub).

---

### Task 1: FocusPlotHeader component + test

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/FocusPlotHeader.tsx`
- Test: `packages/HimalayaUI/frontend/test/FocusPlotHeader.test.tsx`

`FocusPlotHeader` props:
```ts
interface FocusPlotHeaderProps {
  sampleName: string;        // serif title (display_name || name || "Sample N")
  sampleCode: string | null; // mono subline "smp" slug (sample.name)
  beamtime: string | null;   // experiment name
  exposureLabel: string | null; // representative exposure label
}
```

Renders:
- `data-testid="focus-plot-header"` wrapper
- kicker `data-testid="focus-plot-kicker"` text "Integration", `text-xs font-semibold uppercase tracking-wide text-print-accent`
- serif title `data-testid="focus-plot-title"` = sampleName, `text-display leading-tight`
- mono subline `data-testid="focus-plot-sub"` = join of present(`sampleCode`, `beamtime`, exposure phrase) with " · ", `font-mono text-xs text-ink-faint`. Exposure phrase = `representative exposure <exposureLabel>` when exposureLabel present.

- [ ] **Step 1: Write failing test** asserting kicker text, serif title text, and subline composition (all three segments; and degraded case where only sampleName present).
- [ ] **Step 2: Run, verify fail** (module not found).
- [ ] **Step 3: Implement FocusPlotHeader.**
- [ ] **Step 4: Run, verify pass.**

### Task 2: PlotCard headerSlot prop + test

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/PlotCard.tsx`
- Test: `packages/HimalayaUI/frontend/test/PlotCard.headerSlot.test.tsx`

Add `headerSlot?: JSX.Element` to `PlotCard` signature; thread to `TitleStrip`. In `TitleStrip`, when `headerSlot` is provided, render `{headerSlot}` on the left (no button, no `onTitleClick`); otherwise render the existing picker button unchanged. Keep the right-side controls cluster outside the branch.

- [ ] **Step 1: Write failing test** — render PlotCard wrapped in QueryClient/providers (use test-utils) with `headerSlot={<div data-testid="slot">X</div>}`; assert `slot` present and `plot-title` (picker button) absent. Second test: prop-less PlotCard still renders `plot-title` picker.
- [ ] **Step 2: Run, verify fail.**
- [ ] **Step 3: Implement headerSlot threading in PlotCard + TitleStrip.**
- [ ] **Step 4: Run, verify pass.**

### Task 3: Wire FocusPlotHeader into FocusWorkspaceLayout + test

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/FocusWorkspaceLayout.tsx`
- Test: `packages/HimalayaUI/frontend/test/FocusWorkspaceLayout.header.test.tsx`

In `FocusWorkspaceLayout`: add `useExperiment(experimentId ?? 0)` (gated, returns beamtime name) and read `activeExposureId` from app state + `useExposures(activeSampleId)` to resolve the representative exposure filename/label. Build a `FocusPlotHeader` from `corpusSample` (serif name = display_name||name||`Sample N`; sampleCode = name) and pass as `headerSlot` to `<PlotCard headerSlot={...} />`. Guard: only pass the slot when `corpusSample` is loaded (else pass nothing → PlotCard falls back gracefully, but in focus context the sample is always present once corpus loads).

- [ ] **Step 1: Write failing test** — mock queries, render FocusWorkspaceLayout at a sample route; assert `focus-plot-header` rendered with sample name, and the picker `plot-title` is NOT present.
- [ ] **Step 2: Run, verify fail.**
- [ ] **Step 3: Implement wiring.**
- [ ] **Step 4: Run, verify pass.**

### Task 4: Full verification

- [ ] `npm run build` (tsc --noEmit + vite) green.
- [ ] `npm test` green.
- [ ] Live screenshot of `/sample/2` vs mockup `.plate-head`.
- [ ] Commit per task; open PR.
