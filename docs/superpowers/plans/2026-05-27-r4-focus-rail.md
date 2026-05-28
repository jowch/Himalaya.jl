# R4 — Focus rail as output Implementation Plan

> **For agentic workers:** Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rework the focus-workspace rail into the mockup's "output" surface — a serif "Phase call" block above multi-select "Candidate indexings", with per-phase series ratios, phase-colour hover-preview that dims losing peaks on a swap, and no legacy Miller inset.

**Architecture:** A new `seriesRatio` pure helper derives the √N : √N : … string from a phase + the claimed `ratio_position`s. `PhaseCallBlock` renders the output cards; `PhasePanel` is reworked from +/− `IndexCard` buttons into multi-select `.cand` rows. The rail's `IndicesCard` drops the `MillerPlot` inset. Hover-preview already flows through `hoveredIndexId` → TraceViewer phase-coloured ticks (L-11 trace side is satisfied); we add rail-side phase-colour styling and extend the dim-losing-peaks behaviour so a candidate that would *evict* an active phase dims those orphaned peaks while hovered.

**Tech Stack:** React 18, TypeScript strict, TailwindCSS 4 (Print tokens), Zustand, TanStack Query, Vitest.

---

## File map

- Create `src/lib/seriesRatio.ts` — phase → ratio-series string from claimed positions.
- Create `test/seriesRatio.test.ts` — unit tests for the helper.
- Modify `src/state.ts` — add `previewIndexId` + `setPreviewIndex` (hover-preview channel distinct from the existing `hoveredIndexId`; also used to compute losing peaks).
- Rewrite `src/components/PhasePanel.tsx` — Phase-call block + multi-select candidate rows.
- Modify `src/components/IndicesCard.tsx` — drop the MillerPlot inset.
- Modify `src/components/TraceViewer.tsx` — add "losing peak" dim when a previewed candidate would evict an active phase.
- Modify `src/components/PlotCard.tsx` — pass preview + losing info to TraceViewer.
- Modify `test/PhasePanel.test.tsx` — new assertions (multi-select, phase-call, ratio).
- Modify `test/IndicesCard*`/layout tests if they assert the Miller panel.

---

### Task 1: `seriesRatio` helper

**Files:**
- Create: `src/lib/seriesRatio.ts`
- Test: `test/seriesRatio.test.ts`

The mockup shows `√2 : √3 : √4` (Pn3m), `1 : 2 : 3` (Lamellar). The canonical radicands per phase mirror `src/phase.jl` `phaseratios`. The index's claimed positions come from `peaks[].ratio_position` (1-based). We format the radicands for the claimed positions (deduped, sorted), capped at 4 terms with an ellipsis.

- [ ] **Step 1: Write the failing test**

```ts
import { describe, it, expect } from "vitest";
import { seriesRatio } from "../src/lib/seriesRatio";

describe("seriesRatio", () => {
  it("formats cubic Pn3m radicands as √N terms", () => {
    expect(seriesRatio("Pn3m", [1, 2, 3])).toBe("√2 : √3 : √4");
  });
  it("formats lamellar as bare integers", () => {
    expect(seriesRatio("Lamellar", [1, 2, 3])).toBe("1 : 2 : 3");
  });
  it("formats hexagonal mixed integer/radical", () => {
    expect(seriesRatio("Hexagonal", [1, 2, 3])).toBe("1 : √3 : 2");
  });
  it("collapses perfect squares to integers (√4 → 2)", () => {
    expect(seriesRatio("Im3m", [1, 2])).toBe("√2 : 2");
  });
  it("dedupes and sorts positions", () => {
    expect(seriesRatio("Lamellar", [3, 1, 1, 2])).toBe("1 : 2 : 3");
  });
  it("caps at four terms with an ellipsis", () => {
    expect(seriesRatio("Lamellar", [1, 2, 3, 4, 5])).toBe("1 : 2 : 3 : 4 …");
  });
  it("returns empty string for no positions or unknown phase", () => {
    expect(seriesRatio("Lamellar", [])).toBe("");
    expect(seriesRatio("Nonsense", [1])).toBe("");
  });
  it("ignores positions beyond the known series length", () => {
    expect(seriesRatio("Pn3m", [1, 999])).toBe("√2");
  });
});
```

- [ ] **Step 2: Run it; expect FAIL (module missing).**

Run: `npm test -- seriesRatio`

- [ ] **Step 3: Implement**

```ts
/**
 * seriesRatio — format a phase's defining reflection ratio series as
 * "√2 : √3 : √4" (cubic) or "1 : 2 : 3" (lamellar), restricted to the
 * positions an index actually claims. Mirrors `src/phase.jl phaseratios`:
 * the per-phase array of radicands h²+k²+l² (the value UNDER the radical).
 *
 * R4 finding L-9/L-10: the rail surfaces this ratio per phase + per candidate.
 */

// Radicands (values under √) per phase, position-indexed (1-based via the
// `ratio_position` field). Source of truth: src/phase.jl.
const RADICANDS: Record<string, number[]> = {
  Lamellar:  [1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121],
  Hexagonal: [1, 3, 4, 7, 9, 11, 12, 13, 16, 19],
  Square:    [1, 2, 4, 5, 8, 9, 10, 13, 16, 17, 18, 20],
  Pn3m:      [2, 3, 4, 6, 8, 9, 10, 11, 12, 14, 16],
  Im3m:      [2, 4, 6, 8, 10, 12, 14, 16, 18, 20],
  Ia3d:      [6, 8, 14, 16, 20, 22, 24, 26],
  Fm3m:      [3, 4, 8, 11, 12],
  Fd3m:      [3, 8, 11, 12, 16, 19, 24, 27, 32, 35, 36],
};

const MAX_TERMS = 4;

function term(radicand: number): string {
  const root = Math.sqrt(radicand);
  // Perfect square → integer (√4 → 2); else √N.
  return Number.isInteger(root) ? String(root) : `√${radicand}`;
}

export function seriesRatio(phase: string, positions: number[]): string {
  const radicands = RADICANDS[phase];
  if (!radicands) return "";
  const claimed = [...new Set(positions)]
    .filter((p) => p >= 1 && p <= radicands.length)
    .sort((a, b) => a - b);
  if (claimed.length === 0) return "";
  const trimmed = claimed.slice(0, MAX_TERMS);
  const out = trimmed.map((p) => term(radicands[p - 1]!)).join(" : ");
  return claimed.length > MAX_TERMS ? `${out} …` : out;
}
```

- [ ] **Step 4: Run it; expect PASS.** `npm test -- seriesRatio`
- [ ] **Step 5: Commit.** `feat(focus-rail): seriesRatio helper for phase ratio series (#227 L-9/L-10)`

---

### Task 2: `previewIndexId` state channel

**Files:**
- Modify: `src/state.ts`

The existing `hoveredIndexId` drives the trace's phase-coloured strong-tick preview (L-11 trace side already correct). We add a parallel `previewIndexId` set by the candidate rows so the rail's own hover-preview + the "losing peaks" dim are computed from a single channel. We set BOTH on hover (keep `hoveredIndexId` so the trace strong-tick keeps working) — `previewIndexId` is the one that also computes evicted peaks.

Decision: reuse `hoveredIndexId` rather than add a second id — they always move together. So **no new state field**; instead the losing-peaks computation lives in `PlotCard` (which already resolves `hoveredIndex`). This task is a no-op placeholder folded into Task 4/5. Skip.

---

### Task 3: Rewrite PhasePanel into the output rail

**Files:**
- Modify: `src/components/PhasePanel.tsx`
- Create (inline in same file): `PhaseCallBlock`, `CandidateRow`
- Test: `test/PhasePanel.test.tsx`

The rail has two sections:
1. **Phase call** (`.phasecall`): a card with one `.pc-block` per active-set phase — serif name + mono score + lattice + score bar + series ratio; a `Coexistence · N phases` tag header when >1.
2. **Candidate indexings** (`.cand` rows): every auto index as a multi-select checkbox row (`.c-mark` box, `.c-name`, `explains N peaks`, score + mini bar). Clicking toggles membership (add/remove from the active group). Hovering sets `hoveredIndexId`. Active rows show the checked box + `· in the call`. Speculative indices keep their existing add/remove/delete handling but render as the same `.cand` rows with a delete affordance, plus the "+ Add speculative" button.

- [ ] **Step 1: Write failing tests** — replace the alternatives `+ button` assertions with multi-select checkbox + phase-call + ratio assertions.

```tsx
// in test/PhasePanel.test.tsx, ADD a new describe block:
describe("<PhasePanel> — phase-call output block (R4 L-9/L-10)", () => {
  it("renders a Phase call block per active-set phase with serif name, score and series ratio", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 0.89,
          r_squared: 0.998, lattice_d: 197, ngc: -1.5, status: "candidate", kind: "auto",
          predicted_q: [0.045, 0.055, 0.064],
          peaks: [
            { peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.045 },
            { peak_id: 2, ratio_position: 2, residual: 0, q_observed: 0.055 },
            { peak_id: 3, ratio_position: 3, residual: 0, q_observed: 0.064 },
          ] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const block = await screen.findByTestId("phase-call-block-10");
    expect(block).toHaveTextContent("Pn3m");
    expect(block).toHaveTextContent("0.89");
    expect(block).toHaveTextContent("√2 : √3 : √4");
  });

  it("shows a Coexistence header when two phases are in the call", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 0.89,
          r_squared: 0.99, lattice_d: 197, ngc: -1.5, status: "candidate", kind: "auto",
          predicted_q: [0.045], peaks: [{ peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.045 }] },
        { id: 11, exposure_id: 42, phase: "Lamellar", basis: 0.3, score: 0.96,
          r_squared: 0.99, lattice_d: 61, ngc: null, status: "candidate", kind: "auto",
          predicted_q: [0.103], peaks: [{ peak_id: 4, ratio_position: 1, residual: 0, q_observed: 0.103 }] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10, 11] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    await screen.findByTestId("phase-call-block-10");
    expect(screen.getByTestId("coexistence-tag")).toHaveTextContent(/Coexistence.*2 phases/i);
  });

  it("shows an empty phase-call message when no index is in the call", async () => {
    mockAll(
      [
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
          r_squared: 0.71, lattice_d: 9.1, ngc: null, status: "candidate", kind: "auto",
          predicted_q: [0.4], peaks: [{ peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.4 }] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    expect(await screen.findByTestId("phase-call-empty")).toBeInTheDocument();
  });
});

describe("<PhasePanel> — candidate multi-select (R4 L-10)", () => {
  it("renders candidates as checkboxes reflecting active-set membership", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 0.89,
          r_squared: 0.99, lattice_d: 197, ngc: -1.5, status: "candidate", kind: "auto",
          predicted_q: [0.045], peaks: [{ peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.045 }] },
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
          r_squared: 0.71, lattice_d: 9.1, ngc: null, status: "candidate", kind: "auto",
          predicted_q: [0.4], peaks: [{ peak_id: 2, ratio_position: 1, residual: 0, q_observed: 0.4 }] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const inCall = await screen.findByRole("checkbox", { name: /Pn3m/i });
    expect(inCall).toBeChecked();
    const candidate = screen.getByRole("checkbox", { name: /Im3m/i });
    expect(candidate).not.toBeChecked();
  });

  it("toggling an unchecked candidate posts add-to-group", async () => {
    vi.spyOn(global, "fetch").mockImplementation(async (input) => {
      const u = typeof input === "string" ? input : (input as Request).url;
      if (u.endsWith("/indices")) return new Response(JSON.stringify([
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 0.89, r_squared: 0.99,
          lattice_d: 197, ngc: -1.5, status: "candidate", kind: "auto", predicted_q: [0.045],
          peaks: [{ peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.045 }] },
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6, r_squared: 0.71,
          lattice_d: 9.1, ngc: null, status: "candidate", kind: "auto", predicted_q: [0.4],
          peaks: [{ peak_id: 2, ratio_position: 1, residual: 0, q_observed: 0.4 }] },
      ]), { status: 200, headers: { "Content-Type": "application/json" } });
      if (u.endsWith("/groups")) return new Response(JSON.stringify([
        { id: 2, exposure_id: 42, kind: "custom", active: true, members: [10] },
      ]), { status: 200, headers: { "Content-Type": "application/json" } });
      if (u.endsWith("/api/groups/2/members")) return new Response(JSON.stringify({
        id: 2, exposure_id: 42, kind: "custom", active: true, members: [10, 11],
      }), { status: 200, headers: { "Content-Type": "application/json" } });
      return new Response("not found", { status: 404 });
    });
    renderWithProviders(<PhasePanel exposureId={42} />);
    const candidate = await screen.findByRole("checkbox", { name: /Im3m/i });
    fireEvent.click(candidate);
    await waitFor(() => {
      const spy = global.fetch as unknown as { mock: { calls: unknown[][] } };
      const urls = spy.mock.calls.map((c) => typeof c[0] === "string" ? c[0] : (c[0] as Request).url);
      expect(urls).toContain("/api/groups/2/members");
    });
  });

  it("hovering a candidate sets hoveredIndexId; leaving clears it", async () => {
    useAppState.setState({ hoveredIndexId: undefined });
    mockAll(
      [
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6, r_squared: 0.71,
          lattice_d: 9.1, ngc: null, status: "candidate", kind: "auto", predicted_q: [0.4],
          peaks: [{ peak_id: 2, ratio_position: 1, residual: 0, q_observed: 0.4 }] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const row = await screen.findByTestId("candidate-row-11");
    fireEvent.mouseEnter(row);
    expect(useAppState.getState().hoveredIndexId).toBe(11);
    fireEvent.mouseLeave(row);
    expect(useAppState.getState().hoveredIndexId).toBeUndefined();
  });
});
```

Delete the now-obsolete `+ button` / `−`-based tests ("renders alternative indices with a + button", "renders each active-group index … remove button", "clicking remove calls DELETE …", "clicking + on an alternative posts …", "hovering an alternative sets hoveredIndexId", "dims alternatives with r_squared below 0.98", "does not dim alternatives with null r_squared") — replaced by the multi-select equivalents above. Keep the κ test (still valid: candidate row may not show κ, so rework that test to check the phase-call block shows lattice instead). Keep the "no exposure" test.

- [ ] **Step 2: Run; expect FAIL.** `npm test -- PhasePanel`

- [ ] **Step 3: Rewrite `PhasePanel.tsx`** (full new file below).

- [ ] **Step 4: Run; expect PASS.** `npm test -- PhasePanel`

- [ ] **Step 5: Commit.** `feat(focus-rail): phase-call block + multi-select candidates (#227 L-9/L-10)`

---

### Task 4: Drop the Miller inset from the rail

**Files:**
- Modify: `src/components/IndicesCard.tsx`
- Test: layout/IndicesCard tests that assert `miller-panel`.

- [ ] **Step 1:** grep tests for `miller-panel` / `MillerPlot` in the focus layout; relax/remove assertions that the rail contains a Miller panel. (MillerPlot component + its own test stay — only the rail usage is dropped.)
- [ ] **Step 2:** Simplify `IndicesCard` to render only `<PhasePanel/>` (no MillerPlot, no `useMemo` of activeGroupIndices/hoveredIndex).
- [ ] **Step 3: Run** `npm test -- IndicesCard FocusWorkspace` — expect PASS.
- [ ] **Step 4: Commit.** `feat(focus-rail): drop legacy Miller inset from the Print rail (#227 L-12)`

---

### Task 5: Phase-colour hover-preview + dim losing peaks on swap

**Files:**
- Modify: `src/components/PlotCard.tsx`, `src/components/TraceViewer.tsx`
- Test: `test/TraceViewer*` or a focused new test.

The trace strong-tick for the hovered index is already phase-coloured (`indexTicks` uses `phaseColor`). What's missing (L-11): when hovering a candidate that competes with an active-set phase (claims an overlapping peak), the peaks the swap would *orphan* should dim. We compute, in `PlotCard`, the set of peak ids that would be evicted if the hovered candidate were committed, and pass it to TraceViewer to render those peak triangles dimmed.

- [ ] **Step 1: Write failing test** — a `TraceViewer` test asserting that peaks in `losingPeakIds` get reduced opacity. (Use the existing TraceViewer test harness pattern; if none, add `data-losing` attr to peak nodes and assert.)
- [ ] **Step 2: Run; FAIL.**
- [ ] **Step 3: Implement** — add `losingPeakIds?: Set<number>` prop to TraceViewer; in the peak-triangle loop, if `losingPeakIds?.has(peak.id)` set `opacity` low + `data-losing`. Compute the set in PlotCard from `hoveredIndex` vs active members (overlap → evicted active peaks not in hovered claim).
- [ ] **Step 4: Run; PASS.**
- [ ] **Step 5: Commit.** `feat(focus-rail): dim losing peaks on candidate swap preview (#227 L-11)`

---

### Task 6: Verify

- [ ] `npm run build` (tsc --noEmit + vite build) green.
- [ ] `npm test` green (full suite).
- [ ] Live screenshots per the live-verification runbook.
- [ ] Commit any test fixups.
