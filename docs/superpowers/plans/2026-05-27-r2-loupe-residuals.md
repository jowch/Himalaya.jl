# R2 — Loupe Fidelity Residuals Implementation Plan

> **For agentic workers:** Steps use checkbox (`- [ ]`) syntax for tracking. TDD throughout.

**Goal:** Close the four loupe affordance residuals from "The Print" fidelity audit — signal-strength meter (M-8), contact-sheet/loupe segmented view switch (M-9), subtitle "exposure N of M" (L-9), and bare tag tokens (L-10).

**Architecture:** All four are presentation-layer fixes in three frontend files. M-8 derives an honest signal proxy from the active exposure's detected peak count (more reflections = stronger ordering signal) via the existing `usePeaks` hook, rendering a 5-bar meter matching the mockup `.signal-bars`. M-9 adds a route-aware segmented `<Link>` pair in `CorpusTopbar` (the two existing routes `/samples` ↔ `/samples/loupe/:id` presented as a toggle; on `/samples` the loupe segment is disabled since no sample is selected). L-9 rewrites the loupe subtitle to `sample-id · exposure N of M`. L-10 drops the `key: value` rendering to a bare value token.

**Tech Stack:** React 18, TypeScript strict, TailwindCSS 4, Vitest + RTL, react-router-dom v6.

---

## File Structure

- `packages/HimalayaUI/frontend/src/components/LoupeSidebar.tsx` — add the signal-strength meter row (M-8) + bare tag tokens (L-10). Receives a `signalLevel: number` (0–5) prop computed by the parent.
- `packages/HimalayaUI/frontend/src/pages/LoupePage.tsx` — compute the signal level from `usePeaks(activeId)`, pass it down (M-8); rewrite the subtitle to `name · exposure N of M` (L-9).
- `packages/HimalayaUI/frontend/src/components/CorpusTopbar.tsx` — add the "Contact sheet | Loupe" segmented switch (M-9).
- Tests: `test/LoupeSidebar.test.tsx`, `test/LoupePage.test.tsx`, `test/CorpusTopbar.test.tsx`.

**Signal proxy decision (honest labeling):** The Exposure API carries no true signal metric. Peaks do carry `prominence`/`intensity`, and the *count* of detected peaks is a genuine, defensible proxy for ordering-signal strength (a strongly-ordered phase produces more resolvable Bragg reflections). The meter shows `min(peakCount, 5)` filled bars of 5, labeled `signal` exactly as the mockup. When peaks are still loading or absent, the meter shows 0 filled bars. This matches the mockup's 5-bar `.signal-bars` widget without inventing a fake metric.

---

## Task 1: Signal-strength meter in LoupeSidebar (M-8)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/LoupeSidebar.tsx`
- Test: `packages/HimalayaUI/frontend/test/LoupeSidebar.test.tsx`

- [ ] **Step 1: Write failing tests** — append a `describe` block to `LoupeSidebar.test.tsx`. Add `signalLevel` to `defaultProps()` (default `0`).

```tsx
describe("LoupeSidebar — signal meter (M-8)", () => {
  it("renders a 5-bar signal meter in the meta list", () => {
    render(<LoupeSidebar {...defaultProps()} signalLevel={3} />);
    const meter = screen.getByTestId("loupe-meta-signal");
    expect(meter).toBeInTheDocument();
    expect(meter.querySelectorAll("[data-bar]")).toHaveLength(5);
  });

  it("fills exactly signalLevel bars", () => {
    render(<LoupeSidebar {...defaultProps()} signalLevel={3} />);
    const meter = screen.getByTestId("loupe-meta-signal");
    expect(meter.querySelectorAll('[data-bar="on"]')).toHaveLength(3);
    expect(meter.querySelectorAll('[data-bar="off"]')).toHaveLength(2);
  });

  it("clamps signalLevel into the 0..5 range", () => {
    const { rerender } = render(
      <LoupeSidebar {...defaultProps()} signalLevel={9} />,
    );
    expect(
      screen.getByTestId("loupe-meta-signal").querySelectorAll('[data-bar="on"]'),
    ).toHaveLength(5);
    rerender(<LoupeSidebar {...defaultProps()} signalLevel={-2} />);
    expect(
      screen.getByTestId("loupe-meta-signal").querySelectorAll('[data-bar="on"]'),
    ).toHaveLength(0);
  });
});
```

Also update the shared `defaultProps()` to include `signalLevel: 0`.

- [ ] **Step 2: Run, verify fail** — `npm test -- LoupeSidebar` → FAIL (prop/testid missing).

- [ ] **Step 3: Implement** — add `signalLevel: number` to `Props`, destructure it, add a `SignalMeter` sub-component and a new `MetaRow`-style row. Bars: 5 `<i>` with `data-bar="on"|"off"`, on = `bg-ink-soft`, off = `bg-hair-strong`, sized `w-[5px] h-[11px] rounded-[1px]` per mockup `.signal-bars i`.

- [ ] **Step 4: Run, verify pass** — `npm test -- LoupeSidebar` → PASS.

- [ ] **Step 5: Commit** — `M-8: signal-strength meter in loupe sidebar`.

---

## Task 2: Derive signal level + pass it down (M-8 wiring)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/LoupePage.tsx`
- Test: `packages/HimalayaUI/frontend/test/LoupePage.test.tsx`

- [ ] **Step 1: Write failing test** — mock `usePeaks` in the LoupePage test (add to the `vi.mock("../src/queries", …)` factory). Add a `peaksQ` holder. Assert that with 3 peaks for the active exposure, the rendered meter has 3 filled bars.

```tsx
// in vi.hoisted holder: peaksQ: {} as { data?: { id: number }[] }
// in vi.mock factory: usePeaks: () => h.peaksQ,
// in each beforeEach: h.peaksQ = { data: [] };

describe("LoupePage — signal meter (M-8)", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    h.corpusQ = { data: [sample()], isLoading: false };
    h.exposuresQ = { data: [exposure({ id: 100 })], isLoading: false };
    h.experimentQ = { data: { id: 3, name: "Beamtime March", path: "/x" } };
    h.peaksQ = { data: [{ id: 1 }, { id: 2 }, { id: 3 }] };
  });

  it("fills the meter from the active exposure's peak count", () => {
    renderAt("/samples/loupe/7");
    expect(
      screen.getByTestId("loupe-meta-signal").querySelectorAll('[data-bar="on"]'),
    ).toHaveLength(3);
  });
});
```

- [ ] **Step 2: Run, verify fail** — `npm test -- LoupePage` → FAIL.

- [ ] **Step 3: Implement** — import `usePeaks`, call `usePeaks(activeId)`, compute `const signalLevel = peaksQ.data?.length ?? 0;`, pass `signalLevel={signalLevel}` to `<LoupeSidebar>`. Also pass `signalLevel={0}` in the `LOUPE_FIXTURE` so the boneyard fixture still type-checks.

- [ ] **Step 4: Run, verify pass** — `npm test -- LoupePage` → PASS.

- [ ] **Step 5: Commit** — `M-8: derive loupe signal level from peak count`.

---

## Task 3: Subtitle "exposure N of M" (L-9)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/LoupePage.tsx`
- Test: `packages/HimalayaUI/frontend/test/LoupePage.test.tsx`

- [ ] **Step 1: Write failing test**

```tsx
describe("LoupePage — subtitle (L-9)", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    h.corpusQ = { data: [sample()], isLoading: false };
    h.exposuresQ = {
      data: [exposure({ id: 100 }), exposure({ id: 101 })],
      isLoading: false,
    };
    h.experimentQ = { data: { id: 3, name: "Beamtime March", path: "/x" } };
    h.peaksQ = { data: [] };
  });

  it("shows sample-id · exposure N of M", () => {
    renderAt("/samples/loupe/7");
    const sub = screen.getByTestId("loupe-subtitle");
    expect(sub).toHaveTextContent("JC042");
    expect(sub).toHaveTextContent("exposure 1 of 2");
  });
});
```

- [ ] **Step 2: Run, verify fail** — FAIL (testid + text).

- [ ] **Step 3: Implement** — give the subtitle `<p>` `data-testid="loupe-subtitle"`. Compute `frameIndex` / `total` in LoupePage from `exposures` and `activeId`, render `{sample.name ?? "—"} · exposure {frameIndex + 1} of {total}`. Drop the `experimentName` from the subtitle (it stays computed only if still otherwise referenced; remove if now dead).

- [ ] **Step 4: Run, verify pass** — PASS. Also re-run the `identity` test which previously asserted `/Beamtime March/`; update it to assert the new subtitle (`exposure 1 of 1`).

- [ ] **Step 5: Commit** — `L-9: loupe subtitle reads sample-id · exposure N of M`.

---

## Task 4: Bare tag tokens (L-10)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/LoupeSidebar.tsx`
- Test: `packages/HimalayaUI/frontend/test/LoupeSidebar.test.tsx`

- [ ] **Step 1: Update the existing tag test** — the current test asserts `"lipid: DOPE"`. Change it to assert the bare value and the absence of the key form.

```tsx
it("renders tags as bare value tokens (L-10)", () => {
  const props = defaultProps();
  render(
    <LoupeSidebar
      {...props}
      sample={sample({
        tags: [{ id: 9, key: "lipid", value: "DOPE", source: "manual" }],
      })}
    />,
  );
  const tags = screen.getByTestId("loupe-tags");
  expect(tags).toHaveTextContent("DOPE");
  expect(tags).not.toHaveTextContent("lipid: DOPE");
  expect(screen.queryByLabelText("Remove lipid tag")).not.toBeInTheDocument();
});
```

- [ ] **Step 2: Run, verify fail** — FAIL (still renders `lipid: DOPE`).

- [ ] **Step 3: Implement** — change `{tag.key}: {tag.value}` to `{tag.value}`.

- [ ] **Step 4: Run, verify pass** — PASS.

- [ ] **Step 5: Commit** — `L-10: loupe tags render as bare value tokens`.

---

## Task 5: Segmented view switch in CorpusTopbar (M-9)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/CorpusTopbar.tsx`
- Test: `packages/HimalayaUI/frontend/test/CorpusTopbar.test.tsx`

**Behaviour:** Two segments, "Contact sheet" → `/samples` and "Loupe" → the loupe. The loupe route is sample-scoped (`/samples/loupe/:id`); the topbar has no sample id when on `/samples`. Resolution: when on a loupe route, both segments are active links and the loupe segment marks itself `on`; when on `/samples` (or any non-loupe route), the "Contact sheet" segment is `on` and the "Loupe" segment renders disabled (no sample selected to open). This is the route-based presentation the finding explicitly allows. The `?beamtime=` query is preserved on the Contact-sheet link so toggling back keeps the filter.

- [ ] **Step 1: Write failing tests** — append to `CorpusTopbar.test.tsx`.

```tsx
describe("CorpusTopbar — view switch (M-9)", () => {
  it("shows a Contact sheet | Loupe segmented switch", () => {
    mockExperiments();
    renderTopbar("/samples");
    const seg = screen.getByTestId("view-seg");
    expect(seg).toHaveTextContent("Contact sheet");
    expect(seg).toHaveTextContent("Loupe");
  });

  it("marks Contact sheet active on /samples and disables Loupe", () => {
    mockExperiments();
    renderTopbar("/samples");
    expect(screen.getByTestId("view-seg-sheet")).toHaveAttribute("data-active", "true");
    expect(screen.getByTestId("view-seg-loupe")).toBeDisabled();
  });

  it("marks Loupe active on a loupe route and links sheet back to /samples", () => {
    mockExperiments();
    renderTopbar("/samples/loupe/2");
    expect(screen.getByTestId("view-seg-loupe")).toHaveAttribute("data-active", "true");
    const sheet = screen.getByTestId("view-seg-sheet");
    expect(sheet).toHaveAttribute("href", "/samples");
  });

  it("preserves ?beamtime= on the Contact sheet link", () => {
    mockExperiments();
    renderTopbar("/samples/loupe/2?beamtime=1");
    expect(screen.getByTestId("view-seg-sheet")).toHaveAttribute(
      "href",
      "/samples?beamtime=1",
    );
  });
});
```

- [ ] **Step 2: Run, verify fail** — `npm test -- CorpusTopbar` → FAIL.

- [ ] **Step 3: Implement** — replace the trailing `<span className="flex-1" />` with `<span className="flex-1" />` followed by the segmented control. Detect loupe via `pathname.startsWith("/samples/loupe")`. Build the sheet href as `/samples` + preserved beamtime query. The loupe segment: when on a loupe route it's a non-link `on` button (already there); when on `/samples` it's a disabled button. Style per mockup `.seg`: container `flex bg-plate border border-hair-strong rounded-md overflow-hidden`; buttons `text-[11.5px] font-semibold px-3 py-1.5`; active `bg-ink text-paper`, inactive `text-ink-faint`.

- [ ] **Step 4: Run, verify pass** — `npm test -- CorpusTopbar` → PASS.

- [ ] **Step 5: Commit** — `M-9: contact-sheet/loupe segmented view switch in topbar`.

---

## Task 6: Full verification

- [ ] **Step 1:** `npm run build` → tsc + vite clean.
- [ ] **Step 2:** `npm test` → all green.
- [ ] **Step 3:** Live screenshots vs `sample-table.html` loupe (signal meter, segmented switch toggling /samples↔loupe, subtitle, bare tags). Save `/tmp/r2-*.png`.
- [ ] **Step 4:** request-pr-review → commit → push `r2-loupe-residuals` → open PR vs `main`.

---

## Self-Review

- **Spec coverage:** M-8 (Tasks 1+2), M-9 (Task 5), L-9 (Task 3), L-10 (Task 4). All four #225 scope items covered.
- **Out of scope honored:** no tag editing (#207), no contact-sheet/theme/phase changes, no big-frame ✕ (M-10 belongs to R2's finding list but issue #225 acceptance lists only meter/switch/subtitle/tag-value; M-10 is the grease-pencil X — NOT in #225's "Scope" section per the agent kit, which enumerates M-8/M-9/L-9/L-10 only).
- **Type consistency:** `signalLevel: number` prop named identically in LoupeSidebar Props, LoupePage call site, and the fixture.
