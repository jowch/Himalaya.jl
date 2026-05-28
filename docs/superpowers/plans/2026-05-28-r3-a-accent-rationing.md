# R3-A — Accent Rationing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restore terracotta accent rationing in both directions — repaint auto-peak triangles in the phase color of the index that claims them (terracotta only for the transient q-link `.hot` state), and flip the folio's "+ New series" button from ink to `button-accent`.

**Architecture:** Three surgical edits across three frontend files. (1) In `TraceViewer.tsx`, build a peak->phase-color map from the already-available `activeGroupIndices` prop, and use it to resolve each peak triangle's base color; reserve `var(--color-accent)` for the q-link hot state only. (2) In `PlotCard.tsx`, delete the terracotta "auto peak" legend swatch. (3) In `SeriesFolioPage.tsx`, swap the "+ New series" button's `bg-ink border-ink` for `bg-print-accent border-print-accent`. The TraceViewer change is the load-bearing one and is covered by two new Vitest cases asserting the SVG `fill` attribute (a data/SVG-attribute assertion, never a Tailwind class-string test).

**Tech Stack:** React + TypeScript, Vitest + Testing Library (JSDOM), Observable Plot SVG overlay, Tailwind with Print design tokens.

---

## Background — read before implementing

- **Issue:** `gh issue view 254 --repo jowch/Himalaya.jl` (the human-approved spec).
- **Findings:** `docs/2026-05-28-the-print-round3-findings.md` — R3-F01 (peaks over-use accent) + R3-Y01 ("+ New series" under-uses accent) + cross-surface pattern section 1.
- **Design rule:** `DESIGN.md` section 2 (terracotta = the grease pencil, hue 38, small fraction of any screen, for interaction) + section 5 (names "+ New series" as the canonical `button-accent` example).
- **Mockup:** `docs/redesign-mockups/focus-workspace.html:740-753` — peak markers (`.pk`) painted in the **phase color** of the claiming index via `claimOf(p.i)`, falling back to `var(--unindexed)` when unclaimed; `--accent` reserved for `.pk.hot`.

### Key facts about the existing code (verified against `main` @ d4f2df1)

- `TraceViewer.tsx` already receives `activeGroupIndices: IndexEntry[]` as a prop (`TraceViewerProps`, line 12) and it is already in the peak-drawing effect's dependency array (line 617). **No new props need threading** — this keeps the edit surgical and avoids touching the `PlotCard.tsx` call site of `<TraceViewer>`.
- `IndexEntry` (`api.ts:268-282`) has `peaks: IndexPeakRef[]`. `IndexPeakRef` (`api.ts:261`) is `{ peak_id, ratio_position, residual, q_observed }`. So a peak `p` is "claimed" by an index `ix` iff some `ix.peaks[*].peak_id === p.id`.
- `phaseColor(phase: string)` is already imported at `TraceViewer.tsx:4` from `../phases` and is used by `indexTicks` (line 148). It returns a CSS color string for a phase.
- The peak fill is resolved at `TraceViewer.tsx:409-441`. The current code sets `baseColor` from `peak.source` (auto -> `var(--color-accent)`, manual -> `var(--color-peak-manual)`). `baseColor` is reused in three downstream places: the hover halo stroke (line 461), the optimistic-placeholder stroke (line 483), and the final `fill` fallthrough (line 439). **All three must keep using the resolved `baseColor`** so the q-link hot halo still tracks the peak's own color.
- The q-link hot state is `litByQ` (line 450-451): `hoveredQ !== undefined && Math.abs(peak.q - hoveredQ) <= peak.q * Q_LINK_REL_TOL`. This is the ONLY state that should paint terracotta on a peak.
- The "auto peak" legend swatch is `PlotCard.tsx:684`: `<LegendItem symbol={<TriangleSvg color="var(--color-accent)" />} label="auto peak" />`.
- The "+ New series" button is `SeriesFolioPage.tsx:91-98`; the class string to change is on line 95: `className="rounded-md border border-ink bg-ink px-3 py-1.5 text-xs font-semibold text-paper hover:opacity-90"`.
- Existing TraceViewer tests (`test/TraceViewer.test.tsx`) query `[data-role="peak-root"] polygon` and assert on `getAttribute("fill")` / `getAttribute("fill-opacity")` / `getAttribute("data-peak-id")`. New tests follow this exact pattern. See `test/AGENTS.md`: data/SVG-attribute assertions only, never Tailwind class-string tests.

### Color-resolution decision (the load-bearing design choice)

Replace the current `baseColor` derivation. New precedence for a peak's `baseColor`:

1. **q-link hot** (`litByQ === true`) -> `var(--color-accent)` (terracotta). This is the only terracotta on a peak, and it is transient.
2. **Claimed by an active index** -> `phaseColor(claimingPhase)`, where `claimingPhase` is the phase of the first index in `activeGroupIndices` whose `peaks` contains this peak's id.
3. **Manual, unclaimed** -> `var(--color-peak-manual)` (preserve manual-peak identity; manual peaks are user marks, and the mockup's unclaimed fallback is for auto peaks the analysis didn't index).
4. **Auto, unclaimed** -> `var(--color-ink-faint)` (the mockup's `var(--unindexed)` equivalent; token exists at `styles.css:57`).

Rationale: the issue's Done-when says "Auto-peak triangles render in the phase color of their claiming index (or `--color-ink-faint` when unindexed); terracotta is reserved for the `.hot` q-link state." Manual peaks are out of the R3-F01 scope (the finding names only the auto-peak triangles + the "auto peak" legend swatch), so manual peaks keep `--color-peak-manual` to avoid scope creep. The `excluded`/`losing`/`dimOthers` branches (lines 419-437) keep their existing override behavior; they just consume the new `baseColor` where they currently consume the old one.

---

## File Structure

- Modify: `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx` — color-resolution helper + `baseColor` rewrite (peak fill block, ~lines 409-441; add a small `claimColorByPeakId` map inside the effect before the peak loop).
- Modify: `packages/HimalayaUI/frontend/src/components/PlotCard.tsx` — delete the "auto peak" legend swatch (line 684).
- Modify: `packages/HimalayaUI/frontend/src/pages/SeriesFolioPage.tsx` — flip "+ New series" button class (line 95).
- Test: `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx` — two new cases (claimed -> phase color; hot -> terracotta).

---

## Task 1: Repaint auto-peak triangles in phase color (TraceViewer)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx:409-441`
- Test: `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx`

- [ ] **Step 1: Write the failing tests**

Add these two cases to `test/TraceViewer.test.tsx` (inside the existing top-level `describe`, alongside the other peak-rendering tests). Reuse the existing render-helper conventions in that file — pass the same props the neighbouring tests pass (`trace`, `peaks`, `activeGroupIndices`, plus the default required props). Look at the existing "optimistic placeholder" test (~lines 207-254) for the exact prop shape; mirror it. The phase color is asserted by comparing against `phaseColor(...)` imported from `../src/phases`, NOT a hard-coded hex.

```tsx
import { phaseColor } from "../src/phases";

it("paints an auto peak in the phase color of the index that claims it", () => {
	const peaks = [
		{ id: 1, q: 0.1, source: "auto", excluded: false },
		// add any other Peak fields the neighbouring tests include
	];
	const claimingIndex = {
		// minimal IndexEntry: mirror the shape used by the tick tests (~line 357)
		id: 10,
		exposure_id: 1,
		phase: "Pn3m",
		basis: 0.1,
		score: null,
		r_squared: null,
		lattice_d: null,
		ngc: null,
		status: "candidate",
		kind: "auto",
		inputs_hash: null,
		peaks: [{ peak_id: 1, ratio_position: 0, residual: 0, q_observed: 0.1 }],
		predicted_q: [0.1],
	};
	const { container } = render(
		<TraceViewer
			/* ...same baseline props as the optimistic-placeholder test... */
			peaks={peaks}
			activeGroupIndices={[claimingIndex]}
		/>,
	);
	const tri = container.querySelector(
		'[data-role="peak-root"] polygon[data-peak-id="1"]',
	);
	expect(tri).not.toBeNull();
	expect(tri!.getAttribute("fill")).toBe(phaseColor("Pn3m"));
});

it("paints a peak in terracotta only when hoveredQ matches its q (q-link hot)", () => {
	const peaks = [
		{ id: 1, q: 0.1, source: "auto", excluded: false },
	];
	const claimingIndex = {
		id: 10, exposure_id: 1, phase: "Pn3m", basis: 0.1, score: null,
		r_squared: null, lattice_d: null, ngc: null, status: "candidate",
		kind: "auto", inputs_hash: null,
		peaks: [{ peak_id: 1, ratio_position: 0, residual: 0, q_observed: 0.1 }],
		predicted_q: [0.1],
	};
	const { container } = render(
		<TraceViewer
			/* ...same baseline props... */
			peaks={peaks}
			activeGroupIndices={[claimingIndex]}
			hoveredQ={0.1}
		/>,
	);
	const tri = container.querySelector(
		'[data-role="peak-root"] polygon[data-peak-id="1"]',
	);
	expect(tri!.getAttribute("fill")).toBe("var(--color-accent)");
});
```

> NOTE for the implementer: before writing, open `test/TraceViewer.test.tsx` and copy the EXACT baseline prop set + `Peak`/`IndexEntry` object shape from the nearest existing peak test (the optimistic-placeholder test ~207-254 and the tick test ~357). Do not invent fields. The two snippets above show only the fields load-bearing for this assertion; fill the rest from the neighbouring tests so the component renders.

- [ ] **Step 2: Run the tests to verify they fail**

Run (from `packages/HimalayaUI/frontend/`):
```bash
npx vitest run test/TraceViewer.test.tsx -t "phase color"
npx vitest run test/TraceViewer.test.tsx -t "terracotta only when hoveredQ"
```
Expected: the first FAILS — current code paints auto peaks `var(--color-accent)`, so `fill` is `var(--color-accent)`, not `phaseColor("Pn3m")`. The second may already pass (terracotta is the current auto color) — that is fine; it is a guard against the regression direction.

- [ ] **Step 3: Implement the color resolution**

In `TraceViewer.tsx`, inside the peak-drawing effect, just before the `for (const p of peaks)` loop at line 400, build a lookup from peak id -> claiming phase color:

```tsx
// R3-A (#254): peaks are painted in the phase color of the index that claims
// them (accent rationing — DESIGN.md section 2). Terracotta is reserved for the
// transient q-link hot state below; it is no longer a peak's resting color.
const claimColorByPeakId = new Map<number, string>();
for (const ix of activeGroupIndices) {
	for (const ref of ix.peaks) {
		if (!claimColorByPeakId.has(ref.peak_id)) {
			claimColorByPeakId.set(ref.peak_id, phaseColor(ix.phase));
		}
	}
}
```

Then replace the current base-color block at lines 410-414:

```tsx
const isAuto = peak.source === "auto";
// Bright/neon for "active workflow": auto = ice blue, manual = magenta.
const baseColor = isAuto
	? "var(--color-accent)"
	: "var(--color-peak-manual)";
```

with:

```tsx
const isAuto = peak.source === "auto";
// R3-A (#254): resting peak color follows the claiming index's phase
// (DESIGN.md section 2 accent rationing). Terracotta (`--color-accent`) is
// rationed to the q-link hot state only — see `litByQ` below. The old
// "bright/neon for active workflow" dark-era coloring is retired.
const litByQ = hoveredQ !== undefined
	&& Math.abs(peak.q - hoveredQ) <= peak.q * Q_LINK_REL_TOL;
const claimColor = claimColorByPeakId.get(peak.id);
const baseColor = litByQ
	? "var(--color-accent)"
	: claimColor
		? claimColor
		: isAuto
			? "var(--color-ink-faint)"
			: "var(--color-peak-manual)";
```

Then DELETE the now-duplicate `litByQ` declaration at lines 450-451 (the halo block still needs `litByQ`; it now reads the one declared above). The halo block (lines 449-465) becomes:

```tsx
const litByPeakId = hoveredPeakId === peak.id;
if (litByPeakId || litByQ) {
```

(i.e. remove the `const litByQ = ...` line at 450-451, keep everything else — `litByPeakId` stays, the `if` keeps both conditions).

Leave the `excluded`/`losing`/`dimOthers` branches (lines 419-437), the halo stroke (`baseColor`, line 461), the optimistic stroke (`baseColor`, line 483), and the final `fill = baseColor` (line 439) untouched — they all consume the new `baseColor` automatically.

- [ ] **Step 4: Run the tests to verify they pass**

Run (from `packages/HimalayaUI/frontend/`):
```bash
npx vitest run test/TraceViewer.test.tsx
```
Expected: PASS — including the two new cases and ALL pre-existing TraceViewer tests (the dimOthers / excluded / losing / optimistic / tick tests must stay green — those branches are unchanged).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/TraceViewer.tsx \
        packages/HimalayaUI/frontend/test/TraceViewer.test.tsx
git commit -m "R3-A (#254): paint auto-peak triangles in claiming-phase color, ration terracotta to q-link hot"
```

---

## Task 2: Remove the terracotta "auto peak" legend swatch (PlotCard)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/PlotCard.tsx:684`

- [ ] **Step 1: Delete the swatch line**

In `PlotCard.tsx`, in `PlotLegend` (lines 678-700), delete line 684:

```tsx
<LegendItem symbol={<TriangleSvg color="var(--color-accent)" />} label="auto peak" />
```

Leave the `manual peak`, `excluded`, and `predicted <phase>` legend items as-is. (The `excluded` swatch at line 689 also uses `var(--color-accent)`; it is NOT named in the issue's Done-when and is out of scope for this PR — do not touch it. Flagging only.)

- [ ] **Step 2: Verify the build typechecks and the bundle compiles**

Run (from `packages/HimalayaUI/frontend/`):
```bash
npm run build
```
Expected: PASS (`tsc --noEmit` clean, vite build succeeds). No legend test should reference the removed "auto peak" item; if one does, the build/test will surface it.

- [ ] **Step 3: Run the frontend unit suite for any legend assertions**

Run (from `packages/HimalayaUI/frontend/`):
```bash
npm test
```
Expected: PASS. If a pre-existing test asserted the "auto peak" legend text, update it in this commit (it is now intentionally removed per the issue Done-when).

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/PlotCard.tsx
git commit -m "R3-A (#254): drop terracotta 'auto peak' legend swatch"
```

---

## Task 3: Flip "+ New series" to button-accent (SeriesFolioPage)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/SeriesFolioPage.tsx:95`

- [ ] **Step 1: Swap the button color classes**

In `SeriesFolioPage.tsx`, on the `+ New series` button (lines 91-98), change line 95 from:

```tsx
className="rounded-md border border-ink bg-ink px-3 py-1.5 text-xs font-semibold text-paper hover:opacity-90"
```

to:

```tsx
className="rounded-md border border-print-accent bg-print-accent px-3 py-1.5 text-xs font-semibold text-paper hover:opacity-90"
```

(`text-paper` stays — the accent fill carries white-paper text per DESIGN.md section 5; only the `border-ink bg-ink` -> `border-print-accent bg-print-accent` swap is in scope. This is the surgical change the issue names; do NOT hoist into a shared `<AccentButton>` in this PR — the issue lists that as an optional alternative and it would widen the diff that wave-2 must rebase on.)

- [ ] **Step 2: Verify the build and unit suite**

Run (from `packages/HimalayaUI/frontend/`):
```bash
npm run build && npm test
```
Expected: PASS. The `series-folio-new` test id is unchanged, so the existing SeriesFolioPage tests keep matching by `data-testid`, not by color class.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/pages/SeriesFolioPage.tsx
git commit -m "R3-A (#254): flip '+ New series' to button-accent (bg-print-accent)"
```

---

## Final verification (the gate before PR)

- [ ] **Frontend gate (REQUIRED):** from `packages/HimalayaUI/frontend/`, run BOTH:
```bash
npm run build
npm test
```
Both must pass with zero failures. Capture the tail of each.

- [ ] **Visual "<=2 terracotta marks" Done-when** (orchestrator verifies in the live harness; this plan documents how it is checked): on a focus-indexed view (a sample with >=1 active index and >=5 auto peaks), the rendered trace plate must show **<=2 terracotta marks** (excluding the transient `.hot` q-link state, which only appears on hover). After this change, a resting focus-indexed trace should show **zero** permanent terracotta peak triangles — every claimed auto peak renders in its phase color, every unclaimed auto peak renders `--color-ink-faint`, and the only terracotta on the surface is the legitimate grease-pencil controls (e.g. reject actions), which are <=2 by the rationing rule. The folio's "+ New series" now carries its one sanctioned piece of accent. Verification method: live Playwright screenshot of the focus-indexed surface + a DOM count of `[data-role="peak-root"] polygon[fill="var(--color-accent)"]` at rest (must be 0) — the orchestrator runs this; the Vitest data-* assertions in Task 1 are the in-repo proxy.

---

## Self-review notes

- **Spec coverage:** Done-when (1) phase-color auto peaks + ink-faint unindexed -> Task 1. (2) terracotta reserved for `.hot` -> Task 1 (`litByQ`). (3) PlotLegend "auto peak" swatch removed -> Task 2. (4) "+ New series" -> `bg-print-accent` -> Task 3. (5) <=2 terracotta per view -> Final verification. (6) Vitest covers (a) claimed->phase color, (b) terracotta only when `hoveredQ===peak.q` -> Task 1 Step 1.
- **Surgical-edit constraint:** Task 1 adds a local map + rewrites one `baseColor` block inside the existing effect; no new props, no change to the `<TraceViewer>` call site in PlotCard. Task 2 deletes one line. Task 3 changes one class string. Wave-2 (#257) rebases cleanly.
- **No Tailwind class-string test:** Task 1 asserts the SVG `fill` attribute (a data/SVG-attribute assertion), per `test/AGENTS.md`.
- **Type consistency:** `claimColorByPeakId: Map<number,string>`, keyed on `IndexPeakRef.peak_id` (number) and `Peak.id` (number); `phaseColor(ix.phase): string`. `litByQ` declared once, consumed by `baseColor` and the halo `if`.
