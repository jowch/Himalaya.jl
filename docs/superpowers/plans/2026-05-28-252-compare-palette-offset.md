# COMPARE_PALETTE_LIGHT phase-offset ≥13° invariant Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Re-tune the two `COMPARE_PALETTE_LIGHT` hues that violate the ≥13° phase-offset floor (263→251, 295→279), factor the shared `angularHueDistance` helper into `lib/color-distance.ts`, and add a `presets.test.ts` that enforces the floor against the light palette — so the next palette nudge can't silently re-introduce a phase/member color conflation in figure export (especially the #208 heatmap fill surface).

**Architecture:** GitHub issue #252, Option B (recommended). The on-screen `COMPARE_PALETTE` (`lib/comparison/coloring.ts`) already holds the ≥13° floor, proven by `test/coloring.test.ts`. Its light-background sibling `COMPARE_PALETTE_LIGHT` (`lib/figure-export/presets.ts`) inherited the *old* hue assignments and two entries sit below the floor: `[6] 263` is 1° from Lamellar (264) and `[7] 295` is 5° from Ia3d (300). We promote the test-local `angularHueDistance` helper into production `lib/color-distance.ts`, point `coloring.test.ts` at it (no behavior change, just dedup), add `test/presets.test.ts` asserting the same floor on `COMPARE_PALETTE_LIGHT`, then re-tune the two hues to make it green. The replacement hues (251, 279) mirror the on-screen palette's same-slot hue layout so the on-screen→export hue mapping stays stable.

**Tech Stack:** TypeScript, Vitest (JSDOM), OKLCH color space. Frontend at `packages/HimalayaUI/frontend/`.

---

## Background facts (verified against `main` in this worktree)

`PHASE_PALETTE` hues (`src/phases.ts:32-41`): Pn3m 58, Im3m 162, Ia3d 300, Fm3m 18, Fd3m 318, Hexagonal 350, Lamellar 264, Square 132.

Current `COMPARE_PALETTE_LIGHT` hues (`src/lib/figure-export/presets.ts:58-71`):
`33, 77, 147, 175, 200, 220, 263, 295, 285, 333, 3, 105`.

Angular distance of each to its nearest phase hue (shortest-arc, same formula the test uses):

| idx | hue | nearest phase | distance | status |
|-----|-----|---------------|----------|--------|
| 0 | 33 | Fm3m 18 | 15° | ok |
| 1 | 77 | Pn3m 58 | 19° | ok |
| 2 | 147 | Square 132 | 15° | ok |
| 3 | 175 | Im3m 162 | 13° | ok (at floor) |
| 4 | 200 | Im3m 162 | 38° | ok |
| 5 | 220 | Lamellar 264 | 44° | ok |
| 6 | **263** | **Lamellar 264** | **1°** | **VIOLATION** |
| 7 | **295** | **Ia3d 300** | **5°** | **VIOLATION** |
| 8 | 285 | Ia3d 300 | 15° | ok |
| 9 | 333 | Fd3m 318 | 15° | ok |
| 10 | 3 | Fm3m 18 | 15° | ok |
| 11 | 105 | Square 132 | 27° | ok |

Replacement hues (verified):
- `[6] 263 → 251`: vs Lamellar 264 = **13°** (at floor, ok); vs Im3m 162 = 89°; vs Ia3d 300 = 49°. All ≥13°.
- `[7] 295 → 279`: vs Ia3d 300 = **21°**; vs Lamellar 264 = 15°; vs Fd3m 318 = 39°. All ≥13°.

**IMPORTANT — divergence from the issue text to resolve during implementation.** The issue's Option A says "Move `[7] 295 → 285` (already at 285 in the on-screen palette per #208 r1 fix)". But on-screen `COMPARE_PALETTE` (`coloring.ts:69-82`) actually has slot **[7] = 279** and slot **[8] = 285**. A literal "mirror the on-screen palette" sets light `[7] → 279`, not 285. Setting light `[7] → 285` makes light `[7]` and `[8]` an *identical string* (`oklch(0.50 0.14 285)` both), collapsing two adjacent legend swatches into one. **This plan uses `[7] → 279` to truly mirror the on-screen hue layout and keep slots [7]/[8] distinct, satisfying both the issue's "keep on-screen→export hue mapping stable" intent and the "all sit ≥13°" acceptance criterion.** Flag to the reviewer if 285 is specifically required; the test only enforces the floor, so either choice is test-safe.

So the two edits to `COMPARE_PALETTE_LIGHT`:
- `[6]` hue `263 → 251`
- `[7]` hue `295 → 279`

Both ≥13° from every phase hue; slots [7]=279 / [8]=285 stay distinct, matching the on-screen palette's hue layout.

`angularHueDistance` currently lives only in `test/coloring.test.ts:117-121` (test-local). No production module exports it. `parseOklch` (`:111-115`) is NOT in scope for extraction — leave it test-local; the new `color-distance.ts` takes plain numbers.

`COMPARE_PALETTE_LIGHT` importers (no changes needed — pure value retune): `src/lib/figure-export/adapters/multiTraceAdapter.ts:5,57`; referenced by `src/lib/figure-export/marks/multiTraceExportMarks.ts:100`.

Test location convention (`test/AGENTS.md`): unit tests live flat in `frontend/test/`, mirroring `src/`. New file → `test/presets.test.ts`. Never assert on Tailwind class strings (N/A here — pure color math).

Single-file Vitest run: `node_modules/.bin/vitest run test/<file>` from `packages/HimalayaUI/frontend/`.

---

## File Structure

- **Create** `src/lib/color-distance.ts` — one tiny pure function `angularHueDistance(a, b)`: shortest angular distance between two hues in degrees. Production module so both `coloring.test.ts` and the new `presets.test.ts` import one proven implementation (DRY; the issue mandates this factoring).
- **Create** `test/color-distance.test.ts` — pins the helper's behavior (wraparound, symmetry, identity, known regression anchors).
- **Create** `test/presets.test.ts` — asserts every `COMPARE_PALETTE_LIGHT` hue sits ≥13° from every `PHASE_PALETTE` hue (mirrors the on-screen block in `coloring.test.ts:77-97`), plus a string-equality non-collision check.
- **Modify** `test/coloring.test.ts:117-121` — delete the test-local `angularHueDistance`, import it from `../src/lib/color-distance` (dedup; no assertion change).
- **Modify** `src/lib/figure-export/presets.ts:58-71` — retune slot `[6] 263→251` and slot `[7] 295→279`; update the docstring (`:52-56`) to drop the "tracked but out of scope" disclaimer and state the floor now holds + is tested.

---

## Task 1: Extract `angularHueDistance` into a production helper

**Files:**
- Create: `src/lib/color-distance.ts`
- Test: `test/color-distance.test.ts`

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/frontend/test/color-distance.test.ts`:

```ts
/**
 * Unit tests for the shared angular-hue-distance helper extracted from
 * coloring.test.ts (issue #252). Load-bearing for the phase-offset ≥13°
 * invariant asserted in both coloring.test.ts (on-screen COMPARE_PALETTE) and
 * presets.test.ts (export COMPARE_PALETTE_LIGHT).
 */
import { describe, it, expect } from "vitest";
import { angularHueDistance } from "../src/lib/color-distance";

describe("angularHueDistance", () => {
  it("is zero for identical hues", () => {
    expect(angularHueDistance(120, 120)).toBe(0);
  });

  it("is symmetric", () => {
    expect(angularHueDistance(10, 200)).toBe(angularHueDistance(200, 10));
  });

  it("takes the shortest arc across the 0/360 wraparound", () => {
    // 350 -> 10 is 20deg the short way, not 340deg the long way.
    expect(angularHueDistance(350, 10)).toBe(20);
    expect(angularHueDistance(10, 350)).toBe(20);
  });

  it("never exceeds 180deg", () => {
    expect(angularHueDistance(0, 270)).toBe(90);
    expect(angularHueDistance(0, 180)).toBe(180);
  });

  it("matches the known phase-offset distances (regression anchors)", () => {
    // From issue #252: light palette violations vs phase hues.
    expect(angularHueDistance(263, 264)).toBe(1);   // [6] vs Lamellar
    expect(angularHueDistance(295, 300)).toBe(5);   // [7] vs Ia3d
    // ...and the chosen replacements clear the floor.
    expect(angularHueDistance(251, 264)).toBe(13);  // [6]->251 vs Lamellar
    expect(angularHueDistance(279, 300)).toBe(21);  // [7]->279 vs Ia3d
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`):
```bash
node_modules/.bin/vitest run test/color-distance.test.ts
```
Expected: FAIL — module `../src/lib/color-distance` does not exist (resolve/import error).

- [ ] **Step 3: Write minimal implementation**

Create `packages/HimalayaUI/frontend/src/lib/color-distance.ts`:

```ts
/**
 * Shortest angular distance between two hues, in degrees (mod 360, range
 * [0, 180]). Extracted from `test/coloring.test.ts` (issue #252) so the
 * phase-offset ≥13° invariant can run against both the on-screen
 * `COMPARE_PALETTE` (`lib/comparison/coloring.ts`) and the export
 * `COMPARE_PALETTE_LIGHT` (`lib/figure-export/presets.ts`) from one proven
 * implementation. Phase/member colors visually conflate when their hues sit
 * too close; this is the load-bearing perceptual check behind that floor.
 */
export function angularHueDistance(a: number, b: number): number {
  const d = (Math.abs(((a - b) % 360) + 540) % 360) - 180;
  return Math.abs(d);
}
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```bash
node_modules/.bin/vitest run test/color-distance.test.ts
```
Expected: PASS (5 tests).

- [ ] **Step 5: Commit**

```bash
git add src/lib/color-distance.ts test/color-distance.test.ts
git commit -m "Extract angularHueDistance into lib/color-distance (#252)"
```

---

## Task 2: Point `coloring.test.ts` at the shared helper (dedup)

**Files:**
- Modify: `test/coloring.test.ts`

- [ ] **Step 1: Replace the test-local helper with an import**

In `packages/HimalayaUI/frontend/test/coloring.test.ts`, add the import (after the existing `phases` import, currently around `:21`):

```ts
import { angularHueDistance } from "../src/lib/color-distance";
```

Then delete the test-local definition (currently `:117-121`):

```ts
/** Shortest angular distance between two hues in degrees, mod 360. */
function angularHueDistance(a: number, b: number): number {
  const d = Math.abs(((a - b) % 360) + 540) % 360 - 180;
  return Math.abs(d);
}
```

Leave `parseOklch`, `oklchToLinearSrgb`, `relativeLuminance`, `contrastRatio` exactly as they are — only `angularHueDistance` moves. The ≥13° assertion block (`:77-97`) is unchanged; it now calls the imported helper.

- [ ] **Step 2: Run the test to verify it still passes (no behavior change)**

Run:
```bash
node_modules/.bin/vitest run test/coloring.test.ts
```
Expected: PASS (22 tests) — identical count to baseline; pure refactor.

- [ ] **Step 3: Commit**

```bash
git add test/coloring.test.ts
git commit -m "coloring.test: use shared angularHueDistance helper (#252)"
```

---

## Task 3: Add the failing `presets.test.ts` invariant (RED)

**Files:**
- Test: `test/presets.test.ts`

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/frontend/test/presets.test.ts`:

```ts
/**
 * Figure-export palette invariant (issue #252).
 *
 * `COMPARE_PALETTE_LIGHT` is the light-background sibling of the on-screen
 * `COMPARE_PALETTE`. Since PR #208 it also drives heatmap fill colors in
 * figure export, so a member trace/fill must never read as perceptually
 * identical to a phase color — the same ≥13° hue-offset floor the on-screen
 * palette holds (see test/coloring.test.ts) applies here.
 *
 * Floor rationale: at L≈0.55 C≈0.13 a 13deg hue shift is dE2000 ≈ 2.5 — the
 * smallest gap that keeps byPhase vs bySample/distinct from conflating. The
 * eight phase hues pack the warm + purple sectors tightly enough that >=15deg
 * everywhere AND twelve mutually-distinct entries is infeasible, so 13deg (not
 * 15deg) is the floor. Keep this floor numeric and in lockstep with both
 * palettes' docstrings if hues are re-laid.
 */
import { describe, it, expect } from "vitest";
import { COMPARE_PALETTE_LIGHT } from "../src/lib/figure-export/presets";
import { PHASE_PALETTE } from "../src/phases";
import { angularHueDistance } from "../src/lib/color-distance";

/** Parse the hue channel out of an `oklch(L C h)` string. */
function oklchHue(s: string): number {
  const m = /oklch\(\s*[\d.]+\s+[\d.]+\s+([\d.-]+)\s*\)/i.exec(s);
  if (!m) throw new Error(`not an oklch() string: ${s}`);
  return parseFloat(m[1]!);
}

const PHASE_OFFSET_FLOOR_DEG = 13;

describe("COMPARE_PALETTE_LIGHT — phase-offset invariant (#252)", () => {
  it("exports a palette of 10–12 colors", () => {
    expect(COMPARE_PALETTE_LIGHT.length).toBeGreaterThanOrEqual(10);
    expect(COMPARE_PALETTE_LIGHT.length).toBeLessThanOrEqual(12);
  });

  it("no entry exactly equals a phase color (string-equality first-line check)", () => {
    const phaseColors = new Set(Object.values(PHASE_PALETTE));
    for (const c of COMPARE_PALETTE_LIGHT) {
      expect(phaseColors.has(c), `${c} collides with a PHASE_PALETTE entry`).toBe(false);
    }
  });

  it("every entry sits >=13deg from every PHASE_PALETTE hue", () => {
    const phaseHues = Object.values(PHASE_PALETTE).map(oklchHue);
    for (const c of COMPARE_PALETTE_LIGHT) {
      const ch = oklchHue(c);
      for (const ph of phaseHues) {
        const dist = angularHueDistance(ch, ph);
        expect(
          dist,
          `LIGHT ${c} (hue ${ch}) vs PHASE hue ${ph}`,
        ).toBeGreaterThanOrEqual(PHASE_OFFSET_FLOOR_DEG);
      }
    }
  });
});
```

- [ ] **Step 2: Run the test to verify it fails on the two known violations**

Run:
```bash
node_modules/.bin/vitest run test/presets.test.ts
```
Expected: FAIL on the ">=13deg" case — the message names `LIGHT oklch(0.55 0.12 263) (hue 263) vs PHASE hue 264` (distance 1) and `LIGHT oklch(0.50 0.14 295) (hue 295) vs PHASE hue 300` (distance 5). The first two `it` blocks pass.

- [ ] **Step 3: Commit the RED test**

```bash
git add test/presets.test.ts
git commit -m "test: assert COMPARE_PALETTE_LIGHT phase-offset floor (RED) (#252)"
```

---

## Task 4: Re-tune the two offending hues (GREEN)

**Files:**
- Modify: `src/lib/figure-export/presets.ts:58-71` (palette) and `:42-57` (docstring)

- [ ] **Step 1: Retune the two palette entries**

In `packages/HimalayaUI/frontend/src/lib/figure-export/presets.ts`, change slot `[6]`:

```ts
  "oklch(0.55 0.12 263)", // lavender
```
to:
```ts
  "oklch(0.55 0.12 251)", // lavender (was 263 = 1deg from Lamellar 264; #252)
```

and change slot `[7]`:

```ts
  "oklch(0.50 0.14 295)", // purple
```
to:
```ts
  "oklch(0.50 0.14 279)", // purple (was 295 = 5deg from Ia3d 300; #252)
```

(Keep the L/C channels exactly as-is — only the hue moves. Slot [8] stays `oklch(0.50 0.14 285)`, so [7]=279 / [8]=285 mirror the on-screen palette's hue layout and stay visually distinct.)

- [ ] **Step 2: Update the docstring to drop the out-of-scope disclaimer**

In the same file, replace the phase-offset paragraph in the `COMPARE_PALETTE_LIGHT` docstring (currently `:52-56`):

```
 * **Phase-offset parity (PR #251 r1).** Entry [8] was relocated from hue 315
 * to 285 to keep parity with `COMPARE_PALETTE`'s same-index fix (3° collision
 * with Fd3m 318). Several other entries in this light variant — 175, 263,
 * 295, 3 — sit closer to phase hues than the on-screen ≥13° floor; tracked
 * but out of scope for #208 (issue covers the on-screen render-core).
```

with:

```
 * **Phase-offset invariant (#252).** Every entry sits >=13deg from every
 * `PHASE_PALETTE` hue, matching the on-screen `COMPARE_PALETTE` floor — since
 * PR #208 this palette also drives heatmap fill colors, so a member trace/fill
 * must never read as perceptually identical to a phase color. Entry [8] was
 * relocated 315->285 (3deg from Fd3m 318) in PR #251 r1; entries [6] 263->251
 * (1deg from Lamellar 264) and [7] 295->279 (5deg from Ia3d 300) were re-tuned
 * in #252. `test/presets.test.ts` pins the floor — keep it green if you
 * re-tune hues here.
```

- [ ] **Step 3: Run the invariant test to verify it now passes (GREEN)**

Run:
```bash
node_modules/.bin/vitest run test/presets.test.ts
```
Expected: PASS (3 tests).

- [ ] **Step 4: Commit**

```bash
git add src/lib/figure-export/presets.ts
git commit -m "Re-tune COMPARE_PALETTE_LIGHT hues to clear >=13deg phase floor (#252)"
```

---

## Task 5: Full-suite verification + prod build gate

**Files:** none (verification only)

- [ ] **Step 1: Run the full Vitest suite**

Run (from `packages/HimalayaUI/frontend/`):
```bash
npm test > /tmp/vitest-252.out 2>&1
grep -E "FAIL|passed|failed|Test Files" /tmp/vitest-252.out
```
Expected: all test files pass; no FAIL lines. `color-distance.test.ts`, `coloring.test.ts`, and `presets.test.ts` all green.

- [ ] **Step 2: Run the production build gate (REQUIRED — Vitest alone misses Vite prod-build breakage)**

Run:
```bash
npm run build > /tmp/build-252.out 2>&1; echo "exit=$?"; tail -20 /tmp/build-252.out
```
Expected: `exit=0`. `npm run build` is `tsc --noEmit && vite build`; both the type-check (catches the new `color-distance.ts` import wiring) and the bundle must succeed.

- [ ] **Step 3: Confirm no remaining out-of-scope disclaimer or stale helper**

Run:
```bash
grep -rn "tracked but out of scope" src/lib/figure-export/presets.ts; echo "---"
grep -rn "function angularHueDistance" test/ src/
```
Expected: first grep returns nothing (disclaimer removed). Second grep returns NOTHING (the helper is now an exported arrow/function in `src/lib/color-distance.ts` declared as `export function angularHueDistance` — confirm only that single production definition exists and no test-local copy remains; adjust the grep to `grep -rn "angularHueDistance" test/ src/` to eyeball the import sites).

---

## Manual visual confirmation (acceptance criterion)

The issue's final acceptance bullet asks to re-export a representative figure with at least one Lamellar and one Ia3d trace + heatmap mode and confirm no member-vs-phase visual conflation. This is a manual/human-in-the-loop step performed at PR review against a running app (`bin/himalaya serve` + Compare page -> figure export -> heatmap). NOT automatable in this TDD loop; flag it in the PR description for the human reviewer to spot-check, backed by the now-green numeric floor in `presets.test.ts`. The numeric invariant is the durable guard; the visual check is a one-time confirmation.

---

## Self-Review notes

- **Spec coverage:** (1) "all entries >=13deg from every phase hue" -> Task 4 retune, proven by Task 3 test. (2) "a test asserts this floor, factored to share with on-screen" -> Tasks 1-3 (shared `angularHueDistance`, both `presets.test.ts` and `coloring.test.ts` import it). (3) "docstring updated to remove out-of-scope disclaimer" -> Task 4 Step 2 + Task 5 Step 3 grep. (4) "re-export visual confirmation" -> Manual section. All four acceptance bullets covered.
- **Floors not counts:** palette-length check uses a 10-12 range (matching `coloring.test.ts`); the invariant is a `>=13` floor, not an exact-distance assertion — a future re-tune that keeps the floor stays green.
- **Frontend gate:** Task 5 includes `npm run build`, not Vitest alone, per the orchestrator constraint.
- **No Tailwind class assertions:** N/A — pure OKLCH math.
- **Open question flagged for reviewer:** light `[7]` set to **279** (true mirror of on-screen `COMPARE_PALETTE[7]`) rather than the issue's literal **285** (which would duplicate slot [8] as an identical string). Either clears the floor; 279 keeps slots [7]/[8] distinct. See the "IMPORTANT — divergence" note above. If the reviewer requires 285, change one literal in Task 4 Step 1 and the test still passes.
