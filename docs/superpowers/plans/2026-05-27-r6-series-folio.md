# R6 — Series folio: live mini-waterfall + per-sample phase strip — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:test-driven-development per task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Turn folio cards into real miniatures — a live per-series mini-waterfall SVG + a per-sample phase strip with transition caption — and add the missing folio chrome (recipe/draft pills, filter chips, 3-way sort, "+ New series" action, provenance footer, serif titles).

**Architecture:** The folio listing (`SeriesSummary`) carries only top-3 phases, not an ordered per-sample sequence, so each card fetches its full `Series` detail via `useSeries(id)` and reads `members[].snapshot` (ordered by `display_order`). The mini-waterfall is a pure read-only SVG generator that reuses the mockup `figSVG()` geometry (log-x mapping, stacked baselines, fill+line+ticks) fed by each member's `snapshot.effective_peaks` and `confirmed_index` — NOT the heavy interactive `MultiTracePlot` (one Observable-Plot instance per card across a masonry is the wrong tool). The phase strip + caption derive from the ordered `confirmed_index.phase` per member. Filter/sort logic is pure functions on the listing.

**Tech Stack:** React 18, TypeScript strict, TanStack Query, Vitest + RTL, SVG, Tailwind v4 "The Print" tokens.

---

## File structure

- **Create** `src/lib/series/folioFigure.ts` — pure geometry: `buildMiniWaterfall(members) -> {viewBox, traces[]}` and `buildPhaseStrip(members) -> {segments[], caption}`. Unit-tested in isolation (no DOM).
- **Create** `src/lib/series/folioFilter.ts` — pure `filterSort(series, {search, filter, sort}) -> SeriesSummary[]`.
- **Create** `src/components/SeriesMiniWaterfall.tsx` — read-only SVG component rendering `buildMiniWaterfall` output.
- **Create** `src/components/SeriesPhaseStrip.tsx` — cell-per-sample strip + caption from `buildPhaseStrip`.
- **Modify** `src/components/SeriesFolioCard.tsx` — fetch detail, compose figure + strip + pills + serif title + footer.
- **Modify** `src/pages/SeriesFolioPage.tsx` — filter chips, 3-way sort, "+ New series" action, serif headline.
- **Modify/add tests:** `test/folioFigure.test.ts`, `test/folioFilter.test.ts`, `test/SeriesMiniWaterfall.test.tsx`, `test/SeriesPhaseStrip.test.tsx`, `test/SeriesFolioCard.test.tsx`, `test/SeriesFolioPage.test.tsx`.

**Data note:** "cross-experiment" is not derivable from the available data (no sample→experiment join on the folio path). The filter chip is present per F-D but matches on a `crossExperiment` flag that is currently always false (graceful: the chip yields an empty state on this corpus). Provenance footer (F-H) shows member count + ordering variable + edited timestamp + author; the literal beamtime string is not available, so we show ordering-variable provenance, which IS available.

---

## Task 1: Phase-strip + waterfall geometry (pure)

**Files:**
- Create: `src/lib/series/folioFigure.ts`
- Test: `test/folioFigure.test.ts`

Types (define in `folioFigure.ts`):

```ts
import type { SeriesMember } from "../../api";

export interface StripSegment {
  /** dominant phase, or null = unindexed */
  phase: string | null;
  /** when two phases coexist, the second phase for a 2-stop gradient */
  coexistWith: string | null;
}
export interface PhaseStripModel {
  segments: StripSegment[];
  /** "transition" (first ≠ last), "throughout" (single phase), "none" (no calls) */
  kind: "transition" | "throughout" | "none";
  first: string | null;
  last: string | null;
}

/** Dominant phase of a member = its confirmed_index.phase (null if unindexed). */
export function memberPhase(m: SeriesMember): string | null {
  return m.snapshot?.confirmed_index?.phase ?? null;
}

export function buildPhaseStrip(members: SeriesMember[]): PhaseStripModel { /* … */ }
```

Coexistence: the snapshot exposes a single `confirmed_index`, so true 2-phase coexistence isn't representable from one member. We model coexistence as: when a member's `snapshot.effective_peaks` includes peaks NOT in `confirmed_index.peak_ids` AND the neighbours differ — too speculative. **Decision (YAGNI):** drive coexistence purely off data we have — a member coexists when it has a confirmed index AND its immediate predecessor or successor has a *different* confirmed phase is a transition, not coexistence. Since one member = one confirmed phase, we render each segment as a single phase; the 2-stop gradient is reserved for a future multi-index snapshot and is exercised by a unit test against a synthetic member whose snapshot we extend with a `coexist_phase` only in tests. To keep the gradient code path real and tested without inventing backend fields, `buildPhaseStrip` accepts an optional `coexistResolver(m) -> string | null`; production passes none (always single), the test passes one. This keeps the gradient rendering wired (mockup F-B) and covered, without faking data in production.

- [ ] **Step 1: Write failing tests**

```ts
import { describe, it, expect } from "vitest";
import { buildPhaseStrip, buildMiniWaterfall } from "../src/lib/series/folioFigure";
import type { SeriesMember } from "../src/api";

function member(phase: string | null, peaks: number[] = [], order = 0): SeriesMember {
  return {
    id: order + 1, series_id: 1, exposure_id: order + 100, display_order: order,
    band_height: 1, y_offset: 0, normalization: "max", color_override: null,
    label_override: null, q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: phase === null && peaks.length === 0 ? null : {
      effective_peaks: peaks.map((q, i) => ({ id: i + 1, q, intensity: 1, sharpness: 1, source: "auto" as const })),
      confirmed_index: phase === null ? null : {
        id: order + 1, phase, lattice_d: 100, r_squared: 0.99, ngc: 3, peak_ids: peaks.map((_, i) => i + 1),
      },
      analysis_inputs_hash: "h",
    },
    is_stale: false, created_by: 1, created_at: null,
  };
}

describe("buildPhaseStrip", () => {
  it("emits one segment per member in order", () => {
    const s = buildPhaseStrip([member("Pn3m", [0.04], 0), member("Lamellar", [0.1], 1)]);
    expect(s.segments).toHaveLength(2);
    expect(s.segments[0]!.phase).toBe("Pn3m");
    expect(s.segments[1]!.phase).toBe("Lamellar");
  });
  it("classifies a transition (first ≠ last)", () => {
    const s = buildPhaseStrip([member("Pn3m", [0.04], 0), member("Lamellar", [0.1], 1)]);
    expect(s.kind).toBe("transition");
    expect(s.first).toBe("Pn3m"); expect(s.last).toBe("Lamellar");
  });
  it("classifies a single phase throughout", () => {
    const s = buildPhaseStrip([member("Pn3m", [0.04], 0), member("Pn3m", [0.041], 1)]);
    expect(s.kind).toBe("throughout");
  });
  it("classifies no clear phase when no member is indexed", () => {
    const s = buildPhaseStrip([member(null, [], 0), member(null, [], 1)]);
    expect(s.kind).toBe("none");
    expect(s.segments[0]!.phase).toBeNull();
  });
  it("uses first/last INDEXED member for the caption, skipping unindexed ends", () => {
    const s = buildPhaseStrip([member(null, [], 0), member("Pn3m", [0.04], 1), member("Lamellar", [0.1], 2), member(null, [], 3)]);
    expect(s.first).toBe("Pn3m"); expect(s.last).toBe("Lamellar"); expect(s.kind).toBe("transition");
  });
  it("renders a coexistence gradient when a resolver supplies a second phase", () => {
    const s = buildPhaseStrip([member("Pn3m", [0.04, 0.1], 0)], (m) => "Lamellar");
    expect(s.segments[0]!.coexistWith).toBe("Lamellar");
  });
});

describe("buildMiniWaterfall", () => {
  it("emits one trace path per member with a baseline y stacked top→bottom", () => {
    const wf = buildMiniWaterfall([member("Pn3m", [0.04], 0), member("Lamellar", [0.1], 1)]);
    expect(wf.traces).toHaveLength(2);
    // baselines descend down the SVG: later members sit lower (larger y)
    expect(wf.traces[1]!.baselineY).toBeGreaterThan(wf.traces[0]!.baselineY);
    // each trace has a non-empty line path and a color
    expect(wf.traces[0]!.linePath.startsWith("M")).toBe(true);
    expect(wf.traces[0]!.color).toMatch(/oklch|var/);
  });
  it("places peak ticks at log-mapped x within the viewBox width", () => {
    const wf = buildMiniWaterfall([member("Pn3m", [0.04, 0.057], 0)]);
    const t = wf.traces[0]!;
    expect(t.ticks.length).toBeGreaterThan(0);
    for (const tk of t.ticks) {
      expect(tk.x).toBeGreaterThanOrEqual(0);
      expect(tk.x).toBeLessThanOrEqual(wf.width);
    }
  });
  it("scales viewBox height with member count", () => {
    const a = buildMiniWaterfall([member("Pn3m", [0.04], 0)]);
    const b = buildMiniWaterfall([member("Pn3m", [0.04], 0), member("Pn3m", [0.04], 1), member("Pn3m", [0.04], 2)]);
    expect(b.height).toBeGreaterThan(a.height);
  });
  it("handles an empty member list without throwing", () => {
    const wf = buildMiniWaterfall([]);
    expect(wf.traces).toHaveLength(0);
    expect(wf.height).toBeGreaterThan(0);
  });
});
```

- [ ] **Step 2: Run, expect fail** — `npm test -- folioFigure` → module not found.

- [ ] **Step 3: Implement `folioFigure.ts`.** Geometry mirrors mockup `figSVG()`: `W=340, stepY=31, amp=24, padT=padB=12, padL=padR=14`; log-x map `xOf(q)=padL+(ln q-ln Qmin)/(ln Qmax-ln Qmin)*PW` with `Qmin=0.02, Qmax=0.42`; intensity model `intensity(q,peaks)=0.26*(0.02/q)^0.9 + Σ amp_i*0.6^i*gauss`. Use `effective_peaks` (q + intensity, normalized to the member's max peak intensity → amp) as the Bragg comps; color = `phaseColor(confirmed_index.phase)` or `var(--color-ink-faint)` if unindexed. Build `linePath`, `fillPath` (line closed to baseline), and `ticks` (first ≤3 peaks). Output `{width, height, traces:[{baselineY, linePath, fillPath, color, ticks:[{x, h, color, opacity}]}]}`.

- [ ] **Step 4: Run, expect pass** — `npm test -- folioFigure`.

- [ ] **Step 5: Commit** — `feat(folio): pure mini-waterfall + phase-strip geometry (R6 F-A/F-B #229)`.

## Task 2: Filter/sort logic (pure)

**Files:** Create `src/lib/series/folioFilter.ts`; Test `test/folioFilter.test.ts`.

```ts
export type FolioFilter = "all" | "transition" | "cross";
export type FolioSort = "recent" | "variable" | "size";
export interface FolioControls { search: string; filter: FolioFilter; sort: FolioSort; }
export interface FolioRow { /* SeriesSummary + derived */ }
export function filterSort(series: SeriesSummary[], c: FolioControls): SeriesSummary[];
```

- "transition" filter: `member_phase_count > 1` (proxy for "has a transition" usable from the listing — a series spanning >1 distinct phase). "cross": `false` for now (no data; documented). "recent": preserve backend order. "variable": no ordering var on summary → sort by title as a stable proxy is wrong; instead sort by `member_count` desc is "size". For "variable" we need the ordering variable, absent on the summary. **Decision:** "variable" sorts by title `localeCompare` (the only stable text key on the summary), labeled "Variable" per mockup; documented as a listing-data limitation. "size" sorts `member_count` desc.

- [ ] **Step 1: failing tests** covering: search title match (case-insensitive); "transition" keeps only `member_phase_count>1`; "cross" yields empty; "size" orders by member_count desc; "recent" preserves input order; empty search returns all.
- [ ] **Step 2: run, fail.**
- [ ] **Step 3: implement.**
- [ ] **Step 4: run, pass.**
- [ ] **Step 5: commit** — `feat(folio): filter chips + 3-way sort logic (R6 F-D #229)`.

## Task 3: SeriesMiniWaterfall + SeriesPhaseStrip components

**Files:** Create `src/components/SeriesMiniWaterfall.tsx`, `src/components/SeriesPhaseStrip.tsx`; Tests `test/SeriesMiniWaterfall.test.tsx`, `test/SeriesPhaseStrip.test.tsx`.

- `SeriesMiniWaterfall({members})`: renders `buildMiniWaterfall`'s SVG. `data-testid="series-mini-waterfall"`, `data-trace-count`. Classes from mockup: `.wf-base` dashed hair baseline, `.wf-fill` opacity 0.09, `.wf-line` stroke 1.6, `.wf-tick`. Use inline styles / Tailwind utilities (no asserting on classes in tests — assert on `data-*` + element counts).
- `SeriesPhaseStrip({members})`: renders `buildPhaseStrip`. `data-testid="series-phase-strip"`, one `[data-testid="ps-seg"]` per member, gradient via inline `background`. Caption `[data-testid="ps-cap"]` shows "Pn3m → Lamellar" / "Pn3m throughout" / "No clear phase". Phase names colored via `phaseColor`.

- [ ] **Step 1: failing tests** — waterfall: renders an `<svg>`, `data-trace-count` matches member count, N `path.wf-line` (query by tag/data-testid). strip: N `ps-seg`, caption text for transition/throughout/none, gradient segment has `linear-gradient` in style when coexist resolver fires (component test can render a member array; coexistence in component is single-phase so test the throughout/transition/none captions + seg count, and unit-test gradient in folioFigure Task 1).
- [ ] **Step 2: run, fail.**
- [ ] **Step 3: implement both components.**
- [ ] **Step 4: run, pass.**
- [ ] **Step 5: commit** — `feat(folio): SeriesMiniWaterfall + SeriesPhaseStrip read-only components (R6 F-A/F-B #229)`.

## Task 4: Rework SeriesFolioCard (detail fetch, pills, serif, footer)

**Files:** Modify `src/components/SeriesFolioCard.tsx`; Test `test/SeriesFolioCard.test.tsx`.

Card now: fetches `useSeries(series.id)`; while loading or on a draft with zero members, falls back to a phase-swatch strip from `series.member_phases` (keeps the cold/no-data path graceful). When detail loads with ≥1 member, renders `SeriesMiniWaterfall` + `SeriesPhaseStrip`.

Chrome:
- `card-kick`: `fig-n` kicker = "Recipe" if draft else `Fig. {n}` (n = a stable per-card ordinal passed as a prop `figNumber?`). Pill: draft → `pill-draft` "Draft"; else `newMatches>0` → `pill-new` "+N new match". `newMatches` is not on the summary → **prop `newMatches?: number` defaulting 0** (no data source yet; structurally present per F-C, always 0 in production, exercised in tests). The draft pill replaces the bare lowercase "draft" text.
- Title: serif via `.text-headline` role (or a `card-title` sized serif). `data-testid` unchanged.
- Meta: `{member_count} samples · by {ordering_variable ?? "—"}` (ordering var from detail).
- Footer (`card-foot`): rule + provenance (member count / ordering var) on the left, `edited {relative} · {author}` on the right. Relative time from `updated_at ?? last_event_at`.

Add a small `relativeDays(iso)` + `formatEdited(days)` helper inline (mockup mapping: 0→"just now",1→"yesterday",<7→"N days ago",<14→"1 week ago",else "2 weeks ago"). Keep `mockNow` injectable for deterministic tests via an optional `now?: Date` prop.

- [ ] **Step 1: failing tests** — extend `SeriesFolioCard.test.tsx` (mock `useSeries`):
  - draft (`content_hash:""`) → `pill-draft` present (testid `series-card-:id-pill`), data-draft true.
  - `newMatches={2}` non-draft → pill text "+2 new match".
  - detail with 2 indexed members → `series-mini-waterfall` + `series-phase-strip` present, NOT the swatch fallback.
  - no detail / loading → swatch fallback (`series-card-:id-swatches`) present.
  - footer shows author + an edited string (inject `now`).
  - title still renders; onOpen still fires.
- [ ] **Step 2: run, fail.**
- [ ] **Step 3: implement.** Mock `useSeries` in the test via `vi.mock("../src/queries", …)` returning `{ data, isLoading }`.
- [ ] **Step 4: run, pass.**
- [ ] **Step 5: commit** — `feat(folio): card miniature + pills + serif title + provenance footer (R6 F-A/B/C/H,X-1 #229)`.

## Task 5: Rework SeriesFolioPage (chips, sort, +New series, serif headline)

**Files:** Modify `src/pages/SeriesFolioPage.tsx`; Test `test/SeriesFolioPage.test.tsx`.

- Header kicker "Folio" + serif `<h1>` "Saved series" (`.text-headline`) + sub sentence. Count uses serif numeral.
- Controls row: search (kept) + filter chips `All / Has transition / Cross-experiment` (`data-testid="series-folio-chip-all|transition|cross"`, `aria-pressed`) + 3-way sort seg `Recent / Variable / Largest` (`series-folio-sort-recent|variable|size`). Replace the old `recency|title` toggle. Drive via `filterSort` from Task 2.
- Primary action: header-right `+ New series` button (`btn-solid` style: `bg-ink text-paper`) → `navigate("/series/new")`, `data-testid="series-folio-new"`.
- Pass a `figNumber` (1-based index among non-draft visible cards) to each card.

- [ ] **Step 1: failing tests** — update `SeriesFolioPage.test.tsx`:
  - chips render + clicking "Has transition" filters to `member_phase_count>1`.
  - sort "Largest" orders by member_count desc.
  - "+ New series" navigates to `/series/new` (route marker).
  - serif headline present (`series-folio-heading` testid).
  - existing search + empty + error + card-per-series tests still pass (adjust the removed `sort-title` testid references).
  - Mock `useSeries` (used by the card) so cards render the fallback without network.
- [ ] **Step 2: run, fail.**
- [ ] **Step 3: implement.**
- [ ] **Step 4: run, pass.**
- [ ] **Step 5: commit** — `feat(folio): filter chips, 3-way sort, +New series, serif header (R6 F-D/F-G,X-1 #229)`.

## Task 6: Full suite + build + live verification

- [ ] `npm test` green (whole frontend suite).
- [ ] `npm run build` green (tsc --noEmit + vite build).
- [ ] Live: `VITE_API_PORT=8091 npm run dev --host 127.0.0.1 --port 5193`; Playwright screenshot `/series`; compare to `series-folio.html`. Verify live waterfall (not grey box), per-sample strip + caption, pills, chips + sorts, "+ New series", serif titles, footer.
- [ ] Tear down dev server; commit any screenshot-driven fixes.

## Self-review notes
- Spec coverage: F-A (Task1/3/4), F-B (Task1/3/4), F-C (Task4), F-D (Task2/5), F-G (Task5), F-H (Task4), X-1 (Task4/5). All covered.
- Data limits documented: cross-experiment chip + newMatches pill + "Variable" sort are structurally present but data-starved on this corpus — wired, defaulted, and tested with injected data.
