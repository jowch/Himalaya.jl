# R7 — Series Scoping Proposal UX Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the bare manual-entry list at `/series/new` with the machine-proposes / human-confirms worksheet from `series-scoping.html` — a centred worksheet plate with an autogroup summary, a parsed "Ordered by" field, per-row trace sparklines, amber "check the read" confidence flags with confirm-not-fill-out ink values, a "Himalaya also found" loose-match section, a preview phase strip, reorder/undo affordances, and a narrative gate footer.

**Architecture:** Keep the existing mutation surface untouched (`useScopeSeries.mutate({ key, tags })` and `scopeSeriesMutator`). All new behaviour is (a) pure derivation modules under `src/lib/scoping/` (parse a numeric sort key from a value string; split loose-matches off the proposal) that are unit-tested in isolation, and (b) a rewritten `SeriesScopingPage` + new presentational components that *surface* `proposeOrdering`'s output plus per-member trace/index data fetched through the EXISTING `useStableQueryMap` plumbing (`useMemberTraces` + a new sibling `useMemberIndices`). No new backend routes. No change to `proposeOrdering`'s core key-selection logic.

**Tech Stack:** React 18 + TypeScript strict, TanStack Query, Zustand, Tailwind v4 "The Print" tokens, Vitest + Testing Library, boneyard-js skeletons, Observable-Plot-free hand-rolled SVG sparklines.

---

## Data-source decisions (read before starting)

The picker projection (`PickerSampleRow`: `{ sample, indexing_exposure_id, all_exposures }`) carries **no** trace or phase data. The mockup's sparklines and phase strip are driven by per-sample phase composition. Resolution — reuse existing plumbing, add no routes:

1. **Sparklines** — fetch the live `(q,I)` trace per row via the existing `useMemberTraces(exposureIds)` helper, keyed on each row's `indexing_exposure_id`. Rows with `indexing_exposure_id === null` render an empty sparkline frame.
2. **Dominant phase (sparkline colour + preview strip)** — fetch indices per row via a NEW sibling hook `useMemberIndices(exposureIds)` that mirrors `useMemberTraces` exactly (same `useStableQueryMap` pattern, `queryKeys.indices(id)` / `api.listIndices`). Dominant phase = the index with the highest `score`; coexistence = top two distinct phases. Rows with no indices → `unindexed` grey.
3. **Numeric ordering ("low to high")** — `proposeOrdering` returns value strings only. A NEW pure `parseSortKey(value: string): number | null` extracts a numeric sort key (trailing number, or the divisor of an `a : b` ratio). Used ONLY to render rows low→high; never mutates row identity or the batch payload. Unparseable values sort last (stable).
4. **"Ordered by" label** — render the raw `orderingKey`, humanised (underscores → spaces) via a pure `humanizeKey`. The key is the real data; humanising is presentation only.
5. **"Himalaya also found" loose matches** — a NEW pure `splitProposal(proposal)` partitions `proposal.rows` into `members` (rows with a parsed value for the ordering key) and `looseMatches` (rows that lack the key entirely — `value === ""`). Loose matches start EXCLUDED; the existing inline-flag/confirm flow is reserved for members whose value parsed but is uncertain. On a cold corpus (no `orderingKey`) every row is a loose match and the page shows the cold fallback. This is a deliberate, tested behaviour change that fixes the cold→warm transition (S-A).

**Build-gate invariant (unchanged contract):** `canBuild` still requires `orderingKey !== undefined`, ≥1 included non-flagged row, and no included row flagged. The batch payload still maps included rows to `{ sampleId, value }`. Existing `scoping.test.tsx` behavioural assertions (gating, no-silent-loss, D5 exclude) MUST stay green — update only the DOM-query specifics they touch (test-ids preserved where possible).

---

## File structure

- Create `src/lib/scoping/parseSortKey.ts` — pure: value string → numeric sort key | null.
- Create `src/lib/scoping/splitProposal.ts` — pure: OrderingProposal → `{ members, looseMatches }`; also `humanizeKey`.
- Create `src/lib/scoping/dominantPhase.ts` — pure: `IndexEntry[]` → `{ dominant: string | null; coexist: string | null }`.
- Create `src/lib/plot/sparkline.ts` — pure: `(trace, opts)` → SVG path `d` string + viewBox dims (no React).
- Create `src/components/ScopingSparkline.tsx` — presentational SVG sparkline (consumes `sparkline.ts` + `phaseColor`).
- Create `src/components/ScopingRow.tsx` — one member row: grip, sparkline, name/id, confirm-not-fill-out value cell.
- Create `src/components/ScopingValueCell.tsx` — the confirm-not-fill-out value (ink text ⇄ inline edit; flagged = amber).
- Create `src/components/ScopingAutogroupCard.tsx` — the autogroup summary card.
- Create `src/components/ScopingOrderField.tsx` — the "Ordered by" field + note.
- Create `src/components/ScopingLooseMatches.tsx` — "Himalaya also found" section.
- Create `src/components/ScopingPhaseStrip.tsx` — preview phase strip + caption.
- Create `src/components/ScopingFoot.tsx` — narrative gate footer (state line + note + build button).
- Modify `src/queries.ts` — add `useMemberIndices` sibling hook.
- Rewrite `src/pages/SeriesScopingPage.tsx` — compose the worksheet plate; in-session undo.
- Modify `src/components/ScopingConfirmModal.tsx` — Print build-button tokens (S-B/S-C).
- Tests: `test/parseSortKey.test.ts`, `test/splitProposal.test.ts`, `test/dominantPhase.test.ts`, `test/sparkline.test.ts`, `test/ScopingSparkline.test.tsx`, `test/ScopingValueCell.test.tsx`, `test/ScopingRow.test.tsx`, `test/ScopingLooseMatches.test.tsx`, `test/ScopingPhaseStrip.test.tsx`, `test/ScopingFoot.test.tsx`, `test/queries-member-indices.test.tsx`, and updated `test/scoping.test.tsx`.

---

## Task 1: `parseSortKey` — numeric sort key from a value string

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/scoping/parseSortKey.ts`
- Test: `packages/HimalayaUI/frontend/test/parseSortKey.test.ts`

- [ ] **Step 1: Write the failing test**

```ts
import { describe, it, expect } from "vitest";
import { parseSortKey } from "../src/lib/scoping/parseSortKey";

describe("parseSortKey", () => {
  it("parses a plain number", () => {
    expect(parseSortKey("2.5")).toBe(2.5);
    expect(parseSortKey("0")).toBe(0);
  });
  it("parses the divisor of an a : b ratio (low→high by the second term)", () => {
    expect(parseSortKey("1 : 0")).toBe(0);
    expect(parseSortKey("1 : 0.25")).toBe(0.25);
    expect(parseSortKey("1:4")).toBe(4);
  });
  it("parses a trailing number embedded in text", () => {
    expect(parseSortKey("dose 30mM")).toBe(30);
  });
  it("returns null for an unparseable value", () => {
    expect(parseSortKey("DOPC")).toBeNull();
    expect(parseSortKey("")).toBeNull();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd packages/HimalayaUI/frontend && npx vitest run test/parseSortKey.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Write minimal implementation**

```ts
/**
 * Derive a numeric sort key from a parsed ordering value so the worksheet can
 * render members low→high (series-scoping.html: SERIES.sort by `key`). Pure;
 * never mutates row identity or the batch payload. Unparseable → null (those
 * rows sort last, stable).
 *
 * Rules, in order:
 *  - `a : b` ratio → the second term `b` (LL37 : lipid titration sorts by the
 *    lipid side, matching the mockup's `1 : 0 … 1 : 4`).
 *  - otherwise the first standalone number anywhere in the string.
 */
export function parseSortKey(value: string): number | null {
  const ratio = value.match(/^\s*-?\d*\.?\d+\s*:\s*(-?\d*\.?\d+)\s*$/);
  if (ratio) return Number(ratio[1]);
  const num = value.match(/-?\d*\.?\d+/);
  return num ? Number(num[0]) : null;
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run test/parseSortKey.test.ts`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/scoping/parseSortKey.ts packages/HimalayaUI/frontend/test/parseSortKey.test.ts
git commit -m "feat(scoping): parseSortKey — numeric ordering key from a value string (R7 #230)"
```

---

## Task 2: `splitProposal` + `humanizeKey` — members vs loose matches

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/scoping/splitProposal.ts`
- Test: `packages/HimalayaUI/frontend/test/splitProposal.test.ts`

- [ ] **Step 1: Write the failing test**

```ts
import { describe, it, expect } from "vitest";
import { splitProposal, humanizeKey } from "../src/lib/scoping/splitProposal";
import type { OrderingProposal } from "../src/lib/scoping/proposeOrdering";

const proposal = (rows: OrderingProposal["rows"], key: string | undefined): OrderingProposal =>
  ({ orderingKey: key, rows });

describe("splitProposal", () => {
  it("members carry a parsed value; loose matches lack the key (value empty)", () => {
    const p = proposal([
      { sampleId: 10, sampleName: "A", value: "1:1", flagged: false, include: true },
      { sampleId: 11, sampleName: "B", value: "", flagged: true, include: true },
    ], "ratio");
    const { members, looseMatches } = splitProposal(p);
    expect(members.map((r) => r.sampleId)).toEqual([10]);
    expect(looseMatches.map((r) => r.sampleId)).toEqual([11]);
  });
  it("loose matches start excluded (include=false)", () => {
    const p = proposal([
      { sampleId: 11, sampleName: "B", value: "", flagged: true, include: true },
    ], "ratio");
    expect(splitProposal(p).looseMatches[0].include).toBe(false);
  });
  it("a cold corpus (no orderingKey) puts every row in looseMatches", () => {
    const p = proposal([
      { sampleId: 10, sampleName: "A", value: "", flagged: true, include: true },
    ], undefined);
    const { members, looseMatches } = splitProposal(p);
    expect(members).toEqual([]);
    expect(looseMatches).toHaveLength(1);
  });
});

describe("humanizeKey", () => {
  it("turns snake/kebab keys into spaced words", () => {
    expect(humanizeKey("ll37_lipid_ratio")).toBe("ll37 lipid ratio");
    expect(humanizeKey("ratio")).toBe("ratio");
    expect(humanizeKey(undefined)).toBe("—");
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run test/splitProposal.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Write minimal implementation**

```ts
import type { OrderingProposal, OrderingRow } from "./proposeOrdering";

export interface SplitProposal {
  members: OrderingRow[];      // rows that have a parsed value for the ordering key
  looseMatches: OrderingRow[]; // rows that lack the key entirely — surfaced, not assumed
}

/**
 * Partition a proposal into the members the machine WOULD assume into the
 * series (a parsed value for the ordering key) and the loose matches it
 * NOTICED but will not assume (no value for the key). Loose matches start
 * excluded — the human adds them deliberately. Cold corpus (no orderingKey)
 * → everything is a loose match (the page shows its cold fallback). Pure:
 * surfaces proposeOrdering's output, does not change its key selection.
 */
export function splitProposal(p: OrderingProposal): SplitProposal {
  if (p.orderingKey === undefined) {
    return { members: [], looseMatches: p.rows.map((r) => ({ ...r, include: false })) };
  }
  const members: OrderingRow[] = [];
  const looseMatches: OrderingRow[] = [];
  for (const r of p.rows) {
    if (r.value !== "") members.push(r);
    else looseMatches.push({ ...r, include: false });
  }
  return { members, looseMatches };
}

/** Present a tag key as a human label (underscores/hyphens → spaces). */
export function humanizeKey(key: string | undefined): string {
  if (key === undefined || key === "") return "—";
  return key.replace(/[_-]+/g, " ");
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run test/splitProposal.test.ts`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/scoping/splitProposal.ts packages/HimalayaUI/frontend/test/splitProposal.test.ts
git commit -m "feat(scoping): splitProposal + humanizeKey — members vs Himalaya-also-found (R7 #230)"
```

---

## Task 3: `dominantPhase` — phase from indices

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/scoping/dominantPhase.ts`
- Test: `packages/HimalayaUI/frontend/test/dominantPhase.test.ts`

- [ ] **Step 1: Write the failing test**

```ts
import { describe, it, expect } from "vitest";
import { dominantPhase } from "../src/lib/scoping/dominantPhase";
import type { IndexEntry } from "../src/api";

const idx = (phase: string, score: number): IndexEntry => ({
  id: 1, exposure_id: 1, phase, basis: 1, score, r_squared: null, lattice_d: null,
  ngc: null, status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [],
});

describe("dominantPhase", () => {
  it("picks the highest-scored phase as dominant", () => {
    expect(dominantPhase([idx("Lamellar", 0.4), idx("Pn3m", 0.9)]).dominant).toBe("Pn3m");
  });
  it("reports the second distinct phase as coexist", () => {
    const r = dominantPhase([idx("Pn3m", 0.9), idx("Lamellar", 0.6)]);
    expect(r.dominant).toBe("Pn3m");
    expect(r.coexist).toBe("Lamellar");
  });
  it("no coexist when only one phase is present", () => {
    expect(dominantPhase([idx("Pn3m", 0.9)]).coexist).toBeNull();
  });
  it("returns null dominant for an unindexed exposure", () => {
    expect(dominantPhase([]).dominant).toBeNull();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run test/dominantPhase.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Write minimal implementation**

```ts
import type { IndexEntry } from "../../api";

export interface PhaseRead {
  dominant: string | null; // highest-scored phase
  coexist: string | null;  // second distinct phase, if any (two-phase coexistence)
}

/**
 * Read a sample's dominant phase from its candidate indices: highest `score`
 * wins; the next distinct phase (if any) is the coexistence partner that
 * drives the preview strip's two-phase gradient. Pure; null score sorts last.
 */
export function dominantPhase(indices: IndexEntry[]): PhaseRead {
  const ranked = [...indices].sort((a, b) => (b.score ?? -Infinity) - (a.score ?? -Infinity));
  const dominant = ranked[0]?.phase ?? null;
  const coexist = ranked.find((i) => i.phase !== dominant)?.phase ?? null;
  return { dominant, coexist };
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run test/dominantPhase.test.ts`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/scoping/dominantPhase.ts packages/HimalayaUI/frontend/test/dominantPhase.test.ts
git commit -m "feat(scoping): dominantPhase — phase read from candidate indices (R7 #230)"
```

---

## Task 4: `sparkline` — pure SVG path from a trace

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/plot/sparkline.ts`
- Test: `packages/HimalayaUI/frontend/test/sparkline.test.ts`

- [ ] **Step 1: Write the failing test**

```ts
import { describe, it, expect } from "vitest";
import { sparklinePath, SPARK_W, SPARK_H } from "../src/lib/plot/sparkline";

describe("sparklinePath", () => {
  it("returns null for an empty trace", () => {
    expect(sparklinePath({ q: [], I: [], sigma: [] })).toBeNull();
  });
  it("produces a path that starts with M and stays within the viewBox", () => {
    const q = Array.from({ length: 50 }, (_, i) => 0.02 + i * 0.008);
    const I = q.map((x) => Math.exp(-((x - 0.1) ** 2) / 0.0002));
    const p = sparklinePath({ q, I, sigma: q.map(() => 0) });
    expect(p).not.toBeNull();
    expect(p!.startsWith("M")).toBe(true);
    const ys = [...p!.matchAll(/[ML]\s*[\d.]+\s+([\d.]+)/g)].map((m) => Number(m[1]));
    expect(Math.min(...ys)).toBeGreaterThanOrEqual(0);
    expect(Math.max(...ys)).toBeLessThanOrEqual(SPARK_H);
  });
  it("ignores non-finite / non-positive intensities without throwing", () => {
    const p = sparklinePath({ q: [0.02, 0.04, 0.06], I: [NaN, 0, 5], sigma: [0, 0, 0] });
    expect(p).not.toBeNull();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run test/sparkline.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Write minimal implementation**

```ts
import type { Trace } from "../../api";

export const SPARK_W = 76;
export const SPARK_H = 28;
const PAD_X = 4;
const PAD_B = 4;
const AMP = 17;

/**
 * Hand-rolled mini-trace path for the scoping worksheet sparkline
 * (series-scoping.html `sparkline()`): log-x, linear-y, peak-relative scale.
 * Pure — returns just the SVG path `d` (the React wrapper supplies stroke/fill
 * and the phase colour). Returns null when there is nothing to draw.
 */
export function sparklinePath(trace: Trace): string | null {
  const pts: Array<[number, number]> = [];
  for (let i = 0; i < trace.q.length; i++) {
    const q = trace.q[i], I = trace.I[i];
    if (Number.isFinite(q) && q > 0 && Number.isFinite(I)) pts.push([q, Math.max(0, I)]);
  }
  if (pts.length < 2) return null;
  const qMin = pts[0][0], qMax = pts[pts.length - 1][0];
  if (qMax <= qMin) return null;
  const lnMin = Math.log(qMin), lnSpan = Math.log(qMax) - lnMin;
  const xOf = (q: number) => PAD_X + ((Math.log(q) - lnMin) / lnSpan) * (SPARK_W - 2 * PAD_X);
  const peak = Math.max(...pts.map((p) => p[1]), 1e-6);
  const hS = AMP / peak;
  const y0 = SPARK_H - PAD_B;
  let d = "";
  pts.forEach(([q, I], i) => {
    d += (i ? "L" : "M") + xOf(q).toFixed(1) + " " + (y0 - I * hS).toFixed(1) + " ";
  });
  return d.trim();
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run test/sparkline.test.ts`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/plot/sparkline.ts packages/HimalayaUI/frontend/test/sparkline.test.ts
git commit -m "feat(plot): sparklinePath — pure mini-trace SVG path for scoping (R7 #230)"
```

---

## Task 5: `useMemberIndices` — parallel indices fetch sibling

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/queries.ts` (add after `useMemberTraces`, ~line 247)
- Test: `packages/HimalayaUI/frontend/test/queries-member-indices.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { useMemberIndices } from "../src/queries";

function wrap(client = makeClient()) {
  return ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
}

describe("useMemberIndices", () => {
  beforeEach(() => vi.restoreAllMocks());
  it("fetches indices for each exposure id and maps them by id", async () => {
    vi.spyOn(global, "fetch").mockImplementation((input) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url.endsWith("/api/exposures/7/indices"))
        return Promise.resolve(new Response(JSON.stringify([
          { id: 1, exposure_id: 7, phase: "Pn3m", basis: 1, score: 0.9, r_squared: null,
            lattice_d: null, ngc: null, status: "candidate", kind: "auto",
            inputs_hash: null, peaks: [], predicted_q: [] },
        ]), { status: 200, headers: { "Content-Type": "application/json" } }));
      return Promise.resolve(new Response("[]", { status: 200, headers: { "Content-Type": "application/json" } }));
    });
    const { result } = renderHook(() => useMemberIndices([7]), { wrapper: wrap() });
    await waitFor(() => expect(result.current.get(7)?.[0]?.phase).toBe("Pn3m"));
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run test/queries-member-indices.test.tsx`
Expected: FAIL — `useMemberIndices` is not exported.

- [ ] **Step 3: Write minimal implementation** (insert after `useMemberTraces`)

```ts
/**
 * Sibling of `useMemberTraces`: fetch candidate indices for a variable list of
 * exposure ids in parallel for the scoping worksheet's phase reads (sparkline
 * colour + preview strip). Cache keys mirror `useIndices` so single-exposure
 * pages and scoping share one cache row. Returns `Map<exposure_id, IndexEntry[]>`.
 */
export function useMemberIndices(exposureIds: number[]): Map<number, api.IndexEntry[]> {
  return useStableQueryMap(exposureIds, (id) => ({
    queryKey: queryKeys.indices(id),
    queryFn: () => api.listIndices(id),
  }));
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run test/queries-member-indices.test.tsx`
Expected: PASS (1 test).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/queries.ts packages/HimalayaUI/frontend/test/queries-member-indices.test.tsx
git commit -m "feat(queries): useMemberIndices — parallel per-member index fetch (R7 #230)"
```

---

## Task 6: `ScopingSparkline` component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ScopingSparkline.tsx`
- Test: `packages/HimalayaUI/frontend/test/ScopingSparkline.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
import { describe, it, expect } from "vitest";
import { render } from "@testing-library/react";
import { ScopingSparkline } from "../src/components/ScopingSparkline";

const trace = (() => {
  const q = Array.from({ length: 40 }, (_, i) => 0.02 + i * 0.01);
  return { q, I: q.map((x) => Math.exp(-((x - 0.1) ** 2) / 0.0002)), sigma: q.map(() => 0) };
})();

describe("ScopingSparkline", () => {
  it("renders an svg path when a trace is present", () => {
    const { container } = render(<ScopingSparkline trace={trace} phase="Pn3m" />);
    expect(container.querySelector("path")).not.toBeNull();
  });
  it("renders an empty frame (no path) when the trace is undefined", () => {
    const { container, getByTestId } = render(<ScopingSparkline trace={undefined} phase={null} />);
    expect(getByTestId("scoping-sparkline-empty")).toBeInTheDocument();
    expect(container.querySelector("path")).toBeNull();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run test/ScopingSparkline.test.tsx`
Expected: FAIL — module not found.

- [ ] **Step 3: Write minimal implementation**

```tsx
import type { Trace } from "../api";
import { phaseColor } from "../phases";
import { sparklinePath, SPARK_W, SPARK_H } from "../lib/plot/sparkline";

interface Props {
  trace: Trace | undefined;
  phase: string | null;
}

/**
 * The per-row mini-trace on the scoping worksheet (series-scoping.html `.spark`).
 * Stroked + 10%-fill in the row's dominant phase colour; unindexed → neutral.
 * An empty frame stands in while the trace is loading or the sample is
 * unindexed, so rows keep a stable 76×28 rhythm.
 */
export function ScopingSparkline({ trace, phase }: Props): JSX.Element {
  const d = trace ? sparklinePath(trace) : null;
  const color = phase ? phaseColor(phase) : "var(--color-ink-faint)";
  if (!d) {
    return (
      <span
        data-testid="scoping-sparkline-empty"
        className="h-7 w-[76px] shrink-0 overflow-hidden rounded-sm border border-hair bg-paper-sunk"
        aria-hidden
      />
    );
  }
  const y0 = SPARK_H - 4;
  return (
    <span className="h-7 w-[76px] shrink-0 overflow-hidden rounded-sm border border-hair bg-paper-sunk">
      <svg viewBox={`0 0 ${SPARK_W} ${SPARK_H}`} className="block h-full w-full" aria-hidden>
        <line x1={4} y1={y0} x2={SPARK_W - 4} y2={y0} stroke="var(--color-hair)" strokeWidth={1} />
        <path d={`${d} L${SPARK_W - 4} ${y0} L4 ${y0} Z`} fill={color} opacity={0.1} />
        <path d={d} fill="none" stroke={color} strokeWidth={1.4} strokeLinejoin="round" />
      </svg>
    </span>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run test/ScopingSparkline.test.tsx`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ScopingSparkline.tsx packages/HimalayaUI/frontend/test/ScopingSparkline.test.tsx
git commit -m "feat(scoping): ScopingSparkline — per-row mini trace (R7 #230 S-A)"
```

---

## Task 7: `ScopingValueCell` — confirm-not-fill-out value (S-E)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ScopingValueCell.tsx`
- Test: `packages/HimalayaUI/frontend/test/ScopingValueCell.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ScopingValueCell } from "../src/components/ScopingValueCell";

describe("ScopingValueCell", () => {
  it("a confident value renders as ink text, not an input", () => {
    render(<ScopingValueCell sampleId={10} sampleName="A" value="1:1" flagged={false}
      onChange={() => {}} onToggleFlag={() => {}} />);
    expect(screen.getByTestId("scoping-value-10")).toHaveTextContent("1:1");
    expect(screen.queryByRole("textbox")).toBeNull();
  });
  it("a flagged value shows the amber 'check the read' affordance", () => {
    render(<ScopingValueCell sampleId={11} sampleName="B" value="1 : 0" flagged
      onChange={() => {}} onToggleFlag={() => {}} />);
    const cell = screen.getByTestId("scoping-value-11");
    expect(cell).toHaveAttribute("data-flagged", "true");
    expect(cell).toHaveTextContent(/check the read/i);
  });
  it("clicking a confident value re-opens it as an input", () => {
    render(<ScopingValueCell sampleId={10} sampleName="A" value="1:1" flagged={false}
      onChange={() => {}} onToggleFlag={() => {}} />);
    fireEvent.click(screen.getByTestId("scoping-value-10"));
    expect(screen.getByTestId("scoping-value-input-10")).toHaveValue("1:1");
  });
  it("committing an edit calls onChange with the new value", () => {
    const onChange = vi.fn();
    render(<ScopingValueCell sampleId={10} sampleName="A" value="1:1" flagged={false}
      onChange={onChange} onToggleFlag={() => {}} />);
    fireEvent.click(screen.getByTestId("scoping-value-10"));
    const input = screen.getByTestId("scoping-value-input-10");
    fireEvent.change(input, { target: { value: "2:1" } });
    fireEvent.blur(input);
    expect(onChange).toHaveBeenCalledWith("2:1");
  });
  it("clicking a flagged value accepts the read (onToggleFlag)", () => {
    const onToggleFlag = vi.fn();
    render(<ScopingValueCell sampleId={11} sampleName="B" value="1 : 0" flagged
      onChange={() => {}} onToggleFlag={onToggleFlag} />);
    fireEvent.click(screen.getByTestId("scoping-value-11"));
    expect(onToggleFlag).toHaveBeenCalledTimes(1);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run test/ScopingValueCell.test.tsx`
Expected: FAIL — module not found.

- [ ] **Step 3: Write minimal implementation**

```tsx
import { useState } from "react";

interface Props {
  sampleId: number;
  sampleName: string;
  value: string;
  flagged: boolean;
  onChange: (value: string) => void;     // commit a new value (also clears the flag)
  onToggleFlag: () => void;              // accept the read / re-flag (undo-tracked by parent)
}

/**
 * Confirm-not-fill-out value cell (series-scoping.html `.s-val`, finding S-E).
 *  - Confident: ink text; click re-opens it as an input (every value stays
 *    re-openable, including one resolved by mistake).
 *  - Flagged: amber dashed underline + "▸ check the read"; click accepts the
 *    read (clears the flag) — never a permanent text input.
 *  - Editing: a transient input; blur/Enter commits via onChange, Escape cancels.
 */
export function ScopingValueCell({
  sampleId, sampleName, value, flagged, onChange, onToggleFlag,
}: Props): JSX.Element {
  const [editing, setEditing] = useState(false);
  const [draft, setDraft] = useState(value);

  function open(): void {
    if (flagged) { onToggleFlag(); return; } // accept the read in one click
    setDraft(value);
    setEditing(true);
  }
  function commit(): void {
    setEditing(false);
    if (draft !== value) onChange(draft);
  }

  if (editing) {
    return (
      <input
        type="text"
        data-testid={`scoping-value-input-${sampleId}`}
        aria-label={`Value for ${sampleName}`}
        value={draft}
        autoFocus
        onChange={(e) => setDraft(e.target.value)}
        onBlur={commit}
        onKeyDown={(e) => {
          if (e.key === "Enter") commit();
          if (e.key === "Escape") { setDraft(value); setEditing(false); }
        }}
        className="w-24 shrink-0 rounded border border-hair-strong bg-plate px-2 py-1 text-right font-mono text-[13px] text-ink"
      />
    );
  }

  return (
    <button
      type="button"
      data-testid={`scoping-value-${sampleId}`}
      data-flagged={flagged ? "true" : undefined}
      onClick={open}
      title={flagged ? "Click to accept this read" : "Click to re-open this value"}
      className={`group shrink-0 text-right ${flagged ? "text-print-accent" : "text-ink"}`}
    >
      <span
        className={`font-mono text-[13px] font-bold ${
          flagged
            ? "border-b-[1.5px] border-dashed border-print-accent/60 pb-px"
            : "border-b-[1.5px] border-dotted border-transparent pb-px group-hover:border-hair-strong"
        }`}
      >
        {value || "—"}
      </span>
      {flagged ? (
        <span className="mt-0.5 block text-[9px] font-bold uppercase tracking-wide text-print-accent">
          ▸ check the read
        </span>
      ) : null}
    </button>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run test/ScopingValueCell.test.tsx`
Expected: PASS (5 tests).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ScopingValueCell.tsx packages/HimalayaUI/frontend/test/ScopingValueCell.test.tsx
git commit -m "feat(scoping): ScopingValueCell — confirm-not-fill-out value (R7 #230 S-E)"
```

---

## Task 8: `ScopingRow` — one member row (grip + sparkline + name + value)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ScopingRow.tsx`
- Test: `packages/HimalayaUI/frontend/test/ScopingRow.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ScopingRow } from "../src/components/ScopingRow";

const baseRow = { sampleId: 10, sampleName: "Lipid 1-2 + LL37 1:1", value: "1:1", flagged: false, include: true };

describe("ScopingRow", () => {
  it("renders the sample name, a grip and a sparkline frame", () => {
    render(<ScopingRow row={baseRow} trace={undefined} phase={null}
      onChangeValue={() => {}} onToggleFlag={() => {}} />);
    const row = screen.getByTestId("scoping-row-10");
    expect(row).toHaveTextContent("Lipid 1-2 + LL37 1:1");
    expect(screen.getByTestId("scoping-grip-10")).toBeInTheDocument();
    expect(screen.getByTestId("scoping-sparkline-empty")).toBeInTheDocument();
  });
  it("marks the row data-flagged when flagged", () => {
    render(<ScopingRow row={{ ...baseRow, flagged: true }} trace={undefined} phase={null}
      onChangeValue={() => {}} onToggleFlag={() => {}} />);
    expect(screen.getByTestId("scoping-row-10")).toHaveAttribute("data-flagged", "true");
  });
  it("threads value edits up via onChangeValue", () => {
    const onChangeValue = vi.fn();
    render(<ScopingRow row={baseRow} trace={undefined} phase={null}
      onChangeValue={onChangeValue} onToggleFlag={() => {}} />);
    fireEvent.click(screen.getByTestId("scoping-value-10"));
    const input = screen.getByTestId("scoping-value-input-10");
    fireEvent.change(input, { target: { value: "3:1" } });
    fireEvent.blur(input);
    expect(onChangeValue).toHaveBeenCalledWith("3:1");
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run test/ScopingRow.test.tsx`
Expected: FAIL — module not found.

- [ ] **Step 3: Write minimal implementation**

```tsx
import type { Trace } from "../api";
import type { OrderingRow } from "../lib/scoping/proposeOrdering";
import { ScopingSparkline } from "./ScopingSparkline";
import { ScopingValueCell } from "./ScopingValueCell";

interface Props {
  row: OrderingRow;
  trace: Trace | undefined;
  phase: string | null;
  onChangeValue: (value: string) => void;
  onToggleFlag: () => void;
}

/**
 * One member row of the scoping worksheet (series-scoping.html `.srow`):
 * a drag grip, the trace sparkline, the sample name + id, and the
 * confirm-not-fill-out value cell. Flagged rows get an amber row tint.
 */
export function ScopingRow({ row, trace, phase, onChangeValue, onToggleFlag }: Props): JSX.Element {
  return (
    <div
      data-testid={`scoping-row-${row.sampleId}`}
      data-flagged={row.flagged ? "true" : undefined}
      className={`flex items-center gap-3 border-b border-hair px-2 py-2 last:border-b-0 ${
        row.flagged ? "bg-print-accent/5" : ""
      }`}
    >
      <span
        data-testid={`scoping-grip-${row.sampleId}`}
        aria-hidden
        className="shrink-0 cursor-grab select-none leading-none tracking-tighter text-hair-strong"
      >
        ⠿
      </span>
      <ScopingSparkline trace={trace} phase={phase} />
      <span className="min-w-0 flex-1">
        <span className="block truncate text-[13px] font-semibold text-ink">{row.sampleName}</span>
        <span className="block font-mono text-[10.5px] text-ink-faint">smp_{row.sampleId}</span>
      </span>
      <ScopingValueCell
        sampleId={row.sampleId}
        sampleName={row.sampleName}
        value={row.value}
        flagged={row.flagged}
        onChange={onChangeValue}
        onToggleFlag={onToggleFlag}
      />
    </div>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run test/ScopingRow.test.tsx`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ScopingRow.tsx packages/HimalayaUI/frontend/test/ScopingRow.test.tsx
git commit -m "feat(scoping): ScopingRow — worksheet member row (R7 #230 S-A/S-F)"
```

---

## Task 9: `ScopingAutogroupCard` + `ScopingOrderField`

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ScopingAutogroupCard.tsx`
- Create: `packages/HimalayaUI/frontend/src/components/ScopingOrderField.tsx`
- Test: `packages/HimalayaUI/frontend/test/ScopingAutogroup.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { ScopingAutogroupCard } from "../src/components/ScopingAutogroupCard";
import { ScopingOrderField } from "../src/components/ScopingOrderField";

describe("ScopingAutogroupCard", () => {
  it("summarises the grouping: member count, key label, flag count", () => {
    render(<ScopingAutogroupCard memberCount={6} keyLabel="ratio" flagCount={1} />);
    const card = screen.getByTestId("scoping-autogroup");
    expect(card).toHaveTextContent("6 samples");
    expect(card).toHaveTextContent("ratio");
    expect(card).toHaveTextContent(/one needs a look/i);
  });
  it("reads cleanly when nothing is flagged", () => {
    render(<ScopingAutogroupCard memberCount={5} keyLabel="ratio" flagCount={0} />);
    expect(screen.getByTestId("scoping-autogroup")).toHaveTextContent(/parsed cleanly/i);
  });
});

describe("ScopingOrderField", () => {
  it("shows the humanised ordering key as the field value", () => {
    render(<ScopingOrderField keyLabel="ll37 lipid ratio" />);
    expect(screen.getByTestId("scoping-order-field")).toHaveTextContent("ll37 lipid ratio");
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run test/ScopingAutogroup.test.tsx`
Expected: FAIL — modules not found.

- [ ] **Step 3: Write minimal implementations**

`ScopingAutogroupCard.tsx`:

```tsx
interface Props {
  memberCount: number;
  keyLabel: string;
  flagCount: number;
}

/**
 * The machine's grouping summary (series-scoping.html `.autogroup`): one
 * confident breath naming what Himalaya did — grouped N samples, read the
 * order from <key>, M parsed cleanly, K need a look.
 */
export function ScopingAutogroupCard({ memberCount, keyLabel, flagCount }: Props): JSX.Element {
  const clean = memberCount - flagCount;
  const flagPhrase =
    flagCount === 0
      ? "all parsed cleanly"
      : flagCount === 1
        ? `${clean} parsed cleanly, one needs a look`
        : `${clean} parsed cleanly, ${flagCount} need a look`;
  return (
    <div
      data-testid="scoping-autogroup"
      className="mt-4 flex items-start gap-2.5 rounded-md border border-hair bg-paper-sunk px-3 py-3"
    >
      <svg viewBox="0 0 16 16" className="mt-px h-[15px] w-[15px] shrink-0" aria-hidden>
        <path
          d="M8 1.4l1.7 4.2 4.5.3-3.5 2.9 1.2 4.4L8 10.9 4.1 13.2l1.2-4.4L1.8 5.9l4.5-.3z"
          fill="var(--color-print-accent)"
        />
      </svg>
      <p className="text-xs leading-relaxed text-ink-soft">
        You selected <b className="font-semibold text-ink">{memberCount} samples</b>. Himalaya
        grouped them from their names and read the order from{" "}
        <b className="font-semibold text-ink">{keyLabel}</b> — {flagPhrase}.
      </p>
    </div>
  );
}
```

`ScopingOrderField.tsx`:

```tsx
interface Props {
  keyLabel: string;
}

/**
 * The "Ordered by" field (series-scoping.html `.order-field`): the one real
 * decision on this surface. Read-only in this milestone — re-selecting the
 * ordering variable is a follow-up; the field still presents as a decision.
 */
export function ScopingOrderField({ keyLabel }: Props): JSX.Element {
  return (
    <>
      <div className="mb-1.5 mt-5 text-[10.5px] font-bold uppercase tracking-wider text-ink-faint">
        Ordered by
      </div>
      <div
        data-testid="scoping-order-field"
        className="flex items-center justify-between rounded-md border border-hair-strong bg-plate px-3.5 py-3"
      >
        <span className="text-[15px] font-semibold text-ink">{keyLabel}</span>
        <span className="text-xs text-ink-faint" aria-hidden>▾</span>
      </div>
      <p className="mt-1.5 text-[11px] text-ink-faint">
        Read from the sample names.
      </p>
    </>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run test/ScopingAutogroup.test.tsx`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ScopingAutogroupCard.tsx packages/HimalayaUI/frontend/src/components/ScopingOrderField.tsx packages/HimalayaUI/frontend/test/ScopingAutogroup.test.tsx
git commit -m "feat(scoping): autogroup summary + Ordered-by field (R7 #230 S-A/S-G)"
```

---

## Task 10: `ScopingLooseMatches` — "Himalaya also found"

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ScopingLooseMatches.tsx`
- Test: `packages/HimalayaUI/frontend/test/ScopingLooseMatches.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ScopingLooseMatches } from "../src/components/ScopingLooseMatches";

const loose = { sampleId: 21, sampleName: "Lipid 1-1 + LL37 1:1", value: "", flagged: true, include: false };

describe("ScopingLooseMatches", () => {
  it("lists each loose match with an Add affordance", () => {
    render(<ScopingLooseMatches rows={[loose]} traces={new Map()} phases={new Map()}
      onAdd={() => {}} />);
    expect(screen.getByTestId("scoping-loose-21")).toHaveTextContent("Lipid 1-1 + LL37 1:1");
    expect(screen.getByTestId("scoping-loose-add-21")).toBeInTheDocument();
  });
  it("clicking Add threads the sample id up", () => {
    const onAdd = vi.fn();
    render(<ScopingLooseMatches rows={[loose]} traces={new Map()} phases={new Map()} onAdd={onAdd} />);
    fireEvent.click(screen.getByTestId("scoping-loose-add-21"));
    expect(onAdd).toHaveBeenCalledWith(21);
  });
  it("shows the empty note when nothing else matches", () => {
    render(<ScopingLooseMatches rows={[]} traces={new Map()} phases={new Map()} onAdd={() => {}} />);
    expect(screen.getByTestId("scoping-loose-empty")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run test/ScopingLooseMatches.test.tsx`
Expected: FAIL — module not found.

- [ ] **Step 3: Write minimal implementation**

```tsx
import type { Trace } from "../api";
import type { OrderingRow } from "../lib/scoping/proposeOrdering";
import { ScopingSparkline } from "./ScopingSparkline";

interface Props {
  rows: OrderingRow[];
  traces: Map<number, Trace>;
  phases: Map<number, string | null>;
  onAdd: (sampleId: number) => void; // fold a loose match into the series
}

/**
 * "Himalaya also found" (series-scoping.html `.candidates`): samples the machine
 * noticed but would not assume into the series — here, corpus samples that lack
 * a value for the ordering variable. Each is addable in one click; coherence is
 * the human's call.
 */
export function ScopingLooseMatches({ rows, traces, phases, onAdd }: Props): JSX.Element {
  return (
    <div className="mt-5 border-t border-hair pt-4">
      <div className="mb-2 text-[10.5px] font-bold uppercase tracking-wider text-ink-faint">
        Himalaya also found
      </div>
      {rows.length === 0 ? (
        <div data-testid="scoping-loose-empty" className="text-[11.5px] italic text-ink-faint">
          Nothing else in the corpus matches this grouping.
        </div>
      ) : (
        <div className="flex flex-col gap-2">
          {rows.map((r) => (
            <div
              key={r.sampleId}
              data-testid={`scoping-loose-${r.sampleId}`}
              className="flex items-center gap-3 rounded-md border border-dashed border-hair-strong px-2 py-2"
            >
              <span className="opacity-70">
                <ScopingSparkline trace={traces.get(r.sampleId)} phase={phases.get(r.sampleId) ?? null} />
              </span>
              <div className="min-w-0 flex-1">
                <div className="truncate text-[12.5px] font-semibold text-ink-soft">{r.sampleName}</div>
                <div className="text-[11px] text-ink-faint">
                  No <span className="font-semibold text-print-accent">value</span> for the ordering variable — add it to include.
                </div>
              </div>
              <button
                type="button"
                data-testid={`scoping-loose-add-${r.sampleId}`}
                onClick={() => onAdd(r.sampleId)}
                className="shrink-0 rounded-md border border-hair-strong bg-plate px-2.5 py-1.5 text-[11.5px] font-semibold text-ink hover:bg-paper-sunk"
              >
                + Add to series
              </button>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run test/ScopingLooseMatches.test.tsx`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ScopingLooseMatches.tsx packages/HimalayaUI/frontend/test/ScopingLooseMatches.test.tsx
git commit -m "feat(scoping): ScopingLooseMatches — Himalaya also found (R7 #230 S-A)"
```

---

## Task 11: `ScopingPhaseStrip` — preview phase strip

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ScopingPhaseStrip.tsx`
- Test: `packages/HimalayaUI/frontend/test/ScopingPhaseStrip.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { ScopingPhaseStrip } from "../src/components/ScopingPhaseStrip";

describe("ScopingPhaseStrip", () => {
  it("renders one segment per member", () => {
    render(<ScopingPhaseStrip reads={[
      { dominant: "Pn3m", coexist: null },
      { dominant: "Lamellar", coexist: null },
    ]} />);
    expect(screen.getAllByTestId("scoping-ps-seg")).toHaveLength(2);
  });
  it("captions a transition first → last", () => {
    render(<ScopingPhaseStrip reads={[
      { dominant: "Pn3m", coexist: null },
      { dominant: "Lamellar", coexist: null },
    ]} />);
    const cap = screen.getByTestId("scoping-ps-cap");
    expect(cap).toHaveTextContent("Pn3m");
    expect(cap).toHaveTextContent("Lamellar");
  });
  it("captions 'throughout' when first and last agree", () => {
    render(<ScopingPhaseStrip reads={[
      { dominant: "Pn3m", coexist: null },
      { dominant: "Pn3m", coexist: null },
    ]} />);
    expect(screen.getByTestId("scoping-ps-cap")).toHaveTextContent(/throughout/i);
  });
  it("renders nothing meaningful when no members are indexed", () => {
    render(<ScopingPhaseStrip reads={[{ dominant: null, coexist: null }]} />);
    expect(screen.getByTestId("scoping-ps-cap")).toHaveTextContent(/not yet indexed/i);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run test/ScopingPhaseStrip.test.tsx`
Expected: FAIL — module not found.

- [ ] **Step 3: Write minimal implementation**

```tsx
import { phaseColor } from "../phases";
import type { PhaseRead } from "../lib/scoping/dominantPhase";

interface Props {
  reads: PhaseRead[]; // one per member, in display (low→high) order
}

function segBackground(r: PhaseRead): string {
  if (r.dominant === null) return "var(--color-hair)";
  if (r.coexist) {
    return `linear-gradient(100deg, ${phaseColor(r.dominant)} 42%, ${phaseColor(r.coexist)} 58%)`;
  }
  return phaseColor(r.dominant);
}

/**
 * Preview phase strip (series-scoping.html `.preview`): the phase story this
 * order will print — one cell per member (two-phase coexistence = gradient),
 * captioned first→last (or "throughout" when they agree).
 */
export function ScopingPhaseStrip({ reads }: Props): JSX.Element {
  const first = reads.find((r) => r.dominant !== null)?.dominant ?? null;
  const last = [...reads].reverse().find((r) => r.dominant !== null)?.dominant ?? null;
  return (
    <div className="mt-5">
      <div className="mb-1.5 text-[10.5px] font-bold uppercase tracking-wider text-ink-faint">
        Preview — phase across the series
      </div>
      <div className="flex h-2 gap-0.5">
        {reads.map((r, i) => (
          <div
            key={i}
            data-testid="scoping-ps-seg"
            className="flex-1 rounded-[1.5px]"
            style={{ background: segBackground(r) }}
          />
        ))}
      </div>
      <div data-testid="scoping-ps-cap" className="mt-2 flex items-center gap-1.5 text-[11.5px] text-ink-soft">
        {first === null ? (
          <span className="text-ink-faint">Members not yet indexed — phase preview unavailable.</span>
        ) : first === last ? (
          <>
            <span className="font-semibold" style={{ color: phaseColor(first) }}>{first}</span>
            <span className="text-ink-faint">throughout</span>
          </>
        ) : (
          <>
            <span className="font-semibold" style={{ color: phaseColor(first) }}>{first}</span>
            <span className="text-ink-faint">→</span>
            <span className="font-semibold" style={{ color: phaseColor(last!) }}>{last}</span>
          </>
        )}
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run test/ScopingPhaseStrip.test.tsx`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ScopingPhaseStrip.tsx packages/HimalayaUI/frontend/test/ScopingPhaseStrip.test.tsx
git commit -m "feat(scoping): ScopingPhaseStrip — preview phase strip (R7 #230 S-A)"
```

---

## Task 12: `ScopingFoot` — narrative gate footer (S-A/S-B/S-C)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ScopingFoot.tsx`
- Test: `packages/HimalayaUI/frontend/test/ScopingFoot.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ScopingFoot } from "../src/components/ScopingFoot";

describe("ScopingFoot", () => {
  it("warns and disables the build when values are flagged", () => {
    render(<ScopingFoot flagCount={1} memberCount={6} keyLabel="ratio" canBuild={false} onBuild={() => {}} />);
    expect(screen.getByTestId("scoping-foot-state")).toHaveTextContent(/to check before you can build/i);
    expect(screen.getByTestId("scoping-open-confirm")).toBeDisabled();
  });
  it("reads ready and enables the build when clear", () => {
    render(<ScopingFoot flagCount={0} memberCount={6} keyLabel="ratio" canBuild onBuild={() => {}} />);
    expect(screen.getByTestId("scoping-foot-state")).toHaveTextContent(/ready to build/i);
    expect(screen.getByTestId("scoping-open-confirm")).not.toBeDisabled();
  });
  it("fires onBuild when the enabled button is clicked", () => {
    const onBuild = vi.fn();
    render(<ScopingFoot flagCount={0} memberCount={1} keyLabel="ratio" canBuild onBuild={onBuild} />);
    fireEvent.click(screen.getByTestId("scoping-open-confirm"));
    expect(onBuild).toHaveBeenCalledTimes(1);
  });
  it("the build button carries Print ink tokens, not ice-blue accent", () => {
    render(<ScopingFoot flagCount={0} memberCount={1} keyLabel="ratio" canBuild onBuild={() => {}} />);
    const btn = screen.getByTestId("scoping-open-confirm");
    expect(btn.className).toContain("bg-ink");
    expect(btn.className).toContain("text-paper");
    expect(btn.className).not.toContain("bg-accent");
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run test/ScopingFoot.test.tsx`
Expected: FAIL — module not found.

> Note: asserting on `bg-ink`/`text-paper` class strings is the explicit S-B/S-C
> remediation check (the finding is literally "wrong token"); this is the one
> sanctioned exception to the no-Tailwind-in-tests rule, scoped to the token fix.

- [ ] **Step 3: Write minimal implementation**

```tsx
interface Props {
  flagCount: number;
  memberCount: number;
  keyLabel: string;
  canBuild: boolean;
  onBuild: () => void;
}

/**
 * The narrative gate footer (series-scoping.html `.scope-foot`): a confirmation
 * state line (amber "N to check" / sage "ready"), the metadata-as-byproduct
 * note, and the single durable action. Build button = Print ink (S-B), greyed
 * ink when gated (S-C) — never the ice-blue accent.
 */
export function ScopingFoot({ flagCount, memberCount, keyLabel, canBuild, onBuild }: Props): JSX.Element {
  const ready = flagCount === 0;
  const stateText = ready
    ? `All ${memberCount} values confirmed — ready to build`
    : `${flagCount} value${flagCount === 1 ? "" : "s"} to check before you can build`;
  return (
    <div className="mt-6 flex items-center justify-between gap-5 border-t border-hair pt-4">
      <div className="flex flex-col gap-1">
        <div
          data-testid="scoping-foot-state"
          className={`flex items-center gap-2 text-[12.5px] font-semibold ${ready ? "text-ink" : "text-print-accent"}`}
        >
          <span
            className="h-2 w-2 shrink-0 rounded-full"
            style={{ background: ready ? phaseReady : "var(--color-print-accent)" }}
          />
          {stateText}
        </div>
        <div className="max-w-[42ch] text-[10.5px] text-ink-faint">
          Confirming records the {keyLabel} on every sample — the next series that needs it already knows.
        </div>
      </div>
      <button
        type="button"
        data-testid="scoping-open-confirm"
        disabled={!canBuild}
        onClick={onBuild}
        title={canBuild ? undefined : "Check the flagged values above before building"}
        className="shrink-0 rounded-md border border-ink bg-ink px-4.5 py-2.5 text-[13px] font-semibold text-paper transition-colors hover:bg-ink/85 disabled:cursor-not-allowed disabled:border-hair-strong disabled:bg-paper-sunk disabled:text-ink-faint"
      >
        Confirm &amp; build →
      </button>
    </div>
  );
}

const phaseReady = "oklch(0.520 0.120 162)"; // sage — the ready dot (matches Im3m hue)
```

> If `px-4.5` is not in the Tailwind scale, substitute `px-[18px] py-[11px]`.

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run test/ScopingFoot.test.tsx`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ScopingFoot.tsx packages/HimalayaUI/frontend/test/ScopingFoot.test.tsx
git commit -m "feat(scoping): ScopingFoot — narrative gate footer w/ Print build button (R7 #230 S-A/S-B/S-C)"
```

---

## Task 13: Fix `ScopingConfirmModal` build button tokens (S-B/S-C)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ScopingConfirmModal.tsx:51-55`
- Test: append to `packages/HimalayaUI/frontend/test/scoping.test.tsx` (ScopingConfirmModal describe block)

- [ ] **Step 1: Write the failing test** (add inside the existing `describe("ScopingConfirmModal", …)`)

```tsx
it("the confirm-build button uses Print ink tokens, not the accent (S-B)", () => {
  render(<ScopingConfirmModal open orderingKey="ratio" count={2}
    onConfirm={() => {}} onClose={() => {}} />);
  const btn = screen.getByTestId("scoping-confirm-build");
  expect(btn.className).toContain("bg-ink");
  expect(btn.className).toContain("text-paper");
  expect(btn.className).not.toContain("bg-accent");
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run test/scoping.test.tsx -t "Print ink tokens"`
Expected: FAIL — current class is `bg-accent`.

- [ ] **Step 3: Edit the button class**

Change `ScopingConfirmModal.tsx` line 53 from:

```tsx
className="rounded bg-accent px-3 py-1.5 text-sm font-semibold text-paper">
```

to:

```tsx
className="rounded border border-ink bg-ink px-3 py-1.5 text-sm font-semibold text-paper hover:bg-ink/85">
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run test/scoping.test.tsx -t "Print ink tokens"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ScopingConfirmModal.tsx packages/HimalayaUI/frontend/test/scoping.test.tsx
git commit -m "fix(scoping): confirm-modal build button inherits Print ink (R7 #230 S-B)"
```

---

## Task 14: Rewrite `SeriesScopingPage` — compose the worksheet plate

This is the integration task. It threads the pure modules + per-member fetches + in-session undo into the worksheet. The existing mutation surface (`useScopeSeries`, `pendingBuildRef`, navigate-on-success, error banner) is preserved verbatim.

**Files:**
- Rewrite: `packages/HimalayaUI/frontend/src/pages/SeriesScopingPage.tsx`
- Test: update `packages/HimalayaUI/frontend/test/scoping.test.tsx` (SeriesScopingPage block) + add worksheet + undo specs.

- [ ] **Step 1: Update the existing SeriesScopingPage specs to the worksheet DOM**

In `test/scoping.test.tsx`, the `SeriesScopingPage` describe block: the proposal test currently asserts `scoping-value-10` is an `<input>` with `.toHaveValue("1:1")`. Under confirm-not-fill-out it is ink text. Replace that assertion:

Change:
```tsx
expect(screen.getByTestId("scoping-value-10")).toHaveValue("1:1");
```
to:
```tsx
expect(screen.getByTestId("scoping-value-10")).toHaveTextContent("1:1");
```

The gating / no-loss / D5 tests use `scoping-open-confirm`, `scoping-include-*`, `scoping-confirm-build`, `scoping-error-banner`, `folio-stub` — preserve those test-ids in the rewrite so they keep passing. NOTE: the D5 test unchecks `scoping-include-11`; in the worksheet, member include is the loose-match Add/remove model, so update D5 — see Step 2.

- [ ] **Step 2: Add worksheet + undo + loose-match specs** (append to the `SeriesScopingPage` describe)

```tsx
it("renders the worksheet plate with autogroup summary + serif title + sparklines", async () => {
  vi.spyOn(global, "fetch").mockImplementation((input) => {
    const url = typeof input === "string" ? input : (input as Request).url;
    if (url.endsWith("/api/sample-tags")) return Promise.resolve(jsonRes([{ key: "ratio", value: "1:1" }]));
    if (url.endsWith("/api/picker-samples")) return Promise.resolve(jsonRes([pickerRow(10, "A", "1:1"), pickerRow(11, "B", "2:1")]));
    if (url.includes("/indices")) return Promise.resolve(jsonRes([]));
    if (url.includes("/trace")) return Promise.resolve(jsonRes({ q: [0.02, 0.1, 0.4], I: [1, 5, 1], sigma: [0, 0, 0] }));
    return Promise.resolve(jsonRes([]));
  });
  renderScoping();
  await waitFor(() => expect(screen.getByTestId("scoping-autogroup")).toBeInTheDocument());
  expect(screen.getByTestId("scoping-plate")).toBeInTheDocument();
  expect(screen.getByTestId("scoping-title")).toBeInTheDocument();
  await waitFor(() => expect(screen.getByTestId("scoping-row-10")).toBeInTheDocument());
});

it("a sample missing the ordering key falls into Himalaya-also-found (excluded), not the member list", async () => {
  vi.spyOn(global, "fetch").mockImplementation((input) => {
    const url = typeof input === "string" ? input : (input as Request).url;
    if (url.endsWith("/api/sample-tags")) return Promise.resolve(jsonRes([{ key: "ratio", value: "1:1" }]));
    if (url.endsWith("/api/picker-samples")) return Promise.resolve(jsonRes([pickerRow(10, "A", "1:1"), pickerRow(11, "B")]));
    if (url.includes("/indices")) return Promise.resolve(jsonRes([]));
    if (url.includes("/trace")) return Promise.resolve(jsonRes({ q: [0.02, 0.4], I: [1, 2], sigma: [0, 0] }));
    return Promise.resolve(jsonRes([]));
  });
  renderScoping();
  await waitFor(() => expect(screen.getByTestId("scoping-row-10")).toBeInTheDocument());
  // sample 11 lacks the key → loose match, not a member row
  expect(screen.queryByTestId("scoping-row-11")).toBeNull();
  expect(screen.getByTestId("scoping-loose-11")).toBeInTheDocument();
  // build is allowed: the one member (10) is non-flagged
  await waitFor(() => expect(screen.getByTestId("scoping-open-confirm")).not.toBeDisabled());
});

it("re-opening a confident value then Undo restores it (in-session ⌘Z)", async () => {
  vi.spyOn(global, "fetch").mockImplementation((input) => {
    const url = typeof input === "string" ? input : (input as Request).url;
    if (url.endsWith("/api/sample-tags")) return Promise.resolve(jsonRes([{ key: "ratio", value: "1:1" }]));
    if (url.endsWith("/api/picker-samples")) return Promise.resolve(jsonRes([pickerRow(10, "A", "1:1")]));
    if (url.includes("/indices")) return Promise.resolve(jsonRes([]));
    if (url.includes("/trace")) return Promise.resolve(jsonRes({ q: [0.02, 0.4], I: [1, 2], sigma: [0, 0] }));
    return Promise.resolve(jsonRes([]));
  });
  renderScoping();
  // confident value (10) → click flags it open as needing a look (toggle)
  const cell = await screen.findByTestId("scoping-value-10");
  fireEvent.click(cell);
  await waitFor(() => expect(screen.getByTestId("scoping-row-10")).toHaveAttribute("data-flagged", "true"));
  // undo restores the confirmed state
  fireEvent.click(screen.getByTestId("scoping-undo"));
  await waitFor(() => expect(screen.getByTestId("scoping-row-10")).not.toHaveAttribute("data-flagged"));
});

it("discard navigates back to the folio", async () => {
  vi.spyOn(global, "fetch").mockImplementation(() =>
    Promise.resolve(jsonRes([])));
  renderScoping();
  fireEvent.click(await screen.findByTestId("scoping-discard"));
  await waitFor(() => expect(screen.getByTestId("folio-stub")).toBeInTheDocument());
});
```

Update the existing **D5** test: instead of unchecking `scoping-include-11`, it now asserts that sample 11 (which lacks the ordering key in the existing fixture? No — it has `2:1`). Keep the D5 contract by Adding the loose match and then removing it, OR: since both fixture samples have values, the D5 "exclude a member" affordance becomes "a member is removed". Add a per-member remove control: in the rewrite each member row's value cell stays, and exclusion is via Undo of an Add — but members are auto-included. To preserve the D5 batch-exclusion contract with minimal surface, KEEP `scoping-include-{id}` as a small per-row checkbox is NOT in the mockup. Instead, change D5 to: give sample 11 NO ordering-key value so it is a loose match (excluded by default), then confirm only sample 10 is in the batch:

```tsx
it("D5: a sample without the ordering value is omitted from the batch (loose, excluded)", async () => {
  const posted: any[] = [];
  vi.spyOn(global, "fetch").mockImplementation((input, init) => {
    const url = typeof input === "string" ? input : (input as Request).url;
    if (url.endsWith("/api/sample-tags")) return Promise.resolve(jsonRes([{ key: "ratio", value: "1:1" }]));
    if (url.endsWith("/api/picker-samples"))
      return Promise.resolve(jsonRes([pickerRow(10, "A", "1:1"), pickerRow(11, "B")])); // 11: no ratio
    if (url.includes("/indices")) return Promise.resolve(jsonRes([]));
    if (url.includes("/trace")) return Promise.resolve(jsonRes({ q: [0.02, 0.4], I: [1, 2], sigma: [0, 0] }));
    if (url.endsWith("/api/samples/tags/batch")) {
      posted.push(init?.body ? JSON.parse(String(init.body)) : null);
      return Promise.resolve(jsonRes([], 201));
    }
    return Promise.resolve(jsonRes([]));
  });
  renderScoping();
  await waitFor(() => expect(screen.getByTestId("scoping-open-confirm")).not.toBeDisabled());
  fireEvent.click(screen.getByTestId("scoping-open-confirm"));
  fireEvent.click(await screen.findByTestId("scoping-confirm-build"));
  await waitFor(() => expect(posted[0]).toMatchObject({ key: "ratio", source: "scoping", tags: [{ sample_id: 10, value: "1:1" }] }));
  expect((posted[0].tags as { sample_id: number }[]).some((t) => t.sample_id === 11)).toBe(false);
});
```

Delete the now-obsolete original D5 test (the `scoping-include-11` checkbox one).

- [ ] **Step 3: Run the updated specs to verify they fail**

Run: `npx vitest run test/scoping.test.tsx`
Expected: FAIL — worksheet test-ids (`scoping-plate`, `scoping-autogroup`, `scoping-undo`, `scoping-discard`, `scoping-loose-11`) don't exist yet.

- [ ] **Step 4: Rewrite `SeriesScopingPage.tsx`**

```tsx
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useNavigate } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import {
  useCorpusSampleTags,
  useCorpusPickerSamples,
  useScopeSeries,
  useMemberTraces,
  useMemberIndices,
} from "../queries";
import { proposeOrdering, type OrderingRow } from "../lib/scoping/proposeOrdering";
import { splitProposal, humanizeKey } from "../lib/scoping/splitProposal";
import { parseSortKey } from "../lib/scoping/parseSortKey";
import { dominantPhase } from "../lib/scoping/dominantPhase";
import { ScopingAutogroupCard } from "../components/ScopingAutogroupCard";
import { ScopingOrderField } from "../components/ScopingOrderField";
import { ScopingRow } from "../components/ScopingRow";
import { ScopingLooseMatches } from "../components/ScopingLooseMatches";
import { ScopingPhaseStrip } from "../components/ScopingPhaseStrip";
import { ScopingFoot } from "../components/ScopingFoot";
import { ScopingConfirmModal } from "../components/ScopingConfirmModal";

const SCOPING_FIXTURE = (
  <div className="space-y-2">
    {[0, 1, 2, 3].map((i) => (
      <div key={i} className="flex items-center gap-3 rounded border border-hair p-3">
        <div className="h-4 w-4 rounded bg-paper-sunk" />
        <div className="h-4 w-1/3 rounded bg-paper-sunk" />
        <div className="h-4 w-1/4 rounded bg-paper-sunk" />
      </div>
    ))}
  </div>
);

interface UndoEntry {
  rows: OrderingRow[];
  loose: OrderingRow[];
  label: string;
}

export function SeriesScopingPage(): JSX.Element {
  const navigate = useNavigate();
  const tagsQ = useCorpusSampleTags();
  const pickerQ = useCorpusPickerSamples();
  const scopeSeries = useScopeSeries();

  const isLoading = tagsQ.isLoading || pickerQ.isLoading;
  const isError = tagsQ.isError || pickerQ.isError;

  const proposal = useMemo(
    () => proposeOrdering(tagsQ.data ?? [], pickerQ.data ?? []),
    [tagsQ.data, pickerQ.data],
  );
  const split = useMemo(() => splitProposal(proposal), [proposal]);

  // Local editable copies — seeded from the split once queries resolve, then
  // user-owned. `rows` = members; `loose` = Himalaya-also-found.
  const [rows, setRows] = useState<OrderingRow[]>([]);
  const [loose, setLoose] = useState<OrderingRow[]>([]);
  useEffect(() => {
    setRows(split.members);
    setLoose(split.looseMatches);
  }, [split.members, split.looseMatches]);

  // In-session undo stack (S-F). Each mutating action snapshots prior state.
  const [history, setHistory] = useState<UndoEntry[]>([]);
  const snapshot = useCallback((label: string) => {
    setHistory((h) => [...h, { rows, loose, label }]);
  }, [rows, loose]);
  const undo = useCallback(() => {
    setHistory((h) => {
      if (h.length === 0) return h;
      const last = h[h.length - 1];
      setRows(last.rows);
      setLoose(last.loose);
      return h.slice(0, -1);
    });
  }, []);

  useEffect(() => {
    function onKey(e: KeyboardEvent): void {
      if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === "z") {
        e.preventDefault();
        undo();
      }
    }
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [undo]);

  // Member ordering: render low→high by parsed sort key; unparseable last,
  // stable by sampleId. Pure surfacing — does NOT reorder the batch payload.
  const orderedRows = useMemo(() => {
    return [...rows].sort((a, b) => {
      const ka = parseSortKey(a.value), kb = parseSortKey(b.value);
      if (ka === null && kb === null) return a.sampleId - b.sampleId;
      if (ka === null) return 1;
      if (kb === null) return -1;
      return ka - kb || a.sampleId - b.sampleId;
    });
  }, [rows]);

  // Per-member trace + phase fetches keyed on the indexing exposure id.
  const pickerById = useMemo(() => {
    const m = new Map<number, number | null>();
    for (const s of pickerQ.data ?? []) m.set(s.sample.id, s.indexing_exposure_id);
    return m;
  }, [pickerQ.data]);
  const allSampleIds = useMemo(() => [...rows, ...loose].map((r) => r.sampleId), [rows, loose]);
  const exposureIds = useMemo(
    () => allSampleIds.map((id) => pickerById.get(id)).filter((e): e is number => e != null),
    [allSampleIds, pickerById],
  );
  const tracesByExposure = useMemoTraces(exposureIds);
  const indicesByExposure = useMemoIndices(exposureIds);

  // Re-key trace/phase maps from exposure id → sample id for the row components.
  const traceBySample = useMemo(() => {
    const m = new Map<number, import("../api").Trace>();
    for (const id of allSampleIds) {
      const eid = pickerById.get(id);
      if (eid != null && tracesByExposure.has(eid)) m.set(id, tracesByExposure.get(eid)!);
    }
    return m;
  }, [allSampleIds, pickerById, tracesByExposure]);
  const phaseBySample = useMemo(() => {
    const m = new Map<number, string | null>();
    for (const id of allSampleIds) {
      const eid = pickerById.get(id);
      const idx = eid != null ? indicesByExposure.get(eid) : undefined;
      m.set(id, idx ? dominantPhase(idx).dominant : null);
    }
    return m;
  }, [allSampleIds, pickerById, indicesByExposure]);

  const phaseReads = useMemo(
    () => orderedRows.map((r) => {
      const eid = pickerById.get(r.sampleId);
      const idx = eid != null ? indicesByExposure.get(eid) : undefined;
      return idx ? dominantPhase(idx) : { dominant: null, coexist: null };
    }),
    [orderedRows, pickerById, indicesByExposure],
  );

  const [confirmOpen, setConfirmOpen] = useState(false);

  const flagCount = rows.filter((r) => r.flagged).length;
  const included = rows.filter((r) => r.include && !r.flagged);
  const canBuild =
    proposal.orderingKey !== undefined &&
    included.length > 0 &&
    rows.every((r) => !r.include || !r.flagged);

  function changeValue(sampleId: number, value: string): void {
    snapshot(`edited smp_${sampleId}`);
    setRows((prev) => prev.map((r) =>
      r.sampleId === sampleId ? { ...r, value, flagged: value === "" } : r));
  }
  function toggleFlag(sampleId: number): void {
    snapshot(`toggled smp_${sampleId}`);
    setRows((prev) => prev.map((r) =>
      r.sampleId === sampleId ? { ...r, flagged: !r.flagged } : r));
  }
  function addLoose(sampleId: number): void {
    snapshot(`added smp_${sampleId}`);
    setLoose((prev) => prev.filter((r) => r.sampleId !== sampleId));
    setRows((prev) => {
      const found = loose.find((r) => r.sampleId === sampleId);
      return found ? [...prev, { ...found, include: true }] : prev;
    });
  }

  const pendingBuildRef = useRef(false);

  function handleConfirm(): void {
    if (proposal.orderingKey === undefined) return;
    pendingBuildRef.current = true;
    scopeSeries.mutate({
      key: proposal.orderingKey,
      tags: included.map((r) => ({ sampleId: r.sampleId, value: r.value })),
    });
    setConfirmOpen(false);
  }

  useEffect(() => {
    if (!scopeSeries.isSuccess || !pendingBuildRef.current) return;
    pendingBuildRef.current = false;
    navigate("/series");
  }, [scopeSeries.isSuccess, navigate]);

  useEffect(() => {
    if (scopeSeries.error) pendingBuildRef.current = false;
  }, [scopeSeries.error]);

  const keyLabel = humanizeKey(proposal.orderingKey);

  return (
    <div data-testid="scoping-page" className="flex flex-1 flex-col items-center overflow-auto px-10 py-10">
      <div className="flex w-full max-w-[760px] items-center justify-end pb-3">
        <button
          type="button"
          data-testid="scoping-discard"
          onClick={() => navigate("/series")}
          className="text-xs font-semibold text-ink-faint hover:text-ink"
        >
          Discard
        </button>
      </div>

      <div
        data-testid="scoping-plate"
        className="w-full max-w-[760px] rounded-md border border-hair bg-plate px-8 py-7 shadow-[0_1px_1px_rgba(60,52,40,.04),0_18px_40px_-22px_rgba(60,52,40,.22)]"
      >
        <div className="mb-1.5 text-[11px] font-bold uppercase tracking-[0.14em] text-print-accent">
          New series
        </div>
        <h1 data-testid="scoping-title" className="text-display text-ink">
          {proposal.orderingKey ? `Series by ${keyLabel}` : "New series"}
        </h1>

        {scopeSeries.error ? (
          <div
            data-testid="scoping-error-banner"
            role="alert"
            className="mt-4 rounded border border-print-accent bg-paper-sunk px-4 py-2 text-sm text-print-accent"
          >
            Could not write the scoping tags. Nothing was saved — adjust and try Confirm &amp; build again.
          </div>
        ) : null}

        {isError ? (
          <div data-testid="scoping-error" className="px-4 py-8 text-sm text-ink-soft">
            Could not load corpus tags. Try reloading the page.
          </div>
        ) : (
          <Skeleton name="scoping" className="w-full" loading={isLoading} stagger={50} transition={200} fixture={SCOPING_FIXTURE}>
            {proposal.orderingKey === undefined ? (
              <div data-testid="scoping-empty" className="px-4 py-12 text-center text-sm text-ink-faint">
                {rows.length === 0 && loose.length === 0
                  ? "No samples in the corpus to scope."
                  : "No shared ordering variable yet — tag these samples to propose a series."}
              </div>
            ) : (
              <>
                <ScopingAutogroupCard memberCount={rows.length} keyLabel={keyLabel} flagCount={flagCount} />
                <ScopingOrderField keyLabel={keyLabel} />

                <div className="mb-1 mt-6 flex items-baseline justify-between">
                  <span className="text-[10.5px] font-bold uppercase tracking-wider text-ink-faint">The series</span>
                  <span className="flex items-baseline gap-3.5">
                    {history.length > 0 ? (
                      <button
                        type="button"
                        data-testid="scoping-undo"
                        onClick={undo}
                        title={`Step back: ${history[history.length - 1].label}`}
                        className="text-[11px] font-semibold text-print-accent hover:underline"
                      >
                        ↺ Undo last change
                      </button>
                    ) : null}
                    <span className="font-mono text-[10.5px] text-ink-faint">
                      {rows.length} samples · low to high
                    </span>
                  </span>
                </div>

                <div>
                  {orderedRows.map((r) => (
                    <ScopingRow
                      key={r.sampleId}
                      row={r}
                      trace={traceBySample.get(r.sampleId)}
                      phase={phaseBySample.get(r.sampleId) ?? null}
                      onChangeValue={(v) => changeValue(r.sampleId, v)}
                      onToggleFlag={() => toggleFlag(r.sampleId)}
                    />
                  ))}
                </div>

                <ScopingLooseMatches rows={loose} traces={traceBySample} phases={phaseBySample} onAdd={addLoose} />
                <ScopingPhaseStrip reads={phaseReads} />
                <ScopingFoot
                  flagCount={flagCount}
                  memberCount={rows.length}
                  keyLabel={keyLabel}
                  canBuild={canBuild}
                  onBuild={() => setConfirmOpen(true)}
                />
              </>
            )}
          </Skeleton>
        )}
      </div>

      <ScopingConfirmModal
        open={confirmOpen}
        orderingKey={proposal.orderingKey}
        count={included.length}
        onConfirm={handleConfirm}
        onClose={() => setConfirmOpen(false)}
      />
    </div>
  );
}

// `useMemberTraces` / `useMemberIndices` are stable-mapped over the exposure-id
// list; aliasing here keeps the JSX body readable.
function useMemoTraces(ids: number[]) { return useMemberTraces(ids); }
function useMemoIndices(ids: number[]) { return useMemberIndices(ids); }
```

> Drop the two trivial alias wrappers if the reviewer prefers direct calls; they
> exist only to keep imports adjacent. Direct `useMemberTraces(exposureIds)` /
> `useMemberIndices(exposureIds)` in the body is equivalent.

- [ ] **Step 5: Run the full scoping suite to verify it passes**

Run: `npx vitest run test/scoping.test.tsx test/proposeOrdering.test.ts`
Expected: PASS (all member/gating/undo/loose/discard specs).

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/pages/SeriesScopingPage.tsx packages/HimalayaUI/frontend/test/scoping.test.tsx
git commit -m "feat(scoping): worksheet plate — autogroup/sparklines/loose/preview/gate (R7 #230 S-A/S-D..S-H)"
```

---

## Task 15: Full build + suite gate

**Files:** none (verification only).

- [ ] **Step 1: Type-check + build**

Run: `cd packages/HimalayaUI/frontend && npm run build`
Expected: `tsc --noEmit` clean, `vite build` succeeds.

- [ ] **Step 2: Full unit suite**

Run: `npm test`
Expected: all green (new scoping suites + untouched suites).

- [ ] **Step 3: Fix any fallout** (e.g. a renamed test-id referenced by `queue/scopeSeries.contract.test.ts` or `queue/applyRemoteToCache.scoping.test.ts` — those test the mutator/cache, NOT the page DOM, so they should be unaffected; confirm).

- [ ] **Step 4: Commit any fixes**

```bash
git commit -am "test(scoping): suite + build green after worksheet rewrite (R7 #230)"
```

---

## Task 16: Live verification (warm + cold) — screenshots

**Files:** none (verification only). See agent runbook in the task brief.

- [ ] **Step 1: Seed a private DB copy with a numeric tag key** so `proposeOrdering` has an ordering variable (cold shared DB has 0 sample_tags).

```bash
cp /Users/me/.claude/jobs/01683885/himalaya-dev.db /tmp/r7-dev.db
sqlite3 /tmp/r7-dev.db '.schema sample_tags'
# INSERT (key='concentration', value=<varying>, source='manual') for ~5 samples in one experiment
```

- [ ] **Step 2: Start a private backend on 8092**, wait for `/api/samples`.

- [ ] **Step 3: Warm path** — Vite on 5194 with `VITE_API_PORT=8092`; screenshot `http://127.0.0.1:5194/series/new`; compare to `series-scoping.html`. Save `/tmp/r7-warm.png`.

- [ ] **Step 4: Cold path** — re-point at :8091 (0 tags); confirm the cold fallback ("No shared ordering variable yet…"). Save `/tmp/r7-cold.png`.

- [ ] **Step 5: Tear down** — `lsof -ti:5194 | xargs -r kill`; `pkill -f "port 8092"`. Leave :8091 running.

---

## Self-review notes (spec coverage)

| Finding | Task |
|---|---|
| S-A autogroup summary | 9, 14 |
| S-A "Ordered by" parsed field | 2 (humanizeKey), 9, 14 |
| S-A per-row trace sparklines | 4, 6, 8, 14 |
| S-A amber check-the-read flags | 7, 8 |
| S-A "Himalaya also found" | 2 (split), 10, 14 |
| S-A preview phase strip | 3, 5, 11, 14 |
| S-A narrative gate footer | 12, 14 |
| S-B build button Print ink | 12, 13 |
| S-C greyed-ink disabled | 12 |
| S-D 760px worksheet plate | 14 |
| S-E confirm-not-fill-out | 7 |
| S-F reorder grip + low→high + Undo/⌘Z | 8 (grip), 14 (ordering + undo) |
| S-G serif series name + autogroup sentence + footnote | 9, 14 |
| S-H Discard affordance | 14 |

**Deferred (noted, not gold-plated):** drag-to-reorder is surfaced as a grip handle + automatic low→high sort (full DnD reorder of the durable key is out of scope — ordering is derived from the value, not a manual permutation). Re-selecting the ordering variable from the "Ordered by" field is presented as a decision affordance but is read-only this milestone (no UI exists to pick a different corpus key; follow-up).
