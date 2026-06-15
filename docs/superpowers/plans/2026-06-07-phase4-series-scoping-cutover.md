# Phase-4 Series-Scoping Cutover Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the legacy `/series/new` scoping page body with a greenfield `src/print/pages/SeriesScopingPage.tsx` assembled from the already-built `src/print` composites (`ScopePlate` / `ScopeSampleRow` / `ScopeCandidateRow`), reusing the CARRIED scoping logic untouched and deleting the OLD presentation.

**Architecture:** This is a "wiring, not building" cutover — every greenfield composite already landed in Batch 12. The new page wires CARRIED queries (`useCorpusSampleTags` / `useCorpusPickerSamples` / `useMemberTraces` / `useMemberIndices`), CARRIED lib (`proposeOrdering` / `splitProposal` / `parseSortKey` / `dominantPhase`), and the CARRIED mutation (`useScopeSeries`) into the greenfield `ScopePlate`, plus a thin tested page-derivation module (`scopingDerive.ts`) for the foot state / build gate / payload / preview. The greenfield `ScopingAssembly.stories.tsx` page-simulation is the structural template; the legacy page is consulted only for the CARRIED-logic wiring (which hook, what the mutation posts), never for presentation. The cutover then repoints the route and grep-deletes the OLD page + components + tests.

**Tech Stack:** React 18 + TypeScript (strict, `exactOptionalPropertyTypes: true`), TanStack Query, React Router, boneyard-js Skeleton, Tailwind (token-only classes outside `src/print/ui/**` — `lint:design` enforced), Vitest + RTL, Playwright.

---

## Provenance ledger (NEW / OLD / CARRIED)

**NEW** (this plan creates, under `src/print/**`):
- `src/print/pages/scopingDerive.ts` — pure page derivations (foot state, build gate, payload, preview segments).
- `src/print/pages/SeriesScopingPage.tsx` — the greenfield page.
- Tests: `test/print-pages/scopingDerive.test.ts`, `test/print-pages/SeriesScopingPage.test.tsx`.

**CARRIED** (reuse, do NOT touch):
- Queries: `useCorpusSampleTags`, `useCorpusPickerSamples`, `useMemberTraces`, `useMemberIndices`, `useScopeSeries` (`src/queries.ts`).
- Lib: `src/lib/scoping/proposeOrdering.ts` (`proposeOrdering`, `OrderingRow`), `splitProposal.ts` (`splitProposal`, `humanizeKey`), `parseSortKey.ts` (`parseSortKey`), `dominantPhase.ts` (`dominantPhase`, `PhaseRead`).
- Greenfield composites (built, Batch 12): `src/print/components/ScopePlate.tsx`, `ScopeSampleRow.tsx`, `ScopeCandidateRow.tsx`, `AutoGroup.tsx`, `useDragReorder.ts` (`useDragReorder`, `reorder`), `src/print/plot/Sparkline.tsx`.
- Primitives: `src/print/ui/**` (`PhaseStrip`/`PhaseSegment`, `Field`, `FlagButton`, `GripHandle`, `Button`, `Kicker`, `Card`, `Dot`, `EmptyState`), `src/print/components/PageFrame.tsx` (`width="scoping"` = `max-w-[760px]`).

**OLD** (delete on cutover — enumerate by grep, do NOT trust a fixed list):
- `src/pages/SeriesScopingPage.tsx` + every `src/components/Scoping*.tsx` that only it consumes.
- Their unit tests under `test/Scoping*.test.tsx` + `test/scoping.test.tsx`.

---

## Scope decisions (Rule-3 / controls-don't-lie — settled; do NOT re-litigate mid-task)

These are the honesty calls the research surfaced. The mockup is a sketch; where it would make a control lie, the honest behaviour wins.

| # | Question | Decision | Why |
|---|---|---|---|
| 1 | Mockup says "You selected N samples on the contact sheet." | **Honest grouping copy** — describe what actually happens: "Himalaya grouped N samples … by their {key}, read from the sample names." | No contact-sheet→scoping selection is plumbed; `proposeOrdering` groups the whole corpus picker projection by the most-frequent tag key. Claiming a selection that isn't wired is a lie. |
| 2 | Legacy `ScopingRow` had an inline value-edit `<input>`. | **No inline editing.** Greenfield `ScopeSampleRow` exposes only `onToggleFlag`; the page must not reintroduce a text input. | The greenfield model is "machine parses the value → human *confirms the read* by clicking the flag," not "edit the number." Re-adding editing is exactly the legacy-affordance leak the no-legacy-port rule forbids. |
| 3 | Mockup shows a ▾ dropdown on the "Ordered by" field. | **Display-only order field** — pass `orderedBy={keyLabel}` static; omit `orderOptions`/`onOrderSelect`. | `proposeOrdering` auto-selects the ordering key and has no "force this key" path. A dropdown that re-labels the field without re-grouping would lie (pick "Time" → nothing reorders). Re-grouping by an arbitrary variable is a future feature. The page-sim's `ORDER_OPTIONS` are illustrative only. |
| 4 | Legacy opened a `ScopingConfirmModal` before writing. | **No confirm modal** — "Confirm & build" writes directly + navigates. | Mockup + page-sim build directly; the foot note already carries the consequence copy ("Confirming records the {key} on every sample…"). Tags are additive, non-destructive metadata. |
| 5 | Mockup shows row grip handles (drag-reorder). | **Keep drag-reorder, but display-only.** Wire `useDragReorder`; it rewrites the displayed `order` + the preview strip. It does NOT change the written `(key,value)` tags. | In the mockup + the page-sim proves it with the existing hook. Scoping writes *tags*, not a series (series creation is the builder, Plan 6), so order is never persisted here — and no copy claims it is. Keyboard-reorder a11y is a documented follow-up (grips are `aria-hidden`). |
| 6 | Which load/error states must survive? | **All four honest states:** corpus-load error, loading skeleton, no-shared-ordering-variable empty (+ contact-sheet CTA), and write-error banner. | The folio cutover taught that a swallowed `isError` masquerades as an empty result. Each state is distinct and load-bearing. |

---

## Revision during implementation — Option A: honest commit-gate

The code-quality review surfaced that the plan's flag + candidate-add model was wrong against the real data contract, and the user chose **Option A** (2026-06-07). The load-bearing fact, confirmed in `proposeOrdering.ts`/`splitProposal.ts`:

> A row's `flagged` ⟺ `value === ""` ⟺ it is a **loose match** (`include:false`). So **members always have a value and are never machine-flagged**, and **candidates always have `value === ""`**. There is no "uncertain parse" middle state. Scoping's only durable effect is **writing the parsed `(key,value)` reads onto member samples** — it does not create a series.

Option A reconciles the page to this:
- **Flag = "skip this read from the write"** (excluded, *never* blocking) — not "confirm/fix a read." Members start kept; clicking skips. `canScopeBuild` = key exists ∧ ≥1 kept member (`include && !flagged && value !== ""`). `buildScopePayload` writes only kept members and **guards `value !== ""`** (an empty value can never reach the backend → no corruption).
- **No inline candidate "+ Add"** — a loose match has no value and we don't edit values, so it can never become a committable member. Candidates render as **informational discovery** (a buttonless list: dimmed sparkline + name + "lacks the {key} — tag it on the contact sheet if it belongs"). `ScopeCandidateRow` renders its add button unconditionally, so the page does NOT use that composite for candidates (avoids mutating a shared primitive); the composite is left intact for a future force-add feature.
- **Foot/grouping copy** reframed honestly ("Confirm the reads and build — skip any it misread"; "records the {key} on every **kept** sample"). Supersedes the plan's grouping-copy and `buildFootState(flagCount, memberCount)` (now `buildFootState(keptCount, skippedCount)`).
- Deferred: a distinct "skipped" row treatment (currently reuses the flagged accent wash).

Implemented in commit `b73da2c` (on top of `cfe17a4`). The Task 1/2 code blocks below are the pre-revision form — follow this section where they conflict.

## Task 1: Page-derivation module (`scopingDerive.ts`)

Pure, tested page derivations extracted from the legacy page's inline logic so the build gate and foot state are unit-testable independent of React.

**Files:**
- Create: `src/print/pages/scopingDerive.ts`
- Test: `test/print-pages/scopingDerive.test.ts`

- [ ] **Step 1: Write the failing test**

```ts
// test/print-pages/scopingDerive.test.ts
import { describe, it, expect } from "vitest";
import {
  buildFootState,
  canScopeBuild,
  buildScopePayload,
  toPreviewSegments,
} from "../../src/print/pages/scopingDerive";
import type { OrderingRow } from "../../src/lib/scoping/proposeOrdering";
import type { PhaseRead } from "../../src/lib/scoping/dominantPhase";

const row = (over: Partial<OrderingRow>): OrderingRow => ({
  sampleId: 1,
  sampleName: "s",
  value: "1 : 0.5",
  flagged: false,
  include: true,
  ...over,
});

describe("buildFootState", () => {
  it("is ready when no flags, with the member count", () => {
    expect(buildFootState(0, 6)).toEqual({
      kind: "ready",
      text: "All 6 values confirmed — ready to build",
    });
  });
  it("warns with singular/plural flag wording", () => {
    expect(buildFootState(1, 6)).toEqual({
      kind: "warn",
      text: "1 value to check before you can build",
    });
    expect(buildFootState(2, 6).text).toBe("2 values to check before you can build");
  });
});

describe("canScopeBuild", () => {
  it("is false without an ordering key", () => {
    expect(canScopeBuild([row({})], undefined)).toBe(false);
  });
  it("is false with no included members", () => {
    expect(canScopeBuild([row({ include: false })], "ratio")).toBe(false);
  });
  it("is false when an included member is flagged", () => {
    expect(canScopeBuild([row({ flagged: true })], "ratio")).toBe(false);
  });
  it("is true when all included members are unflagged (flagged loose ignored)", () => {
    expect(
      canScopeBuild([row({ sampleId: 1 }), row({ sampleId: 2, include: false, flagged: true })], "ratio"),
    ).toBe(true);
  });
});

describe("buildScopePayload", () => {
  it("writes only included, unflagged members as {sampleId, value}", () => {
    const rows = [
      row({ sampleId: 1, value: "1 : 0" }),
      row({ sampleId: 2, flagged: true }),
      row({ sampleId: 3, include: false }),
    ];
    expect(buildScopePayload(rows)).toEqual([{ sampleId: 1, value: "1 : 0" }]);
  });
});

describe("toPreviewSegments", () => {
  it("omits coexistWith when there is no second phase", () => {
    const reads: PhaseRead[] = [{ dominant: "Pn3m", coexist: null }];
    expect(toPreviewSegments(reads)).toEqual([{ phase: "Pn3m" }]);
  });
  it("wraps the coexist partner in an array", () => {
    const reads: PhaseRead[] = [{ dominant: "Pn3m", coexist: "Lamellar" }];
    expect(toPreviewSegments(reads)).toEqual([{ phase: "Pn3m", coexistWith: ["Lamellar"] }]);
  });
  it("passes a null dominant through as an unindexed cell", () => {
    expect(toPreviewSegments([{ dominant: null, coexist: null }])).toEqual([{ phase: null }]);
  });
});
```

- [ ] **Step 2: Run it; verify it fails** — `npx vitest run test/print-pages/scopingDerive.test.ts` → FAIL (module not found).

- [ ] **Step 3: Implement**

```ts
// src/print/pages/scopingDerive.ts
import type { OrderingRow } from "../../lib/scoping/proposeOrdering";
import type { PhaseRead } from "../../lib/scoping/dominantPhase";
import type { PhaseSegment } from "../ui";

export interface FootState {
  kind: "warn" | "ready";
  text: string;
}

/** The confirm-gate state line. Warn while any value still needs a look; ready
 *  once all are confirmed. `memberCount` is the member total (the denominator
 *  in the "All N values confirmed" copy). */
export function buildFootState(flagCount: number, memberCount: number): FootState {
  if (flagCount > 0) {
    return {
      kind: "warn",
      text: `${flagCount} value${flagCount === 1 ? "" : "s"} to check before you can build`,
    };
  }
  return { kind: "ready", text: `All ${memberCount} values confirmed — ready to build` };
}

/** Build gate (CARRIED contract from the legacy page): an ordering key must
 *  exist, at least one member must be included, and no included member may be
 *  flagged. The batch route 400s on an empty array. */
export function canScopeBuild(rows: OrderingRow[], orderingKey: string | undefined): boolean {
  if (orderingKey === undefined) return false;
  const included = rows.filter((r) => r.include && !r.flagged);
  return included.length > 0 && rows.every((r) => !r.include || !r.flagged);
}

/** The batch payload: only included, non-flagged members are written. */
export function buildScopePayload(rows: OrderingRow[]): { sampleId: number; value: string }[] {
  return rows
    .filter((r) => r.include && !r.flagged)
    .map((r) => ({ sampleId: r.sampleId, value: r.value }));
}

/** Map per-member phase reads (dominant + optional coexist partner) onto the
 *  PhaseStrip preview segments. A null dominant stays a null (unindexed) cell;
 *  a coexist partner becomes the single-element `coexistWith` array (the
 *  greenfield strip takes `string[]`, not the legacy `string | null`). */
export function toPreviewSegments(reads: PhaseRead[]): PhaseSegment[] {
  return reads.map((r) =>
    r.coexist ? { phase: r.dominant, coexistWith: [r.coexist] } : { phase: r.dominant },
  );
}
```

- [ ] **Step 4: Run it; verify it passes** — `npx vitest run test/print-pages/scopingDerive.test.ts` → PASS.

- [ ] **Step 5: Gate + commit**

```bash
npx tsc --noEmit -p tsconfig.build.json   # 0 errors
npm run lint:design                        # 0 violations
git add src/print/pages/scopingDerive.ts test/print-pages/scopingDerive.test.ts
git commit -m "$(cat <<'EOF'
Greenfield scoping: pure page-derivation module (foot state / build gate / payload / preview)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: The greenfield page (`SeriesScopingPage.tsx`)

Wire the CARRIED data + Task-1 derivations into `ScopePlate`, mirroring the page-sim's local state (flag toggle, candidate fold, undo + ⌘Z, display-only drag-reorder) and the legacy's CARRIED-logic wiring (trace/phase maps, scope mutation + navigate). Honour all six scope decisions.

**Files:**
- Create: `src/print/pages/SeriesScopingPage.tsx`
- Test: `test/print-pages/SeriesScopingPage.test.tsx`
- Reference (read, do not edit): `src/print/components/ScopingAssembly.stories.tsx` (structural template), `src/print/pages/SeriesFolioPage.tsx` (sibling page pattern), `test/print-pages/SeriesFolioPage.test.tsx` (sibling test / query-mock pattern), `packages/HimalayaUI/frontend/test/AGENTS.md` (Vitest/JSDOM rules).

- [ ] **Step 1: Write the failing test (core assertions; add the rest after the page exists)**

Mirror the query-mocking harness in `test/print-pages/SeriesFolioPage.test.tsx`. Mock `../../src/queries` so `useCorpusSampleTags` / `useCorpusPickerSamples` return resolved data, `useMemberTraces` / `useMemberIndices` return `Map`s, and `useScopeSeries` returns a controllable `{ mutate, isSuccess, error }`. Render inside a `MemoryRouter`.

```tsx
// test/print-pages/SeriesScopingPage.test.tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";

const mutate = vi.fn();
let scopeState = { mutate, isSuccess: false, error: null as Error | null };
let tagsState = { data: [{ key: "ratio", value: "1 : 0.5" }], isLoading: false, isError: false };
let pickerState = {
  data: [
    { sample: { id: 1, name: "A", display_name: "A", tags: [{ key: "ratio", value: "1 : 0" }], experiment_id: 1, notes: "" }, indexing_exposure_id: 37, all_exposures: [] },
    { sample: { id: 2, name: "B", display_name: "B", tags: [{ key: "ratio", value: "1 : 0.5" }], experiment_id: 1, notes: "" }, indexing_exposure_id: 65, all_exposures: [] },
  ],
  isLoading: false,
  isError: false,
};

vi.mock("../../src/queries", () => ({
  useCorpusSampleTags: () => tagsState,
  useCorpusPickerSamples: () => pickerState,
  useScopeSeries: () => scopeState,
  useMemberTraces: () => new Map([[37, { q: [0.1, 0.2], I: [1, 2], sigma: [0, 0] }], [65, { q: [0.1, 0.2], I: [2, 1], sigma: [0, 0] }]]),
  useMemberIndices: () => new Map([[37, [{ id: 1, exposure_id: 37, phase: "Pn3m", score: 0.9 }]], [65, [{ id: 2, exposure_id: 65, phase: "Pn3m", score: 0.8 }]]]),
}));

import { SeriesScopingPage } from "../../src/print/pages/SeriesScopingPage";

function renderPage(): void {
  const qc = new QueryClient();
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/series/new"]}>
        <SeriesScopingPage />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

beforeEach(() => {
  mutate.mockReset();
  scopeState = { mutate, isSuccess: false, error: null };
  tagsState = { data: [{ key: "ratio", value: "1 : 0.5" }], isLoading: false, isError: false };
});

describe("SeriesScopingPage", () => {
  it("renders a scope-sample-row per member", () => {
    renderPage();
    expect(screen.getAllByTestId("scope-sample-row")).toHaveLength(2);
  });

  it("surfaces a corpus-load error distinctly (not as an empty state)", () => {
    tagsState = { data: [], isLoading: false, isError: true };
    renderPage();
    expect(screen.getByText(/couldn't load the corpus/i)).toBeInTheDocument();
    expect(screen.queryByTestId("scope-sample-row")).not.toBeInTheDocument();
  });

  it("writes only confirmed members and navigates on build success", () => {
    renderPage();
    fireEvent.click(screen.getByRole("button", { name: /confirm & build/i }));
    expect(mutate).toHaveBeenCalledWith(
      expect.objectContaining({ key: "ratio", tags: expect.arrayContaining([{ sampleId: 2, value: "1 : 0.5" }]) }),
    );
  });
});
```

> NOTE for the implementer: align the mock object shapes (`SampleTagPair`, `PickerSampleRow`, `Trace`, `IndexEntry`, and the `useScopeSeries` return) with the real exports in `src/queries.ts` / `src/api.ts`. The shapes above are representative — read the real types and the sibling `SeriesFolioPage.test.tsx` harness before finalising. After the page exists, add assertions for: the **warn vs ready** foot line, **build disabled** while a member is flagged, **flag toggle** flipping a row, the **no-ordering-variable empty state** + contact-sheet CTA, and the **write-error banner** (`scopeState.error` set).

- [ ] **Step 2: Run it; verify it fails** — `npx vitest run test/print-pages/SeriesScopingPage.test.tsx` → FAIL (module not found).

- [ ] **Step 3: Implement the page**

```tsx
// src/print/pages/SeriesScopingPage.tsx
import { useEffect, useMemo, useRef, useState } from "react";
import { useNavigate } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { PageFrame } from "../components/PageFrame";
import { ScopePlate } from "../components/ScopePlate";
import { ScopeSampleRow } from "../components/ScopeSampleRow";
import { ScopeCandidateRow } from "../components/ScopeCandidateRow";
import { useDragReorder, reorder } from "../components/useDragReorder";
import { EmptyState } from "../ui";
import type { Trace } from "../../api";
import {
  useCorpusSampleTags,
  useCorpusPickerSamples,
  useScopeSeries,
  useMemberTraces,
  useMemberIndices,
} from "../../queries";
import { proposeOrdering, type OrderingRow } from "../../lib/scoping/proposeOrdering";
import { splitProposal, humanizeKey } from "../../lib/scoping/splitProposal";
import { parseSortKey } from "../../lib/scoping/parseSortKey";
import { dominantPhase } from "../../lib/scoping/dominantPhase";
import { buildFootState, canScopeBuild, buildScopePayload, toPreviewSegments } from "./scopingDerive";

const EMPTY_TRACE: Trace = { q: [], I: [], sigma: [] };

// Token-only skeleton fixture (no inline appearance literals — design-guard clean).
// No scoping.bones.json capture exists yet (deferred, needs the data volume).
const SCOPING_FIXTURE = (
  <div className="space-y-2">
    {[0, 1, 2, 3].map((i) => (
      <div key={i} className="flex items-center gap-3 rounded border border-hair p-3">
        <div className="h-4 w-4 rounded bg-paper-sunk" />
        <div className="h-4 w-1/3 rounded bg-paper-sunk" />
        <div className="ml-auto h-4 w-1/5 rounded bg-paper-sunk" />
      </div>
    ))}
  </div>
);

type HistoryEntry =
  | { type: "flag"; id: number; prev: boolean; label: string }
  | { type: "add"; id: number; label: string };

/**
 * SeriesScopingPage (greenfield) — the machine-proposes / human-confirms
 * scoping worksheet at /series/new. The confirm-and-build GATE that *writes*
 * the structured (key,value) sample_tags — NOT series creation (that is the
 * builder). Assembled from src/print composites + the series-scoping mockup;
 * carried logic only (proposeOrdering/splitProposal/dominantPhase + the
 * useScopeSeries batch write), no legacy presentation. See the plan's scope
 * table for the six honesty decisions (display-only order field, no confirm
 * modal, display-only drag-reorder, honest grouping copy, flag-not-edit).
 */
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
  const keyLabel = humanizeKey(proposal.orderingKey);

  // Default value-sorted member order (low → high; unparseable last, stable by id).
  const seededOrder = useMemo(
    () =>
      [...split.members]
        .sort((a, b) => {
          const ka = parseSortKey(a.value);
          const kb = parseSortKey(b.value);
          if (ka === null && kb === null) return a.sampleId - b.sampleId;
          if (ka === null) return 1;
          if (kb === null) return -1;
          return ka - kb || a.sampleId - b.sampleId;
        })
        .map((r) => r.sampleId),
    [split.members],
  );

  // Local, user-owned copies seeded once the proposal resolves.
  const [rows, setRows] = useState<OrderingRow[]>([]);
  const [loose, setLoose] = useState<OrderingRow[]>([]);
  const [order, setOrder] = useState<number[]>([]);
  const [history, setHistory] = useState<HistoryEntry[]>([]);
  useEffect(() => {
    setRows(split.members);
    setLoose(split.looseMatches);
    setOrder(seededOrder);
    setHistory([]);
  }, [split.members, split.looseMatches, seededOrder]);

  // Display-only manual reorder (scope decision #5): rewrites `order` + the
  // preview; never touches the written (key,value) payload.
  const { dragItemProps, dropEdge } = useDragReorder((from, to) =>
    setOrder((o) => reorder(o, from, to)),
  );

  const byId = useMemo(() => new Map(rows.map((r) => [r.sampleId, r])), [rows]);
  const sorted = useMemo(
    () => order.map((id) => byId.get(id)).filter((r): r is OrderingRow => r != null),
    [order, byId],
  );

  // Clicking a value confirms / re-opens the machine's read (flag-not-edit,
  // scope decision #2). Recorded so Undo / ⌘Z steps it back.
  const toggleFlag = (id: number): void => {
    const m = rows.find((r) => r.sampleId === id);
    if (!m) return;
    setHistory((h) => [...h, { type: "flag", id, prev: m.flagged, label: `smp_${id}` }]);
    setRows((cur) => cur.map((r) => (r.sampleId === id ? { ...r, flagged: !r.flagged } : r)));
  };

  // Folding a candidate in: include it and APPEND to the displayed order
  // (manual order wins — never re-sort).
  const addCandidate = (id: number): void => {
    const c = loose.find((r) => r.sampleId === id);
    if (!c) return;
    setLoose((cur) => cur.filter((r) => r.sampleId !== id));
    setRows((cur) => [...cur, { ...c, include: true }]);
    setOrder((o) => [...o, id]);
    setHistory((h) => [...h, { type: "add", id, label: `smp_${id}` }]);
  };

  const undo = (): void => {
    setHistory((h) => {
      const e = h[h.length - 1];
      if (!e) return h;
      if (e.type === "flag") {
        setRows((cur) => cur.map((r) => (r.sampleId === e.id ? { ...r, flagged: e.prev } : r)));
      } else {
        setRows((cur) => {
          const m = cur.find((r) => r.sampleId === e.id);
          if (m) setLoose((ls) => [...ls, { ...m, include: false }]);
          return cur.filter((r) => r.sampleId !== e.id);
        });
        setOrder((o) => o.filter((id) => id !== e.id));
      }
      return h.slice(0, -1);
    });
  };

  useEffect(() => {
    function onKey(e: KeyboardEvent): void {
      if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === "z") {
        e.preventDefault();
        undo();
      }
    }
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  });

  // Trace + dominant-phase maps, keyed exposure → sample (CARRIED wiring).
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
  const tracesByExposure = useMemberTraces(exposureIds);
  const indicesByExposure = useMemberIndices(exposureIds);

  const traceBySample = useMemo(() => {
    const m = new Map<number, Trace>();
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

  // Preview strip follows the displayed order.
  const preview = useMemo(
    () =>
      toPreviewSegments(
        sorted.map((r) => {
          const eid = pickerById.get(r.sampleId);
          const idx = eid != null ? indicesByExposure.get(eid) : undefined;
          return idx ? dominantPhase(idx) : { dominant: null, coexist: null };
        }),
      ),
    [sorted, pickerById, indicesByExposure],
  );

  const flagCount = rows.filter((r) => r.flagged).length;
  const footState = buildFootState(flagCount, rows.length);
  const canBuild = canScopeBuild(rows, proposal.orderingKey);
  const lastLabel = history.length ? history[history.length - 1]!.label : undefined;

  // Defer navigation until the batch write actually succeeds (pending ref).
  const pendingBuildRef = useRef(false);
  const handleBuild = (): void => {
    if (proposal.orderingKey === undefined) return;
    pendingBuildRef.current = true;
    scopeSeries.mutate({ key: proposal.orderingKey, tags: buildScopePayload(rows) });
  };
  useEffect(() => {
    if (!scopeSeries.isSuccess || !pendingBuildRef.current) return;
    pendingBuildRef.current = false;
    navigate("/series");
  }, [scopeSeries.isSuccess, navigate]);
  useEffect(() => {
    if (scopeSeries.error) pendingBuildRef.current = false;
  }, [scopeSeries.error]);

  // ── State 1: corpus load failed (distinct from an empty result). ──────────
  if (isError) {
    return (
      <PageFrame width="scoping" className="px-10 py-10">
        <EmptyState
          title="Couldn't load the corpus"
          body="The sample tags failed to load. Try reloading the page."
        />
      </PageFrame>
    );
  }

  return (
    <PageFrame width="scoping" className="px-10 py-10">
      <div data-testid="scoping-page">
        {/* Discard (mockup top-right, outside the plate). */}
        <div className="flex items-center justify-end pb-3">
          <button
            type="button"
            data-testid="scoping-discard"
            onClick={() => navigate("/series")}
            className="px-1 py-1.5 text-xs font-semibold text-ink-faint hover:text-ink"
          >
            Discard
          </button>
        </div>

        {/* State 4: the batch write failed. */}
        {scopeSeries.error ? (
          <div
            data-testid="scoping-error-banner"
            role="alert"
            className="mb-4 rounded border border-print-accent bg-paper-sunk px-4 py-2 text-sm text-print-accent"
          >
            Could not write the scoping tags. Nothing was saved. Adjust and try Confirm &amp; build again.
          </div>
        ) : null}

        <Skeleton
          name="scoping"
          className="block w-full"
          loading={isLoading}
          stagger={50}
          transition={200}
          fixture={SCOPING_FIXTURE}
          fallback={<div className="p-8 text-sm text-ink-faint">Loading the worksheet…</div>}
        >
          {proposal.orderingKey === undefined ? (
            /* State 3: nothing shares an ordering variable yet. */
            <EmptyState
              title={rows.length === 0 && loose.length === 0 ? "Nothing to scope yet" : "No shared ordering variable"}
              body={
                rows.length === 0 && loose.length === 0
                  ? "No samples in the corpus to scope."
                  : "These samples share no ordering variable yet. Tag them on the contact sheet to propose a series."
              }
              action={{ label: "Open the contact sheet", onClick: () => navigate("/samples") }}
            />
          ) : (
            <ScopePlate
              seriesName={`Series by ${keyLabel}`}
              grouping={
                <>
                  Himalaya grouped <strong>{rows.length} samples</strong> by their{" "}
                  <strong>{keyLabel}</strong>, read from the sample names.
                  {flagCount > 0 ? (
                    <>
                      {" "}
                      {rows.length - flagCount} parsed cleanly, {flagCount}{" "}
                      {flagCount === 1 ? "needs" : "need"} a look.
                    </>
                  ) : (
                    <> All {rows.length} parsed cleanly.</>
                  )}
                </>
              }
              orderedBy={keyLabel}
              orderNote="Read from the sample names."
              count={`${rows.length} samples · low to high`}
              {...(history.length
                ? { onUndo: undo, ...(lastLabel ? { undoLabel: `Step back: ${lastLabel}` } : {}) }
                : {})}
              rows={sorted.map((r, i) => {
                const dprops = dragItemProps(i);
                const edge = dropEdge(i);
                return (
                  <div
                    key={r.sampleId}
                    {...dprops}
                    className={`relative cursor-grab${dprops["data-dragging"] ? " opacity-50" : ""}`}
                  >
                    {edge && (
                      <span
                        aria-hidden="true"
                        className={`pointer-events-none absolute left-0 right-0 z-10 h-0.5 rounded-full bg-accent ${edge === "top" ? "-top-px" : "-bottom-px"}`}
                      />
                    )}
                    <ScopeSampleRow
                      name={r.sampleName}
                      sampleId={`smp_${r.sampleId}`}
                      trace={traceBySample.get(r.sampleId) ?? EMPTY_TRACE}
                      phase={phaseBySample.get(r.sampleId) ?? null}
                      value={r.value}
                      {...(r.flagged ? { flagged: true } : {})}
                      onToggleFlag={() => toggleFlag(r.sampleId)}
                    />
                  </div>
                );
              })}
              candidates={
                loose.length ? (
                  loose.map((c) => (
                    <ScopeCandidateRow
                      key={c.sampleId}
                      name={c.sampleName}
                      why={
                        <>
                          lacks the <strong className="text-accent font-semibold">{keyLabel}</strong> — add
                          it if it belongs.
                        </>
                      }
                      trace={traceBySample.get(c.sampleId) ?? EMPTY_TRACE}
                      phase={phaseBySample.get(c.sampleId) ?? null}
                      onAdd={() => addCandidate(c.sampleId)}
                    />
                  ))
                ) : (
                  <div className="text-meta text-ink-faint italic">
                    Nothing else in the corpus matches this grouping.
                  </div>
                )
              }
              preview={preview}
              footState={footState}
              footNote={
                <>
                  Confirming records the {keyLabel} on every sample — the next series that needs it already
                  knows.
                </>
              }
              {...(!canBuild ? { buildDisabled: true } : {})}
              onBuild={handleBuild}
            />
          )}
        </Skeleton>
      </div>
    </PageFrame>
  );
}
```

> IMPLEMENTER NOTES:
> - **Verify primitive contracts before trusting this code.** Confirm `EmptyState` accepts an `action={{ label, onClick }}` prop (read `src/print/ui/EmptyState.tsx`); if its action API differs (e.g. `ctaLabel`/`onCta`, or a `ReactNode` slot), adapt the call — do not invent a prop. Confirm `Trace` includes `sigma` (adjust `EMPTY_TRACE` if not). Confirm `Sparkline` renders with an empty `{q:[],I:[]}` without throwing; if it can throw on empty input, guard by only rendering the row's trace when present, or render the EmptyState's neutral sparkline path.
> - `bg-accent` / `text-accent` / `text-print-accent` / `bg-paper-sunk` / `border-hair` are semantic token classes (allowed outside `ui/**`); only bracket-literals (`text-[…]`, `rounded-[…]`), raw hex, and side-stripes fail `lint:design`. The drag insertion line copies the page-sim's exact classes verbatim — keep them identical.
> - Do **not** add `orderOptions`/`onOrderSelect`/`onChangeOrder` to `ScopePlate` (scope decision #3 — display-only order field).

- [ ] **Step 4: Run the test; verify it passes** — `npx vitest run test/print-pages/SeriesScopingPage.test.tsx`. Add the deferred assertions (foot warn/ready, build-disabled, flag toggle, empty state, write-error banner) and make them pass.

- [ ] **Step 5: Gate + commit**

```bash
npx tsc --noEmit -p tsconfig.build.json
npm run lint:design
git add src/print/pages/SeriesScopingPage.tsx test/print-pages/SeriesScopingPage.test.tsx
git commit -m "$(cat <<'EOF'
Greenfield scoping: SeriesScopingPage assembled from ScopePlate + carried logic

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Cutover — repoint the route, grep-delete the OLD page + components + tests

Swap the route to the greenfield page and remove the legacy presentation. **Enumerate the OLD set by grep, not by a memorised list** (the folio cutover was bitten by a stale list + a narrow test filter).

**Files:**
- Modify: `src/components/AppRoutes.tsx`
- Delete: `src/pages/SeriesScopingPage.tsx` + the `src/components/Scoping*.tsx` it solely consumes + their tests.

- [ ] **Step 1: Repoint the route**

In `src/components/AppRoutes.tsx`, change the import from the legacy page to the greenfield one:
```diff
- import { SeriesScopingPage } from "../pages/SeriesScopingPage";
+ import { SeriesScopingPage } from "../print/pages/SeriesScopingPage";
```
Leave the `<Route path="/series/new" element={<SeriesScopingPage />} />` line unchanged.

- [ ] **Step 2: Enumerate the OLD set**

```bash
# Components the legacy page imported (the delete candidates):
grep -nE 'from "\.\./components/Scoping' src/pages/SeriesScopingPage.tsx
# Confirm nothing OUTSIDE the legacy page still imports each candidate before deleting it:
for f in ScopingAutogroupCard ScopingOrderField ScopingRow ScopingLooseMatches ScopingFoot ScopingConfirmModal ScopingSparkline ScopingValueCell ScopingPhaseStrip; do
  echo "== $f =="; grep -rln "components/$f\"" src --include=*.tsx | grep -v "src/pages/SeriesScopingPage.tsx"
done
```
Only delete a component when the second grep shows no remaining non-legacy consumer. (Anything still imported elsewhere stays — note it as a surprise and surface it.)

- [ ] **Step 3: Delete the OLD page, its sole-consumed components, and their tests**

```bash
git rm src/pages/SeriesScopingPage.tsx
# Delete each confirmed-orphan component + its test (adjust the list to the grep result):
git rm src/components/ScopingAutogroupCard.tsx src/components/ScopingOrderField.tsx \
       src/components/ScopingRow.tsx src/components/ScopingLooseMatches.tsx \
       src/components/ScopingFoot.tsx src/components/ScopingConfirmModal.tsx
git rm test/scoping.test.tsx test/ScopingAutogroup.test.tsx test/ScopingFoot.test.tsx \
       test/ScopingLooseMatches.test.tsx test/ScopingPhaseStrip.test.tsx test/ScopingRow.test.tsx \
       test/ScopingSparkline.test.tsx test/ScopingValueCell.test.tsx
```

- [ ] **Step 4: Guard — no greenfield→OLD import leaked, OLD fully unreferenced**

```bash
# src/print/pages must import NEW + CARRIED only, never src/pages or src/components:
grep -rnE 'from "(\.\./)+(pages|components)/' src/print/pages/SeriesScopingPage.tsx && echo "LEAK" || echo "clean"
# No dangling references to any deleted symbol:
grep -rn "ScopingConfirmModal\|ScopingAutogroupCard\|ScopingLooseMatches\|ScopingFoot\|ScopingOrderField" src test && echo "DANGLING" || echo "clean"
```
Both must report `clean`.

- [ ] **Step 5: Full gate (run the WHOLE suite, not a filtered subset)**

```bash
npx tsc --noEmit -p tsconfig.build.json          # 0
npm run lint:design                               # 0
npm test > /tmp/vitest-scoping.out 2>&1; tail -5 /tmp/vitest-scoping.out   # all green
```
If any test still references a deleted testid/symbol (e.g. a route/smoke test asserting `scoping-page` semantics), migrate its assertion to the greenfield DOM (`scoping-page` wrapper testid is preserved; rows are `scope-sample-row`). Do not leave the suite red.

- [ ] **Step 6: Commit**

```bash
git add -A src/components/AppRoutes.tsx
git commit -m "$(cat <<'EOF'
Greenfield scoping: repoint /series/new to the print page, delete legacy scoping presentation

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```
> `git rm` already staged the deletions; `git add` the modified route file. Do NOT `git add -A` the whole tree — stage only the named paths.

---

## Task 4: Migrate the scoping e2e spec to the greenfield DOM

`e2e/series-scoping.spec.ts` exists and asserts against the legacy DOM. Re-point its selectors at the greenfield page (mirroring how `e2e/series-folio.spec.ts` was migrated in Plan 4b).

**Files:**
- Modify: `e2e/series-scoping.spec.ts`
- Reference: `e2e/series-folio.spec.ts` (migrated sibling), `e2e/AGENTS.md` (route-mock + selector patterns).

- [ ] **Step 1: Read the current spec + the sibling** — note which legacy selectors it uses (`scoping-row`, legacy value inputs, the confirm modal) and map each to the greenfield equivalent (`scope-sample-row`, `scope-candidate-row`, the FlagButton value control, the direct "Confirm & build" button — there is no confirm modal now).

- [ ] **Step 2: Rewrite selectors + flow** — the greenfield flow has no inline value editing and no confirm modal: a member is confirmed by clicking its flagged value (FlagButton), and build is a single "Confirm & build" click that navigates to `/series`. Keep the route mocks (`**/api/sample-tags`, `**/api/picker-samples`, `**/api/exposures/*/trace`, `**/api/exposures/*/indices`, `**/api/sample-tags/batch`) and the identity-seed init script from `e2e/AGENTS.md` / the folio spec. Assert: rows render, a flagged value blocks build, clicking it enables build, build POSTs the batch and lands on `/series`.

- [ ] **Step 3: Run it** — `npx playwright test e2e/series-scoping.spec.ts` (auto-starts Vite) → green. If a shared `e2e/smoke.spec.ts` touches `/series/new`, update its assertion to the greenfield DOM too.

- [ ] **Step 4: Commit**

```bash
git add e2e/series-scoping.spec.ts
git commit -m "$(cat <<'EOF'
Greenfield scoping: migrate the series-scoping e2e spec to the print DOM

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Final review (after all four tasks)

- [ ] Dispatch a final whole-implementation code review (the project `frontend-reviewer` agent) over the four commits.
- [ ] Confirm the full gate once more: `npx tsc --noEmit -p tsconfig.build.json` · `npm run lint:design` · `npm test` (whole suite) · `npm run build` · `npx playwright test e2e/series-scoping.spec.ts e2e/smoke.spec.ts`.
- [ ] Do **not** run `finishing-a-development-branch` — the greenfield branch stays unmerged until the whole rebuild lands. Stop after the review + green gate and report.

## Deferred (out of scope for this plan)

- **Live pass** against `/Volumes/data` (real corpus tags, real picker projection, the N-member sparkline fan-out under real latency) — blocked on the volume mount; use the volume-free Playwright harness for now (`feedback_live_audit_harness.md`).
- **`scoping.bones.json` capture** — the page ships with the inline `SCOPING_FIXTURE`; a real boneyard capture waits on the live pass.
- **Re-grouping by an arbitrary ordering variable** (functional order-field dropdown) — needs a "force ordering key" path in `proposeOrdering`; scope decision #3 keeps the field display-only until then.
- **Keyboard-reorder a11y** for the drag list (pointer-only today; grips are `aria-hidden`).
