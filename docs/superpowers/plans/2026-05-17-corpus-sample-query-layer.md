# Corpus Sample Query Layer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `useCorpusSamples()` TanStack Query hook backed by the corpus-wide `GET /api/samples` route, with a distinct `["corpus","samples"]` query-key namespace.

**Architecture:** Pure frontend data-layer addition. A `CorpusSample` type and `listCorpusSamples` fetcher in `api.ts`; a `corpusSamples` key and `useCorpusSamples` hook in `queries.ts`; a new Vitest file. All edits are append-only and touch no existing hook, type, or test.

**Tech Stack:** TypeScript, React, TanStack Query v5, Vitest, `@testing-library/react`.

---

## Background — required reading

The corpus route `GET /api/samples` already exists (`packages/HimalayaUI/src/routes_samples.jl:27`). It returns every sample across all experiments. Each sample row carries:

- the same fields as the per-experiment `Sample` (`id`, `experiment_id`, `name`, `display_name`, `notes`, `tags`), **plus**
- an extra `q_units: string` field resolved from the owning experiment's config (defaults to `"A-1"` when the experiment is missing).

The per-experiment `Sample` interface (`api.ts:28-35`) has no `q_units`. The corpus hook must preserve it — the redesign Phase 3 reads `q_units` for cross-experiment normalization. We model this with a `CorpusSample` interface that `extends Sample`, leaving `Sample` untouched.

The design spec is `docs/superpowers/specs/2026-05-17-corpus-sample-query-layer-design.md`.

All commands below run from `packages/HimalayaUI/frontend/`.

## File structure

| File | Change | Responsibility |
|---|---|---|
| `packages/HimalayaUI/frontend/src/api.ts` | Modify | Add `CorpusSample` interface + `listCorpusSamples` fetcher |
| `packages/HimalayaUI/frontend/src/queries.ts` | Modify | Add `queryKeys.corpusSamples` key + `useCorpusSamples` hook |
| `packages/HimalayaUI/frontend/test/queries-corpus-samples.test.tsx` | Create | Vitest covering the hook fetch + key shape |

---

### Task 1: `CorpusSample` type and `listCorpusSamples` fetcher

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/api.ts` (after the `Sample` interface, line 35; and in the `// Samples` section after `listSamples`, line 105)

This task has no standalone test — the fetcher is exercised through the hook test in Task 2. It is committed separately so Task 2 starts from a compiling tree.

- [ ] **Step 1: Add the `CorpusSample` interface**

In `api.ts`, immediately after the closing `}` of the `Sample` interface (line 35), add:

```ts
// Corpus samples carry q_units (resolved from the owning experiment's
// config) — the per-experiment Sample does not. Phase 3 normalization
// reads this field. Returned by the corpus-wide GET /api/samples route.
export interface CorpusSample extends Sample {
  q_units: string;
}
```

- [ ] **Step 2: Add the `listCorpusSamples` fetcher**

In the `// Samples` section of `api.ts`, immediately after the `listSamples` declaration (it ends at line 105 with the `request<Sample[]>(...)` call), add:

```ts
export const listCorpusSamples = (): Promise<CorpusSample[]> =>
  request<CorpusSample[]>("GET", "/api/samples");
```

- [ ] **Step 3: Verify the project still type-checks**

Run: `npx tsc --noEmit`
Expected: no errors (exit code 0). The new exports are unused so far — that is fine; TypeScript does not flag unused exports.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/api.ts
git commit -m "feat: add CorpusSample type and listCorpusSamples fetcher (#156)"
```

---

### Task 2: `corpusSamples` query key and `useCorpusSamples` hook (TDD)

**Files:**
- Create: `packages/HimalayaUI/frontend/test/queries-corpus-samples.test.tsx`
- Modify: `packages/HimalayaUI/frontend/src/queries.ts` (add key after line 98 `comparisonPins`; add hook after `useSamples`, which ends at line 124)

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/frontend/test/queries-corpus-samples.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { useCorpusSamples, queryKeys } from "../src/queries";

function withClient() {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, wrapper };
}

function mockOnce(status: number, body: unknown): void {
  vi.spyOn(global, "fetch").mockResolvedValueOnce(
    new Response(JSON.stringify(body), {
      status, headers: { "Content-Type": "application/json" },
    }),
  );
}

const CORPUS_ROWS = [
  { id: 1, experiment_id: 1, name: "s1", display_name: null, notes: "",
    tags: [], q_units: "A-1" },
  { id: 2, experiment_id: 2, name: "s2", display_name: null, notes: "",
    tags: [{ id: 9, key: "lipid", value: "DOPC", source: "manual" }],
    q_units: "nm-1" },
];

describe("queries — useCorpusSamples", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("fetches GET /api/samples and returns the corpus list with q_units", async () => {
    const { wrapper } = withClient();
    mockOnce(200, CORPUS_ROWS);
    const { result } = renderHook(() => useCorpusSamples(), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));

    expect(global.fetch).toHaveBeenCalledWith(
      "/api/samples", expect.objectContaining({ method: "GET" }),
    );
    expect(result.current.data).toHaveLength(2);
    expect(result.current.data?.[0].q_units).toBe("A-1");
    expect(result.current.data?.[1].q_units).toBe("nm-1");
    expect(result.current.data?.[1].tags[0].key).toBe("lipid");
  });

  it("uses the ['corpus','samples'] key, distinct from per-experiment samples", () => {
    expect(queryKeys.corpusSamples).toEqual(["corpus", "samples"]);
    // The per-experiment key is experiment-scoped for every id — the corpus
    // key must never deep-equal it, so cache entries cannot collide.
    for (const id of [0, 1, 42]) {
      expect(queryKeys.corpusSamples).not.toEqual(queryKeys.samples(id));
    }
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npm test -- queries-corpus-samples`
Expected: FAIL — `useCorpusSamples` is not exported from `../src/queries` (and `queryKeys.corpusSamples` is undefined).

- [ ] **Step 3: Add the `corpusSamples` query key**

In `queries.ts`, inside the `queryKeys` object, immediately after the `comparisonPins` entry (line 98) and before the closing `};` (line 99), add:

```ts
  // Corpus-wide sample listing (redesign Phase 1, #156). Distinct
  // ["corpus", ...] namespace so corpus add_tag (#159) patches and
  // invalidations never collide with the per-experiment
  // ["experiment", id, "samples"] entries that `samples` produces.
  corpusSamples: ["corpus", "samples"] as const,
```

- [ ] **Step 4: Add the `useCorpusSamples` hook**

In `queries.ts`, immediately after the `useSamples` function (its closing `}` is line 124), add:

```ts
/**
 * Corpus-wide sample list — every sample across all experiments, each with
 * its tags and a q_units from the owning experiment's config. Backs the
 * Phase 1 contact sheet (#160) and sample loupe (#161). Distinct from the
 * per-experiment useSamples(experimentId): no parameter, so no `enabled`
 * gate — the fetch runs unconditionally.
 */
export function useCorpusSamples() {
  return useQuery({
    queryKey: queryKeys.corpusSamples,
    queryFn: () => api.listCorpusSamples(),
  });
}
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `npm test -- queries-corpus-samples`
Expected: PASS — both tests green.

- [ ] **Step 6: Run the build to confirm no type regressions**

Run: `npm run build`
Expected: PASS — `tsc --noEmit` clean and `vite build` succeeds.

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/frontend/src/queries.ts packages/HimalayaUI/frontend/test/queries-corpus-samples.test.tsx
git commit -m "feat: add useCorpusSamples hook and corpus query key (#156)"
```

---

## Self-review notes

- **Spec coverage:** AC1 (fetch + `q_units`) → Task 2 Step 1 test; AC2 (distinct key) → Task 2 Step 1 key-shape test; AC3 (`useSamples` unchanged) → no task touches `useSamples` or `queries-samples.test.tsx`; AC4 (Vitest) → Task 2 file; AC5 (`npm run build`) → Task 2 Step 6.
- **Type consistency:** `CorpusSample` (Task 1) is the return type of `listCorpusSamples` (Task 1) and flows through `useCorpusSamples` (Task 2). `queryKeys.corpusSamples` is defined in Task 2 Step 3 and consumed in Step 4 and the Step 1 test.
- **No placeholders:** every code step shows complete code; every run step shows the command and expected result.
