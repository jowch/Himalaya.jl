# HimalayaUI Plan 5 — Phase Panel, Miller Plot, Active-Group Overlay Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete the phase-assignment UX — render the Miller-index plot, the Phase panel (Active / Alternatives), active-group q-position overlay on the trace, and hover-preview of alternative indices.

**Architecture:** Extend the backend indices endpoint with predicted-q arrays (derived via `Himalaya.phaseratios`). Add TanStack Query hooks for groups + mutations. A client-only `hoveredIndexId` field in Zustand drives the TraceViewer's hover preview. Two new React components — `MillerPlot` and `PhasePanel` — fill the previously-placeholder right-column panes. The TraceViewer gains two additional Observable Plot marks: active-group ticks (solid, phase-colored) and hovered-index ticks (faint, phase-colored).

**Tech Stack:** Julia 1.12, Oxygen.jl (backend), `Himalaya.phaseratios` (core), React 18, TypeScript strict (`exactOptionalPropertyTypes: true`), Vite, Zustand, TanStack Query v5, Observable Plot, Vitest + RTL, Playwright.

**Assumptions about the starting state** (verified at plan-write time):
- `main` is at `8ebdbf0` — Plan 4 fully merged. 170 backend tests, 43 unit tests, 5 E2E tests all pass.
- Backend endpoints exist: `GET /api/exposures/:id/indices` (returns `[{id, exposure_id, phase, basis, score, r_squared, lattice_d, status, peaks: [{peak_id, ratio_position, residual}]}]`), `GET /api/exposures/:id/groups` (returns groups with `members: number[]`), `POST /api/groups/:id/members` and `DELETE /api/groups/:id/members/:index_id`.
- Core `Himalaya.phaseratios(P::Type; normalize=true)` returns a sparse vector of √(sum-of-squares) ratios indexed by `ratio_position`. Phases are named by `string(P)` — e.g., `"Pn3m"`, `"Im3m"`, `"Lamellar"`. The inverse is `getfield(Himalaya, Symbol(name))`.
- Frontend file layout matches Plan 4 result: `src/{api.ts,queries.ts,state.ts,App.tsx,main.tsx,ErrorBoundary.tsx,components/{Navbar,Layout,SampleList,UserModal,ExposureList,TraceViewer,StaleIndicesBanner}.tsx}`.
- Zustand store has `username`, `activeSampleId`, `activeExposureId` (persisted); TanStack Query is the provider for `useExperiment`, `useSamples`, `useExposures`, `useTrace`, `usePeaks`, `useIndices`, `useAddPeak`, `useRemovePeak`, `useReanalyzeExposure`.
- `api.ts` has types `Peak`, `IndexEntry { id, exposure_id, phase, basis, score, r_squared, lattice_d, status }` (no `peaks` field or `predicted_q` yet — Plan 5 adds them).
- E2E spec at `packages/HimalayaUI/frontend/e2e/smoke.spec.ts` has 5 tests including `clicking the trace posts a manual peak`.

**Out of scope (deferred to future plans):**
- Phase panel's **Recent** section (spec calls it "not interactive in v1" — can be a quick follow-up).
- Lamellar/Hexagonal/Square special handling for Miller plot axis label (fine to use generic "ratio" label).
- Group swapping / multiple custom groups per exposure.

**Run commands from:**
- Backend tests: repo root: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'`
- Frontend commands: `packages/HimalayaUI/frontend/`

---

## File Map

**Backend:**
- Modify: `packages/HimalayaUI/src/routes_analysis.jl` (extend `GET /api/exposures/:id/indices` with `predicted_q` and enrich `peaks[]` with `q_observed`)
- Modify: `packages/HimalayaUI/test/test_routes_analysis.jl` (assert new fields)

**Frontend:**
- Modify: `packages/HimalayaUI/frontend/src/api.ts` (enrich `IndexEntry`; add `GroupEntry` + groups fetchers)
- Modify: `packages/HimalayaUI/frontend/src/queries.ts` (add `useGroups`, `useAddIndexToGroup`, `useRemoveIndexFromGroup`; extend `queryKeys.groups`)
- Modify: `packages/HimalayaUI/frontend/src/state.ts` (add `hoveredIndexId`)
- Create: `packages/HimalayaUI/frontend/src/phases.ts` (phase → color lookup — pure data module)
- Create: `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx`
- Create: `packages/HimalayaUI/frontend/src/components/MillerPlot.tsx`
- Modify: `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx` (render active-group + hover-preview overlay marks)
- Modify: `packages/HimalayaUI/frontend/src/App.tsx` (wire `MillerPlot` into `rightTop`, `PhasePanel` into `rightBottom`, pass `activeGroupIndices` and `hoveredIndex` to TraceViewer)
- Modify: `packages/HimalayaUI/frontend/test/smoke.test.tsx` (extend mocks for new endpoints)
- Create tests: `test/{phases,queries-groups,PhasePanel,MillerPlot,TraceViewer}.test.tsx` (TraceViewer test extended in-place)
- Modify: `packages/HimalayaUI/frontend/e2e/smoke.spec.ts` (hover + confirm/exclude E2E)

---

## Task 1: Backend — extend indices endpoint with `predicted_q` + `q_observed`

**Motivation:** The frontend needs the full predicted-q series per index (for overlay ticks + Miller scatter) and the observed q for each matched peak (so the Miller plot doesn't have to cross-reference the peaks endpoint). Compute both server-side from `Himalaya.phaseratios(P; normalize=true)`.

**Files:**
- Modify: `packages/HimalayaUI/src/routes_analysis.jl`
- Modify: `packages/HimalayaUI/test/test_routes_analysis.jl`

- [ ] **Step 1: Write failing tests**

Open `packages/HimalayaUI/test/test_routes_analysis.jl`. Find the existing `GET /api/exposures/{id}/indices` testset (search for `/indices"`). Append these assertions inside the same testset, right after the existing 200 check that validates the current response shape:

```julia
# New in Plan 5: each index has predicted_q (basis × normalized phaseratios)
# and each matched peak has q_observed attached.
list = JSON3.read(String(r.body))
for entry in list
    @test haskey(entry, :predicted_q)
    @test length(entry.predicted_q) > 0
    @test all(q -> q > 0, entry.predicted_q)
    for p in entry.peaks
        @test haskey(p, :q_observed)
        @test p.q_observed > 0
    end
end
```

(The local variable `r` is the existing `HTTP.get("$base/api/exposures/$e_id/indices")` response — reuse it.)

- [ ] **Step 2: Run test to verify failure**

Run from repo root: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'`

Expected: the new assertions fail because neither `predicted_q` nor `q_observed` is present in the current response. Other tests pass.

- [ ] **Step 3: Add a helper + extend the route**

In `packages/HimalayaUI/src/routes_analysis.jl`, replace the `@get "/api/exposures/{id}/indices"` handler body with:

```julia
@get "/api/exposures/{id}/indices" function(req::HTTP.Request, id::Int)
    db = current_db()
    indices = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM indices WHERE exposure_id = ? ORDER BY score DESC", [id]))
    out = map(indices) do ix
        peak_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT ip.peak_id, ip.ratio_position, ip.residual, p.q AS q_observed
             FROM index_peaks ip JOIN peaks p ON p.id = ip.peak_id
             WHERE ip.index_id = ? ORDER BY ip.ratio_position",
            [Int(ix.id)]))
        predicted = predicted_q_for_phase(String(ix.phase), Float64(ix.basis))
        d = row_to_json(ix)
        d[:peaks]        = rows_to_json(peak_rows)
        d[:predicted_q]  = predicted
        d
    end
    HTTP.Response(200, ["Content-Type" => "application/json"],
        JSON3.write(out))
end
```

Add this helper above `register_analysis_routes!()` in the same file:

```julia
"""
    predicted_q_for_phase(phase_name, basis) -> Vector{Float64}

Return predicted q positions (basis × normalized phase ratios) for a phase
named by its string form (e.g., "Pn3m"). Returns empty vector if the phase
is unknown to Himalaya.
"""
function predicted_q_for_phase(phase_name::AbstractString, basis::Float64)::Vector{Float64}
    P = try
        getfield(Himalaya, Symbol(phase_name))
    catch
        return Float64[]
    end
    (P isa Type && P <: Himalaya.Phase) || return Float64[]
    ratios = Himalaya.phaseratios(P; normalize=true)
    [basis * r for r in ratios]
end
```

Make sure the `using Himalaya` line appears at the top of the file (it may already — search for `using Himalaya`; add it if missing).

- [ ] **Step 4: Run tests to verify pass**

Run: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'`

Expected: all 170+ tests pass, including the new assertions.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_analysis.jl \
        packages/HimalayaUI/test/test_routes_analysis.jl
git commit -m "feat(api): /api/exposures/:id/indices includes predicted_q + q_observed"
```

---

## Task 2: Frontend — enrich `IndexEntry` type; add Group types and fetchers

**Motivation:** Teach the typed client about the new index shape and the existing groups endpoints.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/api.ts`
- Modify: `packages/HimalayaUI/frontend/test/api.test.ts`

- [ ] **Step 1: Write failing tests**

Append to `packages/HimalayaUI/frontend/test/api.test.ts`, inside the existing `describe("api", () => { ... })` block:

```ts
it("listIndices returns indices with predicted_q and enriched peaks", async () => {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify([
      {
        id: 1, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 1.0,
        r_squared: 0.99, lattice_d: 12.5, status: "candidate",
        predicted_q: [0.7071, 0.866, 1.0],
        peaks: [{ peak_id: 10, ratio_position: 1, residual: 0.001, q_observed: 0.71 }],
      },
    ]), { status: 200 }),
  );
  const indices = await api.listIndices(42);
  expect(indices).toHaveLength(1);
  expect(indices[0]!.predicted_q).toEqual([0.7071, 0.866, 1.0]);
  expect(indices[0]!.peaks[0]!.q_observed).toBeCloseTo(0.71);
});

it("listGroups fetches groups for exposure", async () => {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify([
      { id: 1, exposure_id: 42, kind: "auto",   active: false, members: [10] },
      { id: 2, exposure_id: 42, kind: "custom", active: true,  members: [10, 11] },
    ]), { status: 200 }),
  );
  const groups = await api.listGroups(42);
  expect(groups).toHaveLength(2);
  expect(groups[1]!.kind).toBe("custom");
  expect(groups[1]!.active).toBe(true);
});

it("addIndexToGroup posts {index_id} with X-Username", async () => {
  const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify({
      id: 2, exposure_id: 42, kind: "custom", active: true, members: [10, 11],
    }), { status: 200 }),
  );
  const g = await api.addIndexToGroup(2, 11, { username: "alice" });
  expect(g.members).toEqual([10, 11]);
  const [, init] = fetchSpy.mock.calls[0]! as [string, RequestInit];
  expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
  expect(init.body).toBe(JSON.stringify({ index_id: 11 }));
});

it("removeIndexFromGroup sends DELETE with X-Username", async () => {
  const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify({
      id: 2, exposure_id: 42, kind: "custom", active: true, members: [10],
    }), { status: 200 }),
  );
  await api.removeIndexFromGroup(2, 11, { username: "alice" });
  const [url, init] = fetchSpy.mock.calls[0]! as [string, RequestInit];
  expect(url).toBe("/api/groups/2/members/11");
  expect(init.method).toBe("DELETE");
  expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
});
```

- [ ] **Step 2: Run to verify failure**

Run from `packages/HimalayaUI/frontend/`: `npm test -- api.test`
Expected: 4 new tests fail — `api.listGroups is not a function`, etc., plus TS errors on `predicted_q`.

- [ ] **Step 3: Extend `api.ts`**

In `packages/HimalayaUI/frontend/src/api.ts`, REPLACE the existing `IndexEntry` interface and `listIndices` export with the extended versions; and APPEND the group section.

Replace:

```ts
// Indices
export interface IndexEntry {
  id: number;
  exposure_id: number;
  phase: string;
  basis: number;
  score: number | null;
  r_squared: number | null;
  lattice_d: number | null;
  status: "candidate" | "stale";
}

export const listIndices = (exposure_id: number) =>
  request<IndexEntry[]>("GET", `/api/exposures/${exposure_id}/indices`);
```

with:

```ts
// Indices
export interface IndexPeakRef {
  peak_id: number;
  ratio_position: number;
  residual: number;
  q_observed: number;
}

export interface IndexEntry {
  id: number;
  exposure_id: number;
  phase: string;
  basis: number;
  score: number | null;
  r_squared: number | null;
  lattice_d: number | null;
  status: "candidate" | "stale";
  peaks: IndexPeakRef[];
  predicted_q: number[];
}

export const listIndices = (exposure_id: number) =>
  request<IndexEntry[]>("GET", `/api/exposures/${exposure_id}/indices`);

// Groups
export interface GroupEntry {
  id: number;
  exposure_id: number;
  kind: "auto" | "custom";
  active: boolean;
  members: number[];
}

export const listGroups = (exposure_id: number) =>
  request<GroupEntry[]>("GET", `/api/exposures/${exposure_id}/groups`);
export const addIndexToGroup = (group_id: number, index_id: number, opts?: AuthOpts) =>
  request<GroupEntry>("POST", `/api/groups/${group_id}/members`, { index_id }, opts);
export const removeIndexFromGroup = (group_id: number, index_id: number, opts?: AuthOpts) =>
  request<GroupEntry>("DELETE", `/api/groups/${group_id}/members/${index_id}`, undefined, opts);
```

(Keep the existing `reanalyzeExposure` export at its current position.)

- [ ] **Step 4: Run to verify pass**

Run: `npm test -- api.test && npm run build`
Expected: all api tests pass, tsc clean.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/api.ts \
        packages/HimalayaUI/frontend/test/api.test.ts
git commit -m "feat(ui): index/group types + group fetchers in api.ts"
```

---

## Task 3: Frontend — Query hooks for groups + mutations

**Motivation:** Surface groups through TanStack Query with correct cache invalidation (mutations touch groups, which downstream affects the trace overlay).

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/queries.ts`
- Create: `packages/HimalayaUI/frontend/test/queries-groups.test.tsx`

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/queries-groups.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor, act } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useGroups, useAddIndexToGroup, useRemoveIndexFromGroup, queryKeys,
} from "../src/queries";

function withClient() {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, wrapper };
}

function mockOnce(status: number, body: unknown): void {
  vi.spyOn(global, "fetch").mockResolvedValueOnce(
    new Response(status === 204 ? null : JSON.stringify(body), {
      status, headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("queries — groups", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("useGroups fetches when exposureId is provided", async () => {
    mockOnce(200, [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [] }]);
    const { wrapper } = withClient();
    const { result } = renderHook(() => useGroups(42), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(result.current.data).toHaveLength(1);
  });

  it("useGroups disabled when exposureId undefined", () => {
    const spy = vi.spyOn(global, "fetch");
    const { wrapper } = withClient();
    const { result } = renderHook(() => useGroups(undefined), { wrapper });
    expect(result.current.fetchStatus).toBe("idle");
    expect(spy).not.toHaveBeenCalled();
  });

  it("useAddIndexToGroup invalidates groups for exposure", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(200, { id: 2, exposure_id: 42, kind: "custom", active: true, members: [10, 11] });
    const { result } = renderHook(() => useAddIndexToGroup(42, 2), { wrapper });
    await act(async () => { await result.current.mutateAsync(11); });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.groups(42) });
  });

  it("useRemoveIndexFromGroup invalidates groups for exposure", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(200, { id: 2, exposure_id: 42, kind: "custom", active: true, members: [10] });
    const { result } = renderHook(() => useRemoveIndexFromGroup(42, 2), { wrapper });
    await act(async () => { await result.current.mutateAsync(11); });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.groups(42) });
  });
});
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- queries-groups`
Expected: 4 tests fail.

- [ ] **Step 3: Extend `queries.ts`**

In `packages/HimalayaUI/frontend/src/queries.ts`, extend the `queryKeys` object by adding one line (after `indices`):

```ts
  groups:     (exposureId: number) => ["exposure", exposureId, "groups"] as const,
```

Then append these exports (below the existing mutation hooks). Reuse the existing `authOpts` helper that was introduced in Plan 4 Task 3.

```ts
export function useGroups(exposureId: number | undefined) {
  return useQuery({
    queryKey: ["exposure", exposureId ?? "none", "groups"] as const,
    queryFn: () => api.listGroups(exposureId as number),
    enabled: exposureId !== undefined,
  });
}

export function useAddIndexToGroup(exposureId: number, groupId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (indexId: number) =>
      api.addIndexToGroup(groupId, indexId, authOpts(username)),
    onSuccess: () => qc.invalidateQueries({ queryKey: queryKeys.groups(exposureId) }),
  });
}

export function useRemoveIndexFromGroup(exposureId: number, groupId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (indexId: number) =>
      api.removeIndexFromGroup(groupId, indexId, authOpts(username)),
    onSuccess: () => qc.invalidateQueries({ queryKey: queryKeys.groups(exposureId) }),
  });
}
```

- [ ] **Step 4: Run to verify pass**

Run: `npm test && npm run build`
Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/queries.ts \
        packages/HimalayaUI/frontend/test/queries-groups.test.tsx
git commit -m "feat(ui): useGroups, useAddIndexToGroup, useRemoveIndexFromGroup hooks"
```

---

## Task 4: Frontend — `hoveredIndexId` in Zustand

**Motivation:** The hover preview needs a single shared piece of UI state the PhasePanel writes and the TraceViewer reads. Zustand is the client-state home; do NOT persist `hoveredIndexId`.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/state.ts`
- Modify: `packages/HimalayaUI/frontend/test/state.test.ts`

- [ ] **Step 1: Write failing tests**

Append to `packages/HimalayaUI/frontend/test/state.test.ts`, inside the existing `describe` block:

```ts
it("hoveredIndexId starts undefined and can be set/cleared", () => {
  useAppState.setState({ hoveredIndexId: undefined });
  expect(useAppState.getState().hoveredIndexId).toBeUndefined();
  useAppState.getState().setHoveredIndex(7);
  expect(useAppState.getState().hoveredIndexId).toBe(7);
  useAppState.getState().setHoveredIndex(undefined);
  expect(useAppState.getState().hoveredIndexId).toBeUndefined();
});

it("hoveredIndexId is NOT in the persisted partition", () => {
  useAppState.setState({ hoveredIndexId: 42 });
  const raw = localStorage.getItem("himalaya-ui:state");
  // After any write, the persist middleware has flushed. Confirm hoveredIndexId
  // does not appear in the stringified storage payload.
  expect(raw ?? "").not.toContain("hoveredIndexId");
});
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- state.test`
Expected: `setHoveredIndex` not defined.

- [ ] **Step 3: Extend `state.ts`**

In `packages/HimalayaUI/frontend/src/state.ts`:

- Add `hoveredIndexId: number | undefined;` to the `AppState` interface.
- Add `setHoveredIndex: (id: number | undefined) => void;` to `AppState`.
- In `create<AppState>()(persist(...))` initializer, add `hoveredIndexId: undefined,` and `setHoveredIndex: (hoveredIndexId) => set({ hoveredIndexId }),`.
- **Important**: leave the `partialize` argument unchanged — it should still only include `{ username, activeSampleId, activeExposureId }`. That keeps `hoveredIndexId` in memory only, never in localStorage.

Full replacement file:

```ts
import { create } from "zustand";
import { persist } from "zustand/middleware";

export const LS_KEY = "himalaya-ui:state";

export interface AppState {
  username: string | undefined;
  activeSampleId: number | undefined;
  activeExposureId: number | undefined;
  hoveredIndexId: number | undefined;
  setUsername: (name: string) => void;
  setActiveSample: (id: number | undefined) => void;
  setActiveExposure: (id: number | undefined) => void;
  setHoveredIndex: (id: number | undefined) => void;
}

export const useAppState = create<AppState>()(
  persist(
    (set) => ({
      username: undefined,
      activeSampleId: undefined,
      activeExposureId: undefined,
      hoveredIndexId: undefined,
      setUsername: (username) => set({ username }),
      setActiveSample: (activeSampleId) =>
        set({ activeSampleId, activeExposureId: undefined }),
      setActiveExposure: (activeExposureId) => set({ activeExposureId }),
      setHoveredIndex: (hoveredIndexId) => set({ hoveredIndexId }),
    }),
    {
      name: LS_KEY,
      version: 1,
      partialize: (s) => ({
        username: s.username,
        activeSampleId: s.activeSampleId,
        activeExposureId: s.activeExposureId,
      }),
    },
  ),
);
```

- [ ] **Step 4: Run to verify pass**

Run: `npm test && npm run build`
Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/state.ts \
        packages/HimalayaUI/frontend/test/state.test.ts
git commit -m "feat(ui): hoveredIndexId in Zustand (client-only, non-persisted)"
```

---

## Task 5: Phase → color lookup module

**Motivation:** Phases need a stable color mapping shared by the TraceViewer overlay, PhasePanel swatches, and MillerPlot series colors. One-line data file.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/phases.ts`
- Create: `packages/HimalayaUI/frontend/test/phases.test.ts`

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/phases.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import { phaseColor, KNOWN_PHASES } from "../src/phases";

describe("phases", () => {
  it("returns a distinct color for each known phase", () => {
    const colors = new Set(KNOWN_PHASES.map((p) => phaseColor(p)));
    expect(colors.size).toBe(KNOWN_PHASES.length);
    for (const c of colors) expect(c).toMatch(/^#[0-9a-f]{6}$/i);
  });

  it("returns a fallback color for unknown phases", () => {
    expect(phaseColor("Unknown")).toMatch(/^#[0-9a-f]{6}$/i);
  });
});
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- phases`

- [ ] **Step 3: Implement `phases.ts`**

Create `packages/HimalayaUI/frontend/src/phases.ts`:

```ts
export const KNOWN_PHASES = [
  "Pn3m", "Im3m", "Ia3d", "Fm3m", "Fd3m",
  "Hexagonal", "Lamellar", "Square",
] as const;

export type KnownPhase = typeof KNOWN_PHASES[number];

const PALETTE: Record<KnownPhase, string> = {
  Pn3m:      "#d97706",  // accent orange
  Im3m:      "#22c55e",  // green
  Ia3d:      "#3b82f6",  // blue
  Fm3m:      "#a855f7",  // purple
  Fd3m:      "#ec4899",  // pink
  Hexagonal: "#eab308",  // yellow
  Lamellar:  "#14b8a6",  // teal
  Square:    "#f97316",  // different orange
};

const FALLBACK = "#9a9894"; // matches --color-fg-muted

export function phaseColor(phase: string): string {
  return (PALETTE as Record<string, string>)[phase] ?? FALLBACK;
}
```

- [ ] **Step 4: Run to verify pass**

Run: `npm test -- phases`
Expected: 2/2 pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/phases.ts \
        packages/HimalayaUI/frontend/test/phases.test.ts
git commit -m "feat(ui): phaseColor lookup for Miller/phase/overlay coloring"
```

---

## Task 6: PhasePanel — active group section

**Motivation:** Start with the simpler half of the panel — show each index currently in the active group (phase, lattice parameter, R²) with a `−` button. Alternatives + hover come in Task 7.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx`
- Create: `packages/HimalayaUI/frontend/test/PhasePanel.test.tsx`

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/PhasePanel.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, waitFor, fireEvent } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { PhasePanel } from "../src/components/PhasePanel";

beforeEach(() => { vi.restoreAllMocks(); });

function mockAll(indices: unknown[], groups: unknown[]): void {
  vi.spyOn(global, "fetch").mockImplementation(async (input) => {
    const u = typeof input === "string" ? input : (input as Request).url;
    if (u.endsWith("/indices")) {
      return new Response(JSON.stringify(indices),
        { status: 200, headers: { "Content-Type": "application/json" } });
    }
    if (u.endsWith("/groups")) {
      return new Response(JSON.stringify(groups),
        { status: 200, headers: { "Content-Type": "application/json" } });
    }
    if (u.match(/\/api\/groups\/\d+\/members\/\d+$/)) {
      return new Response(JSON.stringify({
        id: 2, exposure_id: 42, kind: "custom", active: true, members: [],
      }), { status: 200, headers: { "Content-Type": "application/json" } });
    }
    return new Response("not found", { status: 404 });
  });
}

describe("<PhasePanel> — active group", () => {
  it("renders a hint when no exposure is active", () => {
    renderWithProviders(<PhasePanel exposureId={undefined} />);
    expect(screen.getByText(/no exposure selected/i)).toBeInTheDocument();
  });

  it("renders each active-group index with phase, lattice, R² and a remove button", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 1.0,
          r_squared: 0.998, lattice_d: 12.5, status: "candidate",
          predicted_q: [0.7, 0.9], peaks: [] },
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
          r_squared: 0.71, lattice_d: 9.1, status: "candidate",
          predicted_q: [0.4, 0.6], peaks: [] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    await waitFor(() => expect(screen.getByText("Pn3m")).toBeInTheDocument());
    expect(screen.getByText(/12\.50/)).toBeInTheDocument();
    expect(screen.getByText(/0\.998/)).toBeInTheDocument();
    // The active-group row has a remove button
    const rm = screen.getByRole("button", { name: /remove index 10/i });
    expect(rm).toBeInTheDocument();
  });

  it("clicking remove calls DELETE /api/groups/:gid/members/:indexId", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 1.0,
          r_squared: 0.998, lattice_d: 12.5, status: "candidate",
          predicted_q: [0.7, 0.9], peaks: [] },
      ],
      [{ id: 2, exposure_id: 42, kind: "custom", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const rm = await screen.findByRole("button", { name: /remove index 10/i });
    fireEvent.click(rm);
    await waitFor(() => {
      const spy = global.fetch as unknown as { mock: { calls: unknown[][] } };
      const urls = spy.mock.calls.map((c) =>
        typeof c[0] === "string" ? c[0] : (c[0] as Request).url);
      expect(urls.some((u) => u === "/api/groups/2/members/10")).toBe(true);
    });
  });
});
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- PhasePanel`
Expected: tests fail (`PhasePanel` not defined).

- [ ] **Step 3: Implement `PhasePanel.tsx` (active-group section only)**

Create `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx`:

```tsx
import { useIndices, useGroups, useRemoveIndexFromGroup } from "../queries";
import { phaseColor } from "../phases";
import type { IndexEntry, GroupEntry } from "../api";

export interface PhasePanelProps {
  exposureId: number | undefined;
}

function activeGroup(groups: GroupEntry[]): GroupEntry | undefined {
  return groups.find((g) => g.active);
}

function formatLattice(d: number | null): string {
  return d != null ? d.toFixed(2) : "—";
}

function formatR2(r: number | null): string {
  return r != null ? r.toFixed(3) : "—";
}

export function PhasePanel({ exposureId }: PhasePanelProps): JSX.Element {
  const indicesQ = useIndices(exposureId);
  const groupsQ  = useGroups(exposureId);
  const active   = (groupsQ.data && activeGroup(groupsQ.data)) ?? undefined;
  const removeMember = useRemoveIndexFromGroup(exposureId ?? 0, active?.id ?? 0);

  if (exposureId === undefined) {
    return <p className="text-fg-muted italic">No exposure selected.</p>;
  }
  if (indicesQ.isPending || groupsQ.isPending) {
    return <p className="text-fg-muted">Loading phase assignments…</p>;
  }

  const indices = indicesQ.data ?? [];
  const activeMembers = (active?.members ?? [])
    .map((id) => indices.find((i) => i.id === id))
    .filter((i): i is IndexEntry => i != null);

  return (
    <div className="flex flex-col gap-3">
      <section>
        <h3 className="text-fg-muted text-[12px] uppercase tracking-wide mb-1">
          Active group
        </h3>
        {activeMembers.length === 0 ? (
          <p className="text-fg-muted italic text-[13px]">No indices in the active group.</p>
        ) : (
          <ul className="flex flex-col gap-1">
            {activeMembers.map((ix) => (
              <li
                key={ix.id}
                data-index-id={ix.id}
                className="flex items-center gap-2 px-2 py-1 rounded-md bg-bg-elevated"
              >
                <span
                  className="w-2 h-2 rounded-full"
                  style={{ background: phaseColor(ix.phase) }}
                  aria-hidden
                />
                <span className="font-medium">{ix.phase}</span>
                <span className="text-fg-muted text-[13px]">
                  a={formatLattice(ix.lattice_d)} nm · R²={formatR2(ix.r_squared)}
                </span>
                <button
                  className="ml-auto text-fg-muted hover:text-error focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent rounded-md px-1"
                  aria-label={`Remove index ${ix.id}`}
                  onClick={() => {
                    if (active) removeMember.mutate(ix.id);
                  }}
                >
                  −
                </button>
              </li>
            ))}
          </ul>
        )}
      </section>
    </div>
  );
}
```

- [ ] **Step 4: Run to verify pass**

Run: `npm test -- PhasePanel`
Expected: 3/3 pass.

- [ ] **Step 5: Full suite + build**

Run: `npm test && npm run build`

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/PhasePanel.tsx \
        packages/HimalayaUI/frontend/test/PhasePanel.test.tsx
git commit -m "feat(ui): PhasePanel active-group section with remove button"
```

---

## Task 7: PhasePanel — alternatives section + hover wiring

**Motivation:** Show candidate indices not in the active group, ranked by score. `+` button adds to group. Hover sets/clears `hoveredIndexId` in Zustand, which the TraceViewer will consume in Task 9.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx`
- Modify: `packages/HimalayaUI/frontend/test/PhasePanel.test.tsx`

- [ ] **Step 1: Write failing tests**

Append to `packages/HimalayaUI/frontend/test/PhasePanel.test.tsx` (inside the top-level `describe("<PhasePanel> — active group"` block, or open a new `describe` below it — either works):

```tsx
import { useAppState } from "../src/state";

describe("<PhasePanel> — alternatives", () => {
  it("renders alternative indices with a + button", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 1.0,
          r_squared: 0.998, lattice_d: 12.5, status: "candidate",
          predicted_q: [0.7, 0.9], peaks: [] },
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
          r_squared: 0.71, lattice_d: 9.1, status: "candidate",
          predicted_q: [0.4, 0.6], peaks: [] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    await waitFor(() => expect(screen.getByText("Im3m")).toBeInTheDocument());
    expect(screen.getByRole("button", { name: /add index 11/i })).toBeInTheDocument();
  });

  it("hovering an alternative sets hoveredIndexId; leaving clears it", async () => {
    useAppState.setState({ hoveredIndexId: undefined });
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 1.0,
          r_squared: 0.998, lattice_d: 12.5, status: "candidate",
          predicted_q: [0.7], peaks: [] },
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
          r_squared: 0.71, lattice_d: 9.1, status: "candidate",
          predicted_q: [0.4], peaks: [] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const row = await waitFor(() =>
      document.querySelector<HTMLElement>('[data-alternative-id="11"]') ?? (() => { throw new Error("not found"); })(),
    );
    fireEvent.mouseEnter(row);
    expect(useAppState.getState().hoveredIndexId).toBe(11);
    fireEvent.mouseLeave(row);
    expect(useAppState.getState().hoveredIndexId).toBeUndefined();
  });

  it("clicking + on an alternative posts to /api/groups/:gid/members", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 1.0,
          r_squared: 0.998, lattice_d: 12.5, status: "candidate",
          predicted_q: [0.7], peaks: [] },
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
          r_squared: 0.71, lattice_d: 9.1, status: "candidate",
          predicted_q: [0.4], peaks: [] },
      ],
      [{ id: 2, exposure_id: 42, kind: "custom", active: true, members: [10] }],
    );
    // Intercept POST too
    (global.fetch as unknown as { mockImplementation: (fn: (i: RequestInfo) => Promise<Response>) => void })
      .mockImplementation(async (input: RequestInfo | URL) => {
        const u = typeof input === "string" ? input : (input as Request).url;
        if (u.endsWith("/indices")) {
          return new Response(JSON.stringify([
            { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 1.0,
              r_squared: 0.998, lattice_d: 12.5, status: "candidate",
              predicted_q: [0.7], peaks: [] },
            { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
              r_squared: 0.71, lattice_d: 9.1, status: "candidate",
              predicted_q: [0.4], peaks: [] },
          ]), { status: 200, headers: { "Content-Type": "application/json" } });
        }
        if (u.endsWith("/groups")) {
          return new Response(JSON.stringify([
            { id: 2, exposure_id: 42, kind: "custom", active: true, members: [10] },
          ]), { status: 200, headers: { "Content-Type": "application/json" } });
        }
        if (u.endsWith("/api/groups/2/members")) {
          return new Response(JSON.stringify({
            id: 2, exposure_id: 42, kind: "custom", active: true, members: [10, 11],
          }), { status: 200, headers: { "Content-Type": "application/json" } });
        }
        return new Response("not found", { status: 404 });
      });

    renderWithProviders(<PhasePanel exposureId={42} />);
    const addBtn = await screen.findByRole("button", { name: /add index 11/i });
    fireEvent.click(addBtn);
    await waitFor(() => {
      const spy = global.fetch as unknown as { mock: { calls: unknown[][] } };
      const urls = spy.mock.calls.map((c) =>
        typeof c[0] === "string" ? c[0] : (c[0] as Request).url);
      expect(urls).toContain("/api/groups/2/members");
    });
  });
});
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- PhasePanel`
Expected: 3 new tests fail (no alternatives section yet).

- [ ] **Step 3: Extend `PhasePanel.tsx`**

Replace the whole file with the expanded version:

```tsx
import { useIndices, useGroups, useAddIndexToGroup, useRemoveIndexFromGroup } from "../queries";
import { useAppState } from "../state";
import { phaseColor } from "../phases";
import type { IndexEntry, GroupEntry } from "../api";

export interface PhasePanelProps {
  exposureId: number | undefined;
}

function activeGroup(groups: GroupEntry[]): GroupEntry | undefined {
  return groups.find((g) => g.active);
}

function formatLattice(d: number | null): string {
  return d != null ? d.toFixed(2) : "—";
}

function formatR2(r: number | null): string {
  return r != null ? r.toFixed(3) : "—";
}

export function PhasePanel({ exposureId }: PhasePanelProps): JSX.Element {
  const indicesQ = useIndices(exposureId);
  const groupsQ  = useGroups(exposureId);
  const setHoveredIndex = useAppState((s) => s.setHoveredIndex);
  const active = (groupsQ.data && activeGroup(groupsQ.data)) ?? undefined;
  const addMember    = useAddIndexToGroup(exposureId ?? 0, active?.id ?? 0);
  const removeMember = useRemoveIndexFromGroup(exposureId ?? 0, active?.id ?? 0);

  if (exposureId === undefined) {
    return <p className="text-fg-muted italic">No exposure selected.</p>;
  }
  if (indicesQ.isPending || groupsQ.isPending) {
    return <p className="text-fg-muted">Loading phase assignments…</p>;
  }

  const indices = (indicesQ.data ?? []).slice().sort(
    (a, b) => (b.score ?? 0) - (a.score ?? 0),
  );
  const memberIds = new Set(active?.members ?? []);
  const activeMembers = indices.filter((ix) => memberIds.has(ix.id));
  const alternatives  = indices.filter((ix) => !memberIds.has(ix.id));

  return (
    <div className="flex flex-col gap-3">
      <section>
        <h3 className="text-fg-muted text-[12px] uppercase tracking-wide mb-1">Active group</h3>
        {activeMembers.length === 0 ? (
          <p className="text-fg-muted italic text-[13px]">No indices in the active group.</p>
        ) : (
          <ul className="flex flex-col gap-1">
            {activeMembers.map((ix) => (
              <li
                key={ix.id}
                data-index-id={ix.id}
                className="flex items-center gap-2 px-2 py-1 rounded-md bg-bg-elevated"
              >
                <span
                  className="w-2 h-2 rounded-full"
                  style={{ background: phaseColor(ix.phase) }}
                  aria-hidden
                />
                <span className="font-medium">{ix.phase}</span>
                <span className="text-fg-muted text-[13px]">
                  a={formatLattice(ix.lattice_d)} nm · R²={formatR2(ix.r_squared)}
                </span>
                <button
                  className="ml-auto text-fg-muted hover:text-error focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent rounded-md px-1"
                  aria-label={`Remove index ${ix.id}`}
                  onClick={() => { if (active) removeMember.mutate(ix.id); }}
                >
                  −
                </button>
              </li>
            ))}
          </ul>
        )}
      </section>

      <section>
        <h3 className="text-fg-muted text-[12px] uppercase tracking-wide mb-1">
          Alternatives ({alternatives.length})
        </h3>
        {alternatives.length === 0 ? (
          <p className="text-fg-muted italic text-[13px]">No alternatives.</p>
        ) : (
          <ul className="flex flex-col gap-1">
            {alternatives.map((ix) => (
              <li
                key={ix.id}
                data-alternative-id={ix.id}
                className="flex items-center gap-2 px-2 py-1 rounded-md hover:bg-bg-hover cursor-default"
                onMouseEnter={() => setHoveredIndex(ix.id)}
                onMouseLeave={() => setHoveredIndex(undefined)}
              >
                <span
                  className="w-2 h-2 rounded-full"
                  style={{ background: phaseColor(ix.phase) }}
                  aria-hidden
                />
                <span className="font-medium">{ix.phase}</span>
                <span className="text-fg-muted text-[13px]">
                  a={formatLattice(ix.lattice_d)} nm · R²={formatR2(ix.r_squared)}
                </span>
                <button
                  className="ml-auto text-fg-muted hover:text-success focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent rounded-md px-1"
                  aria-label={`Add index ${ix.id}`}
                  onClick={() => { if (active) addMember.mutate(ix.id); }}
                  disabled={active === undefined}
                >
                  +
                </button>
              </li>
            ))}
          </ul>
        )}
      </section>
    </div>
  );
}
```

- [ ] **Step 4: Run to verify pass**

Run: `npm test -- PhasePanel`
Expected: all 6 tests pass (3 original + 3 new).

- [ ] **Step 5: Full suite + build**

Run: `npm test && npm run build`

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/PhasePanel.tsx \
        packages/HimalayaUI/frontend/test/PhasePanel.test.tsx
git commit -m "feat(ui): PhasePanel alternatives section with hover + add"
```

---

## Task 8: TraceViewer — active-group overlay ticks

**Motivation:** Draw vertical tick marks at each predicted q of every index in the active group, colored by phase. This gives the user instant visual feedback: do the predicted positions match observed peaks?

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx`
- Modify: `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx`

- [ ] **Step 1: Write failing tests**

Append to `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx`:

```tsx
import type { IndexEntry } from "../src/api";

describe("<TraceViewer> — overlays", () => {
  it("calls Plot.ruleX with predicted_q when activeGroupIndices is non-empty", async () => {
    const Plot = await import("@observablehq/plot");
    const ruleX = vi.fn(() => ({ _kind: "ruleX" }));
    // Augment the module mock — add ruleX if not present
    (Plot as unknown as { ruleX: typeof ruleX }).ruleX = ruleX;

    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 30], sigma: [1, 1, 1] };
    const indices: IndexEntry[] = [{
      id: 1, exposure_id: 1, phase: "Pn3m", basis: 0.5, score: 1,
      r_squared: 0.99, lattice_d: 12, status: "candidate",
      predicted_q: [0.7, 0.9],
      peaks: [],
    }];
    render(
      <TraceViewer
        trace={trace}
        peaks={[]}
        activeGroupIndices={indices}
        hoveredIndex={undefined}
        onAddPeak={() => {}}
        onRemovePeak={() => {}}
      />,
    );
    expect(ruleX).toHaveBeenCalled();
  });
});
```

Also update existing TraceViewer test renders (`render(<TraceViewer ... />)` calls) to pass the new required props:

```tsx
activeGroupIndices={[]}
hoveredIndex={undefined}
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- TraceViewer`
Expected: new test fails (and possibly old tests fail with TS errors until props added).

- [ ] **Step 3: Extend `TraceViewer.tsx`**

Replace `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx`:

```tsx
import { useEffect, useRef } from "react";
import * as Plot from "@observablehq/plot";
import type { Trace, Peak, IndexEntry } from "../api";
import { phaseColor } from "../phases";

export interface TraceViewerProps {
  trace: Trace;
  peaks: Peak[];
  activeGroupIndices: IndexEntry[];
  hoveredIndex: IndexEntry | undefined;
  onAddPeak: (q: number) => void;
  onRemovePeak: (peakId: number) => void;
}

export function snapToLocalMax(q: number[], I: number[], qClick: number, K = 3): number {
  let nearest = 0;
  for (let i = 1; i < q.length; i++) {
    if (Math.abs(q[i]! - qClick) < Math.abs(q[nearest]! - qClick)) nearest = i;
  }
  const lo = Math.max(0, nearest - K);
  const hi = Math.min(q.length - 1, nearest + K);
  let best = lo;
  for (let i = lo + 1; i <= hi; i++) {
    if (I[i]! > I[best]!) best = i;
  }
  return q[best]!;
}

export function findNearestPeak(
  peaks: Peak[], qClick: number, tolerance: number,
): Peak | null {
  let best: Peak | null = null;
  let bestDist = Infinity;
  for (const p of peaks) {
    const d = Math.abs(p.q - qClick);
    if (d < bestDist) { best = p; bestDist = d; }
  }
  return best && bestDist <= tolerance ? best : null;
}

function withIntensity(peaks: Peak[], trace: Trace): Array<{ q: number; I: number }> {
  return peaks.map((p) => ({ q: p.q, I: interpolateI(p.q, trace) }));
}

function interpolateI(q: number, trace: Trace): number {
  let nearest = 0;
  for (let i = 1; i < trace.q.length; i++) {
    if (Math.abs(trace.q[i]! - q) < Math.abs(trace.q[nearest]! - q)) nearest = i;
  }
  return trace.I[nearest]!;
}

function indexOverlayData(indices: IndexEntry[]): Array<{ q: number; phase: string; color: string }> {
  return indices.flatMap((ix) =>
    ix.predicted_q.map((q) => ({ q, phase: ix.phase, color: phaseColor(ix.phase) })),
  );
}

export function TraceViewer({
  trace, peaks, activeGroupIndices, hoveredIndex, onAddPeak, onRemovePeak,
}: TraceViewerProps): JSX.Element {
  const hostRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const host = hostRef.current;
    if (!host) return;

    const data = trace.q.map((q, i) => ({
      q, I: trace.I[i]!,
      lo: Math.max(1e-12, trace.I[i]! - trace.sigma[i]!),
      hi: trace.I[i]! + trace.sigma[i]!,
    }));

    const activeOverlay  = indexOverlayData(activeGroupIndices);
    const hoveredOverlay = hoveredIndex ? indexOverlayData([hoveredIndex]) : [];

    const el = Plot.plot({
      width: host.clientWidth || 400,
      height: host.clientHeight || 300,
      marginLeft: 50, marginBottom: 40,
      x: { type: "log", label: "q (nm⁻¹)" },
      y: { type: "log", label: "I (a.u.)" },
      marks: [
        Plot.areaY(data, { x: "q", y1: "lo", y2: "hi",
          fill: "var(--color-accent)", fillOpacity: 0.15 }),
        Plot.line(data, { x: "q", y: "I",
          stroke: "var(--color-fg)", strokeWidth: 1 }),
        Plot.ruleX(activeOverlay, {
          x: "q",
          stroke: (d: { color: string }) => d.color,
          strokeWidth: 1.5,
          strokeOpacity: 0.9,
        }),
        Plot.ruleX(hoveredOverlay, {
          x: "q",
          stroke: (d: { color: string }) => d.color,
          strokeWidth: 1.5,
          strokeOpacity: 0.45,
          strokeDasharray: "3,3",
        }),
        Plot.dot(withIntensity(peaks.filter((p) => p.source === "auto"), trace),
          { x: "q", y: "I", fill: "var(--color-accent)", r: 4 }),
        Plot.dot(withIntensity(peaks.filter((p) => p.source === "manual"), trace),
          { x: "q", y: "I", stroke: "var(--color-warning)", strokeWidth: 2, fill: "none", r: 5 }),
      ],
    });

    host.replaceChildren(el);

    function handleClick(ev: MouseEvent): void {
      const xScale = (el as unknown as { scale: (name: string) => { invert?: (v: number) => number } | undefined }).scale("x");
      if (!xScale?.invert) return;
      const rect = el.getBoundingClientRect();
      const qClick = xScale.invert(ev.clientX - rect.left);
      const tolerance = Math.max(qClick * 0.02, 1e-6);
      const existing = findNearestPeak(peaks, qClick, tolerance);
      if (existing) onRemovePeak(existing.id);
      else          onAddPeak(snapToLocalMax(trace.q, trace.I, qClick));
    }
    (el as unknown as EventTarget).addEventListener("click", handleClick);

    return () => {
      (el as unknown as EventTarget).removeEventListener("click", handleClick);
      host.replaceChildren();
    };
  }, [trace, peaks, activeGroupIndices, hoveredIndex, onAddPeak, onRemovePeak]);

  return <div ref={hostRef} className="w-full h-full" data-testid="trace-viewer" />;
}
```

Also update the mocked `@observablehq/plot` module at the top of `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx` to include `ruleX`:

```tsx
vi.mock("@observablehq/plot", () => ({
  plot: vi.fn(() => {
    const el = document.createElement("div");
    el.setAttribute("data-testid", "plot-svg");
    return el;
  }),
  areaY: vi.fn(() => ({ _kind: "areaY" })),
  line:  vi.fn(() => ({ _kind: "line" })),
  dot:   vi.fn(() => ({ _kind: "dot" })),
  ruleX: vi.fn(() => ({ _kind: "ruleX" })),
}));
```

- [ ] **Step 4: Run to verify pass**

Run: `npm test -- TraceViewer`
Expected: all tests pass.

- [ ] **Step 5: Full suite + build**

Run: `npm test && npm run build`

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/TraceViewer.tsx \
        packages/HimalayaUI/frontend/test/TraceViewer.test.tsx
git commit -m "feat(ui): TraceViewer renders active-group overlay (ruleX)"
```

---

## Task 9: TraceViewer hover preview — already wired via `hoveredIndex` prop

Task 8's implementation already renders the hovered-index's predicted_q via the `hoveredOverlay` Plot.ruleX with dashed stroke. No separate code change needed — but App.tsx needs to compute the hovered IndexEntry from Zustand's `hoveredIndexId` and pass it through. That wiring happens in Task 11. For this task, lock in the overlay-data helper with a unit test.

**Files:**
- Modify: `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx` (new test for hover branch)

- [ ] **Step 1: Write failing test**

Append to `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx` inside the `describe("<TraceViewer> — overlays"` block:

```tsx
it("calls Plot.ruleX a second time (with dashed opts) when hoveredIndex is set", async () => {
  const Plot = await import("@observablehq/plot");
  const ruleX = Plot.ruleX as unknown as ReturnType<typeof vi.fn>;
  ruleX.mockClear();

  const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 30], sigma: [1, 1, 1] };
  const hovered: IndexEntry = {
    id: 11, exposure_id: 1, phase: "Im3m", basis: 0.3, score: 0.5,
    r_squared: 0.7, lattice_d: 9.0, status: "candidate",
    predicted_q: [0.42, 0.6],
    peaks: [],
  };
  render(
    <TraceViewer
      trace={trace}
      peaks={[]}
      activeGroupIndices={[]}
      hoveredIndex={hovered}
      onAddPeak={() => {}}
      onRemovePeak={() => {}}
    />,
  );
  // Two ruleX calls: one for activeGroupIndices (empty -> still called), one for hovered.
  expect(ruleX.mock.calls.length).toBe(2);
  // Second call's data has 2 rows (hovered.predicted_q.length)
  const [data] = ruleX.mock.calls[1]! as [Array<{ q: number }>];
  expect(data.length).toBe(2);
});
```

- [ ] **Step 2: Run**

Run: `npm test -- TraceViewer`
Expected: PASS. (The production code already supports this; the test just locks it in.)

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/test/TraceViewer.test.tsx
git commit -m "test(ui): lock in TraceViewer hover-preview overlay"
```

---

## Task 10: MillerPlot component

**Motivation:** Show a scatter of observed q vs predicted ratio, with a linear fit line per phase. Slope is the lattice parameter d, and R² is annotated. Uses Observable Plot's `linearRegressionY`.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/MillerPlot.tsx`
- Create: `packages/HimalayaUI/frontend/test/MillerPlot.test.tsx`

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/MillerPlot.test.tsx`:

```tsx
import { describe, it, expect, vi } from "vitest";
import { render } from "@testing-library/react";
import type { IndexEntry } from "../src/api";
import { MillerPlot, toScatterData } from "../src/components/MillerPlot";

vi.mock("@observablehq/plot", () => ({
  plot: vi.fn(() => {
    const el = document.createElement("div");
    el.setAttribute("data-testid", "miller-svg");
    return el;
  }),
  dot:   vi.fn(() => ({ _kind: "dot" })),
  linearRegressionY: vi.fn(() => ({ _kind: "linreg" })),
  text:  vi.fn(() => ({ _kind: "text" })),
}));

describe("toScatterData", () => {
  it("maps each peak to (ratio, q_observed, phase)", () => {
    const ix: IndexEntry = {
      id: 1, exposure_id: 1, phase: "Pn3m", basis: 0.5, score: 1,
      r_squared: 0.99, lattice_d: 12, status: "candidate",
      predicted_q: [0.7071, 0.866, 1.0],
      peaks: [
        { peak_id: 10, ratio_position: 1, residual: 0.001, q_observed: 0.71 },
        { peak_id: 11, ratio_position: 3, residual: 0.002, q_observed: 1.01 },
      ],
    };
    const rows = toScatterData([ix]);
    expect(rows).toHaveLength(2);
    expect(rows[0]).toMatchObject({ ratio: 0.7071, q: 0.71, phase: "Pn3m" });
    expect(rows[1]).toMatchObject({ ratio: 1.0, q: 1.01, phase: "Pn3m" });
  });

  it("skips peaks whose ratio_position is out of bounds", () => {
    const ix: IndexEntry = {
      id: 1, exposure_id: 1, phase: "Pn3m", basis: 0.5, score: 1,
      r_squared: 0.99, lattice_d: 12, status: "candidate",
      predicted_q: [0.7],
      peaks: [{ peak_id: 10, ratio_position: 99, residual: 0, q_observed: 1.0 }],
    };
    expect(toScatterData([ix])).toEqual([]);
  });
});

describe("<MillerPlot>", () => {
  it("renders a container with testid 'miller-plot'", () => {
    const { container } = render(<MillerPlot indices={[]} />);
    expect(container.querySelector('[data-testid="miller-plot"]')).not.toBeNull();
  });

  it("invokes Plot.dot and Plot.linearRegressionY when indices have peaks", async () => {
    const Plot = await import("@observablehq/plot");
    (Plot.dot as unknown as { mockClear: () => void }).mockClear();
    (Plot.linearRegressionY as unknown as { mockClear: () => void }).mockClear();
    const ix: IndexEntry = {
      id: 1, exposure_id: 1, phase: "Pn3m", basis: 0.5, score: 1,
      r_squared: 0.99, lattice_d: 12, status: "candidate",
      predicted_q: [0.7, 1.0],
      peaks: [
        { peak_id: 10, ratio_position: 1, residual: 0, q_observed: 0.71 },
        { peak_id: 11, ratio_position: 2, residual: 0, q_observed: 1.0 },
      ],
    };
    render(<MillerPlot indices={[ix]} />);
    expect(Plot.dot).toHaveBeenCalled();
    expect(Plot.linearRegressionY).toHaveBeenCalled();
  });
});
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- MillerPlot`

- [ ] **Step 3: Implement `MillerPlot.tsx`**

Create `packages/HimalayaUI/frontend/src/components/MillerPlot.tsx`:

```tsx
import { useEffect, useRef } from "react";
import * as Plot from "@observablehq/plot";
import type { IndexEntry } from "../api";
import { phaseColor } from "../phases";

export interface MillerPlotProps {
  indices: IndexEntry[];
}

export interface ScatterRow { ratio: number; q: number; phase: string; color: string; }

export function toScatterData(indices: IndexEntry[]): ScatterRow[] {
  const rows: ScatterRow[] = [];
  for (const ix of indices) {
    if (ix.basis <= 0) continue;
    for (const pk of ix.peaks) {
      const pred = ix.predicted_q[pk.ratio_position - 1];
      if (pred == null) continue;
      const ratio = pred / ix.basis; // normalized √(sum of squares)
      rows.push({ ratio, q: pk.q_observed, phase: ix.phase, color: phaseColor(ix.phase) });
    }
  }
  return rows;
}

export function MillerPlot({ indices }: MillerPlotProps): JSX.Element {
  const hostRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const host = hostRef.current;
    if (!host) return;

    const data = toScatterData(indices);

    const el = Plot.plot({
      width:  host.clientWidth  || 360,
      height: host.clientHeight || 260,
      marginLeft: 50, marginBottom: 40,
      x: { label: "√(h²+k²+l²) (normalized)" },
      y: { label: "q (nm⁻¹)" },
      marks: data.length === 0 ? [] : [
        Plot.linearRegressionY(data, {
          x: "ratio", y: "q",
          stroke: "var(--color-fg-muted)",
          strokeDasharray: "4,3",
        }),
        Plot.dot(data, {
          x: "ratio", y: "q",
          fill: (d: ScatterRow) => d.color,
          stroke: "var(--color-bg)",
          strokeWidth: 1,
          r: 4,
          title: (d: ScatterRow) => `${d.phase}\nratio ${d.ratio.toFixed(3)}\nq ${d.q.toFixed(4)}`,
        }),
      ],
    });

    host.replaceChildren(el);
    return () => { host.replaceChildren(); };
  }, [indices]);

  return <div ref={hostRef} className="w-full h-full" data-testid="miller-plot" />;
}
```

- [ ] **Step 4: Run to verify pass**

Run: `npm test -- MillerPlot`
Expected: 4/4 pass.

- [ ] **Step 5: Full suite + build**

Run: `npm test && npm run build`

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/MillerPlot.tsx \
        packages/HimalayaUI/frontend/test/MillerPlot.test.tsx
git commit -m "feat(ui): MillerPlot scatter with linear fit per index"
```

---

## Task 11: Wire MillerPlot, PhasePanel, and new TraceViewer props into App.tsx

**Motivation:** Plug the pieces together. `rightTop` = MillerPlot, `rightBottom` = PhasePanel. Pass `activeGroupIndices` and resolved `hoveredIndex` into TraceViewer.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/App.tsx`
- Modify: `packages/HimalayaUI/frontend/test/smoke.test.tsx`

- [ ] **Step 1: Update `App.tsx`**

Replace `packages/HimalayaUI/frontend/src/App.tsx`:

```tsx
import { useEffect, useState, useMemo } from "react";
import "./styles.css";
import { useAppState } from "./state";
import {
  useExperiment, useSamples, useExposures, useTrace, usePeaks, useIndices, useGroups,
  useAddPeak, useRemovePeak,
} from "./queries";
import { Navbar } from "./components/Navbar";
import { Layout } from "./components/Layout";
import { SampleList } from "./components/SampleList";
import { UserModal } from "./components/UserModal";
import { ExposureList } from "./components/ExposureList";
import { TraceViewer } from "./components/TraceViewer";
import { StaleIndicesBanner } from "./components/StaleIndicesBanner";
import { MillerPlot } from "./components/MillerPlot";
import { PhasePanel } from "./components/PhasePanel";

const EXPERIMENT_ID = 1;

export function App(): JSX.Element {
  const username          = useAppState((s) => s.username);
  const activeSampleId    = useAppState((s) => s.activeSampleId);
  const activeExposureId  = useAppState((s) => s.activeExposureId);
  const hoveredIndexId    = useAppState((s) => s.hoveredIndexId);
  const setUsername       = useAppState((s) => s.setUsername);
  const setActiveSample   = useAppState((s) => s.setActiveSample);
  const setActiveExposure = useAppState((s) => s.setActiveExposure);

  const [modalOpen, setModalOpen] = useState<boolean>(!username);

  const experimentQ = useExperiment(EXPERIMENT_ID);
  const samplesQ    = useSamples(EXPERIMENT_ID);
  const exposuresQ  = useExposures(activeSampleId);
  const traceQ      = useTrace(activeExposureId);
  const peaksQ      = usePeaks(activeExposureId);
  const indicesQ    = useIndices(activeExposureId);
  const groupsQ     = useGroups(activeExposureId);

  const addPeak    = useAddPeak(activeExposureId ?? 0);
  const removePeak = useRemovePeak(activeExposureId ?? 0);

  useEffect(() => {
    const exposures = exposuresQ.data ?? [];
    if (exposures.length === 0) return;
    const stillValid = exposures.some((e) => e.id === activeExposureId);
    if (!stillValid) setActiveExposure(exposures[0]!.id);
  }, [exposuresQ.data, activeExposureId, setActiveExposure]);

  const indices = indicesQ.data ?? [];
  const activeGroup = (groupsQ.data ?? []).find((g) => g.active);
  const activeGroupIndices = useMemo(
    () => (activeGroup?.members ?? [])
      .map((id) => indices.find((i) => i.id === id))
      .filter((i): i is NonNullable<typeof i> => i != null),
    [activeGroup, indices],
  );
  const hoveredIndex = hoveredIndexId != null
    ? indices.find((i) => i.id === hoveredIndexId)
    : undefined;

  const samples      = samplesQ.data ?? [];
  const activeSample = samples.find((s) => s.id === activeSampleId);
  const bootError    = experimentQ.error ?? samplesQ.error;

  const breadcrumb = bootError
    ? `Error: ${(bootError as Error).message}`
    : (experimentQ.data?.name ?? "experiment")
        + (activeSample
          ? ` › ${activeSample.label ?? ""} ${activeSample.name ?? ""}`.trimEnd()
          : "");

  return (
    <>
      <Navbar
        breadcrumb={breadcrumb}
        username={username}
        onUserClick={() => setModalOpen(true)}
      />
      <Layout
        left={
          <SampleList
            samples={samples}
            activeId={activeSampleId}
            onSelect={setActiveSample}
          />
        }
        centerTop={
          <div className="flex flex-col flex-1 min-h-0">
            <StaleIndicesBanner exposureId={activeExposureId} />
            {traceQ.data && peaksQ.data && activeExposureId !== undefined ? (
              <TraceViewer
                trace={traceQ.data}
                peaks={peaksQ.data}
                activeGroupIndices={activeGroupIndices}
                hoveredIndex={hoveredIndex}
                onAddPeak={(q) => addPeak.mutate(q)}
                onRemovePeak={(peakId) => removePeak.mutate(peakId)}
              />
            ) : (
              <p className="text-fg-muted italic flex-1 flex items-center justify-center">
                Select an exposure to view its trace.
              </p>
            )}
          </div>
        }
        centerBottom={<ExposureList />}
        rightTop={<MillerPlot indices={activeGroupIndices} />}
        rightBottom={<PhasePanel exposureId={activeExposureId} />}
      />
      <UserModal
        open={modalOpen}
        onSelect={(name) => {
          setUsername(name);
          setModalOpen(false);
        }}
        onClose={() => setModalOpen(false)}
      />
    </>
  );
}
```

- [ ] **Step 2: Update smoke test mocks**

Replace the `mockFetch({...})` call in `packages/HimalayaUI/frontend/test/smoke.test.tsx` with:

```tsx
mockFetch({
  "/api/users": [],
  "/api/experiments/1": {
    id: 1, name: "demo", path: "/x", data_dir: "/x/data",
    analysis_dir: "/x/analysis", manifest_path: null,
    created_at: "2026-04-22T00:00:00Z",
  },
  "/api/experiments/1/samples": [
    { id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null, tags: [] },
  ],
  "/api/samples/10/exposures": [],
});
```

No change from Plan 4's mock set — the new queries (`useIndices`, `useGroups`) are disabled until an exposure is active, which this smoke test deliberately avoids.

- [ ] **Step 3: Run all**

Run: `npm test && npm run build`
Expected: all green.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/App.tsx \
        packages/HimalayaUI/frontend/test/smoke.test.tsx
git commit -m "feat(ui): wire MillerPlot + PhasePanel + overlay into App"
```

---

## Task 12: E2E — hover preview on an alternative index

**Motivation:** One E2E that exercises the whole assembled flow: boot with an active exposure, hover an alternative in the PhasePanel, confirm the POST wiring for `+`.

**Files:**
- Modify: `packages/HimalayaUI/frontend/e2e/smoke.spec.ts`

- [ ] **Step 1: Add the test**

Append to `packages/HimalayaUI/frontend/e2e/smoke.spec.ts` at top level (after the last existing test):

```ts
test("hovering an alternative and clicking + posts to groups", async ({ page }) => {
  let posted: { index_id: number } | null = null;

  await page.route("**/api/users", (r) => r.fulfill({ status: 200, body: "[]" }));
  await page.route("**/api/experiments/1", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify({
      id: 1, name: "demo", path: "/x", data_dir: "/x/data",
      analysis_dir: "/x/analysis", manifest_path: null,
      created_at: "2026-04-22T00:00:00Z",
    }),
  }));
  await page.route("**/api/experiments/1/samples", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify([
      { id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null, tags: [] },
    ]),
  }));
  await page.route("**/api/samples/10/exposures", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify([
      { id: 100, sample_id: 10, filename: "demo", kind: "file",
        selected: true, tags: [], sources: [] },
    ]),
  }));
  await page.route("**/api/exposures/100/trace", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify({
      q: [0.05, 0.1, 0.15, 0.2], I: [100, 80, 90, 70], sigma: [5, 4, 4, 3],
    }),
  }));
  await page.route("**/api/exposures/100/peaks", (r) =>
    r.fulfill({ status: 200, body: "[]" }));
  await page.route("**/api/exposures/100/indices", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify([
      { id: 1, exposure_id: 100, phase: "Pn3m", basis: 0.5, score: 1.0,
        r_squared: 0.99, lattice_d: 12.5, status: "candidate",
        predicted_q: [0.7, 0.86], peaks: [] },
      { id: 2, exposure_id: 100, phase: "Im3m", basis: 0.3, score: 0.6,
        r_squared: 0.71, lattice_d: 9.1, status: "candidate",
        predicted_q: [0.42, 0.6], peaks: [] },
    ]),
  }));
  await page.route("**/api/exposures/100/groups", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify([
      { id: 1, exposure_id: 100, kind: "auto", active: true, members: [1] },
    ]),
  }));
  await page.route("**/api/groups/1/members", async (r) => {
    if (r.request().method() === "POST") {
      posted = r.request().postDataJSON() as { index_id: number };
      return r.fulfill({
        status: 200,
        body: JSON.stringify({
          id: 2, exposure_id: 100, kind: "custom", active: true, members: [1, posted.index_id],
        }),
      });
    }
    return r.fulfill({ status: 404, body: "nope" });
  });

  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: { username: "alice", activeSampleId: 10, activeExposureId: 100 },
      version: 1,
    }));
  });

  await page.goto("/");
  const altRow = page.locator('[data-alternative-id="2"]');
  await altRow.waitFor();
  await altRow.hover();
  const addBtn = altRow.getByRole("button", { name: /add index 2/i });
  await addBtn.click();
  await expect.poll(() => posted?.index_id, { timeout: 2000 }).toBe(2);
});
```

- [ ] **Step 2: Build + run E2E**

Run: `npm run build && npm run e2e`
Expected: all tests pass (existing 5 + new 1).

If the hover doesn't register, check that `useAppState` selectors are preventing the `+ button` re-render — the button's `onClick` just calls `addMember.mutate(ix.id)`, which doesn't depend on hover. The hover test is really covered by the PhasePanel unit test.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/e2e/smoke.spec.ts
git commit -m "test(e2e): hover alternative and add to group"
```

---

## Final verification

From repo root:

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: 170+ tests pass.

From `packages/HimalayaUI/frontend/`:

```bash
npm test && npm run build && npm run e2e
```

Expected: all unit tests pass (43 existing + ~14 new), build clean, E2E all pass (5 existing + 1 new).

Hand off to Plan 6 (Properties panel tabs: Peaks table / Tags / Notes editor) and Plan 6+ (Recent section, export button, per-user audit view).
