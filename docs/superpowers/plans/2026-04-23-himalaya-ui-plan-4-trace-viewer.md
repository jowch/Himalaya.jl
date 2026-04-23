# HimalayaUI Plan 4 — Trace Viewer & Peak Editing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Render the SAXS trace for the active exposure, let the user add/remove peaks by clicking, and prompt for re-analysis when indices go stale.

**Architecture:** A new backend endpoint parses the raw `.dat` file on demand and returns `{q, I, sigma}` arrays. The frontend gains Observable-Plot-based rendering inside a React component that uses TanStack Query for data and mutations. Exposure selection lives in Zustand (already wired). This plan covers trace rendering, peak markers, click-to-add/remove, and a stale-index banner — but not the phase panel or Miller plot (Plan 5).

**Tech Stack:** Julia 1.12, Oxygen.jl, HimalayaUI (existing backend patterns), React 18, TypeScript strict (`exactOptionalPropertyTypes: true`), Vite, Zustand, **TanStack Query v5**, **Observable Plot** (new dep), Vitest + RTL, Playwright.

**Assumptions about the starting state** (verified at plan-write time):
- `main` is at commit `e303e3a`. All 164 backend tests pass. Frontend: 15 unit + 4 E2E tests pass.
- Existing backend endpoints (all landed in Plan 2): `GET /api/samples/{id}/exposures`, `GET /api/exposures/{id}/peaks`, `POST /api/exposures/{id}/peaks`, `DELETE /api/peaks/{id}`, `GET /api/exposures/{id}/indices`, `POST /api/exposures/{id}/analyze`. `POST /api/exposures/{id}/peaks` marks all indices for the exposure `stale` and returns `stale_indices` count.
- Route files live at `packages/HimalayaUI/src/routes_*.jl`; each exports a `register_*_routes!()` that's called from `server.jl`. Tests use the harness at `packages/HimalayaUI/test/test_http.jl` (`with_test_server(db) do port, base ... end`).
- `load_dat(path)` in `packages/HimalayaUI/src/datfile.jl` returns `(q, I, σ)` from the three-column file.
- Test fixture `test/data/example_tot.dat` (~900 rows of a real SAXS trace) is already used by `test_routes_peaks.jl`; Plan 4 reuses it.
- Frontend file layout follows the Plan 3.5 shape: `src/{api.ts,queries.ts,state.ts,App.tsx,main.tsx,ErrorBoundary.tsx,components/*}`. `queries.ts` currently exports `queryKeys`, `useExperiment`, `useSamples`.
- Zustand store (`state.ts`) already has `activeSampleId`, `activeExposureId`, `setActiveSample`, `setActiveExposure`.
- `api.ts` already exposes `AuthOpts`, `ApiError`, `request<T>`, and the existing CRUD functions. Mutations accept an optional `opts?: AuthOpts`.

**Backend commands run from:** repo root (`/Users/me/projects/Himalaya.jl`).
**Frontend commands run from:** `packages/HimalayaUI/frontend/`.

---

## File Map

**Backend (new or modified):**
- Create: `packages/HimalayaUI/src/routes_trace.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl` (add `include("routes_trace.jl")`)
- Modify: `packages/HimalayaUI/src/server.jl` (call `register_trace_routes!()`)
- Create: `packages/HimalayaUI/test/test_routes_trace.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl` (include the new test file)

**Frontend (new or modified):**
- Modify: `packages/HimalayaUI/frontend/package.json` (add `@observablehq/plot`)
- Modify: `packages/HimalayaUI/frontend/src/api.ts` (Exposure/Peak/Index/Trace types + fetchers + mutations)
- Modify: `packages/HimalayaUI/frontend/src/queries.ts` (useExposures, useTrace, usePeaks, useIndices + useAddPeak, useRemovePeak, useReanalyzeExposure hooks)
- Create: `packages/HimalayaUI/frontend/src/components/ExposureList.tsx`
- Create: `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx`
- Create: `packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx`
- Modify: `packages/HimalayaUI/frontend/src/App.tsx` (wire new components into Layout + auto-select active exposure)
- Create tests: `test/ExposureList.test.tsx`, `test/TraceViewer.test.tsx`, `test/StaleIndicesBanner.test.tsx`, `test/queries.test.tsx`
- Modify: `packages/HimalayaUI/frontend/test/smoke.test.tsx` (extend mocks to cover new endpoints)
- Modify: `packages/HimalayaUI/frontend/e2e/smoke.spec.ts` (new E2E for add-peak happy path)

---

## Task 1: Backend trace endpoint

**Motivation:** The raw `(q, I, σ)` arrays aren't in the DB — the analysis pipeline reads them transiently from `.dat` files. Surface them as an API endpoint so the frontend can render the trace.

**Files:**
- Create: `packages/HimalayaUI/src/routes_trace.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl`
- Modify: `packages/HimalayaUI/src/server.jl`
- Create: `packages/HimalayaUI/test/test_routes_trace.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/test_routes_trace.jl`:

```julia
using Test, HTTP, JSON3

@testset "trace route" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")

    with_test_server(db) do port, base
        r = HTTP.get("$base/api/exposures/$e_id/trace")
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test haskey(body, :q) && haskey(body, :I) && haskey(body, :sigma)
        @test length(body.q) == length(body.I) == length(body.sigma)
        @test length(body.q) > 100
        @test all(q -> q > 0, body.q)

        # 404 for unknown exposure
        r = HTTP.get("$base/api/exposures/99999/trace"; status_exception = false)
        @test r.status == 404
    end
end
```

Also add the include line to `packages/HimalayaUI/test/runtests.jl`. Look at the file — tests are listed one-per-line with `include("test_*.jl")`. Add `include("test_routes_trace.jl")` next to the other `test_routes_*` lines (after `test_routes_peaks.jl`).

- [ ] **Step 2: Run test to verify failure**

Run from repo root: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'`

Expected: the new testset fails with 404 (no `/api/exposures/{id}/trace` route registered). All other tests still pass.

- [ ] **Step 3: Create `routes_trace.jl`**

Create `packages/HimalayaUI/src/routes_trace.jl`:

```julia
using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_trace_routes!()
    @get "/api/exposures/{id}/trace" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT e.filename, e.kind, x.analysis_dir
             FROM exposures e JOIN samples s ON s.id = e.sample_id
             JOIN experiments x ON x.id = s.experiment_id
             WHERE e.id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure not found")))
        row = rows[1]
        String(row.kind) == "file" || return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "trace not available for derived exposures")))
        row.filename === missing && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure has no filename")))

        path = joinpath(String(row.analysis_dir), String(row.filename) * ".dat")
        isfile(path) || return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => ".dat file not found: $path")))

        q, I, σ = load_dat(path)
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:q => q, :I => I, :sigma => σ)))
    end
end
```

- [ ] **Step 4: Wire the route file**

In `packages/HimalayaUI/src/HimalayaUI.jl`, add after the existing `include("routes_peaks.jl")`:

```julia
include("routes_trace.jl")
```

In `packages/HimalayaUI/src/server.jl`, find the block that calls `register_*_routes!()` (search for `register_peaks_routes!()`) and add after it:

```julia
    register_trace_routes!()
```

Keep the same indentation as the neighbouring lines.

- [ ] **Step 5: Run test to verify pass**

Run from repo root: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'`

Expected: all tests pass, including the new trace route testset (2 assertions per branch — the 200 and 404 paths).

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/routes_trace.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/src/server.jl \
        packages/HimalayaUI/test/test_routes_trace.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "feat(api): GET /api/exposures/:id/trace returns q/I/sigma arrays"
```

---

## Task 2: Frontend API types and fetchers

**Motivation:** Before we can wire queries and components, add the typed fetchers for exposures, trace, peaks, indices, and the peak/analyze mutations. All changes are in one file with no React coupling — easy to test in isolation.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/api.ts`
- Modify: `packages/HimalayaUI/frontend/test/api.test.ts`

- [ ] **Step 1: Write failing tests**

Append these to `packages/HimalayaUI/frontend/test/api.test.ts` inside the existing `describe("api", () => { ... })` block:

```ts
it("getTrace returns parsed q/I/sigma arrays", async () => {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify({ q: [0.1, 0.2], I: [10, 20], sigma: [1, 2] }),
      { status: 200 }),
  );
  const t = await api.getTrace(42);
  expect(t.q).toEqual([0.1, 0.2]);
  expect(t.I).toEqual([10, 20]);
  expect(t.sigma).toEqual([1, 2]);
});

it("addPeak posts {q} with X-Username and returns parsed peak", async () => {
  const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify({
      id: 7, exposure_id: 42, q: 0.15, source: "manual", stale_indices: 3,
    }), { status: 201 }),
  );
  const p = await api.addPeak(42, 0.15, { username: "alice" });
  expect(p.id).toBe(7);
  expect(p.source).toBe("manual");
  expect(p.stale_indices).toBe(3);
  const init = fetchSpy.mock.calls[0]![1] as RequestInit;
  expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
  expect(init.body).toBe(JSON.stringify({ q: 0.15 }));
});

it("removePeak sends DELETE with X-Username", async () => {
  const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(null, { status: 204 }),
  );
  await api.removePeak(7, { username: "alice" });
  const [url, init] = fetchSpy.mock.calls[0]! as [string, RequestInit];
  expect(url).toBe("/api/peaks/7");
  expect(init.method).toBe("DELETE");
  expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
});

it("reanalyzeExposure posts empty body and returns {id, analyzed}", async () => {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify({ id: 42, analyzed: true }), { status: 200 }),
  );
  const r = await api.reanalyzeExposure(42, { username: "alice" });
  expect(r.id).toBe(42);
  expect(r.analyzed).toBe(true);
});
```

- [ ] **Step 2: Run tests to verify failure**

Run from `packages/HimalayaUI/frontend/`: `npm test -- api.test`

Expected: four new tests fail with `api.getTrace is not a function` etc.

- [ ] **Step 3: Extend `api.ts`**

Append to `packages/HimalayaUI/frontend/src/api.ts` (after the existing Sample exports):

```ts
// Exposures
export interface Exposure {
  id: number;
  sample_id: number;
  filename: string | null;
  kind: "file" | "averaged" | "background_subtracted";
  selected: boolean;
  tags: unknown[];
  sources: unknown[];
}

export const listExposures = (sample_id: number) =>
  request<Exposure[]>("GET", `/api/samples/${sample_id}/exposures`);

// Trace
export interface Trace {
  q: number[];
  I: number[];
  sigma: number[];
}

export const getTrace = (exposure_id: number) =>
  request<Trace>("GET", `/api/exposures/${exposure_id}/trace`);

// Peaks
export interface Peak {
  id: number;
  exposure_id: number;
  q: number;
  intensity: number | null;
  prominence: number | null;
  sharpness: number | null;
  source: "auto" | "manual";
}

export interface PeakCreated extends Peak { stale_indices: number }

export const listPeaks = (exposure_id: number) =>
  request<Peak[]>("GET", `/api/exposures/${exposure_id}/peaks`);
export const addPeak = (exposure_id: number, q: number, opts?: AuthOpts) =>
  request<PeakCreated>("POST", `/api/exposures/${exposure_id}/peaks`, { q }, opts);
export const removePeak = (peak_id: number, opts?: AuthOpts) =>
  request<void>("DELETE", `/api/peaks/${peak_id}`, undefined, opts);

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

// Analysis
export const reanalyzeExposure = (exposure_id: number, opts?: AuthOpts) =>
  request<{ id: number; analyzed: boolean }>("POST", `/api/exposures/${exposure_id}/analyze`, {}, opts);
```

- [ ] **Step 4: Run tests to verify pass**

Run: `npm test -- api.test`
Expected: 4 old + 4 new tests pass (8 total).

- [ ] **Step 5: Run full suite and build**

Run: `npm test && npm run build`
Expected: all green.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/api.ts \
        packages/HimalayaUI/frontend/test/api.test.ts
git commit -m "feat(ui): add api fetchers for exposures, trace, peaks, indices, analyze"
```

---

## Task 3: TanStack Query hooks for exposures/trace/peaks/indices + mutations

**Motivation:** Expose the new fetchers through React Query hooks with proper cache keys and invalidation. This is the glue layer between `api.ts` and components.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/queries.ts`
- Create: `packages/HimalayaUI/frontend/test/queries.test.tsx`

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/queries.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor, act } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useExposures, useTrace, usePeaks, useIndices,
  useAddPeak, useRemovePeak, useReanalyzeExposure,
  queryKeys,
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
    new Response(JSON.stringify(body), {
      status, headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("queries", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("useExposures fetches when sampleId is provided", async () => {
    mockOnce(200, [{ id: 1, sample_id: 10, filename: "f", kind: "file",
                     selected: false, tags: [], sources: [] }]);
    const { wrapper } = withClient();
    const { result } = renderHook(() => useExposures(10), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(result.current.data).toHaveLength(1);
  });

  it("useExposures is disabled when sampleId is undefined", () => {
    const fetchSpy = vi.spyOn(global, "fetch");
    const { wrapper } = withClient();
    const { result } = renderHook(() => useExposures(undefined), { wrapper });
    expect(result.current.fetchStatus).toBe("idle");
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it("useTrace fetches for a given exposureId", async () => {
    mockOnce(200, { q: [0.1], I: [10], sigma: [1] });
    const { wrapper } = withClient();
    const { result } = renderHook(() => useTrace(42), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(result.current.data?.q).toEqual([0.1]);
  });

  it("usePeaks fetches for a given exposureId", async () => {
    mockOnce(200, []);
    const { wrapper } = withClient();
    const { result } = renderHook(() => usePeaks(42), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
  });

  it("useIndices fetches for a given exposureId", async () => {
    mockOnce(200, []);
    const { wrapper } = withClient();
    const { result } = renderHook(() => useIndices(42), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
  });

  it("useAddPeak invalidates peaks and indices for the exposure on success", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData(queryKeys.peaks(42), []);
    client.setQueryData(queryKeys.indices(42), []);
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(201, { id: 1, exposure_id: 42, q: 0.1, intensity: null,
                    prominence: null, sharpness: null, source: "manual",
                    stale_indices: 0 });
    const { result } = renderHook(() => useAddPeak(42), { wrapper });
    await act(async () => { await result.current.mutateAsync(0.1); });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.peaks(42) });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.indices(42) });
  });

  it("useRemovePeak invalidates peaks and indices for the exposure on success", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(204, null);
    const { result } = renderHook(() => useRemovePeak(42), { wrapper });
    await act(async () => { await result.current.mutateAsync(7); });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.peaks(42) });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.indices(42) });
  });

  it("useReanalyzeExposure invalidates peaks and indices on success", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(200, { id: 42, analyzed: true });
    const { result } = renderHook(() => useReanalyzeExposure(42), { wrapper });
    await act(async () => { await result.current.mutateAsync(); });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.peaks(42) });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.indices(42) });
  });
});
```

- [ ] **Step 2: Run tests to verify failure**

Run: `npm test -- queries.test`
Expected: all 8 tests fail with "useExposures is not a function" / "queryKeys.peaks is not a function" etc.

- [ ] **Step 3: Extend `queries.ts`**

Replace `packages/HimalayaUI/frontend/src/queries.ts`:

```ts
import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";
import * as api from "./api";
import { useAppState } from "./state";

export const queryKeys = {
  experiment: (id: number) => ["experiment", id] as const,
  samples:    (experimentId: number) => ["experiment", experimentId, "samples"] as const,
  exposures:  (sampleId: number) => ["sample", sampleId, "exposures"] as const,
  trace:      (exposureId: number) => ["exposure", exposureId, "trace"] as const,
  peaks:      (exposureId: number) => ["exposure", exposureId, "peaks"] as const,
  indices:    (exposureId: number) => ["exposure", exposureId, "indices"] as const,
};

export function useExperiment(id: number) {
  return useQuery({
    queryKey: queryKeys.experiment(id),
    queryFn: () => api.getExperiment(id),
  });
}

export function useSamples(experimentId: number) {
  return useQuery({
    queryKey: queryKeys.samples(experimentId),
    queryFn: () => api.listSamples(experimentId),
  });
}

export function useExposures(sampleId: number | undefined) {
  return useQuery({
    queryKey: ["sample", sampleId ?? "none", "exposures"] as const,
    queryFn: () => api.listExposures(sampleId as number),
    enabled: sampleId !== undefined,
  });
}

export function useTrace(exposureId: number | undefined) {
  return useQuery({
    queryKey: ["exposure", exposureId ?? "none", "trace"] as const,
    queryFn: () => api.getTrace(exposureId as number),
    enabled: exposureId !== undefined,
  });
}

export function usePeaks(exposureId: number | undefined) {
  return useQuery({
    queryKey: ["exposure", exposureId ?? "none", "peaks"] as const,
    queryFn: () => api.listPeaks(exposureId as number),
    enabled: exposureId !== undefined,
  });
}

export function useIndices(exposureId: number | undefined) {
  return useQuery({
    queryKey: ["exposure", exposureId ?? "none", "indices"] as const,
    queryFn: () => api.listIndices(exposureId as number),
    enabled: exposureId !== undefined,
  });
}

function invalidateExposure(qc: ReturnType<typeof useQueryClient>, exposureId: number): void {
  qc.invalidateQueries({ queryKey: queryKeys.peaks(exposureId) });
  qc.invalidateQueries({ queryKey: queryKeys.indices(exposureId) });
}

export function useAddPeak(exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (q: number) => api.addPeak(exposureId, q, { username }),
    onSuccess: () => invalidateExposure(qc, exposureId),
  });
}

export function useRemovePeak(exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (peakId: number) => api.removePeak(peakId, { username }),
    onSuccess: () => invalidateExposure(qc, exposureId),
  });
}

export function useReanalyzeExposure(exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: () => api.reanalyzeExposure(exposureId, { username }),
    onSuccess: () => invalidateExposure(qc, exposureId),
  });
}
```

- [ ] **Step 4: Run tests to verify pass**

Run: `npm test -- queries.test`
Expected: all 8 pass.

- [ ] **Step 5: Full suite + build**

Run: `npm test && npm run build`
Expected: clean.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/queries.ts \
        packages/HimalayaUI/frontend/test/queries.test.tsx
git commit -m "feat(ui): TanStack Query hooks for exposures, trace, peaks, indices, mutations"
```

---

## Task 4: ExposureList component

**Motivation:** Give the user a way to select which exposure to view. This is the first piece of the bottom-center "Properties panel". In Plan 4 it's the whole panel; Plans 5–6 will add tabs (Peaks/Tags/Notes).

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ExposureList.tsx`
- Create: `packages/HimalayaUI/frontend/test/ExposureList.test.tsx`

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/ExposureList.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, fireEvent, waitFor } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { ExposureList } from "../src/components/ExposureList";
import { useAppState } from "../src/state";

beforeEach(() => {
  vi.restoreAllMocks();
  useAppState.setState({ activeSampleId: undefined, activeExposureId: undefined });
});

function mockExposures(items: unknown[]): void {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify(items), {
      status: 200, headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("<ExposureList>", () => {
  it("shows a hint when no sample is active", () => {
    renderWithProviders(<ExposureList />);
    expect(screen.getByText(/no sample selected/i)).toBeInTheDocument();
  });

  it("renders one row per exposure and marks the active one", async () => {
    useAppState.setState({ activeSampleId: 1, activeExposureId: 11 });
    mockExposures([
      { id: 10, sample_id: 1, filename: "a", kind: "file",
        selected: false, tags: [], sources: [] },
      { id: 11, sample_id: 1, filename: "b", kind: "file",
        selected: true,  tags: [], sources: [] },
    ]);
    renderWithProviders(<ExposureList />);
    await waitFor(() => expect(screen.getByText("a")).toBeInTheDocument());
    expect(screen.getByText("b")).toBeInTheDocument();
    const active = document.querySelector('[data-exposure-id="11"][data-active="true"]');
    expect(active).not.toBeNull();
  });

  it("clicking a row updates activeExposureId in Zustand", async () => {
    useAppState.setState({ activeSampleId: 1, activeExposureId: undefined });
    mockExposures([
      { id: 10, sample_id: 1, filename: "a", kind: "file",
        selected: false, tags: [], sources: [] },
    ]);
    renderWithProviders(<ExposureList />);
    await waitFor(() => expect(screen.getByText("a")).toBeInTheDocument());
    const row = document.querySelector<HTMLElement>('[data-exposure-id="10"]')!;
    fireEvent.click(row);
    expect(useAppState.getState().activeExposureId).toBe(10);
  });
});
```

- [ ] **Step 2: Run tests to verify failure**

Run: `npm test -- ExposureList`
Expected: 3 tests fail (`ExposureList` not defined).

- [ ] **Step 3: Implement `ExposureList.tsx`**

Create `packages/HimalayaUI/frontend/src/components/ExposureList.tsx`:

```tsx
import { useAppState } from "../state";
import { useExposures } from "../queries";

export function ExposureList(): JSX.Element {
  const activeSampleId   = useAppState((s) => s.activeSampleId);
  const activeExposureId = useAppState((s) => s.activeExposureId);
  const setActiveExposure = useAppState((s) => s.setActiveExposure);

  if (activeSampleId === undefined) {
    return <p className="text-fg-muted italic">No sample selected.</p>;
  }

  const q = useExposures(activeSampleId);

  if (q.isPending) return <p className="text-fg-muted">Loading exposures…</p>;
  if (q.error)     return <p className="text-error">Error: {(q.error as Error).message}</p>;

  const exposures = q.data ?? [];
  if (exposures.length === 0) {
    return <p className="text-fg-muted italic">No exposures for this sample.</p>;
  }

  return (
    <ul className="list-none flex flex-col gap-1">
      {exposures.map((e) => {
        const active = e.id === activeExposureId;
        return (
          <li
            key={e.id}
            data-exposure-id={e.id}
            data-active={active}
            className={
              "flex items-center gap-2 px-2 py-1 rounded-md cursor-pointer border-l-2 " +
              (active ? "bg-bg-elevated border-accent" : "border-transparent hover:bg-bg-hover")
            }
            onClick={() => setActiveExposure(e.id)}
          >
            <span className="font-mono text-[13px]">{e.filename ?? "(derived)"}</span>
            <span className="text-fg-muted text-[12px]">{e.kind}</span>
          </li>
        );
      })}
    </ul>
  );
}
```

- [ ] **Step 4: Run tests to verify pass**

Run: `npm test -- ExposureList`
Expected: 3/3 pass.

- [ ] **Step 5: Full suite + build**

Run: `npm test && npm run build`
Expected: clean.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ExposureList.tsx \
        packages/HimalayaUI/frontend/test/ExposureList.test.tsx
git commit -m "feat(ui): ExposureList in properties panel"
```

---

## Task 5: Install Observable Plot; TraceViewer scaffold

**Motivation:** Get the chart library in place and render the basic trace (log-log I(q) with σ ribbon) for the active exposure. No peaks yet — that's Task 6.

**Files:**
- Modify: `packages/HimalayaUI/frontend/package.json`
- Create: `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx`
- Create: `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx`

- [ ] **Step 1: Install `@observablehq/plot`**

Run from `packages/HimalayaUI/frontend/`:

```bash
npm install @observablehq/plot@^0.6.17
```

Expected: `package.json` dependencies gain `@observablehq/plot`.

- [ ] **Step 2: Write failing tests**

Create `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx`:

```tsx
import { describe, it, expect, vi } from "vitest";
import { render } from "@testing-library/react";
import { TraceViewer } from "../src/components/TraceViewer";

vi.mock("@observablehq/plot", () => ({
  plot: vi.fn(() => {
    const el = document.createElement("div");
    el.setAttribute("data-testid", "plot-svg");
    return el;
  }),
  areaY: vi.fn(() => ({ _kind: "areaY" })),
  line:  vi.fn(() => ({ _kind: "line" })),
  dot:   vi.fn(() => ({ _kind: "dot" })),
}));

describe("<TraceViewer>", () => {
  it("renders a container with testid 'trace-viewer'", () => {
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 30], sigma: [1, 1, 1] };
    const { container } = render(
      <TraceViewer
        trace={trace}
        peaks={[]}
        onAddPeak={() => {}}
        onRemovePeak={() => {}}
      />,
    );
    expect(container.querySelector('[data-testid="trace-viewer"]')).not.toBeNull();
  });

  it("invokes Plot.plot with marks including line and areaY for trace", async () => {
    const Plot = await import("@observablehq/plot");
    const trace = { q: [0.1, 0.2], I: [10, 20], sigma: [1, 1] };
    render(
      <TraceViewer trace={trace} peaks={[]} onAddPeak={() => {}} onRemovePeak={() => {}} />,
    );
    expect(Plot.plot).toHaveBeenCalled();
    expect(Plot.areaY).toHaveBeenCalled();
    expect(Plot.line).toHaveBeenCalled();
  });
});
```

- [ ] **Step 3: Run tests to verify failure**

Run: `npm test -- TraceViewer`
Expected: 2 tests fail (`TraceViewer` not defined).

- [ ] **Step 4: Implement `TraceViewer.tsx`**

Create `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx`:

```tsx
import { useEffect, useRef } from "react";
import * as Plot from "@observablehq/plot";
import type { Trace, Peak } from "../api";

export interface TraceViewerProps {
  trace: Trace;
  peaks: Peak[];
  onAddPeak: (q: number) => void;
  onRemovePeak: (peakId: number) => void;
}

export function TraceViewer({
  trace, peaks, onAddPeak: _onAddPeak, onRemovePeak: _onRemovePeak,
}: TraceViewerProps): JSX.Element {
  // onAddPeak/onRemovePeak are wired in Tasks 7–8; parameters are aliased to
  // underscore-prefixed locals here to silence TS6133 "declared but never used".
  void _onAddPeak; void _onRemovePeak;
  const hostRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const host = hostRef.current;
    if (!host) return;

    const data = trace.q.map((q, i) => ({
      q, I: trace.I[i]!,
      lo: Math.max(1e-12, trace.I[i]! - trace.sigma[i]!),
      hi: trace.I[i]! + trace.sigma[i]!,
    }));

    const el = Plot.plot({
      width: host.clientWidth || 400,
      height: host.clientHeight || 300,
      marginLeft: 50, marginBottom: 40,
      x: { type: "log", label: "q (nm⁻¹)" },
      y: { type: "log", label: "I (a.u.)" },
      marks: [
        Plot.areaY(data, { x: "q", y1: "lo", y2: "hi",
          fill: "var(--color-accent)", fillOpacity: 0.15 }),
        Plot.line(data,  { x: "q", y: "I",
          stroke: "var(--color-fg)", strokeWidth: 1 }),
      ],
    });

    host.replaceChildren(el);
    return () => { host.replaceChildren(); };
  }, [trace, peaks]);

  return <div ref={hostRef} className="w-full h-full" data-testid="trace-viewer" />;
}
```

- [ ] **Step 5: Run tests to verify pass**

Run: `npm test -- TraceViewer`
Expected: 2/2 pass.

- [ ] **Step 6: Full suite + build**

Run: `npm test && npm run build`
Expected: clean.

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/frontend/package.json \
        packages/HimalayaUI/frontend/package-lock.json \
        packages/HimalayaUI/frontend/src/components/TraceViewer.tsx \
        packages/HimalayaUI/frontend/test/TraceViewer.test.tsx
git commit -m "feat(ui): TraceViewer scaffold renders log-log I(q) with sigma ribbon"
```

---

## Task 6: Render peak markers on the trace

**Motivation:** Show where the auto-detected and manual peaks land. Auto peaks are filled; manual peaks are outlined in a warning color.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx`
- Modify: `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx`

- [ ] **Step 1: Write a failing test for peak rendering**

Add to `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx`:

```tsx
it("passes auto peaks and manual peaks to separate dot marks", async () => {
  const Plot = await import("@observablehq/plot");
  (Plot.dot as unknown as { mockClear: () => void }).mockClear();
  const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 30], sigma: [1, 1, 1] };
  const peaks = [
    { id: 1, exposure_id: 1, q: 0.1, intensity: null, prominence: null,
      sharpness: null, source: "auto"   as const },
    { id: 2, exposure_id: 1, q: 0.2, intensity: null, prominence: null,
      sharpness: null, source: "manual" as const },
  ];
  render(
    <TraceViewer trace={trace} peaks={peaks}
      onAddPeak={() => {}} onRemovePeak={() => {}} />,
  );
  // One call with auto peaks, one call with manual peaks
  expect((Plot.dot as unknown as { mock: { calls: unknown[] } }).mock.calls.length).toBe(2);
});
```

- [ ] **Step 2: Run test to verify failure**

Run: `npm test -- TraceViewer`
Expected: the new test fails — current `TraceViewer` calls `Plot.dot` zero times.

- [ ] **Step 3: Extend `TraceViewer.tsx`**

Replace the `marks` array inside `Plot.plot({...})` with:

```tsx
marks: [
  Plot.areaY(data, { x: "q", y1: "lo", y2: "hi",
    fill: "var(--color-accent)", fillOpacity: 0.15 }),
  Plot.line(data, { x: "q", y: "I",
    stroke: "var(--color-fg)", strokeWidth: 1 }),
  Plot.dot(withIntensity(peaks.filter((p) => p.source === "auto"), trace),
    { x: "q", y: "I",
      fill: "var(--color-accent)", r: 4 }),
  Plot.dot(withIntensity(peaks.filter((p) => p.source === "manual"), trace),
    { x: "q", y: "I",
      stroke: "var(--color-warning)", strokeWidth: 2, fill: "none", r: 5 }),
],
```

Above the component (module scope), add:

```tsx
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
```

- [ ] **Step 4: Run test to verify pass**

Run: `npm test -- TraceViewer`
Expected: all 3 tests pass.

- [ ] **Step 5: Full suite + build**

Run: `npm test && npm run build`
Expected: clean.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/TraceViewer.tsx \
        packages/HimalayaUI/frontend/test/TraceViewer.test.tsx
git commit -m "feat(ui): TraceViewer renders auto + manual peak markers"
```

---

## Task 7: Click-to-add peak (with snap-to-local-max)

**Motivation:** Let the user add a manual peak by clicking on the trace. Snap the click to the nearest local maximum within a small window so clicks don't land on noise.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx`
- Modify: `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx`

- [ ] **Step 1: Write a failing test for the snap helper**

The snap function is pure and easy to unit-test without mocking Plot. Export it from `TraceViewer.tsx` and test directly.

Add to `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx`:

```tsx
import { snapToLocalMax } from "../src/components/TraceViewer";

describe("snapToLocalMax", () => {
  const q = [0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16];
  const I = [  5,    6,   15,   14,   13,   20,    4];

  it("snaps to the local max within ±3 indices of the click", () => {
    // click at 0.11 -> window covers indices 0..4; local max is at 0.12 (I=15)
    expect(snapToLocalMax(q, I, 0.11, 3)).toBeCloseTo(0.12);
  });

  it("defaults K=3 when K is not provided", () => {
    expect(snapToLocalMax(q, I, 0.11)).toBeCloseTo(0.12);
  });

  it("clips to array bounds", () => {
    expect(snapToLocalMax(q, I, 0.99, 3)).toBeCloseTo(q[q.length - 1]!);
  });
});
```

- [ ] **Step 2: Run test to verify failure**

Run: `npm test -- TraceViewer`
Expected: the 3 new snap tests fail (`snapToLocalMax` not exported).

- [ ] **Step 3: Add `snapToLocalMax` and wire the click handler**

Add at module scope in `TraceViewer.tsx`:

```tsx
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
```

Update the destructure at the top of the component — drop the underscore aliases and the `void` discard lines added in Task 5. Final destructure:

```tsx
export function TraceViewer({
  trace, peaks, onAddPeak, onRemovePeak,
}: TraceViewerProps): JSX.Element {
```

(`onRemovePeak` is still unused until Task 8 — leave it in the destructure; TS doesn't error on unused function parameters, only unused locals. Include it in the effect deps array `[trace, peaks, onAddPeak, onRemovePeak]`.)

Inside the `useEffect`, after `host.replaceChildren(el);` and before the cleanup `return`, add:

```tsx
function handleClick(ev: MouseEvent): void {
  const xScale = (el as unknown as { scale: (name: string) => { invert?: (v: number) => number } | undefined }).scale("x");
  if (!xScale?.invert) return;
  const rect = el.getBoundingClientRect();
  const qClick = xScale.invert(ev.clientX - rect.left);
  const qSnapped = snapToLocalMax(trace.q, trace.I, qClick);
  onAddPeak(qSnapped);
}
(el as unknown as EventTarget).addEventListener("click", handleClick);
```

Update the cleanup return to also remove the listener:

```tsx
return () => {
  (el as unknown as EventTarget).removeEventListener("click", handleClick);
  host.replaceChildren();
};
```

- [ ] **Step 4: Run tests to verify pass**

Run: `npm test -- TraceViewer`
Expected: all 6 tests pass (the 3 existing + 3 new snap tests). The click-handler integration itself is covered by E2E in Task 10 — the unit test for Plot's click DOM event is brittle.

- [ ] **Step 5: Wire `useAddPeak` mutation into `App.tsx`**

Replace `packages/HimalayaUI/frontend/src/App.tsx`:

```tsx
import { useEffect, useState } from "react";
import "./styles.css";
import { useAppState } from "./state";
import {
  useExperiment, useSamples, useExposures, useTrace, usePeaks,
  useAddPeak, useRemovePeak,
} from "./queries";
import { Navbar } from "./components/Navbar";
import { Layout } from "./components/Layout";
import { SampleList } from "./components/SampleList";
import { UserModal } from "./components/UserModal";
import { ExposureList } from "./components/ExposureList";
import { TraceViewer } from "./components/TraceViewer";

const EXPERIMENT_ID = 1;

export function App(): JSX.Element {
  const username          = useAppState((s) => s.username);
  const activeSampleId    = useAppState((s) => s.activeSampleId);
  const activeExposureId  = useAppState((s) => s.activeExposureId);
  const setUsername       = useAppState((s) => s.setUsername);
  const setActiveSample   = useAppState((s) => s.setActiveSample);
  const setActiveExposure = useAppState((s) => s.setActiveExposure);

  const [modalOpen, setModalOpen] = useState<boolean>(!username);

  const experimentQ = useExperiment(EXPERIMENT_ID);
  const samplesQ    = useSamples(EXPERIMENT_ID);
  const exposuresQ  = useExposures(activeSampleId);
  const traceQ      = useTrace(activeExposureId);
  const peaksQ      = usePeaks(activeExposureId);

  const addPeak    = useAddPeak(activeExposureId ?? 0);
  const removePeak = useRemovePeak(activeExposureId ?? 0);

  // Auto-select the first exposure when samples or the active sample changes
  useEffect(() => {
    const exposures = exposuresQ.data ?? [];
    if (exposures.length === 0) return;
    const stillValid = exposures.some((e) => e.id === activeExposureId);
    if (!stillValid) setActiveExposure(exposures[0]!.id);
  }, [exposuresQ.data, activeExposureId, setActiveExposure]);

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
          traceQ.data && peaksQ.data && activeExposureId !== undefined ? (
            <TraceViewer
              trace={traceQ.data}
              peaks={peaksQ.data}
              onAddPeak={(q) => addPeak.mutate(q)}
              onRemovePeak={(peakId) => removePeak.mutate(peakId)}
            />
          ) : undefined
        }
        centerBottom={<ExposureList />}
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

- [ ] **Step 6: Update `smoke.test.tsx` mocks**

The smoke test was added in Task 3 of Plan 3.5 and mocks 3 endpoints. Extend it so boot doesn't fail with the new queries. Replace the `mockFetch({...})` call's argument in `packages/HimalayaUI/frontend/test/smoke.test.tsx` with:

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

(No active exposure → no trace/peaks/indices calls; the `useExposures(undefined)` hook is disabled until a sample is chosen, and with only 10 in the list and no user click the sample is still inactive at boot. If the smoke test previously set an active sample via persisted state, clear it — the test already calls `localStorage.clear()` in beforeEach.)

- [ ] **Step 7: Run all tests**

Run: `npm test`
Expected: all green.

- [ ] **Step 8: Build**

Run: `npm run build`
Expected: clean.

- [ ] **Step 9: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/TraceViewer.tsx \
        packages/HimalayaUI/frontend/test/TraceViewer.test.tsx \
        packages/HimalayaUI/frontend/src/App.tsx \
        packages/HimalayaUI/frontend/test/smoke.test.tsx
git commit -m "feat(ui): click on trace adds manual peak with snap-to-local-max"
```

---

## Task 8: Click a peak marker to remove it

**Motivation:** Close the editing loop — clicking near an existing peak removes it instead of adding a new one. We already have `useRemovePeak` wired through `App.tsx`.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx`
- Modify: `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx`

- [ ] **Step 1: Write failing tests for the nearest-peak helper**

Add to `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx`:

```tsx
import { findNearestPeak } from "../src/components/TraceViewer";

describe("findNearestPeak", () => {
  const peaks = [
    { id: 1, exposure_id: 1, q: 0.10, intensity: null, prominence: null,
      sharpness: null, source: "auto" as const },
    { id: 2, exposure_id: 1, q: 0.30, intensity: null, prominence: null,
      sharpness: null, source: "auto" as const },
  ];

  it("returns the closest peak when within tolerance", () => {
    expect(findNearestPeak(peaks, 0.105, 0.01)?.id).toBe(1);
  });

  it("returns null when outside tolerance", () => {
    expect(findNearestPeak(peaks, 0.20, 0.01)).toBeNull();
  });

  it("returns null for empty list", () => {
    expect(findNearestPeak([], 0.1, 1)).toBeNull();
  });
});
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- TraceViewer`
Expected: 3 new tests fail (`findNearestPeak` not exported).

- [ ] **Step 3: Add `findNearestPeak` and wire it into the click handler**

Add at module scope in `TraceViewer.tsx`:

```tsx
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
```

Replace the body of `handleClick` with:

```tsx
function handleClick(ev: MouseEvent): void {
  const xScale = (el as unknown as { scale: (name: string) => { invert?: (v: number) => number } | undefined }).scale("x");
  if (!xScale?.invert) return;
  const rect = el.getBoundingClientRect();
  const qClick = xScale.invert(ev.clientX - rect.left);

  // Tolerance: 2% of the click q (relative, log-friendly)
  const tolerance = Math.max(qClick * 0.02, 1e-6);
  const existing = findNearestPeak(peaks, qClick, tolerance);
  if (existing) onRemovePeak(existing.id);
  else          onAddPeak(snapToLocalMax(trace.q, trace.I, qClick));
}
```

- [ ] **Step 4: Run tests to verify pass**

Run: `npm test -- TraceViewer`
Expected: all 9 tests pass.

- [ ] **Step 5: Full suite + build**

Run: `npm test && npm run build`
Expected: clean.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/TraceViewer.tsx \
        packages/HimalayaUI/frontend/test/TraceViewer.test.tsx
git commit -m "feat(ui): clicking a peak marker removes it"
```

---

## Task 9: StaleIndicesBanner

**Motivation:** When manual peaks are added or removed, the backend marks affected indices `stale`. Surface that to the user with a banner offering re-analysis.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx`
- Create: `packages/HimalayaUI/frontend/test/StaleIndicesBanner.test.tsx`
- Modify: `packages/HimalayaUI/frontend/src/App.tsx` (render the banner above the trace)

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/StaleIndicesBanner.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, fireEvent, waitFor } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { StaleIndicesBanner } from "../src/components/StaleIndicesBanner";

beforeEach(() => { vi.restoreAllMocks(); });

function mockResp(url: string, status: number, body: unknown): void {
  vi.spyOn(global, "fetch").mockImplementation(async (input) => {
    const u = typeof input === "string" ? input : (input as Request).url;
    if (u.endsWith(url)) {
      return new Response(JSON.stringify(body), {
        status, headers: { "Content-Type": "application/json" },
      });
    }
    return new Response("not found", { status: 404 });
  });
}

describe("<StaleIndicesBanner>", () => {
  it("renders nothing when exposureId is undefined", () => {
    const { container } = renderWithProviders(
      <StaleIndicesBanner exposureId={undefined} />,
    );
    expect(container.textContent).toBe("");
  });

  it("renders nothing when no indices are stale", async () => {
    mockResp("/api/exposures/42/indices", 200, [
      { id: 1, exposure_id: 42, phase: "Pn3m", basis: 0.1, score: 1,
        r_squared: 0.99, lattice_d: 10, status: "candidate" },
    ]);
    const { container } = renderWithProviders(
      <StaleIndicesBanner exposureId={42} />,
    );
    await waitFor(() => expect(container.textContent).toBe(""));
  });

  it("renders a re-analyze button when any index is stale", async () => {
    mockResp("/api/exposures/42/indices", 200, [
      { id: 1, exposure_id: 42, phase: "Pn3m", basis: 0.1, score: 1,
        r_squared: 0.99, lattice_d: 10, status: "stale" },
    ]);
    renderWithProviders(<StaleIndicesBanner exposureId={42} />);
    await waitFor(() =>
      expect(screen.getByRole("button", { name: /re-analyze/i })).toBeInTheDocument(),
    );
  });

  it("clicking re-analyze calls POST /api/exposures/:id/analyze", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(async (input) => {
      const u = typeof input === "string" ? input : (input as Request).url;
      if (u.endsWith("/api/exposures/42/indices")) {
        return new Response(JSON.stringify([
          { id: 1, exposure_id: 42, phase: "Pn3m", basis: 0.1, score: 1,
            r_squared: 0.99, lattice_d: 10, status: "stale" },
        ]), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (u.endsWith("/api/exposures/42/analyze")) {
        return new Response(JSON.stringify({ id: 42, analyzed: true }),
          { status: 200, headers: { "Content-Type": "application/json" } });
      }
      return new Response("not found", { status: 404 });
    });
    renderWithProviders(<StaleIndicesBanner exposureId={42} />);
    const btn = await screen.findByRole("button", { name: /re-analyze/i });
    fireEvent.click(btn);
    await waitFor(() => {
      const urls = fetchSpy.mock.calls.map((c) =>
        typeof c[0] === "string" ? c[0] : (c[0] as Request).url);
      expect(urls).toContain("/api/exposures/42/analyze");
    });
  });
});
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- StaleIndicesBanner`
Expected: 4 tests fail (`StaleIndicesBanner` not defined).

- [ ] **Step 3: Implement `StaleIndicesBanner.tsx`**

Create `packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx`:

```tsx
import { useIndices, useReanalyzeExposure } from "../queries";

export interface StaleIndicesBannerProps {
  exposureId: number | undefined;
}

export function StaleIndicesBanner(
  { exposureId }: StaleIndicesBannerProps,
): JSX.Element | null {
  const q = useIndices(exposureId);
  // Always call the mutation hook so hook order is stable across renders.
  // When exposureId is undefined we never click the button anyway.
  const reanalyze = useReanalyzeExposure(exposureId ?? 0);

  if (exposureId === undefined) return null;
  const indices = q.data ?? [];
  const stale = indices.filter((i) => i.status === "stale");
  if (stale.length === 0) return null;

  return (
    <div
      role="alert"
      className="flex items-center justify-between gap-4 px-3 py-2 mb-2 border border-warning text-fg bg-bg-elevated rounded-md"
    >
      <span>
        {stale.length} {stale.length === 1 ? "index is" : "indices are"} stale.
      </span>
      <button
        className="bg-accent border border-accent text-white rounded-md px-2.5 py-1 hover:brightness-110 focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent"
        disabled={reanalyze.isPending}
        onClick={() => reanalyze.mutate()}
      >
        {reanalyze.isPending ? "Re-analyzing…" : "Re-analyze"}
      </button>
    </div>
  );
}
```

- [ ] **Step 4: Run to verify pass**

Run: `npm test -- StaleIndicesBanner`
Expected: 4/4 pass.

- [ ] **Step 5: Wire the banner into `App.tsx`**

In `packages/HimalayaUI/frontend/src/App.tsx`, import:

```tsx
import { StaleIndicesBanner } from "./components/StaleIndicesBanner";
```

And wrap the `centerTop` prop's value so the banner appears above the viewer:

```tsx
centerTop={
  <div className="flex flex-col flex-1 min-h-0">
    <StaleIndicesBanner exposureId={activeExposureId} />
    {traceQ.data && peaksQ.data && activeExposureId !== undefined ? (
      <TraceViewer
        trace={traceQ.data}
        peaks={peaksQ.data}
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
```

- [ ] **Step 6: Smoke test mock update**

The smoke test triggers `useIndices` whenever an active exposure is set. It currently doesn't set one, so no change is required for the existing smoke test. But to confirm: running `npm test` should still pass.

- [ ] **Step 7: Full suite + build**

Run: `npm test && npm run build`
Expected: all green.

- [ ] **Step 8: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx \
        packages/HimalayaUI/frontend/test/StaleIndicesBanner.test.tsx \
        packages/HimalayaUI/frontend/src/App.tsx
git commit -m "feat(ui): StaleIndicesBanner prompts re-analysis"
```

---

## Task 10: E2E happy path — click to add a manual peak

**Motivation:** Lock in the full flow in Playwright so Plan 5+ refactors can't silently break it.

**Files:**
- Modify: `packages/HimalayaUI/frontend/e2e/smoke.spec.ts`

- [ ] **Step 1: Add the new E2E test**

Append to `packages/HimalayaUI/frontend/e2e/smoke.spec.ts` (after the existing tests, still inside the same `test.describe`, if one exists; otherwise at top level):

```ts
test("clicking the trace posts a manual peak", async ({ page }) => {
  let posted: { q: number } | null = null;

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
      q:     [0.05, 0.10, 0.15, 0.20, 0.25, 0.30],
      I:     [ 100,   80,   90,  500,   85,   70],
      sigma: [  5,    4,    5,   25,    4,    3],
    }),
  }));
  await page.route("**/api/exposures/100/peaks", async (r) => {
    if (r.request().method() === "POST") {
      const body = r.request().postDataJSON() as { q: number };
      posted = body;
      return r.fulfill({
        status: 201,
        body: JSON.stringify({
          id: 999, exposure_id: 100, q: body.q, intensity: null,
          prominence: null, sharpness: null, source: "manual",
          stale_indices: 1,
        }),
      });
    }
    return r.fulfill({ status: 200, body: "[]" });
  });
  await page.route("**/api/exposures/100/indices", (r) =>
    r.fulfill({ status: 200, body: "[]" }));

  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({
        state: { username: "alice", activeSampleId: 10, activeExposureId: 100 },
        version: 1,
      }),
    );
  });

  await page.goto("/");
  const viewer = page.locator('[data-testid="trace-viewer"]');
  await viewer.waitFor();
  const box = await viewer.boundingBox();
  if (!box) throw new Error("trace viewer has no bounding box");

  // Click somewhere in the middle of the plot — exact coords don't matter,
  // the snap will land on one of the data point q values.
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);

  await expect.poll(() => posted?.q, { timeout: 2000 }).toBeGreaterThan(0);
});
```

- [ ] **Step 2: Build the frontend so the test runs against the actual bundle**

Run from `packages/HimalayaUI/frontend/`: `npm run build`
Expected: clean.

- [ ] **Step 3: Run E2E**

Run: `npm run e2e`
Expected: the new test passes along with the existing ones.

If the click doesn't register, it's almost certainly because Observable Plot's SVG intercepts the pointer event differently than a raw `<div>`. Diagnosis: inspect `viewer.boundingBox()` and, if needed, target the inner SVG explicitly with `viewer.locator("svg")`. Don't adjust the click-handling code without first confirming the DOM.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/e2e/smoke.spec.ts
git commit -m "test(e2e): clicking the trace posts a manual peak"
```

---

## Final verification

From repo root:

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: all backend tests pass (existing + 1 new trace-route testset).

From `packages/HimalayaUI/frontend/`:

```bash
npm test && npm run build && npm run e2e
```

Expected: all unit tests pass (old + ~19 new); tsc clean; Vite build clean; Playwright passes all existing E2E + 1 new.

Hand off to Plan 5 (Miller-index plot + Phase panel + active-group overlay + hover preview).
