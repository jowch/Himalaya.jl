# Phase-4 Plan 4a — `GET /api/series/:id/traces` + `useSeriesTraces` (the folio's batch-trace gate)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one backend read route that returns every member trace of a series as a single `exposure_id → Trace` map, plus the frontend `api.getSeriesTraces` + `useSeriesTraces` hook that consumes it — so the greenfield Series-folio page can render `SeriesCard`→`CardFigure` real-curve mini-waterfalls without an O(N×M) per-exposure fetch fan-out.

**Architecture:** A new `@get "/api/series/{id}/traces"` handler in `routes_series.jl` joins `series_members → exposures → samples → experiments`, resolves each member's `.dat` path the same way the single-exposure trace route does, loads it via the shared `load_dat`, and returns a JSON object `{ "<exposure_id>": {q,I,sigma} }` — **skipping** (not failing on) members with no exposure, derived exposures, or a missing `.dat`. The frontend mirrors the existing `getTrace`/`useTrace` pattern; the hook returns `Record<number, Trace>` — the exact shape `toWaterfallRows(members, tracesById)` already indexes, so no client-side conversion.

**Tech Stack:** Julia + Oxygen.jl 1.10 (`@get`, `current_db()`, `JSON3.write`), SQLite.jl/DBInterface, `Test` + in-process `with_test_server`; React 18 + TypeScript (`exactOptionalPropertyTypes: true`) + TanStack Query, Vitest + RTL.

---

## Standing constraints (do not violate)

- **Provenance (from `docs/superpowers/specs/2026-06-03-phase4-cutover-strategy-design.md`):** this plan touches only **CARRIED** infra (`api.ts`, `queries.ts`) + **backend** (`routes_series.jl`, `series.jl`). No `src/print/**` page is built here (that is Plan 4b). No legacy `src/pages/**` or `src/components/**` is imported or edited.
- Commit ONLY specifically-named files. **NEVER** `git add -A` / `git add .`.
- **NEVER** stage `src/bones/registry.ts`, `src/bones/*.bones.json`. (Plan-doc staging policy is unresolved on this branch — see *Plan-durability note* at the foot; default: leave this plan untracked unless told otherwise.)
- Every commit's exact last line: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
- **Frontend** work runs from `packages/HimalayaUI/frontend`. Typecheck: `npx tsc --noEmit -p tsconfig.build.json`. Design guard (`npm run lint:design`) is N/A here — no `.tsx`/appearance touched, but run it anyway for a clean gate (it must stay exit 0).
- **Backend** work runs from repo root with `--project=packages/HimalayaUI`. The Julia suite is slow (5–10 min); see *Backend test loop* below — do **not** re-run it repeatedly with different greps.
- Dispatch implementers SEQUENTIALLY (avoid `.git/index.lock` collisions). Tasks here are ordered backend → frontend; the frontend tasks do not need the backend running (they mock `fetch`).

### Backend test loop (how to run Task 1's test without an 8-minute wait every step)

The canonical gate is the full suite, captured once:
```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|did not pass|fail|traces" /tmp/jl-test.out
```
For a faster inner loop, a standalone driver that includes only the harness + this one file usually works (fall back to the full suite if it errors on a missing `using`):
```bash
julia --project=packages/HimalayaUI -e '
using Test, HimalayaUI, HTTP, JSON3, SQLite, DBInterface, Tables
include("packages/HimalayaUI/test/test_http.jl")          # defines with_test_server
include("packages/HimalayaUI/test/test_routes_series.jl")  # the file under edit
'
```

---

## Reference: verified APIs (mapped from live source 2026-06-06 — confirm line numbers at execution; they drift)

### Backend (`packages/HimalayaUI/src/`)
```julia
# routes_series.jl (415 lines) — existing read-route handler shape (lines 80-85):
@get "/api/series/{id}" function(req::HTTP.Request, id::Int)
    db = current_db()
    out = fetch_series_with_plate(db, id)
    out === nothing && return _json_error(404, "series not found")
    HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
end
# Add the new route AFTER the `/forks` route (~line 91), BEFORE the first `@post` (~line 93).
# _json_error(status, msg) lives in json.jl (~lines 39-46). current_db() is the request-scoped DB.

# routes_trace.jl (lines 4-35) — how ONE exposure's trace path is resolved + loaded:
#   join exposures→samples→experiments for (filename, kind, experiment_id, analysis_dir);
#   require kind == "file" and filename !== missing;
#   cfg = config_from_db(db, experiment_id); path = joinpath(analysis_dir, replace(cfg.integration_pattern, "{name}" => filename));
#   isfile(path); q, I, σ = load_dat(path);  returns Dict(:q=>q, :I=>I, :sigma=>σ)

# datfile.jl (lines 10-14):
load_dat(path::AbstractString) -> (q::Vector{Float64}, I::Vector{Float64}, σ::Vector{Float64})

# config_from_db(db, experiment_id::Int) -> cfg  (cfg.integration_pattern is the "{name}"-templated filename)

# series_members schema (db.jl ~787): columns incl. series_id, exposure_id (nullable), display_order.
```

### Frontend (`packages/HimalayaUI/frontend/src/`)
```ts
// api.ts (186-193) — Trace + the single-trace fetch to mirror:
export interface Trace { q: number[]; I: number[]; sigma: number[]; }
export const getTrace = (exposure_id: number) =>
  request<Trace>("GET", `/api/exposures/${exposure_id}/trace`);
// request<T>(method, path, body?, opts?) is the shared fetch wrapper (api.ts 68-95).

// queries.ts — keys + the single-trace hook to mirror:
trace: (exposureId: number | undefined) => ["exposure", exposureId ?? "none", "trace"] as const,   // line 58-59
series: (id: number | undefined) => ["series", id ?? "none"] as const,                              // existing
export function useTrace(exposureId: number | undefined) {                                          // 244-250
  return useQuery({ queryKey: queryKeys.trace(exposureId), queryFn: () => api.getTrace(exposureId as number), enabled: exposureId !== undefined });
}

// THE CONSUMER CONTRACT (src/print/waterfall/waterfallModel.ts 46-54) — do not change it; match it:
export function toWaterfallRows(members: SeriesMember[], tracesById: Record<number, Trace>): WaterfallRow[]
//   row.trace = (member.exposure_id != null && tracesById[member.exposure_id]) || EMPTY_TRACE;  ← plain object, number index, EMPTY_TRACE fallback
```

### Test patterns
```julia
# test/test_routes_trace.jl (3-27) — the PROVEN .dat fixture (copy this setup verbatim, then add series rows):
analysis_dir = joinpath(tmp, "analysis", "automatic_analysis"); mkpath(analysis_dir)
cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"), joinpath(analysis_dir, "example_tot.dat"))
db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
exp_id = HimalayaUI.init_experiment!(db; path=tmp, data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
with_test_server(db) do port, base
    r = HTTP.get("$base/api/exposures/$e_id/trace"); @test r.status == 200
    body = JSON3.read(String(r.body)); @test haskey(body, :q) && haskey(body, :I) && haskey(body, :sigma)
end
```
```tsx
// frontend test/queries.test.tsx (61-67) — the read-hook test shape to mirror:
it("useTrace fetches for a given exposureId", async () => {
  mockOnce(200, { q: [0.1], I: [10], sigma: [1] });
  const { wrapper } = withClient();
  const { result } = renderHook(() => useTrace(42), { wrapper });
  await waitFor(() => expect(result.current.isSuccess).toBe(true));
  expect(result.current.data?.q).toEqual([0.1]);
});
// withClient()/mockOnce() helpers live at the top of that file; reuse them.
```

---

## File map
- **Modify** `packages/HimalayaUI/src/series.jl` — add the internal helper `series_member_traces(db, series_id) -> Dict{Int,Any}` (resolve + load each member's trace, skip the unresolvable).
- **Modify** `packages/HimalayaUI/src/routes_series.jl` — add `@get "/api/series/{id}/traces"` (existence check → 404, else the helper's map).
- **Modify** `packages/HimalayaUI/test/test_routes_series.jl` — add the route testset (happy path + skip cases + 404).
- **Modify** `packages/HimalayaUI/frontend/src/api.ts` — add `getSeriesTraces`.
- **Modify** `packages/HimalayaUI/frontend/src/queries.ts` — add `queryKeys.seriesTraces` + `useSeriesTraces`.
- **Create** `packages/HimalayaUI/frontend/test/queries-series-traces.test.tsx` — the hook test.

> **DRY note (deliberate non-refactor):** the single-exposure trace route (`routes_trace.jl`) resolves a `.dat` path with near-identical join+config logic. We do **not** refactor it to share a loader in this plan — its 400-vs-404 error semantics differ from the batch route's skip-and-continue semantics, and touching tested route code widens blast radius for no folio benefit. `load_dat` (the actual file parser) **is** the shared primitive. A later dedup is a deferred nicety, noted here so it is not re-discovered.

---

## Task 1: Backend — `series_member_traces` helper + the route (TDD)

**Files:**
- Modify: `packages/HimalayaUI/src/series.jl` (add helper)
- Modify: `packages/HimalayaUI/src/routes_series.jl` (add route ~after line 91)
- Test: `packages/HimalayaUI/test/test_routes_series.jl` (add testset)

- [ ] **Step 1: Write the failing test** — append to `test/test_routes_series.jl`:

```julia
@testset "GET /api/series/{id}/traces" begin
    mktempdir() do tmp
        # Proven .dat fixture (mirrors test_routes_trace.jl): one resolvable exposure.
        analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
        mkpath(analysis_dir)
        cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
           joinpath(analysis_dir, "example_tot.dat"))
        db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
        exp_id = HimalayaUI.init_experiment!(db; path=tmp,
            data_dir=joinpath(tmp, "data"), analysis_dir=analysis_dir)
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
        good   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
        # A second exposure whose .dat does NOT exist on disk → must be SKIPPED, not 500.
        missing_dat = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="nope")

        DBInterface.execute(db, "INSERT INTO series (id, title, state) VALUES (7, 'S7', 'draft')")
        # display_order 0 = good, 1 = missing-dat, 2 = NULL exposure (orphan) → both skipped.
        DBInterface.execute(db, """INSERT INTO series_members (series_id, exposure_id, display_order, created_at)
            VALUES (7, $good, 0, '2026-06-06T00:00:00.000Z')""")
        DBInterface.execute(db, """INSERT INTO series_members (series_id, exposure_id, display_order, created_at)
            VALUES (7, $missing_dat, 1, '2026-06-06T00:00:00.000Z')""")
        DBInterface.execute(db, """INSERT INTO series_members (series_id, exposure_id, display_order, created_at)
            VALUES (7, NULL, 2, '2026-06-06T00:00:00.000Z')""")

        with_test_server(db) do port, base
            # 404 for an unknown series.
            r404 = HTTP.get("$base/api/series/999/traces", ["X-Username" => "alice"];
                            status_exception = false)
            @test r404.status == 404

            r = HTTP.get("$base/api/series/7/traces", ["X-Username" => "alice"])
            @test r.status == 200
            body = JSON3.read(String(r.body), Dict{String, Any})
            # Only the resolvable exposure is present; the missing-dat + NULL members are skipped.
            @test collect(keys(body)) == [string(good)]
            tr = body[string(good)]
            @test haskey(tr, "q") && haskey(tr, "I") && haskey(tr, "sigma")
            @test length(tr["q"]) == length(tr["I"]) == length(tr["sigma"]) > 0

            # An existing series with zero resolvable members → 200 + empty object (not 404).
            DBInterface.execute(db, "INSERT INTO series (id, title, state) VALUES (8, 'S8', 'draft')")
            rEmpty = HTTP.get("$base/api/series/8/traces", ["X-Username" => "alice"])
            @test rEmpty.status == 200
            @test JSON3.read(String(rEmpty.body), Dict{String, Any}) == Dict{String, Any}()
        end
        close(db)
    end
end
```

- [ ] **Step 2: Run → fail** — see *Backend test loop*. Expected: `UndefVarError`/404-on-known-series (route not registered yet) under the new testset.

- [ ] **Step 3: Implement the helper** — append to `packages/HimalayaUI/src/series.jl`:

```julia
"""
    series_member_traces(db, series_id) -> Dict{Int,Any}

Resolve and load every member trace of `series_id`, keyed by exposure_id, in
display order. Members with no exposure, a derived (non-"file") exposure, no
filename, or a missing `.dat` on disk are SKIPPED (omitted from the map) — the
batch route degrades gracefully; the folio's `toWaterfallRows` renders an empty
row for any absent exposure. `config_from_db` is memoized per experiment.
"""
function series_member_traces(db::SQLite.DB, series_id::Integer)
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT e.id AS exposure_id, e.filename, e.kind,
                  x.id AS experiment_id, x.analysis_dir
           FROM series_members sm
           JOIN exposures e   ON e.id = sm.exposure_id
           JOIN samples s     ON s.id = e.sample_id
           JOIN experiments x ON x.id = s.experiment_id
           WHERE sm.series_id = ? AND sm.exposure_id IS NOT NULL
           ORDER BY sm.display_order ASC, sm.id ASC""", [series_id]))

    out = Dict{Int,Any}()
    cfg_cache = Dict{Int,Any}()
    for row in rows
        String(row.kind) == "file" || continue
        row.filename === missing && continue
        eid = Int(row.experiment_id)
        cfg = get!(cfg_cache, eid) do
            config_from_db(db, eid)
        end
        pattern = replace(cfg.integration_pattern, "{name}" => String(row.filename))
        path    = joinpath(String(row.analysis_dir), pattern)
        isfile(path) || continue
        q, I, σ = load_dat(path)
        out[Int(row.exposure_id)] = Dict(:q => q, :I => I, :sigma => σ)
    end
    return out
end
```

> If `series.jl` does not already `using Tables`/`SQLite`/`DBInterface` at module scope, they are imported by `HimalayaUI.jl` for the package — reference them unqualified as the existing `fetch_series_with_plate` (same file) does. Match that file's existing qualification style exactly.

- [ ] **Step 4: Implement the route** — in `packages/HimalayaUI/src/routes_series.jl`, immediately after the `/forks` handler (~line 91):

```julia
@get "/api/series/{id}/traces" function(req::HTTP.Request, id::Int)
    db = current_db()
    exists = !isempty(Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM series WHERE id = ?", [id])))
    exists || return _json_error(404, "series not found")
    out = series_member_traces(db, id)
    HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
end
```

- [ ] **Step 5: Run → pass** — see *Backend test loop*. Expected: the new testset passes (key set `[good]`, q/I/sigma present, empty series → `{}`, unknown series → 404). Confirm no OTHER series testset regressed in `/tmp/jl-test.out`.

- [ ] **Step 6: Commit**
```bash
git add packages/HimalayaUI/src/series.jl packages/HimalayaUI/src/routes_series.jl packages/HimalayaUI/test/test_routes_series.jl
git commit -m "$(printf 'Phase-4: GET /api/series/:id/traces — batch member-trace map for the folio\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 2: Frontend — `api.getSeriesTraces`

**Files:** Modify `packages/HimalayaUI/frontend/src/api.ts` (beside `getTrace`, ~line 193).

- [ ] **Step 1: Implement** — add directly after `getTrace`:

```ts
/** Batch member traces for a series, keyed by exposure_id. Matches the
 *  `toWaterfallRows(members, tracesById)` contract (a plain Record, number index).
 *  Unresolvable members (no exposure / derived / missing .dat) are absent from the map. */
export const getSeriesTraces = (series_id: number) =>
  request<Record<number, Trace>>("GET", `/api/series/${series_id}/traces`);
```

- [ ] **Step 2: Typecheck**
```bash
cd packages/HimalayaUI/frontend && npx tsc --noEmit -p tsconfig.build.json
```
Expected: exit 0 (no callers yet; this just adds an export).

- [ ] **Step 3: Commit**
```bash
git add packages/HimalayaUI/frontend/src/api.ts
git commit -m "$(printf 'Phase-4: api.getSeriesTraces — batch series trace fetch\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 3: Frontend — `useSeriesTraces` hook + queryKey (TDD)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/queries.ts` (add key beside `series`; add hook beside `useTrace`/`useSeries`)
- Test: `packages/HimalayaUI/frontend/test/queries-series-traces.test.tsx`

- [ ] **Step 1: Write the failing test** — create `test/queries-series-traces.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { useSeriesTraces } from "../src/queries";

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
      status,
      headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("queries — useSeriesTraces", () => {
  beforeEach(() => vi.restoreAllMocks());

  it("fetches the exposure_id → Trace map for a series", async () => {
    mockOnce(200, { 1000: { q: [0.1], I: [10], sigma: [1] } });
    const { wrapper } = withClient();
    const { result } = renderHook(() => useSeriesTraces(7), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    // JSON object keys are strings at runtime; number-index access resolves them,
    // exactly as toWaterfallRows does (tracesById[member.exposure_id]).
    expect(result.current.data?.[1000]?.q).toEqual([0.1]);
  });

  it("is disabled (does not fetch) when seriesId is undefined", () => {
    const fetchSpy = vi.spyOn(global, "fetch");
    const { wrapper } = withClient();
    const { result } = renderHook(() => useSeriesTraces(undefined), { wrapper });
    expect(result.current.fetchStatus).toBe("idle");
    expect(fetchSpy).not.toHaveBeenCalled();
  });
});
```

- [ ] **Step 2: Run → fail**
```bash
cd packages/HimalayaUI/frontend && npm test -- queries-series-traces
```
Expected: FAIL — `useSeriesTraces` is not exported from `../src/queries`.

- [ ] **Step 3: Implement** — in `src/queries.ts`, add the key beside `series:` in the `queryKeys` object (mirrors the `trace` nesting under its entity):

```ts
  seriesTraces: (id: number | undefined) =>
    ["series", id ?? "none", "traces"] as const,
```

and add the hook beside `useSeries`:

```ts
/** Batch member traces for a series (one request, `exposure_id → Trace`). Feeds
 *  the greenfield folio's `SeriesCard`→`CardFigure` via `toWaterfallRows`; avoids
 *  the O(N×M) per-exposure fan-out the legacy folio would incur on greenfield curves. */
export function useSeriesTraces(seriesId: number | undefined) {
  return useQuery({
    queryKey: queryKeys.seriesTraces(seriesId),
    queryFn: () => api.getSeriesTraces(seriesId as number),
    enabled: seriesId !== undefined,
  });
}
```

- [ ] **Step 4: Run → pass**
```bash
cd packages/HimalayaUI/frontend && npm test -- queries-series-traces
```
Expected: both tests PASS.

- [ ] **Step 5: Typecheck**
```bash
cd packages/HimalayaUI/frontend && npx tsc --noEmit -p tsconfig.build.json
```
Expected: exit 0.

- [ ] **Step 6: Commit**
```bash
git add packages/HimalayaUI/frontend/src/queries.ts packages/HimalayaUI/frontend/test/queries-series-traces.test.tsx
git commit -m "$(printf 'Phase-4: useSeriesTraces hook + seriesTraces queryKey\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 4: Final gate

- [ ] **Step 1: Frontend gate**
```bash
cd packages/HimalayaUI/frontend
npx tsc --noEmit -p tsconfig.build.json   # exit 0
npm run lint:design                       # exit 0 (no appearance touched)
npm test -- queries-series-traces queries # the new + the sibling trace tests stay green
```

- [ ] **Step 2: Backend gate** — the full HimalayaUI suite green, captured once:
```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-test.out   # the suite passes; new traces testset included
```

- [ ] **Step 3: Manual smoke (optional, if a live backend is handy)** — `bin/himalaya serve <experiment>` then `curl -s localhost:8080/api/series/<id>/traces | head -c 200` returns `{"<exposure_id>":{"q":[…]…}}`. Confirms the route is wired through `register_series_routes!`.

No route-registration edit is needed beyond adding the `@get` inside `routes_series.jl` — `register_series_routes!` already wraps every `@get`/`@post` in that file (it is called from `register_routes!` in `server.jl`). If the new route 404s at runtime, confirm the `@get` is inside the `register_series_routes!` function body, not at module top level.

---

## Self-review checklist (run before declaring done)

- **Spec coverage** (strategy step 4, backend half): the route returns one `exposure_id → Trace` map per series ✅; the `useSeriesTraces` hook + `getSeriesTraces` fetch consume it ✅; folio assembly (the `toWaterfallRows` wiring) is **Plan 4b**, not here. The 409-relax / dead-code / other strategy steps are out of scope.
- **Contract match:** hook returns `Record<number, Trace>` — the exact `toWaterfallRows(members, tracesById)` parameter type (not `Map`); route emits a JSON object with string keys, number-indexed on the client. ✅
- **Skip-don't-fail:** NULL exposure_id, derived (`kind != "file"`), no filename, and missing `.dat` are all `continue`d; an existing series with zero resolvable members returns `200 {}`, only an unknown series id returns `404`. Tested. ✅
- **Type/name consistency:** `getSeriesTraces` (api) ↔ `api.getSeriesTraces` (hook) ↔ `series_member_traces` (Julia helper) ↔ `seriesTraces` (queryKey) — names match across tasks. `Trace` reused, not redefined. ✅
- **No placeholders:** every step has full code + an exact command + expected result. ✅
- **Verify-at-execution:** `routes_series.jl` line numbers (80-85 handler, ~91 insertion point) and `series.jl` qualification style drift — confirm against live source per the reference note. ✅

---

## Plan-durability note (carry to the human)

The Phase-4 plans are inconsistently versioned: `phase4-samples-cutover` was committed (`2485304`) while `phase4-focus`/`phase4-loupe-pilot` exist only as **untracked files inside the gitignored `.claude/worktrees/` tree** — if that worktree is cleaned, those specs vanish. This plan follows the untracked default for now. **Recommendation:** make one deliberate decision — either commit the Phase-4 plan docs (treat them as durable artifacts; this contradicts the legacy "never stage plans" guardrail, so it is a convention change the human should bless) or relocate them outside the disposable worktree. Resolve before the worktree churns.
