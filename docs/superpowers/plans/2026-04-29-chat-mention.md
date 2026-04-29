# Chat @-mention System Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an @-tag mention system to the Himalaya chat that lets users reference peaks, indices, exposures, and samples inline — with hover tooltips, Zustand-wired highlights, and smooth click navigation.

**Architecture:** `[[type:id]]` tokens are stored verbatim in `sample_messages.body TEXT` (no schema change). At render time, `renderMentions()` splits body text on tokens and maps each to a `MentionChip` resolved via `useMentionResolution` (client cache first, lazy fetch fallback). The `MentionCompose` wrapper intercepts `@` in the textarea and shows a `MentionPicker` dropdown.

**Tech Stack:** React 18 + TypeScript strict, TanStack Query v5, Zustand, `@chenglou/pretext`, Oxygen.jl (Julia backend), Vitest + RTL, Playwright E2E.

---

## File Map

**New frontend:**
- `packages/HimalayaUI/frontend/src/lib/renderMentions.tsx` — pure token parser/renderer
- `packages/HimalayaUI/frontend/src/hooks/useMentionResolution.ts` — hybrid resolution hook
- `packages/HimalayaUI/frontend/src/components/MentionChip.tsx` — inline chip + tooltip + Zustand wiring
- `packages/HimalayaUI/frontend/src/components/MentionPicker.tsx` — autocomplete dropdown
- `packages/HimalayaUI/frontend/src/components/MentionCompose.tsx` — compose wrapper with @ trigger

**New tests:**
- `packages/HimalayaUI/frontend/test/renderMentions.test.tsx`
- `packages/HimalayaUI/frontend/test/MentionChip.test.tsx`
- `packages/HimalayaUI/frontend/test/useMentionResolution.test.ts`
- `packages/HimalayaUI/frontend/test/MentionPicker.test.tsx`
- `packages/HimalayaUI/frontend/test/MentionCompose.test.tsx`
- `packages/HimalayaUI/test/test_routes_mentions.jl`

**Modified:**
- `packages/HimalayaUI/frontend/src/state.ts` — add `hoveredPeakId` + `setHoveredPeak`
- `packages/HimalayaUI/frontend/src/api.ts` — add `getPeak`, `getIndex`, `getExposure`, `getSample` fetchers
- `packages/HimalayaUI/frontend/src/queries.ts` — add `queryKeys.peak/index/exposure/sample` + hooks
- `packages/HimalayaUI/frontend/src/components/ChatCard.tsx` — integrate MentionCompose + renderMentions
- `packages/HimalayaUI/frontend/test/ChatCard.test.tsx` — update affected tests
- `packages/HimalayaUI/frontend/test/state.test.ts` — add hoveredPeakId tests
- `packages/HimalayaUI/src/routes_peaks.jl` — add `GET /api/peaks/:id`
- `packages/HimalayaUI/src/routes_exposures.jl` — add `GET /api/exposures/:id`
- `packages/HimalayaUI/src/routes_samples.jl` — add `GET /api/samples/:id`
- `packages/HimalayaUI/src/routes_analysis.jl` — add `GET /api/indices/:id`
- `packages/HimalayaUI/test/runtests.jl` — include new test file
- `packages/HimalayaUI/frontend/package.json` — add `@chenglou/pretext`

**Note on `ngc`:** The `ngc` field is not currently in the `IndexEntry` response. The tooltip implements the slot as a conditional — if `IndexEntry` gains an `ngc: number | null` field in the future, add it to the `GET /api/exposures/:id/indices` query and `IndexEntry` type, and the chip tooltip will display it automatically.

---

## Task 1: Add `hoveredPeakId` to Zustand state

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/state.ts`
- Modify: `packages/HimalayaUI/frontend/test/state.test.ts`

- [ ] **Step 1: Write failing tests**

Add to `packages/HimalayaUI/frontend/test/state.test.ts` inside the `describe("useAppState")` block:

```typescript
it("hoveredPeakId starts undefined and can be set/cleared", () => {
  useAppState.setState({ hoveredPeakId: undefined });
  expect(useAppState.getState().hoveredPeakId).toBeUndefined();
  useAppState.getState().setHoveredPeak(5);
  expect(useAppState.getState().hoveredPeakId).toBe(5);
  useAppState.getState().setHoveredPeak(undefined);
  expect(useAppState.getState().hoveredPeakId).toBeUndefined();
});

it("hoveredPeakId is NOT in the persisted partition", () => {
  useAppState.setState({ hoveredPeakId: 99 });
  const raw = localStorage.getItem(LS_KEY) ?? "";
  expect(raw).not.toContain("hoveredPeakId");
});

it("does NOT persist ephemeral UI fields (navModal*, hoveredIndexId, hoveredPeakId)", () => {
  useAppState.setState({ hoveredIndexId: 3, hoveredPeakId: 5, navModalOpen: true, navModalStep: "sample" });
  const raw = localStorage.getItem(LS_KEY) ?? "";
  expect(raw).not.toContain("hoveredPeakId");
  expect(raw).not.toContain("navModalOpen");
});
```

- [ ] **Step 2: Run tests — expect failure**

```bash
cd packages/HimalayaUI/frontend
node_modules/.bin/vitest run test/state.test.ts
```

Expected: FAIL — `hoveredPeakId` not in state type.

- [ ] **Step 3: Add `hoveredPeakId` to state**

In `packages/HimalayaUI/frontend/src/state.ts`, add to `AppState` interface (in the ephemeral section):

```typescript
hoveredPeakId: number | undefined;
```

Add setter to the interface:

```typescript
setHoveredPeak: (id: number | undefined) => void;
```

Add to the `create` implementation (after `hoveredIndexId: undefined`):

```typescript
hoveredPeakId: undefined,
```

Add setter (after `setHoveredIndex`):

```typescript
setHoveredPeak: (hoveredPeakId) => set({ hoveredPeakId }),
```

The `partialize` block is unchanged — `hoveredPeakId` is ephemeral and must NOT appear there.

- [ ] **Step 4: Run tests — expect pass**

```bash
node_modules/.bin/vitest run test/state.test.ts
```

Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/state.ts packages/HimalayaUI/frontend/test/state.test.ts
git commit -m "feat(ui): add hoveredPeakId ephemeral state for mention chip wiring"
```

---

## Task 2: Backend — four new entity GET routes

**Files:**
- Modify: `packages/HimalayaUI/src/routes_peaks.jl`
- Modify: `packages/HimalayaUI/src/routes_exposures.jl`
- Modify: `packages/HimalayaUI/src/routes_samples.jl`
- Modify: `packages/HimalayaUI/src/routes_analysis.jl`
- Create: `packages/HimalayaUI/test/test_routes_mentions.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl`

The existing routes are all collection or action routes. We need single-entity GET routes for lazy mention resolution.

- [ ] **Step 1: Add test file**

Create `packages/HimalayaUI/test/test_routes_mentions.jl`:

```julia
using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "mention lookup routes" begin
    tmp = mktempdir()
    db  = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="A", name="sampleA")
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="run001")

    # Insert a peak manually
    res = DBInterface.execute(db,
        "INSERT INTO peaks (exposure_id, q, source, excluded) VALUES (?, ?, 'auto', 0)",
        [e_id, 1.223])
    pk_id = Int(DBInterface.lastrowid(res))

    with_test_server(db) do port, base
        @testset "GET /api/peaks/:id" begin
            r = HTTP.get("$base/api/peaks/$pk_id")
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test body.id == pk_id
            @test body.q ≈ 1.223
            @test body.source == "auto"
            @test body.excluded === false

            r404 = HTTP.get("$base/api/peaks/999999"; status_exception=false)
            @test r404.status == 404
        end

        @testset "GET /api/exposures/:id" begin
            r = HTTP.get("$base/api/exposures/$e_id")
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test body.id == e_id
            @test body.filename == "run001"

            r404 = HTTP.get("$base/api/exposures/999999"; status_exception=false)
            @test r404.status == 404
        end

        @testset "GET /api/samples/:id" begin
            r = HTTP.get("$base/api/samples/$s_id")
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test body.id == s_id
            @test body.name == "sampleA"

            r404 = HTTP.get("$base/api/samples/999999"; status_exception=false)
            @test r404.status == 404
        end

        @testset "GET /api/indices/:id" begin
            # No indices until analysis runs — just test 404 path
            r404 = HTTP.get("$base/api/indices/999999"; status_exception=false)
            @test r404.status == 404
        end
    end
end
```

- [ ] **Step 2: Include test in runtests.jl**

In `packages/HimalayaUI/test/runtests.jl`, add:

```julia
include("test_routes_mentions.jl")
```

- [ ] **Step 3: Run tests — expect failure**

```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: FAIL — routes not defined yet.

- [ ] **Step 4: Add `GET /api/peaks/:id`**

In `packages/HimalayaUI/src/routes_peaks.jl`, inside `register_peaks_routes!()`, add before the closing `end`:

```julia
@get "/api/peaks/{id}" function(req::HTTP.Request, id::Int)
    db   = current_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, exposure_id, q, intensity, prominence, sharpness, source, excluded
         FROM peaks WHERE id = ?", [id]))
    isempty(rows) && return HTTP.Response(404,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "peak not found")))
    HTTP.Response(200, ["Content-Type" => "application/json"],
        JSON3.write(row_to_json(rows[1]; bool_keys = (:excluded,))))
end
```

- [ ] **Step 5: Add `GET /api/exposures/:id`**

In `packages/HimalayaUI/src/routes_exposures.jl`, inside `register_exposures_routes!()`, add before the closing `end`:

```julia
@get "/api/exposures/{id}" function(req::HTTP.Request, id::Int)
    db   = current_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM exposures WHERE id = ?", [id]))
    isempty(rows) && return HTTP.Response(404,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "exposure not found")))
    ex   = rows[1]
    tags = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, key, value, source FROM exposure_tags
         WHERE exposure_id = ? ORDER BY id", [id]))
    d        = row_to_json(ex; bool_keys = (:selected,))
    d[:tags] = rows_to_json(tags)
    d[:sources] = []
    HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(d))
end
```

- [ ] **Step 6: Add `GET /api/samples/:id`**

In `packages/HimalayaUI/src/routes_samples.jl`, inside `register_samples_routes!()`, add before the closing `end`:

```julia
@get "/api/samples/{id}" function(req::HTTP.Request, id::Int)
    db   = current_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM samples WHERE id = ?", [id]))
    isempty(rows) && return HTTP.Response(404,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "sample not found")))
    sm   = rows[1]
    tags = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, key, value, source FROM sample_tags
         WHERE sample_id = ? ORDER BY id", [id]))
    d        = row_to_json(sm)
    d[:tags] = rows_to_json(tags)
    HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(d))
end
```

- [ ] **Step 7: Add `GET /api/indices/:id`**

In `packages/HimalayaUI/src/routes_analysis.jl`, inside `register_analysis_routes!()`, add before the closing `end`:

```julia
@get "/api/indices/{id}" function(req::HTTP.Request, id::Int)
    db   = current_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM indices WHERE id = ?", [id]))
    isempty(rows) && return HTTP.Response(404,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "index not found")))
    ix        = rows[1]
    peak_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT ip.peak_id, ip.ratio_position, ip.residual, p.q AS q_observed
         FROM index_peaks ip JOIN peaks p ON p.id = ip.peak_id
         WHERE ip.index_id = ? ORDER BY ip.ratio_position", [id]))
    predicted = predicted_q_for_phase(String(ix.phase), Float64(ix.basis))
    d               = row_to_json(ix)
    d[:peaks]       = rows_to_json(peak_rows)
    d[:predicted_q] = predicted
    HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(d))
end
```

- [ ] **Step 8: Run tests — expect pass**

```bash
julia --project=. -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: all PASS.

- [ ] **Step 9: Commit**

```bash
git add packages/HimalayaUI/src/routes_peaks.jl \
        packages/HimalayaUI/src/routes_exposures.jl \
        packages/HimalayaUI/src/routes_samples.jl \
        packages/HimalayaUI/src/routes_analysis.jl \
        packages/HimalayaUI/test/test_routes_mentions.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "feat(api): add single-entity GET routes for mention resolution"
```

---

## Task 3: Frontend — api.ts fetchers + query hooks

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/api.ts`
- Modify: `packages/HimalayaUI/frontend/src/queries.ts`

- [ ] **Step 1: Add fetchers to api.ts**

Append to `packages/HimalayaUI/frontend/src/api.ts`:

```typescript
// Single-entity fetchers for mention resolution
export const getPeak     = (id: number) => request<Peak>("GET", `/api/peaks/${id}`);
export const getIndex    = (id: number) => request<IndexEntry>("GET", `/api/indices/${id}`);
export const getExposure = (id: number) => request<Exposure>("GET", `/api/exposures/${id}`);
export const getSample   = (id: number) => request<Sample>("GET", `/api/samples/${id}`);
```

- [ ] **Step 2: Add query keys + hooks to queries.ts**

Add to the `queryKeys` object in `packages/HimalayaUI/frontend/src/queries.ts`:

```typescript
peak:     (id: number) => ["peak", id] as const,
index:    (id: number) => ["index", id] as const,
exposure: (id: number) => ["exposure", id] as const,
sample:   (id: number) => ["sample", id] as const,
```

Append the four hooks to `queries.ts`:

```typescript
export function usePeak(id: number | undefined) {
  return useQuery({
    queryKey: id !== undefined ? queryKeys.peak(id) : (["peak", "none"] as const),
    queryFn: () => api.getPeak(id as number),
    enabled: id !== undefined,
    retry: false,
  });
}

export function useIndex(id: number | undefined) {
  return useQuery({
    queryKey: id !== undefined ? queryKeys.index(id) : (["index", "none"] as const),
    queryFn: () => api.getIndex(id as number),
    enabled: id !== undefined,
    retry: false,
  });
}

export function useExposure(id: number | undefined) {
  return useQuery({
    queryKey: id !== undefined ? queryKeys.exposure(id) : (["exposure", "none"] as const),
    queryFn: () => api.getExposure(id as number),
    enabled: id !== undefined,
    retry: false,
  });
}

export function useSampleById(id: number | undefined) {
  return useQuery({
    queryKey: id !== undefined ? queryKeys.sample(id) : (["sample", "none"] as const),
    queryFn: () => api.getSample(id as number),
    enabled: id !== undefined,
    retry: false,
  });
}
```

Note: the hook is named `useSampleById` to avoid shadowing the existing `useSamples` (which takes an `experimentId`).

- [ ] **Step 3: Run tests (no failures expected yet)**

```bash
cd packages/HimalayaUI/frontend
npm run build 2>&1 | tail -20
```

Expected: build passes with no type errors.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/api.ts packages/HimalayaUI/frontend/src/queries.ts
git commit -m "feat(ui): add single-entity fetchers and query hooks for mention resolution"
```

---

## Task 4: `renderMentions` — pure token parser

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/renderMentions.tsx`
- Create: `packages/HimalayaUI/frontend/test/renderMentions.test.tsx`

`renderMentions` is a pure function: it splits a message body string on `[[type:id]]` tokens and returns an array of either plain strings or structured token objects. It does NOT resolve or render — it just parses. The caller decides what to do with each token.

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/renderMentions.test.tsx`:

```typescript
import { describe, it, expect } from "vitest";
import { parseMentions, type MentionToken } from "../src/lib/renderMentions";

describe("parseMentions", () => {
  it("returns single string segment for plain text", () => {
    const result = parseMentions("hello world");
    expect(result).toEqual([{ kind: "text", text: "hello world" }]);
  });

  it("returns a mention token for a single [[type:id]]", () => {
    const result = parseMentions("[[peak:42]]");
    expect(result).toEqual([{ kind: "mention", type: "peak", id: 42 }]);
  });

  it("splits mixed text and tokens correctly", () => {
    const result = parseMentions("shoulder at [[peak:42]] matches [[index:17]]");
    expect(result).toEqual([
      { kind: "text", text: "shoulder at " },
      { kind: "mention", type: "peak", id: 42 },
      { kind: "text", text: " matches " },
      { kind: "mention", type: "index", id: 17 },
    ]);
  });

  it("handles token at start and end", () => {
    const result = parseMentions("[[sample:8]] and [[exposure:33]]");
    expect(result).toEqual([
      { kind: "mention", type: "sample", id: 8 },
      { kind: "text", text: " and " },
      { kind: "mention", type: "exposure", id: 33 },
    ]);
  });

  it("ignores empty text segments between adjacent tokens", () => {
    const result = parseMentions("[[peak:1]][[peak:2]]");
    expect(result).toEqual([
      { kind: "mention", type: "peak", id: 1 },
      { kind: "mention", type: "peak", id: 2 },
    ]);
  });

  it("returns plain text for malformed tokens", () => {
    const result = parseMentions("[[notatype:42]] and [[peak:abc]]");
    // Both are invalid — returned as plain text segments
    expect(result.every((s) => s.kind === "text")).toBe(true);
  });

  it("handles empty string", () => {
    expect(parseMentions("")).toEqual([]);
  });
});
```

- [ ] **Step 2: Run tests — expect failure**

```bash
node_modules/.bin/vitest run test/renderMentions.test.tsx
```

Expected: FAIL — module not found.

- [ ] **Step 3: Implement `parseMentions`**

Create `packages/HimalayaUI/frontend/src/lib/renderMentions.tsx`:

```typescript
const VALID_TYPES = ["peak", "index", "exposure", "sample", "experiment"] as const;
type MentionType = typeof VALID_TYPES[number];

export type TextSegment   = { kind: "text"; text: string };
export type MentionToken  = { kind: "mention"; type: MentionType; id: number };
export type BodySegment   = TextSegment | MentionToken;

const TOKEN_RE = /\[\[(\w+):(\d+)\]\]/g;

export function parseMentions(body: string): BodySegment[] {
  if (!body) return [];
  const segments: BodySegment[] = [];
  let last = 0;

  for (const match of body.matchAll(TOKEN_RE)) {
    const type = match[1];
    const id   = parseInt(match[2], 10);
    if (!(VALID_TYPES as readonly string[]).includes(type) || Number.isNaN(id)) {
      // Treat invalid token as literal text — advance past it
      const textChunk = body.slice(last, match.index! + match[0].length);
      if (textChunk) segments.push({ kind: "text", text: textChunk });
      last = match.index! + match[0].length;
      continue;
    }
    if (match.index! > last) {
      segments.push({ kind: "text", text: body.slice(last, match.index) });
    }
    segments.push({ kind: "mention", type: type as MentionType, id });
    last = match.index! + match[0].length;
  }

  if (last < body.length) {
    segments.push({ kind: "text", text: body.slice(last) });
  }
  return segments;
}
```

- [ ] **Step 4: Run tests — expect pass**

```bash
node_modules/.bin/vitest run test/renderMentions.test.tsx
```

Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/renderMentions.tsx \
        packages/HimalayaUI/frontend/test/renderMentions.test.tsx
git commit -m "feat(ui): add parseMentions — pure [[type:id]] token parser"
```

---

## Task 5: `useMentionResolution` hook

**Files:**
- Create: `packages/HimalayaUI/frontend/src/hooks/useMentionResolution.ts`
- Create: `packages/HimalayaUI/frontend/test/useMentionResolution.test.ts`

This hook accepts an array of `MentionToken` objects and returns a resolution map. It checks TanStack Query cache first; for unknown IDs it fires lazy fetches via the hooks from Task 3.

The hook's return type is `Map<string, ResolvedMention | "loading" | "dead">` where the key is `"${type}:${id}"`.

- [ ] **Step 1: Define `ResolvedMention` type and write failing tests**

Create `packages/HimalayaUI/frontend/test/useMentionResolution.test.ts`:

```typescript
import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { makeClient } from "./test-utils";
import { QueryClientProvider } from "@tanstack/react-query";
import { createElement } from "react";
import * as api from "../src/api";
import { useMentionResolution } from "../src/hooks/useMentionResolution";
import type { MentionToken } from "../src/lib/renderMentions";

const PEAK: api.Peak = {
  id: 42, exposure_id: 1, q: 1.223, intensity: 841.2,
  prominence: 4.2, sharpness: 0.3, source: "auto", excluded: false,
};

function wrapper(client = makeClient()) {
  return ({ children }: { children: React.ReactNode }) =>
    createElement(QueryClientProvider, { client }, children);
}

describe("useMentionResolution", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("returns 'loading' initially when entity is not in cache", async () => {
    vi.spyOn(api, "getPeak").mockResolvedValue(PEAK);
    const tokens: MentionToken[] = [{ kind: "mention", type: "peak", id: 42 }];
    const { result } = renderHook(() => useMentionResolution(tokens), { wrapper: wrapper() });
    expect(result.current.get("peak:42")).toBe("loading");
  });

  it("resolves to entity data after fetch", async () => {
    vi.spyOn(api, "getPeak").mockResolvedValue(PEAK);
    const tokens: MentionToken[] = [{ kind: "mention", type: "peak", id: 42 }];
    const { result } = renderHook(() => useMentionResolution(tokens), { wrapper: wrapper() });
    await waitFor(() => {
      const entry = result.current.get("peak:42");
      expect(entry).not.toBe("loading");
      expect(entry).not.toBe("dead");
    });
    const entry = result.current.get("peak:42");
    expect(entry).toMatchObject({ type: "peak", data: PEAK });
  });

  it("marks entity as 'dead' on 404", async () => {
    vi.spyOn(api, "getPeak").mockRejectedValue(
      new api.ApiError(404, "not found", null)
    );
    const tokens: MentionToken[] = [{ kind: "mention", type: "peak", id: 999 }];
    const { result } = renderHook(() => useMentionResolution(tokens), { wrapper: wrapper() });
    await waitFor(() => result.current.get("peak:999") === "dead");
    expect(result.current.get("peak:999")).toBe("dead");
  });

  it("returns empty map for empty token list", () => {
    const { result } = renderHook(() => useMentionResolution([]), { wrapper: wrapper() });
    expect(result.current.size).toBe(0);
  });
});
```

- [ ] **Step 2: Run tests — expect failure**

```bash
node_modules/.bin/vitest run test/useMentionResolution.test.ts
```

Expected: FAIL — module not found.

- [ ] **Step 3: Implement `useMentionResolution`**

Create `packages/HimalayaUI/frontend/src/hooks/useMentionResolution.ts`:

```typescript
import { useMemo } from "react";
import { useQueries } from "@tanstack/react-query";
import * as api from "../api";
import { queryKeys } from "../queries";
import type { MentionToken } from "../lib/renderMentions";

export type ResolvedMention =
  | { type: "peak";       data: api.Peak }
  | { type: "index";      data: api.IndexEntry }
  | { type: "exposure";   data: api.Exposure }
  | { type: "sample";     data: api.Sample }
  | { type: "experiment"; data: api.Experiment };

type ResolutionEntry = ResolvedMention | "loading" | "dead";

function queryForToken(token: MentionToken) {
  const { type, id } = token;
  switch (type) {
    case "peak":
      return { queryKey: queryKeys.peak(id), queryFn: () => api.getPeak(id), retry: false };
    case "index":
      return { queryKey: queryKeys.index(id), queryFn: () => api.getIndex(id), retry: false };
    case "exposure":
      return { queryKey: queryKeys.exposure(id), queryFn: () => api.getExposure(id), retry: false };
    case "sample":
      return { queryKey: queryKeys.sample(id), queryFn: () => api.getSample(id), retry: false };
    case "experiment":
      return { queryKey: queryKeys.experiment(id), queryFn: () => api.getExperiment(id), retry: false };
  }
}

export function useMentionResolution(tokens: MentionToken[]): Map<string, ResolutionEntry> {
  const queries = useQueries({ queries: tokens.map(queryForToken) });

  return useMemo(() => {
    const map = new Map<string, ResolutionEntry>();
    tokens.forEach((token, i) => {
      const key = `${token.type}:${token.id}`;
      const q   = queries[i];
      if (q.isError) {
        const err = q.error;
        map.set(key, err instanceof api.ApiError && err.status === 404 ? "dead" : "loading");
      } else if (q.isSuccess && q.data !== undefined) {
        map.set(key, { type: token.type, data: q.data } as ResolvedMention);
      } else {
        map.set(key, "loading");
      }
    });
    return map;
  }, [tokens, queries]);
}
```

- [ ] **Step 4: Run tests — expect pass**

```bash
node_modules/.bin/vitest run test/useMentionResolution.test.ts
```

Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/hooks/useMentionResolution.ts \
        packages/HimalayaUI/frontend/test/useMentionResolution.test.ts
git commit -m "feat(ui): add useMentionResolution — hybrid cache-first mention resolution"
```

---

## Task 6: `MentionChip` component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/MentionChip.tsx`
- Create: `packages/HimalayaUI/frontend/test/MentionChip.test.tsx`

Renders a single resolved mention as an inline chip with hover tooltip and Zustand wiring for peak/index highlights.

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/MentionChip.test.tsx`:

```typescript
import { describe, it, expect, beforeEach } from "vitest";
import { screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { MentionChip } from "../src/components/MentionChip";
import { useAppState } from "../src/state";
import type { ResolvedMention } from "../src/hooks/useMentionResolution";
import * as api from "../src/api";

const PEAK: api.Peak = {
  id: 42, exposure_id: 1, q: 1.223, intensity: 841, prominence: 4.2,
  sharpness: 0.3, source: "auto", excluded: false,
};

const INDEX: api.IndexEntry = {
  id: 17, exposure_id: 1, phase: "Pn3m", basis: 1.0, score: 0.91,
  r_squared: 0.998, lattice_d: 5.14, status: "candidate",
  peaks: [], predicted_q: [1.223, 1.414, 1.732],
};

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({ hoveredPeakId: undefined, hoveredIndexId: undefined });
});

describe("<MentionChip> — resolved", () => {
  it("renders peak chip with q value", () => {
    const resolved: ResolvedMention = { type: "peak", data: PEAK };
    renderWithProviders(<MentionChip resolved={resolved} />);
    expect(screen.getByText(/q = 1\.223/)).toBeInTheDocument();
  });

  it("renders index chip with phase and score", () => {
    const resolved: ResolvedMention = { type: "index", data: INDEX };
    renderWithProviders(<MentionChip resolved={resolved} />);
    expect(screen.getByText(/Pn3m/)).toBeInTheDocument();
    expect(screen.getByText(/0\.91/)).toBeInTheDocument();
  });

  it("sets hoveredPeakId on mouse enter for peak chip", async () => {
    const resolved: ResolvedMention = { type: "peak", data: PEAK };
    const user = userEvent.setup();
    renderWithProviders(<MentionChip resolved={resolved} />);
    await user.hover(screen.getByText(/q = 1\.223/));
    expect(useAppState.getState().hoveredPeakId).toBe(42);
  });

  it("clears hoveredPeakId on mouse leave for peak chip", async () => {
    const resolved: ResolvedMention = { type: "peak", data: PEAK };
    const user = userEvent.setup();
    renderWithProviders(<MentionChip resolved={resolved} />);
    await user.hover(screen.getByText(/q = 1\.223/));
    await user.unhover(screen.getByText(/q = 1\.223/));
    expect(useAppState.getState().hoveredPeakId).toBeUndefined();
  });

  it("sets hoveredIndexId on mouse enter for index chip", async () => {
    const resolved: ResolvedMention = { type: "index", data: INDEX };
    const user = userEvent.setup();
    renderWithProviders(<MentionChip resolved={resolved} />);
    await user.hover(screen.getByText(/Pn3m/));
    expect(useAppState.getState().hoveredIndexId).toBe(17);
  });
});

describe("<MentionChip> — loading and dead", () => {
  it("renders nothing visible for 'loading' state", () => {
    renderWithProviders(<MentionChip resolved="loading" originalText="Pn3m · 0.91" />);
    // Loading state shows the original text dimmed
    expect(screen.getByText("Pn3m · 0.91")).toBeInTheDocument();
  });

  it("renders grayed chip for 'dead' state", () => {
    renderWithProviders(<MentionChip resolved="dead" originalText="Pn3m · 0.91" />);
    const chip = screen.getByText("Pn3m · 0.91");
    expect(chip).toBeInTheDocument();
    expect(chip.closest("[data-mention-state='dead']")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run tests — expect failure**

```bash
node_modules/.bin/vitest run test/MentionChip.test.tsx
```

Expected: FAIL — module not found.

- [ ] **Step 3: Implement `MentionChip`**

Create `packages/HimalayaUI/frontend/src/components/MentionChip.tsx`:

```typescript
import { useState, useCallback } from "react";
import { useAppState } from "../state";
import { KNOWN_PHASES } from "../phases";
import type { ResolvedMention } from "../hooks/useMentionResolution";

const CUBIC_PHASES: ReadonlySet<string> = new Set(["Pn3m", "Im3m", "Ia3d", "Fm3m", "Fd3m"]);

interface ChipProps {
  resolved: ResolvedMention | "loading" | "dead";
  originalText?: string;
}

const CHIP_STYLES = {
  peak:     "text-[#7cb8e8] bg-[#15222e] border-[#4a7aaa]",
  index:    "text-[#b5a0d8] bg-[#1e1828] border-[#7a60a8]",
  exposure: "text-[#88c0a8] bg-[#162018] border-[#508070]",
  sample:   "text-[#c0b878] bg-[#201c10] border-[#887840]",
  experiment: "text-[#c0b878] bg-[#201c10] border-[#887840]",
  dead:     "text-[#484848] bg-[#181818] border-[#333333]",
} as const;

const CHIP_HOVER_STYLES = {
  peak:     "hover:text-[#a8d4f8] hover:bg-[#1c3045] hover:border-[#7cb8e8]",
  index:    "hover:text-[#cdb8f0] hover:bg-[#2a2040] hover:border-[#b5a0d8]",
  exposure: "hover:text-[#a8e0c0] hover:bg-[#1c3028] hover:border-[#88c0a8]",
  sample:   "hover:text-[#dcd090] hover:bg-[#282215] hover:border-[#c0b878]",
  experiment: "hover:text-[#dcd090] hover:bg-[#282215] hover:border-[#c0b878]",
  dead:     "hover:text-[#606060] hover:bg-[#1e1e1e] hover:border-[#444444]",
} as const;

function chipLabel(resolved: ResolvedMention): string {
  switch (resolved.type) {
    case "peak":       return `q = ${resolved.data.q.toFixed(3)}`;
    case "index":      return `${resolved.data.phase} · ${(resolved.data.score ?? 0).toFixed(2)}`;
    case "exposure":   return resolved.data.filename ?? `exposure ${resolved.data.id}`;
    case "sample":     return resolved.data.name ?? resolved.data.label ?? `sample ${resolved.data.id}`;
    case "experiment": return resolved.data.name ?? `experiment ${resolved.data.id}`;
  }
}

function TooltipContent({ resolved }: { resolved: ResolvedMention }): JSX.Element | null {
  switch (resolved.type) {
    case "peak":
      return (
        <span>
          source <code>{resolved.data.source}</code>
          {resolved.data.prominence != null && <> · prominence <code>{resolved.data.prominence.toFixed(1)}</code></>}
        </span>
      );
    case "index": {
      const q1 = resolved.data.predicted_q[0];
      const isCubic = CUBIC_PHASES.has(resolved.data.phase);
      return (
        <span>
          {q1 != null && <>q₁ <code>{q1.toFixed(3)}</code> · </>}
          {resolved.data.lattice_d != null && <>d <code>{resolved.data.lattice_d.toFixed(2)} nm</code> · </>}
          R² <code>{(resolved.data.r_squared ?? 0).toFixed(3)}</code>
          {isCubic && (resolved.data as unknown as { ngc?: number }).ngc != null && (
            <> · ngc <code>{(resolved.data as unknown as { ngc: number }).ngc.toFixed(1)}</code></>
          )}
        </span>
      );
    }
    case "exposure":
      return resolved.data.status != null
        ? <span>status <code>{resolved.data.status}</code></span>
        : null;
    case "sample":
    case "experiment":
      return null;
  }
}

export function MentionChip({ resolved, originalText }: ChipProps): JSX.Element {
  const [isHovered, setIsHovered] = useState(false);
  const setHoveredPeak  = useAppState((s) => s.setHoveredPeak);
  const setHoveredIndex = useAppState((s) => s.setHoveredIndex);

  const handleMouseEnter = useCallback(() => {
    setIsHovered(true);
    if (resolved === "loading" || resolved === "dead") return;
    if (resolved.type === "peak")  setHoveredPeak(resolved.data.id);
    if (resolved.type === "index") setHoveredIndex(resolved.data.id);
  }, [resolved, setHoveredPeak, setHoveredIndex]);

  const handleMouseLeave = useCallback(() => {
    setIsHovered(false);
    if (resolved === "loading" || resolved === "dead") return;
    if (resolved.type === "peak")  setHoveredPeak(undefined);
    if (resolved.type === "index") setHoveredIndex(undefined);
  }, [resolved, setHoveredPeak, setHoveredIndex]);

  const state = resolved === "loading" ? "loading" : resolved === "dead" ? "dead" : resolved.type;
  const styleKey = state === "loading" ? "dead" : state;
  const baseStyle = CHIP_STYLES[styleKey as keyof typeof CHIP_STYLES] ?? CHIP_STYLES.dead;
  const hoverStyle = CHIP_HOVER_STYLES[styleKey as keyof typeof CHIP_HOVER_STYLES] ?? CHIP_HOVER_STYLES.dead;

  const label = resolved === "loading" || resolved === "dead"
    ? (originalText ?? "…")
    : chipLabel(resolved);

  const tooltip = isHovered && resolved !== "loading"
    ? (resolved === "dead"
        ? <span className="text-fg-dim">no longer exists</span>
        : <TooltipContent resolved={resolved} />)
    : null;

  return (
    <span
      data-mention-state={state}
      className={`relative inline border-b pb-px px-1 rounded-sm cursor-pointer whitespace-nowrap
                  text-sm transition-colors ${baseStyle} ${hoverStyle}`}
      onMouseEnter={handleMouseEnter}
      onMouseLeave={handleMouseLeave}
    >
      {label}
      {tooltip && (
        <span className="absolute bottom-[calc(100%+6px)] left-1/2 -translate-x-1/2 z-10
                         bg-surface border border-border rounded-md px-2 py-1
                         text-xs text-fg-muted whitespace-nowrap shadow-lg pointer-events-none">
          {tooltip}
          <span className="absolute top-full left-1/2 -translate-x-1/2
                           border-4 border-transparent border-t-border" />
        </span>
      )}
    </span>
  );
}
```

- [ ] **Step 4: Run tests — expect pass**

```bash
node_modules/.bin/vitest run test/MentionChip.test.tsx
```

Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/MentionChip.tsx \
        packages/HimalayaUI/frontend/test/MentionChip.test.tsx
git commit -m "feat(ui): add MentionChip — inline mention chip with hover tooltip and Zustand wiring"
```

---

## Task 7: Wire `parseMentions` + `useMentionResolution` + `MentionChip` into `MessageRow`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ChatCard.tsx`
- Modify: `packages/HimalayaUI/frontend/test/ChatCard.test.tsx`

`MessageRow` currently renders `msg.body` as plain text. We replace the `<p>` content with a rendered mention list. The resolution hook fires at the `MessageRow` level so each message resolves its own tokens independently.

- [ ] **Step 1: Add mention rendering tests to ChatCard.test.tsx**

Add to `packages/HimalayaUI/frontend/test/ChatCard.test.tsx` inside the `describe("<ChatCard>")` block:

```typescript
it("renders a peak mention chip inline", async () => {
  const PEAK: api.Peak = {
    id: 42, exposure_id: 1, q: 1.223, intensity: 841, prominence: 4.2,
    sharpness: 0.3, source: "auto", excluded: false,
  };
  vi.spyOn(api, "listSampleMessages").mockResolvedValue([
    { id: 1, sample_id: 3, author_id: 1, author: "alice",
      body: "see [[peak:42]]", created_at: "2026-04-24 10:00:00" },
  ]);
  vi.spyOn(api, "getPeak").mockResolvedValue(PEAK);
  renderWithProviders(<ChatCard />);
  expect(await screen.findByText(/q = 1\.223/)).toBeInTheDocument();
});

it("renders dead chip when mention entity returns 404", async () => {
  vi.spyOn(api, "listSampleMessages").mockResolvedValue([
    { id: 2, sample_id: 3, author_id: 1, author: "alice",
      body: "old ref [[index:999]]", created_at: "2026-04-24 10:00:00" },
  ]);
  vi.spyOn(api, "getIndex").mockRejectedValue(new api.ApiError(404, "not found", null));
  renderWithProviders(<ChatCard />);
  await screen.findByText("old ref ");
  const chip = document.querySelector("[data-mention-state='dead']");
  expect(chip).not.toBeNull();
});
```

- [ ] **Step 2: Run tests — expect failure**

```bash
node_modules/.bin/vitest run test/ChatCard.test.tsx
```

Expected: new tests FAIL — body still renders as plain text.

- [ ] **Step 3: Update `MessageRow` in ChatCard.tsx**

Replace the `MessageRow` function and add imports in `packages/HimalayaUI/frontend/src/components/ChatCard.tsx`:

```typescript
import { useEffect, useRef, useState, useMemo } from "react";
import { parseMentions } from "../lib/renderMentions";
import { useMentionResolution } from "../hooks/useMentionResolution";
import { MentionChip } from "./MentionChip";
```

Replace `MessageRow`:

```typescript
function MessageRow({ msg }: { msg: SampleMessage }): JSX.Element {
  const authorLabel   = msg.author ?? "deleted user";
  const authorDeleted = msg.author == null;
  const segments      = useMemo(() => parseMentions(msg.body), [msg.body]);
  const mentions      = useMemo(
    () => segments.filter((s): s is import("../lib/renderMentions").MentionToken => s.kind === "mention"),
    [segments],
  );
  const resolutionMap = useMentionResolution(mentions);

  return (
    <div className="flex flex-col gap-0.5 min-w-0" data-testid={`chat-message-${msg.id}`}>
      <div className="flex items-baseline gap-2">
        <span className={authorDeleted ? "text-meta text-fg-dim italic" : "text-meta"}>
          {authorLabel}
        </span>
        <span className="text-fg-dim text-xs">{formatTime(msg.created_at)}</span>
      </div>
      <p className="text-base font-sans text-fg-muted leading-snug break-words whitespace-pre-wrap">
        {segments.map((seg, i) => {
          if (seg.kind === "text") return <span key={i}>{seg.text}</span>;
          const key    = `${seg.type}:${seg.id}`;
          const token  = `[[${seg.type}:${seg.id}]]`;
          const entry  = resolutionMap.get(key) ?? "loading";
          const originalText = token; // fallback display for loading/dead
          return (
            <MentionChip key={i} resolved={entry} originalText={originalText} />
          );
        })}
      </p>
    </div>
  );
}
```

- [ ] **Step 4: Run tests — expect pass**

```bash
node_modules/.bin/vitest run test/ChatCard.test.tsx
```

Expected: all PASS including the two new mention tests.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ChatCard.tsx \
        packages/HimalayaUI/frontend/test/ChatCard.test.tsx
git commit -m "feat(ui): render mention chips inline in MessageRow"
```

---

## Task 8: `MentionPicker` — autocomplete dropdown

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/MentionPicker.tsx`
- Create: `packages/HimalayaUI/frontend/test/MentionPicker.test.tsx`

Fuzzy search across indices, peaks, exposures, and samples loaded from TanStack Query. Grouped and ranked by proximity to active context.

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/MentionPicker.test.tsx`:

```typescript
import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { MentionPicker } from "../src/components/MentionPicker";
import { useAppState } from "../src/state";
import * as api from "../src/api";

const INDICES: api.IndexEntry[] = [
  { id: 17, exposure_id: 1, phase: "Pn3m", basis: 1.0, score: 0.91, r_squared: 0.998,
    lattice_d: 5.14, status: "candidate", peaks: [], predicted_q: [1.223, 1.414] },
];
const PEAKS: api.Peak[] = [
  { id: 42, exposure_id: 1, q: 1.223, intensity: 841, prominence: 4.2,
    sharpness: 0.3, source: "auto", excluded: false },
];
const EXPOSURES: api.Exposure[] = [
  { id: 1, sample_id: 3, filename: "JC001-120.dat", kind: "file",
    selected: true, status: "accepted", image_path: null, image_version: "", tags: [], sources: [] },
];

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({ activeSampleId: 3, activeExposureId: 1, activeExperimentId: 1 });
  vi.spyOn(api, "listIndices").mockResolvedValue(INDICES);
  vi.spyOn(api, "listPeaks").mockResolvedValue(PEAKS);
  vi.spyOn(api, "listExposures").mockResolvedValue(EXPOSURES);
  vi.spyOn(api, "listSamples").mockResolvedValue([]);
});

describe("<MentionPicker>", () => {
  it("shows indices section when query matches phase", async () => {
    renderWithProviders(
      <MentionPicker query="Pn3" onSelect={vi.fn()} onDismiss={vi.fn()} />
    );
    expect(await screen.findByText(/Pn3m/)).toBeInTheDocument();
  });

  it("shows peak section when query matches q value", async () => {
    renderWithProviders(
      <MentionPicker query="1.22" onSelect={vi.fn()} onDismiss={vi.fn()} />
    );
    expect(await screen.findByText(/1\.223/)).toBeInTheDocument();
  });

  it("calls onSelect with [[type:id]] token when row is clicked", async () => {
    const onSelect = vi.fn();
    const user = userEvent.setup();
    renderWithProviders(
      <MentionPicker query="Pn3" onSelect={onSelect} onDismiss={vi.fn()} />
    );
    await user.click(await screen.findByText(/Pn3m/));
    expect(onSelect).toHaveBeenCalledWith("[[index:17]]");
  });

  it("calls onDismiss on Escape key", async () => {
    const onDismiss = vi.fn();
    const user = userEvent.setup();
    renderWithProviders(
      <MentionPicker query="Pn3" onSelect={vi.fn()} onDismiss={onDismiss} />
    );
    await user.keyboard("{Escape}");
    expect(onDismiss).toHaveBeenCalled();
  });

  it("shows empty state when no results match", async () => {
    renderWithProviders(
      <MentionPicker query="zzznothing" onSelect={vi.fn()} onDismiss={vi.fn()} />
    );
    expect(await screen.findByText(/no results/i)).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run tests — expect failure**

```bash
node_modules/.bin/vitest run test/MentionPicker.test.tsx
```

Expected: FAIL — module not found.

- [ ] **Step 3: Implement `MentionPicker`**

Create `packages/HimalayaUI/frontend/src/components/MentionPicker.tsx`:

```typescript
import { useEffect, useRef, useState, useMemo } from "react";
import { useAppState } from "../state";
import { useIndices, usePeaks, useExposures, useSamples } from "../queries";
import { phaseColor } from "../phases";
import type { IndexEntry, Peak, Exposure, Sample } from "../api";

type PickerRow =
  | { kind: "index";    item: IndexEntry }
  | { kind: "peak";     item: Peak }
  | { kind: "exposure"; item: Exposure }
  | { kind: "sample";   item: Sample };

interface MentionPickerProps {
  query: string;
  onSelect: (token: string) => void;
  onDismiss: () => void;
}

function rowToken(row: PickerRow): string {
  switch (row.kind) {
    case "index":    return `[[index:${row.item.id}]]`;
    case "peak":     return `[[peak:${row.item.id}]]`;
    case "exposure": return `[[exposure:${row.item.id}]]`;
    case "sample":   return `[[sample:${row.item.id}]]`;
  }
}

function matchesQuery(text: string, q: string): boolean {
  return text.toLowerCase().includes(q.toLowerCase());
}

// Parse optional type prefix: "@peak:1.22" → { type: "peak", rest: "1.22" }
function parseQuery(raw: string): { type: string | null; rest: string } {
  const m = raw.match(/^(\w+):(.*)$/);
  if (m) return { type: m[1].toLowerCase(), rest: m[2] };
  return { type: null, rest: raw };
}

export function MentionPicker({ query, onSelect, onDismiss }: MentionPickerProps): JSX.Element {
  const activeExposureId  = useAppState((s) => s.activeExposureId);
  const activeSampleId    = useAppState((s) => s.activeSampleId);
  const activeExperimentId = useAppState((s) => s.activeExperimentId);

  const indicesQ  = useIndices(activeExposureId);
  const peaksQ    = usePeaks(activeExposureId);
  const exposuresQ = useExposures(activeSampleId);
  const samplesQ  = useSamples(activeExperimentId ?? 0);

  const [activeIdx, setActiveIdx] = useState(0);
  const listRef = useRef<HTMLDivElement>(null);

  const { type: typeFilter, rest: searchText } = useMemo(() => parseQuery(query), [query]);

  const rows = useMemo((): PickerRow[] => {
    const q = searchText;
    const all: PickerRow[] = [];

    const wantAll   = typeFilter === null;
    const wantIndex = wantAll || typeFilter === "index";
    const wantPeak  = wantAll || typeFilter === "peak";
    const wantExp   = wantAll || typeFilter === "exposure";
    const wantSamp  = wantAll || typeFilter === "sample";

    if (wantIndex) {
      (indicesQ.data ?? [])
        .filter((ix) => !q || matchesQuery(ix.phase, q) ||
          matchesQuery((ix.score ?? 0).toFixed(2), q) ||
          ix.predicted_q.some((pq) => matchesQuery(pq.toFixed(3), q)))
        .forEach((ix) => all.push({ kind: "index", item: ix }));
    }
    if (wantPeak) {
      (peaksQ.data ?? [])
        .filter((pk) => !q || matchesQuery(pk.q.toFixed(3), q))
        .forEach((pk) => all.push({ kind: "peak", item: pk }));
    }
    if (wantExp) {
      (exposuresQ.data ?? [])
        .filter((ex) => !q || matchesQuery(ex.filename ?? "", q))
        .forEach((ex) => all.push({ kind: "exposure", item: ex }));
    }
    if (wantSamp) {
      (samplesQ.data ?? [])
        .filter((sm) => !q || matchesQuery(sm.name ?? sm.label ?? "", q))
        .forEach((sm) => all.push({ kind: "sample", item: sm }));
    }
    return all;
  }, [query, typeFilter, searchText, indicesQ.data, peaksQ.data, exposuresQ.data, samplesQ.data]);

  useEffect(() => { setActiveIdx(0); }, [rows]);

  useEffect(() => {
    function onKey(e: KeyboardEvent) {
      if (e.key === "Escape") { onDismiss(); return; }
      if (e.key === "ArrowDown") { e.preventDefault(); setActiveIdx((i) => Math.min(i + 1, rows.length - 1)); }
      if (e.key === "ArrowUp")   { e.preventDefault(); setActiveIdx((i) => Math.max(i - 1, 0)); }
      if (e.key === "Enter" && rows[activeIdx]) { e.preventDefault(); onSelect(rowToken(rows[activeIdx])); }
    }
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [rows, activeIdx, onSelect, onDismiss]);

  function rowLabel(row: PickerRow): string {
    switch (row.kind) {
      case "index":    return `${row.item.phase} · ${(row.item.score ?? 0).toFixed(2)}`;
      case "peak":     return `q = ${row.item.item?.q.toFixed(3) ?? (row.item as Peak).q.toFixed(3)}`;
      case "exposure": return row.item.filename ?? `exposure ${row.item.id}`;
      case "sample":   return row.item.name ?? row.item.label ?? `sample ${row.item.id}`;
    }
  }

  function rowMeta(row: PickerRow): string | null {
    switch (row.kind) {
      case "index":    return `score ${(row.item.score ?? 0).toFixed(2)}`;
      case "peak":     return `${row.item.source} · prom ${(row.item.prominence ?? 0).toFixed(1)}`;
      case "exposure": return row.item.status ?? null;
      case "sample":   return null;
    }
  }

  return (
    <div
      ref={listRef}
      role="listbox"
      className="absolute bottom-full left-0 right-0 mb-1 z-20
                 bg-surface border border-border rounded-lg overflow-hidden shadow-xl"
    >
      <div className="px-3 py-1.5 border-b border-border bg-bg
                      text-xs text-fg-dim flex justify-between">
        <span>@{query || "…"}</span>
        <span className="text-fg-dim opacity-50">↑↓ navigate · Enter select · Esc dismiss</span>
      </div>
      {rows.length === 0 ? (
        <div className="px-3 py-2 text-xs text-fg-dim">No results</div>
      ) : (
        rows.map((row, i) => {
          const color = row.kind === "index" ? phaseColor(row.item.phase) : undefined;
          return (
            <div
              key={`${row.kind}:${row.item.id}`}
              role="option"
              aria-selected={i === activeIdx}
              onClick={() => onSelect(rowToken(row))}
              className={`px-3 py-1.5 cursor-pointer flex justify-between items-center text-sm
                          ${i === activeIdx ? "bg-surface-hover" : "hover:bg-surface-hover"}`}
            >
              <span style={color ? { color } : undefined}>
                {rowLabel(row)}
              </span>
              {rowMeta(row) && (
                <span className="text-xs text-fg-dim ml-2">{rowMeta(row)}</span>
              )}
            </div>
          );
        })
      )}
    </div>
  );
}
```

Fix the `rowLabel` for peak (copy-paste error):

```typescript
case "peak": return `q = ${(row.item as Peak).q.toFixed(3)}`;
```

- [ ] **Step 4: Run tests — expect pass**

```bash
node_modules/.bin/vitest run test/MentionPicker.test.tsx
```

Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/MentionPicker.tsx \
        packages/HimalayaUI/frontend/test/MentionPicker.test.tsx
git commit -m "feat(ui): add MentionPicker autocomplete dropdown"
```

---

## Task 9: `MentionCompose` — compose wrapper with @ trigger and pretext auto-resize

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/MentionCompose.tsx`
- Create: `packages/HimalayaUI/frontend/test/MentionCompose.test.tsx`
- Modify: `packages/HimalayaUI/frontend/package.json`

- [ ] **Step 1: Install pretext**

```bash
cd packages/HimalayaUI/frontend
npm install @chenglou/pretext
```

Verify install:

```bash
node -e "require('@chenglou/pretext'); console.log('ok')"
```

Expected output: `ok`

- [ ] **Step 2: Write failing tests**

Create `packages/HimalayaUI/frontend/test/MentionCompose.test.tsx`:

```typescript
import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { MentionCompose } from "../src/components/MentionCompose";
import { useAppState } from "../src/state";
import * as api from "../src/api";

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({ activeSampleId: 3, activeExposureId: 1, activeExperimentId: 1 });
  vi.spyOn(api, "listIndices").mockResolvedValue([]);
  vi.spyOn(api, "listPeaks").mockResolvedValue([]);
  vi.spyOn(api, "listExposures").mockResolvedValue([]);
  vi.spyOn(api, "listSamples").mockResolvedValue([]);
});

describe("<MentionCompose>", () => {
  it("shows picker when @ is typed", async () => {
    const user = userEvent.setup();
    renderWithProviders(<MentionCompose disabled={false} onSubmit={vi.fn()} />);
    const ta = screen.getByRole("textbox");
    await user.click(ta);
    await user.keyboard("@");
    expect(await screen.findByRole("listbox")).toBeInTheDocument();
  });

  it("dismisses picker on Escape", async () => {
    const user = userEvent.setup();
    renderWithProviders(<MentionCompose disabled={false} onSubmit={vi.fn()} />);
    const ta = screen.getByRole("textbox");
    await user.click(ta);
    await user.keyboard("@foo");
    await screen.findByRole("listbox");
    await user.keyboard("{Escape}");
    await waitFor(() => {
      expect(screen.queryByRole("listbox")).toBeNull();
    });
  });

  it("inserts [[type:id]] token and closes picker on selection", async () => {
    vi.spyOn(api, "listIndices").mockResolvedValue([
      { id: 17, exposure_id: 1, phase: "Pn3m", basis: 1.0, score: 0.91,
        r_squared: 0.998, lattice_d: 5.14, status: "candidate", peaks: [],
        predicted_q: [1.223] },
    ]);
    const user = userEvent.setup();
    renderWithProviders(<MentionCompose disabled={false} onSubmit={vi.fn()} />);
    const ta = screen.getByRole("textbox");
    await user.click(ta);
    await user.keyboard("see @Pn3");
    const row = await screen.findByText(/Pn3m/);
    await user.click(row);
    expect((ta as HTMLTextAreaElement).value).toBe("see [[index:17]]");
    expect(screen.queryByRole("listbox")).toBeNull();
  });

  it("submits the raw token body on Enter", async () => {
    const onSubmit = vi.fn();
    const user = userEvent.setup();
    renderWithProviders(<MentionCompose disabled={false} onSubmit={onSubmit} />);
    const ta = screen.getByRole("textbox");
    await user.click(ta);
    await user.keyboard("hello");
    await user.keyboard("{Enter}");
    expect(onSubmit).toHaveBeenCalledWith("hello");
  });
});
```

- [ ] **Step 3: Run tests — expect failure**

```bash
node_modules/.bin/vitest run test/MentionCompose.test.tsx
```

Expected: FAIL — module not found.

- [ ] **Step 4: Implement `MentionCompose`**

Create `packages/HimalayaUI/frontend/src/components/MentionCompose.tsx`:

```typescript
import { useState, useRef, useCallback } from "react";
import { prepare, layout } from "@chenglou/pretext";
import { MentionPicker } from "./MentionPicker";

interface MentionComposeProps {
  disabled: boolean;
  onSubmit: (body: string) => void;
}

function computeTextareaHeight(text: string, width: number, font: string): number {
  if (!text || width <= 0) return 40;
  const prepared = prepare(text || " ", font, { whiteSpace: "pre-wrap" });
  const { height } = layout(prepared, width, 20);
  return Math.max(40, Math.min(height + 20, 200));
}

export function MentionCompose({ disabled, onSubmit }: MentionComposeProps): JSX.Element {
  const [text, setText] = useState("");
  const [pickerQuery, setPickerQuery] = useState<string | null>(null);
  const [atStart, setAtStart] = useState(0);
  const ref = useRef<HTMLTextAreaElement>(null);
  const wrapRef = useRef<HTMLDivElement>(null);

  const updateText = useCallback((val: string, cursor: number) => {
    setText(val);

    // Find the closest @ before cursor with no space between
    const before = val.slice(0, cursor);
    const atIdx  = before.lastIndexOf("@");
    if (atIdx !== -1 && !before.slice(atIdx + 1).includes(" ")) {
      setAtStart(atIdx);
      setPickerQuery(before.slice(atIdx + 1));
    } else {
      setPickerQuery(null);
    }

    // Auto-resize via pretext
    if (ref.current && wrapRef.current) {
      const w = wrapRef.current.clientWidth - 20;
      const font = getComputedStyle(ref.current).font;
      ref.current.style.height = `${computeTextareaHeight(val, w, font)}px`;
    }
  }, []);

  const handleChange = (e: React.ChangeEvent<HTMLTextAreaElement>): void => {
    updateText(e.target.value, e.target.selectionStart ?? e.target.value.length);
  };

  const handleSelect = useCallback((token: string) => {
    const before = text.slice(0, atStart);
    const after  = text.slice(ref.current?.selectionStart ?? text.length);
    const next   = before + token + after;
    setText(next);
    setPickerQuery(null);
    ref.current?.focus();
  }, [text, atStart]);

  const trySubmit = useCallback(() => {
    const trimmed = text.trim();
    if (!trimmed || disabled) return;
    onSubmit(trimmed);
    setText("");
    setPickerQuery(null);
    if (ref.current) ref.current.style.height = "40px";
  }, [text, disabled, onSubmit]);

  const onKeyDown = (e: React.KeyboardEvent<HTMLTextAreaElement>): void => {
    if (e.key === "Enter" && !e.shiftKey && pickerQuery === null) {
      e.preventDefault();
      trySubmit();
    }
    // Esc and arrow keys are handled inside MentionPicker via window listener
  };

  return (
    <div ref={wrapRef} className="flex-shrink-0 border-t border-border bg-bg px-2.5 py-2 relative">
      {pickerQuery !== null && (
        <MentionPicker
          query={pickerQuery}
          onSelect={handleSelect}
          onDismiss={() => setPickerQuery(null)}
        />
      )}
      <textarea
        ref={ref}
        value={text}
        onChange={handleChange}
        onKeyDown={onKeyDown}
        placeholder={disabled ? "Sign in to post…" : "Write a note… (@ to mention)"}
        data-testid="chat-compose"
        className="w-full resize-none bg-transparent text-fg text-base font-sans
                   placeholder:text-fg-dim outline-0 border-0"
        style={{ height: "40px" }}
        role="textbox"
      />
      <div className="flex items-center justify-between text-xs text-fg-dim">
        <span>
          <kbd className="border border-border rounded px-1">⏎</kbd> send
          {" · "}
          <kbd className="border border-border rounded px-1">⇧⏎</kbd> newline
          {" · "}
          <kbd className="border border-border rounded px-1">@</kbd> mention
        </span>
      </div>
    </div>
  );
}
```

- [ ] **Step 5: Run tests — expect pass**

```bash
node_modules/.bin/vitest run test/MentionCompose.test.tsx
```

Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/MentionCompose.tsx \
        packages/HimalayaUI/frontend/test/MentionCompose.test.tsx \
        packages/HimalayaUI/frontend/package.json \
        packages/HimalayaUI/frontend/package-lock.json
git commit -m "feat(ui): add MentionCompose with @ trigger, pretext auto-resize, and picker integration"
```

---

## Task 10: Wire `MentionCompose` into `ChatCard`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ChatCard.tsx`

The `Compose` function currently renders a raw textarea. Replace it with `MentionCompose`.

- [ ] **Step 1: Replace `Compose` with `MentionCompose` in ChatCard.tsx**

In `packages/HimalayaUI/frontend/src/components/ChatCard.tsx`:

Add import:
```typescript
import { MentionCompose } from "./MentionCompose";
```

Replace the `<Compose ... />` usage:
```typescript
<MentionCompose
  disabled={username === undefined || postMsg.isPending}
  onSubmit={(body) => postMsg.mutate(body)}
/>
```

Delete the entire `Compose` function and `ComposeProps` interface — they are now replaced by `MentionCompose`.

- [ ] **Step 2: Run full test suite**

```bash
node_modules/.bin/vitest run
```

Expected: all PASS. If existing ChatCard tests that reference `chat-compose` testid fail, verify `MentionCompose` renders `<textarea data-testid="chat-compose" ...>` — it does per the implementation above.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ChatCard.tsx
git commit -m "feat(ui): wire MentionCompose into ChatCard replacing raw Compose"
```

---

## Task 11: Full test suite + build verification

- [ ] **Step 1: Run Julia backend tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: all PASS.

- [ ] **Step 2: Run frontend unit tests**

```bash
cd packages/HimalayaUI/frontend
node_modules/.bin/vitest run
```

Expected: all PASS.

- [ ] **Step 3: Run TypeScript build**

```bash
npm run build
```

Expected: `tsc --noEmit` passes, `vite build` succeeds.

- [ ] **Step 4: Run Playwright E2E**

```bash
npm run e2e
```

Expected: all PASS. If any test flakes on port, run `lsof -ti:5173 | xargs kill -9` and retry.

- [ ] **Step 5: Final commit if any fixes needed**

If any small fixes were made during verification:

```bash
git add -p
git commit -m "fix(ui): mention system post-integration fixes"
```
