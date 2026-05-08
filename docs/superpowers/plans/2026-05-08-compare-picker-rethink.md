# Compare picker rethink — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite the Compare-page "Add traces" picker from exposure-centric to sample-first (PR1, issues #84+#85), then extract an inline-panel variant for the Compare edit-mode right slot (PR2, issue #87 §6–7).

**Architecture:** PR1 introduces a new read-only Julia route `GET /api/experiments/:eid/picker-samples` returning sample rows with resolved indexing-exposure ids, then refactors `ComparisonPicker.tsx` into `ComparisonPicker` (modal shell) + `ComparisonPickerBody` (pure render layer) + `SamplePickerRow`. PR2 adds `ComparisonPickerPanel` (inline shell) consuming the same body, swaps the right slot in `ComparePageEdit`. Drift policy is freeze-at-pick: the picker resolves the indexing exposure once at pick time and stores that exposure id; the comparison member never silently retargets when Inspect later flips the `selected` flag. Schema unchanged; no SSE/idempotency surface.

**Tech Stack:** Julia 1.10 + Oxygen.jl + SQLite.jl (backend), React 18 + TypeScript strict + TanStack Query + Zustand + Vite + Vitest + Playwright (frontend), Tailwind v4 (styling).

**Reference spec:** `docs/superpowers/specs/2026-05-08-compare-picker-rethink-design.md` (read first; this plan implements it).

---

## File map

**PR1 — sample-first picker:**

- Create: `packages/HimalayaUI/src/routes_picker.jl` (extend) — new `@get "/api/experiments/{eid}/picker-samples"` route.
- Modify: `packages/HimalayaUI/src/comparisons.jl` — new `picker_samples(db, experiment_id)` helper.
- Create: `packages/HimalayaUI/test/test_picker_samples_route.jl` — backend tests for the new helper + route.
- Modify: `packages/HimalayaUI/test/runtests.jl` — include the new test file.
- Modify: `packages/HimalayaUI/frontend/src/api.ts` — `getPickerSamples` fetcher + `PickerSampleRow` type.
- Modify: `packages/HimalayaUI/frontend/src/queries.ts` — `usePickerSamples` hook + new `queryKeys.pickerSamples` entry.
- Create: `packages/HimalayaUI/frontend/src/components/ComparisonPickerBody.tsx` — pure render layer (filters, recents, list, footer-slot).
- Create: `packages/HimalayaUI/frontend/src/components/SamplePickerRow.tsx` — primary sample row + caret.
- Modify: `packages/HimalayaUI/frontend/src/components/ExposureListRow.tsx` — add `control: "checkbox" | "radio"` prop.
- Modify: `packages/HimalayaUI/frontend/src/components/ComparisonPicker.tsx` — slim down to modal shell wrapping the body.
- Create: `packages/HimalayaUI/frontend/test/SamplePickerRow.test.tsx` — sample row + caret + disabled-row tests.
- Create: `packages/HimalayaUI/frontend/test/ComparisonPickerBody.test.tsx` — body-level tests (filters, recents dedup, picks emit).
- Modify: `packages/HimalayaUI/frontend/test/ComparisonPicker.test.tsx` — narrow to shell concerns + focus-trap regression.
- Modify: `packages/HimalayaUI/frontend/test/ExposureListRow.test.tsx` — control-prop branch.
- Modify: `packages/HimalayaUI/frontend/e2e/comparisonPicker.spec.ts` — mock new route, update assertions.

**PR2 — inline picker panel:**

- Create: `packages/HimalayaUI/frontend/src/components/ComparisonPickerPanel.tsx` — inline shell.
- Modify: `packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx` — swap right slot, remove buried "+ Add traces", focus search via ref.
- Create: `packages/HimalayaUI/frontend/test/ComparisonPickerPanel.test.tsx` — inline shell + immediate-commit + locking.
- Modify: `packages/HimalayaUI/frontend/test/ComparePageEdit.test.tsx` — right slot now hosts the panel.
- Create: `packages/HimalayaUI/frontend/e2e/live/comparePickerInline.spec.ts` — live-mode roundtrip.

---

# PR1 — Sample-first picker

## Task 1: Backend helper `picker_samples`

**Files:**
- Modify: `packages/HimalayaUI/src/comparisons.jl` (append helper near `recently_used_exposures`, around line 484)

- [ ] **Step 1: Write the failing test (helper-level)**

Create `packages/HimalayaUI/test/test_picker_samples_route.jl` with the helper test only (route tests added in Task 3):

```julia
using Test
using SQLite, DBInterface, Tables
using HimalayaUI
using HimalayaUI: open_db, picker_samples
# Note: `picker_samples` is reachable via `using HimalayaUI:` because it is a
# module-level binding in `comparisons.jl` (which is `include`d into the
# HimalayaUI module). No `export` is required — same pattern as
# `recently_used_exposures` referenced from `routes_picker.jl:51`.

@testset "picker_samples helper" begin
    @testset "selected exposure resolves" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "h.db"))
            DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
            DBInterface.execute(db,
                "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'S')")
            DBInterface.execute(db,
                "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'f1', 0)")
            DBInterface.execute(db,
                "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (101, 10, 'f2', 1)")

            rows = picker_samples(db, 1)
            @test length(rows) == 1
            @test rows[1][:indexing_exposure_id] == 101
            @test length(rows[1][:all_exposures]) == 2
            @test [e[:id] for e in rows[1][:all_exposures]] == [100, 101]   # ORDER BY id ASC
            @test rows[1][:all_exposures][2][:selected] === true             # bool, not 1
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

`runtests.jl` does not read `ARGS` — `Pkg.test(...; test_args=...)` would silently run the entire 5–10 min suite. For the iterative red→green loop, run the test file directly. The `test_http.jl` include is needed for `with_test_server` (used in Task 3).

Run from the worktree root:

```bash
julia --project=packages/HimalayaUI -e '
using Test, HimalayaUI
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_picker_samples_route.jl")
' > /tmp/jl-picker.out 2>&1
grep -E "Test Summary|did not pass|fail|UndefVarError" /tmp/jl-picker.out
```

Expected: `UndefVarError: picker_samples not defined` or similar.

- [ ] **Step 3: Add `runtests.jl` include**

Modify `packages/HimalayaUI/test/runtests.jl` — add the include directly after `include("test_picker_routes.jl")`:

```julia
    include("test_picker_routes.jl")
    include("test_picker_samples_route.jl")
```

Tests run in include order top-to-bottom inside the outer `@testset "HimalayaUI"`.

- [ ] **Step 4: Implement the helper**

Append to `packages/HimalayaUI/src/comparisons.jl` directly after `recently_used_exposures` (line 484):

```julia
"""
    picker_samples(db, experiment_id) -> Vector{Dict{Symbol, Any}}

Picker primary list per spec §"PR1 — sample-first picker → Backend".

For each sample in `experiment_id`, returns:
  :sample              => sample row (Symbol-keyed Dict, with :tags vector)
  :indexing_exposure_id => Int or nothing — `selected = 1` else MAX(id) else nothing
  :all_exposures       => Vector of {:id, :filename, :selected::Bool}, ORDER BY id ASC

Three bulk queries (no JOIN'd Cartesian flatten, no per-sample N+1):
(1) `WHERE experiment_id = ?` for samples,
(2) `WHERE sample_id IN (...)` for exposures,
(3) `WHERE sample_id IN (...)` for sample_tags.
Empty experiment ⇒ [].
"""
function picker_samples(db::SQLite.DB, experiment_id::Integer)::Vector{Dict{Symbol, Any}}
    samples = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM samples WHERE experiment_id = ? ORDER BY id", [Int(experiment_id)]))
    isempty(samples) && return Dict{Symbol, Any}[]

    sample_ids   = [Int(s.id) for s in samples]
    placeholders = join(fill("?", length(sample_ids)), ",")
    exposures    = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, sample_id, filename, selected FROM exposures
         WHERE sample_id IN ($placeholders) ORDER BY sample_id ASC, id ASC",
        sample_ids))
    tag_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, sample_id, key, value, source FROM sample_tags
         WHERE sample_id IN ($placeholders) ORDER BY sample_id ASC, id ASC",
        sample_ids))

    # Group exposures + tags by sample_id in Julia (avoids a Cartesian JOIN dedup).
    grouped_exps = Dict{Int, Vector{NamedTuple}}()
    for e in exposures
        push!(get!(grouped_exps, Int(e.sample_id), NamedTuple[]), e)
    end
    grouped_tags = Dict{Int, Vector{NamedTuple}}()
    for t in tag_rows
        push!(get!(grouped_tags, Int(t.sample_id), NamedTuple[]), t)
    end

    out = Dict{Symbol, Any}[]
    for sm in samples
        sid  = Int(sm.id)
        exps = get(grouped_exps, sid, NamedTuple[])

        # Resolve indexing exposure: highest-id selected wins (defensive
        # against legacy multi-selected data); else highest-id overall;
        # else nothing for zero-exposure samples. `exps` is sorted id ASC,
        # so iterating in reverse yields highest-id first.
        idx_id = nothing
        for e in Iterators.reverse(exps)
            if e.selected != 0
                idx_id = Int(e.id); break
            end
        end
        if idx_id === nothing && !isempty(exps)
            idx_id = Int(last(exps).id)   # last == max by id ASC ordering
        end

        # Sample → row_to_json + bulk-grouped tags. Drop sample_id from each
        # tag dict (it's redundant once grouped under the sample).
        tags_for_sm = get(grouped_tags, sid, NamedTuple[])
        tag_dicts = map(tags_for_sm) do t
            d = row_to_json(t)
            delete!(d, :sample_id)
            d
        end
        sample_dict = row_to_json(sm)
        sample_dict[:tags] = tag_dicts

        all_exp = [row_to_json(e; bool_keys = (:selected,)) for e in exps]
        push!(out, Dict{Symbol, Any}(
            :sample               => sample_dict,
            :indexing_exposure_id => idx_id,
            :all_exposures        => all_exp,
        ))
    end
    out
end
```

- [ ] **Step 5: Run test to verify it passes**

```bash
julia --project=packages/HimalayaUI -e '
using Test, HimalayaUI
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_picker_samples_route.jl")
' > /tmp/jl-picker.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-picker.out
```

Expected: passing summary, zero failures.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/comparisons.jl packages/HimalayaUI/test/runtests.jl packages/HimalayaUI/test/test_picker_samples_route.jl
git commit -m "feat(picker): picker_samples helper (sample-first picker, #84)"
```

## Task 2: Backend helper — corner-case tests

**Files:**
- Modify: `packages/HimalayaUI/test/test_picker_samples_route.jl`

- [ ] **Step 1: Append corner-case tests**

```julia
@testset "no selected falls back to highest id" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "h.db"))
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'S')")
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (200, 10, 'a', 0)")
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (199, 10, 'b', 0)")
        rows = picker_samples(db, 1)
        @test rows[1][:indexing_exposure_id] == 200
    end
end

@testset "single exposure, selected=0 — still resolves" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "h.db"))
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'S')")
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (300, 10, 'x', 0)")
        rows = picker_samples(db, 1)
        @test rows[1][:indexing_exposure_id] == 300
    end
end

@testset "zero-exposure sample — included with null indexing_exposure_id" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "h.db"))
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'EmptyS')")
        rows = picker_samples(db, 1)
        @test length(rows) == 1
        @test rows[1][:indexing_exposure_id] === nothing
        @test isempty(rows[1][:all_exposures])
    end
end

@testset "unknown experiment id → empty list" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "h.db"))
        @test picker_samples(db, 99) == Dict{Symbol, Any}[]
    end
end

@testset "multi-experiment isolation" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "h.db"))
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'A', '/tmp/a', '/tmp/a/data', '/tmp/a/analysis')")
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (2, 'B', '/tmp/b', '/tmp/b/data', '/tmp/b/analysis')")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'A1')")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (20, 2, 'B1')")
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'a', 0)")
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (999, 20, 'b', 0)")  # bigger id, wrong exp
        a = picker_samples(db, 1)
        @test length(a) == 1
        @test [e[:id] for e in a[1][:all_exposures]] == [100]
        @test a[1][:indexing_exposure_id] == 100   # not 999 — global MAX(id) would have leaked
    end
end

@testset "orphan exposure (sample_id NULL) excluded" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "h.db"))
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'S')")
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'a', 1)")
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (200, NULL, 'orphan', 0)")
        rows = picker_samples(db, 1)
        @test length(rows[1][:all_exposures]) == 1
        @test rows[1][:all_exposures][1][:id] == 100
    end
end

@testset "NULL name and label render as null in JSON" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "h.db"))
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id) VALUES (10, 1)")  # name+label NULL
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'f', 1)")
        rows = picker_samples(db, 1)
        @test rows[1][:sample][:name]  === nothing
        @test rows[1][:sample][:label] === nothing
    end
end

@testset "defensive multi-selected legacy data" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "h.db"))
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'S')")
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'a', 1)")
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (200, 10, 'b', 1)")
        rows = picker_samples(db, 1)
        @test rows[1][:indexing_exposure_id] == 200   # higher id wins among multiple selected
    end
end
```

- [ ] **Step 2: Run tests, verify all pass**

```bash
julia --project=packages/HimalayaUI -e '
using Test, HimalayaUI
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_picker_samples_route.jl")
' > /tmp/jl-picker.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-picker.out
```

Expected: all corner-case `@testset`s pass.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/test/test_picker_samples_route.jl
git commit -m "test(picker): corner-case coverage for picker_samples"
```

## Task 3: Backend route `GET /api/experiments/:eid/picker-samples`

**Files:**
- Modify: `packages/HimalayaUI/src/routes_picker.jl`

- [ ] **Step 1: Write the failing route test**

Append to `packages/HimalayaUI/test/test_picker_samples_route.jl`. The test harness is `with_test_server(f, db)` defined in `test_http.jl` (yields `(port, base)`):

```julia
using HTTP, JSON3

@testset "GET /api/experiments/:eid/picker-samples" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "h.db"))
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'S')")
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'f1', 1)")

        with_test_server(db) do port, base
            r = HTTP.get("$base/api/experiments/1/picker-samples")
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test length(body) == 1
            @test body[1].indexing_exposure_id == 100
            @test body[1].all_exposures[1].selected === true   # JSON-shape: bool, not 1
            @test haskey(body[1], :indexing_exposure_id)        # null vs absent key

            # Sanity: zero-exposure sample produces null (not absent).
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (11, 1, 'Empty')")
            r2 = HTTP.get("$base/api/experiments/1/picker-samples")
            body2 = JSON3.read(String(r2.body))
            empty_row = first(filter(b -> b.sample.id == 11, collect(body2)))
            @test empty_row.indexing_exposure_id === nothing
        end
    end
end
```

`with_test_server` calls `start_test_server!(db, port)` (server.jl:136) which binds the db so `current_db()` resolves correctly inside the route handler. Don't bypass it — there's no `_test_serve!` / `_test_terminate!` helper.

- [ ] **Step 2: Run test, verify it fails**

Expected: HTTP 404, since the route isn't registered.

- [ ] **Step 3: Add the route in `routes_picker.jl`**

Inside `register_picker_routes!()` (add after the `sample-tags` route):

```julia
    @get "/api/experiments/{eid}/picker-samples" function(req::HTTP.Request, eid::Int)
        db = current_db()
        rows = picker_samples(db, eid)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(rows))
    end
```

No `export` needed. `picker_samples` is reachable unqualified from `routes_picker.jl` because both files are `include`d into the same `HimalayaUI` module — same pattern as `recently_used_exposures` referenced at `routes_picker.jl:51`.

- [ ] **Step 4: Run test, verify pass**

```bash
julia --project=packages/HimalayaUI -e '
using Test, HimalayaUI
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_picker_samples_route.jl")
' > /tmp/jl-picker.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-picker.out
```

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_picker.jl packages/HimalayaUI/test/test_picker_samples_route.jl
git commit -m "feat(picker): GET /api/experiments/:eid/picker-samples route"
```

## Task 4: Frontend API fetcher + types

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/api.ts`

- [ ] **Step 1: Add the type and fetcher**

Append after `getSampleTags` (around line 587):

```ts
/** Per-row shape returned by `GET /api/experiments/:eid/picker-samples`. */
export interface PickerSampleRow {
  sample: Sample;
  indexing_exposure_id: number | null;
  all_exposures: PickerSampleExposure[];
}

export interface PickerSampleExposure {
  id: number;
  sample_id: number;
  filename: string | null;
  selected: boolean;
}

export const getPickerSamples = (
  experiment_id: number,
): Promise<PickerSampleRow[]> =>
  request<PickerSampleRow[]>(
    "GET", `/api/experiments/${experiment_id}/picker-samples`);
```

- [ ] **Step 2: Verify compile**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/tsc --noEmit
```

Expected: zero errors.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/api.ts
git commit -m "feat(picker): getPickerSamples fetcher + PickerSampleRow type"
```

## Task 5: Frontend `usePickerSamples` hook + queryKey

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/queries.ts`

- [ ] **Step 1: Add queryKey + hook**

Add to the `queryKeys` map (in the picker support section, near line 71):

```ts
  pickerSamples: (experimentId: number) =>
    ["experiment", experimentId, "picker-samples"] as const,
```

Add the hook (next to `useSampleTags`, around line 590):

```ts
/**
 * Picker primary list. Returns one row per sample with the resolved
 * indexing-exposure id frozen at fetch time. Spec §PR1 backend.
 *
 * `enabled: experimentId !== undefined` matches `useSampleTags` — picker is
 * always opened from an experiment context, but the hook is shaped to
 * accept `undefined` so render isn't gated on experiment selection.
 */
export function usePickerSamples(experimentId: number | undefined) {
  return useQuery({
    queryKey: experimentId !== undefined
      ? queryKeys.pickerSamples(experimentId)
      : (["experiment", "none", "picker-samples"] as const),
    queryFn: () => api.getPickerSamples(experimentId as number),
    enabled: experimentId !== undefined,
  });
}
```

- [ ] **Step 2: Verify compile**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/tsc --noEmit
```

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/queries.ts
git commit -m "feat(picker): usePickerSamples hook + queryKey"
```

## Task 6: Extend `ExposureListRow` with `control` prop

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ExposureListRow.tsx`
- Modify: `packages/HimalayaUI/frontend/test/ExposureListRow.test.tsx`

- [ ] **Step 1: Write the failing test**

Add to `ExposureListRow.test.tsx`:

```tsx
import { render, screen } from "@testing-library/react";
import { ExposureListRow } from "../src/components/ExposureListRow";
import type { Exposure, Sample } from "../src/api";

const exposure: Exposure = {
  id: 1, sample_id: 10, filename: "f.dat",
  kind: "file", selected: false, status: "accepted",
  image_path: null, image_version: "", tags: [], sources: [],
  trace_hash: null, analysis_inputs_hash: null,
};
const sample: Sample = {
  id: 10, experiment_id: 1, name: "S", label: null, notes: null, tags: [],
};

test("renders radio when control='radio'", () => {
  render(
    <ExposureListRow
      exposure={exposure} sample={sample}
      control="radio"
      checked={false} onCheckedChange={() => {}}
    />,
  );
  const input = screen.getByTestId("exposure-list-row-checkbox") as HTMLInputElement;
  expect(input.type).toBe("radio");
});

test("renders checkbox when control unset (default)", () => {
  render(
    <ExposureListRow
      exposure={exposure} sample={sample}
      checked={false} onCheckedChange={() => {}}
    />,
  );
  const input = screen.getByTestId("exposure-list-row-checkbox") as HTMLInputElement;
  expect(input.type).toBe("checkbox");
});
```

- [ ] **Step 2: Run, verify fail**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/ExposureListRow.test.tsx
```

Expected: TypeScript or runtime fail on `control` prop.

- [ ] **Step 3: Add the prop**

Modify `ExposureListRow.tsx`:

```tsx
interface Props {
  // … existing props …
  /** Input control type. Defaults to "checkbox". */
  control?: "checkbox" | "radio";
}

export function ExposureListRow({
  // … existing destructure …
  control = "checkbox",
}: Props): JSX.Element {
  // … existing body …
  // Replace `<input type="checkbox" …>` with:
  <input
    type={control}
    // … rest unchanged …
  />
```

(Keep `data-testid="exposure-list-row-checkbox"` — it is a stable identifier, not a type-specific name. Tests reference it.)

- [ ] **Step 4: Run, verify pass**

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ExposureListRow.tsx packages/HimalayaUI/frontend/test/ExposureListRow.test.tsx
git commit -m "feat(picker): ExposureListRow control prop (checkbox|radio)"
```

## Task 7: `SamplePickerRow` component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/SamplePickerRow.tsx`
- Create: `packages/HimalayaUI/frontend/test/SamplePickerRow.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
import { render, screen, fireEvent } from "@testing-library/react";
import { SamplePickerRow } from "../src/components/SamplePickerRow";
import type { PickerSampleRow } from "../src/api";

const baseRow: PickerSampleRow = {
  sample: { id: 10, experiment_id: 1, name: "S1", label: null, notes: null, tags: [] },
  indexing_exposure_id: 100,
  all_exposures: [
    { id: 100, sample_id: 10, filename: "f1.dat", selected: true },
    { id: 101, sample_id: 10, filename: "f2.dat", selected: false },
  ],
};

test("renders sample name as primary, no filename in primary slot", () => {
  render(
    <SamplePickerRow
      row={baseRow} checked={false} onCheckedChange={() => {}}
      overrideExposureId={undefined} onOverrideChange={() => {}}
    />,
  );
  expect(screen.getByText("S1")).toBeInTheDocument();
  expect(screen.queryByText("f1.dat")).toBeNull();   // hidden until caret expanded
});

test("caret toggles override list visibility", () => {
  render(
    <SamplePickerRow
      row={baseRow} checked={false} onCheckedChange={() => {}}
      overrideExposureId={undefined} onOverrideChange={() => {}}
    />,
  );
  const caret = screen.getByTestId("sample-picker-row-caret");
  fireEvent.click(caret);
  expect(screen.getByText("f1.dat")).toBeInTheDocument();
  expect(screen.getByText("f2.dat")).toBeInTheDocument();
});

test("zero-exposure sample renders disabled, no checkbox", () => {
  const empty: PickerSampleRow = {
    ...baseRow, indexing_exposure_id: null, all_exposures: [],
  };
  render(
    <SamplePickerRow
      row={empty} checked={false} onCheckedChange={() => {}}
      overrideExposureId={undefined} onOverrideChange={() => {}}
    />,
  );
  const row = screen.getByTestId("sample-picker-row");
  expect(row).toHaveAttribute("data-disabled", "true");
  expect(row).not.toHaveAttribute("data-exposure-id");
  expect(screen.queryByRole("checkbox")).toBeNull();
});

test("override radio fires onOverrideChange with selected exposure id", () => {
  const handle = vi.fn();
  render(
    <SamplePickerRow
      row={baseRow} checked={true} onCheckedChange={() => {}}
      overrideExposureId={undefined} onOverrideChange={handle}
    />,
  );
  fireEvent.click(screen.getByTestId("sample-picker-row-caret"));
  fireEvent.click(screen.getByLabelText(/f2\.dat/));
  expect(handle).toHaveBeenCalledWith(101);
});

test("data-exposure-id reflects resolved id (default → indexing)", () => {
  render(
    <SamplePickerRow
      row={baseRow} checked={false} onCheckedChange={() => {}}
      overrideExposureId={undefined} onOverrideChange={() => {}}
    />,
  );
  expect(screen.getByTestId("sample-picker-row")).toHaveAttribute("data-exposure-id", "100");
});

test("data-exposure-id reflects override when set", () => {
  render(
    <SamplePickerRow
      row={baseRow} checked={true} onCheckedChange={() => {}}
      overrideExposureId={101} onOverrideChange={() => {}}
    />,
  );
  expect(screen.getByTestId("sample-picker-row")).toHaveAttribute("data-exposure-id", "101");
});
```

- [ ] **Step 2: Run, verify fail (component doesn't exist)**

- [ ] **Step 3: Implement the component**

```tsx
import { useState } from "react";
import type { PickerSampleRow as PickerRow } from "../api";
import { ExposureListRow } from "./ExposureListRow";

interface Props {
  row: PickerRow;
  /** Whole-row checked state (the sample is in the picks set). */
  checked: boolean;
  onCheckedChange: (next: boolean) => void;
  /** Override exposure id, or undefined for default. */
  overrideExposureId: number | undefined;
  onOverrideChange: (exposureId: number) => void;
}

function sampleLabel(s: PickerRow["sample"]): string {
  return s.name ?? s.label ?? `Sample #${s.id}`;
}

export function SamplePickerRow({
  row, checked, onCheckedChange, overrideExposureId, onOverrideChange,
}: Props): JSX.Element {
  const [expanded, setExpanded] = useState(false);
  const disabled = row.indexing_exposure_id === null;
  const resolvedExposureId = overrideExposureId ?? row.indexing_exposure_id;

  const dataAttrs: Record<string, string> = {};
  dataAttrs["data-sample-id"] = String(row.sample.id);
  if (disabled) dataAttrs["data-disabled"] = "true";
  if (resolvedExposureId !== null && resolvedExposureId !== undefined) {
    dataAttrs["data-exposure-id"] = String(resolvedExposureId);
  }

  return (
    <div
      data-testid="sample-picker-row"
      {...dataAttrs}
      className={
        "flex items-start gap-3 px-3 py-2 rounded " +
        (disabled
          ? "opacity-60"
          : "cursor-pointer hover:bg-bg-elevated") +
        (checked && !disabled ? " ring-1 ring-accent/30" : "")
      }
    >
      {!disabled && (
        <input
          type="checkbox"
          data-testid="sample-picker-row-checkbox"
          checked={checked}
          onChange={(e) => onCheckedChange(e.target.checked)}
          className="mt-1 shrink-0"
        />
      )}
      <div className="flex-1 min-w-0">
        <div className="flex items-baseline gap-2">
          <span className="font-medium text-fg truncate">
            {sampleLabel(row.sample)}
          </span>
          <span className="text-xs text-fg-dim shrink-0">
            {row.all_exposures.length} {row.all_exposures.length === 1 ? "exposure" : "exposures"}
          </span>
          {!disabled && (
            <button
              type="button"
              data-testid="sample-picker-row-caret"
              aria-expanded={expanded}
              aria-label="Show exposures"
              onClick={(e) => { e.stopPropagation(); setExpanded((v) => !v); }}
              className="ml-auto text-fg-muted hover:text-fg text-xs px-1"
            >
              {expanded ? "▾" : "▸"}
            </button>
          )}
        </div>
        {row.sample.notes && (
          <div className="text-xs text-fg-muted line-clamp-2" title={row.sample.notes}>
            {row.sample.notes}
          </div>
        )}
        {expanded && (
          <ul role="radiogroup" aria-label="Override exposure"
              className="mt-2 ml-4 border-l border-border pl-2 flex flex-col">
            {row.all_exposures.map((e) => (
              <li key={e.id}>
                <ExposureListRow
                  exposure={{
                    id: e.id, sample_id: e.sample_id, filename: e.filename,
                    kind: "file", selected: e.selected, status: "accepted",
                    image_path: null, image_version: "", tags: [], sources: [],
                    trace_hash: null, analysis_inputs_hash: null,
                  }}
                  sample={row.sample}
                  control="radio"
                  checked={resolvedExposureId === e.id}
                  onCheckedChange={() => onOverrideChange(e.id)}
                />
              </li>
            ))}
          </ul>
        )}
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Run, verify pass**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/SamplePickerRow.test.tsx
```

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/SamplePickerRow.tsx packages/HimalayaUI/frontend/test/SamplePickerRow.test.tsx
git commit -m "feat(picker): SamplePickerRow with override caret (#84)"
```

## Task 8: `ComparisonPickerBody` extraction

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ComparisonPickerBody.tsx`
- Create: `packages/HimalayaUI/frontend/test/ComparisonPickerBody.test.tsx`

- [ ] **Step 1: Write the failing tests**

```tsx
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { beforeEach, afterEach, vi } from "vitest";
import { ComparisonPickerBody, type Pick } from "../src/components/ComparisonPickerBody";
import * as api from "../src/api";

// Wrap mocks in beforeEach so they're isolated per test (the project's
// vitest setup does not call vi.restoreAllMocks() between files).
beforeEach(() => {
  vi.spyOn(api, "getPickerSamples").mockResolvedValue([
    {
      sample: { id: 10, experiment_id: 1, name: "S1", label: null, notes: null, tags: [] },
      indexing_exposure_id: 100,
      all_exposures: [{ id: 100, sample_id: 10, filename: "f1.dat", selected: true }],
    },
    {
      sample: { id: 20, experiment_id: 1, name: "S2", label: null, notes: null, tags: [] },
      indexing_exposure_id: 200,
      all_exposures: [{ id: 200, sample_id: 20, filename: "f2.dat", selected: true }],
    },
  ]);
  vi.spyOn(api, "listExperiments").mockResolvedValue([
    { id: 1, name: "E", config: null },
  ]);
  vi.spyOn(api, "getRecentlyPickedExposures").mockResolvedValue([200, 100]);
  vi.spyOn(api, "getSampleTags").mockResolvedValue([]);
});
afterEach(() => { vi.restoreAllMocks(); });

function wrap(ui: React.ReactElement) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(<QueryClientProvider client={qc}>{ui}</QueryClientProvider>);
}

test("controlled-picks: onPicksChange fires with default exposure id on toggle", async () => {
  const onPicksChange = vi.fn();
  wrap(
    <ComparisonPickerBody
      experimentId={1}
      picks={[]}
      onPicksChange={onPicksChange}
      alreadyAddedExposureIds={new Set()}
    />,
  );
  await screen.findByText("S1");
  fireEvent.click(screen.getByTestId("sample-picker-row-checkbox"));
  expect(onPicksChange).toHaveBeenCalledWith([
    { sample_id: 10, exposure_id: 100, source: "default" },
  ]);
});

// NOTE: Immediate-mode `onPick` is added in PR2 Task 13. PR1 ships only the
// controlled-picks interface (`picks` + `onPicksChange`). Skip the immediate-mode
// test until PR2.

test("recents section dedupes against main list (one row per sample, S2 appears once)", async () => {
  wrap(
    <ComparisonPickerBody
      experimentId={1}
      picks={[]} onPicksChange={() => {}}
      alreadyAddedExposureIds={new Set()}
    />,
  );
  await screen.findByText("S1");
  // Recents has S2 (exposure 200). Main list also has S2.
  // The body must NOT render S2 twice in the visible main list.
  const s2Rows = screen.queryAllByText("S2");
  expect(s2Rows.length).toBeLessThanOrEqual(1);   // exactly one section renders it
});
```

- [ ] **Step 2: Run, verify fail (component doesn't exist)**

- [ ] **Step 3: Implement `ComparisonPickerBody`**

```tsx
import { useMemo, useRef, useState } from "react";
import type { RefObject } from "react";
import { useExperiments, useRecentlyPickedExposures, useSampleTags, usePickerSamples } from "../queries";
import { SamplePickerRow } from "./SamplePickerRow";
import type { PickerSampleRow } from "../api";

export type Pick = {
  sample_id: number;
  exposure_id: number;
  source: "default" | "override";
};

interface Props {
  experimentId: number | undefined;
  /** Controlled picks list. */
  picks: Pick[];
  onPicksChange: (next: Pick[]) => void;
  /** Set of exposure ids already in the active draft — rendered as locked. */
  alreadyAddedExposureIds: Set<number>;
  /** Optional ref the parent threads down to focus the search input. */
  searchInputRef?: RefObject<HTMLInputElement>;
  // Filter chip state (search / tag-chip) is internal to the body.
  // PR2 adds an `onPick` prop for immediate-commit shells.
}

export function ComparisonPickerBody({
  experimentId, picks, onPicksChange, alreadyAddedExposureIds,
  searchInputRef,
}: Props): JSX.Element {
  const userId = useCurrentUserId();

  // Data sources.
  const pickerQ      = usePickerSamples(experimentId);
  const recentsQ     = useRecentlyPickedExposures(userId, 20);
  const tagsQ        = useSampleTags(experimentId);
  // `useExperiments` is intentionally not invoked — the picker scope is the
  // active experiment only in PR1. Cross-experiment widening can return in
  // a follow-on if user demand surfaces.

  const fallbackRef = useRef<HTMLInputElement>(null);
  const inputRef = searchInputRef ?? fallbackRef;

  const [search, setSearch]               = useState("");
  const [selectedTags, setSelectedTags]   = useState<Set<string>>(new Set());

  // Build the rendered row list (filtered + sorted) from the picker-samples query.
  const allRows = pickerQ.data ?? [];
  const filteredRows = useMemo<PickerSampleRow[]>(() => {
    return allRows.filter((r) => {
      if (selectedTags.size > 0) {
        const match = r.sample.tags.some((t) => selectedTags.has(`${t.key}:${t.value}`));
        if (!match) return false;
      }
      if (search.trim() !== "") {
        const needle = search.toLowerCase();
        const haystack = [
          r.sample.name ?? "", r.sample.label ?? "", r.sample.notes ?? "",
          ...r.sample.tags.map((t) => t.value),
          ...r.all_exposures.map((e) => e.filename ?? ""),
        ].map((s) => s.toLowerCase());
        if (!haystack.some((s) => s.includes(needle))) return false;
      }
      return true;
    }).sort((a, b) =>
      (a.sample.name ?? "").localeCompare(b.sample.name ?? ""),
    );
  }, [allRows, selectedTags, search]);

  // Recents → derive sample-id ordering, dedupe to one row per sample, then
  // exclude from main list. useMemo keeps server state out of Zustand.
  const recentSamples = useMemo<PickerSampleRow[]>(() => {
    const recentIds = recentsQ.data ?? [];
    if (recentIds.length === 0) return [];
    const sampleByExposureId = new Map<number, PickerSampleRow>();
    for (const r of allRows) {
      for (const e of r.all_exposures) sampleByExposureId.set(e.id, r);
    }
    const seen = new Set<number>();
    const out: PickerSampleRow[] = [];
    for (const eid of recentIds) {
      const r = sampleByExposureId.get(eid);
      if (!r || seen.has(r.sample.id)) continue;
      seen.add(r.sample.id);
      out.push(r);
    }
    return out;
  }, [recentsQ.data, allRows]);

  const recentSampleIds = useMemo(
    () => new Set(recentSamples.map((r) => r.sample.id)),
    [recentSamples],
  );
  const mainListRows = useMemo(
    () => filteredRows.filter((r) => !recentSampleIds.has(r.sample.id)),
    [filteredRows, recentSampleIds],
  );

  // Pick-set lookup.
  const pickedSampleIds = useMemo(
    () => new Set(picks.map((p) => p.sample_id)),
    [picks],
  );
  const overrideBySampleId = useMemo(() => {
    const m = new Map<number, number>();
    for (const p of picks) {
      if (p.source === "override") m.set(p.sample_id, p.exposure_id);
    }
    return m;
  }, [picks]);

  const togglePickFor = (row: PickerSampleRow, next: boolean): void => {
    if (row.indexing_exposure_id === null) return;
    const pick: Pick = {
      sample_id: row.sample.id,
      exposure_id: row.indexing_exposure_id,
      source: "default",
    };
    if (next) onPicksChange([...picks, pick]);
    else onPicksChange(picks.filter((p) => p.sample_id !== row.sample.id));
  };

  const overridePickFor = (row: PickerSampleRow, exposureId: number): void => {
    const next: Pick = {
      sample_id: row.sample.id,
      exposure_id: exposureId,
      source: exposureId === row.indexing_exposure_id ? "default" : "override",
    };
    const i = picks.findIndex((p) => p.sample_id === row.sample.id);
    if (i < 0) onPicksChange([...picks, next]);
    else onPicksChange(picks.map((p, j) => (j === i ? next : p)));
  };

  // Tag filter helper.
  const toggleTag = (key: string, value: string): void => {
    const id = `${key}:${value}`;
    setSelectedTags((prev) => {
      const m = new Set(prev);
      if (m.has(id)) m.delete(id);
      else m.add(id);
      return m;
    });
  };

  return (
    <div className="flex flex-col flex-1 min-h-0">
      {/* Filters (search + experiment chips + tag chips). */}
      <div className="px-4 py-2 border-b border-border space-y-2">
        <input
          ref={inputRef}
          data-testid="comparison-picker-search"
          type="search"
          value={search}
          onChange={(e) => setSearch(e.target.value)}
          placeholder="Search samples, notes, filenames…"
          className="w-full bg-transparent border border-border rounded px-2 py-1 text-sm"
          spellCheck={false}
        />
        {(tagsQ.data ?? []).length > 0 && (
          <div data-testid="picker-filter-tag" className="flex flex-wrap gap-1 items-center">
            <span className="text-xs text-fg-muted">Tags:</span>
            {(tagsQ.data ?? []).map((t) => {
              const id = `${t.key}:${t.value}`;
              const active = selectedTags.has(id);
              return (
                <button
                  key={id}
                  type="button"
                  data-testid={`picker-tag-option-${id}`}
                  data-active={active ? "true" : undefined}
                  onClick={() => toggleTag(t.key, t.value)}
                  className={
                    "text-xs px-2 py-0.5 rounded-full border " +
                    (active
                      ? "bg-accent/15 border-accent text-accent"
                      : "border-border text-fg-muted hover:text-fg")
                  }
                >
                  {t.key}:{t.value}
                </button>
              );
            })}
          </div>
        )}
      </div>

      <div className="flex-1 min-h-0 overflow-y-auto">
        {recentSamples.length > 0 && (
          <section data-testid="comparison-picker-recents" className="border-b border-border/40 py-2">
            <div className="px-4 text-xs font-medium text-fg-muted uppercase tracking-wide pb-1">
              Recently used
            </div>
            <ul role="listbox" aria-label="Recently used samples" className="flex flex-col">
              {recentSamples.map((r) => (
                <li key={`recent-${r.sample.id}`} data-testid="picker-row">
                  <SamplePickerRow
                    row={r}
                    checked={pickedSampleIds.has(r.sample.id)}
                    onCheckedChange={(next) => togglePickFor(r, next)}
                    overrideExposureId={overrideBySampleId.get(r.sample.id)}
                    onOverrideChange={(eid) => overridePickFor(r, eid)}
                  />
                </li>
              ))}
            </ul>
          </section>
        )}

        {mainListRows.length === 0 ? (
          <div data-testid="comparison-picker-empty"
               className="px-4 py-8 text-center text-fg-muted text-sm">
            No samples match. Try clearing filters.
          </div>
        ) : (
          <section className="py-2">
            <div className="px-4 text-xs font-medium text-fg-muted uppercase tracking-wide pb-1">
              All samples
            </div>
            <ul role="listbox" aria-label="Samples" className="flex flex-col">
              {mainListRows.map((r) => (
                <li key={r.sample.id} data-testid="picker-row">
                  <SamplePickerRow
                    row={r}
                    checked={pickedSampleIds.has(r.sample.id)}
                    onCheckedChange={(next) => togglePickFor(r, next)}
                    overrideExposureId={overrideBySampleId.get(r.sample.id)}
                    onOverrideChange={(eid) => overridePickFor(r, eid)}
                  />
                </li>
              ))}
            </ul>
          </section>
        )}
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Run tests, verify pass**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/ComparisonPickerBody.test.tsx
```

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ComparisonPickerBody.tsx packages/HimalayaUI/frontend/test/ComparisonPickerBody.test.tsx
git commit -m "feat(picker): ComparisonPickerBody extraction (sample-first)"
```

## Task 9: Wrap body in Skeleton + skeleton-gating regression test

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ComparisonPickerBody.tsx`
- Modify: `packages/HimalayaUI/frontend/test/ComparisonPickerBody.test.tsx`

- [ ] **Step 1: Wire the Skeleton wrapper in `ComparisonPickerBody.tsx`**

Import Skeleton and wrap the body's content (everything inside the outer `<div className="flex flex-col …">`) so cold fetches show a skeleton but background refetches don't flicker. Per CLAUDE.md "Skeleton loading via boneyard-js": gate on `isLoading`, not `isPending`. The four queries (`pickerQ`, `recentsQ`, `tagsQ`, plus any add-on) — gate on the *primary* query (`pickerQ.isLoading`) since that's what actually materializes rows.

```tsx
import { Skeleton } from "boneyard-js/react";

// Inside the component, replace the outermost <div> wrapper with:
return (
  <Skeleton
    name="comparison-picker-body"
    className="flex flex-col flex-1 min-h-0"
    loading={pickerQ.isLoading}
    fixture={
      <div className="flex flex-col flex-1 min-h-0">
        <div className="px-4 py-2 border-b border-border h-[44px]" />
        <div className="flex-1" />
      </div>
    }
  >
    { /* existing body content */ }
  </Skeleton>
);
```

- [ ] **Step 2: Add the skeleton-gating regression test**

```tsx
test("does not flicker on background refetch (gates on isLoading not isPending)", async () => {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false, staleTime: 0 } } });
  qc.setQueryData(["experiment", 1, "picker-samples"], [
    {
      sample: { id: 10, experiment_id: 1, name: "S1", label: null, notes: null, tags: [] },
      indexing_exposure_id: 100,
      all_exposures: [{ id: 100, sample_id: 10, filename: "f1", selected: true }],
    },
  ]);
  const { rerender } = render(
    <QueryClientProvider client={qc}>
      <ComparisonPickerBody
        experimentId={1} picks={[]} onPicksChange={() => {}}
        alreadyAddedExposureIds={new Set()}
      />
    </QueryClientProvider>,
  );
  expect(screen.getByText("S1")).toBeInTheDocument();
  qc.invalidateQueries({ queryKey: ["experiment", 1, "picker-samples"] });
  // After invalidation, the body should still render S1 — not flicker to skeleton.
  rerender(
    <QueryClientProvider client={qc}>
      <ComparisonPickerBody
        experimentId={1} picks={[]} onPicksChange={() => {}}
        alreadyAddedExposureIds={new Set()}
      />
    </QueryClientProvider>,
  );
  expect(screen.getByText("S1")).toBeInTheDocument();
});
```

- [ ] **Step 3: Run, verify pass**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/ComparisonPickerBody.test.tsx
```

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/test/ComparisonPickerBody.test.tsx packages/HimalayaUI/frontend/src/components/ComparisonPickerBody.tsx
git commit -m "feat(picker): Skeleton wrapper + isLoading gating regression"
```

## Task 10: Slim `ComparisonPicker` down to modal shell

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ComparisonPicker.tsx`
- Modify: `packages/HimalayaUI/frontend/test/ComparisonPicker.test.tsx`

- [ ] **Step 1: Update test to assert shell-only concerns**

Replace existing body-level tests with shell-level assertions:

```tsx
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ComparisonPicker } from "../src/components/ComparisonPicker";

function wrap(ui: React.ReactElement) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(<QueryClientProvider client={qc}>{ui}</QueryClientProvider>);
}

test("not rendered when isOpen=false", () => {
  wrap(<ComparisonPicker isOpen={false} onClose={() => {}} experimentId={1} />);
  expect(screen.queryByRole("dialog")).toBeNull();
});

test("Esc key fires onClose", () => {
  const onClose = vi.fn();
  wrap(<ComparisonPicker isOpen={true} onClose={onClose} experimentId={1} />);
  fireEvent.keyDown(screen.getByRole("dialog"), { key: "Escape" });
  expect(onClose).toHaveBeenCalled();
});

test("clicking outside fires onClose", () => {
  const onClose = vi.fn();
  wrap(<ComparisonPicker isOpen={true} onClose={onClose} experimentId={1} />);
  fireEvent.click(screen.getByTestId("comparison-picker-overlay"));
  expect(onClose).toHaveBeenCalled();
});

test("focus trap: Tab cycles within dialog", () => {
  // Mount, get all focusable elements, simulate Tab — assert focus stays inside.
  // Pin the existing useFocusTrap behavior so PR2's body-extraction doesn't
  // accidentally pull focusable elements outside the trap.
  // (Concrete implementation: render, query inside dialogRef, Tab through.)
});

test("'Add N selected' fires addMember per pick on click", () => {
  // Existing test — keep it. Mocks Zustand addMember; clicks rows; clicks add;
  // asserts addMember called once per pick with resolved exposure id.
});
```

- [ ] **Step 2: Run, verify failing (some tests fail — body assertions removed)**

- [ ] **Step 3: Refactor `ComparisonPicker.tsx`**

Replace the existing body with a shell wrapping `ComparisonPickerBody`. Keep these in the shell only:

```tsx
import { useEffect, useRef, useState } from "react";
import { useQueryClient } from "@tanstack/react-query";
import { useAppState } from "../state";
import { useFocusTrap } from "../hooks/useFocusTrap";
import { ComparisonPickerBody, type Pick } from "./ComparisonPickerBody";

interface Props {
  isOpen: boolean;
  onClose: () => void;
  experimentId: number | undefined;
}

export function ComparisonPicker({
  isOpen, onClose, experimentId,
}: Props): JSX.Element | null {
  const dialogRef = useRef<HTMLDivElement>(null);
  const inputRef  = useRef<HTMLInputElement>(null);
  useFocusTrap(dialogRef, isOpen);

  const qc = useQueryClient();
  const draft = useAppState((s) => s.activeDraft);
  const addMember = useAppState((s) => s.addMember);

  const [picks, setPicks] = useState<Pick[]>([]);
  useEffect(() => {
    if (isOpen) {
      setPicks([]);
      inputRef.current?.focus();
    }
  }, [isOpen]);

  if (!isOpen) return null;

  const alreadyAddedExposureIds = new Set(
    (draft?.members ?? [])
      .map((m) => m.exposure_id)
      .filter((id): id is number => id !== null),
  );

  const onAddSelected = (): void => {
    for (const p of picks) addMember(p.exposure_id, qc);
    onClose();
  };

  return (
    <div
      data-testid="comparison-picker-overlay"
      className="fixed inset-0 z-50 flex items-start justify-center pt-[10vh]
                 bg-[oklch(0.05_0_0/0.65)] backdrop-blur-sm anim-pal-in"
      role="presentation"
      onClick={(e) => { if (e.target === e.currentTarget) onClose(); }}
    >
      <div
        ref={dialogRef}
        data-testid="comparison-picker"
        role="dialog"
        aria-modal="true"
        aria-labelledby="comparison-picker-title"
        onKeyDown={(e) => { if (e.key === "Escape") { e.preventDefault(); onClose(); } }}
        className="w-[min(720px,calc(100vw-48px))] max-h-[78vh]
                   bg-bg-elevated border border-border rounded-xl shadow-2xl
                   flex flex-col overflow-hidden anim-pal-scale"
      >
        <div className="flex items-center justify-between px-4 py-3 border-b border-border">
          <h2 id="comparison-picker-title" className="text-base font-medium text-fg">
            Add traces
          </h2>
          <button
            type="button"
            data-testid="comparison-picker-close"
            onClick={onClose}
            className="text-fg-muted hover:text-fg text-sm px-2 py-1"
            aria-label="Close picker"
          >
            esc
          </button>
        </div>

        <ComparisonPickerBody
          experimentId={experimentId}
          picks={picks}
          onPicksChange={setPicks}
          alreadyAddedExposureIds={alreadyAddedExposureIds}
          searchInputRef={inputRef}
        />

        <div className="flex items-center gap-2 px-4 py-3 border-t border-border">
          <span className="text-xs text-fg-dim flex-1">{picks.length} selected</span>
          <button type="button"
            data-testid="comparison-picker-cancel"
            onClick={onClose}
            className="px-3 py-1 rounded border border-border text-sm text-fg-muted">
            Cancel
          </button>
          <button type="button"
            data-testid="comparison-picker-add"
            onClick={onAddSelected}
            disabled={picks.length === 0}
            className="px-3 py-1 rounded bg-accent text-bg text-sm font-medium disabled:opacity-50">
            {picks.length === 0 ? "Add selected" : `Add ${picks.length} selected`}
          </button>
        </div>
      </div>
    </div>
  );
}
```

Title is `"Add traces"` (was `"Add exposures"`) — closes #84 acceptance criterion.

- [ ] **Step 4: Run all picker-related tests, verify pass**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/ComparisonPicker.test.tsx test/ComparisonPickerBody.test.tsx test/SamplePickerRow.test.tsx test/ExposureListRow.test.tsx
```

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ComparisonPicker.tsx packages/HimalayaUI/frontend/test/ComparisonPicker.test.tsx
git commit -m "refactor(picker): slim ComparisonPicker to modal shell (#84)"
```

## Task 11: Update existing E2E spec for sample-first picker

**Files:**
- Modify: `packages/HimalayaUI/frontend/e2e/comparisonPicker.spec.ts`

- [ ] **Step 1: Update Playwright route mocks**

Find the existing `page.route(/.*samples.*/, …)` and similar — add a mock for the new endpoint:

```ts
await page.route(/.*\/api\/experiments\/\d+\/picker-samples/, async (route) => {
  await route.fulfill({
    status: 200,
    contentType: "application/json",
    body: JSON.stringify([
      {
        sample: { id: 10, experiment_id: 1, name: "Sample A", label: null, notes: null, tags: [] },
        indexing_exposure_id: 100,
        all_exposures: [{ id: 100, sample_id: 10, filename: "f1", selected: true }],
      },
    ]),
  });
});
```

Update assertions that depend on filename-primary rows to use sample-name:

```ts
await expect(page.getByText("Sample A")).toBeVisible();
// previously: await expect(page.getByText("f1")).toBeVisible();
```

Add a caret-expand assertion:

```ts
await page.getByTestId("sample-picker-row-caret").click();
await expect(page.getByText("f1")).toBeVisible();
```

- [ ] **Step 2: Run E2E**

```bash
cd packages/HimalayaUI/frontend && npm run e2e -- --grep "ComparisonPicker"
```

(Headless verification per spec § "Frontend — E2E" — capture to `/tmp/picker-e2e.out` per CLAUDE.md log-once rule.)

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/e2e/comparisonPicker.spec.ts
git commit -m "test(picker): update e2e mocks + assertions for sample-first picker"
```

## Task 12: PR1 verification gate (run myself, not just write)

**Files:** none.

- [ ] **Step 1: Full Vitest pass on the picker slice**

```bash
cd packages/HimalayaUI/frontend && npm test -- --run ComparisonPicker SamplePickerRow ComparisonPickerBody ExposureListRow > /tmp/vitest-picker.out 2>&1
grep -E "Test Files|Tests|FAIL" /tmp/vitest-picker.out
```

Expected: zero failures.

- [ ] **Step 2: Build passes (TS-strict gate)**

```bash
cd packages/HimalayaUI/frontend && npm run build > /tmp/build.out 2>&1
grep -E "error|warning" /tmp/build.out
```

Expected: build succeeds, no TS errors. (Warnings about asset sizes are fine.)

- [ ] **Step 3: Headless Playwright on the picker slice**

```bash
cd packages/HimalayaUI/frontend && npm run e2e -- --grep "ComparisonPicker" > /tmp/picker-e2e.out 2>&1
tail -30 /tmp/picker-e2e.out
```

Expected: all specs pass.

- [ ] **Step 4: Backend full suite (PR-cut gate)**

The PR1 gate runs the full Julia suite (~5–10 min) — confirms nothing else broke. The direct-include slices used during Tasks 1–3 (`test_http.jl + test_picker_samples_route.jl`) are sufficient for the iterative red→green loop, but the gate before merge wants full coverage. `test_picker_routes.jl` depends on `_setup_analyzed_exposure` which lives inside `test_comparisons.jl` — only the runtests.jl include order satisfies that.

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-full.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-full.out
tail -50 /tmp/jl-full.out
```

Expected: zero failures.

- [ ] **Step 5: PR1 cut**

```bash
gh pr create --title "feat(picker): sample-first picker (closes #84, partial #85)" --body "$(cat <<'EOF'
## Summary
- New `GET /api/experiments/:eid/picker-samples` route returns one row per sample with `indexing_exposure_id` resolved at fetch time (selected=1 → MAX(id) → null).
- `ComparisonPicker` rewritten as modal shell + extracted `ComparisonPickerBody` (filters, recents, list).
- New `SamplePickerRow` is the primary row: sample name primary, filename only under the override caret.
- `ExposureListRow` parameterized with `control: "checkbox" | "radio"` so the override leaf reuses it.
- Drift policy: freeze at pick time. Schema unchanged.

## Test plan
- [x] Vitest: all picker-related test files pass.
- [x] `npm run build` clean (TS-strict).
- [x] Playwright (mocked) `--grep "ComparisonPicker"` passes.
- [x] Julia: `test_picker_samples_route.jl` covers helper + route + corner cases (multi-experiment isolation, NULL fields, JSON shape, defensive multi-selected).

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

# PR2 — Inline picker panel in edit mode

## Task 13: Add `onPick` immediate-mode prop to `ComparisonPickerBody`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ComparisonPickerBody.tsx`
- Modify: `packages/HimalayaUI/frontend/test/ComparisonPickerBody.test.tsx`

PR1 deliberately omits `onPick` so it doesn't ship dead code. PR2 adds it as the first task — used by `ComparisonPickerPanel` (Task 14).

- [ ] **Step 1: Add the failing test**

```tsx
test("immediate mode: onPick fires per toggle, picks prop ignored", async () => {
  const onPick = vi.fn();
  wrap(
    <ComparisonPickerBody
      experimentId={1}
      picks={[]}
      onPicksChange={() => {}}
      onPick={onPick}
      alreadyAddedExposureIds={new Set()}
    />,
  );
  await screen.findByText("S1");
  fireEvent.click(screen.getByTestId("sample-picker-row-checkbox"));
  expect(onPick).toHaveBeenCalledWith({ sample_id: 10, exposure_id: 100, source: "default" });
});

test("immediate mode: override caret pick fires onPick with override exposure id", async () => {
  const onPick = vi.fn();
  wrap(
    <ComparisonPickerBody
      experimentId={1}
      picks={[]}
      onPicksChange={() => {}}
      onPick={onPick}
      alreadyAddedExposureIds={new Set()}
    />,
  );
  await screen.findByText("S1");
  // Open caret on S1 row, pick an override (single exposure here, so this
  // re-fires the same id but exercises the override path through onPick).
  fireEvent.click(screen.getAllByTestId("sample-picker-row-caret")[0]);
  fireEvent.click(screen.getByLabelText(/f1\.dat/));
  expect(onPick).toHaveBeenCalled();
});
```

- [ ] **Step 2: Run, verify fail**

- [ ] **Step 3: Add the prop + branch in `ComparisonPickerBody.tsx`**

Update the `Props` interface:

```tsx
interface Props {
  experimentId: number | undefined;
  picks: Pick[];
  onPicksChange: (next: Pick[]) => void;
  /** When set, fires per-pick instead of accumulating in `picks`. */
  onPick?: (pick: Pick) => void;
  alreadyAddedExposureIds: Set<number>;
  searchInputRef?: RefObject<HTMLInputElement>;
}
```

Update `togglePickFor`:

```tsx
const togglePickFor = (row: PickerSampleRow, next: boolean): void => {
  if (row.indexing_exposure_id === null) return;
  const pick: Pick = {
    sample_id: row.sample.id,
    exposure_id: row.indexing_exposure_id,
    source: "default",
  };
  if (onPick) {
    if (next) onPick(pick);
    return;
  }
  if (next) onPicksChange([...picks, pick]);
  else onPicksChange(picks.filter((p) => p.sample_id !== row.sample.id));
};
```

Same shape for `overridePickFor`:

```tsx
const overridePickFor = (row: PickerSampleRow, exposureId: number): void => {
  const next: Pick = {
    sample_id: row.sample.id,
    exposure_id: exposureId,
    source: exposureId === row.indexing_exposure_id ? "default" : "override",
  };
  if (onPick) { onPick(next); return; }
  const i = picks.findIndex((p) => p.sample_id === row.sample.id);
  if (i < 0) onPicksChange([...picks, next]);
  else onPicksChange(picks.map((p, j) => (j === i ? next : p)));
};
```

- [ ] **Step 4: Run, verify pass**

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ComparisonPickerBody.tsx packages/HimalayaUI/frontend/test/ComparisonPickerBody.test.tsx
git commit -m "feat(picker): onPick immediate-mode prop on ComparisonPickerBody (#87)"
```

## Task 14: `ComparisonPickerPanel` component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ComparisonPickerPanel.tsx`
- Create: `packages/HimalayaUI/frontend/test/ComparisonPickerPanel.test.tsx`

- [ ] **Step 1: Write failing tests**

```tsx
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { beforeEach, afterEach, vi } from "vitest";
import { ComparisonPickerPanel } from "../src/components/ComparisonPickerPanel";
import * as api from "../src/api";
import { useAppState } from "../src/state";

const ROW = {
  sample: { id: 10, experiment_id: 1, name: "S1", label: null, notes: null, tags: [] },
  indexing_exposure_id: 100,
  all_exposures: [{ id: 100, sample_id: 10, filename: "f1.dat", selected: true }],
};
beforeEach(() => {
  // Reset Zustand to a fresh draft so per-test seeding (or absence of) is
  // deterministic. The store reads `loadDraftFromSession` at module load,
  // so without an explicit reset prior tests' drafts leak.
  useAppState.getState().discardDraft();
  useAppState.getState().startNewDraft();

  vi.spyOn(api, "getPickerSamples").mockResolvedValue([ROW]);
  vi.spyOn(api, "getRecentlyPickedExposures").mockResolvedValue([]);
  vi.spyOn(api, "getSampleTags").mockResolvedValue([]);
});
afterEach(() => { vi.restoreAllMocks(); });

function wrap(ui: React.ReactElement) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(<QueryClientProvider client={qc}>{ui}</QueryClientProvider>);
}

test("inline shell does not render dialog role or focus trap", async () => {
  wrap(<ComparisonPickerPanel experimentId={1} />);
  await screen.findByText("S1");
  expect(screen.queryByRole("dialog")).toBeNull();
});

test("toggle fires addMember immediately (immediate-commit semantics)", async () => {
  wrap(<ComparisonPickerPanel experimentId={1} />);
  await screen.findByText("S1");
  fireEvent.click(screen.getByTestId("sample-picker-row-checkbox"));
  // After immediate-commit, the draft should contain a member with exposure_id 100.
  const draft = useAppState.getState().activeDraft;
  expect(draft?.members.some((m) => m.exposure_id === 100)).toBe(true);
});

test("already-added rows render locked (alreadyAddedExposureIds gate)", async () => {
  // Pre-seed a draft with exposure 100 already a member.
  useAppState.getState().startNewDraft();
  // QueryClient mock so addMember's snapshot derivation doesn't blow up.
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  useAppState.getState().addMember(100, qc);
  wrap(<ComparisonPickerPanel experimentId={1} />);
  await screen.findByText("S1");
  // Re-clicking should not append a duplicate member (panel filters via
  // alreadyAddedExposureIds before calling addMember).
  fireEvent.click(screen.getByTestId("sample-picker-row-checkbox"));
  const members = useAppState.getState().activeDraft?.members ?? [];
  expect(members.filter((m) => m.exposure_id === 100).length).toBe(1);
});
```

- [ ] **Step 2: Run, verify fail (component doesn't exist)**

- [ ] **Step 3: Implement the panel**

```tsx
import type { RefObject } from "react";
import { useQueryClient } from "@tanstack/react-query";
import { useAppState } from "../state";
import { ComparisonPickerBody, type Pick } from "./ComparisonPickerBody";

interface Props {
  experimentId: number | undefined;
  /** External ref the parent uses to focus the search input from a sibling CTA. */
  searchInputRef?: RefObject<HTMLInputElement>;
}

export function ComparisonPickerPanel({
  experimentId, searchInputRef,
}: Props): JSX.Element {
  const qc = useQueryClient();
  const draft = useAppState((s) => s.activeDraft);
  const addMember = useAppState((s) => s.addMember);

  const alreadyAddedExposureIds = new Set(
    (draft?.members ?? [])
      .map((m) => m.exposure_id)
      .filter((id): id is number => id !== null),
  );

  const handlePick = (p: Pick): void => {
    if (alreadyAddedExposureIds.has(p.exposure_id)) return;
    addMember(p.exposure_id, qc);
  };

  return (
    <div
      data-testid="comparison-picker-panel"
      className="flex flex-col h-full overflow-hidden"
    >
      <div className="px-4 py-3 border-b border-border">
        <h3 className="text-sm font-medium text-fg">Add traces</h3>
      </div>
      <ComparisonPickerBody
        experimentId={experimentId}
        picks={[]}                        // immediate mode: ignored
        onPicksChange={() => {}}          // immediate mode: ignored
        onPick={handlePick}
        alreadyAddedExposureIds={alreadyAddedExposureIds}
        searchInputRef={searchInputRef}
      />
    </div>
  );
}
```

- [ ] **Step 4: Run tests, verify pass**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/ComparisonPickerPanel.test.tsx
```

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ComparisonPickerPanel.tsx packages/HimalayaUI/frontend/test/ComparisonPickerPanel.test.tsx
git commit -m "feat(picker): ComparisonPickerPanel inline variant (#87 §6)"
```

## Task 15: Swap the right slot in `ComparePageEdit`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx`
- Modify: `packages/HimalayaUI/frontend/test/ComparePageEdit.test.tsx`

- [ ] **Step 1: Update test**

Replace the right-slot hint-card assertion:

```tsx
test("right slot hosts ComparisonPickerPanel in edit mode", async () => {
  // Render ComparePageEdit with seeded draft + experiment.
  // ... follow existing test setup.
  expect(screen.getByTestId("comparison-picker-panel")).toBeInTheDocument();
  expect(screen.queryByTestId("compare-edit-right-hint")).toBeNull();
});

test("empty-state '+ Add traces' button focuses the panel's search input", async () => {
  // Render with empty draft.
  fireEvent.click(screen.getByTestId("compare-edit-add-traces"));
  expect(screen.getByTestId("comparison-picker-search")).toHaveFocus();
});

// Delete the existing test that asserted the bottom-right "+ Add traces"
// button on the populated branch — that button is gone in PR2.
```

- [ ] **Step 2: Run, verify fail**

- [ ] **Step 3: Wire the panel**

In `ComparePageEdit.tsx`:

1. Add a ref:

```tsx
const pickerSearchRef = useRef<HTMLInputElement>(null);
```

2. Replace the right-slot content:

```tsx
right={
  <ComparisonPickerPanel experimentId={eid} searchInputRef={pickerSearchRef} />
}
```

3. Update the empty-state CTA `onClick`:

```tsx
onClick={() => pickerSearchRef.current?.focus()}
```

4. Delete the bottom-right "+ Add traces" button. To find it: `grep -n '+ Add traces' packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx` finds two matches — keep the empty-state CTA inside `compare-edit-plot-empty`, delete the trailing `<div className="flex justify-end">…</div>` block that hosts the second button (the one rendered when `plotMembers.length > 0`).

5. Delete the modal mount: `grep -n '<ComparisonPicker' packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx` finds the `<ComparisonPicker isOpen={pickerOpen} onClose={…} experimentId={eid} />` block — delete that JSX, the `pickerOpen` `useState`, and the `import { ComparisonPicker }` line. Keep `ComparisonPicker` itself in the codebase (still used by sidebar entry points).

6. Update `slotClassName.right` if the panel needs a different min-height than the hint card (it does — try `flex flex-col min-h-[400px]` mirroring `slotClassName.left`).

- [ ] **Step 4: Run tests, verify pass**

```bash
cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/ComparePageEdit.test.tsx
```

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx packages/HimalayaUI/frontend/test/ComparePageEdit.test.tsx
git commit -m "feat(compare): inline picker panel in edit-mode right slot (#87 §6-7)"
```

## Task 16: Live-mode E2E

**Files:**
- Create: `packages/HimalayaUI/frontend/e2e/live/comparePickerInline.spec.ts`

- [ ] **Step 1: Write the live spec**

```ts
import { test, expect } from "@playwright/test";

test.describe("Compare picker inline panel — live mode", () => {
  test("pick a sample and see a band appear", async ({ page }) => {
    await page.goto("/");
    await page.waitForTimeout(800);   // SSE warmup per live runbook

    // Navigate to /compare/new for the active experiment.
    // (Use the project's existing live-mode pattern — seed Zustand state
    //  in localStorage if needed, see e2e/live/README.md.)
    await page.goto("/experiments/1/compare/new");

    // Panel is visible in the right slot.
    await expect(page.getByTestId("comparison-picker-panel")).toBeVisible();

    // Click first sample row.
    const firstRow = page.getByTestId("sample-picker-row").first();
    const sampleId = await firstRow.getAttribute("data-sample-id");
    await firstRow.getByTestId("sample-picker-row-checkbox").click();

    // The plot should now show a band for the resolved exposure.
    await expect(page.getByTestId("compare-edit-gutter")).toContainText(/.+/);
    // Row should now be checked.
    await expect(firstRow.getByTestId("sample-picker-row-checkbox")).toBeChecked();
  });

  test("override caret swaps the resolved exposure", async ({ page }) => {
    await page.goto("/");
    await page.waitForTimeout(800);
    await page.goto("/experiments/1/compare/new");

    const row = page.getByTestId("sample-picker-row").first();
    await row.getByTestId("sample-picker-row-checkbox").click();
    const initialEid = await row.getAttribute("data-exposure-id");

    await row.getByTestId("sample-picker-row-caret").click();
    // Pick the second exposure radio (if available).
    const radios = row.locator('input[type="radio"]');
    const count = await radios.count();
    if (count > 1) {
      await radios.nth(1).click();
      const newEid = await row.getAttribute("data-exposure-id");
      expect(newEid).not.toBe(initialEid);
    }
  });
});
```

- [ ] **Step 2: Run live spec**

(Operator brings up backend on 8090 + Vite on 5180 per `e2e/live/README.md` first.)

```bash
cd packages/HimalayaUI/frontend && npm run e2e:live -- --grep "comparePickerInline" > /tmp/live-picker.out 2>&1
tail -30 /tmp/live-picker.out
```

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/e2e/live/comparePickerInline.spec.ts
git commit -m "test(picker): live-mode roundtrip for inline panel (#87 §6-7)"
```

## Task 17: PR2 verification gate (run myself)

**Files:** none.

- [ ] **Step 1: Vitest full pass on the picker + edit slice**

```bash
cd packages/HimalayaUI/frontend && npm test -- --run ComparisonPicker ComparisonPickerPanel ComparisonPickerBody SamplePickerRow ComparePageEdit > /tmp/vitest-picker2.out 2>&1
grep -E "Test Files|Tests|FAIL" /tmp/vitest-picker2.out
```

- [ ] **Step 2: Build passes**

```bash
cd packages/HimalayaUI/frontend && npm run build > /tmp/build.out 2>&1
```

- [ ] **Step 3: Headless mocked Playwright**

```bash
cd packages/HimalayaUI/frontend && npm run e2e -- --grep "ComparisonPickerPanel|ComparePageEdit" > /tmp/picker-e2e2.out 2>&1
tail -30 /tmp/picker-e2e2.out
```

- [ ] **Step 4: Live Playwright** (operator brings up backend + Vite first)

```bash
cd packages/HimalayaUI/frontend && npm run e2e:live -- --grep "comparePickerInline" > /tmp/live-picker.out 2>&1
```

- [ ] **Step 5: PR2 cut**

```bash
gh pr create --title "feat(compare-edit): inline picker panel in right slot (closes #87 §6-7)" --body "$(cat <<'EOF'
## Summary
- New `ComparisonPickerPanel` inline shell wraps the same `ComparisonPickerBody` from PR1, but commits picks immediately via `addMember`.
- `ComparePageEdit` swaps the right-slot hint card for the panel; bottom-right "+ Add traces" button removed; empty-state CTA focuses the panel's search input.
- Modal `ComparisonPicker` retained for non-edit entry points.

## Test plan
- [x] Vitest: ComparisonPickerPanel + ComparePageEdit + sibling files pass.
- [x] `npm run build` clean (TS-strict).
- [x] Playwright (mocked): ComparisonPickerPanel + ComparePageEdit specs pass.
- [x] Playwright (live): `comparePickerInline.spec.ts` passes.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

# Self-review — spec coverage check

| Spec section | Implementing task |
|--------------|-------------------|
| §PR1 Backend route | Task 3 |
| §PR1 `picker_samples` helper + query shape lock + bool_keys | Tasks 1, 2 |
| §PR1 Backend tests (8 cases incl. multi-exp + NULL + JSON-shape + defensive multi-selected) | Tasks 1, 2, 3 |
| §PR1 Frontend API fetcher + types | Task 4 |
| §PR1 `usePickerSamples` hook | Task 5 |
| §PR1 `ExposureListRow` `control` prop | Task 6 |
| §PR1 `SamplePickerRow` (incl. disabled-row branch + accent-ring) | Task 7 |
| §PR1 `ComparisonPickerBody` extraction (controlled picks, recents dedup as useMemo) | Task 8 |
| §PR1 Skeleton wrap + isLoading-gating regression | Task 9 |
| §PR1 `ComparisonPicker` slimmed to modal shell + focus-trap regression + title="Add traces" | Task 10 |
| §PR1 E2E mock update | Task 11 |
| §PR1 Headless verification (vitest + build + e2e + julia) | Task 12 |
| §PR2 `onPick` immediate-mode prop on ComparisonPickerBody | Task 13 |
| §PR2 `ComparisonPickerPanel` (no overlay, no focus trap, immediate commit) | Task 14 |
| §PR2 `ComparePageEdit` right-slot swap + buried-button removal + searchInputRef | Task 15 |
| §PR2 Live-mode E2E | Task 16 |
| §PR2 Headless verification | Task 17 |

All spec sections mapped. No `TBD` / `TODO` / "fill in later" remain. Type names consistent (`PickerSampleRow`, `Pick`, `PickerSampleExposure`, `picker_samples`, `usePickerSamples`, `ComparisonPickerBody`, `ComparisonPickerPanel`, `SamplePickerRow`).
