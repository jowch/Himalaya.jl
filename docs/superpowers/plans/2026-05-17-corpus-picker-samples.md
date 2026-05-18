# Corpus-wide `GET /api/picker-samples` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a corpus-wide `GET /api/picker-samples` endpoint returning the sample-picker projection across every experiment, as a sibling to the experiment-scoped `GET /api/experiments/{eid}/picker-samples`.

**Architecture:** The existing `picker_samples(db, experiment_id)` helper in `comparisons.jl` is refactored so its experiment-agnostic body (exposure/tag bulk fetch, grouping, indexing-exposure resolution, dict assembly) moves into a private `_picker_samples_projection(db, samples)`. Two thin public methods then dispatch on arity: the 2-arg method keeps the experiment-scoped `samples` query; a new 1-arg method runs an unscoped query. A new route in `routes_picker.jl` calls the 1-arg method.

**Tech Stack:** Julia, Oxygen.jl (`@get`), SQLite.jl / DBInterface, JSON3, stdlib `Test`.

**Spec:** `docs/superpowers/specs/2026-05-17-corpus-picker-samples-design.md`

---

## File Structure

- **Modify** `packages/HimalayaUI/src/comparisons.jl` (function at lines ~490–575) — split `picker_samples` into a shared private helper plus two dispatching public methods.
- **Modify** `packages/HimalayaUI/src/routes_picker.jl` — add the `GET /api/picker-samples` route inside `register_picker_routes!()`.
- **Modify** `packages/HimalayaUI/test/test_picker_samples_route.jl` — append a corpus helper testset and a corpus route testset.

No change to `server.jl` — `register_picker_routes!()` is already wired in.

## Running tests

The full HimalayaUI Julia suite is slow (5–10 min). For fast TDD iteration, run only the two test files involved — `test_http.jl` defines `with_test_server` and must be included first. Run from the repo root:

```bash
julia --project=packages/HimalayaUI -e '
using Test
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_picker_samples_route.jl")
' 2>&1 | tail -40
```

The final verification (Task 3) runs the full suite once.

---

## Task 1: Corpus `picker_samples(db)` method via dispatch refactor

**Files:**
- Modify: `packages/HimalayaUI/src/comparisons.jl` (lines ~490–575)
- Test: `packages/HimalayaUI/test/test_picker_samples_route.jl`

- [ ] **Step 1: Write the failing test**

Append this testset to `packages/HimalayaUI/test/test_picker_samples_route.jl`, immediately after the closing `end` of the existing `@testset "picker_samples helper"` block (before the `@testset "GET /api/experiments/:eid/picker-samples"` block):

```julia
@testset "picker_samples corpus (no experiment_id)" begin
    @testset "empty corpus → empty list" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "h.db"))
            @test picker_samples(db) == Dict{Symbol, Any}[]
        end
    end

    @testset "spans all experiments, ordered by (experiment_id, id)" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "h.db"))
            DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'A', '/tmp/a', '/tmp/a/data', '/tmp/a/analysis')")
            DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (2, 'B', '/tmp/b', '/tmp/b/data', '/tmp/b/analysis')")
            # Insert out of order to prove the ORDER BY, not insertion order.
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (30, 2, 'B-late')")
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'A1')")
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (20, 2, 'B-early')")
            rows = picker_samples(db)
            @test length(rows) == 3
            @test [r[:sample][:id] for r in rows] == [10, 20, 30]            # ORDER BY experiment_id, id
            @test [r[:sample][:experiment_id] for r in rows] == [1, 2, 2]
        end
    end

    @testset "cross-experiment exposure resolution stays correct" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "h.db"))
            DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'A', '/tmp/a', '/tmp/a/data', '/tmp/a/analysis')")
            DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (2, 'B', '/tmp/b', '/tmp/b/data', '/tmp/b/analysis')")
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'A1')")
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (20, 2, 'B1')")
            DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'a', 0)")
            DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (999, 20, 'b', 0)")  # bigger global id, other experiment
            rows = picker_samples(db)
            a = first(filter(r -> r[:sample][:id] == 10, rows))
            @test a[:indexing_exposure_id] == 100                            # not 999 — exposures grouped by sample_id
            @test [e[:id] for e in a[:all_exposures]] == [100]
        end
    end
end
```

No new imports are needed — `test_picker_samples_route.jl` already does `using HimalayaUI: open_db, picker_samples`, and the corpus method shares the `picker_samples` binding.

- [ ] **Step 2: Run the test to verify it fails**

Run:

```bash
julia --project=packages/HimalayaUI -e '
using Test
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_picker_samples_route.jl")
' 2>&1 | tail -40
```

Expected: FAIL inside `picker_samples corpus` — `MethodError: no method matching picker_samples(::SQLite.DB)` (only the 2-arg method exists).

- [ ] **Step 3: Implement the dispatch refactor**

In `packages/HimalayaUI/src/comparisons.jl`, replace the entire docstring + function block (the `"""`-docstring beginning `picker_samples(db, experiment_id)` through the `end` that closes `function picker_samples`, lines ~490–575) with the following:

```julia
"""
    picker_samples(db, experiment_id) -> Vector{Dict{Symbol, Any}}
    picker_samples(db)                -> Vector{Dict{Symbol, Any}}

Picker primary list per spec §"PR1 — sample-first picker → Backend".

The two-argument method is experiment-scoped. The one-argument method returns
the projection for every sample in the corpus, ordered by `(experiment_id, id)`
— used by the corpus-wide `GET /api/picker-samples` route.

For each sample, returns:
  :sample              => sample row (Symbol-keyed Dict, with :tags vector)
  :indexing_exposure_id => Int or nothing — `selected = 1` else MAX(id) else nothing
  :all_exposures       => Vector of {:id, :filename, :selected::Bool}, ORDER BY id ASC

Three bulk queries (no JOIN'd Cartesian flatten, no per-sample N+1):
(1) samples (scoped to one experiment, or the whole corpus),
(2) `WHERE sample_id IN (...)` for exposures,
(3) `WHERE sample_id IN (...)` for sample_tags.
Empty result ⇒ [].
"""
function picker_samples(db::SQLite.DB, experiment_id::Integer)::Vector{Dict{Symbol, Any}}
    # Explicit column list (PR #96 review): keep the picker JSON shape
    # deliberate so a future column added to `samples` doesn't auto-leak
    # into the picker payload.
    samples = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, experiment_id, name, display_name, notes
         FROM samples WHERE experiment_id = ? ORDER BY id",
        [Int(experiment_id)]))
    _picker_samples_projection(db, samples)
end

function picker_samples(db::SQLite.DB)::Vector{Dict{Symbol, Any}}
    # Corpus-wide: every sample, ORDER BY (experiment_id, id) for stable,
    # experiment-grouped output. `experiment_id` is in the column list so a
    # consumer can group client-side. Same explicit column list as the
    # scoped method — the JSON shape must not diverge.
    samples = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, experiment_id, name, display_name, notes
         FROM samples ORDER BY experiment_id, id"))
    _picker_samples_projection(db, samples)
end

# Shared body: everything downstream of the `samples` fetch is driven purely
# by sample_ids and is experiment-agnostic. `samples` is the row table from
# either picker_samples method above.
function _picker_samples_projection(db::SQLite.DB, samples)::Vector{Dict{Symbol, Any}}
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

This preserves the 2-arg method's observable behavior exactly — only the post-`samples` body moved into `_picker_samples_projection`, called identically.

- [ ] **Step 4: Run the tests to verify they pass**

Run:

```bash
julia --project=packages/HimalayaUI -e '
using Test
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_picker_samples_route.jl")
' 2>&1 | tail -40
```

Expected: PASS — the new `picker_samples corpus` testset is green, and the pre-existing `picker_samples helper` testsets (the regression guard for the refactor) remain green.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/comparisons.jl packages/HimalayaUI/test/test_picker_samples_route.jl
git commit -m "feat: add corpus-wide picker_samples(db) via dispatch refactor (#158)"
```

---

## Task 2: `GET /api/picker-samples` route

**Files:**
- Modify: `packages/HimalayaUI/src/routes_picker.jl`
- Test: `packages/HimalayaUI/test/test_picker_samples_route.jl`

- [ ] **Step 1: Write the failing test**

Append this testset to the end of `packages/HimalayaUI/test/test_picker_samples_route.jl` (after the existing `@testset "GET /api/experiments/:eid/picker-samples"` block):

```julia
@testset "GET /api/picker-samples (corpus)" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "h.db"))
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'A', '/tmp/a', '/tmp/a/data', '/tmp/a/analysis')")
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (2, 'B', '/tmp/b', '/tmp/b/data', '/tmp/b/analysis')")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'A1')")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (20, 2, 'B1')")  # zero-exposure
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'f1', 1)")

        with_test_server(db) do port, base
            r = HTTP.get("$base/api/picker-samples")
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test length(body) == 2                          # samples from both experiments
            @test [b.sample.id for b in body] == [10, 20]     # ORDER BY experiment_id, id

            row10 = first(filter(b -> b.sample.id == 10, collect(body)))
            @test row10.indexing_exposure_id == 100
            @test row10.all_exposures[1].selected === true    # JSON-shape: bool, not 1

            row20 = first(filter(b -> b.sample.id == 20, collect(body)))
            @test row20.indexing_exposure_id === nothing      # zero-exposure → null
            @test haskey(row20, :indexing_exposure_id)        # null vs absent key
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:

```bash
julia --project=packages/HimalayaUI -e '
using Test
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_picker_samples_route.jl")
' 2>&1 | tail -40
```

Expected: FAIL inside `GET /api/picker-samples (corpus)` — the route is unregistered, so `HTTP.get` raises `HTTP.Exceptions.StatusError` for a 404 response.

- [ ] **Step 3: Add the route**

In `packages/HimalayaUI/src/routes_picker.jl`, inside `register_picker_routes!()`, add this block immediately after the existing `@get "/api/experiments/{eid}/picker-samples"` block (before the closing `end` of the function):

```julia
    @get "/api/picker-samples" function(req::HTTP.Request)
        # Corpus-wide sibling of the experiment-scoped picker-samples route.
        # Read-only — no with_idempotency, no SSE, no event-log row, matching
        # the other picker routes.
        db = current_db()
        rows = picker_samples(db)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(rows))
    end
```

No change to `server.jl` is needed — `register_picker_routes!()` is already called there.

- [ ] **Step 4: Run the tests to verify they pass**

Run:

```bash
julia --project=packages/HimalayaUI -e '
using Test
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_picker_samples_route.jl")
' 2>&1 | tail -40
```

Expected: PASS — `GET /api/picker-samples (corpus)` is green, and the pre-existing `GET /api/experiments/:eid/picker-samples` route testset is unchanged and green.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_picker.jl packages/HimalayaUI/test/test_picker_samples_route.jl
git commit -m "feat: add corpus-wide GET /api/picker-samples route (#158)"
```

---

## Task 3: Full-suite verification

**Files:** none — verification only.

- [ ] **Step 1: Run the full HimalayaUI Julia suite**

Run from the repo root:

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|did not pass|fail|FAIL" /tmp/jl-test.out
tail -50 /tmp/jl-test.out
```

Expected: every testset passes — in particular `test_picker_samples_route.jl` (both new and pre-existing testsets) and any other suite that touches `comparisons.jl` (e.g. `test_comparisons.jl`, comparison route tests). The suite takes 5–10 min; run it once.

- [ ] **Step 2: Confirm the frontend is unaffected**

No frontend file was touched and `picker-samples` has no frontend consumer. No frontend build or test run is required. Confirm with:

```bash
git diff --name-only main...HEAD
```

Expected: only `comparisons.jl`, `routes_picker.jl`, `test_picker_samples_route.jl`, and the spec/plan docs — no file under `frontend/`.

---

## Self-Review Notes

- **Spec coverage:** corpus route (Task 2), helper refactor (Task 1), Julia tests (Tasks 1 & 2), experiment route unchanged (regression-guarded in Steps 4 + Task 3), append-only / no frontend change (Task 3 Step 2). All spec sections map to a task.
- **Type consistency:** `picker_samples` and `_picker_samples_projection` signatures are consistent across Task 1 and their use in Task 2; the route returns the same `Vector{Dict{Symbol, Any}}` shape as the experiment-scoped route.
- **Sequencing:** the spec notes `routes_picker.jl` is shared with #157 — land #157 first to avoid the append conflict in `register_picker_routes!()`.
