# Corpus-wide `GET /api/samples` Route Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a corpus-wide `GET /api/samples` endpoint that returns every sample across all experiments in one projection — each sample carrying its `tags` and a `q_units` value sourced from the owning experiment's config — with an optional `?experiment_id=` filter.

**Architecture:** A new `@get "/api/samples"` route added to the existing `register_samples_routes!()` function in `routes_samples.jl`. It runs exactly three SQL queries regardless of corpus size (samples, experiment configs, batched tags), then joins them in memory — deliberately avoiding the sibling route's per-sample N+1 tag query. The `q_units` TOML-parse-with-fallback logic, currently inline in `routes_experiments.jl`, is first extracted to a shared `_q_units_from_config` helper so the new route reuses it rather than copying it.

**Tech Stack:** Julia, Oxygen.jl (routing), HTTP.jl, SQLite.jl / DBInterface, JSON3, stdlib TOML, stdlib Test.

---

## Background for the implementing engineer

This is the Himalaya.jl monorepo. The relevant package is `packages/HimalayaUI` — a Julia/Oxygen.jl REST backend. You are working in a git worktree; run all commands from the worktree root.

**Key existing code you will read and extend:**

- `packages/HimalayaUI/src/routes_samples.jl` — registers sample routes inside `register_samples_routes!()`. The existing `@get "/api/experiments/{id}/samples"` (lines 4–18) is the *sibling* of the route you are adding. It is the shape reference but **not** to be modified.
- `packages/HimalayaUI/src/routes_experiments.jl` — its `_experiment_row_to_json` helper (lines 3–20) contains the inline `q_units` extraction logic you will factor out.
- `packages/HimalayaUI/test/test_routes_samples.jl` — the test file you will extend with a new `@testset`.
- `packages/HimalayaUI/test/test_http.jl` — defines the `with_test_server(f, db)` test harness used by all route tests.

**Codebase conventions:**

- Routes return `HTTP.Response(status, ["Content-Type" => "application/json"], JSON3.write(body))`.
- `current_db()` yields the active `SQLite.DB` inside a route handler.
- `row_to_json(row)` (in `src/json.jl`) converts a `NamedTuple` DB row to a mutable `Dict{Symbol,Any}` ready for `JSON3.write`. You can add keys to its result.
- Query the DB with `Tables.rowtable(DBInterface.execute(db, sql, params))`, which returns a `Vector{NamedTuple}`.
- The backend Julia suite is slow (5–10 min full run). During TDD, iterate with the single-file command shown in the tasks below; run the full suite only once at the end.

**Running a single test file fast** (used throughout this plan):

```bash
julia --project=packages/HimalayaUI -e '
using Test, HTTP, JSON3, SQLite, DBInterface, Tables, HimalayaUI
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_routes_samples.jl")'
```

`test_http.jl` must be included first — it defines `with_test_server`, which `test_routes_samples.jl` depends on.

---

## File Structure

- **Modify** `packages/HimalayaUI/src/routes_experiments.jl` — extract `_q_units_from_config` helper; rewire `_experiment_row_to_json` to call it.
- **Modify** `packages/HimalayaUI/src/routes_samples.jl` — add the `@get "/api/samples"` route inside `register_samples_routes!()`.
- **Modify** `packages/HimalayaUI/test/test_routes_samples.jl` — add a `@testset "corpus samples route"`.

No `server.jl` change is needed: `register_samples_routes!()` is already called at `server.jl:122`, so adding a route inside that function registers it automatically.

---

## Task 1: Extract the `_q_units_from_config` helper

This is a pure refactor. The `q_units` extraction logic — parse the config TOML, read `beamline.q_units`, fall back to the ASCII default `"A-1"` on missing or malformed config — currently lives inline inside `_experiment_row_to_json`. Task 2's new route needs the same logic; extract it now so the route reuses it instead of copying it.

There is no new test for this task: the existing `test_routes_experiments.jl` already asserts the `q_units` default and the malformed-TOML fallback (it sets an unbalanced-bracket config and expects `"A-1"`). Those tests staying green *is* the verification.

**Files:**
- Modify: `packages/HimalayaUI/src/routes_experiments.jl:3-20`

- [ ] **Step 1: Read the current helper**

Read `packages/HimalayaUI/src/routes_experiments.jl` lines 1–20. Confirm `_experiment_row_to_json` contains the inline `q_units = if cfg_text isa AbstractString ...` block, and that the file already uses `TOML.parse` (the `TOML` stdlib is in scope for the module).

- [ ] **Step 2: Add the extracted helper and rewire the caller**

Replace the entire `_experiment_row_to_json` function (lines 3–20) with the helper plus a slimmed caller:

```julia
"""
    _q_units_from_config(cfg_text) -> String

Extract `beamline.q_units` from an experiment's stored config TOML.

Accepts anything (a `String`, `missing`, `nothing`): a non-string, an empty
string, or malformed TOML all fall back to the ASCII default `"A-1"` — the UI
prettifies that to "Å⁻¹". Defensive by design: one experiment's malformed
config must never 500 a list endpoint.
"""
function _q_units_from_config(cfg_text)::String
    if cfg_text isa AbstractString && !isempty(cfg_text)
        try
            bl = get(TOML.parse(cfg_text), "beamline", Dict())
            return get(bl, "q_units", "A-1")
        catch
            return "A-1"
        end
    end
    return "A-1"
end

function _experiment_row_to_json(row::NamedTuple)
    d = row_to_json(row)
    d[:q_units] = _q_units_from_config(get(d, :config, nothing))
    d
end
```

- [ ] **Step 3: Verify the experiments route tests still pass**

Run:

```bash
julia --project=packages/HimalayaUI -e '
using Test, HTTP, JSON3, SQLite, DBInterface, Tables, HimalayaUI
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_routes_experiments.jl")'
```

Expected: the `experiments routes` testset passes — in particular the assertions `body2.q_units == "A-1"` (default) and `body4.q_units == "A-1"` (malformed-config fallback). No failures, no errors.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/src/routes_experiments.jl
git commit -m "refactor: extract _q_units_from_config helper from routes_experiments

Factors the config-TOML q_units extraction (with the A-1 fallback) out
of _experiment_row_to_json so the corpus GET /api/samples route can reuse
it instead of duplicating the parse logic.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Add the corpus `GET /api/samples` route

Add the new route via TDD: write the failing test first, watch it fail because the route does not exist, then implement the route, then watch it pass.

**Files:**
- Modify: `packages/HimalayaUI/test/test_routes_samples.jl` (append a new `@testset`)
- Modify: `packages/HimalayaUI/src/routes_samples.jl` (add route inside `register_samples_routes!()`)

- [ ] **Step 1: Write the failing test**

Append this new `@testset` to the end of `packages/HimalayaUI/test/test_routes_samples.jl`, after the existing `@testset "samples routes"` block (after its final `end`):

```julia
@testset "corpus samples route" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)

    e1 = HimalayaUI.init_experiment!(db; name="E1", path="/e1",
        data_dir="/e1/d", analysis_dir="/e1/a")
    e2 = HimalayaUI.init_experiment!(db; name="E2", path="/e2",
        data_dir="/e2/d", analysis_dir="/e2/a")

    # Distinct q_units per experiment, written straight into the config blob.
    # e1 gets an explicit nm-1; e2 is left config-less so it falls back to A-1.
    DBInterface.execute(db,
        "UPDATE experiments SET config = ? WHERE id = ?",
        ["[beamline]\nq_units = \"nm-1\"\n", e1])

    s1 = HimalayaUI.create_sample!(db; experiment_id=e1, name="A1", display_name="UX-A1")
    s2 = HimalayaUI.create_sample!(db; experiment_id=e1, name="A2", display_name="UX-A2")
    s3 = HimalayaUI.create_sample!(db; experiment_id=e2, name="B1", display_name="UX-B1")

    # One tag on s1, so the projection's bundled `tags` array is exercised.
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value, source)
         VALUES (?, 'lipid', 'DOPC', 'manual')", [s1])

    with_test_server(db) do port, base
        # Full corpus: every sample across both experiments.
        r = HTTP.get("$base/api/samples")
        @test r.status == 200
        all = JSON3.read(String(r.body))
        @test length(all) == 3
        @test Set(s.name for s in all) == Set(["A1", "A2", "B1"])

        by_name = Dict(String(s.name) => s for s in all)

        # q_units sourced from each sample's owning experiment.
        @test by_name["A1"].q_units == "nm-1"
        @test by_name["A2"].q_units == "nm-1"
        @test by_name["B1"].q_units == "A-1"   # e2 has no config → default

        # tags bundled in the projection.
        @test length(by_name["A1"].tags) == 1
        @test by_name["A1"].tags[1].key   == "lipid"
        @test by_name["A1"].tags[1].value == "DOPC"
        @test by_name["A2"].tags == []

        # ?experiment_id= filter narrows to one experiment.
        r = HTTP.get("$base/api/samples?experiment_id=$e1")
        @test r.status == 200
        filtered = JSON3.read(String(r.body))
        @test length(filtered) == 2
        @test Set(s.name for s in filtered) == Set(["A1", "A2"])

        # Nonexistent experiment id → empty array (SQL gives this for free).
        r = HTTP.get("$base/api/samples?experiment_id=999999")
        @test r.status == 200
        @test JSON3.read(String(r.body)) == []

        # Malformed experiment_id → 400, not a silent ignore.
        r = HTTP.get("$base/api/samples?experiment_id=abc"; status_exception=false)
        @test r.status == 400
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:

```bash
julia --project=packages/HimalayaUI -e '
using Test, HTTP, JSON3, SQLite, DBInterface, Tables, HimalayaUI
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_routes_samples.jl")'
```

Expected: FAIL. The first `GET /api/samples` request returns `404` (no such route), so `@test r.status == 200` fails. The pre-existing `samples routes` testset still passes.

- [ ] **Step 3: Implement the route**

In `packages/HimalayaUI/src/routes_samples.jl`, add the following route inside `register_samples_routes!()`, immediately after the existing `@get "/api/experiments/{id}/samples"` block (after its closing `end` at line 18, before the `@patch "/api/samples/{id}"` block):

```julia
    # Corpus-wide sample listing — every sample across all experiments, each
    # carrying its `tags` and a `q_units` sourced from the owning experiment's
    # config. Optional `?experiment_id=` narrows to one experiment, so this
    # route subsumes the experiment-scoped `/api/experiments/{id}/samples`.
    #
    # Three queries total, never N+1: samples, experiment configs, and one
    # batched tag query grouped in memory.
    @get "/api/samples" function(req::HTTP.Request)
        db     = current_db()
        params = HTTP.queryparams(HTTP.URI(req.target))

        # Optional ?experiment_id= filter. A non-integer value is a client
        # error, not something to silently ignore.
        exp_filter = nothing
        if haskey(params, "experiment_id")
            exp_filter = tryparse(Int, params["experiment_id"])
            exp_filter === nothing && return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "experiment_id must be an integer")))
        end

        samples = exp_filter === nothing ?
            Tables.rowtable(DBInterface.execute(db,
                "SELECT * FROM samples ORDER BY id")) :
            Tables.rowtable(DBInterface.execute(db,
                "SELECT * FROM samples WHERE experiment_id = ? ORDER BY id",
                [exp_filter]))

        # experiment_id -> q_units, one TOML parse per experiment (not per sample).
        qunits_by_exp = Dict{Int, String}()
        for er in Tables.rowtable(DBInterface.execute(db,
                "SELECT id, config FROM experiments"))
            qunits_by_exp[Int(er.id)] = _q_units_from_config(er.config)
        end

        # One batched tag query, grouped by sample_id. Skipped entirely when
        # there are no samples — an empty `IN ()` is invalid SQL.
        tags_by_sample = Dict{Int, Vector{Any}}()
        if !isempty(samples)
            ids          = [Int(sm.id) for sm in samples]
            placeholders = join(fill("?", length(ids)), ", ")
            tagrows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, sample_id, key, value, source FROM sample_tags
                 WHERE sample_id IN ($placeholders) ORDER BY id", ids))
            for tr in tagrows
                push!(get!(tags_by_sample, Int(tr.sample_id), Any[]),
                      Dict(:id     => Int(tr.id), :key   => tr.key,
                           :value  => tr.value,   :source => tr.source))
            end
        end

        out = map(samples) do sm
            d          = row_to_json(sm)
            d[:tags]   = get(tags_by_sample, Int(sm.id), Any[])
            d[:q_units] = get(qunits_by_exp, Int(sm.experiment_id), "A-1")
            d
        end
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(out))
    end
```

- [ ] **Step 4: Run the test to verify it passes**

Run:

```bash
julia --project=packages/HimalayaUI -e '
using Test, HTTP, JSON3, SQLite, DBInterface, Tables, HimalayaUI
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_routes_samples.jl")'
```

Expected: PASS. Both the pre-existing `samples routes` testset and the new `corpus samples route` testset pass with no failures or errors.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_samples.jl packages/HimalayaUI/test/test_routes_samples.jl
git commit -m "feat: add corpus-wide GET /api/samples route (#154)

Returns every sample across all experiments in one projection, each
carrying its tags and a q_units value from the owning experiment's
config. Optional ?experiment_id= narrows to one experiment, so the
corpus route subsumes the experiment-scoped listing.

Three queries total (samples, experiment configs, batched tags) joined
in memory — no per-sample N+1.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Full-suite verification

The single-file runs above are fast but do not prove the rest of the backend is unaffected. Run the full suite once.

- [ ] **Step 1: Run the full HimalayaUI backend suite**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|did not pass|fail|error" /tmp/jl-test.out
tail -50 /tmp/jl-test.out
```

Expected: every testset passes. The `HimalayaUI` top-level summary reports `0` fails and `0` errors. Pay particular attention to `test_routes_experiments.jl` (Task 1's refactor target) and `test_routes_samples.jl` (Task 2's new testset).

- [ ] **Step 2: If anything fails, fix and re-run**

If the suite reports failures, read the relevant section of `/tmp/jl-test.out`, fix the cause, and re-run Step 1. Do not proceed until the full suite is green. Commit any fix with a descriptive message.

---

## Self-review notes (for the planner — not a task)

**Spec coverage** — every acceptance criterion in `docs/superpowers/specs/2026-05-17-corpus-samples-route-design.md` maps to a task:

- "`GET /api/samples` returns every sample" → Task 2, Step 1 full-corpus assertion + Step 3 route.
- "Each sample carries `q_units` and `tags`" → Task 2 assertions on `by_name[...].q_units` / `.tags`; sourced via `_q_units_from_config` (Task 1).
- "`?experiment_id=` filters; malformed → 400" → Task 2 filter, nonexistent-id, and `experiment_id=abc` assertions.
- "Existing `/api/experiments/{id}/samples` unchanged and green" → not modified; Task 2 keeps the existing `samples routes` testset; Task 3 confirms full suite.
- "A Julia route test exercises corpus, filtered, malformed responses" → Task 2, Step 1.
- "Running frontend unaffected" → route is additive and unconsumed; nothing to verify in code.

**Out-of-scope items** (pagination, the `useCorpusSamples` hook #156, `q_units` normalization, cross-route `/api/sample-tags`) are correctly absent from all tasks.

**Type consistency** — `_q_units_from_config(cfg_text)::String` is defined in Task 1 and called in Task 2 with `er.config` (a `String` or `missing`); the helper's signature is untyped and explicitly handles non-string input. `current_db`, `row_to_json`, `with_test_server`, `HimalayaUI.init_experiment!`, `HimalayaUI.create_sample!`, `HimalayaUI.create_schema!` are all pre-existing and used consistently.
