# Corpus `GET /api/sample-tags` Route (#157 / I0.2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a corpus-wide `GET /api/sample-tags` endpoint that returns every distinct `(key, value)` sample tag across all experiments in the database.

**Architecture:** A read-only `@get` handler added to the existing `register_picker_routes!()` in `routes_picker.jl`. It runs `SELECT DISTINCT key, value FROM sample_tags ORDER BY key, value` — no experiment scoping, so no `JOIN` to `samples`. Append-only: no schema change, no event log row, no SSE, no `with_idempotency`, no `server.jl` change (the `register_picker_routes!()` call already exists at `server.jl:130`).

**Tech Stack:** Julia, Oxygen.jl (routing), SQLite via DBInterface.jl, JSON3.jl. Tests use stdlib `Test` with the HimalayaUI in-process HTTP test server.

**Spec:** Master plan §3.1 (`docs/superpowers/plans/2026-05-17-himalaya-ui-redesign.md`) + issue map I0.2 (`docs/superpowers/plans/2026-05-17-himalaya-ui-redesign-issue-map.md`). No standalone spec doc — this plan's design decisions were settled in brainstorming: distinct `(key, value)` pairs, flat sorted array, response shape `[{key, value}]`. The experiment-scoped sibling `GET /api/experiments/{eid}/sample-tags` is legacy and untouched.

---

## File Structure

- **Modify** `packages/HimalayaUI/src/routes_picker.jl` — add the `@get "/api/sample-tags"` handler inside `register_picker_routes!()`, immediately after the experiment-scoped sample-tags route (after line 69, before the `picker-samples` route at line 71). Also extend the file-header docstring to mention the new corpus route.
- **Modify** `packages/HimalayaUI/test/test_picker_routes.jl` — add two `@testset`s inside the existing `@testset "picker routes"` block, before its closing `end` at line 148.

No other files change. The pre-existing staleness of the file header (it omits the `picker-samples` route) is deliberately left alone — out of scope for this issue.

---

## Task 1: Corpus `GET /api/sample-tags` route

**Files:**
- Modify: `packages/HimalayaUI/src/routes_picker.jl` (handler after line 69; header docstring lines 3-24)
- Test: `packages/HimalayaUI/test/test_picker_routes.jl` (new testsets before line 148)

- [ ] **Step 1: Write the failing tests**

In `packages/HimalayaUI/test/test_picker_routes.jl`, insert the following two testsets immediately before the final `end` on line 148 (the `end` that closes `@testset "picker routes"`). Match the existing 4-space indentation of the sibling testsets:

```julia
    @testset "GET /api/sample-tags: corpus-wide tags across experiments" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            # Second experiment with its own sample + tag. The corpus route
            # MUST include it — this is the inverse of the experiment-scoped
            # route's "only tags in scope" test.
            e2_id = HimalayaUI.init_experiment!(ctx.db; path=tmp * "/e2",
                data_dir=tmp * "/e2/data", analysis_dir=tmp * "/e2/analysis")
            s_other = HimalayaUI.create_sample!(ctx.db; experiment_id=e2_id, name="OTHER")
            DBInterface.execute(ctx.db,
                "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, 'manual')",
                [ctx.sample_id, "lipid", "DOPC"])
            DBInterface.execute(ctx.db,
                "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, 'manual')",
                [s_other, "buffer", "PBS"])
            # A duplicate (key, value) on a third sample in experiment 1 —
            # DISTINCT must collapse it to a single corpus entry.
            s3 = HimalayaUI.create_sample!(ctx.db; experiment_id=ctx.experiment_id, name="D3")
            DBInterface.execute(ctx.db,
                "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, 'manual')",
                [s3, "lipid", "DOPC"])

            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/sample-tags")
                @test r.status == 200
                tags = JSON3.read(String(r.body))
                pairs = Set([(String(t.key), String(t.value)) for t in tags])
                # Tags from BOTH experiments are present.
                @test ("lipid", "DOPC") in pairs
                @test ("buffer", "PBS") in pairs
                # DISTINCT: the duplicate (lipid, DOPC) collapses to one entry.
                @test length(tags) == 2
            end
        end
    end

    @testset "GET /api/sample-tags: empty list when no tags" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/sample-tags")
                @test r.status == 200
                @test JSON3.read(String(r.body)) == []
            end
        end
    end
```

- [ ] **Step 2: Run the tests to verify they fail**

Run (the HimalayaUI Julia suite is slow, 5–10 min — capture once, then grep the file):

```bash
cd /home/jonathanchen/projects/Himalaya.jl/.claude/worktrees/corpus-sample-tags
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -nE "picker routes|sample-tags|Test Summary|did not pass|fail|Error" /tmp/jl-test.out
```

Expected: FAIL/ERROR. The route does not exist yet, so Oxygen returns 404; `HTTP.get` raises `HTTP.StatusError` (status 404) before the `@test r.status == 200` line runs. The `picker routes` testset reports the two new testsets as `Error`.

- [ ] **Step 3: Add the corpus route handler**

In `packages/HimalayaUI/src/routes_picker.jl`, insert this handler inside `register_picker_routes!()`, immediately after the closing `end` of the experiment-scoped `@get "/api/experiments/{eid}/sample-tags"` block (current line 69) and before the `@get "/api/experiments/{eid}/picker-samples"` block (current line 71):

```julia
    @get "/api/sample-tags" function(req::HTTP.Request)
        db = current_db()
        # Corpus-wide sibling of the experiment-scoped route above: every
        # distinct (key, value) tag across the whole database, with no
        # experiment filter — so no JOIN to `samples` is needed. DISTINCT
        # collapses on the (key, value) pair, matching the per-experiment
        # route's contract.
        rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT DISTINCT key, value
               FROM sample_tags
               ORDER BY key, value"""))

        out = [Dict(:key => String(r.key), :value => String(r.value)) for r in rows]
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
    end
```

- [ ] **Step 4: Update the file-header docstring**

In `packages/HimalayaUI/src/routes_picker.jl`, the header comment (lines 3-24) currently documents only the experiment-scoped routes. Add a short paragraph describing the new corpus route. Replace this block:

```julia
#   GET /api/experiments/:eid/sample-tags
#       Distinct (key, value) pairs across every sample in the experiment,
#       used to populate the picker's tag-filter dropdown. CRITICAL: distinct
#       collapses on the (key, value) pair, NOT on value alone — two tags
#       with the same value but different keys must surface as separate
#       options. The picker treats each (key, value) as an independent
#       filter chip per the spec.
```

with:

```julia
#   GET /api/experiments/:eid/sample-tags
#       Distinct (key, value) pairs across every sample in the experiment,
#       used to populate the picker's tag-filter dropdown. CRITICAL: distinct
#       collapses on the (key, value) pair, NOT on value alone — two tags
#       with the same value but different keys must surface as separate
#       options. The picker treats each (key, value) as an independent
#       filter chip per the spec.
#
#   GET /api/sample-tags
#       Corpus-wide sibling of the route above: distinct (key, value) pairs
#       across every sample in the whole database, with no experiment
#       scoping. Feeds the redesign's corpus-level series scoping step
#       (issue map I0.2). Same (key, value) distinct contract.
```

- [ ] **Step 5: Run the tests to verify they pass**

Run:

```bash
cd /home/jonathanchen/projects/Himalaya.jl/.claude/worktrees/corpus-sample-tags
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -nE "picker routes|sample-tags|Test Summary|did not pass|fail" /tmp/jl-test.out
```

Expected: PASS. The `picker routes` testset summary shows all tests passing, including the two new `GET /api/sample-tags` testsets. The pre-existing experiment-scoped `GET /api/experiments/:eid/sample-tags` testsets remain green (acceptance criterion: sibling route unchanged).

- [ ] **Step 6: Commit**

```bash
cd /home/jonathanchen/projects/Himalaya.jl/.claude/worktrees/corpus-sample-tags
git add packages/HimalayaUI/src/routes_picker.jl packages/HimalayaUI/test/test_picker_routes.jl
git commit -m "feat: add corpus-wide GET /api/sample-tags route (#157)

Returns every distinct (key, value) sample tag across all experiments,
as a corpus sibling of the experiment-scoped sample-tags route. Read-only,
append-only: no schema change, no events, no migration. Feeds the
redesign's corpus series-scoping step (issue map I0.2)."
```

---

## Self-Review

**Spec coverage** — issue #157 acceptance criteria:
- "`GET /api/sample-tags` returns every sample tag in the database" → Task 1 Step 3 handler. ✓
- "The experiment-scoped sibling route is unchanged and green" → handler is purely additive; Step 5 grep confirms the sibling testsets stay green. ✓
- "A Julia route test covers the corpus response" → Task 1 Step 1, two testsets (cross-experiment aggregation + distinct collapse; empty list). ✓
- "The running frontend is unaffected" → no frontend files touched; route is additive. ✓
- Issue map I0.2 "register in server.jl" → registration is automatic: the handler lives inside the already-registered `register_picker_routes!()` (`server.jl:130`), so no `server.jl` edit is required. Noted in Architecture and in the commit message.

**Placeholder scan** — no TBD/TODO; every code step shows complete code.

**Type consistency** — handler returns `[{key, value}]` (`Dict(:key, :value)` per row); tests read `t.key`/`t.value` and compare `(String, String)` tuples. Consistent. Empty case returns `[]`, asserted in the second testset. Helper names (`_setup_analyzed_exposure`, `with_test_server`, `HimalayaUI.init_experiment!`, `HimalayaUI.create_sample!`) match existing usage in the same test file (lines 116-118, 78-79).
