# findpeaks fast-skip for zero-peak traces (#300) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Exposures whose trace legitimately yields **zero auto peaks** must fast-skip re-analysis exactly like every other unchanged exposure, instead of re-running `findpeaks` (plus file I/O, plus a durable `analyze_run` row on the curation-route path) on every call.

**Architecture:** Same shape as the merged #298 fix one stage down: the `autopeaks_count > 0` legs of the two `findpeaks_skipped` predicates in `analyze_exposure!` conflate "never analyzed" with "analyzed, found zero peaks". The correct signal is the stored `trace_hash`: it has exactly one writer (`pipeline.jl:1045`), inside the same transaction as `diff_update_auto_peaks!`, so a non-NULL stored hash proves `auto_peaks` (zero rows included) reflects that trace. **The two sites need *different* expressions:** the `trace_known_unchanged` leg never computes a fresh hash (`new_trace_hash` is aliased to the stored one), so it gates on `stored_trace_hash !== nothing`; the slow-path leg compares real hashes and just drops its count leg (`nothing == <String>` is `false`, protecting the never-analyzed case).

**Tech Stack:** Julia (SQLite.jl, stdlib Test), no frontend changes, no schema changes.

## Global Constraints

- Work in the worktree `.claude/worktrees/fix-300-findpeaks-fast-skip` (already created; this plan lives on its branch).
- **Worktree Julia setup MUST dev the local core before instantiate** (bare instantiate pulls registry Himalaya 0.5.3, which lacks `bonnet_lattice` → events bucket errors):
  `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'` (run from worktree root).
- All commands below run from the worktree root unless stated otherwise.
- Never run two Julia backend suites concurrently (port races, shared logs).
- Known environment flakes to NOT misdiagnose: the P99-latency testset (`test_fast_skip.jl:105`) and test-server port-readiness ("did not become ready within 5.0s") fail intermittently regardless of the change.
- TDD: red-test commit and fix commit are separate.
- Commit trailers per session convention (`Co-Authored-By: Claude ...`).

## Verified facts the plan relies on (do not re-derive)

- `Himalaya.findpeaks(q, I, σ)` returns **zero peaks with no error** on constant, constant+noise, smooth-decay, and decay+noise traces (probed 2026-07-17 on core `main`). Result field for count: `length(result.q)`.
- `Himalaya.indexpeaks(Float64[], Float64[])` returns an empty candidate vector, no error.
- `load_dat` = `readdlm(path, Float64)`, requires ≥3 columns; returns columns 1,2,3 as q, I, σ (`packages/HimalayaUI/src/datfile.jl:10-14`).
- `trace_hash` single writer: `packages/HimalayaUI/src/pipeline.jl:1045`, inside `if !findpeaks_skipped` in the one outer tx, adjacent to `diff_update_auto_peaks!`. `auto_peaks` rows are deleted only by diff-update (`pipeline.jl:223`) or with the whole exposure.
- First-ever analyze stamps `trace_hash` AND `analysis_inputs_hash` in the same tx (persist always runs — stored inputs hash NULL never matches), so "trace_hash non-NULL but inputs_hash NULL" cannot exist post-tx.
- Empty persists are battle-tested: prod has 1202 analyzed zero-auto-peak exposures (of 4280 total); pre-fix, the delete-last-add-curation flow already persists an empty set via a fresh empty findpeaks result.
- Prod damage evidence for the PR body: 83 curations sit on zero-auto-peak exposures; worst zero-peak exposure carries 59 durable `analyze_run` rows (vs max 23 for any with-peak exposure; avg 1.2).
- `count_auto_peaks` MUST NOT be deleted — `test_fast_skip.jl` calls `HimalayaUI.count_auto_peaks` in the add-curation testset. Only the now-unused local read inside `analyze_exposure!` goes.
- `test_fast_skip.jl` runs in the `events` GROUP bucket.

---

### Task 1: Failing tests — flat-trace fixture + three red testsets

**Files:**
- Modify: `packages/HimalayaUI/test/test_fast_skip.jl` (append after the `"fast-skip: skipped when exclude collapses indexing to zero indices (#297)"` testset, before `"fast-skip: load_dat REQUIRED when add curation present"`)

**Interfaces:**
- Consumes: `seed_experiment!` (from `test/seed.jl`, already included), `analyze_exposure!`, `HimalayaUI.count_auto_peaks`, `HimalayaUI.read_trace_hash`, `HimalayaUI.read_inputs_hash`, `latest_analyze_run` (file-local helper at `test_fast_skip.jl:43`, returns `(id, payload)`).
- Produces: file-local helpers `setup_flat_analyzed_exposure(tmp) -> (db, exposure_id, analysis_dir, dat_path)` and `analyze_run_count(db, exposure_id) -> Int`.

- [ ] **Step 1: Add the fixture helper and three testsets**

Insert the following block into `packages/HimalayaUI/test/test_fast_skip.jl` immediately after the `end` of the `"...zero indices (#297)"` testset:

```julia
# ──────────────────────────────────────────────────────────────────────────────
# #300: a trace that legitimately yields ZERO auto peaks must fast-skip like
# any other unchanged exposure. Fixture: a flat (featureless) trace — verified
# to produce zero findpeaks results with no error.
# ──────────────────────────────────────────────────────────────────────────────
function setup_flat_analyzed_exposure(tmp::String; name="FlatExp", stem="FLAT01")
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    mkpath(joinpath(tmp, "data"))
    dat_path = joinpath(analysis_dir, stem * ".dat")
    open(dat_path, "w") do io
        for qv in range(0.05, 0.6, length=800)
            println(io, "$(qv) 100.0 1.0")   # q  I(constant)  σ
        end
    end

    db = open_db(joinpath(tmp, "himalaya.db"))
    seeded = seed_experiment!(db, tmp;
        name = name, analysis_dir = analysis_dir,
        stems = [stem], experiment_type = "simple")
    exposure_id = seeded.exposure_ids[1]

    # First analyze: findpeaks finds zero peaks; trace_hash + inputs_hash stamp.
    analyze_exposure!(db, exposure_id, analysis_dir)
    return (db=db, exposure_id=exposure_id, analysis_dir=analysis_dir,
            dat_path=dat_path)
end

analyze_run_count(db, exposure_id) = first(Tables.rowtable(DBInterface.execute(db,
    "SELECT COUNT(*) AS c FROM user_actions WHERE entity_id = ? AND action = 'analyze_run'",
    [exposure_id]))).c

@testset "fast-skip: zero-peak trace fast-skips with zero file I/O (#300)" begin
    mktempdir() do tmp
        ctx = setup_flat_analyzed_exposure(tmp)

        # Preconditions (loud tripwires, like #298's n_indices == 0): the flat
        # fixture really yields zero peaks, and the first analyze stamped both
        # hashes. If core peak-finding changes make this non-zero, this testset
        # no longer covers #300 — adjust the fixture, don't delete the asserts.
        @test HimalayaUI.count_auto_peaks(ctx.db, ctx.exposure_id) == 0
        @test HimalayaUI.read_trace_hash(ctx.db, ctx.exposure_id) !== nothing
        @test HimalayaUI.read_inputs_hash(ctx.db, ctx.exposure_id) !== nothing

        n_before = analyze_run_count(ctx.db, ctx.exposure_id)

        # Hard proof of zero file I/O: delete the .dat. If the fast path fails
        # to engage, analyze_exposure! falls to the slow path and errors
        # ("dat file not found") — a red signal, not a silent pass.
        rm(ctx.dat_path)
        analyze_exposure!(ctx.db, ctx.exposure_id, ctx.analysis_dir;
                          trace_known_unchanged=true)
        @test analyze_run_count(ctx.db, ctx.exposure_id) == n_before
    end
end

@testset "fast-skip: zero-peak trace skips findpeaks on the slow path (#300)" begin
    mktempdir() do tmp
        ctx = setup_flat_analyzed_exposure(tmp)

        # Plain (trace_known_unchanged=false) re-analyze: the file is re-hashed
        # (unavoidable), but findpeaks must NOT re-run — the stored trace_hash
        # matches and auto_peaks (zero rows) already reflect this trace.
        analyze_exposure!(ctx.db, ctx.exposure_id, ctx.analysis_dir)
        run = latest_analyze_run(ctx.db, ctx.exposure_id)
        @test run.payload[:findpeaks_skipped]  === true
        @test run.payload[:indexpeaks_skipped] === true
    end
end

@testset "fast-skip: empty persist after last add curation removed (#300)" begin
    mktempdir() do tmp
        ctx = setup_flat_analyzed_exposure(tmp)
        empty_inputs_hash = HimalayaUI.read_inputs_hash(ctx.db, ctx.exposure_id)

        # Add curation → slow path by design (sharpness sampled from trace).
        DBInterface.execute(ctx.db,
            "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', ?)",
            [ctx.exposure_id, 0.2])
        analyze_exposure!(ctx.db, ctx.exposure_id, ctx.analysis_dir)
        @test HimalayaUI.read_inputs_hash(ctx.db, ctx.exposure_id) != empty_inputs_hash

        # Remove it. The re-analyze must skip findpeaks (trace unchanged, hash
        # stamped) and persist the EMPTY effective set via
        # synthesize_peaks_result — the path #300 makes newly common.
        DBInterface.execute(ctx.db,
            "DELETE FROM peak_curations WHERE exposure_id = ? AND kind = 'add'",
            [ctx.exposure_id])
        analyze_exposure!(ctx.db, ctx.exposure_id, ctx.analysis_dir;
                          trace_known_unchanged=true)
        run = latest_analyze_run(ctx.db, ctx.exposure_id)
        @test run.payload[:findpeaks_skipped] === true
        # Round-trip: back to the exact empty-set hash from the first analyze.
        @test HimalayaUI.read_inputs_hash(ctx.db, ctx.exposure_id) == empty_inputs_hash

        # And the settled state is a full no-op again.
        n_before = analyze_run_count(ctx.db, ctx.exposure_id)
        analyze_exposure!(ctx.db, ctx.exposure_id, ctx.analysis_dir;
                          trace_known_unchanged=true)
        @test analyze_run_count(ctx.db, ctx.exposure_id) == n_before
    end
end
```

- [ ] **Step 2: Run the three testsets to verify they fail for the right reason**

Run (from worktree root):
```bash
cd packages/HimalayaUI && julia --project=. -e \
  'using Test, HimalayaUI; @testset "outer" begin include("test/test_fast_skip.jl") end' \
  2>&1 | grep -E "300\)|Test Failed|Evaluated|dat file not found" | head -20
```

Expected failures (the P99-latency testset may also fail — known flake, ignore):
- Testset 1 errors with `dat file not found` (fast path failed to engage — falls to slow path, file deleted). An *error*, not a Fail, is the expected red here.
- Testset 2 fails `run.payload[:findpeaks_skipped] === true` (evaluates `false === true`).
- Testset 3 fails `run.payload[:findpeaks_skipped] === true` (evaluates `false === true`).

If testset 1's `count_auto_peaks == 0` precondition fails instead, STOP: the flat-fixture assumption broke; do not proceed to Task 2.

- [ ] **Step 3: Commit the red tests**

```bash
git add packages/HimalayaUI/test/test_fast_skip.jl
git commit -m "test(fast-skip): zero-peak traces must fast-skip without findpeaks reruns (#300)"
```

---

### Task 2: The predicate fix

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl:928-929` (docstring), `:963-965` (DB-only reads), `:969-971` (trace_known_unchanged leg), `:994-996` (slow-path leg). Line numbers are as of merge commit `85fc1203`.

**Interfaces:**
- Consumes: `read_trace_hash` (returns `Union{String, Nothing}`, NULL → `nothing`).
- Produces: no signature changes; `findpeaks_skipped::Bool` semantics change as described.

- [ ] **Step 1: Apply the four edits**

Edit 1 — docstring (lines 928-929). Replace:
```julia
Hash-guarded: findpeaks is skipped when `trace_hash` matches the persisted
value AND auto_peaks already exist. indexpeaks is skipped when
```
with:
```julia
Hash-guarded: findpeaks is skipped when `trace_hash` matches the persisted
value — a match proves auto_peaks (possibly zero rows; featureless traces
legitimately yield no peaks, #300) reflects this exact trace, because the
hash's only writer sits in the same transaction as diff_update. indexpeaks
is skipped when
```

Edit 2 — remove the now-unused local read (line 965). Replace:
```julia
    stored_trace_hash  = read_trace_hash(db, exposure_id)
    stored_inputs_hash = read_inputs_hash(db, exposure_id)
    autopeaks_count    = count_auto_peaks(db, exposure_id)
```
with:
```julia
    stored_trace_hash  = read_trace_hash(db, exposure_id)
    stored_inputs_hash = read_inputs_hash(db, exposure_id)
```
(Do NOT delete the `count_auto_peaks` function itself — tests use it.)

Edit 3 — trace_known_unchanged leg (lines 969-971). Replace:
```julia
    if trace_known_unchanged
        new_trace_hash = stored_trace_hash
        findpeaks_skipped = autopeaks_count > 0
    end
```
with:
```julia
    if trace_known_unchanged
        new_trace_hash = stored_trace_hash
        # #300: no fresh hash exists on this path (the caller asserts the file
        # is unchanged), so hash PRESENCE is the "has been analyzed" signal:
        # trace_hash's only writer commits with diff_update, so non-nothing
        # proves auto_peaks — zero rows included — reflect this trace.
        findpeaks_skipped = stored_trace_hash !== nothing
    end
```

Edit 4 — slow-path leg (line 996). Replace:
```julia
        findpeaks_skipped = (stored_trace_hash == new_trace_hash) && (autopeaks_count > 0)
```
with:
```julia
        # #300: hash match alone (mirrors the #297 indexpeaks predicate) —
        # never-analyzed exposures can't match (stored is nothing).
        findpeaks_skipped = stored_trace_hash == new_trace_hash
```

- [ ] **Step 2: Run the three new testsets — verify green**

Same command as Task 1 Step 2. Expected: no `(#300)` testset failures; P99-latency flake permitted.

- [ ] **Step 3: Grep for stale references**

```bash
grep -rn "autopeaks_count" packages/HimalayaUI/src packages/HimalayaUI/test docs
```
Expected: no hits inside `analyze_exposure!`; test references to `count_auto_peaks` remain valid. If any comment still describes the count guard, update it in this commit.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/src/pipeline.jl
git commit -m "fix(analyze): drop autopeaks_count>0 from findpeaks fast-skip (closes #300)"
```

---

### Task 3: Full-suite verification

**Files:** none (verification only).

- [ ] **Step 1: Run the owning bucket serially first**

```bash
GROUP=events HIMALAYA_SUITE_PARALLEL=1 julia --project=packages/HimalayaUI \
  -e 'using Pkg; Pkg.test("HimalayaUI")' > build/test-events-300.log 2>&1; echo "exit=$?"
```
Expected: `exit=0`.

- [ ] **Step 2: Full parallel suite**

```bash
make test-parallel > build/suite-300.log 2>&1; echo "exit=$?"
```
Expected: `exit=0`, or failures limited to the two documented environment flakes (P99 latency; "did not become ready within 5.0s" port races). Re-run any flake-failed bucket serially to confirm.

- [ ] **Step 3: Core suite sanity**

```bash
julia --project=. -e 'using Pkg; Pkg.test()' > build/core-300.log 2>&1; echo "exit=$?"
```
Expected: `exit=0` (no core changes; guards against accidental core edits).

---

### Task 4: Domain reviews + PR

**Files:** none (process).

- [ ] **Step 1: Dispatch `himalaya-reviewer` and `saxs-physics-reviewer`** on the two commits (per the #268/#298 precedent for fast-skip predicate changes). Instruct: static review only, no test runs, return findings in the final message. Scrutiny prompts: (a) any writer/deleter of `trace_hash`/`auto_peaks` outside the diff-update tx; (b) the newly-common `synthesize_peaks_result`-on-empty persist path; (c) fixture soundness of "flat trace → zero peaks".

- [ ] **Step 2: Fix anything they find** (new commit per finding class), re-running Task 3 Step 1 after any src change.

- [ ] **Step 3: Push + draft PR**

```bash
git push -u origin HEAD
gh pr create --draft --title "fix(analyze): fast-skip zero-peak traces on trace-hash presence" --body-file <body>
```
PR body must contain `Closes #300`, the two-different-expressions rationale (hash presence vs hash equality), the prod impact numbers (1202/4280 zero-peak analyzed exposures; 83 curations on peakless exposures; max 59 `analyze_run` rows on one vs 23 for any with-peak exposure), and the verification evidence.

---

## Self-review notes

- Spec coverage: #300's three pointers (trace_known_unchanged leg, slow-path mirror, missing fixture) map to Task 2 Edits 3-4 and Task 1's fixture. The issue's "needs its own decision" is resolved as: hash presence on the assert-path, hash equality on the compute-path.
- The `payload[:findpeaks_skipped]` key exists today (`pipeline.jl:1101`) and `latest_analyze_run` (test helper, `test_fast_skip.jl:43`) returns `(id, payload)` — types consistent with the test code above.
- `open_db`, `seed_experiment!`, `Tables`, `DBInterface` are already imported at the top of `test_fast_skip.jl` — no new imports needed. `analyze_run_count` is a new file-local helper (earlier testsets inline the same query; leave them as-is to keep the diff additive).
