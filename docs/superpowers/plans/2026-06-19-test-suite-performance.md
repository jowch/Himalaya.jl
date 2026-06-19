# HimalayaUI Test-Suite Performance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cut the HimalayaUI Julia backend test suite from ~20 min to under ~3 min by replacing per-test Oxygen HTTP-server boots with in-process `internalrequest` dispatch, amortizing DB construction with a template clone, and sharding the suite across processes — without losing coverage or isolation.

**Architecture:** Four independently-landable milestones gated by a green full-suite run that must match a captured baseline pass count. M0 builds and *proves* the in-process dispatch primitive (a differential harness asserting wire ≡ in-process). M1 mechanically migrates ~25 route-test files to it, keeping a small wire-keeper tier. M2 swaps the `open_db`-backed fixtures to a `VACUUM INTO` template clone. M3 carves the suite into named `GROUP` buckets runnable in parallel.

**Tech Stack:** Julia, stdlib `Test`, Oxygen.jl 1.10.x (`internalrequest`/`serialize=true`), HTTP.jl, SQLite.jl (WAL, `VACUUM INTO`), JSON3.

## Global Constraints

- **No `src/` behavior changes.** This is a test-layer refactor. The only production code touched is test files, `test/runtests.jl`, a new `test/` harness file, and `Makefile`/`scripts` for the parallel runner. (If a src change ever looks necessary, stop and escalate.)
- **No new package dependencies.** Test-time only; everything used already ships.
- **Isolation scrub is mandatory** in every in-process helper `finally`: `HimalayaUI._DB_REF[] = nothing`; `HimalayaUI.SSE_SUBSCRIBERS[] = []` under `HimalayaUI.SSE_LOCK`; `empty!(HimalayaUI.OP_LOCKS)`. **Never** call `Oxygen.resetstate()` in the in-process path (it wipes the router the next file needs).
- **Dispatch with `serialize=true, catch_errors=true`.** `serialize=false` drops `DefaultSerializer`, so thrown handlers (`routes_peaks.jl:167`, FK `SQLiteException`s) would propagate out of `call()` instead of returning 4xx/5xx (spec blocker B2).
- **Route registration is a liveness check, not a sticky boolean** (spec blocker B1): re-`register_routes!()` whenever a probe of `/api/health` is not 200, because a wire keeper's `resetstate()` wipes the router mid-suite under `GROUP=All`.
- **Per-file commits.** Each migrated file is its own commit with green-before (wire) → swap → green-after (in-process), test count unchanged.
- **`GROUP=All` must reproduce the current `runtests.jl` include order exactly** — the serial bisect fallback.
- **Slow-suite protocol** (`packages/HimalayaUI/test/AGENTS.md`): capture the full run once to a file, then grep it. Never re-run with different greps.

Internal (non-exported) symbols are accessed `HimalayaUI.<name>` in tests.

---

## File Structure

- **Create** `packages/HimalayaUI/test/test_inproc.jl` — the `with_inproc_routes` helper + liveness check + scrub. One responsibility: in-process dispatch fixture.
- **Create** `packages/HimalayaUI/test/test_inproc_equivalence.jl` — the differential harness (wire ≡ in-process matrix). Kept permanently.
- **Create** `packages/HimalayaUI/test/test_template_db.jl` — the M2 template builder + `open_prepared_clone`. (No `with_legacy_db` helper: regime-a and legacy-migration testsets are simply left on their existing `SQLite.DB()`/`open_db` fixtures, not routed through an opt-out.)
- **Create** `packages/HimalayaUI/test/test_fixtures.jl` — shared cross-file test helpers extracted in M3.0 so GROUP buckets don't break on missing symbols.
- **Modify** `packages/HimalayaUI/test/test_http.jl` — include `test_inproc.jl` so the helper rides the same load order as `with_test_server`.
- **Modify** ~25 `packages/HimalayaUI/test/test_routes_*.jl` / `test_picker_*.jl` / etc. — mechanical `with_test_server` → `with_inproc_routes` swap (M1).
- **Modify** `packages/HimalayaUI/test/runtests.jl` — per-`@testset` timing (M0); `GROUP` buckets (M3).
- **Modify** `Makefile` — `test-parallel` target (M3).

---

## Milestone M0 — Baseline + proven primitive

### Task M0.1: Capture the regression-floor baseline

**Files:**
- Create: `packages/HimalayaUI/test/BASELINE.md` (committed record)

**Interfaces:**
- Produces: a committed baseline pass count + wall time that every later milestone gate compares against.

- [ ] **Step 1: Run the full suite once, capturing output and wall time**

```bash
cd /Users/me/projects/Himalaya.jl/.claude/worktrees/test-suite-perf
/usr/bin/time -p julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-baseline.out 2>&1; echo "exit=$?"
```
Expected: `exit=0`. (If non-zero, the suite is already red — stop and surface that before any refactor.)

- [ ] **Step 2: Extract the pass count and wall time**

```bash
grep -E "Test Summary|Pass|Fail|Total" /tmp/jl-baseline.out | tail -20
grep -E "^real " /tmp/jl-baseline.out
```
Expected: a top-level `HimalayaUI` Test Summary row with a `Pass` total and `0` fails, plus a `real <seconds>` line.

- [ ] **Step 3: Record the baseline**

Write `packages/HimalayaUI/test/BASELINE.md` with the exact numbers, e.g.:

```markdown
# Test-suite baseline (pre-perf-refactor)

Captured: <date>, commit <sha>
- Total Pass: <N>   Fail: 0   Error: 0   Broken: <B>
- Wall time (real): <S> s
- Command: julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'

Every milestone gate must show Pass == <N> (and Fail==Error==0).
```

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/test/BASELINE.md
git commit -m "test(perf): capture pre-refactor suite baseline (M0.1)"
```

### Task M0.2: Per-`@testset` timing report in runtests.jl

**Files:**
- Modify: `packages/HimalayaUI/test/runtests.jl`

**Interfaces:**
- Produces: a printed, sorted "slowest includes" table after the suite, driving M3 group balancing. No behavior/count change.

- [ ] **Step 1: Wrap each include with `@elapsed` into a timing dict**

Replace the body of the `@testset "HimalayaUI"` block so each `include(f)` is timed. Keep the *exact same file order*. Pattern:

```julia
using Test

const _TEST_TIMES = Dict{String,Float64}()
_timed_include(f) = (_TEST_TIMES[f] = @elapsed include(f))

@testset "HimalayaUI" begin
    _timed_include("test_config.jl")
    _timed_include("test_db.jl")
    # ... every existing include in the SAME order, via _timed_include ...
    _timed_include("test_spa_fallback.jl")
end

# Print slowest-first; harmless under Pkg.test, informative for M3 balancing.
let rows = sort(collect(_TEST_TIMES); by = last, rev = true)
    println("\n── per-file test timing (slowest first) ──")
    for (f, t) in rows
        println(rpad(f, 42), lpad(round(t; digits=1), 8), " s")
    end
end
```

- [ ] **Step 2: Run the suite; confirm same pass count + a timing table**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-m0_2.out 2>&1; echo "exit=$?"
grep -E "per-file test timing|Test Summary" /tmp/jl-m0_2.out | head
```
Expected: `exit=0`, the timing table prints, and the `HimalayaUI` Pass total equals the M0.1 baseline.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/test/runtests.jl
git commit -m "test(perf): per-file @elapsed timing report in runtests (M0.2)"
```

### Task M0.3: Build the `with_inproc_routes` helper

**Files:**
- Create: `packages/HimalayaUI/test/test_inproc.jl`
- Modify: `packages/HimalayaUI/test/test_http.jl`

**Interfaces:**
- Produces: `with_inproc_routes(f, db::SQLite.DB)` — calls `f(call)` where `call(method, target; headers=Pair[], body=UInt8[]) -> HTTP.Response` dispatches in-process. Scrubs `_DB_REF`/`SSE_SUBSCRIBERS`/`OP_LOCKS` in `finally`. Used by M1 migrations and the M0.4 harness.
- Consumes: `HimalayaUI.register_routes!`, `HimalayaUI.bind_db!`, `Oxygen.internalrequest`, `HimalayaUI.SSE_LOCK`, `HimalayaUI.SSE_SUBSCRIBERS`, `HimalayaUI.OP_LOCKS`, `HimalayaUI._DB_REF` (all verified to exist).

- [ ] **Step 1: Write the failing test**

Create the test in a **scratch file** (NOT inline in `test_inproc.jl`). Critical: `test_http.jl` and its transitive includes (`test_inproc.jl`, `test_template_db.jl`, `test_fixtures.jl`) are included once **per GROUP bucket** in M3, so any top-level `@testset` in them would run 3× and break the M3 Pass-sum invariant. These files must define helpers only — no `@testset`. Use the scratch file as the failing test first — create `packages/HimalayaUI/test/_scratch_inproc_test.jl`:

```julia
using Test, HTTP, SQLite
using HimalayaUI
include("test_inproc.jl")

@testset "with_inproc_routes smoke" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        with_inproc_routes(db) do call
            r = call("GET", "/api/health")
            @test r.status == 200
            @test occursin("ok", String(r.body))
        end
        # Scrub ran: no DB bound after the block.
        @test HimalayaUI._DB_REF[] === nothing
    end
end
```

- [ ] **Step 2: Run it to verify it fails (helper not defined yet)**

```bash
cd packages/HimalayaUI && julia --project=. test/_scratch_inproc_test.jl 2>&1 | tail -20
```
Expected: FAIL/ERROR — `test_inproc.jl` does not exist / `with_inproc_routes` not defined.

- [ ] **Step 3: Write the helper**

Create `packages/HimalayaUI/test/test_inproc.jl`:

```julia
# In-process route dispatch fixture — the fast analogue of with_test_server.
# Dispatches HTTP.Requests straight through Oxygen.internalrequest (no socket),
# against the same global Oxygen CONTEXT[] that register_routes! populates.
using HTTP, SQLite
import Oxygen
using HimalayaUI

# Router-liveness check (spec blocker B1). resetstate() — called by BOTH
# start_test_server! and stop_test_server! — swaps CONTEXT[] for a fresh
# empty-router ServerContext, so a wire keeper running before an in-process
# file under GROUP=All wipes the router. Probe a route register_routes! always
# mounts (/api/health) and re-register only when it's gone. Never a sticky bool.
function _ensure_inproc_routes!()
    probe = Oxygen.internalrequest(HTTP.Request("GET", "/api/health");
                                   serialize = true, catch_errors = true)
    probe.status == 200 || HimalayaUI.register_routes!()
    nothing
end

"""
    with_inproc_routes(f, db)

Bind `db` as the singleton, ensure routes are registered against the live
Oxygen context, and pass `f` a `call(method, target; headers, body)` closure
that returns the `HTTP.Response` from `Oxygen.internalrequest` (serialize=true
so thrown handlers become 4xx/5xx exactly as production serve() does).

Mirrors `stop_test_server!`'s scrub in `finally` (minus resetstate()): nulls
`_DB_REF`, clears `SSE_SUBSCRIBERS` under `SSE_LOCK`, empties `OP_LOCKS`.
"""
function with_inproc_routes(f, db::SQLite.DB)
    _ensure_inproc_routes!()
    HimalayaUI.bind_db!(db)
    call = function (method::AbstractString, target::AbstractString;
                     headers = Pair{String,String}[], body = UInt8[])
        req = HTTP.Request(method, target, headers, body)
        return Oxygen.internalrequest(req; serialize = true, catch_errors = true)
    end
    try
        f(call)
    finally
        HimalayaUI._DB_REF[] = nothing
        lock(HimalayaUI.SSE_LOCK) do
            HimalayaUI.SSE_SUBSCRIBERS[] = []
        end
        empty!(HimalayaUI.OP_LOCKS)
    end
end
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
cd packages/HimalayaUI && julia --project=. test/_scratch_inproc_test.jl 2>&1 | tail -20
```
Expected: PASS (`with_inproc_routes smoke`).

- [ ] **Step 5: Wire the helper into the shared load order, remove scratch**

Add to the END of `packages/HimalayaUI/test/test_http.jl` (so `with_inproc_routes` loads alongside `with_test_server`):

```julia
include("test_inproc.jl")
```
Then delete the scratch file:

```bash
rm packages/HimalayaUI/test/_scratch_inproc_test.jl
```

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/test/test_inproc.jl packages/HimalayaUI/test/test_http.jl
git commit -m "test(perf): add with_inproc_routes in-process dispatch helper (M0.3)"
```

### Task M0.4: Differential equivalence harness (wire ≡ in-process)

**Files:**
- Create: `packages/HimalayaUI/test/test_inproc_equivalence.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl` (include it after `test_http.jl`)

**Interfaces:**
- Consumes: `with_test_server`, `with_inproc_routes`, `HimalayaUI.open_db`, fixture builders (`init_experiment!`, `create_sample!`, `create_exposure!`, `analyze_exposure!`).
- Produces: a permanent test asserting that, for a representative request matrix, the wire `HTTP.Response` and the in-process `HTTP.Response` are equivalent (status + body identical; headers identical after dropping a transport-only allowlist).

- [ ] **Step 1: Write the equivalence helper + the matrix as the failing test**

Create `packages/HimalayaUI/test/test_inproc_equivalence.jl`:

```julia
# Differential equivalence harness: proves Oxygen.internalrequest (in-process)
# produces the SAME HTTP.Response as a real socket request, for the shapes most
# likely to diverge. This is the single load-bearing equivalence proof (spec
# §5.1). Keep it green permanently.
using Test, HTTP, SQLite, JSON3, DBInterface
using FileIO, ImageCore, TiffImages   # for the binary-PNG contract row
using HimalayaUI

if !isdefined(@__MODULE__, :with_test_server)
    include("test_http.jl")
end

# Drop only transport-only headers; assert every other header matches so
# contract headers (Content-Type, Cache-Control, X-Image-*) cannot slip.
const _TRANSPORT_HEADERS = Set(lowercase.(
    ["Date", "Server", "Transfer-Encoding", "Connection", "Content-Length"]))

_sig_headers(r::HTTP.Response) = sort([
    lowercase(k) => v for (k, v) in r.headers
    if !(lowercase(k) in _TRANSPORT_HEADERS)])

# Normalize body to a String regardless of representation: the wire path returns
# bytes (Vector{UInt8}) parsed off the socket; the in-process path may leave the
# body as the raw String set by format_response! (misc.jl:402). copy(::String)
# would throw, so branch on the type.
_body(r::HTTP.Response) = r.body isa AbstractString ? String(r.body) : String(copy(r.body))

"Assert a wire response and an in-process response are equivalent."
function _assert_equiv(wire::HTTP.Response, inproc::HTTP.Response; label="")
    @test wire.status == inproc.status
    @test _body(wire) == _body(inproc)
    @test _sig_headers(wire) == _sig_headers(inproc)
end

# A fully-seeded fixture DB (its own temp dir) that both transports run against.
# NOTE: idempotency-replay rows below build their OWN per-transport DB so the
# second transport can't replay the first's cached row.
function _seed(dir)
    analysis_dir = joinpath(dir, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db = HimalayaUI.open_db(joinpath(dir, "h.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=dir,
        data_dir=joinpath(dir, "data"), analysis_dir=analysis_dir)
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
    (; db, exp_id, s_id, e_id)
end

@testset "in-process ≡ wire equivalence" begin

    @testset "GET 200 JSON (list)" begin
        mktempdir() do d
            fx = _seed(d)
            w = with_test_server(fx.db) do port, base
                HTTP.get("$base/api/experiments/$(fx.exp_id)/samples")
            end
            ip = with_inproc_routes(fx.db) do call
                call("GET", "/api/experiments/$(fx.exp_id)/samples")
            end
            _assert_equiv(w, ip; label="list")
        end
    end

    @testset "GET binary PNG image — exercises the X-Image-*/Cache-Control contract (B-2)" begin
        # Must seed a REAL image_path or the route 404s on BOTH transports and the
        # 200 PNG branch (the only code that emits image/png + X-Image-*/Cache-Control)
        # is never reached. TIFF-seeding mirrors test_routes_image.jl:4-19. Use the
        # FULL image (not thumb) so X-Image-Width/Height are present, and so the
        # response is recomputed deterministically each call (same source ⇒ identical
        # PNG bytes across transports).
        mktempdir() do d
            tiff = joinpath(d, "img.tiff")
            save(tiff, Gray.(rand(Float32, 512, 384)))
            db = HimalayaUI.open_db(joinpath(d, "h.db"))   # M0.4 predates M2; use open_db
            exp_id = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
            e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, image_path=tiff)
            path = "/api/exposures/$e_id/image?thumb=0"   # explicit full branch; exercises query parsing
            w  = with_test_server(db) do port, base; HTTP.get("$base$path"; status_exception=false) end
            ip = with_inproc_routes(db) do call; call("GET", path) end
            @test w.status == 200                          # guard: the contract branch was actually hit
            @test HTTP.header(w, "Content-Type") == "image/png"
            @test HTTP.header(w, "X-Image-Width") != ""    # full branch only
            _assert_equiv(w, ip; label="image-png")        # body bytes + X-Image-*/Cache-Control parity
        end
    end

    @testset "uncaught throw → 500 (same status both paths)" begin
        # Drive a route that raises (FK violation → SQLiteException). serialize=true
        # must wrap it into the same status the wire produces.
        mktempdir() do d
            fx = _seed(d)
            hdrs = ["Content-Type" => "application/json", "X-Username" => "alice"]
            bodyj = JSON3.write(Dict(:body => "orphan"))
            w = with_test_server(fx.db) do port, base
                HTTP.post("$base/api/samples/99999/messages"; body=bodyj, headers=hdrs, status_exception=false)
            end
            ip = with_inproc_routes(fx.db) do call
                call("POST", "/api/samples/99999/messages";
                     headers=hdrs, body=Vector{UInt8}(bodyj))
            end
            @test w.status >= 400
            _assert_equiv(w, ip; label="fk-throw")
        end
    end

    @testset "DELETE verb parity (404 on missing id, no body)" begin
        # Exercises the DELETE method + a no-body request + an error-response
        # body, asserting parity. Uses a REAL delete route (/api/peaks/{id},
        # routes_peaks.jl:348) with a nonexistent id so no mutation/setup is
        # needed and both transports see the identical not-found response.
        mktempdir() do d
            fx = _seed(d)
            hdrs = ["X-Username" => "alice"]
            w  = with_test_server(fx.db) do port, base
                HTTP.request("DELETE", "$base/api/peaks/99999", hdrs; status_exception=false)
            end
            ip = with_inproc_routes(fx.db) do call
                call("DELETE", "/api/peaks/99999"; headers=hdrs)
            end
            _assert_equiv(w, ip; label="delete-404")
        end
    end

    @testset "idempotency replay: identical body, per-transport DB" begin
        # Each transport gets its OWN seeded DB + its OWN op_id, so neither
        # replays the other's cached row (idempotency.jl:108-110).
        function _replay(call_with, opid)
            hdrs = ["Content-Type"=>"application/json", "X-Username"=>"alice",
                    "X-Client-Id"=>"tab-1", "X-Client-Op-Id"=>opid]
            bodyj = JSON3.write(Dict(:body => "hello"))
            call_with(hdrs, bodyj)
        end
        wire_bodies = mktempdir() do d
            fx = _seed(d)
            with_test_server(fx.db) do port, base
                r1 = _replay((h,b)->HTTP.post("$base/api/samples/$(fx.s_id)/messages"; body=b, headers=h), "op-wire")
                r2 = _replay((h,b)->HTTP.post("$base/api/samples/$(fx.s_id)/messages"; body=b, headers=h), "op-wire")
                (r1.status, _body(r1), _body(r2))
            end
        end
        inproc_bodies = mktempdir() do d
            fx = _seed(d)
            with_inproc_routes(fx.db) do call
                r1 = _replay((h,b)->call("POST", "/api/samples/$(fx.s_id)/messages"; headers=h, body=Vector{UInt8}(b)), "op-ip")
                r2 = _replay((h,b)->call("POST", "/api/samples/$(fx.s_id)/messages"; headers=h, body=Vector{UInt8}(b)), "op-ip")
                (r1.status, _body(r1), _body(r2))
            end
        end
        @test wire_bodies[1] == inproc_bodies[1]          # same status
        @test wire_bodies[2] == wire_bodies[3]            # wire: replay identical
        @test inproc_bodies[2] == inproc_bodies[3]        # in-process: replay identical
    end

    @testset "GET numeric-array JSON body (trace) — serialization parity" begin
        # The trace route returns numeric arrays (q/I) — the shape most prone to
        # JSON serialization divergence (float precision, key order). routes_trace.jl:4.
        mktempdir() do d
            fx = _seed(d)
            w  = with_test_server(fx.db) do port, base; HTTP.get("$base/api/exposures/$(fx.e_id)/trace"; status_exception=false) end
            ip = with_inproc_routes(fx.db) do call; call("GET", "/api/exposures/$(fx.e_id)/trace") end
            @test w.status == 200
            _assert_equiv(w, ip; label="trace")
        end
    end

    @testset "SSE broadcast fans out identically in-process" begin
        # The wire path and in-process path must both flush the post-commit
        # broadcast to a fake subscriber. Push a (pending=Channel,) sub directly
        # (mirrors test_idempotency_replay_invariant.jl::_capture_sse_during) and
        # assert exactly one frame on each transport. NOTE: in-process flush is
        # synchronous (events.jl _flush_post_commit_broadcasts! runs in-task), so
        # no sleep is needed; the wire path needs a short drain.
        function _count_frames(seed_op)
            pending = Channel{String}(64)
            sub = (pending = pending,)
            lock(HimalayaUI.SSE_LOCK) do; push!(HimalayaUI.SSE_SUBSCRIBERS[], sub); end
            try
                seed_op()
            finally
                lock(HimalayaUI.SSE_LOCK) do
                    filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
                end
                close(pending)
            end
            count(f -> !startswith(f, ":") && occursin("post_message", f), collect(pending))
        end
        hdrs = ["Content-Type"=>"application/json", "X-Username"=>"alice"]
        bodyj = JSON3.write(Dict(:body => "hi"))
        wire_n = mktempdir() do d
            fx = _seed(d)
            _count_frames() do
                with_test_server(fx.db) do port, base
                    HTTP.post("$base/api/samples/$(fx.s_id)/messages"; body=bodyj, headers=hdrs)
                    sleep(0.3)   # wire broadcast fires off-task; allow the drain
                end
            end
        end
        inproc_n = mktempdir() do d
            fx = _seed(d)
            _count_frames() do
                with_inproc_routes(fx.db) do call
                    call("POST", "/api/samples/$(fx.s_id)/messages"; headers=hdrs, body=Vector{UInt8}(bodyj))
                end  # synchronous flush — frame already on the channel
            end
        end
        @test wire_n == 1
        @test inproc_n == 1
    end
end
```

(Optional: the DELETE row above asserts 404 parity. For genuine empty-body 204 coverage, a real 204 route exists at `routes_samples.jl:267` / `routes_exposures.jl:221` — add a row deleting a freshly-created tag per-transport-DB if 204 parity is wanted. Not load-bearing given the other rows.)

- [ ] **Step 2: Run the harness; verify it passes**

```bash
cd packages/HimalayaUI && julia --project=. -e '
  using Test
  include("test/test_inproc_equivalence.jl")' > /tmp/jl-m0_4.out 2>&1; echo "exit=$?"
grep -E "in-process ≡ wire|Test Summary|Fail|Error" /tmp/jl-m0_4.out | tail
```
Expected: `exit=0`, the `in-process ≡ wire equivalence` testset passes with 0 fails. If a row fails, that is a *real* dispatch divergence — investigate before proceeding (do NOT migrate any file until this is green).

- [ ] **Step 3: Add it to runtests.jl right after test_http.jl**

In `runtests.jl`, add `_timed_include("test_inproc_equivalence.jl")` immediately after the `test_http.jl` include (so `with_inproc_routes` is loaded). Keep all other ordering.

- [ ] **Step 4: Run the full suite; confirm green and pass count grew only by the new harness's tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-m0_4b.out 2>&1; echo "exit=$?"
grep -E "Test Summary" /tmp/jl-m0_4b.out | tail -3
```
Expected: `exit=0`; Pass == baseline + (the equivalence testset's assertion count); Fail==Error==0.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/test/test_inproc_equivalence.jl packages/HimalayaUI/test/runtests.jl
git commit -m "test(perf): differential wire≡in-process equivalence harness (M0.4)"
```

**M0 gate:** full suite green; baseline recorded; equivalence harness green. The primitive is proven — only now does migration begin.

---

## Milestone M1 — In-process dispatch migration (the ~80% win)

The mechanical swap recipe (applies to every MIGRATE file):

1. The file already has the `isdefined(@__MODULE__, :with_test_server)` self-include guard (or gains nothing — `with_inproc_routes` rides the same `test_http.jl` include, so no new include is needed).
2. For each block `with_test_server(<db-expr>) do port, base … end`:
   - Replace the header with `with_inproc_routes(<db-expr>) do call`.
   - Replace each `HTTP.get("$base/api/X")` → `call("GET", "/api/X")`.
   - `HTTP.post("$base/api/X"; body=B, headers=H)` → `call("POST", "/api/X"; headers=H, body=Vector{UInt8}(B))`.
   - `HTTP.patch(...)`/`HTTP.delete(...)` likewise (`"PATCH"`/`"DELETE"`).
   - `HTTP.request("PATCH", "$base/api/X", H, B; status_exception=false)` → `call("PATCH", "/api/X"; headers=H, body=Vector{UInt8}(B))`. (In-process `call` never throws on 4xx/5xx, so `status_exception=false` simply drops.)
   - Leave every `@test r.status == …` / `JSON3.read(String(r.body))` assertion **unchanged**.
3. Run the file standalone green-after; commit.

`<db-expr>` is `db` for 18 files and `ctx.db` for the `ctx`-fixture files — copy it verbatim from the existing call.

### Task M1.1: Pilot migration — `test_routes_samples.jl` (22 calls)

**Files:**
- Modify: `packages/HimalayaUI/test/test_routes_samples.jl`

- [ ] **Step 1: Green-before on the wire (record count)**

```bash
cd packages/HimalayaUI && julia --project=. -e 'using Test; include("test/test_http.jl"); include("test/test_routes_samples.jl")' > /tmp/jl-samples-before.out 2>&1; echo "exit=$?"
grep -E "Test Summary|Pass|Fail" /tmp/jl-samples-before.out | tail
```
Expected: `exit=0`; note the Pass count.

- [ ] **Step 2: Apply the swap recipe to all 22 blocks**

Example transformation (from `test_routes_samples.jl:11-26`):

```julia
# BEFORE
with_test_server(db) do port, base
    r = HTTP.get("$base/api/experiments/$exp_id/samples")
    @test r.status == 200
    list = JSON3.read(String(r.body))
    @test length(list) == 2
    r = HTTP.patch("$base/api/samples/$s1";
        body = JSON3.write(Dict(:notes => "hello")),
        headers = ["Content-Type" => "application/json", "X-Username" => "alice"])
    @test r.status == 200
end

# AFTER
with_inproc_routes(db) do call
    r = call("GET", "/api/experiments/$exp_id/samples")
    @test r.status == 200
    list = JSON3.read(String(r.body))
    @test length(list) == 2
    r = call("PATCH", "/api/samples/$s1";
        headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
        body = Vector{UInt8}(JSON3.write(Dict(:notes => "hello"))))
    @test r.status == 200
end
```

And the `HTTP.request` validation block (`:106-122`):

```julia
# BEFORE: HTTP.request("PATCH", "$base/api/samples/$sid", [...headers...], BODY; status_exception=false)
# AFTER:
r = call("PATCH", "/api/samples/$sid";
    headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
    body = Vector{UInt8}(JSON3.write(Dict(:name => "renamed"))))
@test r.status == 400
```

- [ ] **Step 3: Green-after in-process; same count**

```bash
cd packages/HimalayaUI && julia --project=. -e 'using Test; include("test/test_http.jl"); include("test/test_routes_samples.jl")' > /tmp/jl-samples-after.out 2>&1; echo "exit=$?"
grep -E "Test Summary|Pass|Fail" /tmp/jl-samples-after.out | tail
```
Expected: `exit=0`; Pass count **==** the Step 1 count.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/test/test_routes_samples.jl
git commit -m "test(perf): migrate test_routes_samples to in-process dispatch (M1.1)"
```

### Task M1.2: Migrate the remaining MIGRATE files (one commit each)

**Files (each its own commit, applying the M1 recipe + the green-before/after gate from M1.1):**

- [ ] `test_routes_series.jl` (16, `db`)
- [ ] `test_routes_resolve.jl` (15, `ctx.db`)
- [ ] `test_route_response_shapes.jl` (18, mixed `db`/`ctx.db` — copy each verbatim)
- [ ] `test_routes_peaks.jl` (12, `db`)
- [ ] `test_picker_routes.jl` (10, `ctx.db`)
- [ ] `test_assignment_reattach.jl` (5, `db`)
- [ ] `test_assignments.jl` (5, `db`)
- [ ] `test_routes_experiments.jl` (5, `db`)
- [ ] `test_routes_image.jl` (3, `db`)
- [ ] `test_routes_analysis.jl` (3, `db`)
- [ ] `test_speculative.jl` (3, `db`)
- [ ] `test_picker_samples_route.jl` (3, `db`)
- [ ] `test_routes_exposures.jl` (2, `db`)
- [ ] `test_routes_messages.jl` (2, `db`)
- [ ] `test_routes_trace.jl` (2, `db`)
- [ ] `test_route_validation_routing.jl` (1, `ctx.db`)
- [ ] `test_routes_export.jl` (1, `db`)
- [ ] `test_routes_mentions.jl` (1, `db`)
- [ ] `test_routes_users.jl` (1, `db`)

For each: green-before (`include("test/test_http.jl"); include("test/<file>")`), apply recipe, green-after with equal Pass count, commit `test(perf): migrate <file> to in-process dispatch (M1.2)`.

- [ ] **Special-case `test_idempotency_replay_invariant.jl` (4 calls, `db`):** migrate the 4 `with_test_server` blocks AND remove the now-dead `sleep(0.3)` (the in-process `call` flushes the post-commit broadcast synchronously, `events.jl` `_flush_post_commit_broadcasts!`). In `_capture_sse_during` (`:41-65`), delete the `sleep(0.3)` line — the frames are already on the channel by the time `f()` returns in-process. Green-after must still show `length(frames) == 1` for both replay testsets. Commit separately: `test(perf): migrate idempotency replay invariant + drop dead sleep (M1.2)`.

### Task M1.3: M1 milestone gate

- [ ] **Step 1: Confirm keepers untouched**

Verify these still use the wire (no change): `test_sse.jl`, `test_routes_sse_broadcast.jl`, `test_spa_fallback.jl`, `test_routes_health.jl`, `test_routes_status.jl`, `test_idempotency_sse_suppression.jl`. And `test_concurrent_writes.jl` (never used `with_test_server`).

```bash
grep -l "with_test_server" packages/HimalayaUI/test/*.jl
```
Expected: only the keepers above + `test_http.jl` (definition) + `test_inproc_equivalence.jl` (uses it deliberately).

- [ ] **Step 2: Full suite, captured once**

```bash
/usr/bin/time -p julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-m1-gate.out 2>&1; echo "exit=$?"
grep -E "Test Summary|^real " /tmp/jl-m1-gate.out | tail
```
Expected: `exit=0`; Pass **==** baseline (+ the M0.4 harness assertions); Fail==Error==0.

- [ ] **Step 3: Go/no-go threshold**

Read the `real` wall time. **It must be < 600 s (10 min); target ≈300–480 s.** If it is not under 10 min, STOP and revisit the keeper boundary (are too many heavy files still on the wire?) before starting M2. Record the new wall time in `BASELINE.md` under an "M1" heading.

```bash
git add packages/HimalayaUI/test/BASELINE.md
git commit -m "test(perf): record M1 wall time (M1.3 gate)"
```

---

## Milestone M2 — Template DB for the `open_db` family

Scope: only the **FK-ON `open_db`-backed** files. The `:memory:` `SQLite.DB()`+`create_schema!` files (FK OFF) and the `test_db.jl` legacy-migration testsets are left as-is.

### Task M2.1: Template builder + `open_prepared_clone` helper

**Files:**
- Create: `packages/HimalayaUI/test/test_template_db.jl`
- Modify: `packages/HimalayaUI/test/test_http.jl` (include it, so it loads before fixtures)

**Interfaces:**
- Produces:
  - `prepared_template_path() -> String` — path to a process-wide, built-once, fully-migrated, WAL-checkpointed template DB file (zero rows).
  - `open_prepared_clone(dir) -> SQLite.DB` — copies the template into `dir`, opens it, re-applies `PRAGMA foreign_keys=ON` + WAL + `SQLite.finalize_statements!`, returns the connection. Drop-in for `open_db(joinpath(dir, "h.db"))`.
- Consumes: `HimalayaUI.open_db`.

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/_scratch_template_test.jl`:

```julia
using Test, SQLite, DBInterface, Tables
using HimalayaUI
include("test_template_db.jl")

@testset "open_prepared_clone parity with open_db" begin
    mktempdir() do d
        db = open_prepared_clone(d)
        # FK enforcement restored (the whole point):
        fk = Tables.rowtable(DBInterface.execute(db, "PRAGMA foreign_keys"))
        @test fk[1][1] == 1
        # Schema present (a migrated column exists):
        cols = [r.name for r in Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(experiments)"))]
        @test "config" in cols
        # Zero-row template ⇒ first insert id == 1 (sqlite_sequence empty):
        eid = HimalayaUI.create_experiment!(db; name="T", path="/p",
            data_dir="/p/d", analysis_dir="/p/a")
        @test eid == 1
        # FK violation actually throws (proves foreign_keys=ON took effect):
        @test_throws Exception DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, name) VALUES (99999, 'x')")
    end
end
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd packages/HimalayaUI && julia --project=. test/_scratch_template_test.jl 2>&1 | tail -20
```
Expected: FAIL/ERROR — `test_template_db.jl` / `open_prepared_clone` not defined.

- [ ] **Step 3: Write the helper**

Create `packages/HimalayaUI/test/test_template_db.jl`:

```julia
# M2: build a fully-migrated, zero-row template DB ONCE per process, then clone
# it per test instead of re-running create_schema!+migrate_schema! each time.
# Scoped to the FK-ON open_db family only. VACUUM INTO produces a compacted,
# WAL-free single file; the clone re-applies the connection-scoped setup that
# open_db does (foreign_keys=ON, WAL, finalize_statements!) — none of which
# survive a file copy.
using SQLite, DBInterface, Tables
using HimalayaUI

const _TEMPLATE_REF = Ref{Union{String,Nothing}}(nothing)

"Build (once) and return the path to the migrated, checkpointed template DB."
function prepared_template_path()
    _TEMPLATE_REF[] === nothing || return _TEMPLATE_REF[]
    tmpdir = mktempdir(; cleanup = false)         # lives for the process
    src = joinpath(tmpdir, "source.db")
    db = HimalayaUI.open_db(src)                  # full create_schema! + migrate_schema!
    # Checkpoint the WAL into the main file and close, so VACUUM INTO copies a
    # quiescent single file with no -wal sidecar.
    DBInterface.execute(db, "PRAGMA wal_checkpoint(TRUNCATE)")
    tmpl = joinpath(tmpdir, "template.db")
    DBInterface.execute(db, "VACUUM INTO '$(tmpl)'")   # target must NOT pre-exist
    SQLite.close(db)
    _TEMPLATE_REF[] = tmpl
    return tmpl
end

"Clone the template into `dir` and return a connection equivalent to open_db's."
function open_prepared_clone(dir::AbstractString)
    tmpl = prepared_template_path()
    dest = joinpath(dir, "h.db")
    cp(tmpl, dest; force = false)                 # copy ONLY the main file
    db = SQLite.DB(dest)
    SQLite.finalize_statements!(db)
    DBInterface.execute(db, "PRAGMA foreign_keys = ON")
    Tables.rowtable(DBInterface.execute(db, "PRAGMA journal_mode = WAL"))
    SQLite.finalize_statements!(db)
    return db
end
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
cd packages/HimalayaUI && julia --project=. test/_scratch_template_test.jl 2>&1 | tail -20
```
Expected: PASS (FK==1, `config` column present, first id==1, FK violation throws).

- [ ] **Step 5: Wire into load order, remove scratch**

Add to the END of `test/test_http.jl` (after the `test_inproc.jl` include):

```julia
include("test_template_db.jl")
```
```bash
rm packages/HimalayaUI/test/_scratch_template_test.jl
```

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/test/test_template_db.jl packages/HimalayaUI/test/test_http.jl
git commit -m "test(perf): template-DB clone helper for open_db family (M2.1)"
```

### Task M2.2: Swap `open_db` fixtures to `open_prepared_clone` (regime-b files)

The mechanical swap: replace `HimalayaUI.open_db(joinpath(<dir>, "h.db"))` (and `"himalaya.db"`) → `open_prepared_clone(<dir>)` **only** where the fixture is a fresh empty DB the test then populates. Do **not** touch `:memory:` `SQLite.DB()` fixtures, and do **not** touch `test_db.jl`'s legacy-migration testsets (they need `migrate_schema!` to run on a legacy shape).

- [ ] **Step 1: For `test_db.jl`, isolate the legacy-migration testsets behind an opt-out**

These testsets (verbatim list — they call `migrate_schema!`/`open_db` on a legacy-shaped DB and must NOT use the clone): the testsets at `test_db.jl` lines ~178-222, 233-286, 318-386, 388-418, 420-465, 507-537, 563-588, 590-602, 629-642, 644-717, 1029-1134, 1136-1163. Leave all of these exactly as they are (they build `SQLite.DB()` / synthetic legacy schemas directly). Only the *CRUD/behavior* testsets that build a fresh `open_db` fixture are candidates — and `test_db.jl`'s first-insert `id==1` CRUD testset (`:56-84`) uses `SQLite.DB()` + `create_schema!` directly (regime a), so it is **not** swapped either.

- [ ] **Step 2: Find the actual swap candidates (don't trust a static list)**

The candidate set is exactly the files with a fresh-fixture `open_db(joinpath(...))` call — derive it, don't guess:

```bash
cd packages/HimalayaUI/test && grep -cE "open_db\(joinpath" *.jl | grep -v ':0'
```
Swap **only** those `open_db(joinpath(<dir>, "...db"))` calls → `open_prepared_clone(<dir>)`. Notes from the review (apply verbatim):
- **Exclude (no swappable fresh fixture / wrong unit-under-test):** `test_routes_samples.jl`, `test_routes_resolve.jl` (no own `open_db(joinpath)` — they consume the shared `_setup_*` fixtures), and `test_config.jl:599` (its `open_db` *is* the parent-dir-creation unit-under-test — swapping it silently drops `mkpath` coverage and the green gate won't catch it).
- **DO-NOT-SWAP:** `test_migrate_comparisons_to_series.jl` lines 232/248 (the "gated no-op" idempotency testset relies on `migrate_schema!` running on the opened DB — cloning defeats it). Annotate like the `test_db.jl` legacy testsets.
- **Mixed-regime files** (`test_image.jl`, `test_config.jl`, `test_pipeline.jl`): swap ONLY the `open_db(joinpath)` fresh fixtures; leave every bare `SQLite.DB()` (regime-a) call as-is.

- [ ] **Step 3: Swap one commit each (green-before/green-after, equal Pass count)**

For each candidate: green-before (`include test_http.jl; include <file>`), swap, green-after with equal Pass count, commit `test(perf): clone template DB in <file> (M2.2)`.

- [ ] **Special case — the shared fixture (`_setup_analyzed_exposure`/`_setup_for_resolve`) opens its DB once in `test_route_response_shapes.jl:36`.** Swapping that single `open_db` changes the fixture for **three** consumers (`test_route_response_shapes`, `test_routes_resolve`, `test_picker_routes`) that have no own `open_db`. The per-file green-after for `test_route_response_shapes` alone does NOT exercise the other two — run all three consumer files together in the green-after for that commit.

For any file where a fixture asserts a specific autoincrement id > the template's zero baseline OR depends on `migrate_schema!` running, leave that fixture on `open_db` and note it in the commit body. (The green-after gate catches id mismatches; it does NOT catch dropped `mkpath`/migration coverage — hence the explicit exclusions above.)

### Task M2.3: M2 milestone gate

- [ ] **Step 1: Full suite, captured once**

```bash
/usr/bin/time -p julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-m2-gate.out 2>&1; echo "exit=$?"
grep -E "Test Summary|^real " /tmp/jl-m2-gate.out | tail
```
Expected: `exit=0`; Pass **==** the M1 gate count (+0 — template is a fixture swap, not new tests); Fail==Error==0. Record the new wall time in `BASELINE.md` under "M2".

- [ ] **Step 2: Commit the wall-time record**

```bash
git add packages/HimalayaUI/test/BASELINE.md
git commit -m "test(perf): record M2 wall time (M2.3 gate)"
```

---

## Milestone M3 — Named GROUP buckets + parallel runner

> **Blocker B-1 — read first.** The flat serial `runtests.jl` silently lets a helper defined in file A be used by file B (A included earlier). Partitioning into buckets severs that, so co-dependent files in *different* buckets throw `UndefVarError` at include time. Confirmed cross-file helpers: `_setup_analyzed_exposure` (defined **twice, incompatibly** — `test_route_response_shapes.jl:31` 4-field NT, `test_comparisons.jl:59` 5-field-with-`experiment_id`), `_premint_cmp!` (`test_comparisons.jl:19`), `_member_payload` (`test_comparisons.jl:26`), `_setup_for_resolve` (`test_route_response_shapes.jl:50`); consumers span `test_picker_routes`, `test_routes_resolve`, `test_spa_fallback`, `test_routes_series`, `test_comparison_pins`, `test_events`. **Task M3.0 must extract these into a shared file before bucketing**, or `GROUP=routes`/`GROUP=wire` are guaranteed red.

### Task M3.0: Extract cross-file test helpers into a shared fixtures module

**Files:**
- Create: `packages/HimalayaUI/test/test_fixtures.jl`
- Modify: `packages/HimalayaUI/test/test_http.jl` (include it, so every bucket gets it), `test_comparisons.jl`, `test_route_response_shapes.jl` (remove the moved defs; the canonical ones now live in `test_fixtures.jl`).

**Interfaces:**
- Produces: a single definition of each shared helper, included everywhere `test_http.jl` is.

- [ ] **Step 1: Audit every cross-file helper**

```bash
cd packages/HimalayaUI/test
# For each top-level `function _foo(` find which file defines it and which OTHER files use it.
for fn in $(grep -rhoE '^function (_[A-Za-z0-9_!]+)\(' *.jl | sed -E 's/^function (_[A-Za-z0-9_!]+)\(/\1/' | sort -u); do
  defs=$(grep -rl "^function $fn(" *.jl | tr '\n' ' ')
  uses=$(grep -rl "\b$fn\b" *.jl | tr '\n' ' ')
  ndef=$(echo $defs | wc -w); nuse=$(echo $uses | wc -w)
  # flag helpers defined in 1 file but used in >1 file (cross-file), or defined in >1 file (dup)
  if [ "$ndef" -gt 1 ] || { [ "$ndef" -eq 1 ] && [ "$nuse" -gt 1 ]; }; then
    echo "CROSS/DUP: $fn  defs=[$defs] uses=[$uses]"
  fi
done
```
Expected: at minimum `_setup_analyzed_exposure` (DUP — two defs), `_premint_cmp!`, `_member_payload`, `_setup_for_resolve`. Capture the full list — it drives what moves into `test_fixtures.jl`.

- [ ] **Step 2: Resolve the duplicate `_setup_analyzed_exposure`**

Under serial `GROUP=All`, `test_comparisons.jl` (include #44) redefines it *before* `test_picker_routes.jl`/`test_routes_resolve.jl` (#46/#48) consume it, so those consumers actually use the **5-field** (`experiment_id`-bearing) `test_comparisons.jl:59` version. Make `test_fixtures.jl` hold that **5-field** definition (the superset). Read both definitions; confirm the 4-field consumers (`test_route_response_shapes.jl`, `test_spa_fallback.jl`) still work against the 5-field shape (a named-tuple with an *extra* field is safe for callers that ignore it). If any 4-field consumer destructures positionally or asserts the NT's exact fields, keep its local call working — adjust the consumer, not the canonical shape.

- [ ] **Step 3: Create `test_fixtures.jl` with the canonical helpers**

Move the canonical `_setup_analyzed_exposure` (5-field), `_premint_cmp!`, `_member_payload`, `_setup_for_resolve`, and any other CROSS/DUP helper from Step 1 into `packages/HimalayaUI/test/test_fixtures.jl` (verbatim bodies, single definition each). Delete the now-moved definitions from `test_comparisons.jl` and `test_route_response_shapes.jl`.

- [ ] **Step 4: Include it from test_http.jl, before the others**

At the TOP of `test/test_http.jl` body (before the `test_inproc.jl`/`test_template_db.jl` includes):

```julia
if !isdefined(@__MODULE__, :_setup_analyzed_exposure)
    include("test_fixtures.jl")
end
```

- [ ] **Step 5: Green check — full suite still passes (GROUP=All), count unchanged**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-m3_0.out 2>&1; echo "exit=$?"
grep -E "Test Summary" /tmp/jl-m3_0.out | tail -3
```
Expected: `exit=0`; Pass == M2 gate count (pure refactor — same tests, deduped helpers).

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/test/test_fixtures.jl packages/HimalayaUI/test/test_http.jl packages/HimalayaUI/test/test_comparisons.jl packages/HimalayaUI/test/test_route_response_shapes.jl
git commit -m "test(perf): extract shared cross-file test fixtures (M3.0, unblocks bucketing)"
```

### Task M3.1: Carve runtests.jl into GROUP buckets

**Files:**
- Modify: `packages/HimalayaUI/test/runtests.jl`

**Interfaces:**
- Produces: `GROUP=<name>` selects a bucket; `GROUP=All` (default) runs every file in the **exact current order**.

- [ ] **Step 1: Restructure with a GROUP dispatch, preserving order for `All`**

```julia
using Test

const GROUP = get(ENV, "GROUP", "All")
const _TEST_TIMES = Dict{String,Float64}()
_timed_include(f) = (_TEST_TIMES[f] = @elapsed include(f))
_want(g) = GROUP == "All" || GROUP == g

# Buckets balanced from M0.2 timing (heaviest files spread across buckets).
# IMPORTANT: under GROUP=All these run in the SAME total order as before
# (db → pipeline → routes → events → wire → …). Each bucket includes
# test_http.jl first (the isdefined guard makes re-include a no-op).
# Every bucket lists test_http.jl FIRST — it transitively includes test_fixtures.jl
# (M3.0), test_inproc.jl, test_template_db.jl, so any cross-file helper is in scope
# regardless of which bucket a file lands in. (test_http.jl defines no tests, so the
# Pass-sum invariant across buckets holds. The seen-set dedup keeps it included once
# per process.)
const GROUPS = [
    ("db",       ["test_http.jl","test_config.jl","test_db.jl","test_migrate_comparisons_to_series.jl",
                  "test_manifest.jl","test_migrate_toml.jl","test_validate.jl"]),
    ("pipeline", ["test_http.jl","test_datfile.jl","test_hash.jl","test_hash_peak_set_memoization.jl",
                  "test_pipeline.jl","test_auto_group_peak_id_claiming.jl",
                  "test_effective_peaks_sharpness_passthrough.jl","test_fast_skip.jl",
                  "test_json.jl","test_image.jl"]),
    ("routes",   ["test_http.jl","test_routes_users.jl","test_routes_experiments.jl",
                  "test_routes_samples.jl","test_routes_exposures.jl","test_routes_image.jl",
                  "test_routes_peaks.jl","test_routes_messages.jl","test_routes_trace.jl",
                  "test_routes_analysis.jl","test_speculative.jl","test_routes_export.jl",
                  "test_routes_mentions.jl","test_route_response_shapes.jl",
                  "test_route_validation_routing.jl","test_routes_series.jl",
                  "test_picker_routes.jl","test_picker_samples_route.jl","test_routes_resolve.jl",
                  "test_inproc_equivalence.jl"]),
    ("events",   ["test_http.jl","test_actions.jl","test_events.jl","test_assignments.jl",
                  "test_assignment_reattach.jl","test_idempotency.jl",
                  "test_idempotency_replay_invariant.jl","test_concurrent_writes.jl",
                  "test_idempotency_sse_suppression.jl","test_comparisons.jl",
                  "test_comparison_pins.jl"]),
    ("wire",     ["test_http.jl","test_routes_health.jl","test_routes_status.jl",
                  "test_sse.jl","test_routes_sse_broadcast.jl","test_spa_fallback.jl"]),
]

# GROUP=All must reproduce the historical order exactly — assert coverage.
@testset "HimalayaUI" begin
    if GROUP == "All"
        for f in ["test_config.jl","test_db.jl","test_migrate_comparisons_to_series.jl",
                  "test_datfile.jl","test_manifest.jl","test_hash.jl",
                  "test_hash_peak_set_memoization.jl","test_pipeline.jl",
                  "test_auto_group_peak_id_claiming.jl",
                  "test_effective_peaks_sharpness_passthrough.jl","test_fast_skip.jl",
                  "test_json.jl","test_http.jl","test_inproc_equivalence.jl",
                  "test_routes_health.jl","test_routes_users.jl","test_routes_experiments.jl",
                  "test_routes_samples.jl","test_routes_exposures.jl","test_image.jl",
                  "test_routes_image.jl","test_routes_status.jl","test_routes_peaks.jl",
                  "test_routes_messages.jl","test_routes_trace.jl","test_routes_analysis.jl",
                  "test_speculative.jl","test_routes_export.jl","test_routes_mentions.jl",
                  "test_actions.jl","test_events.jl","test_assignments.jl",
                  "test_assignment_reattach.jl","test_sse.jl","test_routes_sse_broadcast.jl",
                  "test_route_response_shapes.jl","test_route_validation_routing.jl",
                  "test_idempotency.jl","test_idempotency_replay_invariant.jl",
                  "test_concurrent_writes.jl","test_idempotency_sse_suppression.jl",
                  "test_comparisons.jl","test_routes_series.jl","test_picker_routes.jl",
                  "test_picker_samples_route.jl","test_routes_resolve.jl",
                  "test_comparison_pins.jl","test_validate.jl","test_migrate_toml.jl",
                  "test_spa_fallback.jl"]
            _timed_include(f)
        end
    else
        seen = Set{String}()
        for (name, files) in GROUPS
            _want(name) || continue
            for f in files
                f in seen && continue   # test_http.jl appears in several buckets
                push!(seen, f)
                _timed_include(f)
            end
        end
    end
end

let rows = sort(collect(_TEST_TIMES); by = last, rev = true)
    println("\n── per-file test timing (slowest first) ──")
    for (f, t) in rows
        println(rpad(f, 42), lpad(round(t; digits=1), 8), " s")
    end
end
```

(Note: `test_inproc_equivalence.jl` includes `test_http.jl` via its own guard, and every bucket lists `test_http.jl` first so the `isdefined` guard makes any in-bucket route file self-sufficient.)

- [ ] **Step 2: GROUP=All reproduces baseline**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-m3-all.out 2>&1; echo "exit=$?"
grep -E "Test Summary" /tmp/jl-m3-all.out | tail -3
```
Expected: `exit=0`; Pass **==** M2 gate count.

- [ ] **Step 3: Each bucket green standalone**

```bash
for g in db pipeline routes events wire; do
  echo "== GROUP=$g =="
  GROUP=$g julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-m3-$g.out 2>&1; echo "exit=$?"
  grep -E "Test Summary" /tmp/jl-m3-$g.out | tail -1
done
```
Expected: each `exit=0`. (Pkg.test passes `GROUP` through the env to the test process.) The five buckets' Pass counts should sum to the `All` count.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/test/runtests.jl
git commit -m "test(perf): named GROUP buckets with GROUP=All serial fallback (M3.1)"
```

### Task M3.2: Parallel runner make-target

**Files:**
- Modify: `Makefile`

- [ ] **Step 1: Add a `test-parallel` target**

Append to `Makefile`:

```makefile
.PHONY: test-parallel
GROUPS := db pipeline routes events wire
test-parallel:
	@echo "Running $(words $(GROUPS)) groups in parallel..."
	@pids=""; rc=0; \
	for g in $(GROUPS); do \
		( GROUP=$$g julia --project=packages/HimalayaUI \
			-e 'using Pkg; Pkg.test("HimalayaUI")' > build/test-$$g.log 2>&1 ) & \
		pids="$$pids $$!"; \
	done; \
	for p in $$pids; do wait $$p || rc=1; done; \
	for g in $(GROUPS); do \
		echo "== $$g =="; grep -E "Test Summary" build/test-$$g.log | tail -1 || true; \
	done; \
	exit $$rc
```

(Optional, deferred: add `--sysimage build/himalaya.so` to each `julia` invocation once a route-dispatch precompile workload is added — the current app sysimage doesn't cover `internalrequest`/`Test`, so re-measure before enabling by default.)

- [ ] **Step 2: Run it; confirm all groups pass and wall time drops**

```bash
mkdir -p build
/usr/bin/time -p make test-parallel > /tmp/jl-m3-parallel.out 2>&1; echo "exit=$?"
grep -E "Test Summary|^real " /tmp/jl-m3-parallel.out | tail
```
Expected: `exit=0`; each group's summary shows 0 fails; `real` wall time is bounded by the slowest single group (target: well under the M2 serial time).

- [ ] **Step 3: Record final wall times + commit**

Add an "M3" section to `BASELINE.md` with the serial (`GROUP=All`) and parallel wall times.

```bash
git add Makefile packages/HimalayaUI/test/BASELINE.md
git commit -m "test(perf): make test-parallel runner over GROUP buckets (M3.2)"
```

### Task M3.3: M3 milestone gate

- [ ] **Step 0: Guard — no `@testset` in the per-bucket-included helper files**

```bash
cd packages/HimalayaUI && grep -l '@testset' test/test_http.jl test/test_inproc.jl test/test_template_db.jl test/test_fixtures.jl 2>/dev/null
```
Expected: NO output. Any hit means that file's tests run once per bucket (3×) and the Pass-sum across buckets won't equal `GROUP=All` — move the test out before proceeding.

- [ ] **Step 1: Final verification — serial fallback + parallel both green, count matches baseline**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-final.out 2>&1; echo "serial exit=$?"
grep -E "Test Summary" /tmp/jl-final.out | tail -3
make test-parallel > /tmp/jl-final-par.out 2>&1; echo "parallel exit=$?"
```
Expected: both `exit=0`; serial Pass **==** baseline (+M0.4 harness); parallel groups all green. Done.

---

## Self-Review

**Spec coverage:** §4 fixture seam → M0.3. §5 equivalence + keepers → M0.4 + M1.1/M1.3. §6 M0–M3 → Milestones M0–M3. §7 isolation hazards → Global Constraints + the `finally` scrub in M0.3/M2.1 + `HIMALAYA_FRONTEND_DIST` (covered by keeping `test_spa_fallback` a wire keeper; the in-process probe of `/api/health` doesn't depend on `dist_dir`). §8 don'ts → Global Constraints (no rollback, no ReTestItems/ParallelTestRunner, no numeric shards, preserve include order). §9 gates → M2.1 (FK/WAL), M0.4 (header allowlist + per-transport DB), M1.2 keeper classification.

**Note on `HIMALAYA_FRONTEND_DIST`:** the spec asks to neutralize it before in-process registration. Because the in-process probe only ever hits `/api/health` (always mounted) and migrated tests never assert the SPA catch-all (that stays a wire keeper), the dist mount being present or absent does not affect migrated tests. If a migrated test ever 404s unexpectedly on a non-`api` path, set `ENV["HIMALAYA_FRONTEND_DIST"]` to a known-absent path at the top of `test_inproc.jl` before first registration. Left as a watch-point, not a forced edit.

**Placeholder scan:** no TBD/TODO; every code step shows real code; the M1.2/M2.2 file lists are explicit (recipe defined once in M1's intro + M2.2 Step 2, applied per enumerated file — mechanical and identical, not a cross-reference to distinct code).

**Type consistency:** `with_inproc_routes(f, db)` → `call(method, target; headers, body)` used identically in M0.4, M1.1, M1.2. `open_prepared_clone(dir)` returns `SQLite.DB`, used as a drop-in for `open_db(...)`. `prepared_template_path()` returns `String`.
