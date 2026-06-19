# HimalayaUI test-suite performance — design spec

**Date:** 2026-06-19
**Status:** Design approved; implementation pending
**Scope:** The HimalayaUI Julia backend test suite (`packages/HimalayaUI/test/`). Core `Himalaya` tests and the frontend Vitest/Playwright suites are out of scope.

## 1. Problem & evidence

The HimalayaUI Julia suite takes **~20 minutes** (measured ~1247s wall). A fan-out investigation + an empirical timed run established that the cost is **per-test fixture/server spin-up, not Julia compile**:

- **~70–80% (~870–1000s)** is the Oxygen HTTP server lifecycle, booted **145 times** across 25 of 50 files via `with_test_server` (`test/test_http.jl:4-12` → `start_test_server!`, `src/server.jl:197-208`). Each call runs `Oxygen.resetstate() → bind_db!(db) → register_routes!() → Oxygen.serve(async=true, parallel=true) → _wait_for_server`, then a full teardown (`stop_test_server!`, `:210-218`).
- Only **one line** varies per test: `bind_db!(db)`. `current_db()` → `_DB_REF[]` (`server.jl:38-47`) is the **sole** DB accessor; the router, port, SSE list, GC timer, and `OP_LOCKS` are all DB-independent singletons.
- `_wait_for_server` (`server.jl:227-`) health-polls `GET /api/health` with `sleep(0.05)` up to a 5s deadline — incurred 145×, the variable tail inflating the average.
- Compile/load (TTFX) is paid **once per process** (~30–90s) — mathematically cannot be 80% of a 20-minute suite.
- Secondary costs (an order of magnitude smaller): the `create_schema!` migration chain (`db.jl:234`) runs ~131× (~2–5s total); `analyze_exposure!` peak-finding runs ~43× in `test_pipeline.jl` (~100–200ms each).

**The suite cannot parallelize today** because Oxygen state + `_DB_REF` are process-global singletons: two servers up in one process corrupt each other.

## 2. Goals & non-goals

**Goals**
- Cut wall time dramatically (target: ~20 min → **< 3 min**) without losing coverage.
- Preserve test isolation and the production-fidelity coverage the suite was built with.
- No new runtime dependencies. Test-time additions limited to what already ships (`SafeTestsets`-style isolation is optional; not required).

**Non-goals**
- No rewrite into `@testitem`/ReTestItems (see §8). No `src/` behavior changes — this is a test-layer refactor.
- Not changing the production threading model or the singleton design itself.

## 3. The four parts (ecosystem-validated)

All four were pressure-tested against how large Julia packages test (Oxygen.jl, HTTP.jl, SQLite.jl, SciML, JuMP, JuliaLang `Base`).

| Part | Verdict | Note |
|---|---|---|
| **#1 In-process dispatch** | ✅ Confirmed — Oxygen's own suite uses `internalrequest` as its default route-test path; real `serve()` is reserved for transport tests (`ssetests.jl`, `streamingtests.jl`). | Middleware-parity gate **closed**: HimalayaUI passes **no** global middleware to `Oxygen.serve` (`server.jl:193,204`); cross-cutting concerns (`with_idempotency`, event logging) are **handler-body wrappers** (`routes_exposures.jl:123`), so `internalrequest`'s default `middleware=[]` drops nothing. **Two corrections from review (see §4):** dispatch with `serialize=true` (preserves throw→500; **B2**), and re-register routes via a **router-liveness check**, not a sticky boolean (`resetstate()` wipes the router; **B1**). |
| **#2 Wire-test tier** | ✅ Confirmed — mirrors Oxygen's `internalrequest`/`ssetests.jl` split. | Conservative keeper set, verified per-file. |
| **#3 Template DB** | ✅ Confirmed — SQLite.jl's own suite copies a prepared fixture per test (`setup_clean_test_db`). | `VACUUM INTO` variant is strictly better (no binary in git, matches migrations, compacted). Rollback-isolation explicitly **rejected**. |
| **#4 Named GROUP sharding** | ✅ Confirmed model (process-level), with amended interface. | SciML named-`GROUP` buckets, **not** numeric shard indices; **not** ReTestItems. |

## 4. Section 1 — the fixture seam (`with_inproc_routes`)

A new helper in `test/test_http.jl`, alongside the untouched `with_test_server`:

```julia
with_inproc_routes(db) do call
    r = call("GET", "/api/samples")          # returns the same HTTP.Response handlers build
    @test r.status == 200
end
```

- **Route registration — a *liveness check*, NOT a sticky boolean (blocker B1).** `with_inproc_routes` must `HimalayaUI.register_routes!()` whenever Oxygen's router is empty/unregistered, *not* trust a once-per-process `Ref{Bool}`. Reason: `stop_test_server!` **and** `start_test_server!` both call `Oxygen.resetstate()`, which replaces `CONTEXT[]` with a fresh empty-router `ServerContext()` (Oxygen `methods.jl:8`). Under `GROUP=All`, a wire keeper that boots+tears down **after** an in-process file ran would wipe the router; a sticky boolean would still read "registered," and the next `call()` would dispatch against an empty router and return a **silent `HTTP.Response(404)`** (not an error), reddening `@test r.status == 200`. The liveness check (re-register if the router is empty) makes any interleaving safe.
- **Per call:** `HimalayaUI.bind_db!(db)`, then `call(method, target; headers, body)` constructs an `HTTP.Request` and returns `Oxygen.internalrequest(req; serialize=true, catch_errors=true)`. **`serialize=true` is load-bearing (blocker B2)** — `setupmiddleware` only includes `DefaultSerializer` (the layer that catches a thrown handler error and returns a 500) when `serialize=true` (Oxygen `core.jl`). With `serialize=false`, a throwing handler (e.g. the bare `error(...)` at `routes_peaks.jl:167`, or an uncaught FK `SQLiteException`) would propagate **out of `call()`** instead of returning a 4xx/5xx Response, breaking tests that assert `r.status >= 400` (`test_routes_messages.jl:80`, `test_routes_samples.jl:279`). `serialize=true`'s `format_response!` passes explicit `HTTP.Response` objects through unchanged *and* wraps throws into 500, matching production `serve()` exactly. (The earlier "`serialize=false` because handlers already return Response" rationale was wrong: not all do — `/api/health` returns a bare `Dict`, and several routes `throw`.)
- **`finally` scrub** (mirrors `stop_test_server!` exactly): `_DB_REF[] = nothing`; `SSE_SUBSCRIBERS[] = []` under `SSE_LOCK`; `empty!(OP_LOCKS)`; **GC timer not started**. (Do **not** call `resetstate()` in the in-process `finally` — that would wipe the router the next in-process file relies on; the liveness check above is what restores it after a wire keeper resets it.)

Migration is a mechanical per-file swap: `with_test_server(db) do port, base … HTTP.get("$base/api/…")` → `with_inproc_routes(db) do call … call("GET", "/api/…")`. The closure returns the identical `HTTP.Response`, so every existing `r.status` / `JSON3.read(String(r.body))` assertion is **unchanged** — keeping the diff reviewable line-by-line. (Invariant verified per-migrated-route, not assumed: each migrated handler body must terminate in an `HTTP.Response`; only `/api/health` returns a bare `Dict` today, and it is a wire keeper. An M1 pre-migration grep guards this.)

**What in-process dispatch deliberately does NOT exercise:** the TCP socket, byte-level HTTP parse, and Oxygen's `async`/`parallel=true` request scheduling. Those belong to the wire tier (§5).

## 5. Section 2 — equivalence guarantee + keeper boundary

Three decoupled safety layers, so no single mechanism carries the whole burden:

1. **Differential equivalence harness** (`test_inproc_equivalence.jl`, built first, kept permanently). This is the **single load-bearing equivalence proof** (layers 2/3 below are migration-correctness + count nets, not independent equivalence proofs), so its matrix must be exhaustive over the response shapes most likely to diverge in `internalrequest`'s request construction. Fires, through **both** a real `with_test_server` socket **and** `with_inproc_routes`:
   - `GET` `200` (JSON) · `POST`/`PATCH` `200` · `DELETE` `204` empty-body.
   - a **query-string GET** (verify `call()` preserves the query in `HTTP.Request.target`; `routes_exposures.jl:48` uses `HTTP.queryparams`).
   - a **numeric-array JSON body** (trace route).
   - a **binary/PNG image** response (asserts `Content-Type: image/png` + the `X-Image-*` contract headers; `routes_exposures.jl:47-60`).
   - a body-validation **4xx** AND an **uncaught-throw→500** row driven by a real throwing handler (e.g. an FK violation), asserting throw→same-status on both transports (guards blocker B2).
   - an **idempotency-replay** pair — using a **fresh clone DB and a distinct `client_op_id` per transport**, else the second transport replays the first's cached row (`idempotency.jl:108-110`) and masks a real difference.
   - an **SSE-broadcast-then-drain**.

   Asserts equivalent `HTTP.Response`: identical status and body, and headers compared by **allowlist** — drop only a fixed transport-only set (`Date`, `Server`, `Transfer-Encoding`, `Connection`, `Content-Length`) and assert **all other** headers match, so contract headers (`Content-Type`, `Cache-Control`, `X-Image-*`) cannot silently slip (see §9). Proves the primitive once; guards against future Oxygen drift. *Answers "is dispatch faithful."*
2. **Per-file green-before/green-after gate.** Each file migrated as its own commit: green-on-wire → mechanical swap → green-in-process, same test count. *Answers "did I migrate this file correctly."* Decoupled from layer 1: a red file with a green harness ⇒ a migration slip, not a dispatch problem.
3. **Full-suite milestone gate.** The slow suite (`> /tmp/jl-test.out`, capture-once per `test/AGENTS.md`) runs at each milestone boundary, asserting pass count **==** the M0 baseline.

**Keeper boundary (conservative, verified per-file against source).** A test stays on `with_test_server` only if it depends on **real transport** (the wire / streaming / static-file serving). Note: **no current test asserts real `parallel=true` HTTP concurrency** — the "concurrency" tests drive `with_idempotency` directly via `@async`/`Threads.@spawn`, not the wire — so the keeper set is narrow. Confirm each per-file:
- SSE over the wire / streaming: `test_sse.jl`, `test_routes_sse_broadcast.jl`, `test_idempotency_sse_suppression.jl` (confirm each — broadcast tests that only push a fake subscriber and drain a `Channel` can move in-process and do NOT need the wire).
- Static-file / env: `test_spa_fallback.jl` (mutates `HIMALAYA_FRONTEND_DIST`, serves `index.html` over the wire), `test_routes_health.jl`, `test_routes_status.jl`.

**Not a wire keeper — scrub-only:** `test_concurrent_writes.jl` never opens a socket (it calls `with_idempotency` directly under `@async`/`Threads.@spawn`). It moves in-process but **must** get the full `_DB_REF`/`OP_LOCKS` `finally`-scrub; do **not** bucket it into the M3 `wire` group.

Wire-tier files boot the server **once per file** (Oxygen `ssetests.jl` pattern), keep the ephemeral `find_free_port()`, and wrap teardown under a timeout so a hung server can't wedge the group.

## 6. Section 3 — milestone sequencing

Four independently landable, independently revertable milestones, with a green gate between each.

**M0 — Baseline + harness** *(no behavior change)*
- Capture the full-suite green baseline (exact pass count + wall time) = the regression floor.
- Add per-`@testset` `@elapsed` timing report to `runtests.jl` (data for M3 group balancing).
- Build the differential equivalence harness (§5.1). Must be green before any migration.
- *Gate:* full suite green · baseline captured · harness green.

**M1 — In-process dispatch (#1, the ~80% win)**
- Add `with_inproc_routes` to `test_http.jl` (router-liveness check + `serialize=true`, per §4).
- Migrate route files one commit each (green-before/green-after, §5.2). Keepers stay on the wire. While migrating `test_idempotency_replay_invariant.jl`, **drop the now-dead `sleep(0.3)` ×4** — under `internalrequest` the post-commit broadcast flush is synchronous in the test task (`events.jl`, `idempotency.jl:189`); replace with a zero-timeout drain.
- *Gate:* each file green in-process, unchanged count · milestone end: full suite green, pass count **==** baseline. **Go/no-go threshold:** M1 must reach **< 10 min** (target ≈5–8 min); if not, revisit the keeper boundary before starting M2 rather than proceeding.

**M2 — Template DB (#3)** — *scoped to the file-backed `open_db` population only.*
- The suite runs **two FK regimes**: ~11 files using in-memory `SQLite.DB()`+`create_schema!` (FK **OFF**) vs ~33 files using `open_db` (FK **ON**). `PRAGMA foreign_keys=ON` is **connection-scoped**, set only inside `open_db` (`db.jl:1907`), and does **not** survive `VACUUM INTO`. So M2 replaces only the **`open_db`-backed** fixtures; leave the `:memory:` tests as-is (a uniform template would silently flip one family's FK enforcement). Enumerate the exact file list in the implementation plan, not "the non-server slice."
- Build the pre-migrated template once at suite start via `VACUUM INTO` (source DB WAL-checkpointed + closed first → a **delete-mode, WAL-free single file**; target path fresh per clone; `cp` copies only the main file, no `-wal`/`-shm` glob; template stays read-only during the run). The template is **schema-only / zero-row** (pure DDL ⇒ empty `sqlite_sequence` ⇒ first-insert id = 1, preserving `test_db.jl:65/72/79`'s `== 1` asserts).
- Define an **"open prepared clone" helper** that opens a clone and re-applies the connection-scoped setup `open_db` does — `PRAGMA foreign_keys=ON`, WAL, `finalize_statements!` — **without** re-running `create_schema!`/`migrate_schema!`. This is a deliberate **fidelity change**, not a transparent 1:1 swap: a migrated+WAL+on-disk template is a *superset* of `create_schema!` output (`migrate_schema!` is additive), so existing column-presence assertions still pass, and the §5.2 green-before/after gate backstops it.
- For the migration-upgrade testsets in `test_db.jl` that build a legacy-shaped DB then assert `migrate_schema!`, prefer a named **`with_legacy_db` opt-out helper** over a skip-list (opt-out, not opt-in); the §5.3 count==baseline gate backstops either way.
- *Gate:* full suite green, pass count **==** baseline, wall time measured; plus a **"FK-violation-throws against a clone"** assertion proving the prepared-clone helper restored `foreign_keys=ON`.

**M3 — Named GROUP buckets + parallel runner (#4)**
- Restructure `runtests.jl` into named `GROUP` buckets (`db`, `routes`, `events`, `wire`, `pipeline` — balanced using M0 timing data). `GROUP=All` (default) reproduces the current serial order exactly = the bisect fallback. The wire tier (§5) is the `wire` group; `test_concurrent_writes.jl` is **not** in it (scrub-only, §5).
- **Each GROUP bucket must include `test_http.jl` first** (it carries `with_test_server`/`with_inproc_routes`); the existing `isdefined(…with_test_server)` self-include guard (e.g. `test_assignment_reattach.jl`) survives M1 unchanged. Verify `GROUP=All` preserves exact serial order.
- Add `make test-parallel`: launches one process per group, aggregates exit codes. **Sysimage caveat:** the app image (`make sysimage` / `precompile_workload.jl`) exercises only `open_db`+CLI — **no `internalrequest`, no `Test`** — so it speeds *package load* but not first-*dispatch*/`Test`-compile TTFX. Optionally add a route-dispatch workload to the precompile file. **Re-measure on M0/M1 timing before enabling parallel-by-default**; local default stays serial `GROUP=All`.
- *Gate:* `GROUP=All` green **==** baseline · each group green standalone · parallel run green · wall time measured.

**Cumulative target:** M1 alone ≈20 min → ≈5–8 min; + M2 trims the migration residual; + M3 over the sysimage → **< 3 min**.

## 7. Isolation hazards (every milestone must handle)

Scrub the shared process-globals between tests exactly as `stop_test_server!` (`server.jl:210-218`) does:

1. **`_DB_REF` (`server.jl:7`)** — sole DB accessor. Re-`bind_db!` per test; null to `nothing` in `finally` so a stray post-test request 500s loudly.
2. **`SSE_SUBSCRIBERS` + `SSE_LOCK` (`server.jl:35-36`)** — reset `SSE_SUBSCRIBERS[] = []` under `SSE_LOCK` in `finally`. SSE-broadcast tests push a fake `(pending=Channel,)` subscriber and drain in-process; verify by draining the channel, not by assuming the wire delivered (queue SSE-wins resolves with a partial).
3. **`OP_LOCKS` (`idempotency.jl:30`)** — `empty!(OP_LOCKS)` in `finally`. Durable idempotency cache lives in each test's own DB. Lock-order invariant (`server.jl:21-28`): `OP_LOCKS_MU` (`idempotency.jl:31`) is **never held across** a `_DB_WRITE_LOCK` acquisition (`with_idempotency` does `OP_LOCKS_MU → release → _DB_WRITE_LOCK`); don't add anything under `OP_LOCKS_MU` that may acquire `_DB_WRITE_LOCK`.
4. **GC timer (`server.jl:144`)** — do not start it in the in-process path.
5. **Router registry (B1)** — register via a **liveness check**, not a sticky boolean: re-`register_routes!()` whenever Oxygen's router is empty. A wire keeper's `resetstate()` wipes the router mid-suite, so under `GROUP=All` an in-process file that runs after a wire keeper must re-register. Never *unconditionally* re-register per test (that re-introduces the cost #1 removes).
6. **`HIMALAYA_FRONTEND_DIST`** — `register_routes!` mounts the `/**` SPA catch-all + `dynamicfiles` **conditionally** on `isdir(dist_dir)` at registration time (`server.jl:53-72`), so registration is environment-dependent. Before the harness registers routes, neutralize `HIMALAYA_FRONTEND_DIST` to a known-absent path so in-process registration is deterministic. (`test_spa_fallback.jl`, which *wants* the SPA mounted, stays a wire keeper. The `/api` 404-masking concern is a non-issue — the explicit `api/` guard at `server.jl:66` handles it; only non-`api` unknown paths flip, owned solely by the keeper.)
7. **Precompile/sysimage** — the M3 parallel runner re-loads `HimalayaUI`+`Oxygen` per process; favor few fat groups + the shared sysimage so the load tax doesn't erase the parallel win (see M3 sysimage caveat).
8. **Process ENV restoration** — tests that mutate `HIMALAYA_*` ENV (`test_config`, `test_routes_image`, `test_routes_experiments`) now share one process with later tests; confirm each restores ENV in a `finally` so a mutation doesn't leak where a fresh server boot previously masked it.

## 8. What NOT to do

- **Don't reuse a DB connection or storage across tests.** Sharing the `_DB_REF` *pointer* (rebound per test) is fine; sharing the underlying SQLite storage is not. The template-clone (#3) must hand each test separate storage over a fresh connection.
- **Don't adopt rollback-isolation** for the DB. `_reingest_inner!` runs its own `BEGIN/COMMIT` (collides with a harness-owned outer transaction), and rollback forces a single shared connection that WAL serializes — contradicting the `parallel=true` write path.
- **Don't rewrite into `@testitem`/ReTestItems or ParallelTestRunner.** ReTestItems' per-worker setup gives the same boot-once win #1 already captures in-process, at the cost of an invasive rewrite of an include-order-dependent suite (its documented worst case). ParallelTestRunner's "remove all includes" file-granularity model breaks our load-bearing include order, fights precompile economics (up to 50 process loads), and its isolation is something `with_inproc_routes`' scrub already provides. Revisit ParallelTestRunner only if M3 timing data shows named groups can't balance the floor.
- **Don't use numeric shard indices.** Named groups are reproducible (`GROUP=events` is always the events files); numeric bins shift when any file is added.
- **Don't break the `runtests.jl` include order.** `test_http.jl` must precede any file using `with_test_server`/`with_inproc_routes`. `GROUP=All` must hold all files in current order as the serial bisect fallback.
- **Don't *unconditionally* re-run `register_routes!()` per test, and don't call `resetstate()` in the in-process `finally`** — re-register only when the router is empty (the B1 liveness check). Unconditional re-registration re-introduces the cost #1 removes; calling `resetstate()` in the in-process path wipes the router the next file needs.

## 9. Open verification gates (resolve during implementation)

- **VACUUM INTO WAL safety:** the migrated source DB is checkpointed + closed before `VACUUM INTO` → a delete-mode WAL-free single file; `cp` copies only the main file (no `-wal`/`-shm` glob); the template stays read-only during the run; each clone uses a fresh path.
- **FK enforcement across the clone:** confirm the "open prepared clone" helper re-applies `PRAGMA foreign_keys=ON` (connection-scoped, lost by `VACUUM INTO`); gate with a FK-violation-throws assertion. M2 touches only the `open_db`-backed files, not the `:memory:` family (§6 M2).
- **Per-file keeper classification:** confirm each SSE file against source before deciding wire-vs-in-process; `test_concurrent_writes.jl` is scrub-only, not wire (§5).
- **Equivalence-harness header comparison (allowlist):** drop a fixed transport-only set (`Date`, `Server`, `Transfer-Encoding`, `Connection`, `Content-Length`) and assert **all other** headers match — an allowlist of drops, not a denylist of checks — so contract headers can't slip. Idempotency-replay row uses a fresh clone DB + distinct `client_op_id` per transport (§5.1).
