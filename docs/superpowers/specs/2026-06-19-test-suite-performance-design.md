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
| **#1 In-process dispatch** | ✅ Confirmed — Oxygen's own suite uses `internalrequest` as its default route-test path; real `serve()` is reserved for transport tests (`ssetests.jl`, `streamingtests.jl`). | Middleware-parity gate **closed**: HimalayaUI passes **no** global middleware to `Oxygen.serve` (`server.jl:193,204`), and all cross-cutting concerns (`with_idempotency`, event logging) are **handler-body wrappers** (e.g. `routes_exposures.jl:123`), so `internalrequest`'s default `middleware=[]` drops nothing. |
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

- **Once per process** (guarded by a module-level `Ref{Bool}`): `Oxygen.resetstate()` + `HimalayaUI.register_routes!()`. Never per test.
- **Per call:** `HimalayaUI.bind_db!(db)`, then `call(method, target; headers, body)` constructs an `HTTP.Request` and returns `Oxygen.internalrequest(req; serialize=false)`. `serialize=false` is **load-bearing** — handlers already return `HTTP.Response`.
- **`finally` scrub** (mirrors `stop_test_server!` exactly): `_DB_REF[] = nothing`; `SSE_SUBSCRIBERS[] = []` under `SSE_LOCK`; `empty!(OP_LOCKS)`; **GC timer not started**.

Migration is a mechanical per-file swap: `with_test_server(db) do port, base … HTTP.get("$base/api/…")` → `with_inproc_routes(db) do call … call("GET", "/api/…")`. The closure returns the identical `HTTP.Response`, so every existing `r.status` / `JSON3.read(String(r.body))` assertion is **unchanged** — keeping the diff reviewable line-by-line.

**What in-process dispatch deliberately does NOT exercise:** the TCP socket, byte-level HTTP parse, and Oxygen's `async`/`parallel=true` request scheduling. Those belong to the wire tier (§5).

## 5. Section 2 — equivalence guarantee + keeper boundary

Three decoupled safety layers, so no single mechanism carries the whole burden:

1. **Differential equivalence harness** (`test_inproc_equivalence.jl`, built first, kept permanently). Fires a representative matrix — `GET/POST/PATCH/DELETE`, `200`/`4xx`/`5xx`, an idempotency-replay pair, an SSE-broadcast-then-drain — through **both** a real `with_test_server` socket **and** `with_inproc_routes`, asserting equivalent `HTTP.Response` — identical status and body, and headers identical after normalizing away transport-only headers (`Date`, `Server`, etc.; see §9). Proves the primitive once; guards against future Oxygen drift. *Answers "is dispatch faithful."*
2. **Per-file green-before/green-after gate.** Each file migrated as its own commit: green-on-wire → mechanical swap → green-in-process, same test count. *Answers "did I migrate this file correctly."* Decoupled from layer 1: a red file with a green harness ⇒ a migration slip, not a dispatch problem.
3. **Full-suite milestone gate.** The slow suite (`> /tmp/jl-test.out`, capture-once per `test/AGENTS.md`) runs at each milestone boundary, asserting pass count **==** the M0 baseline.

**Keeper boundary (conservative, verified per-file against source).** A test stays on `with_test_server` only if it depends on real transport or real concurrent scheduling. Current best estimate, to confirm per-file:
- SSE over the wire / streaming: `test_sse.jl`, `test_routes_sse_broadcast.jl`, `test_idempotency_sse_suppression.jl` (confirm each — broadcast tests that only push a fake subscriber and drain a `Channel` can move in-process).
- Concurrency / `_DB_WRITE_LOCK` races: `test_concurrent_writes.jl` and any idempotency race test that asserts real parallel dispatch.
- Static-file / env: `test_spa_fallback.jl` (mutates `HIMALAYA_FRONTEND_DIST`, serves `index.html` over the wire), `test_routes_health.jl`, `test_routes_status.jl`.

Wire-tier files boot the server **once per file** (Oxygen `ssetests.jl` pattern), keep the ephemeral `find_free_port()`, and wrap teardown under a timeout so a hung server can't wedge the group.

## 6. Section 3 — milestone sequencing

Four independently landable, independently revertable milestones, with a green gate between each.

**M0 — Baseline + harness** *(no behavior change)*
- Capture the full-suite green baseline (exact pass count + wall time) = the regression floor.
- Add per-`@testset` `@elapsed` timing report to `runtests.jl` (data for M3 group balancing).
- Build the differential equivalence harness (§5.1). Must be green before any migration.
- *Gate:* full suite green · baseline captured · harness green.

**M1 — In-process dispatch (#1, the ~80% win)**
- Add `with_inproc_routes` to `test_http.jl`.
- Migrate route files one commit each (green-before/green-after, §5.2). Keepers stay on the wire.
- *Gate:* each file green in-process, unchanged count · milestone end: full suite green, pass count **==** baseline, wall time measured (expect ≈20 min → ≈5–8 min).

**M2 — Template DB (#3)**
- Build a pre-migrated template once at suite start via `VACUUM INTO` (source DB WAL-checkpointed + closed first; target path fresh per clone). Each test gets a fresh `cp` + connection, `close+rm` in `finally`. Replace per-test `open_db`+`create_schema!` in the non-server slice.
- Do **not** convert the migration-upgrade testsets in `test_db.jl` that build a legacy-shaped DB then assert `migrate_schema!` — the template factory is opt-in and must skip them.
- *Gate:* full suite green, pass count **==** baseline, wall time measured.

**M3 — Named GROUP buckets + parallel runner (#4)**
- Restructure `runtests.jl` into named `GROUP` buckets (`db`, `routes`, `events`, `wire`, `pipeline` — balanced using M0 timing data). `GROUP=All` (default) reproduces the current serial order exactly = the bisect fallback. The wire tier (§5) is the `wire` group.
- Add `make test-parallel`: launches one process per group over the sysimage (`julia --sysimage build/himalaya.so`), aggregates exit codes. Local default stays serial `GROUP=All`; parallel-on-by-default gated on post-M1 timing + sysimage making per-process load negligible.
- *Gate:* `GROUP=All` green **==** baseline · each group green standalone · parallel run green · wall time measured.

**Cumulative target:** M1 alone ≈20 min → ≈5–8 min; + M2 trims the migration residual; + M3 over the sysimage → **< 3 min**.

## 7. Isolation hazards (every milestone must handle)

Scrub the shared process-globals between tests exactly as `stop_test_server!` (`server.jl:210-218`) does:

1. **`_DB_REF` (`server.jl:7`)** — sole DB accessor. Re-`bind_db!` per test; null to `nothing` in `finally` so a stray post-test request 500s loudly.
2. **`SSE_SUBSCRIBERS` + `SSE_LOCK` (`server.jl:35-36`)** — reset `SSE_SUBSCRIBERS[] = []` under `SSE_LOCK` in `finally`. SSE-broadcast tests push a fake `(pending=Channel,)` subscriber and drain in-process; verify by draining the channel, not by assuming the wire delivered (queue SSE-wins resolves with a partial).
3. **`OP_LOCKS` (`idempotency.jl:30`)** — `empty!(OP_LOCKS)` in `finally`. Durable idempotency cache lives in each test's own DB. Lock-order invariant (`server.jl:21-28`): `OP_LOCKS_MU` (`idempotency.jl:31`) is **never held across** a `_DB_WRITE_LOCK` acquisition (`with_idempotency` does `OP_LOCKS_MU → release → _DB_WRITE_LOCK`); don't add anything under `OP_LOCKS_MU` that may acquire `_DB_WRITE_LOCK`.
4. **GC timer (`server.jl:144`)** — do not start it in the in-process path.
5. **Router registry** — `register_routes!()` once per process, never per test.
6. **Precompile/sysimage** — the M3 parallel runner re-loads `HimalayaUI`+`Oxygen` per process; favor few fat groups + the shared sysimage so the load tax doesn't erase the parallel win.

## 8. What NOT to do

- **Don't reuse a DB connection or storage across tests.** Sharing the `_DB_REF` *pointer* (rebound per test) is fine; sharing the underlying SQLite storage is not. The template-clone (#3) must hand each test separate storage over a fresh connection.
- **Don't adopt rollback-isolation** for the DB. `_reingest_inner!` runs its own `BEGIN/COMMIT` (collides with a harness-owned outer transaction), and rollback forces a single shared connection that WAL serializes — contradicting the `parallel=true` write path.
- **Don't rewrite into `@testitem`/ReTestItems or ParallelTestRunner.** ReTestItems' per-worker setup gives the same boot-once win #1 already captures in-process, at the cost of an invasive rewrite of an include-order-dependent suite (its documented worst case). ParallelTestRunner's "remove all includes" file-granularity model breaks our load-bearing include order, fights precompile economics (up to 50 process loads), and its isolation is something `with_inproc_routes`' scrub already provides. Revisit ParallelTestRunner only if M3 timing data shows named groups can't balance the floor.
- **Don't use numeric shard indices.** Named groups are reproducible (`GROUP=events` is always the events files); numeric bins shift when any file is added.
- **Don't break the `runtests.jl` include order.** `test_http.jl` must precede any file using `with_test_server`/`with_inproc_routes`. `GROUP=All` must hold all files in current order as the serial bisect fallback.
- **Don't re-run `register_routes!()`/`resetstate()` per test** in the in-process path — that re-introduces the exact cost #1 removes.

## 9. Open verification gates (resolve during implementation)

- **VACUUM INTO WAL safety:** confirm the migrated source DB is checkpointed + closed before the template is built, and each per-test scratch DB starts with no live `-wal` sidecar.
- **Per-file keeper classification:** confirm each SSE/concurrency file against source before deciding wire-vs-in-process.
- **Equivalence-harness header normalization:** decide how to compare headers (sort + drop transport-only headers like `Server`/`Date`) so the byte-identical assertion is meaningful, not flaky.
