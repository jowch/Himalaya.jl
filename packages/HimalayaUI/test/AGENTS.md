# packages/HimalayaUI/test — Backend tests

Julia test suite (stdlib `Test`: `@testset`, `@test`, `@test_throws`). Internal helpers are accessed via `HimalayaUI.<name>`.

## Running the suite (capture once)

The suite is **~3 min serial / ~2 min parallel** (2026-06 perf refactor: in-process route dispatch + template-DB clone + GROUP buckets). Still **don't re-run with different `| grep` filters** — capture once, grep the file.

```bash
# Serial — exact historical order; the bisect fallback. ~3 min.
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-test.out; tail -50 /tmp/jl-test.out

# Parallel — 5 GROUP buckets as separate processes, wall-bounded by the slowest (routes). ~2 min.
make test-parallel                  # logs per bucket → build/test-<g>.log

# One bucket only (db | pipeline | routes | events | wire):
GROUP=routes julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```

`GROUP=All` (default) runs every file in the historical order. The buckets partition that same set; a load-time **drift guard** in `runtests.jl` errors if the bucket lists, `ALL_ORDER`, and the `test_*.jl` files on disk ever diverge — so a new test file can't silently run in no configuration. When you add a `test_*.jl`, add it to **both** `ALL_ORDER` and a bucket in `runtests.jl` (and, for the parallel runner, the bucket already exists in the `Makefile` `GROUPS` line).

## Route tests dispatch in-process — `with_inproc_routes`, not a socket

The default fixture for route tests is `with_inproc_routes` (`test_inproc.jl`), which dispatches `HTTP.Request`s straight through `Oxygen.internalrequest` against the same global Oxygen `CONTEXT[]` — no server boot, no socket, no port. It's the fast analogue of the old `with_test_server`.

```julia
with_inproc_routes(db) do call
    r = call("GET", "/api/experiments/$exp_id/samples")
    @test r.status == 200
    list = JSON3.read(String(r.body))
    # POST/PATCH/DELETE: headers + a Vector{UInt8} body
    r = call("PATCH", "/api/samples/$sid";
             headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
             body = Vector{UInt8}(JSON3.write(Dict(:notes => "hi"))))
end
```

- `call(method, target; headers=Pair[], body=UInt8[])` returns the same `HTTP.Response` the wire path would. `serialize=true, catch_errors=true` is load-bearing: route throws (`error()`, FK `SQLiteException`) become 4xx/5xx exactly as production `serve()` does — `call` never throws on a 4xx/5xx, so drop any `status_exception=false`.
- The `finally` scrub mirrors `stop_test_server!` **minus `resetstate()`**: nulls `_DB_REF`, clears `SSE_SUBSCRIBERS` under `SSE_LOCK`, empties `OP_LOCKS`. Never call `Oxygen.resetstate()` in the in-process path — it wipes the router the next file needs (a liveness probe of `/api/health` re-registers if a wire keeper cleared it).
- **Equivalence is proven, not assumed.** `test_inproc_equivalence.jl` is a differential harness asserting wire `with_test_server` ≡ in-process `with_inproc_routes` (status + body + contract headers) across JSON / binary-PNG / throw-500 / DELETE-404 / idempotency-replay / numeric-array / SSE. Keep it green; if it ever fails, that's a real dispatch divergence — investigate before trusting any migrated test.

### Wire keepers (still use `with_test_server`)

A small tier stays on a real socket because the streaming/SPA/socket behavior **is** the unit under test: `test_routes_health.jl`, `test_routes_status.jl`, `test_sse.jl`, `test_routes_sse_broadcast.jl`, `test_spa_fallback.jl`. New tests that assert on the SPA catch-all, the SSE wire stream, or actual `parallel=true` socket handling belong here; everything else uses `with_inproc_routes`.

## Test-isolation pattern + template-DB clone

Each testset gets its own DB. For a **fresh, empty** fixture, use `open_prepared_clone` (`test_template_db.jl`) — a drop-in for `open_db(joinpath(tmp, "h.db"))` that clones a pre-migrated, zero-row template (built once per process) instead of re-running `create_schema!`+`migrate_schema!` (~37 ms → ~5 ms). It re-applies the connection-scoped `foreign_keys=ON` + WAL the file copy drops.

```julia
mktempdir() do tmp
    db = open_prepared_clone(tmp)   # fresh empty DB, FK ON, WAL; first insert id == 1
    # ...
end
```

**Keep `open_db` (not the clone) where `open_db`/`migrate_schema!` IS the unit under test** — `test_db.jl` and `test_pipeline.jl` build legacy-shaped DBs and assert the migration runs; cloning a pre-migrated template defeats that. `:memory:` `SQLite.DB()` fixtures (FK OFF) also stay as-is.

## Shared cross-file helpers live in `test_fixtures.jl`

Helpers used by more than one file (`_setup_analyzed_exposure`, `_setup_for_resolve`, `_premint_cmp!`, `_member_payload`) have a single definition in `test_fixtures.jl`, included via `test_http.jl`. Don't redefine them per file — under GROUP bucketing a helper defined in file A is not in scope for file B unless it's in `test_fixtures.jl`. **No `@testset` in the include-only files** (`test_http.jl`, `test_fixtures.jl`, `test_inproc.jl`, `test_template_db.jl`) — they're included once per bucket, so a testset there would multiply Pass counts across buckets.

## Oxygen is a process-global singleton

`_DB_REF`, the Oxygen `CONTEXT[]`, `SSE_SUBSCRIBERS`, and `OP_LOCKS` are module-level globals — one per process. So **two route fixtures can't be live at once in a process**, and parallelism is achieved by **sharding across processes** (the GROUP buckets in `make test-parallel`), not by threads within one process. A route test file still can't run standalone (`julia test/test_routes_*.jl`) without first including `test_http.jl` (which pulls in `test_fixtures.jl` + `test_inproc.jl` + `test_template_db.jl`).

**In-memory fixtures need the migration chain too.** A bare `SQLite.DB(); create_schema!(db)` fixture under-provisions the schema: migration-created tables (`series*`, `comparison*`, `speculative_peak_intents`) won't exist, and `_persist_analysis_inner!` reads `speculative_peak_intents` unconditionally. Any fixture that calls `persist_analysis!`/`analyze_exposure!` must run `migrate_schema!(db)` after `create_schema!(db)` — or just use `open_db` on a tmp path.

## In-process SSE subscriber

To assert SSE fanout, register a `(pending = Channel{String}(64),)` directly on `HimalayaUI.SSE_SUBSCRIBERS[]` under `HimalayaUI.SSE_LOCK` instead of opening an HTTP stream. `test_idempotency_replay_invariant.jl::_capture_sse_during` is the canonical pattern.

**In-process dispatch flushes the broadcast synchronously — no sleep needed.** `broadcast_event!` enqueues frames into a task-local post-commit queue (`task_local_storage(POST_COMMIT_BROADCAST_KEY)`), and `with_idempotency` fires `_flush_post_commit_broadcasts!` after its `SQLite.transaction` commits — all in-task. With `with_inproc_routes`, by the time `call(...)` returns the frame is already on the channel, so the migrated invariant test drops the old `sleep(0.3)`. The **wire path** (a `with_test_server` keeper) still needs the short drain, because the socket round-trip lets the off-task flush lag. (`src/events.jl` enqueue/flush; `test_idempotency_replay_invariant.jl` for the in-process pattern, `test_routes_sse_broadcast.jl` for the wire pattern.)

**Julia `do`-block gotcha**: `do` passes the function as the FIRST arg. `_capture_sse_during("kind") do … end` ⇒ `_capture_sse_during(f, "kind")`. Order your signature accordingly.

## SQLite version is pinned — fresh worktrees hang without it

`Manifest.toml` is gitignored, so a fresh `git worktree add` + instantiate re-resolves SQLite. SQLite_jll 3.53.2 and the SQLite.jl wrapper 1.8.1 both hang `migrate_schema!`'s `ALTER TABLE ADD COLUMN … CHECK` forever (a `sqlite3AlterBeginAddColumn` spin), wedging the suite at `test_pipeline.jl` "rowid PKs to AUTOINCREMENT". `Project.toml` caps both (`SQLite = "=1.8.0"`, `SQLite_jll = "~3.51"`); if a fresh checkout still hangs (one julia proc at 100% CPU, no output), check `Pkg.dependencies()` for SQLite drift and `sample <pid>` to confirm the native spin.

## FK-heal regression tests

`_fix_fk_references_after_autoincrement_migration!(db)` should be invoked **directly** from tests rather than through `open_db`, because `open_db` runs `create_schema!` migrations that expect full production schemas — synthetic FK fixtures (`refs_to_samples`, `_migrate_old_*`) break the migration chain. See `test_db.jl` FK-heal regression tests for the pattern.

## Read-only experiment directories

Himalaya must never create, modify, or delete files inside an experiment dir at runtime (except `himalaya config new`). Ingestion is read-only with respect to the experiment dir and now lives in the HTTP scan path (`scan_and_group!` in `ingest.jl`), not the CLI — the manifest-driven `init`/`reingest` CLI ingest was deleted by the ingestion redesign. Tests that need a populated DB without a real on-disk scan use `seed_experiment!` (`test/seed.jl`), which writes experiment/sample/exposure rows directly via the Phase-A writers.

## Curation lifecycle

`analyze_exposure!` synthesises the effective peak set as `auto_peaks − excludes ∪ adds` (see `pipeline.jl::effective_peaks`). Curation-lifecycle tests in `test_pipeline.jl` pin this contract — keep them green when touching the pipeline.
