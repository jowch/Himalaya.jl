# packages/HimalayaUI/test — Backend tests

Julia test suite (stdlib `Test`: `@testset`, `@test`, `@test_throws`). Internal helpers are accessed via `HimalayaUI.<name>`.

## Slow suite, capture once

The HimalayaUI Julia suite takes 5–10 min per invocation (per-test fixture DBs + Oxygen lifecycle dominate). **Do NOT re-run the same suite with different `| grep` filters** — every invocation rebuilds fixtures from scratch.

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-test.out
tail -50 /tmp/jl-test.out
```

## Test-isolation pattern

Each testset opens its own SQLite DB at an explicit path:

```julia
mktempdir() do tmp
    db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    # ...
end
```

This bypasses `default_db_path()` (`HIMALAYA_DB_PATH` or `~/.himalaya/himalaya.db`) and gives every test its own clean DB. The central DB is for production CLI use; tests must not depend on its state.

## Oxygen is a process-global singleton — one test server at a time

`start_test_server!` (`src/server.jl:197-208`) calls `Oxygen.resetstate()`, `bind_db!(db)` (binds the module-level `_DB_REF` singleton), and `register_routes!`. `stop_test_server!` (`:210-218`) must run between servers or state leaks: it terminates Oxygen, nulls `_DB_REF[]`, and clears `SSE_SUBSCRIBERS[]` + `OP_LOCKS`. `with_test_server` (`test_http.jl:4-12`) pairs them — `start_test_server!` then `stop_test_server!` in `finally`.

Two consequences:

- **`runtests.jl` include order is load-bearing.** `with_test_server` lives in `test_http.jl`, which `runtests.jl` includes before any route test. A route test file cannot be run standalone (`julia test/test_routes_*.jl`) without first including `test_http.jl`.
- **The suite cannot be naively parallelized.** One Oxygen state + one `_DB_REF` for the whole process means two servers up at once corrupt each other.

## In-process SSE subscriber

To assert SSE fanout in Julia tests, register a `(pending = Channel{String}(64),)` directly on `HimalayaUI.SSE_SUBSCRIBERS[]` under `HimalayaUI.SSE_LOCK` instead of opening an HTTP streaming connection. Faster, deterministic, no port management. `test_idempotency_replay_invariant.jl::_capture_sse_during` is the canonical pattern.

**Julia `do`-block gotcha**: `do` syntax passes the function as the FIRST arg. `_capture_sse_during("kind") do ... end` ⇒ `_capture_sse_during(f, "kind")`. Order your helper signature accordingly.

**`_capture_sse_during` must sleep before draining frames.** It sleeps ~0.3s after the HTTP call returns, because `broadcast_event!` is NOT synchronous inside the route handler. Peak/curation routes enqueue the broadcast into a task-local post-commit queue (`task_local_storage(POST_COMMIT_BROADCAST_KEY)`); `with_idempotency` fires it via `_flush_post_commit_broadcasts!` only AFTER its `SQLite.transaction` commits. A new SSE-invariant test that omits or shortens the sleep races the flush and sees zero frames. See `test_idempotency_replay_invariant.jl:49-51`; `src/events.jl:206-225` (enqueue), `:258-270` (flush after commit).

## FK-heal regression tests

`_fix_fk_references_after_autoincrement_migration!(db)` should be invoked **directly** from tests rather than through `open_db`, because `open_db` runs `create_schema!` migrations that expect full production schemas — synthetic FK fixtures (`refs_to_samples`, `_migrate_old_*`) break the migration chain. See `test_db.jl` FK-heal regression tests for the pattern.

## Read-only-dir regression

`test_pipeline.jl` snapshots an experiment-directory's contents before/after `cli_init_with_db!` and asserts equality. Keep this test green — Himalaya must never create, modify, or delete files inside an experiment dir at runtime (except `himalaya config new`).

## Curation lifecycle

`analyze_exposure!` synthesises the effective peak set as `auto_peaks − excludes ∪ adds` (see `pipeline.jl::effective_peaks`). Curation-lifecycle tests in `test_pipeline.jl` pin this contract — keep them green when touching the pipeline.
