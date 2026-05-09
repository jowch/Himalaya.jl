# Permalink URLs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Slug-based permalink URLs for every page in HimalayaUI, plus a `/api/resolve` endpoint and a 404 page for stale URLs. Per spec [`docs/superpowers/specs/2026-05-09-permalink-urls-design.md`](../specs/2026-05-09-permalink-urls-design.md).

**Architecture:** A new read-only Julia route `/api/resolve` translates `(experiment, sample, exposure)` slugs → IDs in one round trip; an SPA catch-all (`/**`) serves `index.html` for unknown paths. The frontend integrates with the existing `react-router-dom` setup: `useStateFromUrl` reads `useLocation()` and dispatches Zustand setters; `useUrlFromState` writes via `useNavigate(target, { replace })`. New `<Route>` declarations cover index/inspect URLs; Compare routing is unchanged. A `<StaleUrlPage>` component renders 404s in-place while `recoverFromStaleUrl(opts)` provides one-click recovery.

**Tech Stack:** Julia 1.11 + Oxygen.jl 1.10 + SQLite.jl (backend), React 18 + TypeScript strict + Zustand + TanStack Query 5 + react-router-dom v6 (already in use) + Vitest + Playwright (frontend).

**Architecture amendment** vs. the original spec: the original said "no router library." In fact `BrowserRouter` is already mounted (`main.tsx:8,36`); `<Routes>` is in `AppShell.tsx`; `TabRocker`, `ComparePage`, `ConflictModal`, etc. all use `useNavigate`/`useLocation`/`useParams`. The plan integrates with this rather than rolling our own — `useStateFromUrl` reads `useLocation()`, `useUrlFromState` calls `navigate(target, { replace })`, and AppShell adds new `<Route>` declarations for index/inspect alongside the existing Compare ones. The spec was amended in-place to reflect this.

**Pre-existing dependency:** This branch is built on top of #88 which has already landed at `8ac2bf6`. `samples.name` is stable, convention-enforced, and unique within experiment.

**Key spec sections:**
- §3.1 — `/api/resolve` route + `migrate_exposures_unique_filename!` migration
- §3.2 — SPA catch-all in `server.jl`
- §4.1 — `parseLocation` parser + `ParsedUrl` discriminated union
- §4.2 — `useStateFromUrl` (URL → Zustand) + `<ResolvingFallback>`
- §4.3 — `useUrlFromState` (Zustand → URL) + push/replace policy
- §4.4 — Mounting order + named-action `staleUrlContext` clearing
- §5 — `/` cold-mount redirect via `qc.getQueryData` + resolve-by-id fallback
- §6 — `<StaleUrlPage>` + `recoverFromStaleUrl` + `StaleUrlContext` type + AppShell ladder
- §7 — SSE-driven URL invalidation (no `client_id` filter)
- §8 — Test plan
- §9 — File-level changes

---

## File Structure

**Backend (new):**
- `packages/HimalayaUI/src/routes_resolve.jl` — `register_resolve_routes!` + 200/400/404 handler.
- `packages/HimalayaUI/test/test_routes_resolve.jl` — happy path, 404 variants, 400, tiebreaker.
- `packages/HimalayaUI/test/test_spa_fallback.jl` — `/foo`, `/inspect/exp/sample`, `/api/foo`, asset paths.

**Backend (modified):**
- `packages/HimalayaUI/src/HimalayaUI.jl` — `include("routes_resolve.jl")`.
- `packages/HimalayaUI/src/server.jl` — `register_resolve_routes!()` call + `/**` catch-all.
- `packages/HimalayaUI/src/db.jl` — new `migrate_exposures_unique_filename!(db)` helper called from `migrate_schema!` after `_fix_fk_references_after_autoincrement_migration!`.
- `packages/HimalayaUI/test/test_db.jl` — direct-invocation tests for the migration helper.
- `packages/HimalayaUI/test/test_route_response_shapes.jl` — resolve 200/400/404 shape rows.
- `packages/HimalayaUI/test/runtests.jl` — register the two new test files.

**Frontend (new):**
- `frontend/src/lib/url/parseLocation.ts` — pure parser.
- `frontend/test/parseLocation.test.ts` — round-trip every URL shape.
- `frontend/src/lib/url/emitMode.ts` — push/replace mode flag for the next URL emit (hoisted to break a state.ts ↔ useUrlFromState cycle).
- `frontend/src/hooks/useStateFromUrl.ts` — reads `useLocation()`; dispatches Zustand setters; origin-tagged fetches.
- `frontend/test/useStateFromUrl.test.tsx` — origin-tag race; pre-fetch staleUrlContext clear; / redirect.
- `frontend/src/hooks/useUrlFromState.ts` — calls `useNavigate()`; subscribes to experiments + samples queries; equality guard.
- `frontend/test/useUrlFromState.test.tsx` — push/replace policy + replay-as-rerun no-spurious-emit.
- `frontend/src/components/StaleUrlPage.tsx` — 404 page.
- `frontend/test/StaleUrlPage.test.tsx` — per-variant render + CTA dispatch.
- `frontend/src/components/ResolvingFallback.tsx` — empty layout-bones placeholder.
- `frontend/e2e/permalinks.spec.ts` — Playwright mocked.
- `frontend/e2e/live/permalinks.spec.ts` — Playwright live integration.

**Frontend (modified):**
- `frontend/src/state.ts` — add `staleUrlContext`, `resolving`, `setStaleUrlContext`, `setResolving`, `recoverFromStaleUrl`; `setActive*` clear `staleUrlContext`. Export `StaleUrlContext` and `RecoverOpts` types.
- `frontend/src/api.ts` — add `ResolveSuccess`, `ResolveError404`, `ResolveError400` types and `resolve()` fetcher.
- `frontend/src/components/AppShell.tsx` — render `<ResolvingFallback>` / `<StaleUrlPage>` / `<PageRouter>` ladder under chrome.
- `frontend/src/App.tsx` — mount `useStateFromUrl()` and `useUrlFromState()`.

---

## Task 1: Migration helper `migrate_exposures_unique_filename!`

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (add helper + call from `migrate_schema!`)
- Modify: `packages/HimalayaUI/test/test_db.jl` (add tests using the FK-heal direct-invocation pattern)

- [ ] **Step 1: Write the failing test**

In `packages/HimalayaUI/test/test_db.jl` add a new `@testset` after the existing FK-heal tests:

```julia
@testset "migrate_exposures_unique_filename!" begin
    @testset "clean DB: index exists, no warnings" begin
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
            try
                # open_db already runs the migration; assert the index is in place.
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT name FROM sqlite_master WHERE type='index' AND name='exposures_unique_filename'"))
                @test length(rows) == 1
            finally
                close(db)
            end
        end
    end

    @testset "idempotent re-run: no-op" begin
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
            try
                # Calling the helper directly a second time should be a no-op (no error,
                # no warnings, index still present).
                HimalayaUI.migrate_exposures_unique_filename!(db)
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT name FROM sqlite_master WHERE type='index' AND name='exposures_unique_filename'"))
                @test length(rows) == 1
            finally
                close(db)
            end
        end
    end

    @testset "synthetic duplicates: rename + warn + index created" begin
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
            try
                # Drop the index so we can synthesize duplicates and re-run the helper.
                DBInterface.execute(db, "DROP INDEX IF EXISTS exposures_unique_filename")

                # Seed a sample + two duplicate exposures via raw SQL (bypassing the upsert).
                eid = let res = DBInterface.execute(db,
                          "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e','/p','/d','/a')")
                    Int(DBInterface.lastrowid(res))
                end
                sid = let res = DBInterface.execute(db,
                          "INSERT INTO samples (experiment_id, name) VALUES (?, 'S1')", [eid])
                    Int(DBInterface.lastrowid(res))
                end
                # Two rows with the same (sample_id, filename).
                DBInterface.execute(db,
                    "INSERT INTO exposures (sample_id, filename, kind) VALUES (?, 'JC001-007', 'simple')", [sid])
                DBInterface.execute(db,
                    "INSERT INTO exposures (sample_id, filename, kind) VALUES (?, 'JC001-007', 'simple')", [sid])

                # Run the helper directly (FK-heal pattern).
                @test_logs (:warn, r"Renamed duplicate exposure"i) min_level=Logging.Warn match_mode=:any begin
                    HimalayaUI.migrate_exposures_unique_filename!(db)
                end

                # Oldest id keeps the bare name; second is suffixed.
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, filename FROM exposures WHERE sample_id = ? ORDER BY id ASC", [sid]))
                @test length(rows) == 2
                @test rows[1].filename == "JC001-007"
                @test rows[2].filename == "JC001-007-2"

                # Index now present.
                idx = Tables.rowtable(DBInterface.execute(db,
                    "SELECT name FROM sqlite_master WHERE type='index' AND name='exposures_unique_filename'"))
                @test length(idx) == 1
            finally
                close(db)
            end
        end
    end
end
```

Add `using Logging` to the top of the test file if not already present.

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /opt/Himalaya.jl/.claude/worktrees/permalinks
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|fail|did not pass" /tmp/jl-test.out
```

Expected: failure ("UndefVarError: migrate_exposures_unique_filename! not defined" or "no such index") — the helper doesn't exist yet.

- [ ] **Step 3: Add the helper to `db.jl`**

In `packages/HimalayaUI/src/db.jl`, add a new helper after `migrate_samples_naming!`:

```julia
"""
    migrate_exposures_unique_filename!(db)

Add `UNIQUE INDEX exposures_unique_filename ON exposures(sample_id, filename)`
following the dedupe-then-enforce pattern from `migrate_samples_naming!`.
Renames pre-existing duplicates (oldest id keeps the bare filename;
second-and-later get `<filename>-2`, …) before creating the index, so the
`CREATE UNIQUE INDEX` always succeeds against clean data.

Idempotent on re-run. Wrapped in `SQLite.transaction` so a partial run
cannot leave duplicates renamed without uniqueness enforcement.

Direct-invocation pattern from CLAUDE.md (FK-heal regression tests) —
tests can call this without going through `open_db`.
"""
function migrate_exposures_unique_filename!(db::SQLite.DB)
    SQLite.transaction(db) do
        # Track existing (sample_id, filename) pairs so a row literally named
        # "<x>-2" doesn't collide with our rename target.
        existing = Set{Tuple{Int64,String}}(
            (Int64(r.sample_id), String(r.filename))
            for r in Tables.rowtable(DBInterface.execute(db,
                "SELECT sample_id, filename FROM exposures WHERE sample_id IS NOT NULL AND filename IS NOT NULL")))

        dups = Tables.rowtable(DBInterface.execute(db, """
            SELECT sample_id, filename FROM exposures
            WHERE sample_id IS NOT NULL AND filename IS NOT NULL
            GROUP BY sample_id, filename HAVING COUNT(*) > 1"""))

        for d in dups
            ids = Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM exposures WHERE sample_id = ? AND filename = ? ORDER BY id ASC",
                [d.sample_id, d.filename]))
            for (i, row) in enumerate(ids)
                i == 1 && continue  # oldest keeps the bare name
                suffix_n = i
                new_name = "$(d.filename)-$(suffix_n)"
                while (Int64(d.sample_id), new_name) in existing
                    suffix_n += 1
                    new_name = "$(d.filename)-$(suffix_n)"
                end
                push!(existing, (Int64(d.sample_id), new_name))
                @warn "Renamed duplicate exposure" sample_id=d.sample_id old=d.filename new=new_name id=row.id
                DBInterface.execute(db,
                    "UPDATE exposures SET filename = ? WHERE id = ?",
                    [new_name, row.id])
            end
        end

        DBInterface.execute(db,
            "CREATE UNIQUE INDEX IF NOT EXISTS exposures_unique_filename ON exposures(sample_id, filename)")
    end
    nothing
end
```

Wire the call into `migrate_schema!`. Find the existing block (around `db.jl:312–319`):

```julia
    migrate_samples_naming!(db)
    migrate_pk_to_autoincrement!(db)
    # ...
    _fix_fk_references_after_autoincrement_migration!(db)
```

Add the new call **after** `_fix_fk_references_after_autoincrement_migration!` (so the AUTOINCREMENT rebuild has settled and the index attaches to the rebuilt `exposures` table — placing it earlier would have it dropped along with `_migrate_old_exposures`):

```julia
    _fix_fk_references_after_autoincrement_migration!(db)
    migrate_exposures_unique_filename!(db)
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|fail|did not pass" /tmp/jl-test.out
```

Expected: all assertions pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_db.jl
git commit -m "feat(db): add migrate_exposures_unique_filename! helper

Dedupe-then-enforce migration mirroring migrate_samples_naming! from #88.
Called from migrate_schema! after _fix_fk_references_after_autoincrement_migration!
so the AUTOINCREMENT rebuild has settled before the index attaches."
```

---

## Task 2: Resolve route — name form (200 + 404 happy paths)

**Files:**
- Create: `packages/HimalayaUI/src/routes_resolve.jl`
- Create: `packages/HimalayaUI/test/test_routes_resolve.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl` (add `include`)
- Modify: `packages/HimalayaUI/src/server.jl` (call `register_resolve_routes!()`)
- Modify: `packages/HimalayaUI/test/runtests.jl` (register new test file)

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/test_routes_resolve.jl`:

```julia
using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# Spec §3.1 — read-only `/api/resolve` route. Translates slugs → IDs in
# one round trip. Read-only — no with_idempotency, no SSE, no events.
#
# Helpers `_setup_analyzed_exposure` (defined in test_route_response_shapes.jl)
# and `with_test_server` (defined in test_http.jl) are available because
# runtests.jl includes those files before this one. See Step 3 for the
# correct ordering when adding the new include lines.

# `_setup_for_resolve` is defined in test_route_response_shapes.jl
# (alongside `_setup_analyzed_exposure`) so both Task 2 and Task 4 can
# call it. See Task 2 Step 3 for the helper body.

@testset "GET /api/resolve" begin
    @testset "200: experiment + sample + exposure happy path" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp&sample=S1&exposure=JC001-007")
                @test r.status == 200
                body = JSON3.read(String(r.body))
                @test body.experiment_id == ctx.experiment_id
                @test body.experiment_name == "test-exp"
                @test body.sample_id == ctx.sample_id
                @test body.sample_name == "S1"
                @test body.exposure_id == ctx.exposure_id
                @test body.exposure_filename == "JC001-007"
            end
        end
    end

    @testset "200: experiment-only" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp")
                @test r.status == 200
                body = JSON3.read(String(r.body))
                @test body.experiment_id == ctx.experiment_id
                @test body.experiment_name == "test-exp"
                @test !haskey(body, :sample_id) || body.sample_id === nothing
            end
        end
    end

    @testset "200: id-form lookup returns names" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment_id=$(ctx.experiment_id)&sample_id=$(ctx.sample_id)")
                @test r.status == 200
                body = JSON3.read(String(r.body))
                @test body.experiment_name == "test-exp"
                @test body.sample_name == "S1"
            end
        end
    end

    @testset "404: missing experiment" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=nope"; status_exception=false)
                @test r.status == 404
                body = JSON3.read(String(r.body))
                @test body.error == "not_found"
                @test body.missing == "experiment"
                @test body.missing_value == "nope"
            end
        end
    end

    @testset "404: missing sample (experiment_resolved present)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp&sample=nope"; status_exception=false)
                @test r.status == 404
                body = JSON3.read(String(r.body))
                @test body.error == "not_found"
                @test body.missing == "sample"
                @test body.missing_value == "nope"
                @test body.experiment_resolved.id == ctx.experiment_id
                @test body.experiment_resolved.name == "test-exp"
            end
        end
    end

    @testset "404: missing exposure (experiment_resolved + sample_resolved present)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp&sample=S1&exposure=nope"; status_exception=false)
                @test r.status == 404
                body = JSON3.read(String(r.body))
                @test body.error == "not_found"
                @test body.missing == "exposure"
                @test body.missing_value == "nope"
                @test body.experiment_resolved.id == ctx.experiment_id
                @test body.sample_resolved.id == ctx.sample_id
            end
        end
    end

    @testset "400: malformed numeric param returns 400, not 500" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment_id=abc"; status_exception=false)
                @test r.status == 400
                body = JSON3.read(String(r.body))
                @test body.error == "invalid_id"
            end
        end
    end

    @testset "404: stale-name regression — rename experiment, old name 404s" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            DBInterface.execute(ctx.db, "UPDATE experiments SET name = 'test-exp-renamed' WHERE id = ?",
                                [ctx.experiment_id])
            with_test_server(ctx.db) do port, base
                r1 = HTTP.get("$base/api/resolve?experiment=test-exp"; status_exception=false)
                @test r1.status == 404
                r2 = HTTP.get("$base/api/resolve?experiment=test-exp-renamed")
                @test r2.status == 200
            end
        end
    end

    @testset "404: id-form for deleted sample (cold-mount Zustand staleness)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            # FK enforcement is on (open_db sets PRAGMA foreign_keys = ON);
            # exposures.sample_id REFERENCES samples(id) without ON DELETE.
            # Disable FKs around the synthetic delete to simulate the
            # "sample id no longer exists" cold-mount state without
            # tearing down dependent rows.
            DBInterface.execute(ctx.db, "PRAGMA foreign_keys = OFF")
            DBInterface.execute(ctx.db, "DELETE FROM samples WHERE id = ?", [ctx.sample_id])
            DBInterface.execute(ctx.db, "PRAGMA foreign_keys = ON")
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment_id=$(ctx.experiment_id)&sample_id=$(ctx.sample_id)";
                             status_exception=false)
                @test r.status == 404
                body = JSON3.read(String(r.body))
                @test body.missing == "sample"
            end
        end
    end
end
```

Add the new test file to `packages/HimalayaUI/test/runtests.jl`. The include must come AFTER `test_http.jl` (line 12) and `test_route_response_shapes.jl` (line 33) because the helpers are defined there. A safe placement is alongside `test_picker_routes.jl` (line 38). Add:

```julia
include("test_routes_resolve.jl")
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|fail|did not pass" /tmp/jl-test.out
```

Expected: 404s on all calls (the route doesn't exist).

- [ ] **Step 3a: Hoist a shared `_setup_for_resolve` helper into `test_route_response_shapes.jl`**

In `packages/HimalayaUI/test/test_route_response_shapes.jl`, immediately after the existing `_setup_analyzed_exposure` definition (around line 43), add:

```julia
"""
Like `_setup_analyzed_exposure` but UPDATEs the rows to known slug-resolvable
names ("test-exp" / "S1" / "JC001-007") and captures `experiment_id` for tests
that need it. Used by `test_routes_resolve.jl` and the resolve-shape rows below.
"""
function _setup_for_resolve(tmp::String)
    ctx = _setup_analyzed_exposure(tmp)
    DBInterface.execute(ctx.db, "UPDATE experiments SET name = 'test-exp'")
    DBInterface.execute(ctx.db, "UPDATE samples SET name = 'S1' WHERE id = ?",
                        [ctx.sample_id])
    DBInterface.execute(ctx.db, "UPDATE exposures SET filename = 'JC001-007' WHERE id = ?",
                        [ctx.exposure_id])
    exp_row = Tables.rowtable(DBInterface.execute(ctx.db,
        "SELECT id FROM experiments LIMIT 1"))[1]
    return (db = ctx.db,
            experiment_id = Int(exp_row.id),
            sample_id = ctx.sample_id,
            exposure_id = ctx.exposure_id)
end
```

This hoist is required because `runtests.jl` includes `test_route_response_shapes.jl` (line 33) BEFORE `test_routes_resolve.jl`, and the resolve-shape contract test (Task 4) lives inside `test_route_response_shapes.jl`. Defining the helper there makes it available to both files.

- [ ] **Step 3: Implement the route**

Create `packages/HimalayaUI/src/routes_resolve.jl`:

```julia
using HTTP, JSON3, DBInterface, Tables, Oxygen

# ─────────────────────────────────────────────────────────────────────────────
# Spec §3.1 — `/api/resolve` slug-to-id resolver.
#
# Read-only. No writes; no with_idempotency; no apply_event!; no SSE; no
# client_op_id; no pendingDeferreds participation. Three SELECTs and a
# response. Queue-orthogonal by design — confirmed against
# docs/mutation-queue.md. Future maintainers tempted to extend it with a
# write path: please don't; add a sibling endpoint instead.
#
# Tiebreaker for non-unique experiment names: ORDER BY id ASC LIMIT 1
# (deterministic across runs). samples and exposures are uniqueness-
# enforced post-#88 / Task 1.
# ─────────────────────────────────────────────────────────────────────────────

function _json(status::Int, body)
    HTTP.Response(status, ["Content-Type" => "application/json"], JSON3.write(body))
end

function _has_param(params, name::String)
    haskey(params, name) && !isempty(params[name])
end

function _safe_str(v)
    # Tables.rowtable returns `missing` for SQL NULL. Coerce to "" so the
    # response shape stays string-typed; the caller may post-process.
    ismissing(v) ? "" : String(v)
end

function _parse_id_or_400(s::String, field::String)
    n = tryparse(Int, s)
    n === nothing && return _json(400, Dict(:error => "invalid_id", :field => field))
    return n
end

function register_resolve_routes!()
    @get "/api/resolve" function(req::HTTP.Request)
        db = current_db()
        params = HTTP.queryparams(HTTP.URI(req.target))

        # Mutual-exclusion check: same-entity name+id collision is 400.
        for (n, i) in (("experiment", "experiment_id"),
                       ("sample",     "sample_id"),
                       ("exposure",   "exposure_id"))
            if _has_param(params, n) && _has_param(params, i)
                return _json(400, Dict(:error => "ambiguous_params"))
            end
        end

        # Resolve experiment. NULL-name rows are treated as "no canonical
        # slug" → 404. The frontend has no path to construct a URL for an
        # experiment without a name; rejecting them here keeps the round
        # trip well-formed.
        exp_row = nothing
        if _has_param(params, "experiment")
            name = params["experiment"]
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, name FROM experiments WHERE name = ? ORDER BY id ASC LIMIT 1", [name]))
            isempty(rows) && return _json(404, Dict(
                :error => "not_found", :missing => "experiment", :missing_value => name))
            exp_row = (id=Int(rows[1].id), name=_safe_str(rows[1].name))
        elseif _has_param(params, "experiment_id")
            id_or_resp = _parse_id_or_400(params["experiment_id"], "experiment_id")
            id_or_resp isa HTTP.Response && return id_or_resp
            id = id_or_resp::Int
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, name FROM experiments WHERE id = ? LIMIT 1", [id]))
            isempty(rows) && return _json(404, Dict(
                :error => "not_found", :missing => "experiment",
                :missing_value => string(id)))
            ismissing(rows[1].name) && return _json(404, Dict(
                :error => "not_found", :missing => "experiment",
                :missing_value => string(id),
                :reason => "experiment has no canonical name"))
            exp_row = (id=Int(rows[1].id), name=_safe_str(rows[1].name))
        else
            return _json(400, Dict(:error => "missing_experiment"))
        end

        # Resolve sample (optional).
        sample_row = nothing
        if _has_param(params, "sample") || _has_param(params, "sample_id")
            if _has_param(params, "sample")
                name = params["sample"]
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, name FROM samples WHERE experiment_id = ? AND name = ? LIMIT 1",
                    [exp_row.id, name]))
                isempty(rows) && return _json(404, Dict(
                    :error => "not_found", :missing => "sample", :missing_value => name,
                    :experiment_resolved => Dict(:id => exp_row.id, :name => exp_row.name)))
                sample_row = (id=Int(rows[1].id), name=_safe_str(rows[1].name))
            else
                id_or_resp = _parse_id_or_400(params["sample_id"], "sample_id")
                id_or_resp isa HTTP.Response && return id_or_resp
                id = id_or_resp::Int
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, name FROM samples WHERE id = ? AND experiment_id = ? LIMIT 1",
                    [id, exp_row.id]))
                isempty(rows) && return _json(404, Dict(
                    :error => "not_found", :missing => "sample", :missing_value => string(id),
                    :experiment_resolved => Dict(:id => exp_row.id, :name => exp_row.name)))
                sample_row = (id=Int(rows[1].id), name=_safe_str(rows[1].name))
            end
        end

        # Resolve exposure (optional, requires sample).
        exposure_row = nothing
        if _has_param(params, "exposure") || _has_param(params, "exposure_id")
            sample_row === nothing && return _json(400, Dict(:error => "missing_sample"))
            if _has_param(params, "exposure")
                name = params["exposure"]
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, filename FROM exposures WHERE sample_id = ? AND filename = ? LIMIT 1",
                    [sample_row.id, name]))
                isempty(rows) && return _json(404, Dict(
                    :error => "not_found", :missing => "exposure", :missing_value => name,
                    :experiment_resolved => Dict(:id => exp_row.id, :name => exp_row.name),
                    :sample_resolved     => Dict(:id => sample_row.id, :name => sample_row.name)))
                exposure_row = (id=Int(rows[1].id), filename=_safe_str(rows[1].filename))
            else
                id_or_resp = _parse_id_or_400(params["exposure_id"], "exposure_id")
                id_or_resp isa HTTP.Response && return id_or_resp
                id = id_or_resp::Int
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, filename FROM exposures WHERE id = ? AND sample_id = ? LIMIT 1",
                    [id, sample_row.id]))
                isempty(rows) && return _json(404, Dict(
                    :error => "not_found", :missing => "exposure", :missing_value => string(id),
                    :experiment_resolved => Dict(:id => exp_row.id, :name => exp_row.name),
                    :sample_resolved     => Dict(:id => sample_row.id, :name => sample_row.name)))
                exposure_row = (id=Int(rows[1].id), filename=_safe_str(rows[1].filename))
            end
        end

        # Build success response.
        body = Dict{Symbol,Any}(
            :experiment_id   => exp_row.id,
            :experiment_name => exp_row.name,
        )
        if sample_row !== nothing
            body[:sample_id]   = sample_row.id
            body[:sample_name] = sample_row.name
        end
        if exposure_row !== nothing
            body[:exposure_id]       = exposure_row.id
            body[:exposure_filename] = exposure_row.filename
        end
        _json(200, body)
    end
end
```

Add `include("routes_resolve.jl")` to `packages/HimalayaUI/src/HimalayaUI.jl` near the other route includes.

In `packages/HimalayaUI/src/server.jl`, add `register_resolve_routes!()` to `register_routes!`. Find the existing list:

```julia
    register_picker_routes!()
end
```

Insert before the closing `end`:

```julia
    register_picker_routes!()
    register_resolve_routes!()
end
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|fail|did not pass" /tmp/jl-test.out
```

Expected: all assertions pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_resolve.jl packages/HimalayaUI/src/HimalayaUI.jl packages/HimalayaUI/src/server.jl packages/HimalayaUI/test/test_routes_resolve.jl packages/HimalayaUI/test/runtests.jl
git commit -m "feat(routes): add /api/resolve slug-to-id endpoint

Read-only. Translates (experiment, sample, exposure) slugs to ids in
one round trip. 200 with names+ids; 404 with missing entity + already-
resolved prefix; 400 on same-entity name+id collision."
```

---

## Task 3: Resolve route — tiebreaker + 400 ambiguity test

**Files:**
- Modify: `packages/HimalayaUI/test/test_routes_resolve.jl`

- [ ] **Step 1: Write the failing tests**

Append to the existing `@testset "GET /api/resolve"`:

```julia
    @testset "400: ambiguous params (experiment + experiment_id)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp&experiment_id=$(ctx.experiment_id)";
                             status_exception=false)
                @test r.status == 400
                body = JSON3.read(String(r.body))
                @test body.error == "ambiguous_params"
            end
        end
    end

    @testset "400: ambiguous params (sample + sample_id)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp&sample=S1&sample_id=$(ctx.sample_id)";
                             status_exception=false)
                @test r.status == 400
            end
        end
    end

    @testset "200: mixed name+id across entities is allowed" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                # name-form experiment + id-form sample is fine.
                r = HTTP.get("$base/api/resolve?experiment=test-exp&sample_id=$(ctx.sample_id)")
                @test r.status == 200
                body = JSON3.read(String(r.body))
                @test body.experiment_id == ctx.experiment_id
                @test body.sample_id == ctx.sample_id
            end
        end
    end

    @testset "tiebreaker: duplicate experiment names → lowest id wins" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            # Insert a second experiment with the same name.
            res = DBInterface.execute(ctx.db,
                "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('test-exp','/p','/d','/a')")
            second_id = Int(DBInterface.lastrowid(res))
            @test second_id > ctx.experiment_id

            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp")
                @test r.status == 200
                body = JSON3.read(String(r.body))
                # Lowest id wins deterministically.
                @test body.experiment_id == ctx.experiment_id
                @test body.experiment_id < second_id
            end
        end
    end
```

- [ ] **Step 2: Run tests to verify they pass**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|fail|did not pass" /tmp/jl-test.out
```

Expected: all assertions pass — Task 2's implementation already includes the 400 path and `ORDER BY id ASC LIMIT 1` tiebreaker.

If a test fails, fix Task 2's implementation rather than the test.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/test/test_routes_resolve.jl
git commit -m "test(routes): pin /api/resolve tiebreaker + ambiguity contract

400 on same-entity name+id collision; 200 on cross-entity mixed forms;
deterministic ORDER BY id ASC LIMIT 1 on duplicate experiment names."
```

---

## Task 4: Response-shape contract test rows

**Files:**
- Modify: `packages/HimalayaUI/test/test_route_response_shapes.jl`

- [ ] **Step 1: Open the existing contract file and find the `@testset` block**

```bash
grep -n "@testset" packages/HimalayaUI/test/test_route_response_shapes.jl | head -10
```

The file should have a top-level `@testset "route response shapes"` containing rows that pin the exact key set on each route's response body. Add three new rows.

- [ ] **Step 2: Add the rows**

Inside the existing top-level `@testset`, append. `_setup_for_resolve` was added to `test_route_response_shapes.jl` itself in Task 2 Step 3a, so it's available here without import.

```julia
    @testset "GET /api/resolve 200 (experiment+sample+exposure)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp&sample=S1&exposure=JC001-007")
                body = JSON3.read(String(r.body))
                # Key set frozen.
                @test Set(keys(body)) == Set([
                    :experiment_id, :experiment_name,
                    :sample_id, :sample_name,
                    :exposure_id, :exposure_filename,
                ])
            end
        end
    end

    @testset "GET /api/resolve 404 (missing exposure)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp&sample=S1&exposure=nope";
                             status_exception=false)
                body = JSON3.read(String(r.body))
                @test Set(keys(body)) == Set([
                    :error, :missing, :missing_value,
                    :experiment_resolved, :sample_resolved,
                ])
                @test Set(keys(body.experiment_resolved)) == Set([:id, :name])
                @test Set(keys(body.sample_resolved))     == Set([:id, :name])
            end
        end
    end

    @testset "GET /api/resolve 400 (ambiguous params)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp&experiment_id=$(ctx.experiment_id)";
                             status_exception=false)
                body = JSON3.read(String(r.body))
                @test Set(keys(body)) == Set([:error])
            end
        end
    end
```

- [ ] **Step 3: Run tests to verify they pass**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|fail|did not pass" /tmp/jl-test.out
```

Expected: all assertions pass.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/test/test_route_response_shapes.jl
git commit -m "test(contract): pin /api/resolve 200/400/404 response shapes"
```

---

## Task 5: SPA catch-all in `server.jl`

**Files:**
- Modify: `packages/HimalayaUI/src/server.jl`
- Create: `packages/HimalayaUI/test/test_spa_fallback.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/test_spa_fallback.jl`:

```julia
using Test, HTTP, SQLite, DBInterface
using HimalayaUI

# Helpers are available via runtests.jl include order — see test_routes_resolve.jl.

# Spec §3.2 — SPA catch-all serves index.html for any non-API, non-asset
# path so deep-link URLs like /inspect/exp/sample reach the frontend.
# /api/* always returns 404 to prevent unregistered API typos from being
# masked as 200 HTML responses.

@testset "SPA fallback" begin
    mktempdir() do tmp
        # Synthesize a minimal dist directory.
        dist = joinpath(tmp, "dist")
        mkpath(dist)
        index_html = "<!doctype html><html><body>shell</body></html>"
        write(joinpath(dist, "index.html"), index_html)
        write(joinpath(dist, "asset.png"), b"\x89PNG\r\n\x1a\nfake")

        # Point the server at the synthetic dist.
        prev = get(ENV, "HIMALAYA_FRONTEND_DIST", nothing)
        ENV["HIMALAYA_FRONTEND_DIST"] = dist
        try
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                @testset "single-segment unknown path returns index.html" begin
                    r = HTTP.get("$base/foo"; status_exception=false)
                    @test r.status == 200
                    @test occursin("shell", String(r.body))
                    @test HTTP.header(r, "Cache-Control") == "no-store"
                    @test occursin("text/html", HTTP.header(r, "Content-Type"))
                end

                @testset "multi-segment unknown path returns index.html (pins /** syntax)" begin
                    r = HTTP.get("$base/inspect/exp/sample"; status_exception=false)
                    @test r.status == 200
                    @test occursin("shell", String(r.body))
                end

                @testset "/api/* unregistered returns 404 (does not fall through)" begin
                    r = HTTP.get("$base/api/typo-not-registered"; status_exception=false)
                    @test r.status == 404
                end

                @testset "asset path served by dynamicfiles, not catch-all" begin
                    r = HTTP.get("$base/asset.png"; status_exception=false)
                    @test r.status == 200
                    # dynamicfiles doesn't add no-store; the catch-all does.
                    @test HTTP.header(r, "Cache-Control") != "no-store"
                end
            end
        finally
            if prev === nothing
                delete!(ENV, "HIMALAYA_FRONTEND_DIST")
            else
                ENV["HIMALAYA_FRONTEND_DIST"] = prev
            end
        end
    end
end
```

Add to `packages/HimalayaUI/test/runtests.jl`:

```julia
include("test_spa_fallback.jl")
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|fail|did not pass" /tmp/jl-test.out
```

Expected: failures for the multi-segment path and unknown `/api/*` (without the catch-all, both 404).

- [ ] **Step 3: Add the catch-all to `server.jl`**

In `packages/HimalayaUI/src/server.jl`, find the dist-mount block at `register_routes!` (~line 31–35):

```julia
    dist_dir = get(ENV, "HIMALAYA_FRONTEND_DIST",
                   joinpath(pkgdir(HimalayaUI), "frontend", "dist"))
    if isdir(dist_dir)
        Oxygen.dynamicfiles(dist_dir, "/")
    end
```

Replace with:

```julia
    dist_dir = get(ENV, "HIMALAYA_FRONTEND_DIST",
                   joinpath(pkgdir(HimalayaUI), "frontend", "dist"))
    if isdir(dist_dir)
        Oxygen.dynamicfiles(dist_dir, "/")

        # SPA catch-all (spec §3.2). Doublestar `/**` is HTTP.jl's multi-
        # segment wildcard (Handlers.jl:174,219–227); single-conditional
        # captures like `{rest:.*}` only match one segment. The `api/`
        # guard is critical — without it an unregistered API route would
        # fall through to index.html and mask 404s as 200s.
        @get "/**" function(req::HTTP.Request)
            path = HTTP.URI(req.target).path
            rest = lstrip(path, '/')
            startswith(rest, "api/") && return HTTP.Response(404, "Not found")
            return HTTP.Response(200,
                ["Content-Type" => "text/html; charset=utf-8",
                 "Cache-Control" => "no-store"],
                read(joinpath(dist_dir, "index.html")))
        end
    end
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|fail|did not pass" /tmp/jl-test.out
```

Expected: all assertions pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/server.jl packages/HimalayaUI/test/test_spa_fallback.jl packages/HimalayaUI/test/runtests.jl
git commit -m "feat(server): SPA catch-all for slug URLs

Adds @get \"/**\" inside the dist-mount block. Multi-segment paths reach
index.html; /api/* unregistered returns 404 (no fall-through to HTML)."
```

---

## Task 6: Frontend Zustand state additions

**Files:**
- Modify: `frontend/src/state.ts`
- Create: `frontend/test/staleUrlState.test.ts`

- [ ] **Step 1: Write the failing test**

Create `frontend/test/staleUrlState.test.ts`:

```ts
import { describe, it, expect, beforeEach } from "vitest";
import { useAppState } from "../src/state";

// Spec §4.4 + §6 — state.ts additions for permalink URL handling.

describe("Zustand state — permalink slots", () => {
  beforeEach(() => {
    // Reset relevant slots without nuking the whole store (preserves persisted
    // username so other tests don't fight unexpectedly).
    useAppState.setState({
      staleUrlContext: null,
      resolving: false,
      activeExperimentId: undefined,
      activeSampleId: undefined,
      activeExposureId: undefined,
      navModalOpen: false,
      navModalStep: "experiment",
    });
  });

  it("setStaleUrlContext stores and clears", () => {
    const ctx = {
      kind: "not_found" as const,
      missing: "experiment" as const,
      missing_value: "lipid-typo",
      experiment_resolved: undefined,
      sample_resolved: undefined,
    };
    useAppState.getState().setStaleUrlContext(ctx);
    expect(useAppState.getState().staleUrlContext).toEqual(ctx);
    useAppState.getState().setStaleUrlContext(null);
    expect(useAppState.getState().staleUrlContext).toBeNull();
  });

  it("setResolving toggles", () => {
    expect(useAppState.getState().resolving).toBe(false);
    useAppState.getState().setResolving(true);
    expect(useAppState.getState().resolving).toBe(true);
  });

  it.each([
    ["setActiveExperiment", (id: number) => useAppState.getState().setActiveExperiment(id)],
    ["setActiveSample",     (id: number) => useAppState.getState().setActiveSample(id)],
    ["setActiveExposure",   (id: number) => useAppState.getState().setActiveExposure(id)],
    ["setActivePage",       () => useAppState.getState().setActivePage("inspect")],
  ])("%s clears staleUrlContext", (_label, fn) => {
    useAppState.getState().setStaleUrlContext({
      kind: "unknown_path", raw: "/foo/bar",
    });
    fn(42);
    expect(useAppState.getState().staleUrlContext).toBeNull();
  });

  it("setActive* does NOT clear resolving", () => {
    useAppState.getState().setResolving(true);
    useAppState.getState().setActiveExperiment(5);
    expect(useAppState.getState().resolving).toBe(true);
  });

  it("recoverFromStaleUrl row 1 (experiment): clears stale, opens modal at experiment step", () => {
    useAppState.getState().setStaleUrlContext({
      kind: "not_found", missing: "experiment", missing_value: "x",
      experiment_resolved: undefined, sample_resolved: undefined,
    });
    useAppState.getState().recoverFromStaleUrl({
      step: "experiment",
      experimentId: undefined,
      sampleId: undefined,
    });
    const s = useAppState.getState();
    expect(s.staleUrlContext).toBeNull();
    expect(s.activeSampleId).toBeUndefined();
    expect(s.activeExposureId).toBeUndefined();
    expect(s.navModalOpen).toBe(true);
    expect(s.navModalStep).toBe("experiment");
  });

  it("recoverFromStaleUrl row 2 (sample): sets experimentId, opens modal at sample step", () => {
    useAppState.getState().recoverFromStaleUrl({
      step: "sample",
      experimentId: 17,
      sampleId: undefined,
    });
    const s = useAppState.getState();
    expect(s.staleUrlContext).toBeNull();
    expect(s.activeExperimentId).toBe(17);
    expect(s.activeSampleId).toBeUndefined();
    expect(s.navModalOpen).toBe(true);
    expect(s.navModalStep).toBe("sample");
  });

  it("recoverFromStaleUrl row 3 (exposure): preserves sample, openModal=false", () => {
    useAppState.getState().recoverFromStaleUrl({
      step: "sample",
      experimentId: 17,
      sampleId: 42,
      openModal: false,
    });
    const s = useAppState.getState();
    expect(s.staleUrlContext).toBeNull();
    expect(s.activeExperimentId).toBe(17);
    expect(s.activeSampleId).toBe(42);
    expect(s.activeExposureId).toBeUndefined();
    expect(s.navModalOpen).toBe(false);
    expect(s.navModalStep).toBe("sample");
  });

  it("staleUrlContext and resolving are NOT in the persisted slice", () => {
    // Touch them, then read back the Zustand persist `partialize` output.
    useAppState.getState().setResolving(true);
    useAppState.getState().setStaleUrlContext({ kind: "unknown_path", raw: "/x" });
    const persisted = JSON.parse(localStorage.getItem("himalaya-ui:state") ?? "{}");
    expect(persisted.state?.resolving).toBeUndefined();
    expect(persisted.state?.staleUrlContext).toBeUndefined();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd packages/HimalayaUI/frontend
node_modules/.bin/vitest run test/staleUrlState.test.ts
```

Expected: failures — slots and actions don't exist yet.

- [ ] **Step 3: Add the types and slots to `state.ts`**

Near the top of `frontend/src/state.ts` (after the existing imports, before the `AppState` interface):

```ts
export type StaleUrlContext =
  | { kind: "not_found"; missing: "experiment" | "sample" | "exposure";
      missing_value: string;
      experiment_resolved: { id: number; name: string } | undefined;
      sample_resolved:     { id: number; name: string } | undefined }
  | { kind: "unknown_path"; raw: string };

export type RecoverOpts = {
  step: NavModalStep;
  experimentId: number | undefined;
  sampleId: number | undefined;
  openModal?: boolean;          // default true; row "exposure" passes false
};
```

In the `AppState` interface, add the new fields (alongside the ephemeral block):

```ts
  // ephemeral
  staleUrlContext: StaleUrlContext | null;
  resolving: boolean;
```

And the new actions:

```ts
  setStaleUrlContext: (ctx: StaleUrlContext | null) => void;
  setResolving: (v: boolean) => void;
  recoverFromStaleUrl: (opts: RecoverOpts) => void;
```

In the store definition, add the initial values (alongside the other ephemeral defaults):

```ts
        staleUrlContext: null,
        resolving: false,
```

Add the action implementations:

```ts
        setStaleUrlContext: (staleUrlContext) => set({ staleUrlContext }),
        setResolving: (resolving) => set({ resolving }),
        recoverFromStaleUrl: (opts) =>
          set((s) => ({
            staleUrlContext: null,
            activeExperimentId: opts.experimentId ?? s.activeExperimentId,
            activeSampleId: opts.sampleId ?? undefined,
            activeExposureId: undefined,
            navModalOpen: opts.openModal ?? true,
            navModalStep: opts.step,
          })),
```

Update the existing `setActive*` actions to also clear `staleUrlContext`. Find:

```ts
        setActiveExperiment: (activeExperimentId) =>
          set({
            activeExperimentId,
            activeSampleId: undefined,
            activeExposureId: undefined,
          }),
        setActiveSample: (activeSampleId) =>
          set({ activeSampleId, activeExposureId: undefined }),
        setActiveExposure: (activeExposureId) => set({ activeExposureId }),
        // ...
        setActivePage: (activePage) => set({ activePage }),
```

Replace with:

```ts
        setActiveExperiment: (activeExperimentId) =>
          set({
            activeExperimentId,
            activeSampleId: undefined,
            activeExposureId: undefined,
            staleUrlContext: null,
          }),
        setActiveSample: (activeSampleId) =>
          set({ activeSampleId, activeExposureId: undefined, staleUrlContext: null }),
        setActiveExposure: (activeExposureId) =>
          set({ activeExposureId, staleUrlContext: null }),
        // ...
        setActivePage: (activePage) => set({ activePage, staleUrlContext: null }),
```

Verify `partialize` does NOT include `staleUrlContext` or `resolving` (it shouldn't — both are ephemeral).

- [ ] **Step 4: Run tests to verify they pass**

```bash
node_modules/.bin/vitest run test/staleUrlState.test.ts
```

Expected: all assertions pass.

- [ ] **Step 5: Commit**

```bash
git add frontend/src/state.ts frontend/test/staleUrlState.test.ts
git commit -m "feat(state): add staleUrlContext, resolving, recoverFromStaleUrl

Per spec §4.4 + §6. setActive* actions also clear staleUrlContext.
Both new slots are ephemeral (not persisted)."
```

---

## Task 7: `parseLocation` parser

**Files:**
- Create: `frontend/src/lib/url/parseLocation.ts`
- Create: `frontend/test/parseLocation.test.ts`

- [ ] **Step 1: Write the failing test**

Create `frontend/test/parseLocation.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import { parseLocation } from "../src/lib/url/parseLocation";

// Spec §4.1 — discriminated-union output keyed by `kind`.

describe("parseLocation", () => {
  it("/ → root", () => {
    expect(parseLocation("/", "")).toEqual({ kind: "root" });
    expect(parseLocation("", "")).toEqual({ kind: "root" });
  });

  it("/index → empty index", () => {
    expect(parseLocation("/index", "")).toEqual({
      kind: "index", experiment: undefined, sample: undefined,
    });
  });

  it("/index/<exp> → experiment chosen", () => {
    expect(parseLocation("/index/lipid-screen", "")).toEqual({
      kind: "index", experiment: "lipid-screen", sample: undefined,
    });
  });

  it("/index/<exp>/<sample> → full", () => {
    expect(parseLocation("/index/lipid-screen/JC001", "")).toEqual({
      kind: "index", experiment: "lipid-screen", sample: "JC001",
    });
  });

  it("/inspect with ?exposure=", () => {
    expect(parseLocation("/inspect/lipid/JC001", "?exposure=JC001-007")).toEqual({
      kind: "inspect", experiment: "lipid", sample: "JC001", exposure: "JC001-007",
    });
  });

  it("/inspect ignores ?exposure when sample missing", () => {
    expect(parseLocation("/inspect/lipid", "?exposure=JC001-007")).toEqual({
      kind: "inspect", experiment: "lipid", sample: undefined, exposure: undefined,
    });
  });

  it("/compare (legacy list) → compare:list", () => {
    expect(parseLocation("/compare", "")).toEqual({ kind: "compare", view: "list" });
  });

  it("/compare/all → compare:list", () => {
    expect(parseLocation("/compare/all", "")).toEqual({ kind: "compare", view: "list" });
  });

  it("/compare/all/42 → compare:list (ComparePage handles numeric id)", () => {
    expect(parseLocation("/compare/all/42", "")).toEqual({ kind: "compare", view: "list" });
  });

  it("/compare/all/42/edit → compare:list", () => {
    expect(parseLocation("/compare/all/42/edit", "")).toEqual({ kind: "compare", view: "list" });
  });

  it("/experiments/17/compare → compare:list", () => {
    expect(parseLocation("/experiments/17/compare", "")).toEqual({ kind: "compare", view: "list" });
  });

  it("/experiments/17/compare/42 → compare:list", () => {
    expect(parseLocation("/experiments/17/compare/42", "")).toEqual({ kind: "compare", view: "list" });
  });

  it("/foo/bar → stale", () => {
    const r = parseLocation("/foo/bar", "?x=1");
    expect(r.kind).toBe("stale");
    if (r.kind === "stale") expect(r.raw).toBe("/foo/bar?x=1");
  });

  it("decodes percent-encoded slugs", () => {
    expect(parseLocation("/index/lipid%20screen/JC%20001", "")).toEqual({
      kind: "index", experiment: "lipid screen", sample: "JC 001",
    });
  });

  it("trailing slash tolerant", () => {
    expect(parseLocation("/index/", "")).toEqual({
      kind: "index", experiment: undefined, sample: undefined,
    });
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

```bash
node_modules/.bin/vitest run test/parseLocation.test.ts
```

Expected: module-not-found error.

- [ ] **Step 3: Implement the parser**

Create `frontend/src/lib/url/parseLocation.ts`:

```ts
// Spec §4.1 — pure URL parser. No side effects.

// Compare URLs are owned by ComparePage / ComparePageEdit (which use
// useParams + useNavigate). useStateFromUrl only needs to know "this is
// a compare path so set activePage='compare' and don't 404 it." Numeric
// ids are not tracked in Zustand.
export type ParsedUrl =
  | { kind: "root" }
  | { kind: "index";   experiment: string | undefined; sample: string | undefined }
  | { kind: "inspect"; experiment: string | undefined; sample: string | undefined; exposure: string | undefined }
  | { kind: "compare"; view: "list" }
  | { kind: "stale"; raw: string };

const safeDecode = (s: string): string => {
  try { return decodeURIComponent(s); } catch { return s; }
};

export function parseLocation(pathname: string, search: string): ParsedUrl {
  // Normalize trailing slash and leading slash.
  const segs = pathname.replace(/\/+$/, "").split("/").filter(Boolean);
  if (segs.length === 0) return { kind: "root" };

  const [page, a, b, c] = segs;
  const decode = (s: string | undefined): string | undefined =>
    s === undefined ? undefined : safeDecode(s);

  if (page === "index" && segs.length <= 3) {
    return { kind: "index", experiment: decode(a), sample: decode(b) };
  }

  if (page === "inspect" && segs.length <= 3) {
    const experiment = decode(a);
    const sample = decode(b);
    let exposure: string | undefined = undefined;
    if (experiment !== undefined && sample !== undefined && search) {
      const params = new URLSearchParams(search.startsWith("?") ? search.slice(1) : search);
      const e = params.get("exposure");
      if (e !== null) exposure = e;
    }
    return { kind: "inspect", experiment, sample, exposure };
  }

  // Compare URLs in the actual app live under two roots:
  //   /experiments/:eid/compare(/new|/:id|/:id/edit)
  //   /compare/all(/new|/:id|/:id/edit)
  // The hooks don't track Compare numeric ids in Zustand — ComparePage
  // handles its own URL via useNavigate/useParams. Recognizing the shape
  // here just lets useStateFromUrl set activePage="compare" without
  // falsely flagging the path as "stale".
  const isExperimentCompare = page === "experiments" && segs[2] === "compare";
  const isCompareAll = page === "compare" && a === "all";
  const isLegacyCompareList = page === "compare" && segs.length === 1;

  if (isExperimentCompare || isCompareAll || isLegacyCompareList) {
    return { kind: "compare", view: "list" };
  }

  return { kind: "stale", raw: pathname + (search ?? "") };
}
```

(`c` and the third positional `b/c` are unused for non-`edit` shapes; lint may complain. The shape is intentional — see the destructure-then-check pattern.)

- [ ] **Step 4: Run test to verify it passes**

```bash
node_modules/.bin/vitest run test/parseLocation.test.ts
```

Expected: all assertions pass.

- [ ] **Step 5: Commit**

```bash
git add frontend/src/lib/url/parseLocation.ts frontend/test/parseLocation.test.ts
git commit -m "feat(url): pure parseLocation parser

Discriminated union by kind. Handles all URL shapes from spec §2 plus
percent-encoded slugs and trailing slashes; unknown paths produce
{ kind: 'stale', raw: <full path+search> }."
```

---

## Task 8: `api.ts` — resolve fetcher + types

**Files:**
- Modify: `frontend/src/api.ts`

- [ ] **Step 1: Add the types and fetcher**

In `frontend/src/api.ts`, add near the other interface declarations:

```ts
export interface ResolveSuccess {
  experiment_id: number;
  experiment_name: string;
  sample_id: number | undefined;
  sample_name: string | undefined;
  exposure_id: number | undefined;
  exposure_filename: string | undefined;
}

export interface ResolveError404 {
  error: "not_found";
  missing: "experiment" | "sample" | "exposure";
  missing_value: string;
  experiment_resolved: { id: number; name: string } | undefined;
  sample_resolved:     { id: number; name: string } | undefined;
}

export interface ResolveError400 {
  error: "ambiguous_params" | "missing_experiment" | "missing_sample";
}

export interface ResolveQuery {
  experiment?: string;
  sample?: string;
  exposure?: string;
  experiment_id?: number;
  sample_id?: number;
  exposure_id?: number;
}

/**
 * Calls /api/resolve. Returns either a `ResolveSuccess` (200), a
 * `ResolveError404`, or a `ResolveError400`. Caller distinguishes on
 * `error` field. Read-only — no auth headers.
 *
 * `signal` exposed for AbortController; caller is responsible for
 * the origin-tag race check (§4.2).
 */
export async function resolve(
  q: ResolveQuery,
  signal: AbortSignal | undefined = undefined,
): Promise<ResolveSuccess | ResolveError404 | ResolveError400> {
  const params = new URLSearchParams();
  if (q.experiment !== undefined) params.set("experiment", q.experiment);
  if (q.sample !== undefined) params.set("sample", q.sample);
  if (q.exposure !== undefined) params.set("exposure", q.exposure);
  if (q.experiment_id !== undefined) params.set("experiment_id", String(q.experiment_id));
  if (q.sample_id !== undefined) params.set("sample_id", String(q.sample_id));
  if (q.exposure_id !== undefined) params.set("exposure_id", String(q.exposure_id));

  const res = await fetch(`/api/resolve?${params.toString()}`, { signal });
  const body = await res.json();
  return body as ResolveSuccess | ResolveError404 | ResolveError400;
}
```

- [ ] **Step 2: Verify the project still compiles**

```bash
cd packages/HimalayaUI/frontend
node_modules/.bin/tsc --noEmit
```

Expected: no new errors. (Existing errors unrelated to this change are fine.)

- [ ] **Step 3: Commit**

```bash
git add frontend/src/api.ts
git commit -m "feat(api): typed resolve() fetcher + 200/400/404 types"
```

---

## Task 9: `<ResolvingFallback>` component

**Files:**
- Create: `frontend/src/components/ResolvingFallback.tsx`

- [ ] **Step 1: Implement the component**

Create `frontend/src/components/ResolvingFallback.tsx`:

```tsx
// Spec §4.2 — minimal sub-paint flicker mask. Renders inside <main>
// while a resolve is in flight. No <Skeleton> from boneyard — that
// would itself flash on every navigation.

export function ResolvingFallback(): JSX.Element {
  return (
    <div
      data-testid="resolving"
      className="flex-1 min-h-0 flex flex-col"
      aria-busy="true"
      aria-label="Loading page"
    />
  );
}
```

- [ ] **Step 2: Verify it compiles**

```bash
cd packages/HimalayaUI/frontend
node_modules/.bin/tsc --noEmit
```

Expected: clean.

- [ ] **Step 3: Commit**

```bash
git add frontend/src/components/ResolvingFallback.tsx
git commit -m "feat(component): ResolvingFallback flicker mask"
```

---

## Task 10: `<StaleUrlPage>` component

**Files:**
- Create: `frontend/src/components/StaleUrlPage.tsx`
- Create: `frontend/test/StaleUrlPage.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `frontend/test/StaleUrlPage.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { StaleUrlPage } from "../src/components/StaleUrlPage";
import { useAppState } from "../src/state";
import type { StaleUrlContext } from "../src/state";

describe("StaleUrlPage", () => {
  beforeEach(() => {
    useAppState.setState({
      staleUrlContext: null,
      activeExperimentId: undefined,
      activeSampleId: undefined,
      navModalOpen: false,
      navModalStep: "experiment",
      activePage: "index",
    });
  });

  it("not_found:experiment — renders header, data-missing, CTA", () => {
    const ctx: StaleUrlContext = {
      kind: "not_found", missing: "experiment", missing_value: "lipid-typo",
      experiment_resolved: undefined, sample_resolved: undefined,
    };
    render(<StaleUrlPage staleUrlContext={ctx} />);
    expect(screen.getByTestId("stale-url-page")).toHaveAttribute("data-missing", "experiment");
    expect(screen.getByRole("heading")).toHaveTextContent(/Experiment 'lipid-typo' not found\./);
    fireEvent.click(screen.getByTestId("stale-url-cta"));
    expect(useAppState.getState().navModalOpen).toBe(true);
    expect(useAppState.getState().navModalStep).toBe("experiment");
    expect(useAppState.getState().staleUrlContext).toBeNull();
  });

  it("not_found:sample — preselects experiment and opens modal at sample step", () => {
    const ctx: StaleUrlContext = {
      kind: "not_found", missing: "sample", missing_value: "JC001",
      experiment_resolved: { id: 17, name: "lipid-screen" },
      sample_resolved: undefined,
    };
    render(<StaleUrlPage staleUrlContext={ctx} />);
    expect(screen.getByTestId("stale-url-page")).toHaveAttribute("data-missing", "sample");
    expect(screen.getByRole("heading")).toHaveTextContent(/Sample 'JC001' not found in 'lipid-screen'\./);
    fireEvent.click(screen.getByTestId("stale-url-cta"));
    const s = useAppState.getState();
    expect(s.activeExperimentId).toBe(17);
    expect(s.navModalOpen).toBe(true);
    expect(s.navModalStep).toBe("sample");
  });

  it("not_found:exposure — preserves sample, openModal=false (snap back)", () => {
    const ctx: StaleUrlContext = {
      kind: "not_found", missing: "exposure", missing_value: "JC001-099",
      experiment_resolved: { id: 17, name: "lipid-screen" },
      sample_resolved: { id: 42, name: "JC001" },
    };
    render(<StaleUrlPage staleUrlContext={ctx} />);
    expect(screen.getByTestId("stale-url-page")).toHaveAttribute("data-missing", "exposure");
    fireEvent.click(screen.getByTestId("stale-url-cta"));
    const s = useAppState.getState();
    expect(s.activeExperimentId).toBe(17);
    expect(s.activeSampleId).toBe(42);
    expect(s.navModalOpen).toBe(false);
  });

  it("unknown_path — Page not found, CTA dispatches setActivePage", () => {
    const ctx: StaleUrlContext = { kind: "unknown_path", raw: "/foo/bar" };
    render(<StaleUrlPage staleUrlContext={ctx} />);
    expect(screen.getByTestId("stale-url-page")).toHaveAttribute("data-missing", "path");
    expect(screen.getByRole("heading")).toHaveTextContent(/Page not found\./);
    fireEvent.click(screen.getByTestId("stale-url-cta"));
    expect(useAppState.getState().activePage).toBe("index");
    expect(useAppState.getState().staleUrlContext).toBeNull();
  });

  it("/ keypress triggers CTA", () => {
    const ctx: StaleUrlContext = {
      kind: "not_found", missing: "experiment", missing_value: "x",
      experiment_resolved: undefined, sample_resolved: undefined,
    };
    render(<StaleUrlPage staleUrlContext={ctx} />);
    fireEvent.keyDown(window, { key: "/" });
    expect(useAppState.getState().navModalOpen).toBe(true);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

```bash
node_modules/.bin/vitest run test/StaleUrlPage.test.tsx
```

Expected: module not found.

- [ ] **Step 3: Implement the component**

Create `frontend/src/components/StaleUrlPage.tsx`:

```tsx
import { useEffect } from "react";
import { useAppState } from "../state";
import type { StaleUrlContext } from "../state";

interface Props {
  staleUrlContext: StaleUrlContext;
}

interface VariantUI {
  dataMissing: string;
  header: string;
  ctaLabel: string;
  onPick: () => void;
}

function uiFor(ctx: StaleUrlContext, store: ReturnType<typeof useAppState.getState>): VariantUI {
  if (ctx.kind === "unknown_path") {
    return {
      dataMissing: "path",
      header: "Page not found.",
      ctaLabel: "Go to Index",
      onPick: () => store.setActivePage("index"),
    };
  }
  if (ctx.missing === "experiment") {
    return {
      dataMissing: "experiment",
      header: `Experiment '${ctx.missing_value}' not found.`,
      ctaLabel: "Select an experiment",
      onPick: () => store.recoverFromStaleUrl({
        step: "experiment",
        experimentId: undefined,
        sampleId: undefined,
      }),
    };
  }
  if (ctx.missing === "sample") {
    const expName = ctx.experiment_resolved?.name ?? "?";
    const expId = ctx.experiment_resolved?.id;
    return {
      dataMissing: "sample",
      header: `Sample '${ctx.missing_value}' not found in '${expName}'.`,
      ctaLabel: "Select another sample",
      onPick: () => store.recoverFromStaleUrl({
        step: "sample",
        experimentId: expId,
        sampleId: undefined,
      }),
    };
  }
  // missing === "exposure"
  const sampleName = ctx.sample_resolved?.name ?? "?";
  const expId = ctx.experiment_resolved?.id;
  const sampleId = ctx.sample_resolved?.id;
  return {
    dataMissing: "exposure",
    header: `Exposure '${ctx.missing_value}' not found in '${sampleName}'.`,
    ctaLabel: "Back to sample",
    onPick: () => store.recoverFromStaleUrl({
      step: "sample",
      experimentId: expId,
      sampleId: sampleId,
      openModal: false,
    }),
  };
}

export function StaleUrlPage({ staleUrlContext }: Props): JSX.Element {
  const store = useAppState.getState();
  const ui = uiFor(staleUrlContext, store);

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "/" && !e.metaKey && !e.ctrlKey && !e.altKey) {
        e.preventDefault();
        ui.onPick();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [ui]);

  return (
    <div
      role="alert"
      data-testid="stale-url-page"
      data-missing={ui.dataMissing}
      className="flex-1 min-h-0 flex flex-col items-center justify-center p-8 text-center"
    >
      <h2 className="text-xl mb-2 text-fg">{ui.header}</h2>
      <p className="text-fg-muted mb-6">It may have been renamed or removed.</p>
      <button
        onClick={ui.onPick}
        data-testid="stale-url-cta"
        className="px-4 py-2 rounded border border-accent hover:bg-accent/10"
      >
        {ui.ctaLabel}
        <kbd className="ml-2 text-xs opacity-60">/</kbd>
      </button>
    </div>
  );
}
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
node_modules/.bin/vitest run test/StaleUrlPage.test.tsx
```

Expected: all assertions pass.

- [ ] **Step 5: Commit**

```bash
git add frontend/src/components/StaleUrlPage.tsx frontend/test/StaleUrlPage.test.tsx
git commit -m "feat(component): StaleUrlPage 404 with per-variant CTA"
```

---

## Task 11: AppShell — add `<Route>` declarations for index/inspect + ladder

**Files:**
- Modify: `frontend/src/components/AppShell.tsx`

The existing AppShell uses `<Routes>` for Compare URLs and falls through to a `<ZustandShellPage />` element for everything else. We need to:
1. Add new `<Route>` declarations for `/`, `/index`, `/index/<exp>`, `/index/<exp>/<sample>`, `/inspect`, `/inspect/<exp>`, `/inspect/<exp>/<sample>` so react-router subscribes correctly.
2. Wrap the page-router output in the resolving/stale ladder (only the page region swaps; chrome stays unconditional).
3. Replace `ZustandShellPage` with the new `<PageBody />` ladder.

- [ ] **Step 1: Read current AppShell**

```bash
sed -n '1,140p' frontend/src/components/AppShell.tsx
```

Note: chrome (`<AppHeader />`, the rocker `<div>` with `<TabRocker />`, and `<NavModal />`) is sibling to `<Routes>`. The new ladder lives INSIDE the `<Route path="*" />` element, replacing `ZustandShellPage`.

- [ ] **Step 2: Replace `ZustandShellPage` with `PageBody`**

Find:

```tsx
function ZustandShellPage(): JSX.Element {
  const activePage = useAppState((s) => s.activePage);
  return (
    <>
      {activePage === "index"   && <IndexPage />}
      {activePage === "inspect" && <InspectPage />}
    </>
  );
}
```

Replace with:

```tsx
import { ResolvingFallback } from "./ResolvingFallback";
import { StaleUrlPage } from "./StaleUrlPage";

function PageBody(): JSX.Element {
  const activePage = useAppState((s) => s.activePage);
  const resolving = useAppState((s) => s.resolving);
  const staleUrlContext = useAppState((s) => s.staleUrlContext);

  if (resolving) return <ResolvingFallback />;
  if (staleUrlContext !== null) {
    return <StaleUrlPage staleUrlContext={staleUrlContext} />;
  }
  if (activePage === "index")   return <IndexPage />;
  if (activePage === "inspect") return <InspectPage />;
  // activePage === "compare" never reaches here because compare URLs are
  // matched by their explicit <Route> entries above.
  return <></>;
}
```

- [ ] **Step 3: Add `<Route>` declarations for index/inspect URLs**

Find the existing `<Routes>` block (`AppShell.tsx:112–130`):

```tsx
<Routes>
  <Route path="/experiments/:eid/compare" element={<ComparePage />} />
  ...
  <Route path="*" element={<ZustandShellPage />} />
</Routes>
```

Replace the `path="*"` line with the new declarations and a refined fallback:

```tsx
<Routes>
  <Route path="/experiments/:eid/compare" element={<ComparePage />} />
  <Route path="/experiments/:eid/compare/new" element={<ComparePageEdit />} />
  <Route path="/experiments/:eid/compare/:id" element={<ComparePage />} />
  <Route path="/experiments/:eid/compare/:id/edit" element={<ComparePageEdit />} />
  <Route path="/compare/all" element={<ComparePage />} />
  <Route path="/compare/all/new" element={<ComparePageEdit />} />
  <Route path="/compare/all/:id" element={<ComparePage />} />
  <Route path="/compare/all/:id/edit" element={<ComparePageEdit />} />

  {/* New permalink shapes — all render PageBody, which inspects Zustand
      to decide which page to mount. The URL-sync hooks dispatch state
      based on the matched route, so PageBody only needs to read the
      already-populated Zustand. */}
  <Route path="/" element={<PageBody />} />
  <Route path="/index" element={<PageBody />} />
  <Route path="/index/:experiment" element={<PageBody />} />
  <Route path="/index/:experiment/:sample" element={<PageBody />} />
  <Route path="/inspect" element={<PageBody />} />
  <Route path="/inspect/:experiment" element={<PageBody />} />
  <Route path="/inspect/:experiment/:sample" element={<PageBody />} />

  <Route path="*" element={<PageBody />} />  {/* stale fallback */}
</Routes>
```

- [ ] **Step 4: Verify build**

```bash
cd packages/HimalayaUI/frontend
node_modules/.bin/tsc --noEmit && node_modules/.bin/vitest run
```

Expected: no new TS errors, existing tests still pass.

- [ ] **Step 5: Commit**

```bash
git add frontend/src/components/AppShell.tsx
git commit -m "feat(shell): add index/inspect <Route> + PageBody ladder

PageBody renders <ResolvingFallback>, <StaleUrlPage>, or the page-router
output based on Zustand resolving + staleUrlContext slots. New <Route>
declarations make react-router subscribe to index/inspect URL changes."
```

---

## Task 12: `useStateFromUrl` hook (react-router integrated)

**Files:**
- Create: `frontend/src/hooks/useStateFromUrl.ts`
- Create: `frontend/test/useStateFromUrl.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `frontend/test/useStateFromUrl.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { useStateFromUrl } from "../src/hooks/useStateFromUrl";
import { useAppState } from "../src/state";
import type { ResolveSuccess, ResolveError404 } from "../src/api";

// Spec §4.2

const ok = (body: ResolveSuccess) => Promise.resolve({
  ok: true, status: 200, json: () => Promise.resolve(body),
} as Response);

const notFound = (body: ResolveError404) => Promise.resolve({
  ok: false, status: 404, json: () => Promise.resolve(body),
} as Response);

beforeEach(() => {
  useAppState.setState({
    staleUrlContext: null, resolving: false,
    activeExperimentId: undefined, activeSampleId: undefined, activeExposureId: undefined,
    activePage: "index",
  });
  history.replaceState(null, "", "/");
});

afterEach(() => {
  vi.restoreAllMocks();
});

describe("useStateFromUrl", () => {
  it("on mount with /index/<exp>/<sample>: dispatches setActive*, clears stale", async () => {
    history.replaceState(null, "", "/index/lipid/JC001");
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(() =>
      ok({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    );
    renderHook(() => useStateFromUrl());
    await waitFor(() => {
      expect(useAppState.getState().activeExperimentId).toBe(17);
    });
    expect(useAppState.getState().activeSampleId).toBe(42);
    expect(useAppState.getState().activePage).toBe("index");
    expect(useAppState.getState().resolving).toBe(false);
    expect(fetchSpy).toHaveBeenCalledOnce();
  });

  it("on mount with /<unknown>/path: sets staleUrlContext kind=unknown_path, no fetch", () => {
    history.replaceState(null, "", "/foo/bar");
    const fetchSpy = vi.spyOn(global, "fetch");
    renderHook(() => useStateFromUrl());
    const ctx = useAppState.getState().staleUrlContext;
    expect(ctx?.kind).toBe("unknown_path");
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it("404: sets staleUrlContext from response body", async () => {
    history.replaceState(null, "", "/index/lipid-typo");
    vi.spyOn(global, "fetch").mockImplementation(() =>
      notFound({
        error: "not_found", missing: "experiment", missing_value: "lipid-typo",
        experiment_resolved: undefined, sample_resolved: undefined,
      }),
    );
    renderHook(() => useStateFromUrl());
    await waitFor(() => {
      const ctx = useAppState.getState().staleUrlContext;
      expect(ctx?.kind).toBe("not_found");
    });
  });

  it("origin-tag race: navigation mid-resolve discards the response", async () => {
    history.replaceState(null, "", "/index/lipid/JC001");
    let resolveFetch: ((v: Response) => void) | undefined;
    vi.spyOn(global, "fetch").mockImplementation(() =>
      new Promise<Response>((res) => { resolveFetch = res; }),
    );
    renderHook(() => useStateFromUrl());
    expect(useAppState.getState().resolving).toBe(true);

    // Simulate user navigating mid-flight via TabRocker (Zustand mutation, no popstate).
    history.replaceState(null, "", "/inspect/other/SAMPLE");
    // Now satisfy the original fetch — its origin tag should mismatch.
    resolveFetch?.({
      ok: true, status: 200,
      json: () => Promise.resolve({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    } as Response);
    await new Promise((r) => setTimeout(r, 0));
    // Active state should NOT have been updated by the discarded resolve.
    expect(useAppState.getState().activeExperimentId).toBeUndefined();
  });

  it("pre-fetch clear: staleUrlContext is cleared at start of recognized fetch", async () => {
    useAppState.setState({
      staleUrlContext: { kind: "unknown_path", raw: "/old" },
    });
    history.replaceState(null, "", "/index/lipid/JC001");
    vi.spyOn(global, "fetch").mockImplementation(() =>
      ok({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    );
    renderHook(() => useStateFromUrl());
    // The clear happens synchronously before fetch; assert without waiting.
    expect(useAppState.getState().staleUrlContext).toBeNull();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

```bash
node_modules/.bin/vitest run test/useStateFromUrl.test.tsx
```

Expected: module not found.

- [ ] **Step 3: Implement the hook**

Create `frontend/src/hooks/useStateFromUrl.ts`:

```ts
import { useEffect } from "react";
import { useLocation, useNavigate } from "react-router-dom";
import { useQueryClient } from "@tanstack/react-query";
import { useAppState } from "../state";
import { parseLocation, type ParsedUrl } from "../lib/url/parseLocation";
import * as api from "../api";
import type { ResolveSuccess, ResolveError404, Experiment, Sample } from "../api";
import { queryKeys } from "../queries";
import { emitReplaceNext } from "../lib/url/emitMode";

// Spec §4.2 — URL → Zustand. Reads `useLocation()` so popstate AND
// useNavigate both flow through the hook. Origin-tagged fetches avoid
// the Zustand-mid-flight race (Zustand mutations don't change `location`,
// so AbortController alone is insufficient).

/**
 * Atomic apply of a 200 resolve response. Single setState commit so
 * useUrlFromState recomputes once, no cascading partial URL emits.
 *
 * App-init URL sync (per spec §4.3 trigger table) → replace, so the
 * initial state→URL emit doesn't push a redundant history entry over
 * the URL the user just landed on.
 */
function applySuccess(body: ResolveSuccess, page: "index" | "inspect" | "compare") {
  emitReplaceNext();
  useAppState.setState({
    activePage: page,
    activeExperimentId: body.experiment_id,
    activeSampleId: body.sample_id,
    activeExposureId: body.exposure_id,
    staleUrlContext: null,
    resolving: false,
  });
}

const VALID_PAGES = new Set(["index", "inspect", "compare"]);

export function useStateFromUrl(): void {
  const location = useLocation();
  const navigate = useNavigate();
  const qc = useQueryClient();

  useEffect(() => {
    let cancelled = false;
    const origin = location.pathname + location.search;
    const parsed = parseLocation(location.pathname, location.search);

    if (parsed.kind === "stale") {
      useAppState.setState({
        staleUrlContext: { kind: "unknown_path", raw: parsed.raw },
        resolving: false,
      });
      return;
    }

    if (parsed.kind === "compare") {
      // Compare uses numeric ids resolved by react-router useParams in the
      // ComparePage component itself; just set the active page.
      useAppState.setState({
        activePage: "compare",
        staleUrlContext: null,
        resolving: false,
      });
      return;
    }

    if (parsed.kind === "root") {
      // §5 redirect: build slug URL from persisted Zustand and navigate
      // (replace). All emits go through useNavigate so react-router's
      // location subscription stays consistent.
      const s = useAppState.getState();
      const page = VALID_PAGES.has(s.activePage) ? s.activePage : "index";
      const expId = s.activeExperimentId;
      const sId = s.activeSampleId;
      if (page === "compare") {
        navigate(expId !== undefined
          ? `/experiments/${expId}/compare`
          : "/compare/all", { replace: true });
        return;
      }
      if (expId === undefined) {
        navigate(`/${page}`, { replace: true });
        return;
      }
      // Synchronous cache hit?
      const exps = qc.getQueryData<Experiment[]>(queryKeys.experiments) ?? [];
      const expName = exps.find((e) => e.id === expId)?.name;
      const samples = sId !== undefined
        ? qc.getQueryData<Sample[]>(queryKeys.samples(expId)) ?? []
        : [];
      const sName = sId !== undefined ? samples.find((x) => x.id === sId)?.name : null;
      if (expName !== null && expName !== undefined &&
          (sId === undefined || (sName !== null && sName !== undefined))) {
        const path = sName !== undefined && sName !== null
          ? `/${page}/${encodeURIComponent(expName)}/${encodeURIComponent(sName)}`
          : `/${page}/${encodeURIComponent(expName)}`;
        navigate(path, { replace: true });
        return;
      }
      // Cold-mount fallback: resolve-by-id.
      (async () => {
        const q: api.ResolveQuery = { experiment_id: expId };
        if (sId !== undefined) q.sample_id = sId;
        let body;
        try {
          body = await api.resolve(q);
        } catch {
          if (!cancelled) navigate("/index", { replace: true });
          return;
        }
        if (cancelled) return;
        if ("error" in body) {
          navigate("/index", { replace: true });
          return;
        }
        const path = body.sample_name !== undefined && body.sample_name !== null
          ? `/${page}/${encodeURIComponent(body.experiment_name)}/${encodeURIComponent(body.sample_name)}`
          : `/${page}/${encodeURIComponent(body.experiment_name)}`;
        navigate(path, { replace: true });
      })();
      return;
    }

    // index or inspect — fetch /api/resolve with whichever slugs are present.
    if (parsed.experiment === undefined) {
      useAppState.setState({
        activePage: parsed.kind,
        activeExperimentId: undefined,
        activeSampleId: undefined,
        activeExposureId: undefined,
        staleUrlContext: null,
        resolving: false,
      });
      return;
    }

    // Pre-fetch clear of staleUrlContext + set resolving.
    useAppState.setState({ staleUrlContext: null, resolving: true });

    const ctl = new AbortController();
    const q: api.ResolveQuery = { experiment: parsed.experiment };
    if (parsed.sample !== undefined) q.sample = parsed.sample;
    if (parsed.kind === "inspect" && parsed.exposure !== undefined) {
      q.exposure = parsed.exposure;
    }

    (async () => {
      let body;
      try {
        body = await api.resolve(q, ctl.signal);
      } catch (e) {
        if ((e as Error).name === "AbortError") return;
        if (!cancelled) useAppState.setState({ resolving: false });
        return;
      }
      // Origin-tag check: did the URL change during the fetch?
      if (cancelled || (window.location.pathname + window.location.search) !== origin) {
        return;
      }
      if ("error" in body && body.error === "not_found") {
        useAppState.setState({
          staleUrlContext: {
            kind: "not_found",
            missing: body.missing,
            missing_value: body.missing_value,
            experiment_resolved: body.experiment_resolved,
            sample_resolved: body.sample_resolved,
          },
          resolving: false,
        });
        return;
      }
      if ("error" in body) {
        useAppState.setState({
          staleUrlContext: { kind: "unknown_path", raw: origin },
          resolving: false,
        });
        return;
      }
      applySuccess(body as ResolveSuccess, parsed.kind);
    })();

    return () => {
      cancelled = true;
      ctl.abort();
    };
  }, [location.pathname, location.search, navigate, qc]);
}
```

Key changes vs. the earlier draft:
- **Reads `useLocation()` from react-router** — no manual popstate listener.
- **Effect deps `[location.pathname, location.search, navigate, qc]`** — fires on every URL change (popstate or `useNavigate`).
- **`applySuccess` is one atomic `setState`** — no cascading `setActive*` calls — and calls `emitReplaceNext()` so the URL→state→URL roundtrip on cold mount uses replace mode (matches spec §4.3).
- **Origin-tag race** uses `window.location.pathname + window.location.search` at response time.
- **Root redirect** uses `useNavigate(target, { replace: true })` everywhere — no raw `history.replaceState`, no synthetic `PopStateEvent`.
- **Page validation** — `VALID_PAGES` set guards against an unrecognized persisted `activePage` falling through to `/<garbage>`.

- [ ] **Step 4: Run test to verify it passes**

```bash
node_modules/.bin/vitest run test/useStateFromUrl.test.tsx
```

Expected: all assertions pass.

- [ ] **Step 5: Commit**

```bash
git add frontend/src/hooks/useStateFromUrl.ts frontend/test/useStateFromUrl.test.tsx
git commit -m "feat(hook): useStateFromUrl — URL → Zustand with origin-tag race

Pre-fetch clears staleUrlContext on recognized-kind URLs. Origin tag
prevents Zustand-mutation-mid-resolve from clobbering the user's
navigation choice."
```

---

## Task 13: `useUrlFromState` hook

**Files:**
- Create: `frontend/src/hooks/useUrlFromState.ts`
- Create: `frontend/test/useUrlFromState.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `frontend/test/useUrlFromState.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useUrlFromState } from "../src/hooks/useUrlFromState";
import { useAppState } from "../src/state";
import { queryKeys } from "../src/queries";
import type { Experiment, Sample } from "../src/api";

// Spec §4.3

function makeWrapper() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  qc.setQueryData<Experiment[]>(queryKeys.experiments, [
    { id: 17, name: "lipid-screen", path: "", data_dir: "", analysis_dir: "",
      manifest_path: null, created_at: "", q_units: null },
  ]);
  qc.setQueryData<Sample[]>(queryKeys.samples(17), [
    { id: 42, experiment_id: 17, name: "JC001", display_name: null, notes: null, tags: [] },
  ]);
  const Wrapper = ({ children }: { children: React.ReactNode }) =>
    <QueryClientProvider client={qc}>{children}</QueryClientProvider>;
  return { qc, Wrapper };
}

beforeEach(() => {
  history.replaceState(null, "", "/");
  useAppState.setState({
    activePage: "index",
    activeExperimentId: undefined, activeSampleId: undefined, activeExposureId: undefined,
    staleUrlContext: null, resolving: false,
  });
});

describe("useUrlFromState", () => {
  it("active sample → /index/<exp>/<sample>", () => {
    const { Wrapper } = makeWrapper();
    useAppState.setState({ activePage: "index", activeExperimentId: 17, activeSampleId: 42 });
    renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    expect(location.pathname).toBe("/index/lipid-screen/JC001");
  });

  it("page change → push", () => {
    const { Wrapper } = makeWrapper();
    history.replaceState(null, "", "/index/lipid-screen/JC001");
    useAppState.setState({ activePage: "index", activeExperimentId: 17, activeSampleId: 42 });
    const pushSpy = vi.spyOn(history, "pushState");
    const replaceSpy = vi.spyOn(history, "replaceState");
    const { rerender } = renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    pushSpy.mockClear(); replaceSpy.mockClear();
    useAppState.setState({ activePage: "inspect" });
    rerender();
    expect(pushSpy).toHaveBeenCalledWith(null, "", "/inspect/lipid-screen/JC001");
  });

  it("exposure-only change → replace", () => {
    const { qc, Wrapper } = makeWrapper();
    qc.setQueryData(queryKeys.exposures(42), [
      { id: 100, sample_id: 42, filename: "JC001-007", kind: "simple", selected: 1 },
    ]);
    history.replaceState(null, "", "/inspect/lipid-screen/JC001");
    useAppState.setState({ activePage: "inspect", activeExperimentId: 17, activeSampleId: 42 });
    const pushSpy = vi.spyOn(history, "pushState");
    const replaceSpy = vi.spyOn(history, "replaceState");
    const { rerender } = renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    pushSpy.mockClear(); replaceSpy.mockClear();
    useAppState.setState({ activeExposureId: 100 });
    rerender();
    expect(replaceSpy).toHaveBeenCalled();
    expect(pushSpy).not.toHaveBeenCalled();
    expect(location.search).toBe("?exposure=JC001-007");
  });

  it("equality guard: identical URL does not emit", () => {
    const { Wrapper } = makeWrapper();
    history.replaceState(null, "", "/index/lipid-screen/JC001");
    useAppState.setState({ activePage: "index", activeExperimentId: 17, activeSampleId: 42 });
    const pushSpy = vi.spyOn(history, "pushState");
    const replaceSpy = vi.spyOn(history, "replaceState");
    const { rerender } = renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    pushSpy.mockClear(); replaceSpy.mockClear();
    rerender();
    expect(pushSpy).not.toHaveBeenCalled();
    expect(replaceSpy).not.toHaveBeenCalled();
  });

  it("replay-as-rerun: identical optimistic + confirmed slug → no spurious emit", () => {
    // Simulate the trivial replay case: cache row gets replaced (foreign event)
    // but the same id-name mapping holds. URL recompute should see the same
    // slug and not emit.
    const { qc, Wrapper } = makeWrapper();
    history.replaceState(null, "", "/index/lipid-screen/JC001");
    useAppState.setState({ activePage: "index", activeExperimentId: 17, activeSampleId: 42 });
    const replaceSpy = vi.spyOn(history, "replaceState");
    const { rerender } = renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    replaceSpy.mockClear();
    // Simulate applyRemoteToCache rewriting samples in place.
    qc.setQueryData<Sample[]>(queryKeys.samples(17), [
      { id: 42, experiment_id: 17, name: "JC001", display_name: "JC001 (touched)",
        notes: null, tags: [] },
    ]);
    rerender();
    expect(replaceSpy).not.toHaveBeenCalled();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

```bash
node_modules/.bin/vitest run test/useUrlFromState.test.tsx
```

Expected: module not found.

- [ ] **Step 3: Implement the emit-mode helper (hoisted)**

Create `frontend/src/lib/url/emitMode.ts`:

```ts
// Spec §4.3 — push/replace mode for the next useUrlFromState emit.
// Hoisted into its own module to avoid a state.ts ↔ useUrlFromState cycle
// and to keep test isolation simpler (consume resets to push on read).

let nextEmitMode: "push" | "replace" = "push";

export function emitReplaceNext(): void {
  nextEmitMode = "replace";
}

export function consumeEmitMode(): "push" | "replace" {
  const mode = nextEmitMode;
  nextEmitMode = "push";
  return mode;
}

// For tests: reset between cases.
export function _resetEmitMode(): void {
  nextEmitMode = "push";
}
```

- [ ] **Step 4: Implement the hook**

Create `frontend/src/hooks/useUrlFromState.ts`:

```ts
import { useEffect } from "react";
import { useNavigate, useLocation } from "react-router-dom";
import { useQueryClient, useQuery } from "@tanstack/react-query";
import { useAppState } from "../state";
import { useExperiments, queryKeys } from "../queries";
import * as api from "../api";
import type { Experiment, Sample, Exposure } from "../api";
import { consumeEmitMode } from "../lib/url/emitMode";

// Spec §4.3 — Zustand → URL via react-router useNavigate. Subscribes to
// experiments + samples queries via TanStack so SSE-driven cache rewrites
// trigger a re-render of this hook (§7 invalidation).

function nameForExperiment(list: Experiment[] | undefined, id: number | undefined):
  string | undefined
{
  if (id === undefined || list === undefined) return undefined;
  const found = list.find((e) => e.id === id);
  if (found === undefined) return undefined;
  return found.name === null ? undefined : found.name;
}

function nameForSample(list: Sample[] | undefined, sId: number | undefined):
  string | undefined
{
  if (sId === undefined || list === undefined) return undefined;
  const found = list.find((s) => s.id === sId);
  if (found === undefined) return undefined;
  return found.name === null ? undefined : found.name;
}

function filenameForExposure(qc: ReturnType<typeof useQueryClient>,
                              sId: number | undefined,
                              eId: number | undefined): string | undefined
{
  if (sId === undefined || eId === undefined) return undefined;
  const list = qc.getQueryData<Exposure[]>(queryKeys.exposures(sId)) ?? [];
  const found = list.find((e) => e.id === eId);
  if (found === undefined) return undefined;
  return found.filename === null ? undefined : found.filename;
}

function buildUrl(
  page: "index" | "inspect" | "compare",
  experiment: string | undefined,
  sample: string | undefined,
  exposure: string | undefined,
  current: string,
): string {
  const enc = (s: string) => encodeURIComponent(s);
  if (page === "compare") {
    // Don't try to re-emit a Compare URL — Compare uses numeric ids that
    // useStateFromUrl doesn't track in Zustand. ComparePage handles its
    // own URL via useNavigate. Returning current keeps the equality guard
    // happy and prevents accidental redirect.
    return current;
  }
  const parts = [`/${page}`];
  if (experiment !== undefined) {
    parts.push(enc(experiment));
    if (sample !== undefined) parts.push(enc(sample));
  }
  let url = parts.join("/");
  if (page === "inspect" && exposure !== undefined && sample !== undefined) {
    url += `?exposure=${enc(exposure)}`;
  }
  return url;
}

export function useUrlFromState(): void {
  const navigate = useNavigate();
  const location = useLocation();
  const qc = useQueryClient();

  const activePage = useAppState((s) => s.activePage);
  const activeExperimentId = useAppState((s) => s.activeExperimentId);
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const activeExposureId = useAppState((s) => s.activeExposureId);

  // Subscribe to experiments + (when an experiment is active) samples
  // queries so SSE-driven cache rewrites trigger a re-render. We use
  // `useQuery` directly here (rather than `useSamples`) so we can gate
  // on `enabled: activeExperimentId !== undefined` — `useSamples` doesn't
  // expose `enabled`, and calling it with `0` would fire `GET
  // /api/experiments/0/samples` → 404 with retries on every cold mount.
  const { data: experiments } = useExperiments();
  const samplesQuery = useQuery({
    queryKey: queryKeys.samples(activeExperimentId ?? 0),
    queryFn: () => api.listSamples(activeExperimentId as number),
    enabled: activeExperimentId !== undefined,
  });
  const samples = activeExperimentId !== undefined ? samplesQuery.data : undefined;

  // Track the previous resolved slug pair so we can detect SSE-driven
  // disappearance (a slug went from defined → undefined because the
  // entity was deleted from the cache). Per spec §4.3 + §7, that case
  // should emit replace, not push (otherwise back-button stops at the
  // broken URL).
  const prevSlugsRef = useRef<{ exp: string | undefined; sample: string | undefined }>(
    { exp: undefined, sample: undefined },
  );

  useEffect(() => {
    const expName = nameForExperiment(experiments, activeExperimentId);
    const sampleName = nameForSample(samples, activeSampleId);
    const exposureName = filenameForExposure(qc, activeSampleId, activeExposureId);
    const current = location.pathname + location.search;
    const target = buildUrl(activePage, expName, sampleName, exposureName, current);

    // SSE-driven invalidation detection: a previously-resolvable slug
    // is now undefined (the row vanished from cache). Force replace.
    const prev = prevSlugsRef.current;
    const slugDisappeared =
      (prev.exp !== undefined && expName === undefined && activeExperimentId !== undefined) ||
      (prev.sample !== undefined && sampleName === undefined && activeSampleId !== undefined);
    prevSlugsRef.current = { exp: expName, sample: sampleName };

    if (target === current) return;        // equality guard
    const explicitMode = consumeEmitMode();
    const mode = slugDisappeared ? "replace" : explicitMode;
    navigate(target, { replace: mode === "replace" });
  }, [
    activePage, activeExperimentId, activeSampleId, activeExposureId,
    experiments, samples,
    location.pathname, location.search,
    navigate, qc,
  ]);
}
```

Add `useRef` to the imports:

```ts
import { useEffect, useRef } from "react";
```

- [ ] **Step 5: Wire `emitReplaceNext()` into state.ts call sites**

In `frontend/src/state.ts`, add the import at the top:

```ts
import { emitReplaceNext } from "./lib/url/emitMode";
```

For `setActiveExposure`:
```ts
        setActiveExposure: (activeExposureId) => {
          emitReplaceNext();
          set({ activeExposureId, staleUrlContext: null });
        },
```

For `recoverFromStaleUrl`:
```ts
        recoverFromStaleUrl: (opts) => {
          emitReplaceNext();
          set((s) => ({
            staleUrlContext: null,
            activeExperimentId: opts.experimentId ?? s.activeExperimentId,
            activeSampleId: opts.sampleId ?? undefined,
            activeExposureId: undefined,
            navModalOpen: opts.openModal ?? true,
            navModalStep: opts.step,
          }));
        },
```

The `lib/url/emitMode.ts` module has no imports from `state.ts`, so there's no circular dependency.

- [ ] **Step 4: Run tests to verify they pass**

```bash
node_modules/.bin/vitest run test/useUrlFromState.test.tsx test/staleUrlState.test.ts
```

Expected: all assertions pass for both files. (Re-running staleUrlState.test ensures the `emitReplaceNext` import didn't break it.)

- [ ] **Step 5: If the cycle bites, hoist to `lib/url/emitMode.ts`**

If you see a circular-dep warning or a runtime "ReferenceError: cannot access … before initialization", create `frontend/src/lib/url/emitMode.ts`:

```ts
let nextEmitMode: "push" | "replace" = "push";

export function emitReplaceNext(): void {
  nextEmitMode = "replace";
}

export function consumeEmitMode(): "push" | "replace" {
  const mode = nextEmitMode;
  nextEmitMode = "push";
  return mode;
}
```

Update `useUrlFromState.ts` to import `consumeEmitMode` and call it at the emit site instead of reading `nextEmitMode` directly. Update `state.ts` to import from `./lib/url/emitMode`. Re-run tests.

- [ ] **Step 6: Commit**

```bash
git add frontend/src/hooks/useUrlFromState.ts frontend/src/state.ts frontend/test/useUrlFromState.test.tsx
git commit -m "feat(hook): useUrlFromState — Zustand → URL with push/replace policy

Equality guard suppresses spurious emits during replay-as-rerun. The
emit-mode flag is set by setActiveExposure and recoverFromStaleUrl
to force replace; default is push (TabRocker / NavModal commit)."
```

---

## Task 14: `/` cold-mount redirect — extra test coverage

The redirect logic is implemented inline in Task 12's `useStateFromUrl` (the `parsed.kind === "root"` branch). This task adds the test cases that prove all redirect paths work.

**Files:**
- Modify: `frontend/test/useStateFromUrl.test.tsx`

- [ ] **Step 1: Add the redirect test cases**

Append to `useStateFromUrl.test.tsx`:

```tsx
describe("useStateFromUrl — / redirect", () => {
  it("with persisted experiment+sample: synchronous getQueryData → replaceState to slug URL", async () => {
    history.replaceState(null, "", "/");
    useAppState.setState({
      activePage: "index", activeExperimentId: 17, activeSampleId: 42,
    });
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.experiments, [
      { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "",
        manifest_path: null, created_at: "", q_units: null },
    ]);
    qc.setQueryData(queryKeys.samples(17), [
      { id: 42, experiment_id: 17, name: "JC001", display_name: null, notes: null, tags: [] },
    ]);
    const Wrapper = ({ children }: { children: React.ReactNode }) =>
      <QueryClientProvider client={qc}>{children}</QueryClientProvider>;
    const replaceSpy = vi.spyOn(history, "replaceState");
    renderHook(() => useStateFromUrl(), { wrapper: Wrapper });
    // Expect a synchronous replaceState to /index/lipid/JC001 (no fetch).
    expect(replaceSpy).toHaveBeenCalledWith(null, "", "/index/lipid/JC001");
  });

  it("with persisted ids missing from cache: falls back to /api/resolve?experiment_id=…", async () => {
    history.replaceState(null, "", "/");
    useAppState.setState({
      activePage: "index", activeExperimentId: 17, activeSampleId: 42,
    });
    const qc = new QueryClient();   // empty cache
    const Wrapper = ({ children }: { children: React.ReactNode }) =>
      <QueryClientProvider client={qc}>{children}</QueryClientProvider>;
    vi.spyOn(global, "fetch").mockImplementation(() =>
      ok({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    );
    const replaceSpy = vi.spyOn(history, "replaceState");
    renderHook(() => useStateFromUrl(), { wrapper: Wrapper });
    await waitFor(() => expect(replaceSpy).toHaveBeenCalled());
    expect(replaceSpy).toHaveBeenCalledWith(null, "", "/index/lipid/JC001");
  });

  it("with no persisted state: replaceState to /index", () => {
    history.replaceState(null, "", "/");
    useAppState.setState({
      activePage: "index", activeExperimentId: undefined, activeSampleId: undefined,
    });
    const qc = new QueryClient();
    const Wrapper = ({ children }: { children: React.ReactNode }) =>
      <QueryClientProvider client={qc}>{children}</QueryClientProvider>;
    const replaceSpy = vi.spyOn(history, "replaceState");
    renderHook(() => useStateFromUrl(), { wrapper: Wrapper });
    expect(replaceSpy).toHaveBeenCalledWith(null, "", "/index");
  });
});
```

(`QueryClient` and `QueryClientProvider` already imported in the file from earlier; if not, add the import.)

- [ ] **Step 2: Run tests to verify they pass**

```bash
node_modules/.bin/vitest run test/useStateFromUrl.test.tsx
```

Expected: all assertions pass.

- [ ] **Step 3: Commit**

```bash
git add frontend/test/useStateFromUrl.test.tsx
git commit -m "test(hook): / cold-mount redirect coverage

getQueryData synchronous hit, resolve-by-id fallback, missing-state
fall-through to /index."
```

---

## Task 15: Wire hooks in `App.tsx`

**Files:**
- Modify: `frontend/src/App.tsx`

The hooks need to be called from a component that mounts INSIDE the `<BrowserRouter>` (which is in `main.tsx` wrapping `<App />`), so they can use `useLocation` and `useNavigate`. `App` itself is fine — it's already inside `<BrowserRouter>`.

- [ ] **Step 1: Add the hook calls**

In `frontend/src/App.tsx`, near the other top-of-component hooks:

```tsx
import { useStateFromUrl } from "./hooks/useStateFromUrl";
import { useUrlFromState } from "./hooks/useUrlFromState";
```

Inside the `App` component body, before the existing `useEffect`s:

```tsx
  useStateFromUrl();    // URL → Zustand
  useUrlFromState();    // Zustand → URL
```

Order matters on cold mount — `useStateFromUrl` populates Zustand first; equality guard makes order irrelevant after first render.

- [ ] **Step 2: Verify build + tests still pass**

```bash
node_modules/.bin/tsc --noEmit && node_modules/.bin/vitest run
```

Expected: no new failures.

- [ ] **Step 3: Smoke-test in dev**

```bash
# In one terminal, start the backend on port 8080:
cd /opt/Himalaya.jl/.claude/worktrees/permalinks
bin/himalaya serve /opt/Himalaya.jl/test-experiment --port 8080 &
# In another terminal, start Vite:
cd packages/HimalayaUI/frontend
npm run dev
# Open http://127.0.0.1:5173/index/<your-exp>/<your-sample> and confirm it lands on the right page.
```

Expected: a deep URL pasted into the address bar lands you at the right sample/page. The address bar updates as you click TabRocker, NavModal-commit, etc. Pasting `/foo/bar` shows the StaleUrlPage.

(Skip this step if a real test experiment isn't available; the Playwright tests in Task 18 cover end-to-end without manual setup.)

- [ ] **Step 4: Commit**

```bash
git add frontend/src/App.tsx
git commit -m "feat(app): mount useStateFromUrl + useUrlFromState"
```

---

## Task 16: Playwright mocked spec

**Files:**
- Create: `frontend/e2e/permalinks.spec.ts`

- [ ] **Step 1: Write the spec**

Create `frontend/e2e/permalinks.spec.ts`:

```ts
import { test, expect } from "@playwright/test";

// Spec §8.2 — Playwright mocked. /api/* is intercepted so this runs
// without a backend.

test.beforeEach(async ({ page }) => {
  // Common stubs: list experiments, list samples for the relevant experiment.
  await page.route("**/api/experiments", (route) => {
    route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify([
        { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "",
          manifest_path: null, created_at: "", q_units: null },
      ]),
    });
  });
  await page.route("**/api/experiments/17/samples", (route) => {
    route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify([
        { id: 42, experiment_id: 17, name: "JC001", display_name: null,
          notes: null, tags: [] },
      ]),
    });
  });
});

test("paste deep URL: lands on right page, no flash of wrong content", async ({ page }) => {
  await page.route("**/api/resolve?**", (route) => {
    route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    });
  });
  await page.goto("/index/lipid/JC001");
  // ResolvingFallback should appear briefly (or already past); page must not
  // show a different sample's content.
  await expect(page.locator("[data-testid='resolving']")).toBeAttached();
  // After resolve, page should be at the right sample.
  await expect(page).toHaveURL(/\/index\/lipid\/JC001$/);
});

test("paste stale URL: 404 page → CTA opens NavModal at right step", async ({ page }) => {
  await page.route("**/api/resolve?**", (route) => {
    route.fulfill({
      status: 404, contentType: "application/json",
      body: JSON.stringify({
        error: "not_found", missing: "experiment", missing_value: "lipid-typo",
        experiment_resolved: undefined, sample_resolved: undefined,
      }),
    });
  });
  await page.goto("/index/lipid-typo/JC001");
  await expect(page.locator("[data-testid='stale-url-page']")).toBeVisible();
  await expect(page.locator("[data-testid='stale-url-page']")).toHaveAttribute("data-missing", "experiment");
  await page.locator("[data-testid='stale-url-cta']").click();
  await expect(page.locator("[data-testid='nav-modal']")).toBeVisible();
});

test("TabRocker continuity: switch pages, URL is /<page>/<exp>/<sample>", async ({ page }) => {
  await page.route("**/api/resolve?**", (route) => {
    route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    });
  });
  await page.goto("/index/lipid/JC001");
  await page.locator("[data-testid='tab-inspect']").click();
  await expect(page).toHaveURL(/\/inspect\/lipid\/JC001$/);
  await page.goBack();
  await expect(page).toHaveURL(/\/index\/lipid\/JC001$/);
});

test("/ cold-mount: replaces to last-active slug URL", async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: {
        activePage: "index",
        activeExperimentId: 17,
        activeSampleId: 42,
        username: "test", firstName: undefined, lastName: undefined,
        activeExposureId: undefined,
        tutorialSeen: true, theme: "dark",
      },
      version: 3,
    }));
  });
  await page.goto("/");
  await expect(page).toHaveURL(/\/index\/lipid\/JC001$/);
});
```

- [ ] **Step 2: Run the spec**

```bash
cd packages/HimalayaUI/frontend
npm run e2e -- permalinks.spec.ts
```

Expected: all tests pass. `data-testid="tab-inspect"` is the canonical selector per `TabRocker.tsx:62` (`data-testid={\`tab-${t.id}\`}`); `data-testid="nav-modal"` per `NavModal.tsx:208`.

- [ ] **Step 3: Commit**

```bash
git add frontend/e2e/permalinks.spec.ts
git commit -m "test(e2e): mocked Playwright spec for permalink URLs"
```

---

## Task 17: Playwright live spec

**Files:**
- Create: `frontend/e2e/live/permalinks.spec.ts`

- [ ] **Step 1: Write the spec**

Create `frontend/e2e/live/permalinks.spec.ts`:

```ts
import { test, expect } from "@playwright/test";

// Spec §8.3 — Live integration. Requires backend (port 8090) + Vite (5180)
// per packages/HimalayaUI/frontend/e2e/live/README.md.

test.describe("permalinks — live", () => {
  test("delete sample in tab B → tab A URL invalidates to StaleUrlPage", async ({ browser }) => {
    const ctxA = await browser.newContext();
    const ctxB = await browser.newContext();
    const tabA = await ctxA.newPage();
    const tabB = await ctxB.newPage();

    // Tab A opens an inspect deep link for an existing sample.
    await tabA.goto("/inspect/E2ETest/SeedSample");
    await tabA.waitForTimeout(800);  // SSE handshake (per gotcha)

    // Tab B deletes the sample via direct API.
    await tabB.request.delete("/api/samples/SEED_SAMPLE_ID", {
      headers: { "X-Username": "alice" },
    });

    // Tab A's URL invalidates; StaleUrlPage renders.
    await expect(tabA.locator("[data-testid='stale-url-page']")).toBeVisible({ timeout: 5_000 });
  });

  test("paste URL referencing a deleted sample lands on StaleUrlPage", async ({ page }) => {
    await page.goto("/inspect/E2ETest/AlreadyDeletedSample");
    await expect(page.locator("[data-testid='stale-url-page']")).toBeVisible();
  });

  test("same-user-different-tab: tab A sees tab B's delete (no client_id self-echo)", async ({ browser }) => {
    // Spec §7 — URL invalidation runs on every SSE event regardless of
    // originating client_id, so two tabs of the same username refresh each
    // other's URL state.
    const ctx = await browser.newContext();
    const tabA = await ctx.newPage();
    const tabB = await ctx.newPage();
    // Same browser context = shared sessionStorage for clientId? No — clientId
    // is per-tab via sessionStorage, so the two tabs DO have distinct client_ids
    // even within a single context. The username header is the same on both.
    await tabA.goto("/inspect/E2ETest/SeedSample");
    await tabA.waitForTimeout(800);

    await tabB.request.delete("/api/samples/SEED_SAMPLE_ID", {
      headers: { "X-Username": "alice" },
    });

    await expect(tabA.locator("[data-testid='stale-url-page']")).toBeVisible({ timeout: 5_000 });
  });
});
```

(`SEED_SAMPLE_ID` and `SeedSample` / `E2ETest` reflect the seeded fixture; adapt to the live test runbook in `e2e/live/README.md`.)

- [ ] **Step 2: Run live (manual setup required)**

Per `packages/HimalayaUI/frontend/e2e/live/README.md`, bring up backend + Vite, then:

```bash
npm run e2e:live -- permalinks.spec.ts
```

Expected: tests pass against the seeded dev DB.

(Skip the run if no live environment is set up — the spec file is committed for the runbook to reference.)

- [ ] **Step 3: Commit**

```bash
git add frontend/e2e/live/permalinks.spec.ts
git commit -m "test(e2e-live): permalink live integration — SSE + stale URL"
```

---

## Task 18: Spec compliance + final test sweep

**Files:** None changed. Verification only.

- [ ] **Step 1: Re-read the spec**

```bash
cat docs/superpowers/specs/2026-05-09-permalink-urls-design.md
```

Walk each section:
- §2 grammar — covered by `parseLocation` (Task 7) + `useUrlFromState::buildUrl` (Task 13).
- §3.1 resolve route — Tasks 2, 3, 4.
- §3.2 SPA catch-all — Task 5.
- §4.1 ParsedUrl — Task 7.
- §4.2 useStateFromUrl + ResolvingFallback — Tasks 9, 12, 14.
- §4.3 useUrlFromState — Task 13.
- §4.4 mounting + named-action clearing — Tasks 6, 13, 15.
- §5 / redirect — Task 14.
- §6 StaleUrlPage + recoverFromStaleUrl + AppShell ladder — Tasks 6, 10, 11.
- §7 SSE invalidation — implicit via `useUrlFromState` reading the cache (no separate task; tested in Task 17 live).
- §8 tests — all covered.
- §9 file changes — see File Structure header above.
- §11 risks — no implementation impact.

If any section is missing a task, add one before moving on.

- [ ] **Step 2: Run the full test suites**

```bash
# Frontend
cd packages/HimalayaUI/frontend
npm test > /tmp/vitest.out 2>&1
grep -E "Test Files|Tests" /tmp/vitest.out
tail -20 /tmp/vitest.out

# Backend
cd /opt/Himalaya.jl/.claude/worktrees/permalinks
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-test.out
tail -30 /tmp/jl-test.out

# Frontend Playwright (mocked)
cd packages/HimalayaUI/frontend
npm run e2e > /tmp/pw.out 2>&1
grep -E "passed|failed" /tmp/pw.out
```

Expected: all suites green.

- [ ] **Step 3: Build the frontend bundle**

```bash
cd packages/HimalayaUI/frontend
npm run build
```

Expected: `tsc --noEmit` clean, `vite build` clean. The new TS strict types should compile without `any`.

- [ ] **Step 4: Run pre-merge smoke**

```bash
ls .claude/skills/pre-merge-smoke/
```

If a `pre-merge-smoke` skill exists, follow its checklist. Otherwise:
- Open the dev server and click through Index → Inspect → Compare → back; URL should track.
- Paste a deep URL; should land at the right page.
- Paste `/foo/bar`; should land on StaleUrlPage.

- [ ] **Step 5: Self-review against spec — final check**

For each spec section, confirm a task implemented it. Note any that were merely tested but not implemented; flag for fix.

- [ ] **Step 6: Open PR**

```bash
git push -u origin permalinks
gh pr create --title "feat: slug-based permalink URLs (#89)" --body "$(cat <<'EOF'
## Summary
- New `/api/resolve` route translates slugs → ids in one round trip
- SPA catch-all (`/**`) so deep URLs reach the frontend
- `useStateFromUrl` + `useUrlFromState` URL-sync hooks (no router lib)
- `<StaleUrlPage>` 404 with one-click recovery via `recoverFromStaleUrl`
- `migrate_exposures_unique_filename!` migration (dedupe-then-enforce)

Spec: `docs/superpowers/specs/2026-05-09-permalink-urls-design.md`
Closes #89.

## Test plan
- [ ] `npm test` (frontend Vitest)
- [ ] `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'`
- [ ] `npm run e2e` (mocked Playwright)
- [ ] `npm run e2e:live` (live, manual setup per runbook)
- [ ] Manual: paste deep URL, paste stale URL, TabRocker continuity, `/` cold-mount

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

Expected: PR created. The PR body's test-plan checklist drives the reviewer's smoke pass.

---

## Self-review notes (left in the plan for the implementer's reference)

- All 18 tasks have actual code blocks, not placeholders.
- Type names used in later tasks (`StaleUrlContext`, `RecoverOpts`, `ResolveSuccess`, `ParsedUrl`) are defined in earlier tasks.
- Action signatures consistent: `recoverFromStaleUrl(opts: RecoverOpts)` everywhere, `setStaleUrlContext(ctx | null)` everywhere.
- The `emitReplaceNext()` mechanism is described both inline (Task 13) and as a hoist option if a circular import bites — implementer chooses based on actual runtime behavior.
- Backend tasks 1–5 are independent of the frontend tasks; can land as a prep PR if desired.
- `_setup_analyzed_exposure` and `with_test_server` are existing test helpers (see `test_picker_routes.jl` for the pattern); the new tests use them as-is.

---
