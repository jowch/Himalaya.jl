# App Shell Unification — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Collapse the three coexisting app chromes into one unified two-tier shell (lean TopNav + contextual dock), wire one holistic keyboard set, and complete the two-phase ingest funnel — without re-laying-out the refined surfaces.

**Architecture:** Four independently-shippable milestones in dependency order. **M1 (backend)** adds the `/api/fs/*` endpoints + an `on_progress` ingest callback + pattern parsing — pure Julia, no frontend dep. **M2 (keyboard + sheet)** re-axes the shortcut registry, removes the roving grid in favour of a page-level `{sampleIndex, frameIndex}` cursor, and adds the per-row Restore button. **M3 (shell + dock)** collapses the chromes into one TopNav + a contextual bottom dock. **M4 (funnel)** builds the phase-1 manifest UI on Configuration, the ScanFailedPage, and wires the per-experiment corpus sheet. Each milestone is green-on-its-own.

**Tech Stack:** Julia / Oxygen.jl / SQLite (backend); React 18 / TypeScript (strict, `exactOptionalPropertyTypes`) / Vite / Zustand / TanStack Query / Tailwind v4 `@theme` (frontend). Tests: stdlib `Test` + in-process route dispatch (Julia); Vitest + React Testing Library (frontend); Playwright (e2e).

## Global Constraints

Every task's requirements implicitly include these — copied verbatim from the spec and the project's standing rules:

- **The Print design laws:** appearance lives in `src/print/ui/**` primitives; consumers pass **placement-only** `className`. **NO** inline appearance utilities (`text-[…]`, `rounded-[…]`, raw colour literals, **side-stripe borders**, gradient text) outside the appearance-exempt dirs. Enforced by `npm run lint:design` (a build step) — the build fails otherwise.
- **One radius step** (`rounded.sm` == `rounded.md`, 5px). Reuse `--color-print-accent` (sources `--color-accent`).
- **Never `git add -A`** — stage only the exact files each task names.
- **Branch stays unmerged** — this is `ingestion-redesign`; do not merge to main. Commit per-task; do not push unless asked.
- **TDD by default** — failing test → minimal impl → verify pass → commit. Each task ends green.
- **Backend route tests dispatch in-process** via `with_inproc_routes` (NOT over a socket) — see `packages/HimalayaUI/test/AGENTS.md`. Capture the suite once (`make test-parallel`, ~2 min) and grep; don't re-run per filter.
- **Commit trailers** (every commit):
  ```
  Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01VCzR9pYAh8XgkCf9tCcVDn
  ```
- **`scan_and_group!` keeps its single structural transaction** — `on_progress` fires from the per-exposure analysis loop *outside* the txn; do NOT split the structural write into per-load commits (spec §6.2).
- **Cull verbs are per-surface:** Corpus = `{Drop, Keep, Restore}` (NO `R`); Loupe = `{Drop, Keep, Set-representative, Restore}`.
- **Keyboard axes (app-wide):** `↑/↓` = sample · `←/→` = detail (frame on Corpus/Loupe, candidate on Focus) · `Enter` = Focus · `L` = Loupe · `X/K` = Drop/Keep · `R` = Set-rep (Loupe only) · `Space` = toggle-select · `Shift+←/→` = range · `Backspace` = restore · `Esc` = ladder · `Alt+↑/↓` = reorder · `/`,`⌘K` = find · `?` = help. This **supersedes the 2026-06-13 keyboard-lock**.

**Reference:** spec `docs/superpowers/specs/2026-06-20-app-shell-unification-design.md`; recon dossier facts inline below (file:line from the 8-slice recon, verified against live source).

**Test-helper note:** backend test helpers in this plan are real and verified — `fresh_db()` + `write_prp(path; …)` (`test/test_ingestion_core.jl`), `with_inproc_routes(db) do call … end` where `call(method, target; headers, body)` dispatches via `Oxygen.internalrequest` (`test/test_inproc.jl:40`), `seed_experiment!(db, exp_dir; name, analysis_dir, …)` (`test/seed.jl:31`), `create_experiment!`/`create_sample!`/`create_exposure!` (`db.jl`). **Frontend** test helpers shown in code blocks (`renderSamplesPage`, `mockNavigate`, `activeSampleIndex`, `renderWithRouter`, etc.) are **illustrative** — map them to this project's existing RTL setup (`frontend/test/` + the helpers in `frontend/test/AGENTS.md`); the assertion *intent* is what's normative, not the helper name.

---

## Milestone 1 — Backend ingest API (pure Julia; ships alone)

New file `src/routes_fs.jl` (`validate`, `suggest`, `manifest`), pattern-parsing in the create route, and an `on_progress` callback on `scan_and_group!`. All backend; verified by the `routes`/`pipeline` test buckets.

### Task 1.1: `scan_and_group!` gains an `on_progress` callback

**Files:**
- Modify: `packages/HimalayaUI/src/ingest.jl:45-52` (signature), `:186-194` (analysis loop)
- Test: `packages/HimalayaUI/test/test_ingestion_core.jl` (canonical scan fixture home — defines `fresh_db()` and `write_prp(path; …)`, kwargs at `:20`; many `scan_and_group!` testsets already build a `data_dir` with `touch("$stem.tif")` + `write_prp`)

**Interfaces:**
- Produces: `scan_and_group!(db, experiment_id; analyze=true, tif_pattern, prp_pattern, dat_pattern, on_progress::Union{Function,Nothing}=nothing)`. `on_progress` is called as `on_progress(processed::Int, total::Int)` after each `analyze_exposure!`, where `total == length(new_exposure_ids)`. Fires from the analysis loop **outside** the structural transaction; called once per new exposure even if that exposure's `analyze_exposure!` throws (it's caught). With `analyze=false` it does not fire (nothing to stream).

- [ ] **Step 1: Write the failing test** (mirror an existing `scan_and_group!` fixture in this file; tif+prp only — `analyze_exposure!` warns on the missing `.dat` but `on_progress` still ticks)

```julia
@testset "scan_and_group! on_progress reports per-exposure" begin
    db = fresh_db()
    data_dir = mktempdir()
    for stem in ("e1", "e2", "e3")
        touch(joinpath(data_dir, "$stem.tif"))
        write_prp(joinpath(data_dir, "$stem.prp"))      # default kwargs OK; see :20
    end
    exp_id = create_experiment!(db; name="t", path=data_dir, data_dir=data_dir, analysis_dir=data_dir)
    ticks = Tuple{Int,Int}[]
    scan_and_group!(db, exp_id; analyze=true, on_progress=(p, t) -> push!(ticks, (p, t)))
    @test last(ticks) == (3, 3)           # final tick = all done
    @test issorted(first.(ticks))         # monotonic processed count
    @test all(t -> t[2] == 3, ticks)      # total stable
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `GROUP=pipeline julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' 2>&1 | grep -A3 on_progress`
(If `test_ingestion_core.jl` is in a different bucket, run that bucket — confirm via the bucket map in `test/runtests.jl`.)
Expected: FAIL — `on_progress` is not a keyword of `scan_and_group!`.

- [ ] **Step 3: Add the kwarg and the throttled call**

In `ingest.jl`, add `on_progress::Union{Function,Nothing}=nothing` to the keyword list (after `dat_pattern`). Then replace the analysis loop (`:186-194`):

```julia
    if analyze
        total = length(new_exposure_ids)
        for (i, eid) in enumerate(new_exposure_ids)
            try
                analyze_exposure!(db, eid)
            catch e
                @warn "scan_and_group!: analyze_exposure! failed" exposure_id=eid exception=e
            end
            # Progress fires OUTSIDE the structural txn (already committed above).
            # ponytail: tick every exposure; the SSE 64-cap drops surplus, terminal frame is authoritative.
            on_progress === nothing || on_progress(i, total)
        end
    end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `GROUP=pipeline julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' 2>&1 | grep -E "on_progress|Test Summary"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/ingest.jl packages/HimalayaUI/test/test_ingestion_core.jl
git commit   # message: "feat(ingest): on_progress callback on scan_and_group!" + trailers
```

### Task 1.2: Create route persists patterns + wires `on_progress`

> **Why this shape (verified against live source):** `scan_and_group!` *already* reads the experiment row's `image_pattern`/`metadata_pattern`/`integration_pattern` columns and coalesces them over its kwargs (`ingest.jl` head). So the create route does NOT pass patterns to the scan — it **writes them to the row at create time**, and the existing `scan_and_group!(db, exp_id)` call picks them up. But `create_experiment!` (`db.jl:2103`) does not currently accept or INSERT those columns, so we extend it. The frontend body nests them as `patterns: { image, metadata, integration }` (`api.ts:CreateExperimentBody`), NOT flat fields.

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl:2103` (`create_experiment!` — add 3 kwargs + INSERT columns; columns already exist in schema at `:726-728`)
- Modify: `packages/HimalayaUI/src/routes_experiments.jl:128-174` (parse nested `patterns`, pass to create), `:176-203` (async spawn → wire `on_progress`)
- Test: `packages/HimalayaUI/test/test_ingestion_scan_api.jl` (routes-style scan tests; uses the `with_inproc_routes(db) do call … end` dispatch idiom)

**Interfaces:**
- Consumes: `scan_and_group!(...; on_progress=...)` (Task 1.1); `broadcast_progress!(experiment_id; kind, processed, total)` (`events.jl:1167`, exists).
- Produces: `create_experiment!(db; …, image_pattern::Union{String,Nothing}=nothing, metadata_pattern=nothing, integration_pattern=nothing)`; `POST /api/experiments` reads `body.patterns.{image,metadata,integration}` and persists them to the row; the async scan emits `ingest_progress` frames.

- [ ] **Step 1: Write the failing test** (nested patterns body → persisted column; in-process dispatch idiom)

```julia
@testset "POST /api/experiments persists patterns from the nested body" begin
    db = fresh_db()
    with_inproc_routes(db) do call
        dir = mktempdir(); touch(joinpath(dir, "s.tif"))
        body = Vector{UInt8}(JSON3.write(Dict(
            "path"     => dir,
            "patterns" => Dict("image" => "{name}.tiff"))))
        resp = call("POST", "/api/experiments";
            headers = ["Content-Type" => "application/json"], body = body)
        @test resp.status == 202
        id  = JSON3.read(resp.body).id
        row = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT image_pattern FROM experiments WHERE id = ?", [id])))
        @test row.image_pattern == "{name}.tiff"
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `GROUP=routes julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' 2>&1 | grep -A3 "persists patterns"`
Expected: FAIL — `image_pattern` is `NULL` (route ignores the body field; `create_experiment!` can't store it).

- [ ] **Step 3a: Extend `create_experiment!`**

In `db.jl:2103`, add three kwargs (after `experiment_type`):

```julia
        image_pattern::Union{String,Nothing} = nothing,
        metadata_pattern::Union{String,Nothing} = nothing,
        integration_pattern::Union{String,Nothing} = nothing,
```

Add `image_pattern, metadata_pattern, integration_pattern` to the INSERT column list (after `experiment_type`), add three `?` placeholders to the `VALUES (…)` tuple, and add the three variables to the params array in the same position. (`NULL` when omitted — identical to today's behaviour.)

- [ ] **Step 3b: Parse nested patterns in the create route + pass to create + wire progress**

After the `name`/`analysis_dir` derivation (`routes_experiments.jl:~161`), add:

```julia
        # Patterns are edited on Configuration before Approve; persisted to the row
        # so scan_and_group! (which reads the row + coalesces) uses them. Nested
        # under `patterns` to match the frontend CreateExperimentBody.
        pats = get(body, :patterns, get(body, "patterns", Dict()))
        ppat(k) = (v = get(pats, k, get(pats, string(k), nothing)); v === nothing ? nothing : String(v))
```

Pass `image_pattern = ppat(:image)`, `metadata_pattern = ppat(:metadata)`, `integration_pattern = ppat(:integration)` into the `create_experiment!(db; …)` call (`:165-170`). Then in the async block (`:181`), wire progress onto the existing scan call:

```julia
                scan_and_group!(db, exp_id;
                    on_progress = (p, t) -> broadcast_progress!(exp_id;
                        kind = "ingest_progress", processed = p, total = t))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `GROUP=routes julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' 2>&1 | grep -E "persists patterns|Test Summary"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/src/routes_experiments.jl packages/HimalayaUI/test/test_ingestion_scan_api.jl
git commit   # "feat(routes): persist create patterns + emit ingest_progress" + trailers
```

### Task 1.3: `routes_fs.jl` — `GET /api/fs/suggest` (path autocomplete)

**Files:**
- Create: `packages/HimalayaUI/src/routes_fs.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl:20-34` (add `include("routes_fs.jl")` before `include("server.jl")`), `packages/HimalayaUI/src/server.jl` `register_routes!` (`:49-133`, add `register_fs_routes!()` near the other `register_*` calls)
- Test: `packages/HimalayaUI/test/test_routes_fs.jl` (new; add to the routes bucket include list — confirm where buckets are declared in `test/runtests.jl`)

**Interfaces:**
- Produces: `register_fs_routes!()`; `GET /api/fs/suggest?prefix=<str>` → `{ "suggestions": ["/abs/path/a", "/abs/path/b", ...] }` (directories only, sorted, capped at 20). Matches `api.ts:suggestPaths` (`:240-244`).

- [ ] **Step 1: Write the failing test**

```julia
@testset "GET /api/fs/suggest lists child dirs of a prefix" begin
    db = fresh_db()
    with_inproc_routes(db) do call
        root = mktempdir(); mkpath(joinpath(root, "alpha")); mkpath(joinpath(root, "beta"))
        touch(joinpath(root, "afile.txt"))
        resp = call("GET", "/api/fs/suggest?prefix=$(HTTP.escapeuri(joinpath(root, "a")))")
        @test resp.status == 200
        s = JSON3.read(resp.body).suggestions
        @test any(endswith.(s, "alpha"))
        @test !any(endswith.(s, "afile.txt"))   # files excluded
        @test !any(endswith.(s, "beta"))         # prefix-filtered
    end
end
```

> The in-process harness (`with_inproc_routes` → `_ensure_inproc_routes!`) registers routes via the same `register_routes!` path the server uses — so adding `register_fs_routes!()` there (this task's Step 3) makes `/api/fs/*` dispatchable in-process. Confirm `_ensure_inproc_routes!` calls `register_routes!` (or add the fs registration where it registers the others).

- [ ] **Step 2: Run test to verify it fails**

Run: `GROUP=routes julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' 2>&1 | grep -A3 "fs/suggest"`
Expected: FAIL — 404 (route unregistered).

- [ ] **Step 3: Create `routes_fs.jl` with the suggest route**

Mirror the query-param + `current_db()`/`HTTP.Response`/`JSON3.write` pattern from `routes_picker.jl:33-59` and the error shape from `json.jl:_json_error`.

```julia
using HTTP, JSON3, Oxygen

"""
    register_fs_routes!()

Filesystem probes for the ingest funnel (spec §6.1). Read-only, no DB writes.
- GET /api/fs/suggest?prefix=  → directory autocomplete for the picker.
- GET /api/fs/validate?path=   → cheap picker gate (exists + not-already-an-experiment).
- GET /api/fs/manifest?path=&{image,metadata,integration}_pattern=  → phase-1 file manifest.
"""
function register_fs_routes!()
    @get "/api/fs/suggest" function (req::HTTP.Request)
        q      = HTTP.queryparams(HTTP.URI(req.target))
        prefix = get(q, "prefix", "")
        isempty(prefix) && return _json(200, Dict(:suggestions => String[]))
        dir  = isdir(prefix) ? prefix : dirname(prefix)
        base = isdir(prefix) ? "" : basename(prefix)
        isdir(dir) || return _json(200, Dict(:suggestions => String[]))
        kids = String[]
        for name in readdir(dir; sort = true)
            startswith(name, base) || continue
            full = joinpath(dir, name)
            isdir(full) && push!(kids, full)
            length(kids) >= 20 && break
        end
        _json(200, Dict(:suggestions => kids))
    end
    # validate + manifest added in Tasks 1.4 / 1.5
end
```

> **Reuse `_json(status, body)`** — it already exists at module scope (`routes_resolve.jl`); do NOT define a local one (ponytail review). It serializes + sets the content-type header.

Add `include("routes_fs.jl")` in `HimalayaUI.jl` (before `server.jl`) and call `register_fs_routes!()` inside `register_routes!` in `server.jl`. **Register the new test file in BOTH drift guards** (test/AGENTS.md): add `test_routes_fs.jl` to the `routes` GROUP bucket AND to the `ALL_ORDER` list (and the `All` GROUP entry) in `test/runtests.jl`, or the suite errors at load.

- [ ] **Step 4: Run test to verify it passes**

Run: `GROUP=routes julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' 2>&1 | grep -E "fs/suggest|Test Summary"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_fs.jl packages/HimalayaUI/src/HimalayaUI.jl packages/HimalayaUI/src/server.jl packages/HimalayaUI/test/test_routes_fs.jl packages/HimalayaUI/test/runtests.jl
git commit   # "feat(routes): GET /api/fs/suggest" + trailers
```

### Task 1.4: `GET /api/fs/validate` (cheap picker gate)

**Files:**
- Modify: `packages/HimalayaUI/src/routes_fs.jl` (add the route inside `register_fs_routes!`)
- Test: `packages/HimalayaUI/test/test_routes_fs.jl`

**Interfaces:**
- Produces: `GET /api/fs/validate?path=<dir>` → `{ "ok": bool, "matched": int, "scanned": int, "message": str|null }`. `ok` = `isdir(path) && no experiment already uses it`. Matches `api.ts:validatePath` (`:248-256`). `matched`/`scanned` are a cheap file count (0/0 acceptable pre-pattern; the rich count is `manifest`'s job).

- [ ] **Step 1: Write the failing test**

```julia
@testset "GET /api/fs/validate gates the picker" begin
    db = fresh_db()
    with_inproc_routes(db) do call
        good = mktempdir(); touch(joinpath(good, "x.tif"))
        r1 = call("GET", "/api/fs/validate?path=$(HTTP.escapeuri(good))")
        @test JSON3.read(r1.body).ok == true

        r2 = call("GET", "/api/fs/validate?path=$(HTTP.escapeuri(good * "_nope"))")
        body2 = JSON3.read(r2.body)
        @test body2.ok == false
        @test body2.message !== nothing

        # already-an-experiment → not ok
        create_experiment!(db; name="e", path=good, data_dir=good, analysis_dir=good)
        @test JSON3.read(call("GET", "/api/fs/validate?path=$(HTTP.escapeuri(good))").body).ok == false
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `GROUP=routes julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' 2>&1 | grep -A3 "fs/validate"`
Expected: FAIL — 404.

- [ ] **Step 3: Add the validate route**

Inside `register_fs_routes!`, before the trailing comment:

```julia
    @get "/api/fs/validate" function (req::HTTP.Request)
        q    = HTTP.queryparams(HTTP.URI(req.target))
        path = get(q, "path", "")
        if isempty(path) || !isdir(path)
            return _json(200, Dict(:ok => false, :matched => 0, :scanned => 0,
                              :message => "path does not exist or is not a directory"))
        end
        dup = !isempty(DBInterface.execute(current_db(),
            "SELECT 1 FROM experiments WHERE data_dir = ? LIMIT 1", [path]) |> Tables.rowtable)
        if dup
            return _json(200, Dict(:ok => false, :matched => 0, :scanned => 0,
                              :message => "an experiment already uses this directory"))
        end
        scanned = count(!startswith("."), readdir(path))   # cheap; rich count is /manifest
        _json(200, Dict(:ok => true, :matched => scanned, :scanned => scanned, :message => nothing))
    end
```

(`validate` is kept — it serves the picker's pre-Approve "exists + not-already-an-experiment" gate, a different funnel step than `manifest`, and `api.ts:validatePath` already targets it.) Add `using DBInterface, Tables` to the file head if not already pulled in transitively (match the other routes files).

- [ ] **Step 4: Run test to verify it passes**

Run: `GROUP=routes julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' 2>&1 | grep -E "fs/validate|Test Summary"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_fs.jl packages/HimalayaUI/test/test_routes_fs.jl
git commit   # "feat(routes): GET /api/fs/validate picker gate" + trailers
```

### Task 1.5: `GET /api/fs/manifest` (phase-1 file manifest)

**Files:**
- Modify: `packages/HimalayaUI/src/routes_fs.jl`
- Reuse: `scan_directory` (`grouping.jl:41-84`), `_matches_prefix_with_boundary` (`config.jl:168-189`), `resolve_files` (`config.jl:202-226`)
- Test: `packages/HimalayaUI/test/test_routes_fs.jl`

**Interfaces:**
- Produces: `GET /api/fs/manifest?path=<dir>&image_pattern=&metadata_pattern=&integration_pattern=` → `{ total, matched: {image, metadata, integration}, unmatched: [{file, miss: "metadata"|"integration"}] }`. **Patterns come from query params, NOT the DB** (phase-1 is pre-experiment). No PRP parse, no DB write. (No `nearest` field — YAGNI per ponytail review; the stem + miss-type is the actionable fact. Add a near-miss suggester later only if ScanFailedPage proves it wants one.)

- [ ] **Step 1: Write the failing test**

```julia
@testset "GET /api/fs/manifest counts by type without writing" begin
    db = fresh_db()
    with_inproc_routes(db) do call
        dir = mktempdir()
        touch(joinpath(dir, "s1.tif")); touch(joinpath(dir, "s1.prp")); touch(joinpath(dir, "s1.dat"))
        touch(joinpath(dir, "s2.tif"))   # missing prp + dat
        q = "path=$(HTTP.escapeuri(dir))&image_pattern={name}.tif&metadata_pattern={name}.prp&integration_pattern={name}.dat"
        resp = call("GET", "/api/fs/manifest?$q")
        @test resp.status == 200
        m = JSON3.read(resp.body)
        @test m.matched.image == 2
        @test m.matched.metadata == 1
        @test any(u -> u.miss == "metadata", m.unmatched)   # s2 missing prp
        @test isempty(Tables.rowtable(DBInterface.execute(db, "SELECT 1 FROM experiments")))  # no row created
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `GROUP=routes julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' 2>&1 | grep -A3 "fs/manifest"`
Expected: FAIL — 404.

- [ ] **Step 3: Add the manifest route**

`scan_directory` returns `Vector{ExposureMeta}` but also parses PRP — for phase-1 we want listing only, so do a direct `readdir` + the existing pattern matcher. Add inside `register_fs_routes!`:

```julia
    @get "/api/fs/manifest" function (req::HTTP.Request)
        q    = HTTP.queryparams(HTTP.URI(req.target))
        path = get(q, "path", "")
        isdir(path) || return _json(400, Dict(:error => "path is not a directory"))
        pats = (image       = get(q, "image_pattern", "{name}.tif"),
                metadata    = get(q, "metadata_pattern", "{name}.prp"),
                integration = get(q, "integration_pattern", "{name}.dat"))
        files = filter(!startswith("."), readdir(path))
        # {name}-capture per type, inlined. Do NOT add a module helper — grouping.jl
        # / config.jl already own pattern matching; if the manifest's needs grow,
        # call the shared matcher (config.jl `_matches_prefix_with_boundary` /
        # `resolve_files`) rather than maintaining a second one (ponytail review).
        function stems(pat)
            occursin("{name}", pat) || return Set{String}()
            pre, post = split(pat, "{name}"; limit = 2)
            out = Set{String}()
            for f in files
                (startswith(f, pre) && endswith(f, post) &&
                    length(f) > length(pre) + length(post)) || continue
                push!(out, f[nextind(f, lastindex(pre)):prevind(f, lastindex(f) - length(post) + 1)])
            end
            out
        end
        img, meta, integ = stems(pats.image), stems(pats.metadata), stems(pats.integration)
        unmatched = Dict{String,String}[]
        for s in img, (label, set) in (("metadata", meta), ("integration", integ))
            s in set || push!(unmatched, Dict("file" => s, "miss" => label))
        end
        _json(200, Dict(:total => length(files),
                   :matched => Dict(:image => length(img), :metadata => length(meta), :integration => length(integ)),
                   :unmatched => unmatched))
    end
```

> **No helper functions** — the `stems` closure is local to the route. `_nearest_file` is cut (YAGNI). If pattern matching gets more complex, call `config.jl`'s `_matches_prefix_with_boundary` / `resolve_files` — do not maintain a second matcher.

- [ ] **Step 4: Run test to verify it passes**

Run: `GROUP=routes julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' 2>&1 | grep -E "fs/manifest|Test Summary"`
Expected: PASS.

- [ ] **Step 5: Full backend gate + commit**

Run: `make test-parallel 2>&1 | tail -20` — expect all 5 buckets green.

```bash
git add packages/HimalayaUI/src/routes_fs.jl packages/HimalayaUI/test/test_routes_fs.jl
git commit   # "feat(routes): GET /api/fs/manifest phase-1 file manifest" + trailers
```

**M1 done:** backend exposes the three fs probes, the create route honours patterns, and ingest emits per-exposure progress. The frontend stubs (`validatePath`/`suggestPaths`) now hit live endpoints.

---

## Milestone 2 — Unified keyboard + sheet interaction (frontend; ships alone)

Re-axis the registry, remove the roving grid, add the page-level cursor + per-row Restore. Works on the current `/samples` sheet *before* the shell collapse, so it's independently verifiable.

### Task 2.1: Re-axis the shortcut registry + `?` normalization

> **Real registry shape (verified live, `shortcuts.ts`):** `ShortcutId` is a **union type** (`:10-25`); `SHORTCUTS` is `Record<ShortcutId, ShortcutDef>` where `ShortcutDef = { id, keys: string[], label, group }` — each entry holds a **`keys` array** of normalized combos (e.g. `find` = `["Mod+k", "/"]`), NOT a `.combo` string. `matchShortcut(e)` returns a **`ShortcutId | null` (a string)**, not an object. Consumers bind handlers via **`useShortcuts({ [id]: () => … })`** maps (one per page) — there is NO global window `keydown` switch, and you must NOT add one. The re-axis is therefore mostly a registry edit + an id migration; the page handlers are edited in T2.4/T2.5.

**Files:**
- Modify: `frontend/src/print/shell/shortcuts.ts` — `ShortcutId` union (`:10-25`), `SHORTCUTS` record (`:38-62`), the locked-decision comment (`:6-8`, now false), `eventCombo` (`:67-82`)
- Test: `frontend/test/shortcuts.test.ts` (mirror the existing CapsLock-X normalization test)

**Interfaces — the id migration (old → new):**
| old id | change | new `keys` |
|---|---|---|
| `prevSample` / `nextSample` | keep id, remap keys (was `[`/`]`) | `["ArrowUp"]` / `["ArrowDown"]` |
| `prevCandidate` / `nextCandidate` | **rename** → `prevDetail` / `nextDetail` | `["ArrowLeft"]` / `["ArrowRight"]` |
| `prevExposure` / `nextExposure` | **remove** (exposure-step → filmstrip `onSelect`, T2.5) | — |
| `drop` `keep` `representative` `undo` `redo` `reorderUp` `reorderDown` `dismiss` `find` | unchanged | unchanged |
| — | **add** `openFocus` | `["Enter"]` |
| — | **add** `openLoupe` | `["l"]` |
| — | **add** `toggleSelect` | `[" "]` (Space; `eventCombo` lowercases the single-char `" "`) |
| — | **add** `extendPrev` / `extendNext` | `["Shift+ArrowLeft"]` / `["Shift+ArrowRight"]` |
| — | **add** `selectAll` | `["Mod+a"]` |
| — | **add** `restore` | `["Backspace"]` |
| — | **add** `helpOverlay` | `["?"]` |

`eventCombo` emits the literal token `?` for the `?` key regardless of the Shift bit. "Detail" is page-interpreted (frame on Corpus/Loupe, candidate on Focus) — same id, different per-page handler, which the `useShortcuts({id})` model supports natively.

**Consumer updates (this task makes the registry compile; the handler bindings move in T2.4/T2.5):** every reference to a removed/renamed id must change — `FocusPage.tsx` (binds `prevExposure`/`nextExposure`/`prevCandidate`/`nextCandidate`), `LoupePage.tsx` (binds `prevExposure`/`nextExposure`), `KbdLegend.tsx` (renders the registry — auto-updates, but set `label`/`group` for the new ids), any `aria-keyshortcuts`. End the task with `npx tsc --noEmit` showing **zero** references to `prevExposure`/`nextExposure`/`prevCandidate`/`nextCandidate`.

- [ ] **Step 1: Write the failing tests** (real shape: `Object.values(SHORTCUTS)`, `.keys`, string `matchShortcut`)

```ts
import { eventCombo, SHORTCUTS, matchShortcut } from '../src/print/shell/shortcuts'

const ev = (init: Partial<KeyboardEvent>) => new KeyboardEvent('keydown', init)
const keysOf = (id: string) => Object.values(SHORTCUTS).find(d => d.id === id)?.keys ?? []

test('? normalizes to a stable token regardless of Shift', () => {
  expect(eventCombo(ev({ key: '?', shiftKey: true }))).toBe('?')
})

test('arrows are the sample/detail axis; [ and ] are gone', () => {
  expect(matchShortcut(ev({ key: 'ArrowUp' }))).toBe('prevSample')   // string id, not object
  expect(matchShortcut(ev({ key: 'ArrowLeft' }))).toBe('prevDetail')
  const allKeys = Object.values(SHORTCUTS).flatMap(d => d.keys)
  expect(allKeys).not.toContain('[')
  expect(allKeys).not.toContain(']')
})

test('new verbs are bound', () => {
  expect(keysOf('openFocus')).toContain('Enter')
  expect(keysOf('restore')).toContain('Backspace')
  expect(matchShortcut(ev({ key: ' ' }))).toBe('toggleSelect')
})

test('every combo resolves to exactly one id (no key bound twice)', () => {
  const seen = new Map<string, number>()
  for (const d of Object.values(SHORTCUTS)) for (const k of d.keys) seen.set(k, (seen.get(k) ?? 0) + 1)
  for (const [k, n] of seen) expect(`${k}:${n}`).toBe(`${k}:1`)
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `cd packages/HimalayaUI/frontend && npx vitest run shortcuts.test.ts`
Expected: FAIL — `prevSample` still on `[`; `prevDetail`/`openFocus`/`restore`/`toggleSelect` undefined.

- [ ] **Step 3: Edit the union + record + comment + `eventCombo`**

Apply the migration table to the `ShortcutId` union and the `SHORTCUTS` record (`keys` arrays + `label`/`group` for the new ids). Rewrite the stale "Locked decisions (Jonathan 2026-06-13)" comment (`:6-8`) to state the rev-2 axes and that this supersedes the lock. In `eventCombo`, as the first line of the body:

```ts
  // '?' is Shift+/ on US layouts but layout-variable; emit a stable token so the
  // help binding is layout-robust (mirrors the single-char lowercasing below).
  if (e.key === '?') return '?'
```

- [ ] **Step 4: Run to verify it passes + tsc gate**

Run: `npx vitest run shortcuts.test.ts && npx tsc --noEmit 2>&1 | grep -E "prevExposure|nextExposure|prevCandidate|nextCandidate" || echo "no dangling refs"`
Expected: tests PASS; the grep prints `no dangling refs` (T2.4/T2.5 will rebind; if tsc still flags page handlers, that's expected until those tasks — note which files for the next tasks).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/shell/shortcuts.ts packages/HimalayaUI/frontend/test/shortcuts.test.ts
git commit   # "feat(keys): re-axis registry to sample/detail + new verb ids" + trailers
```

> **Sequencing note:** Steps 3–5 leave `FocusPage`/`LoupePage` referencing the renamed/removed ids (tsc errors there). Do T2.4/T2.5 in the same work session so the tree is green before any push; each of those tasks ends green on its own page.

### Task 2.2: Remove the roving grid

**Files:**
- Remove: `frontend/src/lib/grid/useRovingGrid.ts`, `frontend/src/lib/grid/rovingGrid.ts`
- Modify: `frontend/src/print/components/SheetTable.tsx:40-45,139-153`, `frontend/src/print/components/SampleTableRow.tsx:12-61,217-294`, `frontend/src/print/pages/SamplesPage.tsx:438-444,470-505`
- Test: `frontend/test/SheetTable.test.tsx` (existing; update roving assertions)

**Interfaces:**
- Produces: `SheetTable` no longer accepts `roving`/`dataRowCount`/`focusOnMountRow`; renders a plain `<table>` (role row/cell). `SampleTableRow` drops `rowIndex` + all `useRovingGridContext` wiring. `SAMPLE_TABLE_COLS*` / `sampleTableCols()` (`:81-97`) unchanged (alignment invariant).

- [ ] **Step 1: Update the test to assert plain-table semantics**

```tsx
test('SheetTable renders a plain table (no roving grid)', () => {
  render(<SheetTable head={<tr><th>H</th></tr>}><tr><td>r</td></tr></SheetTable>)
  expect(screen.queryByRole('grid')).toBeNull()
  expect(screen.getByRole('table')).toBeInTheDocument()
})
```

Delete any existing test that asserts `role="grid"` / roving tabindex.

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run SheetTable.test.tsx`
Expected: FAIL — current SheetTable emits `role="grid"`.

- [ ] **Step 3: Delete the modules and strip the wiring**

```bash
git rm packages/HimalayaUI/frontend/src/lib/grid/useRovingGrid.ts packages/HimalayaUI/frontend/src/lib/grid/rovingGrid.ts
```

In `SheetTable.tsx`: remove the `roving`/`dataRowCount`/`focusOnMountRow` props from the interface, delete the `useRovingGrid(...)` call and the `RovingGridProvider` wrap, and render `<table>` (drop `role="grid"` + `onContainerKeyDown`). In `SampleTableRow.tsx`: delete the `:217-294` block (`useRovingGridContext` + `cellRole`/`tab`/`cellRef`/`activatePointer`/`cellKeyDown` + multi-widget inertness), remove the `rowIndex` prop (`:12-61`), and let cells use default `role`. In `SamplesPage.tsx`: drop `roving`/`dataRowCount`/`focusOnMountRow` (`:438-444`) and `rowIndex={i+1}` (`:470-505`); keep `onOpenFocus`/`onSelectExposure`.

- [ ] **Step 4: Run to verify it passes**

Run: `npx vitest run SheetTable.test.tsx SampleTableRow.test.tsx` then `npx tsc --noEmit` and a text sweep for stragglers: `git grep -nE "roving|useRovingGrid|RovingGrid|requestActivate|onContainerKeyDown|focusOnMountRow|dataRowCount" -- 'frontend/src' 'frontend/test'`.
Expected: tests PASS; the grep returns nothing (the 3 files are the *modules*, but props/comments/tests reference the removed API too — clear them all, incl. any `roving`/`rowIndex` props passed by callers and stale comments).

- [ ] **Step 5: Commit**

```bash
git add -u packages/HimalayaUI/frontend/src/print/components/SheetTable.tsx packages/HimalayaUI/frontend/src/print/components/SampleTableRow.tsx packages/HimalayaUI/frontend/src/print/pages/SamplesPage.tsx packages/HimalayaUI/frontend/test
git commit   # "refactor(sheet): remove roving grid → plain table" + trailers
```

### Task 2.3: Per-row Restore button on dropped rows

**Files:**
- Modify: `frontend/src/print/components/KeptCell.tsx:1-45`
- Test: `frontend/test/KeptCell.test.tsx` (new)

**Interfaces:**
- Consumes: an `onRestore?: () => void` prop on `KeptCell` (wired by `SampleTableRow`).
- Produces: when `dropped > 0` and `onRestore` is provided, renders a real `<button>` (use the `Button` primitive, `variant` matching the existing dock-neutral look) labelled `Restore`. This is the pointer target for the `Backspace`=restore key.

- [ ] **Step 1: Write the failing test**

```tsx
test('KeptCell shows a Restore button only when something is dropped', () => {
  const onRestore = vi.fn()
  const { rerender } = render(<KeptCell kept={3} dropped={0} total={3} onRestore={onRestore} />)
  expect(screen.queryByRole('button', { name: /restore/i })).toBeNull()
  rerender(<KeptCell kept={2} dropped={1} total={3} onRestore={onRestore} />)
  fireEvent.click(screen.getByRole('button', { name: /restore/i }))
  expect(onRestore).toHaveBeenCalledOnce()
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run KeptCell.test.tsx`
Expected: FAIL — no Restore button exists.

- [ ] **Step 3: Add the button (appearance via the Button primitive)**

In `KeptCell.tsx`, after the kept/dropped counts, when `dropped > 0 && onRestore`:

```tsx
import { Button } from '../ui/Button'
// ...
{dropped > 0 && onRestore && (
  <Button variant="neutral" className="ml-2" onClick={(e) => { e.stopPropagation(); onRestore() }}>
    Restore
  </Button>
)}
```

(Confirm the exact neutral variant name in `ui/Button.tsx`; the dossier lists Button variants incl. `ghostInverse` — use the plain neutral one, NOT ghostInverse which is dark-surface.) Wire `onRestore` from `SampleTableRow` → the sample's existing restore callback (mirror the verb logic used by `CullBar`'s Restore).

- [ ] **Step 4: Run to verify it passes**

Run: `npx vitest run KeptCell.test.tsx`
Expected: PASS.

- [ ] **Step 5: Lint-design + commit**

Run: `npm run lint:design` (must pass — appearance came from the primitive).

```bash
git add packages/HimalayaUI/frontend/src/print/components/KeptCell.tsx packages/HimalayaUI/frontend/src/print/components/SampleTableRow.tsx packages/HimalayaUI/frontend/test/KeptCell.test.tsx
git commit   # "feat(sheet): per-row Restore on dropped rows" + trailers
```

### Task 2.4: Page-level cursor + Corpus keyboard map (via `useShortcuts`)

> **Architecture (per review):** bind through the existing `useShortcuts({ [id]: handler })` map — do NOT hand-roll a window `keydown` switch (that reintroduces the double-dispatch this redesign kills). No `Alt`-guard is needed: `reorderUp` = `"Alt+ArrowUp"` is a distinct combo from `prevSample` = `"ArrowUp"`, so the registry never confuses them, and `useShortcuts` already suppresses inside inputs/modals.

**Files:**
- Modify: `frontend/src/print/pages/SamplesPage.tsx` — selection-state region (`:200-289`), keyboard block (`:317-339`)
- Test: `frontend/test/SamplesPage.keyboard.test.tsx` (new)

**Interfaces:**
- Consumes: registry ids from T2.1; the `useShortcuts` hook.
- Produces: page-level cursor `{ sampleIndex, frameIndex }` (default `{0,0}`), driven by a `useShortcuts` map: `prevSample`/`nextSample` clamp-step the sample (reset `frameIndex` 0), `prevDetail`/`nextDetail` clamp-step the frame within the active sample, `openFocus` → navigate Focus (`/sample/:id`, flat per spec §6.4), `openLoupe` → `/sample/:id/loupe`, `toggleSelect`/`extendPrev`/`extendNext`/`selectAll` drive the existing selection (`selected`/`anchorRef`/`checkedSamples`), `drop`/`keep`/`restore` apply to the selection-else-active-frame, `dismiss` clears the selection then **returns false** (so the Esc ladder continues). **No `representative` binding.** Pointer click on a row/exposure calls the same `setCursor`.

- [ ] **Step 1: Write the failing tests** (helper names illustrative — map to `frontend/test/` setup)

```tsx
test('Arrows drive the page cursor; Alt+Arrow is a separate combo (no cursor move)', () => {
  renderSamplesPage(/* 3 samples, sample0 has 2 frames */)
  fireEvent.keyDown(window, { key: 'ArrowDown' })
  expect(activeSampleIndex()).toBe(1)
  fireEvent.keyDown(window, { key: 'ArrowUp', altKey: true })   // reorderUp id — not bound here
  expect(activeSampleIndex()).toBe(1)                            // unchanged
})

test('Enter opens Focus for the active sample; r does nothing on Corpus', () => {
  const nav = mockNavigate()
  renderSamplesPage()
  fireEvent.keyDown(window, { key: 'Enter' })
  expect(nav).toHaveBeenCalledWith(expect.stringMatching(/\/sample\//))
  fireEvent.keyDown(window, { key: 'r' })
  expect(setRepresentative).not.toHaveBeenCalled()
})
```

- [ ] **Step 2: Run to verify it fails** — `npx vitest run SamplesPage.keyboard.test.tsx` → FAIL (no cursor; Enter not wired).

- [ ] **Step 3: Add the cursor + the `useShortcuts` map**

```tsx
const [cursor, setCursor] = useState({ sampleIndex: 0, frameIndex: 0 })
const activeSample = samples[cursor.sampleIndex]
const clampSample = (d: number) =>
  setCursor((c) => ({ sampleIndex: clamp(c.sampleIndex + d, 0, samples.length - 1), frameIndex: 0 }))
const clampFrame = (d: number) =>
  setCursor((c) => ({ ...c, frameIndex: clamp(c.frameIndex + d, 0, (samples[c.sampleIndex]?.frames.length ?? 1) - 1) }))

useShortcuts({
  prevSample: () => clampSample(-1),
  nextSample: () => clampSample(1),
  prevDetail: () => clampFrame(-1),
  nextDetail: () => clampFrame(1),
  openFocus:  () => activeSample && navigate(`/sample/${activeSample.id}`),
  openLoupe:  () => activeSample && navigate(`/sample/${activeSample.id}/loupe`),
  toggleSelect: () => toggleInSelection(activeFrame()),
  extendPrev:   () => extendSelection(-1),
  extendNext:   () => extendSelection(1),
  selectAll:    () => selectAllFrames(),
  drop:    () => applyVerb('drop', selectionOrActive()),
  keep:    () => applyVerb('keep', selectionOrActive()),
  restore: () => applyVerb('restore', selectionOrActive()),
  dismiss: () => { if (hasSelection()) { clearSelection(); return undefined } return false },
  // NO representative on Corpus.
})
```

Wire row/exposure `onClick` → `setCursor(...)`. Reuse the existing selection helpers (`selected`/`anchorRef`) — don't add a new store.

- [ ] **Step 4: Run to verify it passes** — `npx vitest run SamplesPage.keyboard.test.tsx` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/pages/SamplesPage.tsx packages/HimalayaUI/frontend/test/SamplesPage.keyboard.test.tsx
git commit   # "feat(corpus): page cursor + useShortcuts keyboard map" + trailers
```

### Task 2.5: Focus + Loupe re-axis (edit their `useShortcuts` maps)

> **No axis swaps in handler bodies — just rename/remove map keys.** `prevSample`/`nextSample` keep their ids (T2.1 already remapped their keys to `↑/↓`), so those handlers don't change. The candidate handlers move from `prevCandidate`/`nextCandidate` → `prevDetail`/`nextDetail` (same bodies, now on `←/→`). The exposure handlers `prevExposure`/`nextExposure` are **deleted** — exposure-stepping becomes filmstrip-only via the existing `ThumbnailGallery` `onSelect` (the `:814-827` region is that gallery, a click control, NOT a keyboard control; the spec's "rail-scoped" exposure stepping = clicking the filmstrip). No `Alt`-guard needed (registry separates `Alt+Arrow`).

**Files:**
- Modify: `frontend/src/print/pages/FocusPage.tsx:485-524` (useShortcuts map), `frontend/src/print/pages/LoupePage.tsx:338-351` (useShortcuts map)
- Test: `frontend/test/FocusPage.keyboard.test.tsx`, `frontend/test/LoupePage.keyboard.test.tsx`

**Interfaces:**
- Produces: Focus map = `prevSample`/`nextSample` (↑/↓, unchanged bodies) · `prevDetail`/`nextDetail` (←/→, the old candidate steppers) · `openLoupe` (→ `/sample/:id/loupe`) · `dismiss` (unchanged Esc ladder; its return-to-corpus target is finalized in T3.2 with the `from=series` marker). **Removed:** `prevExposure`/`nextExposure`. Loupe map = `prevSample`/`nextSample` (↑/↓) · `prevDetail`/`nextDetail` (←/→, the old frame steppers, renamed from `prevExposure`/`nextExposure`) · `representative` (kept) · `drop`/`keep` (kept) · `restore` (new, Backspace).

- [ ] **Step 1: Write the failing tests**

```tsx
// FocusPage
test('Focus: ←/→ steps candidate, ↑/↓ steps sample, exposure keys are gone', () => {
  renderFocus()
  fireEvent.keyDown(window, { key: 'ArrowRight' }); expect(candidateIndex()).toBe(1)
  fireEvent.keyDown(window, { key: 'ArrowDown' });  expect(sampleNav).toHaveBeenCalled()
})
// LoupePage
test('Loupe: arrows = sample/frame; r sets representative; Backspace restores', () => {
  renderLoupe()
  fireEvent.keyDown(window, { key: 'ArrowRight' }); expect(frameIndex()).toBe(1)
  fireEvent.keyDown(window, { key: 'r' });          expect(setRepresentative).toHaveBeenCalled()
  fireEvent.keyDown(window, { key: 'Backspace' });  expect(restore).toHaveBeenCalled()
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run FocusPage.keyboard.test.tsx LoupePage.keyboard.test.tsx`
Expected: FAIL (and these pages currently won't typecheck — T2.1 renamed the ids they bind).

- [ ] **Step 3: Edit the two maps**

`FocusPage.tsx` useShortcuts (`:485-524`): **delete** the `prevExposure`/`nextExposure` entries; **rename** `prevCandidate` → `prevDetail` and `nextCandidate` → `nextDetail` (keep their bodies — `stepInList`/`previewIndexId`); keep `prevSample`/`nextSample`/`dismiss` as-is (the `addArmed` DECLINE in `dismiss` stays); **add** `openLoupe: () => navigate(\`/sample/\${sample.id}/loupe\`)`. Exposure stepping now relies on the `ThumbnailGallery` `onSelect` already wired in the detector panel.

`LoupePage.tsx` useShortcuts (`:338-351`): **rename** `prevExposure` → `prevDetail`, `nextExposure` → `nextDetail` (frame steppers); keep `prevSample`/`nextSample`/`representative`/`drop`/`keep`; **add** `restore: () => restoreVerdict(activeFrame())`.

- [ ] **Step 4: Run to verify it passes**

Run: `npx vitest run FocusPage.keyboard.test.tsx LoupePage.keyboard.test.tsx && npx tsc --noEmit 2>&1 | grep -E "prevExposure|prevCandidate" || echo clean`
Expected: tests PASS; `clean`.

- [ ] **Step 5: `?` overlay + full frontend gate + commit**

Add the `?` help overlay: `frontend/src/print/shell/KbdOverlay.tsx` (reuse `ModalShell` + render `KbdLegend`, which is already registry-driven — `KbdLegend.tsx:12-34`), opened by the `helpOverlay` id. Wire it in `useGlobalShortcuts.ts` alongside the existing `find` special-case (a global, app-wide binding). Test: opening renders the legend; `Esc` closes.

Run: `npm test && npm run build` (lint:design + tsc + vite must pass).

```bash
git add packages/HimalayaUI/frontend/src/print/pages/FocusPage.tsx packages/HimalayaUI/frontend/src/print/pages/LoupePage.tsx packages/HimalayaUI/frontend/src/print/shell/KbdOverlay.tsx packages/HimalayaUI/frontend/src/hooks/useGlobalShortcuts.ts packages/HimalayaUI/frontend/test
git commit   # "feat(keys): Focus/Loupe re-axis + ? overlay" + trailers
```

> **Spec §4.1 reconciliation (minor):** Focus `P` (add-peak) and Series `A`/`⌘Enter` (add/confirm) stay **surface-local button-driven controls** per spec §7 (in-surface command controls untouched) — they are NOT added to the registry, so the `?` overlay lists only the global/navigation keymap. The "one unified set" framing in §4 covers navigation + cull + select; the in-surface verb keys remain with their surfaces.

**M2 done:** one keyboard model across Corpus/Focus/Loupe, the roving grid is gone, per-row Restore exists, and `?` shows the live keymap.

---

## Milestone 3 — Shell collapse + contextual dock (frontend; ships alone)

One `TopNav`, one shell, the routing collapse, and the contextual bottom dock. Bigger blast radius — lands after the keyboard is stable.

### Task 3.1: Unified `TopNav` (replaces both topbars)

**Files:**
- Create: `frontend/src/print/shell/TopNav.tsx`
- Reference: `ExperimentTopNav.tsx:9-48` (2-item nav template), `CorpusTopbar.tsx:28-32,110-119` (wordmark/tabs)
- Test: `frontend/test/TopNav.test.tsx` (new)

**Interfaces:**
- Produces: `<TopNav />` — wordmark `Himalaya · SAXS` → `/experiments`; two section tabs `Experiments` (`/experiments`) and `Series` (`/series`), active-state from the router. **No** Samples tab, **no** beamtime chip, **no** ⚙ (Configuration rides the experiment header, §3.1).

- [ ] **Step 1: Write the failing test**

```tsx
test('TopNav: wordmark→/experiments, two tabs, no Samples/beamtime/gear', () => {
  renderWithRouter(<TopNav />, { route: '/experiments' })
  expect(screen.getByRole('link', { name: /himalaya/i })).toHaveAttribute('href', '/experiments')
  expect(screen.getByRole('link', { name: 'Experiments' })).toBeInTheDocument()
  expect(screen.getByRole('link', { name: 'Series' })).toBeInTheDocument()
  expect(screen.queryByText('Samples')).toBeNull()
  expect(screen.queryByRole('button', { name: /beamtime/i })).toBeNull()
})
```

- [ ] **Step 2: Run to verify it fails** — `npx vitest run TopNav.test.tsx` → FAIL (no file).

- [ ] **Step 3: Build `TopNav` by adapting `ExperimentTopNav`** (copy its 2-item structure; wordmark from `CorpusTopbar:28-32` pointed at `/experiments`; appearance via existing primitives — no inline appearance).

- [ ] **Step 4: Run to verify it passes** — `npx vitest run TopNav.test.tsx` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/shell/TopNav.tsx packages/HimalayaUI/frontend/test/TopNav.test.tsx
git commit   # "feat(shell): unified TopNav" + trailers
```

### Task 3.2: Routing collapse — one shell, `/samples` redirect, `from=series` return

> **Flat Focus/Loupe (spec §6.4, rev 2):** Focus/Loupe stay at `/sample/:id` and `/sample/:id/loupe`. NO `LegacySampleRedirect`, NO experiment-scoped nesting — a resolver would always flash a loading screen. The chrome reads the experiment from the loaded sample's `experiment_id`. The only redirect is `/samples` → `/experiments`.

**Files:**
- Modify: `frontend/src/print/shell/AppRoutes.tsx:94-161` (one shell; `/samples`→`/experiments`)
- Remove: `frontend/src/print/shell/CorpusShell.tsx`, `frontend/src/print/shell/CorpusTopbar.tsx`, `frontend/src/print/shell/ExperimentTopNav.tsx`
- Modify: `frontend/src/print/shell/ExperimentShell.tsx:33-195` (use `TopNav`; header/tabs move into page content per §3.2)
- Modify: `frontend/src/print/pages/FocusPage.tsx` (`dismiss` return target — see Interfaces), and the **Series-member open call site** (add `?from=series`)
- Test: `frontend/test/AppRoutes.test.tsx`, `frontend/test/FocusPage.return.test.tsx`

**Interfaces:**
- Produces: a single shell wrapping all routes (no `CorpusShell`/`CorpusTopbar`/`ExperimentTopNav`, no beamtime chip). `/samples` `<Navigate replace to="/experiments">`. Focus/Loupe stay flat. A Series member opens Focus with `?from=series` (query param). Focus `dismiss` (the Esc-ladder terminal rung that currently does `navigate("/samples")`) becomes: `from=series` → `navigate("/series")`, else → `navigate(\`/experiments/\${sample.experiment_id}/corpus\`)` (expId from the loaded sample). The dock up-link (T3.3) reads the same marker to label `‹ Series` vs `‹ Corpus`.

- [ ] **Step 1: Write the failing tests**

```tsx
test('/samples redirects to /experiments', async () => {
  renderApp({ route: '/samples' })
  await waitFor(() => expect(window.location.pathname).toBe('/experiments'))
})
test('Focus Esc returns to series when opened from series, else to its experiment corpus', () => {
  const nav = mockNavigate()
  renderFocus({ route: '/sample/42?from=series', sample: { id: 42, experiment_id: 7 } })
  fireEvent.keyDown(window, { key: 'Escape' })
  expect(nav).toHaveBeenLastCalledWith('/series')
  renderFocus({ route: '/sample/42', sample: { id: 42, experiment_id: 7 } })
  fireEvent.keyDown(window, { key: 'Escape' })
  expect(nav).toHaveBeenLastCalledWith('/experiments/7/corpus')
})
```

- [ ] **Step 2: Run to verify it fails** — `npx vitest run AppRoutes.test.tsx FocusPage.return.test.tsx` → FAIL.

- [ ] **Step 3: Collapse + wire the marker** — one shell element wrapping `<Outlet/>`; add the `/samples`→`/experiments` `<Navigate replace>`; delete the three shell files + the beamtime chip path; point `ExperimentShell` at `TopNav`. Update `FocusPage` `dismiss` per Interfaces (read `from` via `useSearchParams`; `experiment_id` from the loaded sample). Add `?from=series` at the Series-member open call site (grep the Series builder/folio for the member-open `navigate`).

- [ ] **Step 4: Run to verify it passes** — `npx vitest run AppRoutes.test.tsx FocusPage.return.test.tsx` + `git grep -nE "CorpusShell|CorpusTopbar|ExperimentTopNav|beamtime" -- frontend/src` (expect nothing) → PASS.

- [ ] **Step 5: Commit**

```bash
git rm packages/HimalayaUI/frontend/src/print/shell/CorpusShell.tsx packages/HimalayaUI/frontend/src/print/shell/CorpusTopbar.tsx packages/HimalayaUI/frontend/src/print/shell/ExperimentTopNav.tsx
git add packages/HimalayaUI/frontend/src/print/shell/AppRoutes.tsx packages/HimalayaUI/frontend/src/print/shell/ExperimentShell.tsx packages/HimalayaUI/frontend/src/print/pages/FocusPage.tsx packages/HimalayaUI/frontend/test
git commit   # "feat(shell): collapse to one shell; flat Focus/Loupe + from=series return" + trailers
```

### Task 3.3: Contextual bottom dock (ONE `Dock`, page-composed)

> **One component, no surface enum, no shell state (ponytail review):** build a single `Dock` appearance shell that accepts placement **children**. Each page renders `<Dock>` *inside itself* (below its content), composing its own segments — up-link, steppers, the cull `Button`s it already wires, destination `Button`s — directly from that page's local cursor/callbacks. No `SurfaceDock` switch component, no lifting state up through a new context: the dock lives where the data already is. The §3.3 grammar is just the left-to-right order each page passes its children.

**Files:**
- Create: `frontend/src/print/ui/Dock.tsx` (appearance-only shell)
- Modify: `frontend/src/print/pages/SamplesPage.tsx`, `FocusPage.tsx`, `LoupePage.tsx`, and the Series page(s) — each renders its own `<Dock>` composition
- Reference: `ComposeBar.tsx:45-82` (float API to mirror), `CullBar.tsx:51-105` (verb-button grammar), `floatingDock.ts:1-30` (`centerLaneOccupied` lane coordination)
- Test: `frontend/test/Dock.test.tsx` + per-page dock assertions in the existing page tests

**Interfaces:**
- Produces: `<Dock>{children}</Dock>` — a fixed, **light** bar (plate + hairline + soft upward shadow; all appearance in the `ui/` primitive, placement-only `className`), registering with `floatingDock`'s `centerLaneOccupied` so it doesn't collide with transient floats. Each page composes children in §3.3 order using existing primitives: Corpus = `‹ Experiments · Sample↑↓ · Frame‹› · Drop · Keep · Restore · Loupe · Focus` (**no Set-representative**); Loupe = `‹ Corpus · Sample↑↓ · Frame‹› · Drop · Keep · Set representative · Restore · Focus`; Focus = `‹ Corpus · Sample↑↓ · Loupe`; Series = `‹ All series · Sample↑↓ · Focus`. Button variants: Focus = filled accent + frosted ↵; Loupe = neutral pill; Drop/Keep = coloured outlines; Restore/Set-rep = neutral.

- [ ] **Step 1: Write the failing tests** (Dock primitive + per-page composition)

```tsx
// Dock primitive
test('Dock renders children in a fixed light bar', () => {
  render(<Dock><button>seg</button></Dock>)
  expect(screen.getByRole('button', { name: 'seg' })).toBeInTheDocument()
})
// Corpus composition (in SamplesPage test): no Set-rep
test('Corpus dock = up-link, steppers, Drop/Keep/Restore, destinations — no Set-rep', () => {
  renderSamplesPage()
  expect(screen.getByRole('link', { name: /experiments/i })).toBeInTheDocument()
  for (const v of ['Drop', 'Keep', 'Restore']) expect(screen.getByRole('button', { name: v })).toBeInTheDocument()
  expect(screen.queryByRole('button', { name: /representative/i })).toBeNull()
  expect(screen.getByRole('button', { name: /focus/i })).toBeInTheDocument()
})
// Loupe composition (in LoupePage test): Set-rep AND Restore
test('Loupe dock includes Set representative AND Restore', () => {
  renderLoupe()
  expect(screen.getByRole('button', { name: /set representative/i })).toBeInTheDocument()
  expect(screen.getByRole('button', { name: 'Restore' })).toBeInTheDocument()
})
```

- [ ] **Step 2: Run to verify it fails** — `npx vitest run Dock.test.tsx SamplesPage LoupePage` → FAIL.

- [ ] **Step 3: Build `Dock` + compose per page** — `Dock` = the light fixed bar (appearance in `ui/`, registers `centerLaneOccupied`). In each page, render `<Dock>` with its segments built from `Button` + the page's existing cursor/verb callbacks (the same ones M2 wired). The up-link label on Focus reads the `from=series` marker (T3.2) → `‹ Series` else `‹ Corpus`. No new state, no `SurfaceDock`.

- [ ] **Step 4: Run to verify it passes** — `npx vitest run Dock.test.tsx SamplesPage LoupePage FocusPage` → PASS.

- [ ] **Step 5: lint:design + build + commit**

Run: `npm run build` (lint:design + tsc + vite).

```bash
git add packages/HimalayaUI/frontend/src/print/ui/Dock.tsx packages/HimalayaUI/frontend/src/print/pages/SamplesPage.tsx packages/HimalayaUI/frontend/src/print/pages/FocusPage.tsx packages/HimalayaUI/frontend/src/print/pages/LoupePage.tsx packages/HimalayaUI/frontend/test
git commit   # "feat(shell): one Dock, composed per page" + trailers
```

**M3 done:** one chrome (TopNav + dock), the three-chromes problem resolved, `/samples` + beamtime retired.

---

## Milestone 4 — Ingest funnel surfaces (depends on M1 endpoints + M2 sheet + M3 shell)

The two-phase funnel handoff (picker → draft Configuration → Approve creates), phase-1 manifest UI, the ScanFailedPage, and the per-experiment corpus-sheet wiring.

### Task 4.0: Two-phase funnel handoff — picker commits a path, does NOT create

> **The core of spec §2/§6.1, and it's currently violated:** `NewExperimentPage` today (`:61-71`) has a "Scan and create" button that calls `createExperiment` then navigates to `/corpus` — i.e. it creates on submit. The spec wants the picker to **commit a path only** (no DB row), hand off a **client-side draft** to Configuration, and create **only at Approve**. This task builds that handoff; T4.1 builds the Configuration first-run UI that consumes it.

**Files:**
- Create: `frontend/src/lib/draftExperiment.ts` (a tiny Zustand store: `{ path, patterns, setDraft, clear }` — client-side only, no DB)
- Modify: `frontend/src/print/pages/NewExperimentPage.tsx:23-149` ("Review →" commits path → draft → navigate; remove the create-on-submit)
- Modify: `frontend/src/print/shell/AppRoutes.tsx` (add the draft route `/experiments/new/config` → `ConfigurationPage` in first-run mode)
- Test: `frontend/test/NewExperimentPage.test.tsx`

**Interfaces:**
- Produces: a `useDraftExperiment` store holding `{ path, patterns }`. `NewExperimentPage`'s primary button is **"Review →"**: it runs the `validatePath` gate (T1.4), and on `ok` calls `setDraft({ path })` + `navigate("/experiments/new/config")` — **never `createExperiment`**. The draft route is how first-run Configuration is addressed **without an `:id`** (no DB row yet — resolving the spec's addressing gap). **Cancel** anywhere pre-Approve calls `clear()` + aborts any in-flight phase-1 fetch + `navigate("/experiments")`.

- [ ] **Step 1: Write the failing test**

```tsx
test('Review commits the path to a draft and navigates to draft Configuration WITHOUT creating', async () => {
  const create = vi.spyOn(api, 'createExperiment')
  const nav = mockNavigate()
  mockValidatePath({ ok: true, matched: 5, scanned: 5, message: null })
  renderNewExperiment()
  fireEvent.change(screen.getByRole('textbox'), { target: { value: '/data/run42' } })
  fireEvent.click(screen.getByRole('button', { name: /review/i }))
  await waitFor(() => expect(nav).toHaveBeenCalledWith('/experiments/new/config'))
  expect(create).not.toHaveBeenCalled()
  expect(useDraftExperiment.getState().path).toBe('/data/run42')
})
```

- [ ] **Step 2: Run to verify it fails** — `npx vitest run NewExperimentPage.test.tsx` → FAIL (current page calls createExperiment).

- [ ] **Step 3: Implement the draft store + picker handoff + route** — add `draftExperiment.ts`; change `NewExperimentPage`'s submit to validate → `setDraft` → `navigate("/experiments/new/config")` (delete the `createExperiment` call here); add the `/experiments/new/config` route rendering `ConfigurationPage`. Cancel (a button on the picker + on first-run Configuration) → `clear()` + `navigate("/experiments")`.

- [ ] **Step 4: Run to verify it passes** — `npx vitest run NewExperimentPage.test.tsx` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/draftExperiment.ts packages/HimalayaUI/frontend/src/print/pages/NewExperimentPage.tsx packages/HimalayaUI/frontend/src/print/shell/AppRoutes.tsx packages/HimalayaUI/frontend/test/NewExperimentPage.test.tsx
git commit   # "feat(funnel): picker commits a path to a client-side draft (no create)" + trailers
```

### Task 4.1: Configuration first-run mode — phase-1 manifest + Approve creates

**Files:**
- Modify: `frontend/src/api.ts` (add `fetchManifest`; `validatePath`/`suggestPaths` already target the now-live `/api/fs/*` URLs — no change needed there), `frontend/src/print/pages/ConfigurationPage.tsx:15-30`
- Consumes: `useDraftExperiment` (T4.0)
- Reuse: `ProgressBar`, `NoticePill`, `Field`, `Input` primitives
- Test: `frontend/test/ConfigurationPage.test.tsx`

**Interfaces:**
- Produces: `fetchManifest(path, patterns) → { total, matched: {image,metadata,integration}, unmatched: [{file, miss}] }` (GET `/api/fs/manifest`, query params). **First-run mode is discriminated by the route**: `/experiments/new/config` (no `:id`) → read path+patterns from `useDraftExperiment`; `/experiments/:id/config` → later-edit mode (existing behaviour). First-run shows a `ProgressBar` while the manifest fetch is in flight, then matched-by-type counts + the unmatched list, with **Approve disabled until the manifest resolves**; editing a pattern updates the draft + re-fetches (debounced). **Geometry/sources region is hidden first-run** (no PRP parsed yet, §6.5); shown only in later-edit mode. **Approve** calls `createExperiment({ path, patterns })`, then `clear()`s the draft and `navigate(\`/experiments/\${created.id}/corpus\`)`.

- [ ] **Step 1: Write the failing test** (first-run via the draft route; geometry absent)

```tsx
test('Configuration first-run runs the manifest, hides geometry, gates Approve, then creates', async () => {
  useDraftExperiment.setState({ path: '/data/run42', patterns: {} })
  mockFetchManifest({ total: 4, matched: { image: 2, metadata: 1, integration: 0 }, unmatched: [{ file: 's2', miss: 'metadata' }] })
  const create = vi.spyOn(api, 'createExperiment').mockResolvedValue({ id: 9 } as any)
  const nav = mockNavigate()
  renderConfiguration({ route: '/experiments/new/config' })   // first-run = draft route, no :id
  expect(screen.getByRole('button', { name: /approve/i })).toBeDisabled()        // while indexing
  expect(await screen.findByText(/2 image/i)).toBeInTheDocument()
  expect(screen.queryByText(/geometry/i)).toBeNull()                              // hidden first-run
  expect(screen.getByRole('button', { name: /approve/i })).toBeEnabled()
  fireEvent.click(screen.getByRole('button', { name: /approve/i }))
  await waitFor(() => expect(create).toHaveBeenCalledWith(expect.objectContaining({ path: '/data/run42' })))
  await waitFor(() => expect(nav).toHaveBeenCalledWith('/experiments/9/corpus'))
})
```

- [ ] **Step 2: Run to verify it fails** — `npx vitest run ConfigurationPage.test.tsx` → FAIL (no first-run mode / no fetchManifest).

- [ ] **Step 3: Add `fetchManifest` + first-run mode** — add `fetchManifest(path, patterns)` to `api.ts` (GET `/api/fs/manifest` with query params; mirror `suggestPaths`). In `ConfigurationPage`, branch on `useParams().id` presence: absent → first-run (read `useDraftExperiment`), present → later-edit. First-run: `useQuery` the manifest keyed by `(path, patterns)`; `ProgressBar` while `isFetching`; counts + unmatched (`NoticePill` per miss); **don't render the geometry/sources region**; `disabled={isFetching}` on Approve; Approve → `createExperiment({ path, patterns })` → `clear()` + `navigate(\`/experiments/\${created.id}/corpus\`)`. Later-edit keeps the geometry/sources table.

- [ ] **Step 4: Run to verify it passes** — `npx vitest run ConfigurationPage.test.tsx` → PASS.

- [ ] **Step 5: lint:design + commit**

```bash
git add packages/HimalayaUI/frontend/src/api.ts packages/HimalayaUI/frontend/src/print/pages/ConfigurationPage.tsx packages/HimalayaUI/frontend/test/ConfigurationPage.test.tsx
git commit   # "feat(funnel): Configuration first-run manifest + Approve creates" + trailers
```

### Task 4.2: `ScanFailedPage`

**Files:**
- Create: `frontend/src/print/pages/ScanFailedPage.tsx`
- Modify: `frontend/src/print/pages/ExperimentCorpusPage.tsx` (render it when `ingest_status === 'failed'`)
- Reuse: `Menu`, `Field`, `Button`, `NoticePill`
- Test: `frontend/test/ScanFailedPage.test.tsx`

**Interfaces:**
- Consumes: the manifest `unmatched` shape `[{file, miss}]` (Task 1.5/4.1 — no `nearest`).
- Produces: `<ScanFailedPage experimentId={...} />` — **Open Configuration** (primary); a **scrollable** list of all unmatched files grouped by miss type, each showing its **stem + which sidecar type is missing**; an **adaptive pattern test** (one `Field` per affected type, clearing independently) → "Apply all in Configuration"; "Ingest N that parsed" = a real `Button` with a two-stage in-place confirm.

- [ ] **Step 1: Write the failing test**

```tsx
test('ScanFailedPage groups misses, offers per-type test + ingest-N confirm', () => {
  renderScanFailed({ unmatched: [{ file: 'a', miss: 'metadata' }, { file: 'b', miss: 'integration' }], parsedCount: 5 })
  expect(screen.getByRole('button', { name: /open configuration/i })).toBeInTheDocument()
  expect(screen.getByText('a')).toBeInTheDocument()                     // unmatched stem
  expect(screen.getByText(/metadata/i)).toBeInTheDocument()             // miss-type label
  expect(screen.getAllByRole('textbox').length).toBe(2)                 // one field per affected type
  fireEvent.click(screen.getByRole('button', { name: /ingest 5/i }))
  expect(screen.getByText(/confirm/i)).toBeInTheDocument()              // two-stage
})
```

- [ ] **Step 2: Run to verify it fails** — `npx vitest run ScanFailedPage.test.tsx` → FAIL.

- [ ] **Step 3: Build `ScanFailedPage`** — compose the primitives per the Interfaces; scrollable region for the unmatched list; per-affected-type pattern `Field`s; two-stage confirm via local `armed` state (mirror the existing two-stage confirm pattern used elsewhere — grep for `armed` in pages). Wire into `ExperimentCorpusPage`'s state switch.

- [ ] **Step 4: Run to verify it passes** — `npx vitest run ScanFailedPage.test.tsx` → PASS.

- [ ] **Step 5: lint:design + commit**

```bash
git add packages/HimalayaUI/frontend/src/print/pages/ScanFailedPage.tsx packages/HimalayaUI/frontend/src/print/pages/ExperimentCorpusPage.tsx packages/HimalayaUI/frontend/test/ScanFailedPage.test.tsx
git commit   # "feat(funnel): ScanFailedPage" + trailers
```

### Task 4.3: Wire the per-experiment corpus sheet (`corpus-sheet-slot`)

**Files:**
- Modify: `frontend/src/print/pages/ExperimentCorpusPage.tsx:13-64` (fill the `corpus-sheet-slot` stub, `:56-61`)
- Reuse: the (now plain-table) `SheetTable` + the page cursor from M2; `GroupingReviewPage` for the live-ingest/grouping states
- Test: `frontend/e2e/` (Playwright — the corpus state machine end-to-end)

**Interfaces:**
- Consumes: `useExperiment`/`useLoads` (`queries.ts:53-59,119-139`), `ingestInFlight` (`state.ts:39-95`), the M2 sheet + cursor.
- Produces: `ExperimentCorpusPage` renders the §6.2 state machine — `scanning` → `GroupingReviewPage` (live unfold); `failed` → `ScanFailedPage`; `has-flags` → sheet + review banner; `clean` → sheet; `rescanning` → inline progress. The sheet is scoped to the experiment's samples.

- [ ] **Step 1: Write the failing e2e** (extend `frontend/e2e/`; mock the experiment + loads + SSE ingest frames)

```ts
test('experiment corpus shows the scoped sheet when clean', async ({ page }) => {
  await mockExperiment(page, { id: 7, ingest_status: 'complete', samples: 3 })
  await page.goto('/experiments/7/corpus')
  await expect(page.getByRole('table')).toBeVisible()
  await expect(page.getByRole('row')).toHaveCount(4) // header + 3
})
```

- [ ] **Step 2: Run to verify it fails** — `npm run e2e -- corpus` → FAIL (stub renders nothing).

- [ ] **Step 3: Fill the slot** — replace the `corpus-sheet-slot` stub with the state switch; mount the scoped `SheetTable` (samples from `useLoads`/experiment query) for clean/has-flags; `GroupingReviewPage` for scanning; `ScanFailedPage` for failed; inline progress for rescanning (driven by `ingestInFlight`).

- [ ] **Step 4: Run to verify it passes** — `npm run e2e -- corpus` → PASS.

- [ ] **Step 5: Full gate + commit**

Run: `npm run build && npm test && npm run e2e` (frontend) and `make test-parallel` (backend) — all green.

```bash
git add packages/HimalayaUI/frontend/src/print/pages/ExperimentCorpusPage.tsx packages/HimalayaUI/frontend/e2e
git commit   # "feat(funnel): wire per-experiment corpus sheet + state machine" + trailers
```

**M4 done:** the funnel is complete — picker → Configuration (phase-1 manifest) → Approve → live scan/grouping → Corpus, with the scan-failed branch.

---

## Self-Review (spec coverage) — rev 2

- §1 three-chromes → M3 (T3.1–3.2). §3.1 top tier (no ⚙) → T3.1. §3.2/§3.3 dock + grammar → T3.3 (one `Dock`, page-composed). §4 keyboard set → M2 (T2.1 registry + id-migration, T2.4 Corpus, T2.5 Focus/Loupe). §4.3 roving-grid removal → T2.2; `?` normalization → T2.1; per-row Restore → T2.3. §6.1 two-phase funnel (picker commits path → draft Configuration → Approve creates) → **T4.0 + T4.1**; create-on-approval = existing POST + persisted patterns → T1.2. §6.2 live unfold + `on_progress` → T1.1/T1.2/T4.3. §6.4 routing collapse + `/samples` redirect + **flat Focus/Loupe** + `from=series` return → T3.2. §6.5 surfaces (gallery/picker → T4.0; Configuration → T4.1; Scan-failed → T4.2; Corpus assembly → T4.3). §7 scope (sheet interaction, FocusPage keys) → M2. §8 net-new backend → M1.
- **Deferred (not in this plan, per §8):** banner-vs-dock final placement (polish during T4.3); hashing `scan_signature` (optional follow-up); gallery timeline-rail refinements (gallery reused as-is). Focus `P` / Series `A`/`⌘Enter` stay surface-local per §7 (T2.5 note) — not registry keys.
- **Cross-task type consistency:** the registry ids in T2.1 (`prevSample`/`nextSample`/`prevDetail`/`nextDetail`/`openFocus`/`openLoupe`/`toggleSelect`/`extendPrev`/`extendNext`/`selectAll`/`restore`/`helpOverlay`) are the exact strings bound in T2.4/2.5; `matchShortcut` returns a string id (not an object) everywhere; `fetchManifest`'s `{file, miss}` unmatched shape (T4.1) matches `/api/fs/manifest` (T1.5); `on_progress(p,t)` identical in T1.1/T1.2; `useDraftExperiment` `{path, patterns}` is produced in T4.0 and consumed in T4.1.

## Resolved decisions (rev 2)

- **`/api/fs/validate` ships** — the picker's pre-Approve gate (different funnel step than `manifest`).
- **Flat Focus/Loupe** (`/sample/:id`, `/sample/:id/loupe`) — no `LegacySampleRedirect`, no experiment-scoped nesting (a resolver always flashes a loading screen). The chrome resolves the experiment from the loaded sample; the `from=series` marker (T3.2) handles the return target. Spec §6.4 amended to match.
