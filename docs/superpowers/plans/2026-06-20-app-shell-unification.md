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
        isempty(prefix) && return _json(Dict(:suggestions => String[]))
        dir  = isdir(prefix) ? prefix : dirname(prefix)
        base = isdir(prefix) ? "" : basename(prefix)
        isdir(dir) || return _json(Dict(:suggestions => String[]))
        kids = String[]
        for name in readdir(dir; sort = true)
            startswith(name, base) || continue
            full = joinpath(dir, name)
            isdir(full) && push!(kids, full)
            length(kids) >= 20 && break
        end
        _json(Dict(:suggestions => kids))
    end
    # validate + manifest added in Tasks 1.4 / 1.5
end

# Small local helper mirroring the JSON response idiom used across routes_*.jl.
_json(x) = HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(x))
```

Add `include("routes_fs.jl")` in `HimalayaUI.jl` (before `server.jl`) and call `register_fs_routes!()` inside `register_routes!` in `server.jl`.

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
            return _json(Dict(:ok => false, :matched => 0, :scanned => 0,
                              :message => "path does not exist or is not a directory"))
        end
        dup = !isempty(DBInterface.execute(current_db(),
            "SELECT 1 FROM experiments WHERE data_dir = ? LIMIT 1", [path]) |> Tables.rowtable)
        if dup
            return _json(Dict(:ok => false, :matched => 0, :scanned => 0,
                              :message => "an experiment already uses this directory"))
        end
        scanned = count(!startswith("."), readdir(path))   # cheap; rich count is /manifest
        _json(Dict(:ok => true, :matched => scanned, :scanned => scanned, :message => nothing))
    end
```

Add `using DBInterface, Tables` to the file head if not already pulled in transitively (match the other routes files).

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
- Produces: `GET /api/fs/manifest?path=<dir>&image_pattern=&metadata_pattern=&integration_pattern=` → `{ total, matched: {image, metadata, integration}, unmatched: [{file, miss: "image"|"metadata"|"integration", nearest: str|null}] }`. **Patterns come from query params, NOT the DB** (phase-1 is pre-experiment). No PRP parse, no DB write.

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
        isdir(path) || return HTTP.Response(400, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "path is not a directory")))
        pats = (image       = get(q, "image_pattern", "{name}.tif"),
                metadata    = get(q, "metadata_pattern", "{name}.prp"),
                integration = get(q, "integration_pattern", "{name}.dat"))
        files = filter(!startswith("."), readdir(path))
        # stem set per type via the existing matcher (config.jl). _stem_for returns
        # the {name} capture or nothing.
        stems(pat) = Set(filter(!isnothing, [_stem_for(f, pat) for f in files]))
        img, meta, integ = stems(pats.image), stems(pats.metadata), stems(pats.integration)
        unmatched = Dict{String,Any}[]
        for s in img
            for (label, set) in (("metadata", meta), ("integration", integ))
                if !(s in set)
                    push!(unmatched, Dict(:file => s, :miss => label,
                        :nearest => _nearest_file(s, files)))
                end
            end
        end
        _json(Dict(:total => length(files),
                   :matched => Dict(:image => length(img), :metadata => length(meta), :integration => length(integ)),
                   :unmatched => unmatched))
    end
```

Add two small private helpers to `routes_fs.jl` (or reuse `config.jl` equivalents if `_stem_for`/`_nearest_file` already exist — grep first; the dossier cites `_matches_prefix_with_boundary`/`resolve_files` as the matcher family):

```julia
# Extract the {name} capture from a filename given a "{name}.ext" pattern, else nothing.
function _stem_for(filename::AbstractString, pattern::AbstractString)
    occursin("{name}", pattern) || return nothing
    pre, post = split(pattern, "{name}"; limit = 2)
    (startswith(filename, pre) && endswith(filename, post) &&
        length(filename) > length(pre) + length(post)) || return nothing
    filename[nextind(filename, lastindex(pre)):prevind(filename, lastindex(filename) - length(post) + 1)]
end

# Nearest existing filename to a stem by case-insensitive prefix, for the
# "did you mean…" pairing in the Scan-failed UI. nothing if none.
function _nearest_file(stem::AbstractString, files::Vector{<:AbstractString})
    cands = sort(filter(f -> startswith(lowercase(f), lowercase(stem)), files))
    isempty(cands) ? nothing : first(cands)
end
```

> **Implementer note:** before writing `_stem_for`, grep `config.jl` for an existing stem/pattern extractor (`_matches_prefix_with_boundary`, `resolve_files`). If one already extracts the `{name}` capture, call it instead of duplicating. Keep the matcher single-sourced.

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

**Files:**
- Modify: `frontend/src/print/shell/shortcuts.ts:10-59` (registry), `:67-75` (`eventCombo`)
- Test: `frontend/test/shortcuts.test.ts` (extend the existing combo test near the CapsLock-X case, `shortcuts.ts:57-60` references it)

**Interfaces:**
- Produces: registry IDs (one semantic id per physical key — `matchShortcut` is flat first-wins, so no key maps to two ids):
  `samplePrev`(`ArrowUp`), `sampleNext`(`ArrowDown`), `detailPrev`(`ArrowLeft`), `detailNext`(`ArrowRight`), `openFocus`(`Enter`), `openLoupe`(`l`), `drop`(`x`), `keep`(`k`), `setRep`(`r`), `toggleSelect`(`Space`), `rangePrev`(`Shift+ArrowLeft`), `rangeNext`(`Shift+ArrowRight`), `selectAll`(`Mod+a`), `restore`(`Backspace`), `reorderUp`(`Alt+ArrowUp`), `reorderDown`(`Alt+ArrowDown`), `helpOverlay`(`?`), `undo`(`Mod+z`), `redo`(`Mod+Shift+z`). Drops `[` / `]`. `eventCombo` emits the literal token `?` for the `?` key regardless of the Shift bit.

- [ ] **Step 1: Write the failing tests**

```ts
import { eventCombo, SHORTCUTS, matchShortcut } from '../src/print/shell/shortcuts'

const ev = (init: Partial<KeyboardEvent>) => new KeyboardEvent('keydown', init)

test('? normalizes to a stable token regardless of Shift', () => {
  expect(eventCombo(ev({ key: '?', shiftKey: true }))).toBe('?')
})

test('arrows are the sample/detail axis; [ and ] are gone', () => {
  expect(matchShortcut(ev({ key: 'ArrowUp' }))?.id).toBe('samplePrev')
  expect(matchShortcut(ev({ key: 'ArrowLeft' }))?.id).toBe('detailPrev')
  expect(SHORTCUTS.find(s => s.combo === '[')).toBeUndefined()
})

test('each physical key resolves to exactly one id (flat registry)', () => {
  const byCombo = new Map<string, number>()
  for (const s of SHORTCUTS) byCombo.set(s.combo, (byCombo.get(s.combo) ?? 0) + 1)
  for (const [combo, n] of byCombo) expect(`${combo}:${n}`).toBe(`${combo}:1`)
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `cd packages/HimalayaUI/frontend && npx vitest run shortcuts.test.ts`
Expected: FAIL — old `[`/`]` ids present; `?` not special-cased.

- [ ] **Step 3: Implement the re-axis + `?` special-case**

Rewrite the registry entries per the Interfaces list. In `eventCombo`, before composing modifiers, add:

```ts
  // '?' is Shift+/ on US layouts but layout-variable; emit a stable token so the
  // help binding is layout-robust (mirrors the CapsLock-X normalization above).
  if (e.key === '?') return '?'
```

- [ ] **Step 4: Run to verify it passes**

Run: `npx vitest run shortcuts.test.ts`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/shell/shortcuts.ts packages/HimalayaUI/frontend/test/shortcuts.test.ts
git commit   # "feat(keys): re-axis registry to sample/detail; stable ? token" + trailers
```

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

Run: `npx vitest run SheetTable.test.tsx SampleTableRow.test.tsx` then `npx tsc --noEmit -p tsconfig.json 2>&1 | grep -i grid` (expect no dangling references).
Expected: PASS; no leftover `useRovingGrid` imports.

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
  const { rerender } = render(<KeptCell kept={3} dropped={0} onRestore={onRestore} />)
  expect(screen.queryByRole('button', { name: /restore/i })).toBeNull()
  rerender(<KeptCell kept={2} dropped={1} onRestore={onRestore} />)
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

### Task 2.4: Page-level `{sampleIndex, frameIndex}` cursor + rebound handlers (Corpus)

**Files:**
- Modify: `frontend/src/print/pages/SamplesPage.tsx:200-289` (state), `:317-339` (keyboard)
- Test: `frontend/test/SamplesPage.keyboard.test.tsx` (new or extend)

**Interfaces:**
- Consumes: registry ids from Task 2.1; `useShortcuts` (`:23-43`, DECLINE pattern) + `suppressGlobalKeys` (`lib/keys.ts:39-64`).
- Produces: a page-level cursor `{ sampleIndex, frameIndex }` (default `{0,0}`); `↑/↓` step `sampleIndex` (clamped), `←/→` step `frameIndex` within the active sample, `Enter` opens Focus on the active sample, `L` opens Loupe, `Space` toggles the active frame in the selection, `Shift+←/→` extends, `X/K` drop/keep (selection else active frame), `Backspace` restores, `Esc` clears selection then declines. **No `R` on Corpus.** Pointer clicks set the same cursor. The single window handler DECLINEs when `suppressGlobalKeys(e)` or `Alt` is held (so `Alt+arrow` reorder isn't shadowed).

- [ ] **Step 1: Write the failing tests**

```tsx
test('Arrow keys drive the page cursor; Alt is declined', () => {
  renderSamplesPage(/* 3 samples, sample0 has 2 frames */)
  fireEvent.keyDown(window, { key: 'ArrowDown' })
  expect(activeSampleIndex()).toBe(1)
  fireEvent.keyDown(window, { key: 'ArrowUp', altKey: true })   // reorder, not cursor
  expect(activeSampleIndex()).toBe(1)                            // unchanged
})

test('Enter opens Focus for the active sample; R does nothing on Corpus', () => {
  const nav = mockNavigate()
  renderSamplesPage()
  fireEvent.keyDown(window, { key: 'Enter' })
  expect(nav).toHaveBeenCalledWith(expect.stringMatching(/\/sample\//))
  fireEvent.keyDown(window, { key: 'r' })
  expect(setRepresentative).not.toHaveBeenCalled()
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run SamplesPage.keyboard.test.tsx`
Expected: FAIL — no page cursor / Enter still opens loupe.

- [ ] **Step 3: Implement the cursor + handler**

Add `const [cursor, setCursor] = useState({ sampleIndex: 0, frameIndex: 0 })`. Register one window handler via `useShortcuts` that, after `if (suppressGlobalKeys(e) || e.altKey) return false`, switches on `matchShortcut(e)?.id`: `samplePrev/Next` clamp-step `sampleIndex` (reset `frameIndex` to 0); `detailPrev/Next` clamp-step `frameIndex` within `samples[sampleIndex].frames.length`; `openFocus` → `navigate(focusPath(activeSample))`; `openLoupe` → loupe; `toggleSelect`/`rangePrev`/`rangeNext`/`selectAll` drive the existing selection state (`selected`/`anchorRef`); `drop`/`keep`/`restore` call the existing verb mutators on the selection-or-active-frame; `Escape` clears selection then `return false`. Pointer click handlers (row, exposure) call `setCursor`. Do NOT register `setRep` on this page.

- [ ] **Step 4: Run to verify it passes**

Run: `npx vitest run SamplesPage.keyboard.test.tsx`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/pages/SamplesPage.tsx packages/HimalayaUI/frontend/test/SamplesPage.keyboard.test.tsx
git commit   # "feat(corpus): page-level cursor + rebound keyboard" + trailers
```

### Task 2.5: Focus + Loupe re-axis

**Files:**
- Modify: `frontend/src/print/pages/FocusPage.tsx:472-524`, `frontend/src/print/pages/LoupePage.tsx:197-254,338-351`
- Test: `frontend/test/FocusPage.keyboard.test.tsx`, `frontend/test/LoupePage.keyboard.test.tsx`

**Interfaces:**
- Produces: Focus — `↑/↓` step sample (was candidate), `←/→` step candidate (was `↑/↓`), `Alt`-guarded; `stepInList`/`previewIndexId`/Esc-ladder/`addArmed` DECLINE preserved; exposure-stepping stays in the rail (FO-EXPSKIP, `:814-827`), NOT on bare arrows. Loupe — `↑/↓` step sample, `←/→` step frame (was `[`/`]` for sample), `Backspace` restore, `R` set-representative (Loupe keeps `R`), `Alt`-guarded.

- [ ] **Step 1: Write the failing tests**

```tsx
// FocusPage
test('Focus: up/down steps sample, left/right steps candidate', () => {
  renderFocus()
  fireEvent.keyDown(window, { key: 'ArrowRight' }); expect(candidateIndex()).toBe(1)
  fireEvent.keyDown(window, { key: 'ArrowDown' });  expect(sampleNav).toHaveBeenCalled()
})
// LoupePage
test('Loupe: arrows = sample/frame; R sets representative; Backspace restores', () => {
  renderLoupe()
  fireEvent.keyDown(window, { key: 'ArrowRight' }); expect(frameIndex()).toBe(1)
  fireEvent.keyDown(window, { key: 'r' });          expect(setRepresentative).toHaveBeenCalled()
  fireEvent.keyDown(window, { key: 'Backspace' });  expect(restore).toHaveBeenCalled()
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run FocusPage.keyboard.test.tsx LoupePage.keyboard.test.tsx`
Expected: FAIL — old axes (candidate on ↑/↓; sample on `[`/`]`).

- [ ] **Step 3: Swap the axes in both handlers**

In `FocusPage.tsx:472-524`, repoint candidate prev/next from `ArrowUp/Down` to `detailPrev/detailNext` (`←/→`) and add `samplePrev/sampleNext` (`↑/↓`) → the existing sample-navigation. Add `if (e.altKey) return false` at the top. Leave the FO-EXPSKIP exposure control (`:814-827`) untouched. In `LoupePage.tsx`, move sample-step from `[`/`]` to `↑/↓`, keep frame on `←/→`, add `restore`(Backspace); keep `setRep`(`r`). Add the Alt-guard.

- [ ] **Step 4: Run to verify it passes**

Run: `npx vitest run FocusPage.keyboard.test.tsx LoupePage.keyboard.test.tsx`
Expected: PASS.

- [ ] **Step 5: `?` overlay + full frontend gate + commit**

Add the `?` help overlay: a new `frontend/src/print/shell/KbdOverlay.tsx` modal (reuse `ModalShell` + render `KbdLegend` for the live registry; `KbdLegend.tsx:12-34` is already registry-driven), opened by the `helpOverlay` id in a global handler (`useGlobalShortcuts.ts`). Test it renders the registry and closes on `Esc`.

Run: `npm test && npm run build` (lint:design + tsc + vite must pass).

```bash
git add packages/HimalayaUI/frontend/src/print/pages/FocusPage.tsx packages/HimalayaUI/frontend/src/print/pages/LoupePage.tsx packages/HimalayaUI/frontend/src/print/shell/KbdOverlay.tsx packages/HimalayaUI/frontend/src/hooks/useGlobalShortcuts.ts packages/HimalayaUI/frontend/test
git commit   # "feat(keys): Focus/Loupe re-axis + ? overlay" + trailers
```

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

### Task 3.2: Routing collapse — one shell + redirects

**Files:**
- Modify: `frontend/src/print/shell/AppRoutes.tsx:94-161`
- Remove: `frontend/src/print/shell/CorpusShell.tsx`, `frontend/src/print/shell/CorpusTopbar.tsx`, `frontend/src/print/shell/ExperimentTopNav.tsx`
- Modify: `frontend/src/print/shell/ExperimentShell.tsx:33-195` (use `TopNav`; header/tabs move into page content per §3.2; mount dock after `<Outlet/>`)
- Test: `frontend/test/AppRoutes.test.tsx` (extend)

**Interfaces:**
- Produces: a single shell wrapping all routes. Redirects: `/samples` → `/experiments`; `/sample/:id` → `/experiments/:expId/sample/:id` (resolve `expId` from the sample — use the existing sample→experiment query; if not yet loaded, render a resolver that fetches then `<Navigate replace>`). New nested routes: `/experiments/:id/sample/:sampleId` (Focus), `/experiments/:id/sample/:sampleId/loupe` (Loupe). `/series/*` stays top-level.

- [ ] **Step 1: Write the failing tests**

```tsx
test('/samples redirects to /experiments', async () => {
  renderApp({ route: '/samples' })
  await waitFor(() => expect(window.location.pathname).toBe('/experiments'))
})
test('legacy /sample/:id redirects to the experiment-scoped Focus route', async () => {
  renderApp({ route: '/sample/42', seedSampleExperiment: { 42: 7 } })
  await waitFor(() => expect(window.location.pathname).toBe('/experiments/7/sample/42'))
})
```

- [ ] **Step 2: Run to verify it fails** — `npx vitest run AppRoutes.test.tsx` → FAIL.

- [ ] **Step 3: Rewrite the route table** — one shell element wrapping `<Outlet/>`; add the redirect routes (a small `LegacySampleRedirect` component that reads `:id`, fetches the sample's experiment, and `<Navigate replace to=...>`); delete the three removed files; point `ExperimentShell` at `TopNav`. Delete the beamtime chip code path entirely.

- [ ] **Step 4: Run to verify it passes** — `npx vitest run AppRoutes.test.tsx` + `npx tsc --noEmit` (no references to removed files) → PASS.

- [ ] **Step 5: Commit**

```bash
git rm packages/HimalayaUI/frontend/src/print/shell/CorpusShell.tsx packages/HimalayaUI/frontend/src/print/shell/CorpusTopbar.tsx packages/HimalayaUI/frontend/src/print/shell/ExperimentTopNav.tsx
git add packages/HimalayaUI/frontend/src/print/shell/AppRoutes.tsx packages/HimalayaUI/frontend/src/print/shell/ExperimentShell.tsx packages/HimalayaUI/frontend/test/AppRoutes.test.tsx
git commit   # "feat(shell): collapse to one shell + legacy redirects" + trailers
```

### Task 3.3: Contextual bottom dock

**Files:**
- Create: `frontend/src/print/ui/Dock.tsx` (the light dock primitive), `frontend/src/print/shell/SurfaceDock.tsx` (per-surface population)
- Reference: `ComposeBar.tsx:45-82` (float API), `CullBar.tsx:51-105` (verb grammar), `floatingDock.ts:1-30` (`centerLaneOccupied`)
- Test: `frontend/test/SurfaceDock.test.tsx`

**Interfaces:**
- Consumes: the page cursor + verb callbacks from M2; `floatingDock` store for lane coordination.
- Produces: `<Dock>` (light plate + hairline + soft upward shadow; appearance in the primitive). `<SurfaceDock surface="corpus"|"focus"|"loupe"|"series" ... />` renders the §3.3 grammar: `‹ up-link · cursor steppers · cull (Corpus/Loupe) · destinations`. Buttons are presentational, delegating to the page's existing callbacks (dock owns no state).

- [ ] **Step 1: Write the failing test**

```tsx
test('Corpus dock shows up-link, steppers, cull verbs (no Set-rep), destinations', () => {
  render(<SurfaceDock surface="corpus" {...corpusDockProps} />)
  expect(screen.getByRole('link', { name: /experiments/i })).toBeInTheDocument()
  expect(screen.getByRole('button', { name: 'Drop' })).toBeInTheDocument()
  expect(screen.getByRole('button', { name: 'Keep' })).toBeInTheDocument()
  expect(screen.getByRole('button', { name: 'Restore' })).toBeInTheDocument()
  expect(screen.queryByRole('button', { name: /representative/i })).toBeNull()  // Corpus has no R
  expect(screen.getByRole('button', { name: /focus/i })).toBeInTheDocument()
})
test('Loupe dock includes Set representative AND Restore', () => {
  render(<SurfaceDock surface="loupe" {...loupeDockProps} />)
  expect(screen.getByRole('button', { name: /set representative/i })).toBeInTheDocument()
  expect(screen.getByRole('button', { name: 'Restore' })).toBeInTheDocument()
})
```

- [ ] **Step 2: Run to verify it fails** — `npx vitest run SurfaceDock.test.tsx` → FAIL.

- [ ] **Step 3: Build `Dock` + `SurfaceDock`** — `Dock` is a light fixed bar (appearance only, in `ui/`). `SurfaceDock` switches on `surface` to lay out the grammar, reusing `Button` variants (Focus = filled accent + frosted ↵; Loupe = neutral pill; Drop/Keep = colored outlines; Restore/Set-rep = neutral). Wire `centerLaneOccupied` so it doesn't collide with transient floats. Mount `<SurfaceDock>` in the shell after `<Outlet/>`, fed by the active page's cursor/callbacks (lift via a small context or shell-level state set by each page on mount).

- [ ] **Step 4: Run to verify it passes** — `npx vitest run SurfaceDock.test.tsx` → PASS.

- [ ] **Step 5: lint:design + build + commit**

Run: `npm run build` (lint:design + tsc + vite).

```bash
git add packages/HimalayaUI/frontend/src/print/ui/Dock.tsx packages/HimalayaUI/frontend/src/print/shell/SurfaceDock.tsx packages/HimalayaUI/frontend/src/print/shell/ExperimentShell.tsx packages/HimalayaUI/frontend/test/SurfaceDock.test.tsx
git commit   # "feat(shell): contextual bottom dock" + trailers
```

**M3 done:** one chrome (TopNav + dock), the three-chromes problem resolved, `/samples` + beamtime retired.

---

## Milestone 4 — Ingest funnel surfaces (depends on M1 endpoints + M2 sheet + M3 shell)

Phase-1 manifest UI on Configuration, the ScanFailedPage, and the per-experiment corpus-sheet wiring.

### Task 4.1: Wire `validatePath`/`suggestPaths` to the live endpoints + Configuration phase-1 manifest fetch

**Files:**
- Modify: `frontend/src/api.ts:240-256` (confirm the stubs hit `/api/fs/suggest` + `/api/fs/validate`; add `fetchManifest`), `frontend/src/print/pages/ConfigurationPage.tsx:15-30`
- Reuse: `ProgressBar`, `NoticePill`, `Field`, `Input` primitives
- Test: `frontend/test/ConfigurationPage.test.tsx`

**Interfaces:**
- Produces: `fetchManifest(path, patterns) → { total, matched: {image,metadata,integration}, unmatched: [...] }` (GET `/api/fs/manifest`); ConfigurationPage first-run mode shows a spinner while the fetch is in flight, then the matched-by-type counts + an unmatched list, with **Approve disabled until the manifest resolves**; editing a pattern re-fetches (debounced).

- [ ] **Step 1: Write the failing test**

```tsx
test('Configuration runs the phase-1 manifest and gates Approve', async () => {
  mockFetchManifest({ total: 4, matched: { image: 2, metadata: 1, integration: 0 }, unmatched: [{ file: 's2', miss: 'metadata', nearest: null }] })
  renderConfiguration({ firstRun: true })
  expect(screen.getByRole('button', { name: /approve/i })).toBeDisabled()       // while indexing
  expect(await screen.findByText(/2 image/i)).toBeInTheDocument()
  expect(screen.getByRole('button', { name: /approve/i })).toBeEnabled()        // after manifest
})
```

- [ ] **Step 2: Run to verify it fails** — `npx vitest run ConfigurationPage.test.tsx` → FAIL.

- [ ] **Step 3: Add `fetchManifest` + the first-run UI** — add `fetchManifest` to `api.ts` (mirror `suggestPaths` shape). In `ConfigurationPage` first-run mode, `useQuery` the manifest keyed by `(path, patterns)`; render `ProgressBar` while `isFetching`, then counts + unmatched (reuse `NoticePill` for misses); bind `disabled={isFetching}` on Approve; Approve calls `createExperiment({ path, patterns })` (body already carries patterns per `api.ts:218-226`). Geometry table is NOT shown first-run (no PRP parsed yet, §6.5).

- [ ] **Step 4: Run to verify it passes** — `npx vitest run ConfigurationPage.test.tsx` → PASS.

- [ ] **Step 5: lint:design + commit**

```bash
git add packages/HimalayaUI/frontend/src/api.ts packages/HimalayaUI/frontend/src/print/pages/ConfigurationPage.tsx packages/HimalayaUI/frontend/test/ConfigurationPage.test.tsx
git commit   # "feat(funnel): phase-1 manifest on Configuration" + trailers
```

### Task 4.2: `ScanFailedPage`

**Files:**
- Create: `frontend/src/print/pages/ScanFailedPage.tsx`
- Modify: `frontend/src/print/pages/ExperimentCorpusPage.tsx` (render it when `ingest_status === 'failed'`)
- Reuse: `Menu`, `Field`, `Button`, `NoticePill`
- Test: `frontend/test/ScanFailedPage.test.tsx`

**Interfaces:**
- Consumes: the manifest `unmatched` shape (Task 1.5/4.1).
- Produces: `<ScanFailedPage experimentId={...} />` — **Open Configuration** (primary); a **scrollable** list of all unmatched files grouped by miss type, each paired with its `nearest` file; an **adaptive pattern test** (one `Field` per affected type, clearing independently) → "Apply all in Configuration"; "Ingest N that parsed" = a real `Button` with a two-stage in-place confirm.

- [ ] **Step 1: Write the failing test**

```tsx
test('ScanFailedPage groups misses, offers per-type test + ingest-N confirm', () => {
  renderScanFailed({ unmatched: [{ file: 'a', miss: 'metadata', nearest: 'a.PRP' }, { file: 'b', miss: 'integration', nearest: null }], parsedCount: 5 })
  expect(screen.getByRole('button', { name: /open configuration/i })).toBeInTheDocument()
  expect(screen.getByText('a.PRP')).toBeInTheDocument()                 // nearest pairing
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

## Self-Review (spec coverage)

- §1 three-chromes → M3 (T3.1–3.3). §3.1 top tier (no ⚙) → T3.1. §3.2/§3.3 dock + grammar → T3.3. §4 keyboard set → M2 (T2.1, 2.4, 2.5). §4.3 roving-grid removal → T2.2; `?` normalization → T2.1; per-row Restore → T2.3. §6.1 funnel/phase-1 → T4.1; create-on-approval = existing POST + patterns → T1.2/T4.1. §6.2 live unfold + on_progress → T1.1/T1.2/T4.3. §6.4 routing + redirects → T3.2. §6.5 surfaces (gallery/picker reuse; Configuration → T4.1; Scan-failed → T4.2; Corpus assembly → T4.3). §7 scope (sheet interaction, FocusPage keys) → M2. §8 net-new backend → M1.
- **Deferred (not in this plan, per §8):** banner-vs-dock final placement (visual polish during T4.3); hashing `scan_signature` (optional follow-up); experiment-gallery timeline-rail visual refinements (gallery is reused as-is, T-none).
- **Cross-task type consistency:** the registry ids in T2.1 are the exact strings consumed in T2.4/2.5/3.3; `fetchManifest`'s return shape (T4.1) matches `/api/fs/manifest`'s body (T1.5); `on_progress(p,t)` signature is identical in T1.1 and T1.2.

## Open item to confirm with the user before M3

**Experiment-scoped Focus/Loupe paths** (`/experiments/:id/sample/:sampleId` + `/loupe`) are a routing choice the spec states but the codebase currently serves at flat `/sample/:id`. T3.2 implements the redirect; if you'd rather keep Focus/Loupe at the flat path and only scope Corpus/Config/Grouping, say so and T3.2 drops the nesting (the redirect becomes a no-op). Everything else is unaffected.
