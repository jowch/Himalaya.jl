# HimalayaUI REST API Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an Oxygen.jl REST API on top of the HimalayaUI SQLite pipeline, exposing all endpoints the frontend plans will need, plus a `himalaya serve` CLI command.

**Architecture:** Single-process Julia server. `server.jl` owns Oxygen.jl setup and the `serve(db; port)` entry point. Route handlers are split by resource (`routes_*.jl`) and all receive the live `SQLite.DB` via a module-level `Ref`. A `X-Username` header is parsed once per mutating request and audited to `user_actions`. Integration tests boot a real Oxygen instance against a temp SQLite DB on a free port and drive it through HTTP.jl.

**Tech Stack:** Oxygen.jl 1.x, HTTP.jl, JSON3, existing HimalayaUI pipeline (from Plan 1).

**This is Plan 2 of 6.** Plan 1 (Foundation) is complete: SQLite schema, CRUD, `.dat` parser, manifest parser, `auto_group`, `persist_analysis!`, `init_experiment!`, `analyze_exposure!`, CLI init/analyze/show. Plan 3 starts the Vite/TypeScript frontend.

**Out of scope for this plan (deferred):**
- `POST /api/experiments` — spec lists it, but in v1 deployment the server is started with an already-initialized DB via `himalaya serve <path>`. Creating experiments over HTTP is unused by any v1 workflow. The endpoint is omitted; `himalaya init` remains the only creation path.
- Frontend static file serving is wired in Task 10 (serves `packages/HimalayaUI/frontend/dist/` if present, empty in this plan).

---

## File Map

| File | Responsibility |
|---|---|
| `packages/HimalayaUI/Project.toml` | + Oxygen, HTTP deps |
| `packages/HimalayaUI/src/HimalayaUI.jl` | Module entry; add new includes |
| `packages/HimalayaUI/src/server.jl` | Oxygen app, `serve(db; port)`, DB `Ref`, health endpoint, router reset for tests |
| `packages/HimalayaUI/src/json.jl` | NamedTuple/Row → JSON conversion helpers; consistent key names |
| `packages/HimalayaUI/src/actions.jl` | `get_username(req)`; `log_action!(db, req; action, entity_type, entity_id, note)` |
| `packages/HimalayaUI/src/routes_users.jl` | `/api/users`, `/api/users/:username/actions` |
| `packages/HimalayaUI/src/routes_experiments.jl` | `/api/experiments/:id`, `PATCH`, `POST /analyze` |
| `packages/HimalayaUI/src/routes_samples.jl` | samples + sample_tags |
| `packages/HimalayaUI/src/routes_exposures.jl` | exposures + exposure_tags + select |
| `packages/HimalayaUI/src/routes_peaks.jl` | peaks list / add manual / delete + stale-marking |
| `packages/HimalayaUI/src/routes_analysis.jl` | indices list + groups (GET, POST member w/ custom-group clone, DELETE member) |
| `packages/HimalayaUI/src/routes_export.jl` | CSV + JSON experiment summary |
| `packages/HimalayaUI/src/cli.jl` | add `cli_serve` subcommand |
| `packages/HimalayaUI/test/test_http.jl` | Shared fixture: free-port picker, `with_test_server(f, db)` helper |
| `packages/HimalayaUI/test/test_routes_*.jl` | One per route file |

---

## Oxygen.jl primer (read before starting)

Oxygen 1.x provides singleton (module-level macros) and multi-instance (`Router`) APIs. This plan uses the **singleton API** for simplicity. Handlers are registered with `@get`, `@post`, `@patch`, `@delete`. Path templates use `{param}` syntax. Body JSON is parsed with `Oxygen.json(req)`. Headers with `HTTP.header(req, "X-Username", "")`.

Key calls used in this plan:
- `Oxygen.serve(; host, port, async, show_banner=false)` — start server
- `Oxygen.terminate()` — stop the async server
- `Oxygen.resetstate()` — clear registered routes (essential between test-server starts in the same process)
- Return a `Dict` / `NamedTuple` / `Vector` and it's JSON-encoded automatically.
- For non-200 status: `HTTP.Response(status, body)` or throw a `HTTP.Exceptions.StatusError`.

**If any call below does not match the installed Oxygen version, verify against Oxygen.jl docs via context7 before improvising.** Prefer fixing the call over altering the architecture.

---

## The DB `Ref` pattern

The Oxygen singleton routes are registered at include-time, before we know which DB they'll serve. We stash the live DB in a module-level `Ref{Union{SQLite.DB, Nothing}}`. `serve(db; ...)` sets it; handlers read it via `current_db()`. Tests set the DB the same way, boot an async server on a free port, run assertions, then call `Oxygen.terminate()` + `Oxygen.resetstate()` between cases.

This is a conscious tradeoff: one DB per process, which matches the v1 deployment model ("single Julia process per server, started with one experiment path").

---

### Task 1: Server scaffold + health endpoint + test harness

**Files:**
- Modify: `packages/HimalayaUI/Project.toml`
- Create: `packages/HimalayaUI/src/server.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl`
- Create: `packages/HimalayaUI/test/test_http.jl`
- Create: `packages/HimalayaUI/test/test_routes_health.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Add Oxygen and HTTP to deps**

```bash
cd packages/HimalayaUI
julia --project=. -e 'using Pkg; Pkg.add(["Oxygen", "HTTP"])'
cd ../..
```

This updates `Project.toml`'s `[deps]` and `[compat]`. Verify both entries are present.

- [ ] **Step 2: Write the failing test**

Create `packages/HimalayaUI/test/test_http.jl`:
```julia
using Test, HTTP, Sockets, SQLite, JSON3
using HimalayaUI

"""
    with_test_server(f, db)

Start an async Oxygen server bound to the given db on a free port, pass
(port, base_url) to f, tear the server down afterward. Safe to call
multiple times per test file.
"""
function with_test_server(f, db::SQLite.DB)
    port = _find_free_port()
    HimalayaUI.start_test_server!(db, port)
    try
        f(port, "http://127.0.0.1:$port")
    finally
        HimalayaUI.stop_test_server!()
    end
end

function _find_free_port()
    server = Sockets.listen(Sockets.IPv4(0), 0)
    port   = Sockets.getsockname(server)[2]
    close(server)
    Int(port)
end
```

Create `packages/HimalayaUI/test/test_routes_health.jl`:
```julia
@testset "GET /api/health" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)

    with_test_server(db) do port, base
        resp = HTTP.get("$base/api/health")
        @test resp.status == 200
        body = JSON3.read(String(resp.body))
        @test body.status == "ok"
    end
end
```

Update `packages/HimalayaUI/test/runtests.jl` — add both includes:
```julia
using Test

@testset "HimalayaUI" begin
    include("test_db.jl")
    include("test_datfile.jl")
    include("test_manifest.jl")
    include("test_pipeline.jl")
    include("test_http.jl")
    include("test_routes_health.jl")
end
```

- [ ] **Step 3: Run to verify failure**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: fails — `HimalayaUI.start_test_server!` not defined.

- [ ] **Step 4: Implement `server.jl`**

Create `packages/HimalayaUI/src/server.jl`:
```julia
using Oxygen
using HTTP
using SQLite
using JSON3

# Module-level DB reference. All handlers access the live DB via current_db().
const _DB_REF = Ref{Union{SQLite.DB, Nothing}}(nothing)

"""
    current_db() -> SQLite.DB

Return the SQLite DB the server is currently bound to. Errors if none is set.
"""
function current_db()
    db = _DB_REF[]
    db === nothing && error("no DB bound; call serve(db; ...) or start_test_server!")
    db
end

"""
    bind_db!(db)

Set the module-level DB. Called by serve() and the test harness.
"""
function bind_db!(db::SQLite.DB)
    _DB_REF[] = db
    nothing
end

"""
    register_routes!()

Register all HTTP handlers on Oxygen's singleton router. Idempotent so long as
Oxygen.resetstate() is called first; callers (serve / start_test_server!) are
responsible for resetting state before calling this.
"""
function register_routes!()
    @get "/api/health" function()
        Dict("status" => "ok")
    end
end

"""
    serve(db; host="127.0.0.1", port=8080)

Start the HimalayaUI HTTP server bound to `db`. Blocks the caller.
"""
function serve(db::SQLite.DB; host::String = "127.0.0.1", port::Int = 8080)
    Oxygen.resetstate()
    bind_db!(db)
    register_routes!()
    Oxygen.serve(; host, port, show_banner = false)
end

"""
    start_test_server!(db, port)

Start the server on `port` in async mode. Blocks until ready.
"""
function start_test_server!(db::SQLite.DB, port::Int)
    Oxygen.resetstate()
    bind_db!(db)
    register_routes!()
    Oxygen.serve(; host = "127.0.0.1", port, async = true, show_banner = false)
    _wait_for_server(port)
end

"""
    stop_test_server!()

Tear down the async server and clear module state. Safe to call repeatedly.
"""
function stop_test_server!()
    Oxygen.terminate()
    Oxygen.resetstate()
    _DB_REF[] = nothing
    nothing
end

function _wait_for_server(port::Int; timeout_s = 5.0)
    deadline = time() + timeout_s
    while time() < deadline
        try
            resp = HTTP.get("http://127.0.0.1:$port/api/health";
                            connect_timeout = 1, readtimeout = 1,
                            retry = false, status_exception = false)
            resp.status == 200 && return
        catch
            # connection refused — keep trying
        end
        sleep(0.05)
    end
    error("test server on port $port did not become ready within $(timeout_s)s")
end
```

- [ ] **Step 5: Wire `server.jl` into the module**

Edit `packages/HimalayaUI/src/HimalayaUI.jl` — add `include("server.jl")` after `include("cli.jl")`:
```julia
module HimalayaUI

include("db.jl")
include("datfile.jl")
include("manifest.jl")
include("pipeline.jl")
include("cli.jl")
include("server.jl")

export main

end
```

- [ ] **Step 6: Run tests to verify pass**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all prior tests + 1 new health test pass.

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/Project.toml \
        packages/HimalayaUI/src/server.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/test/test_http.jl \
        packages/HimalayaUI/test/test_routes_health.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "feat(HimalayaUI): Oxygen.jl server scaffold with health endpoint and test harness"
```

---

### Task 2: JSON serialization helpers

**Files:**
- Create: `packages/HimalayaUI/src/json.jl`
- Create: `packages/HimalayaUI/test/test_json.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl`

Handlers will convert `Tables.rowtable` rows (NamedTuples) to JSON-friendly `Dict`s. The conversion is mostly identity for NamedTuples, but we centralize it so we can:
- Strip internal columns if needed.
- Handle `missing` → `nothing` (JSON3 encodes `nothing` as `null`; `missing` throws).
- Convert `Int8` BOOLEAN columns back to `Bool`.

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/test_json.jl`:
```julia
using Test
using HimalayaUI: row_to_json, rows_to_json

@testset "row_to_json" begin
    # NamedTuple round-trip
    nt = (id = 1, name = "foo", active = 1, notes = missing)
    j  = row_to_json(nt)
    @test j[:id]     == 1
    @test j[:name]   == "foo"
    @test j[:active] == 1
    @test j[:notes]  === nothing     # missing → nothing

    # BOOLEAN coercion: only when column name is explicitly listed
    j2 = row_to_json(nt; bool_keys = (:active,))
    @test j2[:active] === true
end

@testset "rows_to_json" begin
    rows = [(id = 1, q = 0.1), (id = 2, q = 0.2)]
    @test rows_to_json(rows) == [Dict(:id => 1, :q => 0.1),
                                  Dict(:id => 2, :q => 0.2)]
end
```

- [ ] **Step 2: Run to verify failure**

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_json.jl")'
```
Expected: `row_to_json` not defined.

- [ ] **Step 3: Implement `json.jl`**

Create `packages/HimalayaUI/src/json.jl`:
```julia
"""
    row_to_json(row; bool_keys = ())

Convert a NamedTuple row to a Symbol-keyed Dict suitable for JSON3 encoding.
Converts `missing` to `nothing`. Coerces any keys listed in `bool_keys` to Bool.
"""
function row_to_json(row::NamedTuple; bool_keys::Tuple = ())
    out = Dict{Symbol, Any}()
    for k in propertynames(row)
        v = getproperty(row, k)
        if v isa Missing
            out[k] = nothing
        elseif k in bool_keys
            out[k] = v != 0
        else
            out[k] = v
        end
    end
    out
end

"""
    rows_to_json(rows; bool_keys = ())

Convert a collection of rows (typically the output of `Tables.rowtable`).
"""
function rows_to_json(rows; bool_keys::Tuple = ())
    [row_to_json(r; bool_keys) for r in rows]
end
```

- [ ] **Step 4: Wire in**

`packages/HimalayaUI/src/HimalayaUI.jl` — add `include("json.jl")` before `server.jl`:
```julia
include("db.jl")
include("datfile.jl")
include("manifest.jl")
include("pipeline.jl")
include("cli.jl")
include("json.jl")
include("server.jl")
```

`packages/HimalayaUI/test/runtests.jl` — add `include("test_json.jl")` after `test_pipeline.jl`.

- [ ] **Step 5: Run tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/json.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/test/test_json.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "feat(HimalayaUI): JSON serialization helpers"
```

---

### Task 3: Users routes + action logging

**Files:**
- Create: `packages/HimalayaUI/src/actions.jl`
- Create: `packages/HimalayaUI/src/routes_users.jl`
- Create: `packages/HimalayaUI/test/test_routes_users.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl`
- Modify: `packages/HimalayaUI/src/server.jl` (call `register_users_routes!`)
- Modify: `packages/HimalayaUI/test/runtests.jl`

Endpoints:
- `GET  /api/users` → list of `{id, username}`
- `POST /api/users` with `{username}` → `{id, username}`, 201 on create, 200 if exists
- `GET  /api/users/:username/actions` → array of user actions (newest first)

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/test_routes_users.jl`:
```julia
@testset "users routes" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)

    with_test_server(db) do port, base
        # Empty list
        r = HTTP.get("$base/api/users")
        @test r.status == 200
        @test JSON3.read(String(r.body)) == []

        # Create
        r = HTTP.post("$base/api/users";
            body = JSON3.write(Dict(:username => "alice")),
            headers = ["Content-Type" => "application/json"])
        @test r.status == 201
        created = JSON3.read(String(r.body))
        @test created.username == "alice"
        @test created.id == 1

        # Idempotent — second create returns existing
        r = HTTP.post("$base/api/users";
            body = JSON3.write(Dict(:username => "alice")),
            headers = ["Content-Type" => "application/json"])
        @test r.status == 200
        @test JSON3.read(String(r.body)).id == 1

        # List now has alice
        r = HTTP.get("$base/api/users")
        users = JSON3.read(String(r.body))
        @test length(users) == 1
        @test users[1].username == "alice"

        # Empty audit trail
        r = HTTP.get("$base/api/users/alice/actions")
        @test r.status == 200
        @test JSON3.read(String(r.body)) == []

        # 404 for unknown user
        r = HTTP.get("$base/api/users/nobody/actions";
                     status_exception = false)
        @test r.status == 404
    end
end

@testset "log_action!" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    HimalayaUI.get_or_create_user!(db, "alice")

    req = HTTP.Request("POST", "/dummy", ["X-Username" => "alice"])
    HimalayaUI.log_action!(db, req;
        action      = "test",
        entity_type = "sample",
        entity_id   = 42,
        note        = "hello")

    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM user_actions ORDER BY id DESC LIMIT 1"))
    @test length(rows) == 1
    @test rows[1].action      == "test"
    @test rows[1].entity_type == "sample"
    @test rows[1].entity_id   == 42
    @test rows[1].note        == "hello"
end
```

Add `using DBInterface, Tables` to the top of `test_routes_users.jl`.

- [ ] **Step 2: Run to verify failure**

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_routes_users.jl")'
```
Expected: `log_action!` / `get_or_create_user!` not defined.

- [ ] **Step 3: Implement `actions.jl`**

Create `packages/HimalayaUI/src/actions.jl`:
```julia
"""
    get_username(req) -> Union{String, Nothing}

Return the `X-Username` header value if present and non-empty, else nothing.
"""
function get_username(req::HTTP.Request)
    v = HTTP.header(req, "X-Username", "")
    isempty(v) ? nothing : String(v)
end

"""
    get_or_create_user!(db, username) -> user_id::Int

Look up the user id by username; create the row on miss. Returns the id.
"""
function get_or_create_user!(db::SQLite.DB, username::String)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM users WHERE username = ?", [username]))
    isempty(rows) || return Int(rows[1].id)
    res = DBInterface.execute(db,
        "INSERT INTO users (username) VALUES (?)", [username])
    Int(DBInterface.lastrowid(res))
end

"""
    log_action!(db, req; action, entity_type, entity_id, note=nothing)

Record an entry in `user_actions`. Username comes from the `X-Username` header;
missing header ⇒ user_id stored as NULL. Safe to call from any mutating
handler.
"""
function log_action!(db::SQLite.DB, req::HTTP.Request;
        action::String,
        entity_type::String,
        entity_id::Integer,
        note::Union{String, Nothing} = nothing)
    username = get_username(req)
    user_id  = username === nothing ? nothing : get_or_create_user!(db, username)
    DBInterface.execute(db,
        "INSERT INTO user_actions (user_id, action, entity_type, entity_id, note)
         VALUES (?, ?, ?, ?, ?)",
        [user_id, action, entity_type, Int(entity_id), note])
    nothing
end
```

Add `using HTTP, DBInterface, Tables, SQLite` to the top.

- [ ] **Step 4: Implement `routes_users.jl`**

Create `packages/HimalayaUI/src/routes_users.jl`:
```julia
function register_users_routes!()
    @get "/api/users" function()
        rows = Tables.rowtable(DBInterface.execute(current_db(),
            "SELECT id, username FROM users ORDER BY id"))
        rows_to_json(rows)
    end

    @post "/api/users" function(req)
        body     = Oxygen.json(req)
        username = String(body.username)
        isempty(username) && return HTTP.Response(400,
            JSON3.write(Dict(:error => "username required")))

        db    = current_db()
        rows  = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, username FROM users WHERE username = ?", [username]))
        if !isempty(rows)
            return HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(row_to_json(rows[1])))
        end

        uid = get_or_create_user!(db, username)
        HTTP.Response(201, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => uid, :username => username)))
    end

    @get "/api/users/{username}/actions" function(req, username::String)
        db = current_db()
        urows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM users WHERE username = ?", [username]))
        isempty(urows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "user not found")))

        uid  = Int(urows[1].id)
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, timestamp, action, entity_type, entity_id, note
             FROM user_actions WHERE user_id = ? ORDER BY id DESC", [uid]))
        rows_to_json(rows)
    end
end
```

Add `using HTTP, JSON3, DBInterface, Tables, Oxygen` to the top.

- [ ] **Step 5: Wire in**

Edit `packages/HimalayaUI/src/HimalayaUI.jl`:
```julia
include("db.jl")
include("datfile.jl")
include("manifest.jl")
include("pipeline.jl")
include("cli.jl")
include("json.jl")
include("actions.jl")
include("routes_users.jl")
include("server.jl")
```

Edit `packages/HimalayaUI/src/server.jl` — inside `register_routes!()`, call the new function:
```julia
function register_routes!()
    @get "/api/health" function()
        Dict("status" => "ok")
    end
    register_users_routes!()
end
```

Edit `packages/HimalayaUI/test/runtests.jl` — add `include("test_routes_users.jl")`.

- [ ] **Step 6: Run tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all pass.

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/src/actions.jl \
        packages/HimalayaUI/src/routes_users.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/src/server.jl \
        packages/HimalayaUI/test/test_routes_users.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "feat(HimalayaUI): users API + user_actions audit logging"
```

---

### Task 4: Experiments routes

**Files:**
- Create: `packages/HimalayaUI/src/routes_experiments.jl`
- Create: `packages/HimalayaUI/test/test_routes_experiments.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl`, `src/server.jl`, `test/runtests.jl`

Endpoints:
- `GET   /api/experiments/{id}` → experiment row
- `PATCH /api/experiments/{id}` with `{name?, data_dir?, analysis_dir?, manifest_path?}` → updated row
- `POST  /api/experiments/{id}/analyze` → runs full batch analysis; returns `{analyzed: N}`

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/test_routes_experiments.jl`:
```julia
using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "experiments routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))

    db = HimalayaUI.open_db(tmp)
    exp_id = HimalayaUI.init_experiment!(db;
        name = "E1", path = tmp,
        data_dir = joinpath(tmp, "data"),
        analysis_dir = analysis_dir)
    s_id = HimalayaUI.create_sample!(db; experiment_id = exp_id,
        label = "D1", name = "UX1")
    HimalayaUI.create_exposure!(db; sample_id = s_id, filename = "example_tot")

    with_test_server(db) do port, base
        # GET
        r = HTTP.get("$base/api/experiments/$exp_id")
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test body.id == exp_id
        @test body.name == "E1"

        # PATCH name
        r = HTTP.patch("$base/api/experiments/$exp_id";
            body = JSON3.write(Dict(:name => "E1-renamed")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 200
        @test JSON3.read(String(r.body)).name == "E1-renamed"

        # POST analyze
        r = HTTP.post("$base/api/experiments/$exp_id/analyze";
            headers = ["X-Username" => "alice"])
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test body.analyzed == 1

        # Analysis ran: peaks exist for exposure
        peak_count = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM peaks"))[1].c
        @test peak_count > 0

        # 404
        r = HTTP.get("$base/api/experiments/999"; status_exception = false)
        @test r.status == 404
    end
end
```

- [ ] **Step 2: Run to verify failure**

```bash
julia --project=packages/HimalayaUI -e '
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_routes_experiments.jl")'
```
Expected: 404 for every route since none are registered.

- [ ] **Step 3: Implement `routes_experiments.jl`**

Create `packages/HimalayaUI/src/routes_experiments.jl`:
```julia
using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_experiments_routes!()
    @get "/api/experiments/{id}" function(req, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(row_to_json(rows[1])))
    end

    @patch "/api/experiments/{id}" function(req, id::Int)
        db   = current_db()
        body = Oxygen.json(req)

        fields = Symbol[]
        vals   = Any[]
        for k in (:name, :data_dir, :analysis_dir, :manifest_path)
            if haskey(body, k)
                push!(fields, k)
                push!(vals, body[k])
            end
        end
        if isempty(fields)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "no updatable fields provided")))
        end

        sets = join(["$(string(f)) = ?" for f in fields], ", ")
        DBInterface.execute(db,
            "UPDATE experiments SET $sets WHERE id = ?", vcat(vals, [id]))

        log_action!(db, req; action = "update_experiment",
            entity_type = "experiment", entity_id = id)

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM experiments WHERE id = ?", [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(row_to_json(rows[1])))
    end

    @post "/api/experiments/{id}/analyze" function(req, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT analysis_dir FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))
        analysis_dir = rows[1].analysis_dir

        samples   = get_samples(db, id)
        analyzed  = 0
        skipped   = String[]
        for sm in samples
            for ex in get_exposures(db, Int(sm.id))
                try
                    analyze_exposure!(db, Int(ex.id), analysis_dir)
                    analyzed += 1
                catch e
                    push!(skipped, "$(sm.label)/$(ex.filename): $(sprint(showerror, e))")
                end
            end
        end

        log_action!(db, req; action = "analyze",
            entity_type = "experiment", entity_id = id,
            note = "analyzed=$analyzed skipped=$(length(skipped))")

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:analyzed => analyzed, :skipped => skipped)))
    end
end
```

- [ ] **Step 4: Wire in**

Edit `packages/HimalayaUI/src/HimalayaUI.jl` — add `include("routes_experiments.jl")` after `routes_users.jl`.

Edit `packages/HimalayaUI/src/server.jl` — inside `register_routes!`, call `register_experiments_routes!()`.

Edit `packages/HimalayaUI/test/runtests.jl` — add `include("test_routes_experiments.jl")`.

- [ ] **Step 5: Run tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```
Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/routes_experiments.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/src/server.jl \
        packages/HimalayaUI/test/test_routes_experiments.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "feat(HimalayaUI): experiments API (GET, PATCH, POST analyze)"
```

---

### Task 5: Samples + sample_tags routes

**Files:**
- Create: `packages/HimalayaUI/src/routes_samples.jl`
- Create: `packages/HimalayaUI/test/test_routes_samples.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl`, `src/server.jl`, `test/runtests.jl`

Endpoints:
- `GET    /api/experiments/{id}/samples` — array of `{id, label, name, notes, tags: [...]}`
- `PATCH  /api/samples/{id}` with `{name?, notes?}` → updated row
- `POST   /api/samples/{id}/tags` with `{key, value}` → created tag
- `DELETE /api/samples/{id}/tags/{tag_id}` → 204

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/test_routes_samples.jl`:
```julia
using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "samples routes" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/t", data_dir="/t/d",
                                             analysis_dir="/t/a")
    s1 = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1", name="UX1")
    s2 = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D2", name="UX2")

    with_test_server(db) do port, base
        # List
        r = HTTP.get("$base/api/experiments/$exp_id/samples")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test length(list) == 2
        @test list[1].label == "D1"
        @test list[1].tags  == []

        # PATCH
        r = HTTP.patch("$base/api/samples/$s1";
            body = JSON3.write(Dict(:notes => "hello")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 200
        @test JSON3.read(String(r.body)).notes == "hello"

        # POST tag
        r = HTTP.post("$base/api/samples/$s1/tags";
            body = JSON3.write(Dict(:key => "lipid", :value => "DOPC")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 201
        tag = JSON3.read(String(r.body))
        @test tag.key   == "lipid"
        @test tag.value == "DOPC"
        tag_id = tag.id

        # Samples list now shows the tag
        r    = HTTP.get("$base/api/experiments/$exp_id/samples")
        list = JSON3.read(String(r.body))
        @test length(list[1].tags) == 1
        @test list[1].tags[1].key == "lipid"

        # DELETE tag
        r = HTTP.delete("$base/api/samples/$s1/tags/$tag_id";
            headers = ["X-Username" => "alice"])
        @test r.status == 204

        r    = HTTP.get("$base/api/experiments/$exp_id/samples")
        list = JSON3.read(String(r.body))
        @test list[1].tags == []
    end
end
```

- [ ] **Step 2: Run to verify failure** — `Pkg.test` or targeted include; expect 404s.

- [ ] **Step 3: Implement `routes_samples.jl`**

```julia
using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_samples_routes!()
    @get "/api/experiments/{id}/samples" function(req, id::Int)
        db      = current_db()
        samples = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM samples WHERE experiment_id = ? ORDER BY id", [id]))
        out = map(samples) do sm
            tags = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, key, value, source FROM sample_tags
                 WHERE sample_id = ? ORDER BY id", [Int(sm.id)]))
            d = row_to_json(sm)
            d[:tags] = rows_to_json(tags)
            d
        end
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(out))
    end

    @patch "/api/samples/{id}" function(req, id::Int)
        db   = current_db()
        body = Oxygen.json(req)
        fields, vals = Symbol[], Any[]
        for k in (:name, :notes)
            if haskey(body, k)
                push!(fields, k); push!(vals, body[k])
            end
        end
        if isempty(fields)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "no updatable fields provided")))
        end
        sets = join(["$(string(f)) = ?" for f in fields], ", ")
        DBInterface.execute(db,
            "UPDATE samples SET $sets WHERE id = ?", vcat(vals, [id]))

        log_action!(db, req; action = "update_sample",
            entity_type = "sample", entity_id = id)

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM samples WHERE id = ?", [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(row_to_json(rows[1])))
    end

    @post "/api/samples/{id}/tags" function(req, id::Int)
        db    = current_db()
        body  = Oxygen.json(req)
        key   = String(body.key)
        value = String(body.value)
        res   = DBInterface.execute(db,
            "INSERT INTO sample_tags (sample_id, key, value, source)
             VALUES (?, ?, ?, 'manual')",
            [id, key, value])
        tag_id = Int(DBInterface.lastrowid(res))

        log_action!(db, req; action = "add_tag",
            entity_type = "sample", entity_id = id,
            note = "$key=$value")

        HTTP.Response(201, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => tag_id, :sample_id => id,
                             :key => key, :value => value, :source => "manual")))
    end

    @delete "/api/samples/{id}/tags/{tag_id}" function(req, id::Int, tag_id::Int)
        db = current_db()
        DBInterface.execute(db,
            "DELETE FROM sample_tags WHERE id = ? AND sample_id = ?",
            [tag_id, id])
        log_action!(db, req; action = "remove_tag",
            entity_type = "sample", entity_id = id,
            note = "tag_id=$tag_id")
        HTTP.Response(204)
    end
end
```

- [ ] **Step 4: Wire in** — include + register_samples_routes!() in server.jl + runtests.jl include.

- [ ] **Step 5: Run tests** — expect all pass.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/routes_samples.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/src/server.jl \
        packages/HimalayaUI/test/test_routes_samples.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "feat(HimalayaUI): samples and sample_tags API"
```

---

### Task 6: Exposures + exposure_tags + select routes

**Files:**
- Create: `packages/HimalayaUI/src/routes_exposures.jl`
- Create: `packages/HimalayaUI/test/test_routes_exposures.jl`
- Modify: module include, server registration, runtests.

Endpoints:
- `GET    /api/samples/{id}/exposures` — list; each has `{id, filename, kind, selected, tags, sources}`
  - `sources` is non-empty for derived exposures (role + source_exposure_id)
- `PATCH  /api/exposures/{id}/select` → marks `selected=TRUE` on this exposure, clears other exposures of the same sample
- `POST   /api/exposures/{id}/tags` with `{key, value}` → created tag
- `DELETE /api/exposures/{id}/tags/{tag_id}` → 204
- `POST   /api/exposures/{id}/analyze` → `{analyzed: true}`

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/test_routes_exposures.jl`:
```julia
using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "exposures routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))

    db     = HimalayaUI.open_db(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id,
        label="D1", name="UX1")
    e1 = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    e2 = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="other")

    with_test_server(db) do port, base
        # List
        r = HTTP.get("$base/api/samples/$s_id/exposures")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test length(list) == 2
        @test list[1].filename == "example_tot"
        @test list[1].selected === false
        @test list[1].tags     == []
        @test list[1].sources  == []

        # Select e2
        r = HTTP.patch("$base/api/exposures/$e2/select";
            headers = ["X-Username" => "alice"])
        @test r.status == 200

        r = HTTP.get("$base/api/samples/$s_id/exposures")
        list = JSON3.read(String(r.body))
        sel  = findfirst(x -> x.id == e2, list)
        @test list[sel].selected === true
        other = findfirst(x -> x.id == e1, list)
        @test list[other].selected === false

        # Add tag to e1
        r = HTTP.post("$base/api/exposures/$e1/tags";
            body = JSON3.write(Dict(:key=>"flag", :value=>"good")),
            headers = ["Content-Type"=>"application/json", "X-Username"=>"alice"])
        @test r.status == 201
        tag_id = JSON3.read(String(r.body)).id

        r = HTTP.get("$base/api/samples/$s_id/exposures")
        list = JSON3.read(String(r.body))
        other = findfirst(x -> x.id == e1, list)
        @test length(list[other].tags) == 1

        # Delete tag
        r = HTTP.delete("$base/api/exposures/$e1/tags/$tag_id";
            headers = ["X-Username" => "alice"])
        @test r.status == 204

        # Analyze single exposure
        r = HTTP.post("$base/api/exposures/$e1/analyze";
            headers = ["X-Username" => "alice"])
        @test r.status == 200
        @test JSON3.read(String(r.body)).analyzed === true

        # 404 path param for nonexistent exposure
        r = HTTP.post("$base/api/exposures/9999/analyze";
            headers = ["X-Username" => "alice"], status_exception = false)
        @test r.status == 404
    end
end
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Implement `routes_exposures.jl`**

```julia
using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_exposures_routes!()
    @get "/api/samples/{id}/exposures" function(req, id::Int)
        db  = current_db()
        exs = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM exposures WHERE sample_id = ? ORDER BY id", [id]))
        out = map(exs) do e
            tags = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, key, value, source FROM exposure_tags
                 WHERE exposure_id = ? ORDER BY id", [Int(e.id)]))
            srcs = Tables.rowtable(DBInterface.execute(db,
                "SELECT source_exposure_id, role FROM exposure_sources
                 WHERE averaged_exposure_id = ?", [Int(e.id)]))
            d = row_to_json(e; bool_keys = (:selected,))
            d[:tags]    = rows_to_json(tags)
            d[:sources] = rows_to_json(srcs)
            d
        end
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(out))
    end

    @patch "/api/exposures/{id}/select" function(req, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id FROM exposures WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure not found")))
        sample_id = Int(rows[1].sample_id)

        DBInterface.execute(db,
            "UPDATE exposures SET selected = 0 WHERE sample_id = ?", [sample_id])
        DBInterface.execute(db,
            "UPDATE exposures SET selected = 1 WHERE id = ?", [id])

        log_action!(db, req; action = "select_exposure",
            entity_type = "exposure", entity_id = id)

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => id, :selected => true)))
    end

    @post "/api/exposures/{id}/tags" function(req, id::Int)
        db    = current_db()
        body  = Oxygen.json(req)
        key   = String(body.key)
        value = String(body.value)
        res   = DBInterface.execute(db,
            "INSERT INTO exposure_tags (exposure_id, key, value, source)
             VALUES (?, ?, ?, 'manual')", [id, key, value])
        tag_id = Int(DBInterface.lastrowid(res))

        log_action!(db, req; action = "add_tag",
            entity_type = "exposure", entity_id = id, note = "$key=$value")

        HTTP.Response(201, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id=>tag_id, :exposure_id=>id,
                             :key=>key, :value=>value, :source=>"manual")))
    end

    @delete "/api/exposures/{id}/tags/{tag_id}" function(req, id::Int, tag_id::Int)
        db = current_db()
        DBInterface.execute(db,
            "DELETE FROM exposure_tags WHERE id = ? AND exposure_id = ?",
            [tag_id, id])
        log_action!(db, req; action = "remove_tag",
            entity_type = "exposure", entity_id = id, note = "tag_id=$tag_id")
        HTTP.Response(204)
    end

    @post "/api/exposures/{id}/analyze" function(req, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT e.id, x.analysis_dir
             FROM exposures e JOIN samples s ON s.id = e.sample_id
             JOIN experiments x ON x.id = s.experiment_id
             WHERE e.id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure not found")))

        analyze_exposure!(db, id, String(rows[1].analysis_dir))

        log_action!(db, req; action = "analyze",
            entity_type = "exposure", entity_id = id)

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => id, :analyzed => true)))
    end
end
```

- [ ] **Step 4: Wire in.**
- [ ] **Step 5: Run tests.**
- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/routes_exposures.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/src/server.jl \
        packages/HimalayaUI/test/test_routes_exposures.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "feat(HimalayaUI): exposures, exposure_tags, select, and per-exposure analyze API"
```

---

### Task 7: Peaks routes + stale-marking helper

**Files:**
- Create: `packages/HimalayaUI/src/routes_peaks.jl`
- Create: `packages/HimalayaUI/test/test_routes_peaks.jl`
- Modify: module, server, runtests.

Endpoints:
- `GET    /api/exposures/{id}/peaks` — auto + manual peaks
- `POST   /api/exposures/{id}/peaks` with `{q}` → creates a manual peak; marks any `indices` that (via `index_peaks`) rely on a peak within a tolerance of this new q as `status = 'stale'`.
  - Simpler semantics implemented here: **any index on this exposure becomes stale** when a manual peak is added. Deferring precise per-peak invalidation is fine — re-analysis reruns the whole pipeline anyway. Document this in the response body as `{stale_indices: N}`.
- `DELETE /api/peaks/{id}` → 204. Marks all indices on that exposure stale.

- [ ] **Step 1: Write failing test**

Create `packages/HimalayaUI/test/test_routes_peaks.jl`:
```julia
using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "peaks routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    with_test_server(db) do port, base
        # GET peaks
        r = HTTP.get("$base/api/exposures/$e_id/peaks")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test length(list) > 0
        @test all(p -> p.source == "auto", list)

        initial_indices = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM indices WHERE exposure_id = ?", [e_id]))
        @test !isempty(initial_indices)

        # POST manual peak
        r = HTTP.post("$base/api/exposures/$e_id/peaks";
            body = JSON3.write(Dict(:q => 0.5)),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 201
        body = JSON3.read(String(r.body))
        @test body.q == 0.5
        @test body.source == "manual"
        peak_id = body.id
        @test body.stale_indices == length(initial_indices)

        # All indices should now be stale
        stale = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM indices
             WHERE exposure_id = ? AND status = 'stale'", [e_id]))
        @test stale[1].c == length(initial_indices)

        # Peaks list now includes the manual peak
        r = HTTP.get("$base/api/exposures/$e_id/peaks")
        list = JSON3.read(String(r.body))
        manual = filter(p -> p.source == "manual", list)
        @test length(manual) == 1

        # DELETE peak
        r = HTTP.delete("$base/api/peaks/$peak_id";
            headers = ["X-Username" => "alice"])
        @test r.status == 204

        r = HTTP.get("$base/api/exposures/$e_id/peaks")
        list = JSON3.read(String(r.body))
        @test !any(p -> p.id == peak_id, list)
    end
end
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Implement `routes_peaks.jl`**

```julia
using HTTP, JSON3, DBInterface, Tables, Oxygen

"""
    mark_all_indices_stale!(db, exposure_id) -> stale_count

Set `status = 'stale'` on every index for this exposure. Returns the count.
"""
function mark_all_indices_stale!(db::SQLite.DB, exposure_id::Int)
    count_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM indices WHERE exposure_id = ?", [exposure_id]))
    DBInterface.execute(db,
        "UPDATE indices SET status = 'stale' WHERE exposure_id = ?", [exposure_id])
    Int(count_rows[1].c)
end

function register_peaks_routes!()
    @get "/api/exposures/{id}/peaks" function(req, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, exposure_id, q, intensity, prominence, sharpness, source
             FROM peaks WHERE exposure_id = ? ORDER BY q", [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(rows_to_json(rows)))
    end

    @post "/api/exposures/{id}/peaks" function(req, id::Int)
        db   = current_db()
        body = Oxygen.json(req)
        q    = Float64(body.q)

        res = DBInterface.execute(db,
            "INSERT INTO peaks (exposure_id, q, source) VALUES (?, ?, 'manual')",
            [id, q])
        peak_id = Int(DBInterface.lastrowid(res))

        stale = mark_all_indices_stale!(db, id)

        log_action!(db, req; action = "add_peak",
            entity_type = "peak", entity_id = peak_id,
            note = "q=$q")

        HTTP.Response(201, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => peak_id, :exposure_id => id,
                             :q => q, :source => "manual",
                             :stale_indices => stale)))
    end

    @delete "/api/peaks/{id}" function(req, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id FROM peaks WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "peak not found")))
        exposure_id = Int(rows[1].exposure_id)

        DBInterface.execute(db,
            "DELETE FROM index_peaks WHERE peak_id = ?", [id])
        DBInterface.execute(db, "DELETE FROM peaks WHERE id = ?", [id])

        mark_all_indices_stale!(db, exposure_id)

        log_action!(db, req; action = "remove_peak",
            entity_type = "peak", entity_id = id)

        HTTP.Response(204)
    end
end
```

- [ ] **Step 4: Wire in.**
- [ ] **Step 5: Run tests.**
- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/routes_peaks.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/src/server.jl \
        packages/HimalayaUI/test/test_routes_peaks.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "feat(HimalayaUI): peaks API with stale-marking on manual edits"
```

---

### Task 8: Indices + Groups routes

**Files:**
- Create: `packages/HimalayaUI/src/routes_analysis.jl`
- Create: `packages/HimalayaUI/test/test_routes_analysis.jl`
- Modify: module, server, runtests.

Endpoints:
- `GET    /api/exposures/{id}/indices` — all candidate indices
  - Each includes a `peaks` array: `[{peak_id, ratio_position, residual}]`
- `GET    /api/exposures/{id}/groups` — all groups (auto + any custom), each with `{id, kind, active, members: [index_id, ...]}`
- `POST   /api/groups/{id}/members` with `{index_id}` — if `{id}` is the auto group and no custom group exists, create a custom group cloned from auto + add this member + activate custom + deactivate auto. If `{id}` is already the custom group, just add the member. Returns the updated custom group.
- `DELETE /api/groups/{id}/members/{index_id}` — same cloning behavior on first edit.

A helper `ensure_custom_group!(db, exposure_id) -> custom_group_id` encapsulates the clone/promote logic.

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/test_routes_analysis.jl`:
```julia
using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "indices + groups routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    with_test_server(db) do port, base
        # Indices
        r = HTTP.get("$base/api/exposures/$e_id/indices")
        @test r.status == 200
        indices = JSON3.read(String(r.body))
        @test length(indices) >= 1
        @test haskey(indices[1], :peaks)

        # Groups — auto only, active
        r = HTTP.get("$base/api/exposures/$e_id/groups")
        @test r.status == 200
        groups = JSON3.read(String(r.body))
        @test length(groups) == 1
        @test groups[1].kind   == "auto"
        @test groups[1].active === true
        auto_gid  = groups[1].id
        auto_mems = groups[1].members

        # Find a candidate index NOT in the auto group
        extra_candidate = nothing
        for ix in indices
            ix.id in auto_mems && continue
            extra_candidate = ix.id
            break
        end

        if extra_candidate !== nothing
            # POST member → clones auto into custom, adds extra_candidate
            r = HTTP.post("$base/api/groups/$auto_gid/members";
                body = JSON3.write(Dict(:index_id => extra_candidate)),
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice"])
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test body.kind   == "custom"
            @test body.active === true
            @test extra_candidate in body.members

            # Auto group now inactive, custom group active
            r = HTTP.get("$base/api/exposures/$e_id/groups")
            groups = JSON3.read(String(r.body))
            @test length(groups) == 2
            auto_g  = first(filter(g -> g.kind == "auto",   groups))
            cust_g  = first(filter(g -> g.kind == "custom", groups))
            @test auto_g.active === false
            @test cust_g.active === true
            @test extra_candidate in cust_g.members
        end

        # DELETE from auto group should also clone (if custom doesn't exist)
        # Reset by wiping groups and re-running analyze to get a fresh auto group
        DBInterface.execute(db,
            "DELETE FROM index_group_members WHERE group_id IN
             (SELECT id FROM index_groups WHERE exposure_id = ?)", [e_id])
        DBInterface.execute(db,
            "DELETE FROM index_groups WHERE exposure_id = ?", [e_id])
        HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

        r = HTTP.get("$base/api/exposures/$e_id/groups")
        groups = JSON3.read(String(r.body))
        @test length(groups) == 1
        auto_gid    = groups[1].id
        auto_mems   = groups[1].members
        if !isempty(auto_mems)
            removed = first(auto_mems)
            r = HTTP.delete("$base/api/groups/$auto_gid/members/$removed";
                headers = ["X-Username" => "alice"])
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test body.kind == "custom"
            @test !(removed in body.members)
        end
    end
end
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Implement `routes_analysis.jl`**

```julia
using HTTP, JSON3, DBInterface, Tables, Oxygen

"""
    ensure_custom_group!(db, exposure_id) -> (group_id, created)

If a custom group exists for this exposure, return its id and `false`.
Otherwise clone the auto group's members into a new custom group, promote
the custom to active and demote auto, return `(new_id, true)`.

Errors if no auto group exists for the exposure.
"""
function ensure_custom_group!(db::SQLite.DB, exposure_id::Int)
    # Existing custom?
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM index_groups
         WHERE exposure_id = ? AND kind = 'custom'", [exposure_id]))
    isempty(rows) || return (Int(rows[1].id), false)

    # Find auto group
    auto_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM index_groups
         WHERE exposure_id = ? AND kind = 'auto'", [exposure_id]))
    isempty(auto_rows) && error("no auto group for exposure $exposure_id")
    auto_id = Int(auto_rows[1].id)

    # Create custom group (active)
    res = DBInterface.execute(db,
        "INSERT INTO index_groups (exposure_id, kind, active)
         VALUES (?, 'custom', 1)", [exposure_id])
    custom_id = Int(DBInterface.lastrowid(res))

    # Clone auto's members into custom
    DBInterface.execute(db,
        "INSERT INTO index_group_members (group_id, index_id)
         SELECT ?, index_id FROM index_group_members WHERE group_id = ?",
        [custom_id, auto_id])

    # Demote auto
    DBInterface.execute(db,
        "UPDATE index_groups SET active = 0 WHERE id = ?", [auto_id])

    (custom_id, true)
end

function _group_with_members(db::SQLite.DB, group_id::Int)
    g = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM index_groups WHERE id = ?", [group_id]))[1]
    members = Tables.rowtable(DBInterface.execute(db,
        "SELECT index_id FROM index_group_members
         WHERE group_id = ? ORDER BY index_id", [group_id]))
    d = row_to_json(g; bool_keys = (:active,))
    d[:members] = [Int(m.index_id) for m in members]
    d
end

function register_analysis_routes!()
    @get "/api/exposures/{id}/indices" function(req, id::Int)
        db = current_db()
        indices = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM indices WHERE exposure_id = ? ORDER BY score DESC", [id]))
        out = map(indices) do ix
            peaks = Tables.rowtable(DBInterface.execute(db,
                "SELECT peak_id, ratio_position, residual
                 FROM index_peaks WHERE index_id = ? ORDER BY ratio_position",
                [Int(ix.id)]))
            d = row_to_json(ix)
            d[:peaks] = rows_to_json(peaks)
            d
        end
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(out))
    end

    @get "/api/exposures/{id}/groups" function(req, id::Int)
        db     = current_db()
        groups = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM index_groups WHERE exposure_id = ? ORDER BY id", [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write([_group_with_members(db, Int(g.id)) for g in groups]))
    end

    @post "/api/groups/{id}/members" function(req, id::Int)
        db   = current_db()
        body = Oxygen.json(req)
        index_id = Int(body.index_id)

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, kind FROM index_groups WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "group not found")))
        exposure_id = Int(rows[1].exposure_id)

        custom_id, _ = ensure_custom_group!(db, exposure_id)

        DBInterface.execute(db,
            "INSERT OR IGNORE INTO index_group_members (group_id, index_id)
             VALUES (?, ?)", [custom_id, index_id])

        log_action!(db, req; action = "confirm_index",
            entity_type = "index", entity_id = index_id)

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(_group_with_members(db, custom_id)))
    end

    @delete "/api/groups/{id}/members/{index_id}" function(req, id::Int, index_id::Int)
        db = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id FROM index_groups WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "group not found")))
        exposure_id = Int(rows[1].exposure_id)

        custom_id, _ = ensure_custom_group!(db, exposure_id)
        DBInterface.execute(db,
            "DELETE FROM index_group_members
             WHERE group_id = ? AND index_id = ?", [custom_id, index_id])

        log_action!(db, req; action = "exclude_index",
            entity_type = "index", entity_id = index_id)

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(_group_with_members(db, custom_id)))
    end
end
```

- [ ] **Step 4: Wire in.**
- [ ] **Step 5: Run tests.**
- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/routes_analysis.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/src/server.jl \
        packages/HimalayaUI/test/test_routes_analysis.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "feat(HimalayaUI): indices and groups API with custom-group clone-on-edit"
```

---

### Task 9: Export endpoint (JSON + CSV)

**Files:**
- Create: `packages/HimalayaUI/src/routes_export.jl`
- Create: `packages/HimalayaUI/test/test_routes_export.jl`
- Modify: module, server, runtests.

Endpoint:
- `GET /api/experiments/{id}/export?format=json|csv` — summary per sample. Default format: `json`.

**JSON schema** (one object per sample):
```json
{
  "label": "D1", "name": "UX1", "notes": null,
  "tags": [{"key": "lipid", "value": "DOPC"}],
  "exposures": [
    {
      "filename": "example_tot", "selected": true,
      "active_group_kind": "custom",
      "indices": [
        {"phase": "Pn3m", "basis": 0.123, "score": 786.25,
         "r_squared": 0.999, "lattice_d": 12.3}
      ]
    }
  ]
}
```

**CSV columns** (one row per sample, flattening the active group): `sample_label, sample_name, exposure_filename, phases, lattice_ds, r_squareds, scores, tags`. Multi-valued columns are semicolon-joined.

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/test_routes_export.jl`:
```julia
using Test, HTTP, JSON3, SQLite, DBInterface, Tables, CSV

@testset "export routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id,
        label="D1", name="UX1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    with_test_server(db) do port, base
        # JSON
        r = HTTP.get("$base/api/experiments/$exp_id/export?format=json")
        @test r.status == 200
        @test occursin("application/json", HTTP.header(r, "Content-Type"))
        body = JSON3.read(String(r.body))
        @test length(body) == 1
        s = body[1]
        @test s.label == "D1"
        @test length(s.exposures) == 1
        @test s.exposures[1].filename == "example_tot"
        @test length(s.exposures[1].indices) >= 1

        # CSV
        r = HTTP.get("$base/api/experiments/$exp_id/export?format=csv")
        @test r.status == 200
        @test occursin("text/csv", HTTP.header(r, "Content-Type"))
        csv_body = String(r.body)
        @test startswith(csv_body, "sample_label,sample_name,exposure_filename,phases")
        @test occursin("D1,UX1,example_tot", csv_body)

        # Default = json
        r = HTTP.get("$base/api/experiments/$exp_id/export")
        @test occursin("application/json", HTTP.header(r, "Content-Type"))

        # Invalid format
        r = HTTP.get("$base/api/experiments/$exp_id/export?format=xml";
                     status_exception = false)
        @test r.status == 400
    end
end
```

Add `CSV` to `[extras]` / `[targets]` in `packages/HimalayaUI/Project.toml` if not already there (it's a regular dep — already added in Plan 1 Task 1; CSV is imported in manifest.jl). No new dep needed.

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Implement `routes_export.jl`**

```julia
using HTTP, JSON3, DBInterface, Tables, Oxygen, CSV

"""
    build_export(db, experiment_id) -> Vector{Dict}

Compose the per-sample summary used by both JSON and CSV exports.
"""
function build_export(db::SQLite.DB, experiment_id::Int)
    samples = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM samples WHERE experiment_id = ? ORDER BY id",
        [experiment_id]))

    out = Dict[]
    for sm in samples
        sid  = Int(sm.id)
        tags = Tables.rowtable(DBInterface.execute(db,
            "SELECT key, value FROM sample_tags WHERE sample_id = ?", [sid]))

        exposures = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM exposures WHERE sample_id = ? ORDER BY id", [sid]))

        ex_out = Dict[]
        for e in exposures
            eid = Int(e.id)
            active = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, kind FROM index_groups
                 WHERE exposure_id = ? AND active = 1
                 ORDER BY id DESC LIMIT 1", [eid]))
            ag_kind = isempty(active) ? nothing : String(active[1].kind)
            active_id = isempty(active) ? nothing : Int(active[1].id)

            indices = if active_id === nothing
                Tables.rowtable(DBInterface.execute(db,
                    "SELECT phase, basis, score, r_squared, lattice_d
                     FROM indices WHERE exposure_id = ? ORDER BY score DESC",
                    [eid]))
            else
                Tables.rowtable(DBInterface.execute(db,
                    "SELECT i.phase, i.basis, i.score, i.r_squared, i.lattice_d
                     FROM indices i
                     JOIN index_group_members m ON m.index_id = i.id
                     WHERE m.group_id = ? ORDER BY i.score DESC", [active_id]))
            end

            push!(ex_out, Dict(
                :filename          => e.filename,
                :selected          => e.selected != 0,
                :active_group_kind => ag_kind,
                :indices           => rows_to_json(indices),
            ))
        end

        push!(out, Dict(
            :label     => sm.label,
            :name      => sm.name,
            :notes     => sm.notes,
            :tags      => rows_to_json(tags),
            :exposures => ex_out,
        ))
    end
    out
end

function _export_csv(summary)
    header = ["sample_label","sample_name","exposure_filename",
              "phases","lattice_ds","r_squareds","scores","tags"]
    rows = Vector{Vector{Any}}()
    for sm in summary
        tag_str = join(["$(t[:key])=$(t[:value])" for t in sm[:tags]], ";")
        if isempty(sm[:exposures])
            push!(rows, [sm[:label], sm[:name], "", "", "", "", "", tag_str])
            continue
        end
        for ex in sm[:exposures]
            phases = join([String(ix[:phase])        for ix in ex[:indices]], ";")
            lds    = join([string(ix[:lattice_d])    for ix in ex[:indices]], ";")
            r2s    = join([string(ix[:r_squared])    for ix in ex[:indices]], ";")
            scs    = join([string(ix[:score])        for ix in ex[:indices]], ";")
            push!(rows, [sm[:label], sm[:name], ex[:filename],
                         phases, lds, r2s, scs, tag_str])
        end
    end
    io = IOBuffer()
    CSV.write(io, (; (Symbol(h) => [r[i] for r in rows] for (i, h) in enumerate(header))...))
    String(take!(io))
end

function register_export_routes!()
    @get "/api/experiments/{id}/export" function(req, id::Int)
        fmt = try
            q = HTTP.queryparams(HTTP.URI(req.target))
            get(q, "format", "json")
        catch
            "json"
        end

        db      = current_db()
        summary = build_export(db, id)

        if fmt == "json"
            return HTTP.Response(200,
                ["Content-Type" => "application/json"],
                JSON3.write(summary))
        elseif fmt == "csv"
            return HTTP.Response(200,
                ["Content-Type" => "text/csv"],
                _export_csv(summary))
        else
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "format must be json or csv")))
        end
    end
end
```

- [ ] **Step 4: Wire in** — include in HimalayaUI.jl, register in server.jl, add test include.

- [ ] **Step 5: Run tests** — expect pass.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/routes_export.jl \
        packages/HimalayaUI/src/HimalayaUI.jl \
        packages/HimalayaUI/src/server.jl \
        packages/HimalayaUI/test/test_routes_export.jl \
        packages/HimalayaUI/test/runtests.jl
git commit -m "feat(HimalayaUI): experiment export API (JSON + CSV)"
```

---

### Task 10: `himalaya serve` CLI + static file serving

**Files:**
- Modify: `packages/HimalayaUI/src/server.jl` (add `dynamicfiles` for `frontend/dist/`)
- Modify: `packages/HimalayaUI/src/cli.jl` (add `cli_serve` subcommand)
- Create: `packages/HimalayaUI/frontend/dist/.gitkeep` (so the dir exists for static serving)

Serve the frontend `dist/` directory at `/`. If the directory is empty (as in this plan), the API is still accessible at `/api/*`. Plan 3+ will populate `frontend/dist/`.

The CLI subcommand opens the DB from `<experiment_path>/himalaya.db` and invokes `HimalayaUI.serve(db; port)`.

- [ ] **Step 1: Add static file serving**

Edit `packages/HimalayaUI/src/server.jl`. At the top of `register_routes!()` (before the API routes), add static serving:

```julia
function register_routes!()
    # Static frontend assets (served only if the directory exists)
    dist_dir = joinpath(pkgdir(HimalayaUI), "frontend", "dist")
    isdir(dist_dir) && Oxygen.dynamicfiles(dist_dir, "/")

    @get "/api/health" function()
        Dict("status" => "ok")
    end
    register_users_routes!()
    register_experiments_routes!()
    register_samples_routes!()
    register_exposures_routes!()
    register_peaks_routes!()
    register_analysis_routes!()
    register_export_routes!()
end
```

Note: `pkgdir(HimalayaUI)` returns the path to the package root. `HimalayaUI` module reference is already available.

Also ensure the `dist` directory exists with at least one tracked file so git preserves it:

```bash
mkdir -p packages/HimalayaUI/frontend/dist
touch packages/HimalayaUI/frontend/dist/.gitkeep
```

- [ ] **Step 2: Add `cli_serve` to `cli.jl`**

Append to `packages/HimalayaUI/src/cli.jl`:

```julia
function cli_serve(args)
    s = ArgParseSettings(prog = "himalaya serve")
    @add_arg_table! s begin
        "experiment_path"
            required = true
        "--port"
            arg_type = Int
            default  = 8080
        "--host"
            default  = "127.0.0.1"
    end
    p = parse_args(args, s; as_symbols = true)

    db_path = joinpath(p[:experiment_path], "himalaya.db")
    isfile(db_path) || error("no himalaya.db at $db_path — run `himalaya init` first")

    db = open_db(p[:experiment_path])
    println("HimalayaUI serving $(p[:experiment_path]) on http://$(p[:host]):$(p[:port])")
    serve(db; host = p[:host], port = p[:port])
end
```

Edit the `main` dispatch in `cli.jl` to add the `serve` branch:

```julia
function main(args = copy(ARGS))
    isempty(args) && (println("Usage: himalaya <command> [args]"); return)
    cmd = popfirst!(args)

    if cmd == "init"
        cli_init(args)
    elseif cmd == "analyze"
        cli_analyze(args)
    elseif cmd == "show"
        cli_show(args)
    elseif cmd == "serve"
        cli_serve(args)
    else
        println("Unknown command: $cmd. Available: init, analyze, show, serve")
    end
end
```

- [ ] **Step 3: Smoke test**

```bash
# Use the experiment dir seeded in Plan 1 smoke tests
julia --project=packages/HimalayaUI -e '
using HimalayaUI
@async main(["serve", "/tmp/himalaya_test", "--port", "8765"])
sleep(1)

using HTTP, JSON3
r = HTTP.get("http://127.0.0.1:8765/api/health")
println("health: ", JSON3.read(String(r.body)))

r = HTTP.get("http://127.0.0.1:8765/api/experiments/1")
println("exp:    ", JSON3.read(String(r.body)).name)

HimalayaUI.stop_test_server!()'
```

Expected: `health: (status = "ok",)` and `exp: TestRun` (or whatever the Plan 1 smoke test set).

If `/tmp/himalaya_test` has been cleared since Plan 1, reseed it:
```bash
rm -rf /tmp/himalaya_test
mkdir -p /tmp/himalaya_test/analysis/automatic_analysis
cp test/data/example_tot.dat /tmp/himalaya_test/analysis/automatic_analysis/example_tot.dat
julia --project=packages/HimalayaUI -e '
using HimalayaUI; main(["init", "/tmp/himalaya_test", "--name", "TestRun"])'
```

- [ ] **Step 4: Run full suite**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
julia --project=. -e 'using Pkg; Pkg.test()'
```
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/server.jl \
        packages/HimalayaUI/src/cli.jl \
        packages/HimalayaUI/frontend/dist/.gitkeep
git commit -m "feat(HimalayaUI): himalaya serve CLI + static frontend mount"
```

---

## Self-Review

**Spec coverage:**

| Spec §4 endpoint | Task |
|---|---|
| `GET /api/users` | 3 |
| `POST /api/users` | 3 |
| `GET /api/users/:username/actions` | 3 |
| `POST /api/experiments` | Out of scope (see header) |
| `GET /api/experiments/:id` | 4 |
| `PATCH /api/experiments/:id` | 4 |
| `POST /api/experiments/:id/analyze` | 4 |
| `GET /api/experiments/:id/export` | 9 |
| `GET /api/experiments/:id/samples` | 5 |
| `PATCH /api/samples/:id` | 5 |
| `POST /api/samples/:id/tags` | 5 |
| `DELETE /api/samples/:id/tags/:tag_id` | 5 |
| `GET /api/samples/:id/exposures` | 6 |
| `PATCH /api/exposures/:id/select` | 6 |
| `POST /api/exposures/:id/tags` | 6 |
| `DELETE /api/exposures/:id/tags/:tag_id` | 6 |
| `POST /api/exposures/:id/analyze` | 6 |
| `GET /api/exposures/:id/peaks` | 7 |
| `POST /api/exposures/:id/peaks` | 7 |
| `DELETE /api/peaks/:id` | 7 |
| `GET /api/exposures/:id/indices` | 8 |
| `GET /api/exposures/:id/groups` | 8 |
| `POST /api/groups/:id/members` | 8 |
| `DELETE /api/groups/:id/members/:index_id` | 8 |
| `X-Username` → `user_actions` | 3 helper, used in 4–8 |
| `himalaya serve` CLI | 10 |
| Static file serving (for future frontend) | 10 |

**Coverage gap (intentional):** `POST /api/experiments` — v1 deployment uses `himalaya init` + `himalaya serve`, no API-driven creation.

**Naming consistency check:**
- `current_db()`, `bind_db!()`, `start_test_server!()`, `stop_test_server!()` — used consistently.
- `log_action!(db, req; action, entity_type, entity_id, note)` — same signature everywhere.
- `row_to_json` / `rows_to_json` — both defined in Task 2, referenced throughout.
- `register_*_routes!()` — consistent naming.
- `ensure_custom_group!(db, exposure_id) -> (id, created)` — Task 8 only consumer.

**No placeholders verified** — every task has complete code blocks, exact paths, and runnable commands.

**Ready for execution.**
