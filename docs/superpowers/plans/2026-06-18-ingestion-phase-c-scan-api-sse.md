# Ingestion Redesign — Phase C: Scan/Rescan API + SSE Progress + Scheduler Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire the scan and rescan endpoints, the `broadcast_progress!` helper, and the per-experiment rescan scheduler into the existing server infrastructure so that `POST /api/experiments` kicks off a first scan asynchronously with live SSE progress, `POST /api/experiments/{id}/scan` runs a cheap change-check then additive ingest, `GET /api/experiments/{id}` returns typed geometry + ingest status, `GET /api/experiments/{id}/loads` returns the Load▸Sample▸Exposure roll-up, `PATCH /api/experiments/{id}` writes human geometry overrides (replacing today's 400 stub), and `DELETE /api/experiments/{id}` tears down the experiment and its rescan timer.

**Architecture:** All new routes live in `routes_experiments.jl` and are registered by `register_experiments_routes!()`. The scheduler is a new `Dict{Int,Timer}` module-level constant in `server.jl`, mirroring `GC_TIMER` but keyed per experiment; its Timer callback `@spawn`s the tick body and does both `trylock` and `unlock` on that spawned task (Julia `ReentrantLock` is task-owned). `broadcast_progress!` calls `_try_put!` directly (from `events.jl`) with no `user_actions` row — it rides the existing `"curation"` SSE event name. **Phase B is complete and lives in the same `HimalayaUI` module**, so this plan calls `scan_and_group!(db, experiment_id)` and `cheap_change_check(db, experiment_id)` **directly**: both resolve the experiment's `data_dir` from the row internally (no positional `root_dir`), `scan_and_group!` has no `additive` kwarg (it is idempotent via its dedup INSERT keys, so first-scan == rescan), and no `isdefined` guards are needed (the functions are defined at `include` time).

**Tech Stack:** Julia, SQLite.jl / DBInterface.jl / Tables.jl, Oxygen.jl, HTTP.jl, JSON3.jl, stdlib `Test`. Backend package at `packages/HimalayaUI/`.

**Spec:** `docs/superpowers/specs/2026-06-15-ingestion-redesign-design.md`
- §9.2 — REST endpoints (`POST /api/experiments`, `POST /api/experiments/{id}/scan`, `GET /api/experiments/{id}` extended, `GET /api/experiments/{id}/loads`, geometry `PATCH`, `DELETE`)
- §9.3 — `broadcast_progress!` helper (`_try_put!` direct, `curation` event name, top-level `kind` + `payload.experiment_id`)
- §9.4 — per-experiment rescan scheduler (`Dict{Int,Timer}`, `ReentrantLock` per experiment, `@spawn`, tiered backoff to DB columns `last_scan_tier`/`consecutive_empty_ticks`)

**Source of truth for current code:** anchors were line-verified 2026-06-18 against `routes_experiments.jl`, `server.jl`, and `events.jl`, but line numbers drift — confirm each anchor with a quick `grep`/read before editing.

---

## File Structure

| File | Responsibility | This plan |
|---|---|---|
| `packages/HimalayaUI/src/routes_experiments.jl` | New routes: POST experiments, POST scan, GET extended, GET loads, PATCH geometry, DELETE; geometry serializer | MODIFY |
| `packages/HimalayaUI/src/server.jl` | Rescan timer registry (`RESCAN_TIMERS`, `RESCAN_LOCKS`); `start_rescan_scheduler!`, `stop_rescan_scheduler!`, `stop_all_rescan_timers!`; teardown in `stop_test_server!` | MODIFY |
| `packages/HimalayaUI/src/events.jl` | `broadcast_progress!` helper | MODIFY |
| `packages/HimalayaUI/test/test_ingestion_scan_api.jl` | New standalone test file for this phase | CREATE |
| `packages/HimalayaUI/test/runtests.jl` | Test registry | MODIFY |

**Out of scope (other phases):**
- Phase A — schema migrations, `create_load!`/`create_sample!`/`create_exposure!` signatures, `analyze_exposure!` rewire. This plan assumes those are already landed.
- Phase B — `prp.jl`, `geometry.jl`, `grouping.jl`, `ingest.jl` (`scan_and_group!`, `cheap_change_check`). **Complete (verified 2026-06-19).** This plan calls `scan_and_group!(db, id)` and `cheap_change_check(db, id)` directly but does not define them.
- Phase D — `exposure_moved`, `sample_renamed`, `sample_created`, `sample_split`, `grouping_flag_dismissed` event kinds (there is **no** `sample_merged` — merge fans out as `exposure_moved`); grouping-review structural edits.
- Phase E — frontend routes, `display_name→name` TypeScript sweep, `applyRemoteToCache` `ingest_*` arms, store wiring.

---

## Conventions for every task

- **Run a single test file** from the repo root:
  `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
  The file ends with its own `@testset`s and runs directly; the full suite via `Pkg.test` stays for CI.
- Every test builds a **temp DB** (`mktempdir`) so tests never touch a real database.
- HTTP-layer tests use `with_test_server` (defined in `test/test_http.jl`; include it at the top of the file).
- SSE-layer tests inject a fake subscriber directly into `SSE_SUBSCRIBERS[]` (pattern from `test_sse.jl`).
- **Commit after each task** once its test passes.
- Migration sentinel names follow `"MIGRATION_<THING>"` — no new migrations are needed in Phase C; all schema columns were added in Phase A.

---

## Task 0: Test harness

**Files:**
- Create: `packages/HimalayaUI/test/test_ingestion_scan_api.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Create the standalone test file**

```julia
# packages/HimalayaUI/test/test_ingestion_scan_api.jl
using Test
using HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

include("test_http.jl")   # provides with_test_server

"""
Build a fresh temp DB with one experiment whose data directory exists on disk.
Returns (db, dir, exp_id) where dir is a tempdir that the experiment points at.
"""
function scan_test_db()
    dir = mktempdir()
    db  = HimalayaUI.open_db(joinpath(dir, "himalaya.db"))
    exp_id = HimalayaUI.create_experiment!(db;
        name         = "TestExp",
        path         = dir,
        data_dir     = dir,
        analysis_dir = joinpath(dir, "analysis"))
    (db, dir, exp_id)
end

@testset "ingestion scan API + SSE + scheduler (Phase C)" begin
    # task subtestsets appended below
end
```

> **Note:** `create_experiment!` is the function modified in Phase A to accept geometry kwargs. Confirm the current signature with `grep "function create_experiment!"` in `db.jl` and adjust the call if the `path` param was renamed or made optional.

- [ ] **Step 2: Register the file in runtests.jl**

Add in `packages/HimalayaUI/test/runtests.jl`, after `include("test_ingestion_schema.jl")` (the Phase A file):

```julia
include("test_ingestion_scan_api.jl")
```

- [ ] **Step 3: Run to confirm the harness loads**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: PASS (empty `@testset`; proves imports and `scan_test_db` compile).

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/test/test_ingestion_scan_api.jl packages/HimalayaUI/test/runtests.jl
git commit -m "test: scaffold Phase C scan API test harness"
```

---

## Task 1: `broadcast_progress!`

The spec (§9.3) requires a helper that calls `_try_put!` directly — no `user_actions` row, no `event_id`, no FK contract. It rides the `"curation"` SSE event name so no new `addEventListener` is needed in the frontend. The payload carries `kind` (one of `ingest_started | ingest_progress | ingest_complete | ingest_failed`) plus a non-zero `experiment_id` and count fields.

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl` (add after `broadcast_event!`, ~line 1131)
- Test: `packages/HimalayaUI/test/test_ingestion_scan_api.jl`

- [ ] **Step 1: Write the failing test**

Append to the `@testset "ingestion scan API + SSE + scheduler (Phase C)"`:

```julia
@testset "broadcast_progress! emits curation frame with no user_actions row" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        HimalayaUI.bind_db!(db)

        exp_id = HimalayaUI.create_experiment!(db;
            name = "BP", path = dir, data_dir = dir, analysis_dir = dir)

        pending = Channel{String}(64)
        sub = (pending = pending,)
        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end

        before_count = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions"))[1].c

        HimalayaUI.broadcast_progress!(exp_id; kind = "ingest_started",
            processed = 0, total = 680)

        after_count = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions"))[1].c

        # No durable row written
        @test after_count == before_count

        # Frame is on the channel
        @test isready(pending)
        frame = take!(pending)
        @test occursin("event: curation", frame)

        data_line = first(filter(l -> startswith(l, "data: "), split(frame, '\n')))
        obj = JSON3.read(replace(data_line, r"^data: " => ""))
        # Frame is payload-wrapped (mirrors broadcast_event!): `kind` top-level,
        # experiment/count fields under `payload` so it parses like every other
        # "curation" frame and the frontend reads `payload.experiment_id`.
        @test obj.kind == "ingest_started"
        @test obj.payload.experiment_id == exp_id
        @test obj.payload.processed == 0
        @test obj.payload.total == 680

        lock(HimalayaUI.SSE_LOCK) do
            filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
        end
        close(pending)
        HimalayaUI.SSE_SUBSCRIBERS[] = []
    end
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: FAIL (`broadcast_progress!` undefined).

- [ ] **Step 3: Implement `broadcast_progress!`**

In `events.jl`, add after `broadcast_event!` (~line 1131):

```julia
"""
    broadcast_progress!(experiment_id; kind, processed, total, kwargs...)

Emit a transient ingest-progress SSE frame WITHOUT writing a `user_actions` row.
Calls `_try_put!` directly (per spec §9.3): no `event_id`, no FK contract.

Rides the `"curation"` SSE event name so the existing frontend subscriber
(`App.tsx` `addEventListener("curation", …)`) receives the frame without a new
event type. The frame is PAYLOAD-WRAPPED, mirroring `broadcast_event!`
(events.jl:1107-1120): `kind` is top-level and the experiment/count fields live
under a `payload` sub-object, so the frame parses through the same `curation`
path as every structural event. The `applyRemoteToCache` four `ingest_*` arms
discriminate on the top-level `kind`; the companion `App.tsx` listener updates
`ingestInFlight` from the same frames.

`kind` must be one of: `ingest_started`, `ingest_progress`, `ingest_complete`,
`ingest_failed`. The `experiment_id` is always a non-zero positive integer (the
real experiments.id); the frontend reads `payload.experiment_id`, never
`remote.entity_id`.

Progress frames are best-effort: the SSE channel cap is 64 (server.jl:92); a
680-exposure scan may drop intermediate `ingest_progress` frames. Treat
`ingest_complete` / `ingest_failed` as the authoritative terminal state and
always broadcast those, tolerating drops of intermediate progress ticks.
"""
function broadcast_progress!(experiment_id::Integer;
                              kind::String,
                              processed::Integer = 0,
                              total::Integer     = 0,
                              kwargs...)
    # Build the wire dict ONCE. The events.jl Dates import is selective
    # (`using Dates: now, UTC, format, @dateformat_str`) — bare `Dates` is NOT
    # bound, so use the same unqualified `format(now(UTC), dateformat"…")`
    # expression broadcast_event! uses (events.jl ~line 1115).
    #
    # Frame shape MIRRORS broadcast_event! (events.jl:1107-1120): `:kind` is
    # TOP-LEVEL (the frontend discriminates on it) and the experiment/count
    # fields live under a `:payload` sub-object, so the frame parses exactly like
    # every other `"curation"` frame and the frontend reads `payload.experiment_id`
    # (never a flat top-level field). `:ts` rides top-level, same as broadcast_event!.
    ts      = format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ")
    payload = Dict{Symbol, Any}(
        :experiment_id => Int(experiment_id),
        :processed     => Int(processed),
        :total         => Int(total),
    )
    merge!(payload, Dict{Symbol, Any}(kwargs))
    fields = Dict{Symbol, Any}(
        :kind    => kind,
        :ts      => ts,
        :payload => payload,
    )
    msg   = JSON3.write(fields)
    frame = "event: curation\ndata: $msg\n\n"
    lock(SSE_LOCK) do
        to_drop = []
        for sub in SSE_SUBSCRIBERS[]
            _try_put!(sub.pending, frame) || push!(to_drop, sub)
        end
        for sub in to_drop
            filter!(x -> x !== sub, SSE_SUBSCRIBERS[])
        end
    end
    nothing
end
```

> **Import note (verified):** `events.jl:2` imports Dates **selectively** — `using Dates: now, UTC, format, @dateformat_str`. Bare `Dates` is therefore UNbound; a qualified `Dates.now(...)`/`Dates.format(...)` is an `UndefVarError`. Use the unqualified `format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ")` form — identical to the timestamp expression `broadcast_event!` already uses (events.jl ~line 1115). Do NOT add `using Dates`; the four needed names are already in scope.

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/events.jl packages/HimalayaUI/test/test_ingestion_scan_api.jl
git commit -m "feat(sse): broadcast_progress! — direct _try_put!, no user_actions row"
```

---

## Task 2: `create_load!` writer

Phase A added the `loads` table DDL. This task adds the `create_load!` writer function to `db.jl` so that `scan_and_group!` (Phase B) and the scheduler (this phase) can persist Load rows. (If Phase B already added `create_load!`, skip the implementation step and only add the test assertion.)

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (add after `create_experiment!`)
- Test: `packages/HimalayaUI/test/test_ingestion_scan_api.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "create_load! persists a load row" begin
    db, dir, exp_id = scan_test_db()
    lid = HimalayaUI.create_load!(db;
        experiment_id = exp_id,
        load_index    = 1,
        start_time    = "2026-04-12T08:00:00",
        end_time      = "2026-04-12T08:45:00",
        frame_count   = 48)
    row = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM loads WHERE id = ?", [lid])))
    @test row.experiment_id == exp_id
    @test row.load_index    == 1
    @test row.frame_count   == 48
    @test row.start_time    == "2026-04-12T08:00:00"
    SQLite.close(db)
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: FAIL (`create_load!` undefined).

- [ ] **Step 3: Implement `create_load!`**

In `db.jl`, add after `create_experiment!` (~line 1826):

```julia
function create_load!(db::SQLite.DB;
        experiment_id::Integer,
        load_index::Integer,
        session_id::Union{Integer, Nothing}    = nothing,
        start_time::Union{String, Nothing}     = nothing,
        end_time::Union{String, Nothing}       = nothing,
        frame_count::Integer                   = 0,
        note::Union{String, Nothing}           = nothing)
    res = DBInterface.execute(db, """
        INSERT INTO loads (experiment_id, load_index, session_id, start_time, end_time, frame_count, note)
        VALUES (?,?,?,?,?,?,?)
    """, [experiment_id, load_index, session_id, start_time, end_time, frame_count, note])
    Int(DBInterface.lastrowid(res))
end
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_ingestion_scan_api.jl
git commit -m "feat(db): create_load! writer"
```

---

## Task 3: `GET /api/experiments/{id}` — extended response (typed geometry + ingest status + loads roll-up summary)

The current `GET /api/experiments/{id}` reads only the raw row and overlays `_beamline_from_config`. Phase C replaces the overlay with a serializer that reads the Phase-A typed columns directly: geometry fields + `*_source` tags + `last_scanned_at` + `ingest_status`. It also adds a `stats` sub-object (total loads/samples/exposures) so the frontend experiment header can populate the stat ledger without a separate request.

**Files:**
- Modify: `packages/HimalayaUI/src/routes_experiments.jl` (`_experiment_row_to_json`, `GET /api/experiments/{id}`)
- Test: `packages/HimalayaUI/test/test_ingestion_scan_api.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "GET /api/experiments/{id} includes typed geometry + ingest_status + stats" begin
    db, dir, exp_id = scan_test_db()
    # Seed typed geometry columns directly (Phase B would do this via scan_and_group!)
    DBInterface.execute(db, """
        UPDATE experiments SET
            beam_center_x       = 421.409,
            beam_center_y       = 836.946,
            beam_center_x_source = 'setup',
            beam_center_y_source = 'setup',
            pixel_size_um       = 172.0,
            pixel_size_um_source = 'prp',
            flight_path_m       = 1.8095,
            flight_path_m_source = 'setup',
            energy_kev          = 9.0,
            energy_kev_source   = 'prp',
            q_units             = 'A^-1',
            ingest_status       = 'complete',
            last_scanned_at     = '2026-04-12T10:00:00Z'
        WHERE id = ?
    """, [exp_id])
    # Seed a load + sample + exposure so stats are non-zero
    lid = HimalayaUI.create_load!(db; experiment_id = exp_id, load_index = 1, frame_count = 4)
    sid = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "HA85 (S01P15)",
        load_id = lid, slot_index = 15)
    HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = sid, filename = "f.tif")

    with_test_server(db) do port, base
        r = HTTP.get("$base/api/experiments/$exp_id")
        @test r.status == 200
        body = JSON3.read(String(r.body))

        # Typed geometry
        @test body.beam_center_x       ≈ 421.409
        @test body.beam_center_x_source == "setup"
        @test body.flight_path_m       ≈ 1.8095
        @test body.flight_path_m_source == "setup"
        @test body.energy_kev          ≈ 9.0
        @test body.energy_kev_source   == "prp"
        @test body.pixel_size_um       ≈ 172.0
        @test body.q_units             == "A^-1"
        @test body.ingest_status       == "complete"
        @test body.last_scanned_at     == "2026-04-12T10:00:00Z"

        # Stats roll-up
        @test haskey(body, :stats)
        @test body.stats.loads     == 1
        @test body.stats.samples   == 1
        @test body.stats.exposures == 1
    end
    SQLite.close(db)
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: FAIL (response missing `:stats`, or `:beam_center_x_source` absent, or reading from TOML blob).

- [ ] **Step 3: Extend `_experiment_row_to_json` and the GET route**

In `routes_experiments.jl`, replace `_experiment_row_to_json` and extend the `GET /api/experiments/{id}` handler. Read the current function at ~line 43 first.

```julia
"""
    _experiment_stats(db, exp_id) -> NamedTuple

Cheap roll-up of counts for the shared experiment header stat ledger.
"""
function _experiment_stats(db::SQLite.DB, exp_id::Integer)
    loads = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM loads WHERE experiment_id = ?", [exp_id]))[1].c
    samples = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM samples WHERE experiment_id = ?", [exp_id]))[1].c
    exposures = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM exposures WHERE experiment_id = ?", [exp_id]))[1].c
    (loads = Int(loads), samples = Int(samples), exposures = Int(exposures))
end

"""
    _experiment_row_to_json(row, db) -> Dict

Serialize an experiments row to the wire format. Reads typed geometry columns
from Phase A directly (no TOML overlay). Falls back to `_beamline_from_config`
for legacy rows that still have their geometry only in the TOML `config` blob
(experiments ingested before Phase A).
"""
function _experiment_row_to_json(row::NamedTuple, db::Union{SQLite.DB, Nothing} = nothing)
    d = row_to_json(row)
    # Prefer typed columns (Phase A); fall back to TOML blob for legacy rows.
    has_typed = !isnothing(get(d, :beam_center_x, nothing)) ||
                !isnothing(get(d, :energy_kev, nothing))
    if !has_typed
        bl = _beamline_from_config(get(d, :config, nothing))
        d[:q_units]              = bl.q_units
        d[:beam_center_x]        = bl.beam_center_x
        d[:beam_center_y]        = bl.beam_center_y
        d[:pixel_size_um]        = bl.pixel_size_um
        # energy_kev / flight_path_m are real columns even pre-Phase-A
        # (live create_experiment! writes them); surface their VALUE keys too so
        # the wire shape is identical to the typed path (a legacy row simply
        # reports whatever those columns hold — possibly nothing — never absent).
        d[:energy_kev]           = get(d, :energy_kev, nothing)
        d[:flight_path_m]        = get(d, :flight_path_m, nothing)
        d[:beam_center_x_source] = "default"
        d[:beam_center_y_source] = "default"
        d[:pixel_size_um_source] = "default"
        d[:energy_kev_source]    = "default"
        d[:flight_path_m_source] = "default"
        d[:q_units_source]       = "default"
    end
    # Add stats roll-up when db is supplied (single-row endpoint).
    if db !== nothing
        exp_id = Int(row.id)
        d[:stats] = _experiment_stats(db, exp_id)
    end
    d
end
```

> **Backward compat note:** the old `_experiment_row_to_json(row)` 1-arg form is called by the `GET /api/experiments` list route. Keep that call compiling by making `db` optional (default `nothing`); the list route skips the per-row stats query (acceptable — the home page has its own lighter queries). Update the single-experiment GET to pass `db`:

```julia
@get "/api/experiments/{id}" function(req::HTTP.Request, id::Int)
    db   = current_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM experiments WHERE id = ?", [id]))
    isempty(rows) && return HTTP.Response(404,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "experiment not found")))
    HTTP.Response(200, ["Content-Type" => "application/json"],
        JSON3.write(_experiment_row_to_json(rows[1], db)))
end
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_experiments.jl packages/HimalayaUI/test/test_ingestion_scan_api.jl
git commit -m "feat(api): GET /api/experiments/{id} typed geometry + ingest_status + stats"
```

---

## Task 4: `GET /api/experiments/{id}/loads` — Load▸Sample▸Exposure roll-up

A dedicated endpoint for the grouping-review surface (spec §9.2). Returns an array of loads, each embedding its samples, each sample embedding its exposures. Keyed `queryKeys.loads(id)` on the frontend — distinct from the flat `queryKeys.samples(experimentId)`.

**Files:**
- Modify: `packages/HimalayaUI/src/routes_experiments.jl` (add route)
- Test: `packages/HimalayaUI/test/test_ingestion_scan_api.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "GET /api/experiments/{id}/loads returns Load▸Sample▸Exposure roll-up" begin
    db, dir, exp_id = scan_test_db()
    lid1 = HimalayaUI.create_load!(db; experiment_id = exp_id, load_index = 1, frame_count = 8)
    lid2 = HimalayaUI.create_load!(db; experiment_id = exp_id, load_index = 2, frame_count = 4)
    sid1 = HimalayaUI.create_sample!(db; experiment_id = exp_id,
        name = "HA85 (S01P01)", load_id = lid1, slot_index = 1)
    sid2 = HimalayaUI.create_sample!(db; experiment_id = exp_id,
        name = "HA85 (S02P01)", load_id = lid2, slot_index = 1)
    xid1 = HimalayaUI.create_exposure!(db; experiment_id = exp_id,
        sample_id = sid1, filename = "f1.tif", horizontal_position = 12.4, frame_no = 1)
    xid2 = HimalayaUI.create_exposure!(db; experiment_id = exp_id,
        sample_id = sid1, filename = "f2.tif", horizontal_position = 12.4, frame_no = 2)
    xid3 = HimalayaUI.create_exposure!(db; experiment_id = exp_id,
        sample_id = sid2, filename = "f3.tif", horizontal_position = 39.5, frame_no = 1)

    with_test_server(db) do port, base
        r = HTTP.get("$base/api/experiments/$exp_id/loads")
        @test r.status == 200
        loads = JSON3.read(String(r.body))
        @test length(loads) == 2

        l1 = first(filter(l -> l.load_index == 1, loads))
        @test length(l1.samples) == 1
        @test l1.samples[1].name == "HA85 (S01P01)"
        @test length(l1.samples[1].exposures) == 2
        exps = l1.samples[1].exposures
        filenames = [e.filename for e in exps]
        @test "f1.tif" in filenames && "f2.tif" in filenames

        l2 = first(filter(l -> l.load_index == 2, loads))
        @test length(l2.samples) == 1
        @test length(l2.samples[1].exposures) == 1

        # 404 for unknown experiment
        r2 = HTTP.get("$base/api/experiments/999999/loads"; status_exception = false)
        @test r2.status == 404
    end
    SQLite.close(db)
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: FAIL (404 — route does not exist).

- [ ] **Step 3: Implement the route**

In `routes_experiments.jl`, add inside `register_experiments_routes!()`:

```julia
@get "/api/experiments/{id}/loads" function(req::HTTP.Request, id::Int)
    db   = current_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM experiments WHERE id = ?", [id]))
    isempty(rows) && return HTTP.Response(404,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "experiment not found")))

    loads_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM loads WHERE experiment_id = ? ORDER BY load_index", [id]))

    result = map(loads_rows) do lr
        samples_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM samples WHERE load_id = ? ORDER BY slot_index", [Int(lr.id)]))

        samples = map(samples_rows) do sr
            exposures_rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, filename, prp_path, timestamp, exposure_time,
                        horizontal_position, scan_id, frame_no, status, selected,
                        image_path, content_fingerprint
                   FROM exposures WHERE sample_id = ? ORDER BY frame_no, id",
                [Int(sr.id)]))
            d = row_to_json(sr)
            # exposures.selected is a SQLite 0/1 int; coerce to a JSON bool so this
            # endpoint matches every other exposure serializer (routes_analysis.jl /
            # comparisons.jl ~493 use `row_to_json(e; bool_keys=(:selected,))`).
            d[:exposures] = [row_to_json(e; bool_keys = (:selected,)) for e in exposures_rows]
            d
        end

        d = row_to_json(lr)
        d[:samples] = samples
        d
    end

    HTTP.Response(200, ["Content-Type" => "application/json"],
        JSON3.write(result))
end
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_experiments.jl packages/HimalayaUI/test/test_ingestion_scan_api.jl
git commit -m "feat(api): GET /api/experiments/{id}/loads roll-up"
```

---

## Task 5: `PATCH /api/experiments/{id}` — geometry override (replaces 400 stub)

Today `PATCH /api/experiments/{id}` returns 400 for every call (`routes_experiments.jl:73-85`). Phase C replaces the stub: a partial `ExperimentGeometryPatch` writes the named typed geometry columns and sets their `*_source` to `"user"`. User-set fields must never be refreshed on rescan (spec §9.2 + §4 never-clobber). The route accepts **only the geometry fields**. `name`/`description` editing is deferred to Phase E1 (it needs a `description` column migration that does not exist yet), so `name` stays read-only here (PATCH name → 400, preserving the current derived-name contract). Path fields (`data_dir`, `analysis_dir`, `manifest_path`) also remain read-only (set by create, must sync with the filesystem).

**Files:**
- Modify: `packages/HimalayaUI/src/routes_experiments.jl` (replace `@patch "/api/experiments/{id}"` body)
- Test: `packages/HimalayaUI/test/test_ingestion_scan_api.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "PATCH /api/experiments/{id} writes geometry overrides, marks source=user" begin
    db, dir, exp_id = scan_test_db()

    with_test_server(db) do port, base
        # Patch a geometry field
        r = HTTP.patch("$base/api/experiments/$exp_id";
            body    = JSON3.write(Dict(
                :flight_path_m  => 1.7500,
                :beam_center_x  => 500.0,
                :beam_center_y  => 800.0,
                :pixel_size_um  => 172.0,
                :energy_kev     => 9.0,
                :q_units        => "nm^-1")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            status_exception = false)
        @test r.status == 200

        # Verify persisted values + source override
        row = Tables.rowtable(DBInterface.execute(db, """
            SELECT flight_path_m, flight_path_m_source,
                   beam_center_x, beam_center_x_source,
                   q_units
              FROM experiments WHERE id = ?
        """, [exp_id]))[1]
        @test row.flight_path_m       ≈ 1.7500
        @test row.flight_path_m_source == "user"
        @test row.beam_center_x       ≈ 500.0
        @test row.beam_center_x_source == "user"
        @test row.q_units             == "nm^-1"

        # Partial patch: only flight_path_m. Others untouched.
        r2 = HTTP.patch("$base/api/experiments/$exp_id";
            body    = JSON3.write(Dict(:energy_kev => 12.0)),
            headers = ["Content-Type" => "application/json"],
            status_exception = false)
        @test r2.status == 200
        row2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT energy_kev, energy_kev_source, flight_path_m FROM experiments WHERE id=?",
            [exp_id]))[1]
        @test row2.energy_kev        ≈ 12.0
        @test row2.energy_kev_source == "user"
        @test row2.flight_path_m     ≈ 1.7500  # previous value preserved

        # 404 for unknown experiment
        r3 = HTTP.patch("$base/api/experiments/999999";
            body    = JSON3.write(Dict(:energy_kev => 1.0)),
            headers = ["Content-Type" => "application/json"],
            status_exception = false)
        @test r3.status == 404

        # Read-only fields rejected
        r4 = HTTP.patch("$base/api/experiments/$exp_id";
            body    = JSON3.write(Dict(:data_dir => "/evil")),
            headers = ["Content-Type" => "application/json"],
            status_exception = false)
        @test r4.status == 400

        # name stays read-only here — rename lands in Phase E1 (derived-name
        # contract preserved until then).
        r5 = HTTP.patch("$base/api/experiments/$exp_id";
            body    = JSON3.write(Dict(:name => "Renamed")),
            headers = ["Content-Type" => "application/json"],
            status_exception = false)
        @test r5.status == 400
    end
    SQLite.close(db)
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: FAIL (PATCH still returns 400 for everything).

- [ ] **Step 3: Replace the PATCH stub**

First, add the two field-list `const`s at **module top level** — NOT inside
`register_experiments_routes!()`. A `const` in function/local scope is a Julia
syntax error, and these lists are module-wide policy. Place them alongside
`_beamline_from_config` / `_experiment_row_to_json` (top of `routes_experiments.jl`,
outside any function):

```julia
# Mutable geometry fields. Each writable field has a companion *_source column
# (set to "user" on override; never refreshed by rescan). name/description and
# the path fields (data_dir, analysis_dir, manifest_path) remain read-only here —
# name/description editing lands in Phase E1 (needs a description-column migration);
# path fields are set at create time and must stay in sync with the filesystem.
const _GEOMETRY_PATCH_FIELDS = [
    "flight_path_m", "beam_center_x", "beam_center_y",
    "pixel_size_um", "energy_kev", "q_units",
]
const _READONLY_FIELDS = ["data_dir", "analysis_dir", "manifest_path", "path",
                           "id", "created_at", "name", "description"]
```

Then, INSIDE `register_experiments_routes!()`, replace the entire
`@patch "/api/experiments/{id}"` block (~lines 73–85) with just the route body:

```julia
@patch "/api/experiments/{id}" function(req::HTTP.Request, id::Int)
    db   = current_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM experiments WHERE id = ?", [id]))
    isempty(rows) && return HTTP.Response(404,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "experiment not found")))

    body = json(req)

    # Reject any attempt to write read-only path/id fields.
    for k in _READONLY_FIELDS
        (haskey(body, Symbol(k)) || haskey(body, k)) &&
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "$k is read-only; change it via create/scan")))
    end

    # Build SET clauses for geometry fields present in the body.
    set_clauses = String[]
    params      = Any[]
    for field in _GEOMETRY_PATCH_FIELDS
        val = get(body, Symbol(field), get(body, field, nothing))
        val === nothing && continue
        push!(set_clauses, "$field = ?")
        push!(params, val)
        push!(set_clauses, "$(field)_source = 'user'")
    end

    isempty(set_clauses) && return HTTP.Response(200,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:id => id, :updated => false)))

    push!(params, id)
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            DBInterface.execute(db,
                "UPDATE experiments SET $(join(set_clauses, ", ")) WHERE id = ?",
                params)
        end
    end

    updated_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM experiments WHERE id = ?", [id]))
    HTTP.Response(200, ["Content-Type" => "application/json"],
        JSON3.write(_experiment_row_to_json(updated_rows[1], db)))
end
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_experiments.jl packages/HimalayaUI/test/test_ingestion_scan_api.jl
git commit -m "feat(api): PATCH /api/experiments/{id} geometry overrides (replaces 400 stub)"
```

---

## Task 6: `DELETE /api/experiments/{id}` with timer teardown

No experiment-delete route exists today (spec §9.4: "includes building that delete path with timer teardown"). **There is NO `ON DELETE CASCADE` to lean on:** the live schema (db.jl:32-152, confirmed 2026-06-18) keys a deep tree of structural tables off `samples`/`exposures`/`indices`, and almost none of those FKs cascade — `samples.experiment_id` (db.jl:34), `exposures.sample_id` (db.jl:42), and the exposure-keyed tables (`exposure_sources`, `exposure_tags`, `indices`, `auto_peaks`, `peak_curations`, `index_groups`, `assignments`, `assignment_members`, plus the indices-keyed `index_peaks`/`index_group_members`) all delete-block. Corrected Plan A does not add cascades (SQLite cannot ALTER one onto an existing column). The delete therefore disables FK enforcement connection-wide for the transaction (the codebase's own migration teardown idiom, db.jl:1620-1647) and deletes the whole tree explicitly in FK order. The rescan timer for the experiment (Tasks 7–8) must be stopped before the DB delete so the timer callback cannot fire against a non-existent experiment.

**Files:**
- Modify: `packages/HimalayaUI/src/routes_experiments.jl` (add `DELETE` route)
- Test: `packages/HimalayaUI/test/test_ingestion_scan_api.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "DELETE /api/experiments/{id} removes experiment + cascade" begin
    db, dir, exp_id = scan_test_db()
    lid = HimalayaUI.create_load!(db; experiment_id = exp_id, load_index = 1, frame_count = 2)
    sid = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "A",
        load_id = lid, slot_index = 1)
    HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = sid, filename = "f.tif")

    with_test_server(db) do port, base
        r = HTTP.delete("$base/api/experiments/$exp_id",
            headers = ["X-Username" => "alice"],
            status_exception = false)
        @test r.status == 200

        # Experiment is gone
        r2 = HTTP.get("$base/api/experiments/$exp_id"; status_exception = false)
        @test r2.status == 404

        # Cascades: loads, samples, exposures removed
        cnt_loads = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM loads WHERE experiment_id = ?", [exp_id]))[1].c
        cnt_samples = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM samples WHERE experiment_id = ?", [exp_id]))[1].c
        @test cnt_loads   == 0
        @test cnt_samples == 0

        # 404 on repeat delete
        r3 = HTTP.delete("$base/api/experiments/$exp_id",
            headers = ["X-Username" => "alice"],
            status_exception = false)
        @test r3.status == 404
    end
    SQLite.close(db)
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: FAIL (405 Method Not Allowed — route does not exist).

- [ ] **Step 3: Implement the DELETE route**

In `routes_experiments.jl`, add inside `register_experiments_routes!()`:

```julia
@delete "/api/experiments/{id}" function(req::HTTP.Request, id::Int)
    db   = current_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM experiments WHERE id = ?", [id]))
    isempty(rows) && return HTTP.Response(404,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "experiment not found")))

    # Stop the rescan timer BEFORE the DB delete so the timer callback cannot
    # fire against a non-existent row. stop_rescan_scheduler! is defined in server.jl
    # and is a no-op when no timer is running for this id.
    stop_rescan_scheduler!(id)

    # Live schema (db.jl:32-152, confirmed 2026-06-18) keys a DEEP tree of
    # structural tables off exposures/samples/indices, and almost none of those
    # FKs declare `ON DELETE CASCADE`: samples.experiment_id (db.jl:34),
    # exposures.sample_id (db.jl:42), and the exposure-keyed tables
    # exposure_sources, exposure_tags, indices, auto_peaks, peak_curations,
    # index_groups, assignments, assignment_members (db.jl:60-152), plus the
    # indices-keyed index_peaks / index_group_members. Corrected Plan A does NOT
    # add cascades (SQLite cannot ALTER a cascade onto an existing column).
    # FK enforcement is ON at the connection level (open_db, db.jl:1907), so a
    # bare `DELETE FROM experiments` would FK-fail and 500.
    #
    # Enumerating ~16 child deletes in FK order is brittle and easy to leave
    # incomplete as the schema grows. Instead follow the codebase's own teardown
    # idiom for cross-FK structural surgery (db.jl:1620-1647 and the other
    # migrations): toggle `PRAGMA foreign_keys = OFF` at the CONNECTION level
    # OUTSIDE the transaction (it is a documented no-op mid-transaction), do the
    # parent + cascade deletes, then restore `ON` in a `finally`. Delete in
    # FK-child→parent order anyway so the row set is internally consistent.
    #
    # Concurrency note: this disables FK checks connection-wide for the duration
    # of the delete on the shared singleton connection (parallel=true). Acceptable
    # because (a) experiment delete is a rare admin action, (b) the whole delete
    # is serialized under `_DB_WRITE_LOCK`, and (c) it mirrors the established
    # migration precedent. Do NOT issue `PRAGMA foreign_keys` INSIDE the
    # transaction — SQLite silently ignores it there.
    lock(_DB_WRITE_LOCK) do
        DBInterface.execute(db, "PRAGMA foreign_keys = OFF")
        try
            SQLite.transaction(db) do
                # exposure-keyed structural tables (children of exposures/indices)
                DBInterface.execute(db, """
                    DELETE FROM index_peaks WHERE index_id IN
                      (SELECT i.id FROM indices i
                         JOIN exposures e ON e.id = i.exposure_id
                        WHERE e.experiment_id = ?)""", [id])
                DBInterface.execute(db, """
                    DELETE FROM index_group_members WHERE index_id IN
                      (SELECT i.id FROM indices i
                         JOIN exposures e ON e.id = i.exposure_id
                        WHERE e.experiment_id = ?)""", [id])
                for tbl in ("assignment_members", "assignments", "index_groups",
                            "indices", "auto_peaks", "peak_curations",
                            "exposure_sources", "exposure_tags")
                    # exposure_sources keys two exposure columns; clean both.
                    if tbl == "exposure_sources"
                        DBInterface.execute(db, """
                            DELETE FROM exposure_sources WHERE averaged_exposure_id IN
                              (SELECT id FROM exposures WHERE experiment_id = ?)
                               OR source_exposure_id IN
                              (SELECT id FROM exposures WHERE experiment_id = ?)""",
                            [id, id])
                    else
                        DBInterface.execute(db, """
                            DELETE FROM $tbl WHERE exposure_id IN
                              (SELECT id FROM exposures WHERE experiment_id = ?)""",
                            [id])
                    end
                end
                # sample-keyed tables, then the core rows.
                DBInterface.execute(db, """
                    DELETE FROM sample_tags WHERE sample_id IN
                      (SELECT id FROM samples WHERE experiment_id = ?)""", [id])
                DBInterface.execute(db, """
                    DELETE FROM sample_messages WHERE sample_id IN
                      (SELECT id FROM samples WHERE experiment_id = ?)""", [id])
                DBInterface.execute(db,
                    "DELETE FROM exposures WHERE experiment_id = ?", [id])
                DBInterface.execute(db,
                    "DELETE FROM samples WHERE experiment_id = ?", [id])
                DBInterface.execute(db,
                    "DELETE FROM loads WHERE experiment_id = ?", [id])
                DBInterface.execute(db,
                    "DELETE FROM experiments WHERE id = ?", [id])
            end
        finally
            DBInterface.execute(db, "PRAGMA foreign_keys = ON")
        end
    end

    log_action!(db, req; action = "experiment_deleted",
        entity_type = "experiment", entity_id = id)

    HTTP.Response(200, ["Content-Type" => "application/json"],
        JSON3.write(Dict(:id => id, :deleted => true)))
end
```

> **Cascade note (verified 2026-06-18, NOT cascading):** `grep -n "REFERENCES experiments\|REFERENCES samples\|REFERENCES exposures" packages/HimalayaUI/src/db.jl` confirms `samples.experiment_id`, `exposures.sample_id`, and every exposure-keyed structural table lack `ON DELETE CASCADE`. Phase A adds cascade ONLY to `exposures.experiment_id` and `loads.experiment_id`; the sample/exposure/index sub-tree does not cascade. That is why this route disables FK enforcement at the connection level (outside the transaction) and deletes the tree explicitly. Keep the delete list in sync if Phase B/D add new exposure- or index-keyed tables.

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_experiments.jl packages/HimalayaUI/test/test_ingestion_scan_api.jl
git commit -m "feat(api): DELETE /api/experiments/{id} with cascade + timer teardown"
```

---

## Task 7: Per-experiment rescan scheduler infrastructure (`server.jl`)

Add the `RESCAN_TIMERS` and `RESCAN_LOCKS` module-level constants and the three lifecycle functions (`start_rescan_scheduler!`, `stop_rescan_scheduler!`, `stop_all_rescan_timers!`). The scheduler itself (the tick logic and tiered backoff) is wired in Task 8; here we only establish the registry and the no-op skeleton that Task 6's `stop_rescan_scheduler!` call already needs.

**Files:**
- Modify: `packages/HimalayaUI/src/server.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_scan_api.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "rescan timer start/stop lifecycle (unit)" begin
    # This test never opens a real DB — it only exercises the timer registry.
    # start_rescan_scheduler! requires a db for the tick body; we pass a temp one.
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db;
            name = "T", path = dir, data_dir = dir, analysis_dir = dir)

        # Real Timers are armed below; guarantee they are all closed and the DB
        # is shut even if an assertion throws mid-test (otherwise a leaked Timer
        # keeps firing against a closed DB and poisons later testsets).
        try
            # Starting a scheduler registers a timer
            HimalayaUI.start_rescan_scheduler!(db, exp_id;
                tick_interval_seconds = 3600)
            @test haskey(HimalayaUI.RESCAN_TIMERS, exp_id)

            # Stopping it removes the entry
            HimalayaUI.stop_rescan_scheduler!(exp_id)
            @test !haskey(HimalayaUI.RESCAN_TIMERS, exp_id)

            # Stopping a non-existent id is a no-op
            @test_nowarn HimalayaUI.stop_rescan_scheduler!(exp_id)

            # stop_all_rescan_timers! clears all entries
            HimalayaUI.start_rescan_scheduler!(db, exp_id; tick_interval_seconds = 3600)
            HimalayaUI.stop_all_rescan_timers!()
            @test isempty(HimalayaUI.RESCAN_TIMERS)
        finally
            HimalayaUI.stop_all_rescan_timers!()
            SQLite.close(db)
        end
    end
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: FAIL (`RESCAN_TIMERS` undefined).

- [ ] **Step 3: Add registry + lifecycle functions to `server.jl`**

In `server.jl`, add after the `GC_TIMER` const and `stop_gc_timer!` (~line 168):

```julia
# Per-experiment rescan timer registry (spec §9.4).
# Keyed by experiment_id (Int). Each entry is a Timer whose tick runs the
# cheap change-check + additive scan on a @spawn'd task so the libuv timer
# thread is never blocked by multi-minute scan work.
# Distinct from GC_TIMER (a single module Ref). Not persisted across restarts;
# on boot, the scheduler is re-armed by the first scan or manual rescan.
const RESCAN_TIMERS = Dict{Int, Timer}()
const RESCAN_TIMERS_MU = ReentrantLock()

# Per-experiment reentrancy guard (spec §9.4: protect with an explicit
# ReentrantLock, not the DB_WRITE_LOCK which carries write-ordering semantics).
# Keyed by experiment_id. A timer tick attempts trylock; if already locked, the
# tick is skipped (a scan is already in flight).
const RESCAN_LOCKS = Dict{Int, ReentrantLock}()
const RESCAN_LOCKS_MU = ReentrantLock()

"""
    _rescan_lock(experiment_id) -> ReentrantLock

Return the per-experiment reentrancy lock, creating it on first access.
"""
function _rescan_lock(experiment_id::Int)
    lock(RESCAN_LOCKS_MU) do
        get!(RESCAN_LOCKS, experiment_id, ReentrantLock())
    end
end

"""
    stop_rescan_scheduler!(experiment_id)

Stop and remove the rescan Timer for `experiment_id`. No-op if no timer is
registered. Called before `DELETE /api/experiments/{id}` so the timer callback
cannot fire against a non-existent row.
"""
function stop_rescan_scheduler!(experiment_id::Int)
    lock(RESCAN_TIMERS_MU) do
        t = get(RESCAN_TIMERS, experiment_id, nothing)
        if t !== nothing
            close(t)
            delete!(RESCAN_TIMERS, experiment_id)
        end
    end
    nothing
end

"""
    stop_all_rescan_timers!()

Close all per-experiment rescan timers. Called on server shutdown (mirror
`stop_gc_timer!`).
"""
function stop_all_rescan_timers!()
    lock(RESCAN_TIMERS_MU) do
        for (_, t) in RESCAN_TIMERS
            try; close(t); catch; end
        end
        empty!(RESCAN_TIMERS)
    end
    nothing
end

"""
    start_rescan_scheduler!(db, experiment_id; tick_interval_seconds, cheap_check_fn)

Arm a per-experiment rescan Timer. On each tick the libuv timer thread does the
MINIMUM possible work — it only `@spawn`s the tick body. EVERYTHING that touches
the per-experiment `ReentrantLock` (both `trylock` AND `unlock`) runs INSIDE the
spawned task:

1. `@spawn` the body so the libuv timer thread is never blocked.
2. The spawned body `trylock`s the per-experiment lock; if already held (a scan
   is in flight) it returns immediately, skipping this tick.
3. With the lock held, it calls `_rescan_tick!`, then `unlock`s in a `finally`.

CRITICAL (Julia `ReentrantLock` is task-owned): the `trylock` and the matching
`unlock` MUST execute on the SAME task. If `trylock` ran on the timer thread and
`unlock` ran on the `@spawn`'d task, the unlock would throw "unlock from wrong
thread", the lock would never release, and every later tick would be silently
skipped. That is why both calls are inside the `@spawn` block below.

`cheap_check_fn` is an optional change-detection seam that defaults to the real
`cheap_change_check` (Phase B). Task 8's backoff test injects a stub here to drive
the tier transitions deterministically without filesystem fixtures; the production
path passes nothing and the real cheap-check runs. There is no `scan_fn` seam —
the tick calls `scan_and_group!` directly (a no-op on a directory with no matching
triplets, so the test's empty temp dir needs no stub).

Idempotent: calling twice closes the previous timer before installing a new one.
"""
function start_rescan_scheduler!(db::SQLite.DB, experiment_id::Int;
                                  tick_interval_seconds::Real = 3600.0,
                                  cheap_check_fn = nothing)
    # Close any existing timer for this experiment.
    stop_rescan_scheduler!(experiment_id)

    timer = Timer(tick_interval_seconds; interval = tick_interval_seconds) do _
        # Do NOTHING that touches the per-experiment lock on the timer thread.
        # Spawn first, then trylock/unlock both on the spawned task (same-task
        # ownership — see docstring).
        Threads.@spawn begin
            lk = _rescan_lock(experiment_id)
            trylock(lk) || return  # scan already in flight; skip this tick
            try
                _rescan_tick!(db, experiment_id; cheap_check_fn = cheap_check_fn)
            catch err
                @warn "rescan tick failed" experiment_id = experiment_id exception = err
            finally
                unlock(lk)
            end
        end
    end

    lock(RESCAN_TIMERS_MU) do
        RESCAN_TIMERS[experiment_id] = timer
    end
    nothing
end
```

> **Phase-B contract (complete, verified 2026-06-19):** both functions exist in
> `ingest.jl` and load with the module, so the tick calls them directly:
> - **`scan_and_group!(db, experiment_id; analyze=true, tif_pattern, prp_pattern, dat_pattern)`**
>   — resolves `data_dir`/`analysis_dir` from the experiment row internally. No
>   positional `root_dir`, no `additive` (idempotent via dedup INSERT keys).
> - **`cheap_change_check(db, experiment_id; image_pattern="{name}.tif")::Bool`**
>   — on-disk tif count vs persisted exposure count; returns `true` for an
>   unknown/unreadable experiment, `false` for a missing dir. The tick calls it
>   directly (the default of the `cheap_check_fn` seam).

Also update `stop_test_server!` (~line 211) to call `stop_all_rescan_timers!()`:

```julia
function stop_test_server!()
    stop_gc_timer!()
    stop_all_rescan_timers!()   # <- add this line
    Oxygen.terminate()
    Oxygen.resetstate()
    _DB_REF[] = nothing
    SSE_SUBSCRIBERS[] = []
    empty!(OP_LOCKS)
    nothing
end
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/server.jl packages/HimalayaUI/test/test_ingestion_scan_api.jl
git commit -m "feat(server): per-experiment rescan timer registry + lifecycle functions"
```

---

## Task 8: `_rescan_tick!` — tiered backoff + DB state

The tick body: runs the cheap change-check, conditionally triggers `scan_and_group!`, updates `consecutive_empty_ticks`, and applies tiered backoff by restarting the Timer with the appropriate interval. Backoff tiers (spec §9.4): `fast` (default `tick_interval_seconds`), `daily` (86400 s), `stopped`. Transitions: N consecutive empty ticks at `fast` → `daily`; M consecutive empty daily ticks → `stopped`. Re-arm to `fast` when a changed directory is detected.

**Files:**
- Modify: `packages/HimalayaUI/src/server.jl` (add `_rescan_tick!`)
- Test: `packages/HimalayaUI/test/test_ingestion_scan_api.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "_rescan_tick! tiered backoff persists to DB" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db;
            name = "BT", path = dir, data_dir = dir, analysis_dir = dir)

      # The "change" branch below calls start_rescan_scheduler!, which arms a REAL
      # Timer. Wrap the whole body so that timer (and any other) is always closed
      # and the DB shut even if an assertion throws — a leaked Timer would keep
      # firing against a closed DB.
      try
        ticks_before_daily = 3   # configurable in the call below

        # Inject the change decision so tier transitions are deterministic without
        # filesystem fixtures. No scan stub is needed: on a "change" the real
        # scan_and_group! runs against the EMPTY temp dir (no matching triplets) —
        # a harmless no-op.
        no_change  = (_, _) -> false
        has_change = (_, _) -> true

        # Run ticks_before_daily empty ticks; should stay in 'fast' tier until threshold
        for _ in 1:ticks_before_daily
            HimalayaUI._rescan_tick!(db, exp_id;
                cheap_check_fn = no_change,
                fast_interval = 3600.0, ticks_before_daily = ticks_before_daily,
                ticks_before_stop = 2)
        end
        row = Tables.rowtable(DBInterface.execute(db,
            "SELECT last_scan_tier, consecutive_empty_ticks FROM experiments WHERE id=?",
            [exp_id]))[1]
        # After exactly ticks_before_daily empty ticks, tier advances to daily
        @test row.last_scan_tier == "daily"
        @test row.consecutive_empty_ticks == 0  # reset on tier transition

        # One more empty daily tick → not yet stopped (need ticks_before_stop=2 more)
        HimalayaUI._rescan_tick!(db, exp_id;
            cheap_check_fn = no_change,
            fast_interval = 3600.0, ticks_before_daily = ticks_before_daily,
            ticks_before_stop = 2)
        row2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT last_scan_tier, consecutive_empty_ticks FROM experiments WHERE id=?",
            [exp_id]))[1]
        @test row2.last_scan_tier == "daily"
        @test row2.consecutive_empty_ticks == 1

        # Simulate a change: re-arms back to fast. cheap_check_fn returns true, then
        # the real scan_and_group! runs against the EMPTY temp dir — a harmless no-op.
        HimalayaUI._rescan_tick!(db, exp_id;
            cheap_check_fn = has_change,
            fast_interval = 3600.0, ticks_before_daily = ticks_before_daily,
            ticks_before_stop = 2)
        row3 = Tables.rowtable(DBInterface.execute(db,
            "SELECT last_scan_tier, consecutive_empty_ticks FROM experiments WHERE id=?",
            [exp_id]))[1]
        @test row3.last_scan_tier == "fast"
        @test row3.consecutive_empty_ticks == 0
      finally
        # The change-branch re-armed a fast-tier Timer; close it before the DB.
        HimalayaUI.stop_rescan_scheduler!(exp_id)
        HimalayaUI.stop_all_rescan_timers!()
        SQLite.close(db)
      end
    end
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: FAIL (`_rescan_tick!` undefined).

- [ ] **Step 3: Implement `_rescan_tick!`**

Add to `server.jl` after the scheduler functions from Task 7:

```julia
"""
    _rescan_tick!(db, experiment_id; cheap_check_fn=nothing, kwargs...)

One tick of the per-experiment rescan scheduler. Called from the `@spawn`'d
body inside the Timer callback so the libuv timer thread is never blocked.

Change detection runs `cheap_check_fn(db, experiment_id)` when a seam is injected
(Task 8's backoff test), otherwise the real `cheap_change_check(db, experiment_id)`
(Phase B). On a detected change it runs `scan_and_group!(db, experiment_id)`. Both
Phase B functions resolve the experiment's `data_dir` from the row themselves, so
this tick threads no path.

Steps:
1. Read `(last_scan_tier, consecutive_empty_ticks)`; if the experiment is gone, return.
2. Determine `changed::Bool`.
3. If changed: run the scan, reset `consecutive_empty_ticks`, re-arm at the `fast` tier.
4. If not changed: increment the counter; advance fast→daily at `ticks_before_daily`,
   daily→stopped at `ticks_before_stop`, resetting the counter on each transition.
5. Persist `last_scan_tier` + `consecutive_empty_ticks` so restarts don't reset
   quiet experiments to the fast tier.

kwargs (tunable thresholds — defaults match spec §9.4; exposed for testing):
- `fast_interval`: seconds between fast-tier ticks (default 3600.0 = 1 h)
- `daily_interval`: seconds between daily-tier ticks (default 86400.0 = 24 h)
- `ticks_before_daily`: consecutive empty fast ticks before advancing to daily (default 6)
- `ticks_before_stop`: consecutive empty daily ticks before stopping (default 3)
"""
function _rescan_tick!(db::SQLite.DB, experiment_id::Int;
                        cheap_check_fn = nothing,
                        fast_interval::Real    = 3600.0,
                        daily_interval::Real   = 86400.0,
                        ticks_before_daily::Int = 6,
                        ticks_before_stop::Int  = 3)

    # Single existence + state read. Empty ⇒ experiment deleted between timer arm
    # and tick (the DELETE route stops the timer first, so this is belt-and-suspenders).
    row = Tables.rowtable(DBInterface.execute(db,
        "SELECT last_scan_tier, consecutive_empty_ticks FROM experiments WHERE id = ?",
        [experiment_id]))
    isempty(row) && return

    # Change detection: injected seam (test) or the real cheap_change_check. Both
    # take (db, experiment_id) and resolve data_dir themselves.
    changed = try
        cheap_check_fn !== nothing ? cheap_check_fn(db, experiment_id) :
                                     cheap_change_check(db, experiment_id)
    catch err
        @warn "cheap change-check failed" experiment_id = experiment_id exception = err
        false
    end

    if changed
        try
            scan_and_group!(db, experiment_id)
        catch err
            @warn "scan failed during rescan tick" experiment_id = experiment_id exception = err
        end
        # Re-arm at fast tier on a detected change.
        lock(_DB_WRITE_LOCK) do
            SQLite.transaction(db) do
                DBInterface.execute(db, """
                    UPDATE experiments
                       SET last_scan_tier         = 'fast',
                           consecutive_empty_ticks = 0
                     WHERE id = ?
                """, [experiment_id])
            end
        end
        start_rescan_scheduler!(db, experiment_id;
            tick_interval_seconds = fast_interval)
        return
    end

    # No change — increment tick counter and check backoff thresholds.
    tier      = String(row[1].last_scan_tier)
    new_ticks = Int(row[1].consecutive_empty_ticks) + 1

    if tier == "fast" && new_ticks >= ticks_before_daily
        # Advance to daily
        lock(_DB_WRITE_LOCK) do
            SQLite.transaction(db) do
                DBInterface.execute(db, """
                    UPDATE experiments
                       SET last_scan_tier         = 'daily',
                           consecutive_empty_ticks = 0
                     WHERE id = ?
                """, [experiment_id])
            end
        end
        start_rescan_scheduler!(db, experiment_id;
            tick_interval_seconds = daily_interval)
    elseif tier == "daily" && new_ticks >= ticks_before_stop
        # Stop the scheduler
        lock(_DB_WRITE_LOCK) do
            SQLite.transaction(db) do
                DBInterface.execute(db, """
                    UPDATE experiments
                       SET last_scan_tier         = 'stopped',
                           consecutive_empty_ticks = 0
                     WHERE id = ?
                """, [experiment_id])
            end
        end
        stop_rescan_scheduler!(experiment_id)
    else
        # Stay in current tier, just increment the counter
        lock(_DB_WRITE_LOCK) do
            SQLite.transaction(db) do
                DBInterface.execute(db, """
                    UPDATE experiments
                       SET consecutive_empty_ticks = ?
                     WHERE id = ?
                """, [new_ticks, experiment_id])
            end
        end
    end
    nothing
end
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/server.jl packages/HimalayaUI/test/test_ingestion_scan_api.jl
git commit -m "feat(server): _rescan_tick! tiered backoff (fast→daily→stopped), persisted to DB"
```

---

## Task 9: `POST /api/experiments/{id}/scan` — rescan (cheap change-check + additive)

The rescan endpoint runs the cheap change-check first; if nothing changed, it skips the scan. If new files appeared, it calls `scan_and_group!(db, exp_id)`, broadcasts progress via `broadcast_progress!`, and re-arms the rescan scheduler at the fast tier (spec §9.2). (The old `POST /api/experiments/{id}/reingest` route was already removed in Phase B.)

**Files:**
- Modify: `packages/HimalayaUI/src/routes_experiments.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_scan_api.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "POST /api/experiments/{id}/scan no-change is idempotent" begin
    db, dir, exp_id = scan_test_db()

    with_test_server(db) do port, base
        # First scan on empty dir: no files, should return 200 + changed=false
        # (since scan_and_group! is Phase B and may not be present in test env,
        #  the route must not crash when the directory has no matching files;
        #  test the HTTP contract, not the full scan logic)
        r = HTTP.post("$base/api/experiments/$exp_id/scan";
            headers = ["X-Username" => "alice"],
            status_exception = false)
        # 200 or 202 are both acceptable (implementation may vary)
        @test r.status in (200, 202)
        body = JSON3.read(String(r.body))
        @test haskey(body, :status)

        # 404 for unknown experiment
        r2 = HTTP.post("$base/api/experiments/999999/scan";
            headers = ["X-Username" => "alice"],
            status_exception = false)
        @test r2.status == 404
    end
    SQLite.close(db)
end
```

> The test is intentionally lightweight here because `scan_and_group!` is Phase B. The route's contract (correct HTTP status, 404 for unknown, no crash) is what this plan can assert without the full grouping engine. The Phase B plan should add a deeper integration test for the scan logic.

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: FAIL (404 — route does not exist).

- [ ] **Step 3: Implement the route**

In `routes_experiments.jl`, add inside `register_experiments_routes!()`:

```julia
@post "/api/experiments/{id}/scan" function(req::HTTP.Request, id::Int)
    db   = current_db()
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM experiments WHERE id = ?", [id]))
    isempty(rows) && return HTTP.Response(404,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "experiment not found")))

    # Mark as scanning immediately so the frontend header can show progress.
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            DBInterface.execute(db,
                "UPDATE experiments SET ingest_status = 'scanning' WHERE id = ?", [id])
        end
    end

    broadcast_progress!(id; kind = "ingest_started", processed = 0, total = 0)

    # Run the cheap change-check + additive scan on a @spawn'd task so this request
    # returns immediately; progress streams over SSE. Both Phase B functions
    # (cheap_change_check, scan_and_group!) resolve the experiment's data_dir from
    # the row themselves, and return gracefully on an empty directory.
    Threads.@spawn begin
        try
            changed = cheap_change_check(db, id)
            if changed
                scan_and_group!(db, id)
                start_rescan_scheduler!(db, id)   # re-arm the fast-tier scheduler
            end

            lock(_DB_WRITE_LOCK) do
                SQLite.transaction(db) do
                    DBInterface.execute(db,
                        "UPDATE experiments SET ingest_status = 'complete', last_scanned_at = ? WHERE id = ?",
                        [format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ"), id])
                end
            end
            broadcast_progress!(id; kind = "ingest_complete",
                processed = 0, total = 0, changed = changed)
        catch err
            @warn "scan failed" experiment_id = id exception = err
            lock(_DB_WRITE_LOCK) do
                SQLite.transaction(db) do
                    DBInterface.execute(db,
                        "UPDATE experiments SET ingest_status = 'failed' WHERE id = ?", [id])
                end
            end
            broadcast_progress!(id; kind = "ingest_failed",
                processed = 0, total = 0,
                error = sprint(showerror, err))
        end
    end

    log_action!(db, req; action = "scan",
        entity_type = "experiment", entity_id = id)

    HTTP.Response(202, ["Content-Type" => "application/json"],
        JSON3.write(Dict(:status => "scanning", :experiment_id => id)))
end
```

> **Import note (verified):** all `src/*.jl` files are `include`d into the single `HimalayaUI` module (no submodules — see `HimalayaUI.jl`), and `events.jl:2` already does `using Dates: now, UTC, format, @dateformat_str`. Those four names are therefore in scope module-wide, including here. Use the **unqualified** `format(now(UTC), dateformat"…")` form — a qualified `Dates.format`/`Dates.now` is an `UndefVarError` because bare `Dates` is never imported. Do NOT add `using Dates`. `cheap_change_check` and `scan_and_group!` are Phase B functions in the same module (no `isdefined` guards needed); call them directly.

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_experiments.jl packages/HimalayaUI/test/test_ingestion_scan_api.jl
git commit -m "feat(api): POST /api/experiments/{id}/scan — async additive rescan + SSE progress"
```

---

## Task 10: `POST /api/experiments` — create-from-directory + async first scan

The new route replaces the pre-create-`experiment.toml`-then-`init` dance. It receives `{ path, name?, patterns? }`, creates the experiment row immediately (using the directory basename as the default name), kicks off the first scan asynchronously with SSE progress, and returns `{ id, status: "scanning" }` (spec §9.2).

**Files:**
- Modify: `packages/HimalayaUI/src/routes_experiments.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_scan_api.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "POST /api/experiments creates experiment + starts async scan" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "himalaya.db"))

        with_test_server(db) do port, base
            r = HTTP.post("$base/api/experiments";
                body    = JSON3.write(Dict(:path => dir, :name => "MyExp")),
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice"],
                status_exception = false)
            @test r.status in (200, 201, 202)
            body = JSON3.read(String(r.body))
            @test haskey(body, :id)
            @test body.id isa Integer
            @test body.id > 0

            # Experiment row is created
            exp_rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, name, data_dir, ingest_status FROM experiments WHERE id = ?",
                [body.id]))
            @test length(exp_rows) == 1
            @test exp_rows[1].name == "MyExp"
            @test exp_rows[1].data_dir == dir

            # ingest_status transitions to scanning (may have already completed in test)
            @test exp_rows[1].ingest_status in ("scanning", "complete", "failed")

            # Missing path → 400
            r2 = HTTP.post("$base/api/experiments";
                body    = JSON3.write(Dict(:name => "NoPath")),
                headers = ["Content-Type" => "application/json"],
                status_exception = false)
            @test r2.status == 400

            # Non-existent path → 400
            r3 = HTTP.post("$base/api/experiments";
                body    = JSON3.write(Dict(:path => "/does/not/exist/xyz123")),
                headers = ["Content-Type" => "application/json"],
                status_exception = false)
            @test r3.status == 400
        end
        SQLite.close(db)
    end
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: FAIL (route does not exist — returns 404 or 405).

- [ ] **Step 3: Implement the route**

In `routes_experiments.jl`, add inside `register_experiments_routes!()` **before** the existing `@get "/api/experiments"` handler:

```julia
@post "/api/experiments" function(req::HTTP.Request)
    db   = current_db()
    body = json(req)

    # Required: path to the data directory.
    path_val = get(body, :path, get(body, "path", nothing))
    path_val === nothing && return HTTP.Response(400,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "path is required")))
    data_dir = String(path_val)

    isdir(data_dir) || return HTTP.Response(400,
        ["Content-Type" => "application/json"],
        JSON3.write(Dict(:error => "path does not exist or is not a directory",
                         :path  => data_dir)))

    # Derive defaults.
    name_val  = get(body, :name, get(body, "name", nothing))
    exp_name  = name_val !== nothing ? String(name_val) : basename(rstrip(data_dir, '/'))
    # analysis_dir convention: look for an `analysis` subdirectory; fall back to data_dir.
    analysis_dir = let ad = joinpath(data_dir, "analysis")
        isdir(ad) ? ad : data_dir
    end

    exp_id = lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            create_experiment!(db;
                name         = exp_name,
                path         = data_dir,
                data_dir     = data_dir,
                analysis_dir = analysis_dir,
                ingest_status = "scanning")
        end
    end

    broadcast_progress!(exp_id; kind = "ingest_started", processed = 0, total = 0)

    # Kick off first scan asynchronously.
    Threads.@spawn begin
        try
            # scan_and_group! (Phase B, ingest.jl) resolves data_dir from the row
            # and is idempotent (dedup INSERT keys), so first-scan == rescan.
            scan_and_group!(db, exp_id)
            lock(_DB_WRITE_LOCK) do
                SQLite.transaction(db) do
                    DBInterface.execute(db,
                        "UPDATE experiments SET ingest_status = 'complete', last_scanned_at = ? WHERE id = ?",
                        [format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ"), exp_id])
                end
            end
            broadcast_progress!(exp_id; kind = "ingest_complete", processed = 0, total = 0)
            # Arm the rescan scheduler after a successful first scan.
            start_rescan_scheduler!(db, exp_id)
        catch err
            @warn "first scan failed" experiment_id = exp_id exception = err
            lock(_DB_WRITE_LOCK) do
                SQLite.transaction(db) do
                    DBInterface.execute(db,
                        "UPDATE experiments SET ingest_status = 'failed' WHERE id = ?", [exp_id])
                end
            end
            broadcast_progress!(exp_id; kind = "ingest_failed",
                processed = 0, total = 0, error = sprint(showerror, err))
        end
    end

    log_action!(db, req; action = "experiment_created",
        entity_type = "experiment", entity_id = exp_id)

    HTTP.Response(202, ["Content-Type" => "application/json"],
        JSON3.write(Dict(:id => exp_id, :status => "scanning",
                         :name => exp_name, :data_dir => data_dir)))
end
```

> **`create_experiment!` call note:** the Phase A plan extended `create_experiment!` to accept geometry kwargs. The current live signature (read `db.jl:1808`) requires `path`, `data_dir`, `analysis_dir` as positional or keyword args. Confirm the exact signature before this commit and adjust the call as needed. The `ingest_status` kwarg was added in Phase A; if Phase A used `ingest_status = "idle"` as the default, passing it explicitly here overrides to `"scanning"`.

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_experiments.jl packages/HimalayaUI/test/test_ingestion_scan_api.jl
git commit -m "feat(api): POST /api/experiments create-from-directory + async first scan"
```

---

## Task 11: Full route registration + load sanity check

Register the new routes in `HimalayaUI.jl` (or confirm they are registered inside `register_experiments_routes!` — which is already called from `register_routes!`), then run a clean module load and the full backend suite to catch any regression.

**Files:**
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl` (add includes if new `.jl` files were created — in this plan, none were; all new code is in existing files)
- Test: full suite

- [ ] **Step 1: Confirm route registration wiring**

Read the current `register_experiments_routes!()` body and confirm every new route (`@post "/api/experiments"`, `@post "/api/experiments/{id}/scan"`, `@get "/api/experiments/{id}/loads"`, `@delete "/api/experiments/{id}"`, and the new `@patch "/api/experiments/{id}"`) is present inside it. All routes in this plan are added inside the existing function and registered transitively via `register_routes!()` → `register_experiments_routes!()` in `server.jl`. No new includes are needed in `HimalayaUI.jl`.

Run: `grep -n "@post\|@get\|@patch\|@delete" packages/HimalayaUI/src/routes_experiments.jl`
Confirm all expected routes appear.

- [ ] **Step 2: Module load check**

Run: `julia --project=packages/HimalayaUI -e 'using HimalayaUI; println("OK")'`
Expected: prints `OK` with no errors or warnings about undefined names.

- [ ] **Step 3: Run the full Phase C test file**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_scan_api.jl`
Expected: all testsets PASS.

- [ ] **Step 4: Run the full backend suite**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test-c.out 2>&1
tail -60 /tmp/jl-test-c.out
```
Expected: all existing tests pass. Investigate any failure; likely causes:
- The `PATCH /api/experiments/{id}` test in `test_routes_experiments.jl` asserts the old 400 response — update it to expect 200 and check the geometry write (or find the exact assertion and update the expectation).
- `stop_test_server!` now calls `stop_all_rescan_timers!()` — any test that leaves a Timer running after `stop_test_server!` may observe a slightly different teardown; confirm no `Timer closed` warnings that mask test failures.

- [ ] **Step 5: Fix any regressions and commit**

```bash
git add packages/HimalayaUI/src/ packages/HimalayaUI/test/
git commit -m "fix(tests): update PATCH experiment expectation + rescan teardown"
```

---

## Self-Review

### Spec coverage

| Spec requirement | Task |
|---|---|
| `POST /api/experiments` (create-from-directory, async first scan) | Task 10 |
| `POST /api/experiments/{id}/scan` (cheap check, additive ingest, SSE progress) | Task 9 |
| `GET /api/experiments/{id}` extended (typed geometry + `*_source` + ingest status + stats) | Task 3 |
| `GET /api/experiments/{id}/loads` (Load▸Sample▸Exposure roll-up, dedicated endpoint) | Task 4 |
| Geometry PATCH replacing the read-only 400 stub (`routes_experiments.jl:73-85`) | Task 5 |
| `DELETE /api/experiments/{id}` with timer teardown | Task 6 |
| `broadcast_progress!` (`_try_put!` direct, `curation` event, `kind` + `experiment_id`, no `user_actions` row) | Task 1 |
| Per-experiment rescan scheduler (`Dict{Int,Timer}`, `ReentrantLock` per experiment, `@spawn`) | Task 7 |
| Tiered backoff persisted to `last_scan_tier`/`consecutive_empty_ticks` | Task 8 |
| `create_load!` writer (prerequisite for all scanner callers) | Task 2 |
| `stop_all_rescan_timers!` on server shutdown + timer teardown on experiment delete | Tasks 7 + 6 |

### Placeholder scan

Every code step shows the complete implementation. There are no `isdefined` guards or placeholders: Phase B is complete and in the same module, so Tasks 7–10 call `scan_and_group!(db, id)` and `cheap_change_check(db, id)` directly.

### Phase-B contract (complete, verified 2026-06-19)

- **`scan_and_group!(db, experiment_id; analyze=true, tif_pattern, prp_pattern, dat_pattern)`** — defined by Phase B in `ingest.jl`. Resolves `data_dir`/`analysis_dir` from the experiment row internally (**no positional `root_dir`**); **no `additive` kwarg** (idempotent via dedup INSERT keys, so first-scan == rescan). Every caller in this phase (scan route, create route, `_rescan_tick!`) calls it as `scan_and_group!(db, id)`.
- **`cheap_change_check(db, experiment_id; image_pattern="{name}.tif")::Bool`** — defined by Phase B in `ingest.jl`. On-disk tif count vs persisted exposure count; `true` for an unknown/unreadable experiment, `false` for a missing dir. Called directly (the default of `_rescan_tick!`'s `cheap_check_fn` seam, which the backoff test overrides).

### Deferred / out of scope (this phase)

- **§9.2 directory-picker backend (path autocomplete + validate-path):** the `DirectoryPickerField` endpoints listed in spec §9.2 as a create-flow prerequisite are intentionally deferred per open question §13. `POST /api/experiments` validates the supplied path with a plain `isdir` check; live filesystem autocomplete + a dedicated validate-path endpoint are a later increment, not part of Phase C.

### Review applied (2026-06-18)

A verify-before-review pass against live source folded in these fixes:
- **P0** rescan scheduler: `trylock`/`unlock` moved BOTH inside the `@spawn`'d task (Julia `ReentrantLock` is task-owned; cross-task unlock throws and wedges the timer) — Task 7.
- **P0** `broadcast_progress!` + Tasks 9/10 timestamps: use unqualified `format(now(UTC), dateformat"…")` (events.jl imports Dates selectively; bare `Dates` is unbound module-wide) — Tasks 1, 9, 10.
- **P0** DELETE route: no `ON DELETE CASCADE` exists on the sample/exposure/index sub-tree; replaced cascade-reliance with an explicit FK-ordered teardown under connection-level `PRAGMA foreign_keys = OFF` (the migration idiom) — Task 6.
- **P1** Task 5 `const`s moved to module top level (a `const` inside `register_experiments_routes!()` is a syntax error).
- **P1** Phase-B handoff unified on `isdefined`-resolved names (`scan_and_group!`, future `cheap_change_check`) with the positional `root_dir`; dropped the non-shadowing `_default_*` stubs — Tasks 7–10.
- **P1** `broadcast_progress!` frame is PAYLOAD-WRAPPED to match `broadcast_event!` (events.jl:1107-1120) and spec §9.3/§9.6: `:kind` (+ `:ts`) top-level, `experiment_id`/`processed`/`total`/`kwargs` under `:payload`. An earlier flat-frame fix put those fields at top level, which would NOT parse through the frontend's `curation` payload path (E1/E2 read `payload.experiment_id`); reverted to wrapped. Docstring signature → `kwargs...` — Task 1.
- **P2** dropped the no-op in-transaction `PRAGMA foreign_keys = ON` from the delete; loads roll-up exposures now serialize with `bool_keys=(:selected,)`; legacy `_experiment_row_to_json` fallback now emits a uniform wire shape; timer-leak try/finally added to the Task 7/8 tests.

### Revision (2026-06-19) — Phase B contract correction + ponytail cuts

Phase B shipped (and was ponytail-trimmed) after this plan was first drafted, invalidating its central assumptions. This revision realigns to the real `ingest.jl`:
- **Stale contract fixed:** `scan_and_group!` has no positional `root_dir` and no `additive` kwarg (both dropped in Phase B); `cheap_change_check` exists (was assumed a "gap"). Every call in Tasks 8–10 and `_rescan_tick!` is now `scan_and_group!(db, id)` / `cheap_change_check(db, id)`.
- **Dropped the `isdefined` apparatus** (Tasks 7–10): Phase B is in the same module and loaded at `include` time, so the by-name resolution + assume-changed-gap fallback are dead. Direct calls.
- **Dropped the `scan_fn` injection seam**; kept one `cheap_check_fn` seam defaulted to the real `cheap_change_check`. The backoff test injects the change-decision and lets the real scan run as a no-op on its empty temp dir — no `scan_fn` stub, no filesystem fixtures.
- **`_rescan_tick!` no longer reads/threads `data_dir`** — both callees resolve it themselves; a single state read at the top doubles as the existence guard.
- Verified `create_load!` is genuinely absent (Task 2 stands); `manifest_path` column still exists (Task 5 readonly-reject entry stands). Removed stale "supersedes `/reingest`" notes (that route was deleted in Phase B).

### Type/name consistency

- `broadcast_progress!(experiment_id::Integer; kind, processed, total, kwargs...)` — emits a payload-wrapped frame; matches the test assertions (`obj.kind` top-level; `obj.payload.experiment_id`, `obj.payload.processed`, `obj.payload.total`).
- `RESCAN_TIMERS::Dict{Int,Timer}`, `RESCAN_LOCKS::Dict{Int,ReentrantLock}` — referenced identically in Task 7 (definition), Task 6 (route calls `stop_rescan_scheduler!`), and Task 11 (registration check).
- `start_rescan_scheduler!(db, experiment_id; tick_interval_seconds, cheap_check_fn)` and `stop_rescan_scheduler!(experiment_id)` — identical signature in definition (Task 7), the tick re-arm call (Task 8), and the scan route calls (Tasks 9, 10).
- `_rescan_tick!(db, experiment_id; cheap_check_fn, fast_interval, daily_interval, ticks_before_daily, ticks_before_stop)` — kwargs match test call and implementation; `cheap_check_fn` accepts a callable (test stub) or `nothing` (production → real `cheap_change_check`). No `scan_fn` seam; the scan call is direct.
- `_experiment_row_to_json(row, db)` — 2-arg form used in the extended GET route (Task 3); 1-arg form (db=nothing) preserved for the list route backward compat.
- `_experiment_stats(db, exp_id)` — defined and used only in Task 3; returns a NamedTuple embedded as `:stats` in the response.

### Known accuracy risks

- `create_experiment!` kwarg names: the corrected Plan A Task 6 signature (`name`, required `path`/`data_dir`/`analysis_dir`, geometry kwargs, `ingest_status = "idle"` default) was confirmed against Plan A; this phase's calls pass `name`/`path`/`data_dir`/`analysis_dir`/`ingest_status`, all valid.
- DELETE route FK teardown: verified that `samples.experiment_id` (db.jl:34), `exposures.sample_id` (db.jl:42), and the exposure/index-keyed sub-tree (db.jl:60-152) do NOT cascade; the route deletes the tree explicitly under connection-level `PRAGMA foreign_keys = OFF` (migration idiom, db.jl:1620-1647). Keep the delete list in sync if Phase B/D add new exposure- or index-keyed tables.
- `Dates`: NOT imported as a bare module anywhere; `events.jl:2` brings `now, UTC, format, @dateformat_str` into module scope (all files share the single `HimalayaUI` module). All timestamp code uses the unqualified forms — verified.
- `@delete` macro: confirmed in live source (`routes_series.jl:266`, `routes_peaks.jl:348`, etc. all use `@delete`). No change needed.
