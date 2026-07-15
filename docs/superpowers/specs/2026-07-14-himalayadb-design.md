# HimalayaDB.jl — design spec

**Date:** 2026-07-14
**Status:** Approved for implementation planning
**Scope:** A standalone, read-only Julia package giving programmatic access to the curated SAXS annotations produced by HimalayaUI.

## Purpose

Let a scientist pull HimalayaUI's curated annotations — the effective peak set, index candidates, and confirmed phase assignments — into a Julia session for downstream **analysis and plotting**, reading the SQLite file **directly and offline** (no running server). The package encapsulates the non-obvious domain logic (chiefly the effective-peaks UNION) that raw SQL users would otherwise reimplement and get wrong.

## Decisions (settled during brainstorming)

1. **Access mode:** reads the SQLite file directly, offline. Not a REST client.
2. **Traces:** DB-only core; trace-curve loading is a separate opt-in function. Annotations always work from just the DB file.
3. **Return shape:** Tables.jl rows by default (analysis/plotting path) + opt-in reconstruction of core `Himalaya.Index{P}` types. Plus a DataFrames weakdep extension for convenience.
4. **Relationship to HimalayaUI:** standalone, read-only. Does **not** extract HimalayaUI's DB code. Duplicates only the small load-bearing query logic, guarded by a contract test. Extraction/de-duplication is explicitly deferred.

## Non-goals (v1)

- **No writes, ever.** No migrations, no `chmod`, no `apply_event!`, no SSE, no idempotency. Those are HimalayaUI application/event concerns and stay there.
- **No HTTP / Oxygen dependency.**
- **Does not own the schema.** It `SELECT`s against tables HimalayaUI creates; it never `CREATE`s or `ALTER`s.
- **Out of v1 read surface:** `series`/`comparisons` (curated multi-trace groupings), `sample_tags`/`exposure_tags`, chat messages, `user_actions` history/attribution. Deferred (YAGNI); add if a real need appears.
- **No extraction of HimalayaUI's DB layer.** Deferred until the read surface proves itself.

## Architecture

A new package at `packages/HimalayaDB/`, sibling to `packages/HimalayaUI/` in the monorepo. Monorepo placement is deliberate: it makes the contract test (below) trivial because both packages are present in one checkout.

```
packages/HimalayaDB/
  Project.toml            # deps + [weakdeps]/[extensions] for DataFrames
  src/HimalayaDB.jl       # module, includes, exports
  src/connect.jl          # read-only DB open + path resolution
  src/queries.jl          # the curated read functions (SQL)
  src/reconstruct.jl      # opt-in Himalaya.Index{P} reconstruction
  src/trace.jl            # opt-in .dat trace-curve loading
  ext/DataFramesExt.jl    # weakdep: dataframe() convenience
  test/
    runtests.jl
    fixture.jl            # build a small deterministic DB for tests
    test_queries.jl       # unit tests on the fixture
    test_contract.jl      # cross-check vs HimalayaUI (test-only dep)
```

### Data source (from the codebase map)

All anchors are `packages/HimalayaUI/src/` unless noted.

- One central SQLite DB. Path resolves from `ENV["HIMALAYA_DB_PATH"]`, else `~/.himalaya/himalaya.db` (`db.jl:1750`, `default_db_path`).
- WAL journal mode (`db.jl:1797`) → read-only connections can coexist with the live writer.
- Schema: `const SCHEMA` (`db.jl:7-193`) + migration functions.
- **Curation is table-based, not a flag:**
  - `auto_peaks` (`db.jl:97-105`) — algorithm-detected peaks.
  - `peak_curations` (`db.jl:107-114`) — human edits, `kind IN ('exclude','add')`, with `created_by`.
  - **Effective (curated) peak set is computed, never stored** = `auto_peaks − excludes ∪ adds`, tolerance-matched on q. Reference implementations: `effective_peaks` (`pipeline.jl:103`) and the read-shaped `get_peaks_for_exposure` (`pipeline.jl:600-618`, UNION tagging `source ∈ {auto,manual}` + `excluded`).
  - `indices` (`db.jl:72-83`) — candidate index-choices; `kind ∈ {auto, speculative}`.
  - Confirmed phase = active/`custom` `index_groups` (`db.jl:121-128`) via `index_group_members` (`db.jl:133-137`), written by `index_confirmed`/`index_unconfirmed` events.
- **Event log is source of truth; view tables are materialized and provably re-derivable** (`rebuild_views_from_log!`, `events.jl:1034`). Consequence: a reader reads the view tables directly — no event replay needed.
- Trace (q, I) is **not** in SQLite. It lives in per-experiment `.dat` files under the experiment's `analysis_dir`; parser is `datfile.jl`. `exposures` (`db.jl:45-55`) holds `filename`, `image_path`, `trace_hash` — the pointers, not the arrays.
- Core `Himalaya` types (repo-root `src/`): `Index{P<:Phase}` (`src/index.jl:7-11`) + the `Phase` hierarchy (`src/phase.jl:36-45`). **No `Peak` struct** — peaks are NamedTuples / rows. Phase short name is stored in `indices.phase` and reconstructed via `getfield(Himalaya, Symbol(name))`, validated `P isa Type && P <: Himalaya.Phase`.

## Components

### 1. `connect.jl` — read-only open

```
connect(path = get(ENV, "HIMALAYA_DB_PATH", default_himalaya_db_path())) -> SQLite.DB
```

- Opens **read-only**: `mode=ro` URI + `PRAGMA query_only = ON`. **No migrations, no `chmod`** — unlike `HimalayaUI.open_db` (`db.jl:1767`), which mutates the file (runs migrations, `chmod` to 0664 at `db.jl:1819-1826`).
- `default_himalaya_db_path()` mirrors HimalayaUI's default (`~/.himalaya/himalaya.db`) but must **not** create directories.

> **Known implementation risk (resolve in the plan):** opening a WAL DB read-only *while HimalayaUI is concurrently writing* requires directory write access to create the `-shm`/`-wal` sidecars; `immutable=1` is **not** safe under a concurrent writer. For the primary offline use case (app not running) a plain `mode=ro` open is fine. The plan must pin the exact open string and cover both the offline and concurrent-reader cases, and decide what `connect` does when only the `.db` file (no sidecars, no dir write access) is present.

### 2. `queries.jl` — the curated read surface

All return Tables.jl-compatible rows (NamedTuple rowtables via `Tables.rowtable(DBInterface.execute(...))`, matching HimalayaUI's convention where SQL `NULL` → `missing`).

Navigation hierarchy:
- `experiments(db)` — all experiments.
- `samples(db; experiment)` — samples, optionally filtered by experiment id.
- `exposures(db; sample)` — exposures, optionally filtered by sample id.

Curated annotations (the headline surface):
- `curated_peaks(db, exposure_id)` — the **effective peak set**: `auto_peaks ∪ adds − excludes`, each row tagged `source ∈ {auto, manual}` and `excluded`. Port of `get_peaks_for_exposure` (`pipeline.jl:600`). This is the ~20-line UNION that is the package's reason to exist.
- `index_candidates(db, exposure_id)` — all `indices` rows (auto + speculative), ordered by score. Port of `get_indices_for_exposure` (`pipeline.jl:654`).
- `confirmed_indices(db, exposure_id)` — the human-confirmed subset (active `index_group_members`). Answers "what phase did the curator actually commit to." Port of the confirmed-group query behind `get_groups_for_exposure` (`pipeline.jl:662`).

### 3. `reconstruct.jl` — opt-in core types

- `reconstruct_index(db, index_id) -> Himalaya.Index{P}` — rebuilds the core type from an `indices` row + its `index_peaks` supporting peaks, giving `predictpeaks` / `missingpeaks` / `phaseratios`. Sole reason HimalayaDB depends on core `Himalaya`. Reconstruction logic mirrors HimalayaUI's `_serialized_indices_for_broadcast` join (`pipeline.jl:991`) but returns the typed `Index{P}` rather than the JSON broadcast shape.

### 4. `trace.jl` — opt-in trace loading

- `load_trace(db, exposure_id) -> (q, I)` — resolves the exposure's `.dat` path from the DB (exposure `filename` + experiment `analysis_dir`), then parses it. Requires a small `.dat` reader — the one piece of duplicated logic beyond SQL (a few lines; port the relevant slice of `datfile.jl`). Errors clearly if the on-disk data dir is absent (offline-with-DB-only case).

### 5. `ext/DataFramesExt.jl` — DataFrames convenience (weakdep)

- `DataFrames` declared under `[weakdeps]` + `[extensions]` in `Project.toml`. The extension provides `dataframe(rows) -> DataFrame` and loads only when the user has `using DataFrames`. Base package stays lightweight; results are already `DataFrame`-ready via Tables.jl, so this is pure ergonomics (a named, tested helper) rather than new capability.

## Data flow

```
user (Julia REPL / notebook, offline)
  └─ HimalayaDB.connect(path)                     read-only SQLite.DB (no migrations)
       └─ curated_peaks(db, exposure_id)          SELECT/UNION over auto_peaks + peak_curations
            └─ Tables.jl rows  ──(optional)──►  DataFrame(...) / HimalayaDB.dataframe(...)
       └─ reconstruct_index(db, index_id)         SELECT indices + index_peaks → Himalaya.Index{P}
       └─ load_trace(db, exposure_id)             DB lookup → .dat path → parse → (q, I)
```

## Error handling

- **Missing DB file / unreadable path:** clear error from `connect`, naming the resolved path and the `HIMALAYA_DB_PATH` env var.
- **WAL sidecar / read-only edge cases:** see the connection risk note; the plan pins behavior.
- **Unknown / dangling ids** (`exposure_id`, `index_id` not present): return empty rows for list queries; error for single-entity reconstructions (`reconstruct_index` on a missing id).
- **Trace data dir absent** (`load_trace` when only the `.db` was copied): explicit, actionable error distinguishing "no such exposure" from "data dir not reachable."
- **Phase reconstruction guard:** validate `P isa Type && P <: Himalaya.Phase` before use (mirrors HimalayaUI), erroring on an unrecognized phase name rather than dispatching on garbage.
- **Read-only enforcement:** `query_only=ON` makes accidental writes fail loudly at the SQLite layer.

## Testing

Uses stdlib `Test`, per repo convention.

1. **Fixture (`test/fixture.jl`):** build a small deterministic DB — a couple of experiments/samples/exposures, `auto_peaks`, at least one `exclude` and one `add` curation, candidate `indices`, and a confirmed `index_group`. Preference: construct it via HimalayaUI's own writers so the fixture cannot drift from the real schema (test-only dependency on HimalayaUI). The plan decides build-in-test vs. committed fixture file.
2. **Unit tests (`test/test_queries.jl`):** assert `curated_peaks` correctly applies excludes and adds (source/excluded tagging), `confirmed_indices` returns only the confirmed subset, `reconstruct_index` yields the right `Index{P}` with correct predicted peaks, `load_trace` returns matching-length q/I. Regression-floor style where against real-ish data, exact where deterministic.
3. **Contract test (`test/test_contract.jl`) — the drift guard:** run HimalayaDB's readers and HimalayaUI's equivalents (`get_peaks_for_exposure`, the confirmed-group query) against the **same fixture DB**; assert identical results. This is what makes "duplicate ~30 lines of SQL" safe: if HimalayaUI's schema or logic drifts, this test fails. Carries a test-only dependency on HimalayaUI.
4. **DataFrames extension test:** with `DataFrames` loaded, `dataframe(curated_peaks(...))` returns a well-formed `DataFrame`.

## Dependencies

- **Runtime:** `SQLite`, `DBInterface`, `Tables`, core `Himalaya` (for `reconstruct_index`).
- **Weakdep:** `DataFrames` (extension only).
- **Test-only:** `HimalayaUI` (fixture construction + contract test), `Test`.

## Deferred / future

- Extracting HimalayaUI's shared read layer into HimalayaDB (single source of truth) — revisit once the read surface is proven.
- `series` / `comparisons` curated groupings.
- `sample_tags` / `exposure_tags`.
- `user_actions` history / attribution / audit reads.
- Live-server (REST) access path.
