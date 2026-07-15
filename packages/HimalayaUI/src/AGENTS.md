# packages/HimalayaUI/src — Backend

Julia/Oxygen.jl REST backend, SQLite DB, event-log dispatcher, SSE multiplayer, and analysis pipeline.

## Where to look

| Task | Location | Notes |
|------|----------|-------|
| DB schema + CRUD | `db.jl` | FK enforcement ON; `AUTOINCREMENT` on mention-target tables |
| Experiment config | `config.jl` | `ExperimentConfig` + `load_config` + `resolve_files` |
| Manifest | `manifest.jl` | `ManifestSample` + `parse_manifest` (config-driven) |
| Pre-ingest validation | `validate.jl` | `ManifestViolation` + `validate_manifest` |
| Analysis pipeline | `pipeline.jl` | `analyze_exposure!`, `auto_group`, `persist_analysis!` |
| Curation contract | `pipeline.jl::effective_peaks` | `auto_peaks − excludes ∪ adds` |
| Event log + SSE | `events.jl` | `apply_event!` is the sole writer to view tables |
| Idempotency | `idempotency.jl` | `with_idempotency` + `InTransaction` sentinel |
| Speculative indices | `speculative.jl` | Build/delete + re-resolve |
| CLI | `cli.jl` | `init`, `analyze`, `show`, `serve`, `config`, `reingest` |
| REST routes | `routes_*.jl` | One file per resource |
| Audit log | `actions.jl` | `X-Username` extraction + `user_actions` writer |
| TIFF → PNG | `image.jl` | `load_and_lognormalize` (Q0f31-aware) |
| Server | `server.jl` | Oxygen app + `serve(db)` + test harness |
| JSON serialization | `json.jl` | row → Dict helpers |

## SQLite.jl

- **`DBInterface.lastrowid` takes the query *result*, not the db**: `res = DBInterface.execute(db, sql, params); id = Int(DBInterface.lastrowid(res))`.
- **Raw rows from `DBInterface.execute` lose values after the query closes.** Materialize with `Tables.rowtable(DBInterface.execute(...))` to get stable `NamedTuple`s (access via `row.name`).
- **`Tables.rowtable` returns `missing` for SQL NULL, not `nothing`.** Use `ismissing(row.field)`; normalize to `nothing` when mixing with literals: `existing = ismissing(row.field) ? nothing : row.field`. `routes_users.jl` is the canonical NULL-fill example.
- **`samples.name` is the stable identifier; `samples.display_name` is editable.** `name` is set at ingest from manifest column 2 and is never UI-mutable. `display_name` (column 3) is user-editable via PATCH; reingest never clobbers it. `UNIQUE INDEX (experiment_id, name)`. Use the `sampleDisplayName(s)` helper everywhere in the frontend.
- **FK enforcement is on.** `open_db` runs `PRAGMA foreign_keys = ON` per connection. FKs that must survive user deletion need `ON DELETE SET NULL` in the schema DDL — at the column, not at call sites. `index_groups.created_by` and `user_actions.user_id` already have this.
- **`PRAGMA schema_version = N+1` requires `writable_schema = ON`** on most SQLite builds. Used as a VACUUM-free schema-cache invalidator; silently no-ops if `writable_schema` is off. The FK-heal fallback toggles it back on for the bump.
- **PKs use `AUTOINCREMENT` on mention-targets** (`experiments`, `samples`, `exposures`, `peaks`, `indices`). Plain `INTEGER PRIMARY KEY` is rowid-aliased and reuses freed ids — chat `@`-mentions of a deleted entity could rebind to a new one. `migrate_pk_to_autoincrement!` heals legacy DBs on next `open_db`.
- **New timestamps use `comparison_now_iso()`** (T-separator, ms, Z). Legacy `CURRENT_TIMESTAMP` (space-sep) survives in `user_actions.timestamp` and `comparison_pins.pinned_at` — don't sort across the two formats as strings (issue #76).
- **`migrate_comparisons_to_series!` (#171) is the ONE migration that writes `user_actions` rows by raw INSERT — never `apply_event!`** (a migration must not broadcast N replay-as-reruns). `client_op_id`/`client_id` are NULL (so no `idempotent_responses` cache rows; the partial unique index does not apply). The sentinel marker row (`schema_migrations.name='comparisons_to_series'`, `MIGRATION_COMPARISONS_TO_SERIES`) is written LAST inside the migration's own `SQLite.transaction` — the gate flips only on a fully-committed copy. View rows are produced by folding the synthesized payloads through `update_view_for_event!`, never by direct table writes. It runs in `migrate_schema!` after `migrate_series!` (tables + sentinel exist) and the two compare migrations (reads their columns). Migrated series NULL their fork lineage (`forked_from_id`/`forked_at_hash`): `comparisons.forked_from_id` and `series.forked_from_id` are different id-spaces, so copying the id would point at an arbitrary/nonexistent series. FK enforcement is ON during this migration (`migrate_compare_relax_nullability!` re-enabled `PRAGMA foreign_keys=ON` in its `finally`, before this runs), so the dangling write would FK-throw anyway; the `PRAGMA foreign_key_check(series)` test confirms the copy leaves no dangling ref.

## Oxygen.jl 1.10.x

- Use the singleton API: `@get "/path/{id}" function(req::HTTP.Request, id::Int) ... end`. Typed function args extract path params.
- Parse JSON body with `json(req)` (unqualified, imported via `using Oxygen`), **not** `Oxygen.json(req)`.
- Test harness pattern: `Oxygen.resetstate()` before `Oxygen.serve(; async=true)`, `Oxygen.terminate()` after. A module-level `Ref{Union{SQLite.DB, Nothing}}` holds the live DB (matches the one-experiment-per-process model).
- Mount static files with `Oxygen.dynamicfiles(dir, "/")` — only if `isdir(dir)`, so empty frontends don't break tests.
- Oxygen emits a harmless warning about OpenAPI schema generation for some routes; ignore it.
- **`parallel = true` is on** (#115). Every request handler dispatches across the worker thread pool — without it, all handlers stickied to tid=1 even with `JULIA_NUM_THREADS > 1`. Routes therefore must not assume cooperative single-threaded execution. Writers serialize at the Julia level via `_DB_WRITE_LOCK` (server.jl), a `ReentrantLock` acquired by every `SQLite.transaction(db)` site on the singleton: `with_idempotency` (both branches), the default `apply_event!`, `persist_analysis!`, `analyze_exposure!`, and the `gc_idempotent_responses!` DELETE. Closes #122 Race 1 (`SQLite.transaction` TOCTOU — loud 500s + silent savepoint nesting) and Race 2 writer-vs-writer (`db.stmt_wrappers` Dict mutation). Reentrant so a route body inside `with_idempotency` calling `analyze_exposure!` doesn't self-deadlock. **Residual #122 limit:** reader-vs-writer `stmt_wrappers` mutation remains possible because reads don't take the lock. Rare; surfaces as intermittent test flakes or sporadic 500s. Per-request reader connections will close it — tracked as a follow-up.

## Pipeline + curation

- **`analyze_exposure!` curation contract.** Before calling `Himalaya.indexpeaks`, it calls `effective_peaks(db, exposure_id, q, I)` which synthesises `auto_peaks − peak_curations(kind='exclude') ∪ peak_curations(kind='add')`, with sharpness for adds sampled from `Himalaya.sharpness(I)`. Without this, exclusions don't affect scoring and manual peaks never land in `IndexEntry.peaks`. Touch only with curation-lifecycle regression tests in `test_pipeline.jl` green. See `docs/event-log.md` §1 for the auto/curation table split.
- **`persist_analysis!` is transactional.** The auto-peak diff-update + index re-resolve sequence is wrapped in `SQLite.transaction`. New write steps must go inside `_persist_analysis_inner!` so they stay atomic. Same pattern applies to `reingest!` in `cli.jl` (`_reingest_inner!`).
- **`_reingest_inner!` return shape.** `NamedTuple{(:status, :added_samples, :added_exposures, :manifest_path)}` where `:status` is `:ok` or `:no_manifest`. The `POST /api/experiments/:id/reingest` route echoes those fields in JSON (HTTP 200 in both cases — `status` is the discriminator).
- **`speculative_peak_intents` is the durable record; `index_peaks` is a resolved view.** Analysis wipes and rebuilds `index_peaks` for speculatives every run, resolving from intents (creation-time ratio→q assignments, frozen, cascade-deleted with the index). Never write intents from discovery results, and never let a failed re-resolution delete them — zero resolved peaks must leave the row `candidate` with empty `index_peaks` (do not write `status='stale'`; it's a dead enum normalized away on open).

## Experiment config

`experiment.toml` is the source of truth. Every experiment dir has one; `himalaya init` slurps it into `experiments.config` so the DB is self-contained. `analyze_exposure!` reads via `config_from_db` (in-memory parse — no tempfile per call), falling back to `simple.toml` defaults when `config IS NULL`.

- **`himalaya config new --dir <path>`** is the ONLY command that writes inside an experiment dir; refuses to overwrite.
- **`layout.exposure_type`** is validated at parse against `VALID_EXPOSURE_TYPES` (currently `("simple",)`) — extend the tuple before adding a new type.
- Malformed TOML throws a wrapped `Invalid TOML in <path>: …` error.
- **`_build_config(::AbstractDict)`** is the shared helper between `config_from_db` and `load_config`.
- **`parse_manifest` has two methods.** `parse_manifest(source)` is a backward-compat wrapper using `simple.toml` defaults. `parse_manifest(cfg::ExperimentConfig, source)` is the config-driven version — use this in new code. Both accept IO and paths via `readlines(source)`.

Read `docs/experiment-config.md` before touching `config.jl`, `manifest.jl`, or the cli init/reingest paths.

## Filesystem + DB layout

- **Read-only experiment directories at runtime.** Himalaya never creates, modifies, or deletes any file inside an experiment dir during `init`, `analyze`, `reingest`, or `serve`. Sole exception: `himalaya config new --dir`. A regression test in `test_pipeline.jl` snapshots dir contents before/after `cli_init_with_db!` — keep it green.
- **Central DB.** All CLI commands open the same DB resolved by `default_db_path()` in `db.jl`: `HIMALAYA_DB_PATH` if set, else `~/.himalaya/himalaya.db` (parent auto-created). One DB stores every experiment; experiment dirs are pure read-only data sources. Tests pass an explicit file path (`open_db(joinpath(tmp, "himalaya.db"))`) for isolation.
- **Filename ↔ exposure association via filesystem prefix scan.** Manifest filename entries are always prefixes: `JC001-004` expands to four, each scanned via `resolve_files(cfg, dir, prefix, cfg.integration_pattern)` against the filesystem. The manifest declares intent; disk decides what exists. Missing files produce a warning, not an error. When debugging "exposures missing after init/reingest," check actual files in `analysis_dir` first.

## Image route

- **Detector TIFFs are Q0f31 fixed-point.** TiffImages loads as `Gray{Q0f31}` (= `Fixed{Int32, 31}`). `Float32.(channelview(raw))` divides raw photon counts by 2³¹ (~2.1e9), making `log1p` a numerical no-op (~4.7e-10 per count). To recover photon counts, use `reinterpret.(Int32, channelview(raw))`. Then `max(., 0)` clips beamstop/dead-pixel negatives, `log1p` compresses, and a p99-of-positives clip prevents the direct beam from crushing diffraction-ring contrast. See `image.jl::load_and_lognormalize`.
- **`Cache-Control: no-store`** on `routes_exposures.jl`'s image route, paired with `cache: "no-store"` on the frontend `fetch()`, stops the browser from serving stale PNGs across analysis re-runs. Don't change to a longer max-age without invalidation tied to exposure id + analysis version.

## Exposure-state quirks

- **`exposures.selected` is sample-scoped LWW.** `PATCH /api/exposures/:id/select` clears `selected = 0` across all exposures in the sample, then sets one. Under multiplayer this is intentional — concurrent selects produce a single resolved value. **Don't add If-Match to this route.**

## Phase-type serialization

- `string(Himalaya.Pn3m)` returns the fully-qualified `"Himalaya.Pn3m"`. When storing phase names in SQLite, use `string(nameof(P))` → `"Pn3m"`.
- The inverse is `getfield(Himalaya, Symbol(name))` (always validate with `P isa Type && P <: Himalaya.Phase` before calling `phaseratios`).

## Stdlib deps

Stdlibs used directly in a package (`Sockets`, `Printf`, `SparseArrays`, `DelimitedFiles`, `TOML`, …) must be listed in `Project.toml`'s `[deps]` — `Pkg.add` them like regular packages.

## Working in worktrees

- `Manifest.toml` is gitignored, so fresh worktrees re-resolve against the registry. The Wong Lab registry publishes v0.5+, so `Pkg.instantiate()` finds a current `findpeaks` and tests pass without copying anything from main.
- To make core changes made *inside* a worktree flow through to HimalayaUI: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.develop(path="../..")'` once. The `worktree-setup` skill documents both paths.

## Anti-patterns

- **No route may INSERT/DELETE into `peak_curations` or `index_group_members` except through `apply_event!`.** It is the sole writer to view tables.
- **From inside `with_idempotency`, call `apply_event!(InTransaction(), …)`** — never the public `apply_event!(db, req; …)`. The public method opens a nested savepoint AND broadcasts immediately, bypassing the post-commit queue. See `docs/event-log.md` §3a + `docs/mutation-queue.md`.
- **Experiment dirs are read-only at runtime** (except `himalaya config new`).
