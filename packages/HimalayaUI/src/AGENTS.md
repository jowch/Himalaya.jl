# packages/HimalayaUI/src — Backend

Julia/Oxygen.jl REST backend, SQLite DB, event-log dispatcher, SSE multiplayer, and analysis pipeline.

## Where to look

| Task | Location | Notes |
|------|----------|-------|
| DB schema + CRUD | `db.jl` | FK enforcement ON; `AUTOINCREMENT` on mention-target tables |
| Experiment config | `config.jl` | `ExperimentConfig` + `load_config` + `resolve_files` |
| Ingestion (auto-scan) | `ingest.jl` | `scan_and_group!` — sole ingestion entry point (HTTP scan); read-only wrt experiment dir |
| Analysis pipeline | `pipeline.jl` | `analyze_exposure!`, `auto_group`, `persist_analysis!` |
| Curation contract | `pipeline.jl::effective_peaks` | `auto_peaks − excludes ∪ adds` |
| Event log + SSE | `events.jl` | `apply_event!` is the sole writer to view tables |
| Idempotency | `idempotency.jl` | `with_idempotency` + `InTransaction` sentinel |
| Speculative indices | `speculative.jl` | Build/delete + re-resolve |
| CLI | `cli.jl` | `analyze`, `show`, `serve`, `config`, `migrate-toml` (ingestion is HTTP-scan-only, no `init`/`reingest`) |
| REST routes | `routes_*.jl` | One file per resource |
| Audit log | `actions.jl` | `X-Username` extraction + `user_actions` writer |
| TIFF → PNG | `image.jl` | `load_and_lognormalize` (Q0f31-aware) |
| Server | `server.jl` | Oxygen app + `serve(db)` + test harness |
| JSON serialization | `json.jl` | row → Dict helpers |

## SQLite.jl

- **`DBInterface.lastrowid` takes the query *result*, not the db**: `res = DBInterface.execute(db, sql, params); id = Int(DBInterface.lastrowid(res))`.
- **Raw rows from `DBInterface.execute` lose values after the query closes.** Materialize with `Tables.rowtable(DBInterface.execute(...))` to get stable `NamedTuple`s (access via `row.name`).
- **`Tables.rowtable` returns `missing` for SQL NULL, not `nothing`.** Use `ismissing(row.field)`; normalize to `nothing` when mixing with literals: `existing = ismissing(row.field) ? nothing : row.field`. `routes_users.jl` is the canonical NULL-fill example.
- **`samples.name` is the stable identifier; `samples.display_name` is editable.** `name` is set at ingest from manifest column 2 and is never UI-mutable. `display_name` (column 3) is user-editable via PATCH, **but reingest refreshes it from the manifest** (`_reingest_inner!` runs `UPDATE samples SET display_name = ?, notes = ?`, `cli.jl`) — overwriting UI edits unless the manifest carries the same value (see `docs/experiment-config.md`). `UNIQUE INDEX (experiment_id, name)`. Use the `sampleDisplayName(s)` helper everywhere in the frontend.
- **`exposure_sources.role` encodes two derived-exposure use-cases in one table** (`db.jl`): `role='signal'` alone = simple averaging; `role='signal'` + `role='background'` together = background subtraction. The single text-enum column is deliberate — a third provenance type extends it without a schema migration. (Derived-exposure construction itself is still a future feature.)
- **FK enforcement is on.** `open_db` runs `PRAGMA foreign_keys = ON` per connection. FKs that must survive user deletion need `ON DELETE SET NULL` in the schema DDL — at the column, not at call sites. `index_groups.created_by` and `user_actions.user_id` already have this.
- **`PRAGMA schema_version = N+1` requires `writable_schema = ON`** on most SQLite builds. Used as a VACUUM-free schema-cache invalidator; silently no-ops if `writable_schema` is off. The FK-heal fallback toggles it back on for the bump.
- **PKs use `AUTOINCREMENT` on mention-targets** (`experiments`, `samples`, `exposures`, `auto_peaks`, `peak_curations`, `indices`). Plain `INTEGER PRIMARY KEY` is rowid-aliased and reuses freed ids — chat `@`-mentions of a deleted entity could rebind to a new one. `migrate_pk_to_autoincrement!` heals legacy DBs on next `open_db`. (The literal `_AUTOINCREMENT_ENTITIES` tuple still lists the removed `peaks` table purely for legacy-DB healing.)
- **New timestamps use `comparison_now_iso()`** (T-separator, ms, Z). Legacy `CURRENT_TIMESTAMP` (space-sep) survives in `user_actions.timestamp` and `comparison_pins.pinned_at` — don't sort across the two formats as strings (issue #76).
- **`open_db` enables WAL + finalizes statements twice.** Every file-based DB gets `PRAGMA journal_mode = WAL` (skipped for `:memory:` and `file:` URIs) — a one-time persistent file-header change enabling concurrent readers + one writer, load-bearing for `parallel = true` request handling (#115). `open_db` calls `SQLite.finalize_statements!` twice: once after the DDL migrations (`migrate_schema!`), once after the WAL PRAGMA. Without them a cached prepared statement from a migration transaction can reference a dropped `_migrate_old_*` table and fail the first real query. Test harnesses that bypass `open_db` can hit intermittent "database is locked" / "no such table: main._migrate_old_*". See `db.jl::open_db`.
- **`_MIGRATION_TEMP_ENTITIES` carries a substring-uniqueness invariant.** `_fix_fk_references_after_autoincrement_migration!` heals FK refs via an UNANCHORED substring `REPLACE` on `sqlite_master.sql` for each `_migrate_old_<entity>` name in the tuple — safe ONLY because no entity name is a substring of another (adding `peak` while `peaks` exists would corrupt every `peaks` ref). The constant (`= (_AUTOINCREMENT_ENTITIES..., "index_peaks")`) must stay in lockstep with `migrate_pk_to_autoincrement!` and `migrate_r2_widen_index_peaks_pk!`; if you add an entity, re-check the invariant or anchor the replace. See `db.jl` (constant + heal loop).
- **`migrate_comparisons_to_series!` (#171) is the ONE migration that writes `user_actions` rows by raw INSERT — never `apply_event!`** (a migration must not broadcast N replay-as-reruns). `client_op_id`/`client_id` are NULL (so no `idempotent_responses` cache rows; the partial unique index does not apply). The sentinel marker row (`schema_migrations.name='comparisons_to_series'`, `MIGRATION_COMPARISONS_TO_SERIES`) is written LAST inside the migration's own `SQLite.transaction` — the gate flips only on a fully-committed copy. View rows are produced by folding the synthesized payloads through `update_view_for_event!`, never by direct table writes. It runs in `migrate_schema!` after `migrate_series!` (tables + sentinel exist) and the two compare migrations (reads their columns). Migrated series NULL their fork lineage (`forked_from_id`/`forked_at_hash`): `comparisons.forked_from_id` and `series.forked_from_id` are different id-spaces, so copying the id would point at an arbitrary/nonexistent series. FK enforcement is ON during this migration (`migrate_compare_relax_nullability!` re-enabled `PRAGMA foreign_keys=ON` in its `finally`, before this runs), so the dangling write would FK-throw anyway; the `PRAGMA foreign_key_check(series)` test confirms the copy leaves no dangling ref.

## Oxygen.jl 1.10.x

- Use the singleton API: `@get "/path/{id}" function(req::HTTP.Request, id::Int) ... end`. Typed function args extract path params.
- Parse JSON body with `json(req)` (unqualified, imported via `using Oxygen`), **not** `Oxygen.json(req)`.
- Test harness pattern: `Oxygen.resetstate()` before `Oxygen.serve(; async=true)`, `Oxygen.terminate()` after. A module-level `Ref{Union{SQLite.DB, Nothing}}` holds the live DB (matches the one-experiment-per-process model).
- Mount static files with `Oxygen.dynamicfiles(dir, "/")` — only if `isdir(dir)`, so empty frontends don't break tests.
- Oxygen emits a harmless warning about OpenAPI schema generation for some routes; ignore it.
- **Custom response headers are dropped on idempotency cache replay.** A route wrapped in `with_idempotency` that sets custom headers (`Cache-Control`, CORS, `X-Image-*`, …) replays only `Content-Type: application/json` on a cache HIT — the `idempotent_responses` row stores `status_code` + `body`, not headers. Silent, not an error. A new idempotent route needing custom headers must persist them in the cache schema or accept they vanish on retry. See `idempotency.jl::_lookup_cached_response`.
- **`parallel = true` is on** (#115). Every request handler dispatches across the worker thread pool — without it, all handlers stickied to tid=1 even with `JULIA_NUM_THREADS > 1`. Routes therefore must not assume cooperative single-threaded execution. Writers serialize at the Julia level via `_DB_WRITE_LOCK` (server.jl), a `ReentrantLock` acquired by every `SQLite.transaction(db)` site on the singleton: `with_idempotency` (both branches), the default `apply_event!`, `persist_analysis!`, `analyze_exposure!`, and the `gc_idempotent_responses!` DELETE. Closes #122 Race 1 (`SQLite.transaction` TOCTOU — loud 500s + silent savepoint nesting) and Race 2 writer-vs-writer (`db.stmt_wrappers` Dict mutation). Reentrant so a route body inside `with_idempotency` calling `analyze_exposure!` doesn't self-deadlock. **Residual #122 limit:** reader-vs-writer `stmt_wrappers` mutation remains possible because reads don't take the lock. Rare; surfaces as intermittent test flakes or sporadic 500s. Per-request reader connections will close it — tracked as a follow-up.

## Pipeline + curation

- **`analyze_exposure!` curation contract.** Before calling `Himalaya.indexpeaks`, it calls `effective_peaks(db, exposure_id, q, I)` which synthesises `auto_peaks − peak_curations(kind='exclude') ∪ peak_curations(kind='add')`, with sharpness for adds sampled from `Himalaya.sharpness(I)`. Without this, exclusions don't affect scoring and manual peaks never land in `IndexEntry.peaks`. Touch only with curation-lifecycle regression tests in `test_pipeline.jl` green. See `docs/event-log.md` §1 for the auto/curation table split.
- **`persist_analysis!` is transactional.** The auto-peak diff-update + index re-resolve sequence is wrapped in `SQLite.transaction`. New write steps must go inside `_persist_analysis_inner!` so they stay atomic. The ingestion path (`scan_and_group!` in `ingest.jl`) follows the same inner-transaction pattern.
- **`index_peaks.ratio_position` is the Julia↔DB round-trip contract for the index data model.** A core `Index{P}.peaks` is a `SparseVector` (`src/index.jl`); the slot index *is* the ratio-series position. Persist via `SparseArrays.findnz(idx.peaks)` → `(ratio_positions, peak_qvals)` written to `index_peaks` (`pipeline.jl`); the read path rebuilds the `SparseVector` from sorted `ratio_position`s. Any migration or rewrite of the persist/read path must preserve this slot identity or reconstruction breaks silently.
- **Route-local `_ngc_for_phase` intentionally diverges from `Himalaya.ngc` for unit honesty.** `routes_analysis.jl::_ngc_for_phase` uses the bare Gauss–Bonnet form (`κ = -2π·(χ/A₀)/a²`) so κ carries the same length unit as the experiment's lattice parameter (`Å⁻²` for q in `A-1`, `nm⁻²` for q in `nm-1`). `Himalaya.ngc` hard-codes a `(10/a)²` factor assuming `a` in Å and emitting `nm⁻²` — wrong for `nm-1`-q experiments. Don't "DRY" the route back onto `Himalaya.ngc`.
- **`speculative_peak_intents` is the durable record; `index_peaks` is a resolved view.** Analysis wipes and rebuilds `index_peaks` for speculatives every run, resolving from intents (creation-time ratio→q assignments, frozen, cascade-deleted with the index). Never write intents from discovery results **unless the index is `basis_locked`** (see below), and never let a failed re-resolution delete them — zero resolved peaks must leave the row `candidate` with empty `index_peaks` (do not write `status='stale'`; it's a dead enum normalized away on open).
- **`indices.basis_locked = 1` is what makes scan-derived intents legitimate.** The no-discovery-intents rule above exists for one reason: resolved intents feed `_persist_analysis_inner!`'s least-squares refit (`pipeline.jl`), so intents invented by a tolerance scan can drag an index's lattice somewhere the user never put it. `insert_custom_index!` (the custom-index modal) writes exactly such scan-derived intents — legitimately, because it also sets `basis_locked = 1`, and the pipeline skips the refit for locked rows: it re-resolves their peaks and recomputes score/R²/residuals **against the stored basis**, never solving for a new one. The invariant is therefore *intents may be scan-derived iff they cannot move the comb*. If you add another writer of scan-derived intents, it must lock the basis too — or the guarantee that a user-chosen lattice survives reanalysis silently dies.

  **The pipeline also skips its auto-discovery pass entirely for locked rows**, and that skip — not any persisted column — is what freezes a locked index's claim set. The commit bounds its claim to the reflections the modal *drew* (a ratio set, not a count: see `compute_snap` on Hexagonal); discovery scans the whole core series, so re-enabling it for locked rows would let a reanalysis claim orders the user was never shown and make the rail's "explains N peaks" exceed the modal's "N of M land" again. It would also search against a comb refit from the resolved peaks — the very comb the locked row refuses to store. Re-resolution from intents can only ever reproduce a subset of what was drawn; nothing may re-widen it.

## Experiment config

`experiment.toml` is the source of truth. Every experiment dir has one; `himalaya init` slurps it into `experiments.config` so the DB is self-contained. `analyze_exposure!` reads via `config_from_db` (in-memory parse — no tempfile per call), falling back to `simple.toml` defaults when `config IS NULL`.

- **`himalaya config new --dir <path>`** is the ONLY command that writes inside an experiment dir; refuses to overwrite.
- **`layout.exposure_type`** is validated at parse against `VALID_EXPOSURE_TYPES` (currently `("simple",)`) — extend the tuple before adding a new type.
- Malformed TOML throws a wrapped `Invalid TOML in <path>: …` error.
- **`_build_config(::AbstractDict)`** is the shared helper between `config_from_db` and `load_config`.

Read `docs/experiment-config.md` before touching `config.jl` or the ingestion (`ingest.jl` / HTTP scan) path.

## Filesystem + DB layout

- **Read-only experiment directories at runtime.** Himalaya never creates, modifies, or deletes any file inside an experiment dir during ingestion (`scan_and_group!` / HTTP scan), `analyze`, or `serve`. Sole exception: `himalaya config new --dir`.
- **Central DB.** All CLI commands open the same DB resolved by `default_db_path()` in `db.jl`: `HIMALAYA_DB_PATH` if set, else `~/.himalaya/himalaya.db` (parent auto-created). One DB stores every experiment; experiment dirs are pure read-only data sources. Tests pass an explicit file path (`open_db(joinpath(tmp, "himalaya.db"))`) for isolation.
- **Filename ↔ exposure association via filesystem prefix scan.** Manifest filename entries are always prefixes: `JC001-004` expands to four, each scanned via `resolve_files(cfg, dir, prefix, cfg.integration_pattern)` against the filesystem. The manifest declares intent; disk decides what exists. Missing files produce a warning, not an error. When debugging "exposures missing after init/reingest," check actual files in `analysis_dir` first.

## Image route

- **Detector TIFFs are Q0f31 fixed-point.** TiffImages loads as `Gray{Q0f31}` (= `Fixed{Int32, 31}`). `Float32.(channelview(raw))` divides raw photon counts by 2³¹ (~2.1e9), making `log1p` a numerical no-op (~4.7e-10 per count). To recover photon counts, use `reinterpret.(Int32, channelview(raw))`. Then `max(., 0)` clips beamstop/dead-pixel negatives, `log1p` compresses, and a p99-of-positives clip prevents the direct beam from crushing diffraction-ring contrast. See `image.jl::load_and_lognormalize`.
- **`Cache-Control: private, max-age=31536000, immutable`** on `routes_exposures.jl`'s image route. The frontend appends `?v=<image_version_token>` (= `IMAGE_PROCESSING_VERSION` + TIFF mtime), so the URL itself is the cache key — a TIFF rewrite or a `IMAGE_PROCESSING_VERSION` bump changes the URL and forces a refetch. Bump the version const (image.jl) whenever rendered bytes change; never lengthen the max-age without keeping the `?v=` invalidation intact.
- **Thumb disk cache (issue #261).** `?thumb=1` requests go through `ensure_thumb_cached(db, exposure_id, path)` — a `<db_dir>/cache/thumb-128/{id}-{token}.png` cache, where `token = image_version_token(path)` so a version bump or TIFF rewrite naturally re-keys (stale entries are left, never read). Guarded on `db.file == ":memory:"` (renders fresh, no write) so in-memory test DBs work unchanged. The thumb path uses `load_and_lognormalize_thumb` (downscale-BEFORE-lognormalize) — a triage preview, NOT the science view; a faint narrow ring may under-render. `prewarm_thumbnails!(db; threads, overwrite, experiment_id)` populates the cache from `scan_and_group!` (after the ingest tx + analyze, gated on the `analyze` flag so test runs skip it); workers are FS-only — never touch the DB connection inside the `@threads` loop, and a per-row `try`/`catch` keeps one undecodable TIFF from aborting the whole `@threads` batch. **The ingest call is scoped `experiment_id=` and passes `overwrite=false`.** Unscoped it walked every exposure in every experiment on every scan (O(corpus) per ingest, and a long silent tail with the progress bar pinned at 100%); scoping by *experiment* rather than by the ids just inserted is deliberate — it preserves the repair path, so a thumbnail that failed to warm gets another attempt next scan. The tradeoff: `overwrite=true` (which defeats whole-second mtime granularity on a re-scan) costs a full re-render of every exposure every time, so the ingest path gives it up and leaves only the same-second-rewrite hole open. The migration path (`regroup_experiment!`) intentionally does NOT prewarm: it is a one-time operator step where lazy first-view generation is acceptable, and prewarming a full corpus over SMB would add minutes.

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

- **An evicted SSE subscriber's channel MUST be closed, not just dropped from `SSE_SUBSCRIBERS`.** All fan-out goes through `_fanout_frame!` (events.jl); `_try_put!` returning false means the subscriber is closed or saturated, and `_fanout_frame!` then both un-registers it *and* `close(sub.pending)`. Dropping it alone leaves the `/api/events` handler parked on `for frame in pending` (server.jl) forever — the HTTP response never finishes, so the browser's EventSource sees a healthy socket, never fires `onerror`, and never auto-reconnects; the client goes permanently deaf to *every* curation event until a reload. `docs/event-log.md` §"SSE" already specified this behavior ("a pruned subscriber reconnects on the EventSource side and gets a fresh subscription") — the doc was right and the code was wrong. An ingest scan is the reliable trigger (it out-runs the 64-slot channel). Any new broadcaster must call `_fanout_frame!` rather than hand-rolling the subscriber loop.

- **No route may INSERT/DELETE into `peak_curations` or `index_group_members` except through `apply_event!`.** It is the sole writer to view tables.
- **From inside `with_idempotency`, call `apply_event!(InTransaction(), …)`** — never the public `apply_event!(db, req; …)`. The public method opens a nested savepoint AND broadcasts immediately, bypassing the post-commit queue. See `docs/event-log.md` §3a + `docs/mutation-queue.md`.
- **Experiment dirs are read-only at runtime** (except `himalaya config new`).
