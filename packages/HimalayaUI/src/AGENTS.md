# packages/HimalayaUI/src — Backend

## OVERVIEW
Julia/Oxygen.jl REST backend, SQLite DB, event-log dispatcher, SSE multiplayer, and analysis pipeline.

## WHERE TO LOOK
| Task | Location | Notes |
|------|----------|-------|
| DB schema + CRUD | `db.jl` | FK enforcement ON; `AUTOINCREMENT` on mention-target tables |
| Experiment config | `config.jl` | ExperimentConfig + load_config + resolve_files |
| Manifest | `manifest.jl` | ManifestSample + parse_manifest (config-driven) |
| Analysis pipeline | `pipeline.jl` | analyze_exposure!, auto_group, persist_analysis! |
| Event log + SSE | `events.jl` | `apply_event!` is sole writer to view tables |
| Idempotency | `idempotency.jl` | Stripe-style idempotency keys for queue mutations |
| CLI | `cli.jl` | Subcommands: init, analyze, show, serve, config, reingest |
| REST routes | `routes_*.jl` | One file per resource (exposures, peaks, indices, etc.) |
| Server | `server.jl` | Oxygen.jl static asset mount + route registration |

## CONVENTIONS
- `DBInterface.lastrowid(result)` — takes the **result**, not the db
- `Tables.rowtable(DBInterface.execute(...))` — materialize rows before query closes
- SQL NULL → `missing` (not `nothing`): use `ismissing(row.field)`
- `apply_event!(InTransaction(), ...)` inside `with_idempotency` — never the public `apply_event!(db, req; ...)` from within a queue mutation
- `persist_analysis!` and `reingest!` are transactional — add new writes inside `_persist_analysis_inner!` / `_reingest_inner!`
- New timestamps use `comparison_now_iso()` (T-separator, ms, Z); legacy `CURRENT_TIMESTAMP` (space-sep) survives in `user_actions.timestamp` and `comparison_pins.pinned_at`. Don't sort across the two formats as strings (issue #76).

## ANTI-PATTERNS
- No route may INSERT/DELETE into `peak_curations` or `index_group_members` except through `apply_event!`
- Experiment dirs are read-only at runtime (except `himalaya config new`)
