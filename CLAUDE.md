# Himalaya.jl — Claude working notes

## What this is

A Julia package for **indexing SAXS diffraction patterns** — given a 1D integration trace, find Bragg peaks and identify the liquid-crystalline phase (Pn3m, Im3m, Ia3d, Fm3m, Fd3m, Hexagonal, Lamellar, Square) by fitting peak q-values to known phase-ratio series.

## Code layout

Monorepo: the core `Himalaya` package lives at the root; sub-packages live under `packages/`.

```
src/                         # core Himalaya package
  Himalaya.jl                # module entry; exports public API
  peakfinding.jl             # findpeaks (persistence + sharpness + kneedle)
  persistence.jl             # topological persistence helper
  sharpness.jl               # Savitzky-Golay / CWT curvature
  threshold.jl               # kneedle elbow finder
  phase.jl                   # Phase abstract type hierarchy + phaseratios
  index.jl                   # Index struct, indexpeaks, score, fit
  util.jl
packages/
  HimalayaUI/                # web-app sub-package
    src/
      db.jl                  # SQLite schema + CRUD
      datfile.jl             # three-column .dat parser
      manifest.jl            # Google Sheets CSV → ManifestSample
      pipeline.jl            # analyze_exposure!, auto_group, persist_analysis!
      cli.jl                 # himalaya init/analyze/show/serve
      json.jl                # row → Dict serialization
      actions.jl             # X-Username extraction + user_actions logger
      routes_*.jl            # one per REST resource
      server.jl              # Oxygen.jl app + serve(db) + test harness
    frontend/dist/           # served as static at / (empty until Plan 3+)
    test/
docs/
  peak-finding.md            # narrative design notes for findpeaks
  superpowers/               # specs and plans
test/                        # core Himalaya tests
examples/                    # scripts using Himalaya (not part of the package)
scratch/                     # gitignored — exploratory scripts and trace data
```

## Running tests

```bash
# Core Himalaya
julia --project=. -e 'using Pkg; Pkg.test()'

# HimalayaUI sub-package
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'

# One test file in isolation
julia --project=. -e 'using Himalaya, Test; include("test/foo.jl")'
```

Tests use stdlib `Test` (`@testset`, `@test`, `@test_throws`). Internal (non-exported) helpers are accessed via `Himalaya.<name>` in tests.

## Running the app

```bash
# From a fresh experiment dir:
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  init /path/to/experiment --manifest manifest.csv
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  analyze /path/to/experiment
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  serve /path/to/experiment --port 8080
```

`serve` blocks. Frontend is served from `packages/HimalayaUI/frontend/dist/` if present.

## Conventions

- **TDD by default.** Failing test → minimal implementation → verify pass → commit. Each step is its own commit.
- **One responsibility per file.** When a `src/` file accumulates multiple concepts, split it.
- **Regression floors, not hard-coded counts.** For tests against real-data fixtures, prefer `recall ≥ floor` / `spurious ≤ ceiling` assertions over exact counts. Raising a floor is a deliberate commit.
- **Worktrees for feature branches.** For multi-step rewrites, `git worktree add ../Himalaya-<topic> -b <topic>`. Keeps main clean.
- **`Manifest.toml` is gitignored** — Julia library convention. Consumers re-resolve.

## HimalayaUI gotchas

**SQLite.jl:**
- `DBInterface.lastrowid` takes the query **result**, not the db: `res = DBInterface.execute(db, sql, params); id = Int(DBInterface.lastrowid(res))`.
- Raw rows from `DBInterface.execute` lose their values after the query closes. Materialize with `Tables.rowtable(DBInterface.execute(...))` to get stable `NamedTuple`s (access fields via `row.name`).

**Oxygen.jl 1.10.x:**
- Use the singleton API: `@get "/path/{id}" function(req::HTTP.Request, id::Int) ... end`. Typed function args extract path params.
- Parse JSON body with `json(req)` (unqualified, imported via `using Oxygen`), **not** `Oxygen.json(req)`.
- Test harness pattern: `Oxygen.resetstate()` before `Oxygen.serve(; async=true)`, `Oxygen.terminate()` after. A module-level `Ref{Union{SQLite.DB, Nothing}}` holds the live DB (matches the one-experiment-per-process deployment model).
- Mount static files with `Oxygen.dynamicfiles(dir, "/")` — only if `isdir(dir)`, so empty frontends don't break tests.
- Oxygen emits a harmless warning about OpenAPI schema generation for some routes; ignore it.

**Stdlib deps must be explicit.** Stdlibs used directly in a package (`Sockets`, `Printf`, `SparseArrays`, `DelimitedFiles`, etc.) must be listed in `Project.toml`'s `[deps]` — `Pkg.add` them like regular packages.

## Current state

- Core Himalaya: `v0.5.0` on `main` — v2 peak-finding (persistence + sharpness + kneedle) is merged.
- HimalayaUI: Plans 1 & 2 of 6 complete on `main` — SQLite pipeline, full REST API (Oxygen.jl), and CLI (`himalaya init/analyze/show/serve`). Plans 3–6 cover the TypeScript/Vite frontend; see [docs/superpowers/plans/](docs/superpowers/plans/).

## Further reading

- [docs/peak-finding.md](docs/peak-finding.md) — narrative design notes, non-obvious defaults, out-of-scope decisions.
- [docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md](docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md) — web app design spec (schema, API, UI layout). Load-bearing for all HimalayaUI work.
- [docs/superpowers/plans/](docs/superpowers/plans/) — implementation plans (one per sub-project).
- [docs/future-feature-ideas.md](docs/future-feature-ideas.md) — intentionally-deferred features.
