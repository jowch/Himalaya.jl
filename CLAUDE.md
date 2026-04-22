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
  HimalayaUI/                # web-app sub-package: SQLite, pipeline, CLI, REST (planned)
    src/{db,datfile,manifest,pipeline,cli}.jl
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

## Conventions

- **TDD by default.** Failing test → minimal implementation → verify pass → commit. Each step is its own commit.
- **One responsibility per file.** When a `src/` file accumulates multiple concepts, split it.
- **Regression floors, not hard-coded counts.** For tests against real-data fixtures, prefer `recall ≥ floor` / `spurious ≤ ceiling` assertions over exact counts. Raising a floor is a deliberate commit.
- **Worktrees for feature branches.** For multi-step rewrites, `git worktree add ../Himalaya-<topic> -b <topic>`. Keeps main clean.
- **`Manifest.toml` is gitignored** — Julia library convention. Consumers re-resolve.

## HimalayaUI gotchas (SQLite.jl)

- `DBInterface.lastrowid` takes the query **result**, not the db: `res = DBInterface.execute(db, sql, params); id = Int(DBInterface.lastrowid(res))`.
- Raw rows from `DBInterface.execute` lose their values after the query closes. Materialize with `Tables.rowtable(DBInterface.execute(...))` to get stable `NamedTuple`s (access fields via `row.name`).

## Current state

- Core Himalaya: `v0.5.0` on `main` — v2 peak-finding (persistence + sharpness + kneedle) is merged.
- HimalayaUI: Plan 1 of 6 (Foundation) complete on `main` — SQLite-backed pipeline and CLI (`himalaya init / analyze / show`). Plans 2–6 cover the REST API and frontend; see [docs/superpowers/plans/](docs/superpowers/plans/).

## Further reading

- [docs/peak-finding.md](docs/peak-finding.md) — narrative design notes, non-obvious defaults, out-of-scope decisions.
- [docs/superpowers/specs/](docs/superpowers/specs/) — formal design specifications.
- [docs/superpowers/plans/](docs/superpowers/plans/) — implementation plans.
