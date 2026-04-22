# Himalaya.jl — Claude working notes

## What this is

A Julia package for **indexing SAXS diffraction patterns** — given a 1D integration trace, find Bragg peaks and identify the liquid-crystalline phase (Pn3m, Im3m, Ia3d, Fm3m, Fd3m, Hexagonal, Lamellar, Square) by fitting peak q-values to known phase-ratio series.

## Code layout

```
src/
  Himalaya.jl       # module entry; exports public API
  peakfinding.jl    # findpeaks — being rewritten; see `peakfinding-rewrite` branch
  phase.jl          # Phase abstract type hierarchy + phaseratios
  index.jl          # Index struct, indexpeaks, score, fit
  util.jl           # small helpers
docs/
  peak-finding.md   # narrative design notes for findpeaks
  superpowers/      # specs and plans
test/
  runtests.jl       # orchestrator
  index.jl          # indexing tests
examples/           # scripts using Himalaya (not part of the package)
scratch/            # gitignored — exploratory scripts and trace data
```

## Running tests

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'                       # first time
julia --project=. -e 'using Pkg; Pkg.test()'                              # full suite
julia --project=. -e 'using Himalaya, Test; include("test/foo.jl")'       # one file
```

Tests use stdlib `Test` (`@testset`, `@test`, `@test_throws`). Internal (non-exported) helpers are accessed via `Himalaya.<name>` in tests.

## Conventions

- **TDD by default.** Failing test → minimal implementation → verify pass → commit. Each step is its own commit.
- **One responsibility per file.** When a `src/` file accumulates multiple concepts, split it.
- **Regression floors, not hard-coded counts.** For tests against real-data fixtures, prefer `recall ≥ floor` / `spurious ≤ ceiling` assertions over exact counts. Raising a floor is a deliberate commit.
- **Worktrees for feature branches.** For multi-step rewrites, `git worktree add ../Himalaya-<topic> -b <topic>`. Keeps main clean.
- **`Manifest.toml` is gitignored** — Julia library convention. Consumers re-resolve.

## Current state

`v0.4.5` on `main`. A ground-up peak-finding rewrite (persistent homology + sharpness + kneedle thresholds) lives on the `peakfinding-rewrite` branch, not yet merged.

## Further reading

- [docs/peak-finding.md](docs/peak-finding.md) — narrative design notes, non-obvious defaults, out-of-scope decisions.
- [docs/superpowers/specs/](docs/superpowers/specs/) — formal design specifications.
- [docs/superpowers/plans/](docs/superpowers/plans/) — implementation plans.
