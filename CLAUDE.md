# Himalaya.jl — Claude working notes

## What this is

A Julia monorepo for **indexing SAXS diffraction patterns**. The core `Himalaya` package finds Bragg peaks in a 1D integration trace and identifies the liquid-crystalline phase (Pn3m, Im3m, Ia3d, Fm3m, Fd3m, Hexagonal, Lamellar, Square) by fitting peak q-values to known phase-ratio series. `HimalayaUI` (under `packages/`) is a full-stack web app — Julia/Oxygen.jl REST backend + React/Vite frontend — a thinking tool for forming and testing phase hypotheses about SAXS diffraction patterns, both deep within a single pattern and across several.

## How to use this file

This file is the **top-level index**. Module-specific gotchas live in nested `AGENTS.md` files — always read the AGENTS.md closest to the code you're touching. Design docs live in `docs/`.

```
CLAUDE.md                                            # this file (index)
AGENTS.md                                            # quick reference (you may already have read it)
src/AGENTS.md                                        # core peak-finding/indexing
packages/HimalayaUI/src/AGENTS.md                    # backend (SQLite, Oxygen, pipeline)
packages/HimalayaUI/test/AGENTS.md                   # Julia test patterns
packages/HimalayaUI/frontend/src/AGENTS.md           # frontend conventions
packages/HimalayaUI/frontend/src/shell/AGENTS.md  # app shell / routing / modal gotchas
packages/HimalayaUI/frontend/src/lib/queue/AGENTS.md  # mutation queue internals
packages/HimalayaUI/frontend/test/AGENTS.md          # Vitest/JSDOM patterns
packages/HimalayaUI/frontend/e2e/AGENTS.md           # Playwright patterns
```

## Read first

If this is your first session on this repo, skim these in order before touching code:

1. [docs/peak-finding.md](docs/peak-finding.md) — why `findpeaks` is the way it is. Load-bearing.
2. [docs/experiment-config.md](docs/experiment-config.md) — required if touching the New-experiment funnel (`routes_fs.jl`), experiment create/scan (`routes_experiments.jl`, `ingest.jl`), or `config.jl`.
3. [docs/scoring.md](docs/scoring.md) — required if touching `score`, `auto_group`, or `remove_subsets`.
4. [docs/event-log.md](docs/event-log.md) — required if touching `events.jl`, `hash.jl`, the `apply_event!` call sites, or the SSE handler.
5. [docs/mutation-queue.md](docs/mutation-queue.md) — required if touching `lib/queue/`, `idempotency.jl`, `with_idempotency`, or `applyRemoteToCache.ts`.
6. [docs/contract-testing.md](docs/contract-testing.md) — six-layer testing rule. Required before fixing a queue/SSE/cache reconciliation bug.

## Code layout

```
src/                                # core Himalaya package — see src/AGENTS.md
test/                               # core tests
packages/HimalayaPeakModel/         # peak candidates, features, labels, swappable learner/feature-set/recipe, CV (depends on Himalaya, DecisionTree, SHA)
tools/peak-classifier/              # CLI over HimalayaPeakModel + HimalayaUI.labeled_exposures
packages/HimalayaDB/                # read-only Julia reader of the central DB for analysis scripts; contract-tested against HimalayaUI
packages/HimalayaUI/
  src/                              # backend — see packages/HimalayaUI/src/AGENTS.md
  test/                             # backend tests — see packages/HimalayaUI/test/AGENTS.md
  frontend/
    src/                            # frontend — see frontend/src/AGENTS.md
      ui/                           # every appearance-owning component (primitives and workspace components)
      shell/                        # app shell, routing, modals — see shell/AGENTS.md
      pages/ components/            # experiments home, New-experiment funnel, their composites
      plot/ detector/               # render layers
      workspace/                    # the three-pane triage workspace: WorkspacePage and its logic
      hooks/                        # custom React hooks
      lib/                          # url, plot, peak-ref, custom-index, label and bridge helpers
        queue/                      # mutation queue — see lib/queue/AGENTS.md
    test/                           # Vitest unit tests — see frontend/test/AGENTS.md
    e2e/                            # Playwright mocked — see e2e/AGENTS.md
    dist/                           # vite build output; served by Oxygen.jl in prod
docs/                               # design docs (peak-finding, scoring, event-log, …)
examples/                           # scripts using Himalaya (not part of the package)
scratch/                            # gitignored — exploratory scripts and trace data
.claude/
  skills/                           # review-pr, worktree-setup, new-route, new-event-kind, …
  agents/                           # frontend-reviewer, himalaya-reviewer, queue-reviewer, saxs-physics-reviewer
  settings.json                     # hooks: Vitest --run flag, pre-tool-use guards
```

## Running tests

```bash
# First-time setup (also run after `git worktree add`). `make dev-setup` runs the same
# steps (with `npm ci`) and undoes Pkg.develop's rewrite of HimalayaUI/Project.toml:
# HimalayaUI must `develop` the local core BEFORE instantiating — Manifest.toml is
# gitignored and there is no [sources], so a bare instantiate resolves Himalaya from
# the REGISTRY (the published core, not this checkout). `[compat] Himalaya = "0.6"`
# only catches a mismatch across VERSIONS: an unbumped core edit resolves cleanly
# and silently disagrees with this checkout. `develop` is the fix, not a workaround.
# For development use `make dev-setup`: it first writes packages/HimalayaUI/LocalPreferences.toml
# (gitignored) with `precompile_workload = false`, so a src/ edit recompiles in ~5 s, not
# ~100 s. `make setup` is the deploy path and keeps the route warm-up. Delete that file (or
# set it to true) to build the warm-up in dev, e.g. to time a cold start. The manual steps
# below don't write it. See packages/HimalayaUI/src/AGENTS.md.
(cd packages/HimalayaUI/frontend && npm install)
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

# Core Himalaya
julia --project=. -e 'using Pkg; Pkg.test()'

# HimalayaDB (read-only DB reader; contract tests run HimalayaUI's own readers on its fixture)
julia --project=packages/HimalayaDB -e 'using Pkg; Pkg.test()'

# HimalayaUI backend (Julia) — ~3 min serial. Capture once. See packages/HimalayaUI/test/AGENTS.md
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
# Faster (~2 min): run the 7 GROUP buckets as parallel processes:
make test-parallel
# One bucket only (db|pipeline|routes|routes2|events|wire|migration):
GROUP=routes julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'

# HimalayaUI frontend — from packages/HimalayaUI/frontend/
npm test              # Vitest unit tests (one-shot)
npm run test:watch    # Vitest watch mode
npm run e2e           # Playwright E2E (mocked; builds and serves the bundle; runs in CI)
npm run build         # lint:design + lint:ratchet + eslint + tsc --noEmit + vite build (must pass before PR)
npm run typecheck:test  # tsc over src + test + e2e (Vitest and Playwright skip types; CI runs this)
```

Tests use stdlib `Test` (`@testset`, `@test`, `@test_throws`). Internal (non-exported) helpers are accessed via `Himalaya.<name>` / `HimalayaUI.<name>` in tests.

**The Julia backend suite runs ~3 min serial (~2 min via `make test-parallel`).** Still capture output once and grep the file rather than re-running with different filters. Same for `npm test`. Route tests dispatch in-process (`with_inproc_routes`), not over a socket — see `packages/HimalayaUI/test/AGENTS.md` for the dispatch pattern, GROUP buckets, and the wire-keeper boundary.

## Running the app

```bash
# Fast path via compiled sysimage (build once, ~5 min):
make sysimage          # creates build/himalaya.so
# serve has no path arg — it serves the central DB (HIMALAYA_DB_PATH / default_db_path),
# creating an empty one if none exists. Create experiments in the UI (New experiment).
bin/himalaya serve --port 8080
bin/himalaya analyze -e <id|name|path>       # or --all; show takes -e too
bin/himalaya merge <exp> <exp>... [--apply]   # one-time: fold per-setup experiments into one beamtime (dry run without --apply)

# Without sysimage (slower cold start, no build required):
julia --threads=4,1 --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  serve --port 8080
```

`serve` blocks. Frontend is served from `packages/HimalayaUI/frontend/dist/` if present.

**Env vars** (see `packages/HimalayaUI/.env.example`): `HIMALAYA_DB_PATH`, `HIMALAYA_HOST`, `HIMALAYA_PORT`, `HIMALAYA_FRONTEND_DIST`.

Frontend dev loop (live backend, port cleanup, side-by-side prod-DB sessions): [docs/frontend-dev-loop.md](docs/frontend-dev-loop.md).

**No lab data handy (or in a cloud session)?** `make dev-setup && make demo` generates a synthetic SSRL-style experiment into gitignored `scratch/`, ingests it and starts backend + Vite; `make demo-empty` starts on an empty DB to walk the ingest funnel. Never commit the generated data. See [docs/demo.md](docs/demo.md), including the cloud-environment requirements.

**Verifying a change in the real app:** use the `verify-himalaya` skill. `bin/verify doctor`, `bin/verify seed`, then drive the browser (`open`, `snapshot`, `click`, `press`, `screenshot`) and run `bin/verify check-truth` against the demo's `TRUTH.csv`. Its Feature Map (`.claude/skills/verify-himalaya/references/features/README.md`) says how to reach and drive each feature. CLI tests: `npm run test:verify` in `packages/HimalayaUI/frontend`.

## Conventions

- **TDD by default.** Failing test → minimal implementation → verify pass → commit. Each step is its own commit.
- **One responsibility per file.** When a `src/` file accumulates multiple concepts, split it.
- **Regression floors, not hard-coded counts.** For tests against real-data fixtures, prefer `recall ≥ floor` / `spurious ≤ ceiling` assertions over exact counts. Raising a floor is a deliberate commit.
- **Worktrees for feature branches.** For multi-step rewrites, `git worktree add .claude/worktrees/<topic> -b <topic>`. `.claude/worktrees/` is gitignored.
- **`Manifest.toml` is gitignored** — Julia library convention. Consumers re-resolve.

## Where to find gotchas

Module-specific conventions and anti-patterns live in the AGENTS.md file nearest the code you're editing. The deeper the file, the more specific the rules:

| Editing… | Read… |
|---|---|
| Peak finding, scoring, phase types | [src/AGENTS.md](src/AGENTS.md) |
| Routes, SQLite, pipeline, events, image, idempotency | [packages/HimalayaUI/src/AGENTS.md](packages/HimalayaUI/src/AGENTS.md) |
| Backend tests, in-process SSE, FK-heal fixtures | [packages/HimalayaUI/test/AGENTS.md](packages/HimalayaUI/test/AGENTS.md) |
| Zustand, TanStack Query, SSE wiring, Tailwind | [packages/HimalayaUI/frontend/src/AGENTS.md](packages/HimalayaUI/frontend/src/AGENTS.md) |
| App shell / routing / modal gotchas (AppRoutes, ExperimentEntry, ShellPalette, OnboardingFlow) | [packages/HimalayaUI/frontend/src/shell/AGENTS.md](packages/HimalayaUI/frontend/src/shell/AGENTS.md) |
| Mutation queue internals | [packages/HimalayaUI/frontend/src/lib/queue/AGENTS.md](packages/HimalayaUI/frontend/src/lib/queue/AGENTS.md) |
| Vitest / JSDOM / RTL patterns | [packages/HimalayaUI/frontend/test/AGENTS.md](packages/HimalayaUI/frontend/test/AGENTS.md) |
| Playwright selectors, port binding, live-mode timing | [packages/HimalayaUI/frontend/e2e/AGENTS.md](packages/HimalayaUI/frontend/e2e/AGENTS.md) |
| Running a dev session against a live backend (port cleanup, side-by-side prod-DB, `VITE_API_PORT`) | [docs/frontend-dev-loop.md](docs/frontend-dev-loop.md) |

## Current state

- **Core Himalaya** `v0.6`: v2 peak-finding (persistence + sharpness + kneedle). See `docs/peak-finding.md` and `CHANGELOG.md`.
- **Backend:** transactional SQLite pipeline, the HTTP ingest funnel (`/api/fs/*` → `POST /api/experiments` → `scan_and_group!`), REST API (Oxygen.jl), CLI (`analyze`, `show`, `serve`, `merge`), TIFF→PNG image route with Q0f31-aware lognormalize. One experiment per beamtime folder, one `setups` row per setup subfolder; rescans are insert-only and add new setups. There is no `experiment.toml`. See `docs/experiment-config.md`.
- **Events and multiplayer:** every edit goes through `apply_event!` into the `user_actions` log, with content-hash memoization and SSE at `GET /api/events`. Multiplayer is last-write-wins with no conflict UI; undo reads the server before it restores anything. See `docs/event-log.md`.
- **Mutation queue:** a per-mutation `client_op_id` keys both the backend `with_idempotency` cache and the frontend `pendingDeferreds` registry. See `docs/mutation-queue.md`.
- **Frontend:** `index.html → src/main.tsx → App`. The experiments home and the one-screen New-experiment funnel live in `src/pages/`, the shell and onboarding in `src/shell/`. Triage happens in the three-pane workspace at `/experiments/:id/workspace/:sampleId` (`src/workspace/`): sample list, trace plot with peak editing and custom index, detector overlay, phase calls, comparisons and subtraction, bridges, labels, undo, and the ⌘K palette.
- **Design system:** closed-look components in `src/ui/`. Consumers pass placement-only `className`, and `scripts/check-design.mjs` (`lint:design`) fails the build on inline appearance utilities elsewhere. Storybook (`npm run storybook`) is the visual catalog.
- **Removed features get a fresh design if they return.** Focus, Loupe, the Index chat, the old Series/Compare page, sample chat, tags, the corpus page, grouping review and the configuration page are gone. `migrate_drop_retired_tables!` drops their tables on `open_db` and the log rows stay.
- **Deferred:** see [docs/future-feature-ideas.md](docs/future-feature-ideas.md).

## Further reading

**Design docs (`docs/`):**

- [peak-finding.md](docs/peak-finding.md) — narrative design notes, non-obvious defaults, out-of-scope decisions.
- [scoring.md](docs/scoring.md) — how and why of the index scoring formula (coverage × consistency).
- [experiment-config.md](docs/experiment-config.md) — New-experiment funnel, experiment configuration fields, file patterns, geometry, scan and grouping, rescans, read-only contract, CLI reference.
- [frontend-dev-loop.md](docs/frontend-dev-loop.md) — running a dev session against a live backend.
- [demo.md](docs/demo.md) — synthetic demo experiment, `make demo`, and running the app in a cloud session.
- [event-log.md](docs/event-log.md) — `apply_event!` dispatcher contract, hash memoization invariants, SSE multiplayer semantics.
- [mutation-queue.md](docs/mutation-queue.md) — server-reconciliation queue architecture, `client_op_id` lifecycle, replay-as-rerun, `with_idempotency`.
- [contract-testing.md](docs/contract-testing.md) — six-layer rule and canonical paired test files.
- [future-feature-ideas.md](docs/future-feature-ideas.md) — intentionally-deferred features.

**Package-specific:**

- [packages/HimalayaUI/.env.example](packages/HimalayaUI/.env.example) — deployment env vars.

