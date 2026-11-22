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
packages/HimalayaUI/frontend/src/print/shell/AGENTS.md  # app shell / routing / modal gotchas
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
packages/HimalayaDB/                # read-only Julia reader of the central DB for analysis scripts; contract-tested against HimalayaUI
packages/HimalayaUI/
  src/                              # backend — see packages/HimalayaUI/src/AGENTS.md
  test/                             # backend tests — see packages/HimalayaUI/test/AGENTS.md
  frontend/
    src/                            # frontend — see frontend/src/AGENTS.md
      workspace/                    # the three-pane triage workspace (fresh build); ui/ holds its primitives
      print/                        # all presentation: shell/, ui/, components/, pages/, render layers
        shell/                      # app shell, routing, modals — see print/shell/AGENTS.md
      hooks/                        # custom React hooks
      lib/                          # url, plot, peak-ref, custom-index, label and bridge helpers
        queue/                      # mutation queue — see lib/queue/AGENTS.md
      bones/                        # committed *.bones.json skeleton captures
    test/                           # Vitest unit tests — see frontend/test/AGENTS.md
    e2e/                            # Playwright mocked — see e2e/AGENTS.md
    e2e/live/                       # Playwright live integration — see e2e/live/README.md
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
# the REGISTRY (different physics from this checkout). `[compat] Himalaya = "0.6"`
# makes that fail loudly instead of silently; `develop` is the fix, not a workaround.
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
npm run e2e:live      # Live integration tests (requires backend + Vite running)
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
| Zustand, TanStack Query, SSE wiring, Tailwind, boneyard | [packages/HimalayaUI/frontend/src/AGENTS.md](packages/HimalayaUI/frontend/src/AGENTS.md) |
| App shell / routing / modal gotchas (CorpusShell, CorpusTopbar, AppRoutes, NavModal, OnboardingFlow) | [packages/HimalayaUI/frontend/src/print/shell/AGENTS.md](packages/HimalayaUI/frontend/src/print/shell/AGENTS.md) |
| Mutation queue internals | [packages/HimalayaUI/frontend/src/lib/queue/AGENTS.md](packages/HimalayaUI/frontend/src/lib/queue/AGENTS.md) |
| Vitest / JSDOM / RTL patterns | [packages/HimalayaUI/frontend/test/AGENTS.md](packages/HimalayaUI/frontend/test/AGENTS.md) |
| Playwright selectors, port binding, live-mode timing | [packages/HimalayaUI/frontend/e2e/AGENTS.md](packages/HimalayaUI/frontend/e2e/AGENTS.md) |
| Running a dev session against a live backend (port cleanup, side-by-side prod-DB, `VITE_API_PORT`) | [docs/frontend-dev-loop.md](docs/frontend-dev-loop.md) |

## Current state

- **Greenfield cutover MERGED to `main` (PR #281, merge-commit `dcac451`, 2026-06-15):** "The Print" is the sole app. The real app lives in `packages/HimalayaUI/frontend/src/print/App.tsx` (`PrintApp`); `index.html → src/print/main.tsx` is the single entry (mounts `#app`); the legacy `src/App.tsx`/`src/main.tsx`/`print.html` were deleted. `npm run build` emits `dist/index.html` containing the real app. The `worktree-greenfield-ui-rebuild` branch has been deleted. Follow-ups on `main`: #282 (dropped the dead `/api/comparisons` frontend client) and #283 (beam-center calibration).
- **Three-pane workspace replaces Focus, Loupe and Series (2026-11-20):** triage happens at `/experiments/:id/workspace/:sampleId` (`src/workspace/`). The Focus, Loupe and Series pages, their components and their tests were deleted; `/sample/:id`, `/sample/:id/loupe` and `/index/:exp/:sample` redirect to the workspace. The experiments home and ingest funnel remain. **The corpus page was deleted 2026-11-21:** opening an experiment (`/experiments/:id`, `ExperimentEntry`) goes straight to the workspace on the first unreviewed sample; Rescan and "Show ingest problems" are ⌘K palette commands, and a notice above the sample list shows ingest progress, a failed ingest and exposures without a trace or with a failed analysis. **The grouping review page, its merge/split/move/rename routes, the configuration page and the scan-failed page were deleted 2026-11-21** (ingest redesign PR 3); a mistaken split is fixed with a bridge. Some of the history below describes the deleted pages.
- Core Himalaya: `v0.5.3` on `main` — v2 peak-finding (persistence + sharpness + kneedle).
- HimalayaUI — Plans 1–8 + Focus/Index workspace + Inspect page + Series (folio/scoping/builder) + experiment-config system + skeleton loading + multiplayer + instrumentation foundation + mutation queue + slug permalinks + figure export + functional-redesign sweep (M1–M3) + component-library extraction (enforced design system) complete:
  - **Backend:** transactional SQLite pipeline, HTTP ingest funnel (`/api/fs/*` → `POST /api/experiments` → `scan_and_group!`), FK enforcement, REST API (Oxygen.jl), CLI (`analyze`, `show`, `serve`), TIFF→PNG image route with Q0f31-aware lognormalize, env-driven deployment.
  - **Ingest:** the New-experiment funnel infers `data_dir`/`analysis_dir`/file patterns/geometry from the directory tree; `scan_and_group!` groups exposures into loads and samples by time and position, insert-only on rescan. There is no `experiment.toml`. See `docs/experiment-config.md`.
  - **Frontend:** experiments home, the one-screen New-experiment funnel (every setup under a picked folder), and the three-pane triage workspace (`src/workspace/`: sample list, trace plot with peak editing and custom index, detector overlay, phase calls, compare and bridges, labels). OnboardingFlow + NavModal with focus trapping; skeleton loading on data-driven surfaces. The chat Index, @-mentions, Focus, Loupe, the corpus page and Series are gone.
  - **Plan 7 — Multiplayer + Instrumentation:** Auto/curation peak split, diff-update preserves auto peak IDs, content-hash memoization, structured `user_actions` log via `apply_event!`, SSE multiplayer at `GET /api/events`. R5b (If-Match conflict resolution) **cancelled 2026-06-03** — no conflict UI; multiplayer stays last-write-wins, replaced by edit-tracking → undo/redo → versioning (designed in Layer 4). See `docs/event-log.md` §"Conflict resolution".
  - **Plan 8 — Mutation queue + idempotency:** Per-mutation `client_op_id` keys both the backend `with_idempotency` cache and the frontend `pendingDeferreds` registry. Frontend `useQueueMutation` + `handleRemoteEvent` implement own-op confirmation and foreign-event replay-as-rerun. `analyze_run` no-op fast path suppresses both the SSE frame and the durable `user_actions` row.
  - **Series and comparisons deleted (2026-11-21):** the Compare page was folded into Series 2026-05-29, Series was hidden 2026-11-20 and then deleted on both sides 2026-11-21: frontend, `/api/series*` and picker routes, the series/comparison event kinds, and their nine tables (`migrate_drop_retired_tables!` drops them on `open_db`; the log rows stay). Sample chat (`/api/samples/{id}/messages`, `sample_messages`, the `post_message` kind) was deleted the same way. `/series*` and `/compare*` are now stale URLs. Anything similar gets a fresh design. Multiplayer is last-write-wins (no conflict UI — see the Plan 7 note above). Slug-based permalink URLs round-trip through `useStateFromUrl`.
  - **Component library + enforced design system (2026-05-29):** The Print's recurring patterns extracted into a closed-look primitive library under `src/print/ui/` (Button, Card, PhaseChip, Input, Popover, Tooltip, …; unused primitives were deleted 2026-11-22) — consumers pass **placement-only** `className`; appearance lives in the primitives (the closed-look/open-placement contract). Enforced by `scripts/check-design.mjs`, a **pure-absolute** `lint:design` build step (+ a warn-only PostToolUse hook) that fails the build on any inline appearance utility (`text-[…]`, `rounded-[…]`, raw colour literals, side-stripes) outside the appearance-exempt dirs (`src/print/ui/**`, the new workspace's `src/workspace/ui/**`, plus the `print/{plot,detector}/` render layers) (rules #3/#5 also allowlist the colour-authoring files `phases.ts` and `print/main.tsx`). Radius collapsed to one 5px step (`rounded.sm` == `rounded.md`); `--color-print-accent` sources from `--color-accent`; the live visual catalog is Storybook (`npm run storybook` from `packages/HimalayaUI/frontend`) — the hand-maintained `docs/design-system.html` mock was retired 2026-06-21 (it drifted from the living primitives; the `.stories.tsx` files are the source of truth).
  - **Test coverage:** ~1000 Julia (HimalayaUI) · ~100 Julia (core) · ~279 Vitest files · 11 Playwright E2E spec files (mocked) + 4 Playwright live-integration specs.
- Deferred: holistic trace-plot-card / peak-move redesign (M4 — gated on rethinking the `auto_peaks`/`peak_curations` curation model), Phase panel Recent section, export UI, per-user audit view, derived-exposure construction. See [docs/future-feature-ideas.md](docs/future-feature-ideas.md).

## Further reading

**Design docs (`docs/`):**

- [peak-finding.md](docs/peak-finding.md) — narrative design notes, non-obvious defaults, out-of-scope decisions.
- [scoring.md](docs/scoring.md) — how and why of the index scoring formula (coverage × consistency).
- [experiment-config.md](docs/experiment-config.md) — New-experiment funnel, experiment configuration fields, file patterns, geometry, scan and grouping, rescans, read-only contract, legacy TOML, CLI reference.
- [frontend-dev-loop.md](docs/frontend-dev-loop.md) — running a dev session against a live backend.
- [demo.md](docs/demo.md) — synthetic demo experiment, `make demo`, and running the app in a cloud session.
- [event-log.md](docs/event-log.md) — `apply_event!` dispatcher contract, hash memoization invariants, SSE multiplayer semantics.
- [mutation-queue.md](docs/mutation-queue.md) — server-reconciliation queue architecture, `client_op_id` lifecycle, replay-as-rerun, `with_idempotency`.
- [contract-testing.md](docs/contract-testing.md) — six-layer rule and canonical paired test files.
- [future-feature-ideas.md](docs/future-feature-ideas.md) — intentionally-deferred features.

**Package-specific:**

- [packages/HimalayaUI/docs/boneyard.md](packages/HimalayaUI/docs/boneyard.md) — skeleton loading reference.
- [packages/HimalayaUI/frontend/e2e/live/README.md](packages/HimalayaUI/frontend/e2e/live/README.md) — runbook for live-integration Playwright tests.
- [packages/HimalayaUI/.env.example](packages/HimalayaUI/.env.example) — deployment env vars.

> The per-sub-project design specs and implementation plans that used to live under `docs/superpowers/` were consolidated into the living docs above and removed (2026-06-18); their full history is preserved in git.
