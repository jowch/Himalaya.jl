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
2. [docs/experiment-config.md](docs/experiment-config.md) — required if touching `config.jl`, `manifest.jl`, or cli init/reingest.
3. [docs/scoring.md](docs/scoring.md) — required if touching `score`, `auto_group`, or `remove_subsets`.
4. [docs/event-log.md](docs/event-log.md) — required if touching `events.jl`, `hash.jl`, the `apply_event!` call sites, or the SSE handler.
5. [docs/mutation-queue.md](docs/mutation-queue.md) — required if touching `lib/queue/`, `idempotency.jl`, `with_idempotency`, or `applyRemoteToCache.ts`.
6. [docs/contract-testing.md](docs/contract-testing.md) — six-layer testing rule. Required before fixing a queue/SSE/cache reconciliation bug.

## Code layout

```
src/                                # core Himalaya package — see src/AGENTS.md
test/                               # core tests
packages/HimalayaUI/
  src/                              # backend — see packages/HimalayaUI/src/AGENTS.md
  test/                             # backend tests — see packages/HimalayaUI/test/AGENTS.md
  configs/                          # built-in experiment.toml templates (simple.toml)
  frontend/
    src/                            # frontend — see frontend/src/AGENTS.md
      print/                        # all presentation: shell/, ui/, components/, pages/, render layers
        shell/                      # app shell, routing, modals — see print/shell/AGENTS.md
      hooks/                        # custom React hooks
      lib/                          # url, plot, comparison, figure-export helpers
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
# First-time setup (also run after `git worktree add`):
(cd packages/HimalayaUI/frontend && npm install)
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.instantiate()'

# Core Himalaya
julia --project=. -e 'using Pkg; Pkg.test()'

# HimalayaUI backend (Julia) — slow, capture once. See packages/HimalayaUI/test/AGENTS.md
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1

# HimalayaUI frontend — from packages/HimalayaUI/frontend/
npm test              # Vitest unit tests (one-shot)
npm run test:watch    # Vitest watch mode
npm run e2e           # Playwright E2E (mocked; auto-starts Vite)
npm run e2e:live      # Live integration tests (requires backend + Vite running)
npm run build         # lint:design + tsc --noEmit + vite build (must pass before PR)
```

Tests use stdlib `Test` (`@testset`, `@test`, `@test_throws`). Internal (non-exported) helpers are accessed via `Himalaya.<name>` / `HimalayaUI.<name>` in tests.

**The Julia backend suite is slow** (5–10 min). Capture output once, grep the file. Same for `npm test`. Detailed slow-suite guidance in `packages/HimalayaUI/test/AGENTS.md`.

## Running the app

```bash
# Fast path via compiled sysimage (build once, ~5 min):
make sysimage          # creates build/himalaya.so
bin/himalaya config new --type simple --dir /path/to/experiment
# edit /path/to/experiment/experiment.toml to set name + column mappings
bin/himalaya init /path/to/experiment
bin/himalaya analyze /path/to/experiment
bin/himalaya serve /path/to/experiment --port 8080
# After editing manifest.csv or experiment.toml:
bin/himalaya reingest /path/to/experiment

# Without sysimage (slower cold start, no build required):
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  serve /path/to/experiment --port 8080
```

`serve` blocks. Frontend is served from `packages/HimalayaUI/frontend/dist/` if present.

**Env vars** (see `packages/HimalayaUI/.env.example`): `HIMALAYA_DB_PATH`, `HIMALAYA_CONFIGS_DIR`, `HIMALAYA_HOST`, `HIMALAYA_PORT`, `HIMALAYA_FRONTEND_DIST`.

Frontend dev loop (live backend, port cleanup, side-by-side prod-DB sessions): [docs/frontend-dev-loop.md](docs/frontend-dev-loop.md).

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

- **Greenfield cutover MERGED to `main` (PR #281, merge-commit `dcac451`, 2026-06-15):** "The Print" is the sole app. The real app lives in `packages/HimalayaUI/frontend/src/print/App.tsx` (`PrintApp`); `index.html → src/print/main.tsx` is the single entry (mounts `#app`); the legacy `src/App.tsx`/`src/main.tsx`/`print.html` were deleted. `npm run build` emits `dist/index.html` containing the real app. The `worktree-greenfield-ui-rebuild` branch has been deleted. Follow-ups on `main`: #282 (dropped the dead `/api/comparisons` frontend client) and #283 (beam-center calibration). Plan: `docs/superpowers/plans/2026-06-14-greenfield-main-cutover.md`.
- Core Himalaya: `v0.5.3` on `main` — v2 peak-finding (persistence + sharpness + kneedle).
- HimalayaUI — Plans 1–8 + Focus/Index workspace + Inspect page + Series (folio/scoping/builder) + experiment-config system + skeleton loading + multiplayer + instrumentation foundation + mutation queue + slug permalinks + figure export + functional-redesign sweep (M1–M3) + component-library extraction (enforced design system) complete:
  - **Backend:** transactional SQLite pipeline (incl. `_reingest_inner!`), FK enforcement, REST API (Oxygen.jl), CLI (`config new/list`, `init`, `analyze`, `reingest`, `show`, `serve`), TIFF→PNG image route with Q0f31-aware lognormalize, env-driven deployment.
  - **Adapter-driven I/O:** `experiment.toml` per experiment, positional or named columns, configurable file patterns, prefix-based filesystem discovery.
  - **Frontend:** corpus contact sheet → loupe → Focus workspace (`/sample/:id`: trace plate on the d3 `TracePlot` engine, detector + combs panels, phase-call + assignment rail), with clickable corpus→sample doors, peak editing + auto-fit + log/linear toggle, custom-index modal, Loupe/Inspect page (detector image + thumbnail gallery + reject-reason chips + sample metadata), OnboardingFlow + NavModal with focus trapping, skeleton loading on all data-driven surfaces. The earlier three-card *chat* Index and the @-mention subsystem were **retired 2026-05-29** (presentation deleted, message data plane parked).
  - **Plan 7 — Multiplayer + Instrumentation:** Auto/curation peak split, diff-update preserves auto peak IDs, content-hash memoization, structured `user_actions` log via `apply_event!`, SSE multiplayer at `GET /api/events`. R5b (If-Match conflict resolution) **cancelled 2026-06-03** — no conflict UI; multiplayer stays last-write-wins, replaced by edit-tracking → undo/redo → versioning (designed in Layer 4). See `docs/event-log.md` §"Conflict resolution".
  - **Plan 8 — Mutation queue + idempotency:** Per-mutation `client_op_id` keys both the backend `with_idempotency` cache and the frontend `pendingDeferreds` registry. Frontend `useQueueMutation` + `handleRemoteEvent` implement own-op confirmation and foreign-event replay-as-rerun. `analyze_run` no-op fast path suppresses both the SSE frame and the durable `user_actions` row.
  - **Series + picker + figure export + permalinks:** the standalone Compare page was **folded into Series 2026-05-29** (`/compare/*` redirects to `/series`); the Series folio/scoping/builder renders multi-trace overlays on the bespoke d3 plot engine (`src/print/plot/`) with sample-first picker and PNG/SVG copy/save. Multiplayer is last-write-wins (no conflict UI — see the Plan 7 note above). Slug-based permalink URLs round-trip through `useStateFromUrl`.
  - **Component library + enforced design system (2026-05-29):** The Print's recurring patterns extracted into a closed-look primitive library under `src/print/ui/` (~50 primitives — Button, Card, SegmentedControl, PhaseChip, Input, Field, Menu, Tooltip, … — list illustrative) — consumers pass **placement-only** `className`; appearance lives in the primitives (the closed-look/open-placement contract). Enforced by `scripts/check-design.mjs`, a **pure-absolute** `lint:design` build step (+ a warn-only PostToolUse hook) that fails the build on any inline appearance utility (`text-[…]`, `rounded-[…]`, raw colour literals, side-stripes) outside `src/print/ui/**` (rules #3/#5 allowlist the colour-authoring files: `phases.ts`, `lib/comparison/coloring.ts`, `lib/figure-export/**`, the detector/heatmap layers, `print/main.tsx`). Radius collapsed to one 5px step (`rounded.sm` == `rounded.md`); `--color-print-accent` sources from `--color-accent`; static catalog at `docs/design-system.html`. Plan: [docs/superpowers/plans/2026-05-29-component-library-extraction.md](docs/superpowers/plans/2026-05-29-component-library-extraction.md).
  - **Test coverage:** ~1000 Julia (HimalayaUI) · ~100 Julia (core) · ~279 Vitest files · 11 Playwright E2E spec files (mocked) + 4 Playwright live-integration specs.
- Deferred: holistic trace-plot-card / peak-move redesign (M4 — gated on rethinking the `auto_peaks`/`peak_curations` curation model), Phase panel Recent section, export UI, per-user audit view, derived-exposure construction. See [docs/future-feature-ideas.md](docs/future-feature-ideas.md).

## Further reading

**Design docs (`docs/`):**

- [peak-finding.md](docs/peak-finding.md) — narrative design notes, non-obvious defaults, out-of-scope decisions.
- [scoring.md](docs/scoring.md) — how and why of the index scoring formula (coverage × consistency).
- [experiment-config.md](docs/experiment-config.md) — `experiment.toml` schema, read-only contract, filename association, CLI reference.
- [frontend-dev-loop.md](docs/frontend-dev-loop.md) — running a dev session against a live backend.
- [event-log.md](docs/event-log.md) — `apply_event!` dispatcher contract, hash memoization invariants, SSE multiplayer semantics.
- [mutation-queue.md](docs/mutation-queue.md) — server-reconciliation queue architecture, `client_op_id` lifecycle, replay-as-rerun, `with_idempotency`.
- [contract-testing.md](docs/contract-testing.md) — six-layer rule and canonical paired test files.
- [future-feature-ideas.md](docs/future-feature-ideas.md) — intentionally-deferred features.

**Package-specific:**

- [packages/HimalayaUI/docs/boneyard.md](packages/HimalayaUI/docs/boneyard.md) — skeleton loading reference.
- [packages/HimalayaUI/frontend/e2e/live/README.md](packages/HimalayaUI/frontend/e2e/live/README.md) — runbook for live-integration Playwright tests.
- [packages/HimalayaUI/.env.example](packages/HimalayaUI/.env.example) — deployment env vars.

**Specs and plans:**

- [docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md](docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md) — web-app design spec (schema, API, UI).
- [docs/superpowers/specs/2026-04-28-experiment-config-design.md](docs/superpowers/specs/2026-04-28-experiment-config-design.md) — config system design.
- [docs/superpowers/plans/](docs/superpowers/plans/) — implementation plans (one per sub-project).
