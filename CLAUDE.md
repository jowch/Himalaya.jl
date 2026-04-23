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
      routes_*.jl            # one per REST resource (users, experiments,
                             #   samples, exposures, peaks, analysis, trace, export)
      server.jl              # Oxygen.jl app + serve(db) + test harness
    test/
    frontend/                # React 18 + Vite + TS strict
      src/
        main.tsx             # entry: StrictMode > ErrorBoundary > QueryClientProvider > App
        App.tsx              # composition root; Zustand selectors + TanStack Query
        api.ts               # typed fetchers (AuthOpts per-call for mutations)
        state.ts             # Zustand client state (activeSampleId, hoveredIndexId, …)
        queries.ts           # TanStack Query hooks + queryKeys
        phases.ts            # phase → color palette
        styles.css           # Tailwind v4 + @theme tokens
        components/          # Navbar, Layout, SampleList, UserModal,
                             #   ExposureList, TraceViewer, MillerPlot,
                             #   PhasePanel, PropertiesPanel (+ tabs), …
      test/                  # Vitest + React Testing Library
      e2e/                   # Playwright (mocks /api via page.route)
      dist/                  # Vite build output; served by Oxygen.jl in prod
docs/
  peak-finding.md            # narrative design notes for findpeaks
  superpowers/               # specs and plans
test/                        # core Himalaya tests
examples/                    # scripts using Himalaya (not part of the package)
scratch/                     # gitignored — exploratory scripts and trace data
```

## Running tests

```bash
# First-time setup (also run after `git worktree add`):
(cd packages/HimalayaUI/frontend && npm install)
# Worktrees only: copy Manifest.toml from main so Himalaya core resolves to
# the local v0.5.0 — see "Himalaya core resolution in worktrees" gotcha below.

# Core Himalaya
julia --project=. -e 'using Pkg; Pkg.test()'

# HimalayaUI backend (Julia)
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'

# HimalayaUI frontend — from packages/HimalayaUI/frontend/
npm test              # Vitest unit tests (one-shot)
npm run test:watch    # Vitest watch mode
npm run e2e           # Playwright E2E (auto-starts Vite via playwright.config.ts)
npm run build         # tsc --noEmit + vite build (must pass before PR)

# One test file in isolation (core)
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

**Frontend dev loop:** run `himalaya serve` (backend on :8080) in one terminal and `npm run dev` (Vite on :5173) in another — Vite proxies `/api/*` to :8080 (see `vite.config.ts`).

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

**Phase-type serialization:** `string(Himalaya.Pn3m)` returns the fully-qualified `"Himalaya.Pn3m"`. When storing phase names in SQLite, use `string(nameof(P))` → `"Pn3m"`. The inverse is `getfield(Himalaya, Symbol(name))` (always validate with `P isa Type && P <: Himalaya.Phase` before calling `phaseratios`).

**`Himalaya` core resolution in worktrees:** `packages/HimalayaUI/Manifest.toml` pins `Himalaya = "c5c84187..." path = "../.."` to the local v0.5.0. `Manifest.toml` is gitignored, so fresh worktrees re-resolve against the registry and pull the older published v0.4.5 (which has a different `findpeaks` signature). After `git worktree add`, copy `Manifest.toml` from main before running `Pkg.test`.

## HimalayaUI frontend gotchas

**TypeScript strict + `exactOptionalPropertyTypes: true`.** `set({ username: undefined })` fails — optional fields declared as `string | undefined` rather than `username?: string` keep this ergonomic. For passing optional values through (e.g., `AuthOpts`), use the `authOpts(username)` helper in `queries.ts` which returns `{}` or `{ username }` — never `{ username: undefined }`.

**State split (load-bearing):** Zustand owns *client* state (active sample/exposure, hoveredIndexId, username). TanStack Query owns *server* state (experiments, samples, exposures, peaks, indices, groups). Mutations invalidate scoped query keys (`queryKeys.peaks(id)`, `queryKeys.groups(id)`) — don't mix the two concerns in the same hook.

**Observable Plot inside React:** the plot element has a runtime `.scale(name).invert(px)` method that isn't in DOM types; cast with `(el as unknown as { scale: ... })`. Used by TraceViewer to translate click pixel coords to q values.

**E2E selectors:** Playwright tests use `data-testid`, `role`, or stable `data-*` attributes (`data-sample-id`, `data-exposure-id`, `data-alternative-id`, `data-active`). Never assert on Tailwind class strings — they change when styling evolves.

**Tailwind v4 theming:** the dark palette is defined once in `styles.css` via `@theme { --color-* ... }`. Component files use utility classes (`bg-bg`, `text-fg-muted`, `border-accent`). If you need a new color, add it to `@theme` first.

## Current state

- Core Himalaya: `v0.5.0` on `main` — v2 peak-finding (persistence + sharpness + kneedle).
- HimalayaUI: **Plans 1–6 complete.** Backend: SQLite pipeline, REST API (Oxygen.jl), CLI. Frontend: full single-page analysis UI (sample list, trace viewer with peak editing, Miller plot, phase panel, tabbed properties panel). The web-app design spec (§1–10) is substantively implemented; see [docs/superpowers/plans/](docs/superpowers/plans/) for per-plan narrative.
- Deferred for later: Phase panel Recent section, export UI, per-user audit view, beamline-config editor, derived-exposure construction. See [docs/future-feature-ideas.md](docs/future-feature-ideas.md).

## Further reading

- [docs/peak-finding.md](docs/peak-finding.md) — narrative design notes, non-obvious defaults, out-of-scope decisions.
- [docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md](docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md) — web app design spec (schema, API, UI layout). Load-bearing for all HimalayaUI work.
- [docs/superpowers/plans/](docs/superpowers/plans/) — implementation plans (one per sub-project).
- [docs/future-feature-ideas.md](docs/future-feature-ideas.md) — intentionally-deferred features.
