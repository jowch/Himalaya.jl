# Himalaya.jl — Claude working notes

## What this is

A Julia monorepo for **indexing SAXS diffraction patterns**. The core `Himalaya` package finds Bragg peaks in a 1D integration trace and identifies the liquid-crystalline phase (Pn3m, Im3m, Ia3d, Fm3m, Fd3m, Hexagonal, Lamellar, Square) by fitting peak q-values to known phase-ratio series. `HimalayaUI` (under `packages/`) is a full-stack web app — Julia/Oxygen.jl REST backend + React/Vite frontend — for running and curating analyses on a batch of SAXS exposures.

## Read first

If this is your first session on this repo, skim these in order before touching code:

1. [docs/peak-finding.md](docs/peak-finding.md) — why findpeaks is the way it is. Load-bearing.
2. [docs/experiment-config.md](docs/experiment-config.md) — required if touching `config.jl`, `manifest.jl`, or cli init/reingest.
3. [docs/scoring.md](docs/scoring.md) — required if touching `score`, `auto_group`, or `remove_subsets`.

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
  index.jl                   # Index struct, indexpeaks, score
  util.jl
packages/
  HimalayaUI/                # web-app sub-package
    .env.example             # documented env vars (HIMALAYA_DB_PATH etc.)
    configs/                 # built-in experiment.toml templates (simple.toml)
    src/
      db.jl                  # SQLite schema + CRUD
      datfile.jl             # three-column .dat parser
      config.jl              # ExperimentConfig + load_config + resolve_files
      manifest.jl            # ManifestSample + parse_manifest (config-driven)
      pipeline.jl            # analyze_exposure!, auto_group, persist_analysis!
      cli.jl                 # himalaya config/init/analyze/reingest/show/serve
      json.jl                # row → Dict serialization
      actions.jl             # X-Username extraction + user_actions logger
      image.jl               # TIFF load + log-normalize + PNG encode for /image route
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
        components/          # AppShell, AppHeader, TabRocker, TitleButton,
                             #   OnboardingFlow, NavModal, UtilityCluster,
                             #   ChatCard, MentionChip, MentionCompose,
                             #   PlotCard, TraceViewer,
                             #   IndicesCard, PhasePanel, StaleIndicesBanner,
                             #   MillerPlot, Pn3mIcon, ui/…
                             # Inspect: DetectorImage, DetectorImageCard,
                             #   ThumbnailGallery, SampleMetadataCard
                             # Mentions: MentionChip, MentionCompose, MentionPicker
        hooks/               # useFocusTrap, useMentionResolution
        lib/                 # renderMentions.tsx (parseMentions tokenizer)
        bones/               # Committed boneyard skeleton captures (*.bones.json)
                             #   + auto-generated registry.ts
        pages/               # IndexPage (three-card workspace),
                             #   InspectPage (curate exposures), ComparePage
      test/                  # Vitest + React Testing Library
      e2e/                   # Playwright (mocks /api via page.route)
      dist/                  # Vite build output; served by Oxygen.jl in prod
docs/
  peak-finding.md            # findpeaks design (persistence + sharpness + kneedle)
  scoring.md                 # index scoring formula rationale
  experiment-config.md       # experiment.toml format + read-only contract
  future-feature-ideas.md    # intentionally-deferred features
  superpowers/               # specs and plans
test/                        # core Himalaya tests
examples/                    # scripts using Himalaya (not part of the package)
scratch/                     # gitignored — exploratory scripts and trace data
.claude/
  skills/                    # project-specific Claude Code skills:
                             #   review-pr, worktree-setup, new-route,
                             #   review-response, e2e-mock-mode, seed-test-state
  agents/                    # frontend-reviewer (project-specific review agent)
  settings.json              # hooks: Vitest --run flag, pre-tool-use guards
```

## Running tests

```bash
# First-time setup (also run after `git worktree add`):
(cd packages/HimalayaUI/frontend && npm install)
# Worktrees only: copy Manifest.toml files from main so Himalaya core resolves to
# the local v0.5.0 — see "Himalaya core resolution in worktrees" gotcha below.
# Also copy (or instantiate) scripts/Manifest.toml before running `make sysimage`.

# Core Himalaya
julia --project=. -e 'using Pkg; Pkg.test()'

# HimalayaUI backend (Julia)
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'

# HimalayaUI frontend — from packages/HimalayaUI/frontend/
npm test              # Vitest unit tests (one-shot)
npm run test:watch    # Vitest watch mode
npm run e2e           # Playwright E2E (auto-starts Vite via playwright.config.ts)
npm run build         # tsc --noEmit + vite build (must pass before PR)

# Single frontend test file / single E2E test by name
node_modules/.bin/vitest run test/DetectorImage.test.tsx
node_modules/.bin/playwright test --grep "Reject → Other"

# One test file in isolation (core)
julia --project=. -e 'using Himalaya, Test; include("test/foo.jl")'
```

Tests use stdlib `Test` (`@testset`, `@test`, `@test_throws`). Internal (non-exported) helpers are accessed via `Himalaya.<name>` in tests.

## Running the app

```bash
# Fast path via compiled sysimage (build once, ~5 min):
make sysimage          # creates scratch/himalaya.so
bin/himalaya config new --type simple --dir /path/to/experiment
# edit /path/to/experiment/experiment.toml to set name + column mappings
bin/himalaya init /path/to/experiment
bin/himalaya analyze /path/to/experiment
bin/himalaya serve /path/to/experiment --port 8080
# After editing manifest.csv or experiment.toml:
bin/himalaya reingest /path/to/experiment

# Without sysimage (slower cold start, no build required):
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  config new --type simple --dir /path/to/experiment
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  init /path/to/experiment
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  serve /path/to/experiment --port 8080
```

`serve` blocks. Frontend is served from `packages/HimalayaUI/frontend/dist/` if present.

**Env vars** (see `packages/HimalayaUI/.env.example`): `HIMALAYA_DB_PATH` overrides the per-experiment DB path for `/opt`-style centralised deployment; `HIMALAYA_CONFIGS_DIR`, `HIMALAYA_HOST`, `HIMALAYA_PORT`, `HIMALAYA_FRONTEND_DIST` override the corresponding defaults.

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
- **`Tables.rowtable` returns `missing` for SQL NULL, not `nothing`.** When comparing nullable columns in Julia code (e.g. "is this name field unset?"), use `ismissing(row.field)` and normalize to `nothing` if you need to mix with literals: `existing = ismissing(row.field) ? nothing : row.field`. The `routes_users.jl` NULL-fill enrichment path is the canonical example.
- **FK enforcement is on.** `open_db` runs `PRAGMA foreign_keys = ON` on every connection. Any FK column that references `users(id)` and must survive user deletion needs `ON DELETE SET NULL` in the schema DDL — add it there, not at call sites. `index_groups.created_by` and `user_actions.user_id` already have this.
- **PKs use `AUTOINCREMENT` on mention-targets** (`experiments`, `samples`, `exposures`, `peaks`, `indices`). Plain `INTEGER PRIMARY KEY` is rowid-aliased and reuses freed ids on deletion — chat `@`-mentions of a deleted entity could silently rebind to a new entity that took the same id. `migrate_pk_to_autoincrement!` in `db.jl` rebuilds these tables on existing DBs (rename → `create_schema!` → copy shared columns → drop) so legacy DBs heal on next `open_db`.
- **`analyze_exposure!` curation contract.** Before calling `Himalaya.indexpeaks`, the function (1) drops auto peaks the user previously excluded and (2) unions in current manual peaks with sharpness sampled from `Himalaya.sharpness(I)`. Without these, a user's "this is noise" exclusion has no effect on candidate scoring, and a user-marked manual peak at a phase's predicted ratio position never lands in `IndexEntry.peaks`. Touch this only with regression tests in `test_pipeline.jl` (the "incorporates manual peaks" / "ignores excluded auto peaks" testsets) green.
- **`persist_analysis!` is transactional.** The delete-then-reinsert sequence in `pipeline.jl` is wrapped in `SQLite.transaction`. If you add new write steps to that function, put them inside `_persist_analysis_inner!` so they stay atomic. Same pattern applies to `reingest!` in `cli.jl` (`_reingest_inner!`) — wrap any new multi-write CLI operations the same way. `_reingest_inner!` returns a `NamedTuple{(:status, :added_samples, :added_exposures, :manifest_path)}` where `:status` is `:ok` or `:no_manifest`; the route at `POST /api/experiments/:id/reingest` echoes those fields in JSON (HTTP 200 in both cases — `status` is the discriminator).

**Oxygen.jl 1.10.x:**
- Use the singleton API: `@get "/path/{id}" function(req::HTTP.Request, id::Int) ... end`. Typed function args extract path params.
- Parse JSON body with `json(req)` (unqualified, imported via `using Oxygen`), **not** `Oxygen.json(req)`.
- Test harness pattern: `Oxygen.resetstate()` before `Oxygen.serve(; async=true)`, `Oxygen.terminate()` after. A module-level `Ref{Union{SQLite.DB, Nothing}}` holds the live DB (matches the one-experiment-per-process deployment model).
- Mount static files with `Oxygen.dynamicfiles(dir, "/")` — only if `isdir(dir)`, so empty frontends don't break tests.
- Oxygen emits a harmless warning about OpenAPI schema generation for some routes; ignore it.

**Stdlib deps must be explicit.** Stdlibs used directly in a package (`Sockets`, `Printf`, `SparseArrays`, `DelimitedFiles`, `TOML`, etc.) must be listed in `Project.toml`'s `[deps]` — `Pkg.add` them like regular packages.

**Experiment config (`experiment.toml`) is the source of truth.** Every experiment directory has an `experiment.toml` describing manifest column layout, file patterns, and beamline params. Generate one with `himalaya config new --dir <path>` (the only command that writes inside an experiment directory; refuses to overwrite). `himalaya init` reads it and stores the full TOML blob in `experiments.config` so the DB is self-contained. `analyze_exposure!` reads the integration pattern via `config_from_db`, falling back to `simple.toml` defaults when `config IS NULL` (legacy experiments). `config_from_db` and `load_config` share a `_build_config(::AbstractDict)` helper — the DB blob is parsed in-memory (no tempfile) so `analyze_exposure!` doesn't pay disk I/O per call. `layout.exposure_type` is validated at parse time against `VALID_EXPOSURE_TYPES` (currently `("simple",)`); extend that tuple before introducing a new exposure type. Malformed `experiment.toml` produces a wrapped `Invalid TOML in <path>: …` error. To change the layout or column mapping, edit `experiment.toml` then run `himalaya reingest <path>` — preserves curation (status, manual peaks). Read [docs/experiment-config.md](docs/experiment-config.md) before touching `config.jl`, `manifest.jl`, or the cli init/reingest paths.

**Read-only experiment directories at runtime.** Himalaya never creates, modifies, or deletes any file inside an experiment directory during `init`, `analyze`, `reingest`, or `serve`. The sole exception is `himalaya config new --dir`, which writes `experiment.toml` once during setup. A regression test in `test_pipeline.jl` snapshots the directory contents before/after `cli_init_with_db!` — keep it green.

**Central DB.** All CLI commands open the same DB resolved by `default_db_path()` in `db.jl`: `HIMALAYA_DB_PATH` if set, else `~/.himalaya/himalaya.db` (parent dir auto-created). One DB stores every experiment ever registered; experiment dirs are pure read-only data sources. Tests pass an explicit file path (`open_db(joinpath(tmp, "himalaya.db"))`) to keep each testset isolated.

**Filename ↔ exposure association via filesystem prefix scan.** Manifest filename entries are always treated as prefixes: `JC001-004` expands to four prefixes, each scanned via `resolve_files(cfg, dir, prefix, cfg.integration_pattern)` against the filesystem. The manifest declares intent; disk decides what exists. Missing files produce a warning, not an error. When debugging "exposures missing after init/reingest," check the actual files in `analysis_dir` first — the manifest is a hint.

**`parse_manifest` has two methods.** `parse_manifest(source)` is a backward-compat wrapper using `simple.toml` defaults. `parse_manifest(cfg::ExperimentConfig, source)` is the config-driven version — use this in new code. Both accept IO and paths via `readlines(source)`.

**Index scoring:** `score(index)` returns a value in `[0, 1]` — product of `coverage` (harmonic-weighted fraction of expected peaks found, `1/rank` weight per position) and `consistency` (`1/(1+CV)` of peak sharpnesses). `totalprom` and the `prom` field on `Index` no longer exist — the struct now has `sharpness::SparseVector`. Guard `cv` against zero mean before dividing (all-zero sharpness is valid and should score as consistent). `auto_group` and `remove_subsets` in `pipeline.jl` both depend on `score` ordering — correctness of auto-analysis flows from score quality. R² is stored per index but is NOT part of the score; it is a UI hard gate (threshold 0.98 in `PhasePanel`).

**Phase-type serialization:** `string(Himalaya.Pn3m)` returns the fully-qualified `"Himalaya.Pn3m"`. When storing phase names in SQLite, use `string(nameof(P))` → `"Pn3m"`. The inverse is `getfield(Himalaya, Symbol(name))` (always validate with `P isa Type && P <: Himalaya.Phase` before calling `phaseratios`).

**Detector TIFFs are Q0f31 fixed-point.** TiffImages loads as `Gray{Q0f31}` (= `Fixed{Int32, 31}`). `Float32.(channelview(raw))` divides raw photon counts by 2³¹ (~2.1e9), making `log1p` numerically a no-op (~4.7e-10 per count). To recover photon counts, use `reinterpret.(Int32, channelview(raw))`. Then `max(., 0)` clips beamstop/dead-pixel negatives, `log1p` compresses, and a p99-of-positives clip prevents the direct beam from crushing diffraction-ring contrast. See `image.jl::load_and_lognormalize`.

**Image route uses `Cache-Control: no-store`** (`routes_exposures.jl`). Combined with `cache: "no-store"` on the frontend `fetch()`, this stops the browser from serving stale PNGs across analysis re-runs. Don't change to a longer max-age without invalidation tied to exposure id + analysis version.

**`Himalaya` core resolution in worktrees:** `packages/HimalayaUI/Manifest.toml` pins `Himalaya = "c5c84187..." path = "../.."` to the local v0.5.0. `Manifest.toml` is gitignored, so fresh worktrees re-resolve against the registry and pull the older published v0.4.5 (which has a different `findpeaks` signature). After `git worktree add`, copy `Manifest.toml` from main before running `Pkg.test`.

## HimalayaUI frontend gotchas

**TypeScript strict + `exactOptionalPropertyTypes: true`.** `set({ username: undefined })` fails — optional fields declared as `string | undefined` rather than `username?: string` keep this ergonomic. For passing optional values through (e.g., `AuthOpts`), use the `authOpts(username)` helper in `queries.ts` which returns `{}` or `{ username }` — never `{ username: undefined }`.

**Zustand: use named actions, not `setState`.** The store exposes specific actions (`clearUsername`, `setTheme`, `openNavModal`, etc.). Prefer these over `useAppState.setState({ ... })` — direct setState bypasses encapsulation and triggers lint warnings. If you need a new state transition, add a named action to `state.ts`.

**State split (load-bearing):** Zustand owns *client* state (active sample/exposure, hoveredIndexId, username). TanStack Query owns *server* state (experiments, samples, exposures, peaks, indices, groups). Mutations invalidate scoped query keys (`queryKeys.peaks(id)`, `queryKeys.groups(id)`) — don't mix the two concerns in the same hook.

**Observable Plot inside React:** the plot element has a runtime `.scale(name).invert(px)` method that isn't in DOM types; cast with `(el as unknown as { scale: ... })`. Used by TraceViewer to translate click pixel coords to q values.

**E2E selectors:** Playwright tests use `data-testid`, `role`, or stable `data-*` attributes (`data-sample-id`, `data-exposure-id`, `data-alternative-id`, `data-active`, `data-low-r2`). Never assert on Tailwind class strings — they change when styling evolves. For Vitest/RTL tests, use `screen.getByText("X").closest("li")` + `toHaveAttribute` rather than `document.querySelector` — the latter bypasses RTL's async-aware retry logic.

**Playwright port binding:** `playwright.config.ts` expects the dev server on `http://127.0.0.1:5173`, not `localhost`. If another process has that port, tests hang for 60 s then fail. Kill with `lsof -ti:5173 | xargs kill -9`, then re-run. If starting Vite separately before `npm run e2e`, bind it explicitly: `npm run dev -- --host 127.0.0.1`.

**Focus trapping in modals.** `src/hooks/useFocusTrap.ts` exports `useFocusTrap(containerRef, active)`. Call it inside any modal or overlay that should keep Tab focus within its bounds. It intercepts Tab/Shift+Tab to cycle among focusable children and restores the previously-focused element on cleanup. NavModal and OnboardingFlow already use it.

**`QNumInput` is exported from `PlotCard.tsx`** for unit testing. It implements a focus-gated controlled input: external `value` prop changes are synced to draft state only when the input is not focused, preventing wheel-zoom events from interrupting mid-edit. Follow this pattern for any numeric input that can be updated by external events.

**`StaleIndicesBanner` is mounted in `PhasePanel`.** When any index has `status: "stale"`, the banner renders above the index list with a Re-analyze button that posts to `/api/exposures/:id/analyze`. If you add new routes that mark indices stale, no further UI wiring is needed — the banner appears automatically.

**Imperative render functions in effects: use `useCallback`.** Wrap any function that is both defined inside a component and used as a `useEffect` dependency in `useCallback` with its true deps. The effect then depends on `[theCallback]` alone — no redundant dep list, no eslint-disable. `TraceViewer`'s overlay renderer follows this pattern.

**`Fm3m` missing from `indexpeaks` dispatch** — the all-phases loop in `src/index.jl` omits `Fm3m`. The phase is defined and `minpeaks`/`phaseratios` exist, but `indexpeaks` can never return an `Fm3m` index. Known pre-existing gap, not something to fix opportunistically.

**Tailwind v4 theming:** the dark palette is defined once in `styles.css` via `@theme { --color-* ... }`. Component files use utility classes (`bg-bg`, `text-fg-muted`, `border-accent`). If you need a new color, add it to `@theme` first.

**`ImageBitmap.close()` neuters width/height** to 0 per the Web spec. Capture dims **before** closing: `const { width, height } = bitmap; bitmap.close();`. There's a regression test in `test/DetectorImage.test.tsx` using getter-based mocks that simulates the neutering — keep it green if you touch `DetectorImage.tsx`.

**`DetectorImage` auto-rotates to landscape** via a ResizeObserver on the wrapper div: when `containerAspect > imageAspect * 1.25`, rotate the canvas 90° and JS-set `maxWidth`/`maxHeight` to swap the layout box (CSS-only doesn't cut it because `transform: rotate` doesn't change a canvas's bounding box). Re-evaluates inside `renderImage` after each new image so swapping exposures with different aspects re-checks. JSDOM lacks `ResizeObserver` — the stub in `test/setup.ts` keeps unit tests honest.

**TraceViewer auto-fit is floor-only.** `PlotCard::computeFit` sets `yDomain = [max(p01·0.5, fullMax/1e5), fullMax·1.2]` — bottom is the 1st percentile of *positive* in-window intensities scaled down (suppresses dead-pixel zeros while keeping the low-signal tail visible), clamped at `fullMax/1e5` so a single near-zero pixel can't blow the y-range past five log decades. Top is the *full* trace max (so peaks-vs-beam relative scale stays visible without resetting). When peaks exist, x is also tightened to `[firstPeak·0.7, lastPeak·1.3]`. Auto-fires on `activeExposureId` change. Double-click → `onReset` clears both axes.

**Skeleton loading via boneyard-js.** Each load-gated card/list (PlotCard, PhasePanel, ChatCard, DetectorImageCard, SampleMetadataCard, NavModal) wraps its content in `<Skeleton>` from `boneyard-js/react`. Several non-obvious rules:
- **Use `loading={query.isLoading}`, not `isPending`.** `isLoading = isPending && isFetching` — disabled queries (`undefined` selection) and background refetches over cached data both stay skeleton-free; only true cold fetches trigger the animation. Wrong gating causes skeleton flicker on every refetch.
- **`fixture` prop is `ReactNode`** (JSX rendered during boneyard's headless capture so the CLI can measure layout), NOT raw data. Pass real components with mock props — e.g. `<TraceViewer trace={…} peaks={…} … />` with no-op handlers.
- **Always set `fallback`** to mirror the original italic HintText. When bones aren't yet captured for that skeleton name, the runtime renders `fallback` during loading; without it the area is blank.
- **`className` on `<Skeleton>` is load-bearing** — boneyard wraps children in two extra divs which would otherwise break parent flex chains (e.g. ChatCard's message list collapsed to 60px). Pass `flex-1 min-h-0 flex flex-col` (or `h-full w-full`) so the outer wrapper inherits the original child's layout role. Companion CSS in `styles.css`: `[data-boneyard-content] { display: contents }` makes the inner wrapper transparent to layout so the children's flex behaviour passes through.
- **`configureBoneyard()` lives in `main.tsx`,** NOT in `bones/registry.ts`. The Vite HMR plugin regenerates `registry.ts` on every capture and would wipe any config call there. The values must mirror `boneyard.config.json` (which the capture CLI reads but the runtime does not) — both files have to be updated together when the card background colour or animation changes.
- **Bones are committed,** not gitignored. `src/bones/*.bones.json` + the auto-generated `registry.ts` are required for prod builds to render skeletons; without them, `<Skeleton>` falls through to `fallback`. Capture organically during dev (the Vite plugin re-captures on every HMR update) and commit deliberately to widen prod coverage.
- **JSDOM lacks `window.matchMedia`** — boneyard's dark-mode detection calls it on mount. The stub in `test/setup.ts` keeps unit tests honest; same pattern as the `ResizeObserver` stub above.

## Current state

- Core Himalaya: `v0.5.0` on `main` — v2 peak-finding (persistence + sharpness + kneedle).
- HimalayaUI — Plans 1–6 + three-card Index redesign + Inspect page + experiment-config system + skeleton loading + chat @-mentions + build infrastructure complete:
  - **Backend:** transactional SQLite pipeline (incl. `_reingest_inner!`), FK enforcement, REST API (Oxygen.jl), CLI (`config new/list`, `init`, `analyze`, `reingest`, `show`, `serve`), TIFF→PNG image route with Q0f31-aware lognormalize, env-driven deployment (`HIMALAYA_DB_PATH`, `HIMALAYA_CONFIGS_DIR`).
  - **Adapter-driven I/O:** `experiment.toml` per experiment, positional or named columns, configurable file patterns, prefix-based filesystem discovery.
  - **Frontend:** three-card Index workspace (chat | trace plot | index choices), Inspect page (detector image + thumbnail filmstrip + reject-reason chips + sample metadata), trace viewer with peak editing + auto-fit y-floor + log/linear x toggle, auto-rotating detector canvas, Miller plot, PhasePanel with curate + stale-indices reanalyze, OnboardingFlow + NavModal with focus trapping. Skeleton loading screens via boneyard-js on all major data-driven cards. Chat @-mention system (`@peak`, `@index`, `@exposure`, `@sample`) via `MentionChip` / `MentionCompose` / `useMentionResolution`.
  - **Test coverage:** 452 Julia (HimalayaUI) · 90 Julia (core) · 174 Vitest · 16 Playwright E2E (7 inspect + 9 smoke).
- Deferred for later: Phase panel Recent section, export UI, per-user audit view, derived-exposure construction (raw / aggregated / background-subtracted exposure types — schema reserves `exposure_type` field), additional config templates beyond `simple.toml`. See [docs/future-feature-ideas.md](docs/future-feature-ideas.md).

## Further reading

- [docs/peak-finding.md](docs/peak-finding.md) — narrative design notes, non-obvious defaults, out-of-scope decisions.
- [docs/scoring.md](docs/scoring.md) — how and why of the index scoring formula (coverage × consistency).
- [docs/experiment-config.md](docs/experiment-config.md) — `experiment.toml` schema, read-only contract, filename association, CLI reference. Required reading before touching `config.jl`, `manifest.jl`, or the cli init/reingest paths.
- [docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md](docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md) — web app design spec (schema, API, UI layout). Load-bearing for all HimalayaUI work.
- [docs/superpowers/specs/2026-04-28-experiment-config-design.md](docs/superpowers/specs/2026-04-28-experiment-config-design.md) — config system design spec.
- [docs/superpowers/plans/](docs/superpowers/plans/) — implementation plans (one per sub-project).
- [docs/future-feature-ideas.md](docs/future-feature-ideas.md) — intentionally-deferred features.
- [packages/HimalayaUI/.env.example](packages/HimalayaUI/.env.example) — deployment env vars.
