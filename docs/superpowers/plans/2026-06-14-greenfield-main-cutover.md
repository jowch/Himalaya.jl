# Greenfield → main cutover plan

**Date:** 2026-06-14 (revised 2026-06-15 after red-team `w5brxg3pj`)
**Branch:** `worktree-greenfield-ui-rebuild` (stays UNMERGED until Jonathan's go-ahead)
**Status:** PLAN, de-risked — red-team verdict **go-with-fixes**, fixes folded in below.
Load-bearing claims independently verified against live source (see "Evidence").

---

## TL;DR

The cutover is a **promotion**, landed as a **merge-commit**. The production build
entry `PrintApp` (`src/print/App.tsx`) is a **15-line stub**; the real app
(`AppRoutes` + SSE + mutation queue) lives in `src/App.tsx`, reachable only from
the legacy `index.html → src/main.tsx`. So `npm run build` ships a blank page
today. The fix:

1. **Promote** the real app body into `PrintApp`.
2. Make **`index.html` the single greenfield entry** (repoint its script to
   `/src/print/main.tsx`, keep `#app`), and **delete `print.html`** + the legacy
   `src/main.tsx` + `src/App.tsx`.
3. Gate on build + e2e + backend suite proving the bundle is the real app, not the
   stub.

**Why `index.html`-canonical (not `print.html`):** the backend static-serve
(`server.jl:70`), the SPA-fallback test (`test_spa_fallback.jl`), the CSS height
rule (`styles.css:134 #app`), and the Vite dev-root all already assume
`index.html`. Keeping it as the entry makes all four work unchanged;
`print.html`-canonical fought all four (red-team `w5brxg3pj`). `print.html` was
rebuild-era scaffolding.

## Evidence (verified)

| Claim | Status | Evidence |
|---|---|---|
| `PrintApp` (`src/print/App.tsx`) is a 15-line stub | **VERIFIED** | renders only `<main data-testid="print-shell"><h1>The Print · greenfield shell</h1></main>` |
| `src/App.tsx` is the real app, legacy-entry-only | **VERIFIED** | imports `AppRoutes` + `EventSource('/api/events')` + queue `attachPersistence`/`rehydrate` + shell siblings |
| `index.html` mounts `#app` via `/src/main.tsx`; `print.html` mounts `#print-root` via `/src/print/main.tsx` | **VERIFIED** | both files read |
| `styles.css:134` height chain keys on `#app` | **VERIFIED** | `html, body, #app { height: 100%; }` — so the entry MUST mount `#app` or the app collapses to 0px |
| `server.jl:70` + `test_spa_fallback.jl:17` both use `index.html` | **VERIFIED** | server reads `joinpath(dist_dir,"index.html")`; test writes `dist/index.html` and asserts it is served |
| `lib/queue/**` is LIVE | **VERIFIED** | imported by `src/App.tsx` (the real app) + ~15 `test/queue/**` specs |
| Backend `index_groups`/`migrate_assignments!`/`active_group_kind` are KEEP | per survey `w70n20q1o` (high) | migration-load-bearing / unconsumed-but-public; cutover does NOT touch them |

## Decisions (confirmed by Jonathan)

- **Merge strategy: merge-commit (`--no-ff`).** The branch evolved frontend AND the
  Julia backend/data-model; the backend migration history is disaster-recovery
  load-bearing and must stay bisectable. Not squash, not fresh-from-main.
- **Entry strategy: single entry — `index.html` canonical** (flipped 2026-06-15
  after the red-team showed `print.html`-canonical fights the serve/test/CSS/dev-root
  stack). One html (`index.html`), one main (`src/print/main.tsx`), one app
  (`PrintApp`).

## Manifest

### PROMOTE (the load-bearing edits)
- **`src/print/App.tsx`** — replace the stub body with `src/App.tsx`'s body:
  `AppRoutes` + `OnboardingFlow` + `NavModal` + `ToastContainer`/`LiveRegion` +
  `InfrastructureBanner`, plus the 4 effects (DEV-gated `exposeTestHelpers`;
  `EventSource('/api/events')` `curation` → `handleRemoteEvent`;
  `attachPersistence`; `rehydrate` w/ dropped/failed toasts) + the 6-element JSX
  fragment. Keep the export name `PrintApp`. Drop the "Foundation placeholder"
  docstring. **Mirror `src/App.tsx` exactly** — do NOT add a conflict
  modal/bridge or a `useAppState` import (that red-team claim broke under scrutiny).
  Import-path deltas from the new `src/print/` location (verified):
  `./print/shell/*` → `./shell/*`; `./print/ui` → `./ui`; `./lib/queue/*` →
  `../lib/queue/*`; `./lib/toast` → `../lib/toast`; `./styles.css` is already
  loaded by `src/print/main.tsx:1` (omit it here — no double import).
- **`src/print/main.tsx`** — three edits: (a) **mount `#app`** (not `#print-root`),
  so the `styles.css #app` height chain holds; (b) carry the boneyard
  `configureBoneyard({...})` block from `src/main.tsx:25-30`; (c) carry the registry
  import but **rewrite the path** `import "./bones/registry"` → `"../bones/registry"`
  (registry lives at `src/bones/`, not `src/print/bones/`). QueryClient / StrictMode /
  BrowserRouter / ErrorBoundary wrappers already match `src/main.tsx`.

### EDIT
- **`index.html`** — repoint `<script src="/src/main.tsx">` → `"/src/print/main.tsx"`.
  Keep `<div id="app">`. (Optionally update `<title>` to "Himalaya · The Print".)
- **`vite.config.ts`** — remove the `rollupOptions.input = { print: ... }` override
  so Vite uses the default `index.html` entry (build → `dist/index.html`); delete the
  stale "index.html left for dev reference until cutover" comment.
- **`packages/HimalayaUI/frontend/scripts/check-design.mjs:42`** — rename
  `COLOR_AUTHORING_ALLOWLIST` entry `"main.tsx"` → `"print/main.tsx"` (the carried
  `configureBoneyard` oklch literals fail `lint:design` otherwise; `print/main.tsx`
  is not covered by any `isExcluded()` prefix). **Same commit as the boneyard carry.**
- **`test/smoke.test.tsx`** — `import { App } from '../src/App'` →
  `import { PrintApp as App } from '../src/print/App'`. Becomes the real coverage of
  the production root.
- **`test/print-shell.test.tsx`** — its only assertion is the stub `print-shell`
  landmark, which vanishes with promotion. The promoted `PrintApp` calls
  `useQueryClient()` + renders `<AppRoutes>` (`<Routes>`/`useParams`), so a bare
  `render(<PrintApp/>)` throws. **DELETE it** (coverage moves to the repointed
  `smoke.test.tsx`, which already wraps in providers + `MemoryRouter`), or rewrite
  with `renderWithProviders` + a Router. Do NOT just swap the testid string.
- **Docs:** `src/AGENTS.md` (`:9` entry ref `main.tsx`→`print/main.tsx`/`PrintApp`;
  `:63` design-guard allowlist ref `main.tsx`→`print/main.tsx`),
  `docs/superpowers/plans/2026-06-09-frontend-full-migration.md` (mark the
  pages-assembly promotion DONE), `CLAUDE.md` (one-line cutover note).

### DELETE (only AFTER promotion + the build/e2e gate)
- `print.html`, `src/main.tsx`, `src/App.tsx`.
- **KEEP** `test/ErrorBoundary.test.tsx` (ErrorBoundary survives; the inventory's
  "delete" was on a false legacy-only premise).

### KEEP unchanged (red-team-confirmed — these all assume `index.html`)
- `packages/HimalayaUI/src/server.jl:70` — stays `index.html`. **No edit.**
- `packages/HimalayaUI/test/test_spa_fallback.jl` — stays `index.html`. **No edit.**
- `src/styles.css:134` — stays `#app`. **No edit.**
- All of `src/print/**`, `src/api.ts`/`queries.ts`/`state.ts`/`phases.ts`/`styles.css`/
  `ErrorBoundary.tsx`, `src/hooks/**`, `src/lib/**` (incl. **`lib/queue/**` — LIVE**).
- The entire Julia backend incl. `index_groups`, `migrate_assignments!`,
  `seed_assignment_if_absent!`, `persist_analysis!` re-attachment, `active_group_kind`.
  DO NOT delete any backend file/table (the shadowed auto-group write is a
  multi-release deprecation, not this cutover). Optional: add a disaster-recovery
  note to the `migrate_assignments!` docstring.

## Risks (post-flip)

- **P0 — shipping a stub:** if the promotion is incomplete, the build serves a blank
  page. Gate: after `npm run build`, grep `dist/assets/*.js` for a real-app marker
  (e.g. an `AppRoutes` route path / `OnboardingFlow`) AND run the e2e + backend
  static-serve suites. (Under `index.html`-canonical, dev and build both route
  through the promoted `PrintApp`, so a green e2e on the dev server + a real-app
  bundle grep together close this.)
- **P1 — `lib/queue/**` is LIVE:** never delete it or `test/queue/**`.
- **P2 — staging hygiene:** working tree has `bones/registry.ts` +
  `contact-sheet.bones.json` modified and `graphify-out/` untracked. Stage only
  cutover files **by name**. No `git add -A`.
- *(Dissolved by the flip: backend static-serve edit, `test_spa_fallback` break,
  `#print-root` height collapse, Vite dev-root failure — all moot now that
  `index.html`/`#app` stays the entry.)*

## Resolved questions / out-of-scope

- **Backend routes:** nothing removable; `/groups` + `/api/comparisons/*` are
  genuinely absent (survey `w70n20q1o`). The dead `/api/comparisons/*` *frontend*
  client cluster (api + hooks + queue mutators + registry arms, zero `print/`
  consumers) is a **post-cutover follow-up**, not this PR.
- **`active_group_kind`:** unconsumed by the frontend; KEEP (future-deprecation).
- **Pre-flight, human/ops (not codebase-answerable):** does any external CI/deploy
  script expect a build artifact by name? (We now emit `dist/index.html`, the
  conventional name — lower risk than `print.html` would have been.) Does any
  external tool parse the `/api/experiments/{id}/export` JSON? Confirm before merge.

## Folded-in fix (session review `w93ybassx`, `fix-then-ship`)

- **P1 — `addArmed`/`xDomain` leak across same-route sample nav** (`FocusPage.tsx:198-199`).
  React Router doesn't remount FocusPage on `[`/`]` steps, so the "+ Peak" arm + the
  zoom window survive a sample switch — the first click on the next sample silently
  mutates *its* peaks. Fix: one effect keyed on `activeSampleId` resetting both
  (`setAddArmed(false); setXDomain(null)`) + a navigate-then-assert regression test.
- Non-blocking (defer): `dodgeX` right-edge clamp (P3 — do NOT re-add the global
  recenter); `routes_samples.jl` rollup `ORDER BY` tie-break (P3).

## Ordered execution steps

0. **Pre-flight:** confirm `git branch --show-current` = `worktree-greenfield-ui-rebuild`;
   working tree clean except intended changes. Stage by name; no `git add -A`.
1. **Promote** `src/App.tsx` body → `PrintApp` (path deltas above; omit `styles.css`
   import; no conflict modal / `useAppState`).
2. **Reconcile `src/print/main.tsx`** (mount `#app`; carry `configureBoneyard`;
   `import "../bones/registry"`) **and rename `check-design.mjs:42`
   `"main.tsx"`→`"print/main.tsx"`** in the same commit (else `lint:design` fails).
3. **Repoint `index.html`** script → `/src/print/main.tsx`; drop the `vite.config.ts`
   `rollupOptions.input` override + stale comment.
4. **Tests:** repoint `smoke.test.tsx`; **delete** (or provider-rewrite)
   `print-shell.test.tsx`. `npm test` green (esp. `test/queue/**`).
5. **Fold in** the `addArmed`/`xDomain` reset fix + its regression test.
6. **Delete** `print.html`, `src/main.tsx`, `src/App.tsx`. Verify no stray `.js`
   shadows (`find src -name '*.js'` → empty; the build uses `--noEmit`, which can't
   emit — the hazard is an emitting `tsc`, not present here).
7. **Build gate:** `npm run build` (lint:design + `tsc --noEmit` + vite) → emits
   `dist/index.html`; grep `dist/assets` for a real-app marker (NOT the stub `<h1>`).
8. **Docs:** `vite.config.ts` comment, `src/AGENTS.md:9,63`, the frontend-migration
   plan, `CLAUDE.md`. (No `server.jl` edit under `index.html`-canonical.)
9. **e2e gate:** `npm run e2e` (mocked) — `/` resolves the app; routing/SSE/queue work
   through `PrintApp`. This is the load-bearing dev-server check.
10. **Backend suite:** `julia --project=packages/HimalayaUI` HimalayaUI tests (slow;
    capture + grep) — incl. `test_spa_fallback.jl` (still `index.html`, should pass
    untouched).
11. **Commit** as focused commits (Co-Authored-By trailer). Push.
12. **Open the PR** into `main` via `gh`: call out the `PrintApp` promotion, the
    `index.html`-canonical single-entry, the folded-in `addArmed` fix, and **merge via
    merge-commit** to preserve backend migration history.
13. **STOP.** Do not merge — branch stays unmerged pending Jonathan's go-ahead; when
    given, merge `--no-ff`.

**Rollback:** if any gate (7/9/10) goes red, the deletes (step 6) haven't been
committed yet — `git restore` the working tree (or `git reset --hard` to pre-cutover
HEAD on this unmerged branch) and re-run. Commit the deletes only after all gates are
green.
