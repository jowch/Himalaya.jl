# Greenfield → main cutover plan

**Date:** 2026-06-14
**Branch:** `greenfield-ui-rebuild` (stays UNMERGED until Jonathan's explicit go-ahead)
**Status:** PLAN — not yet executed. Produced from the `greenfield-cutover-inventory`
review workflow; the load-bearing claims were independently verified against live
source (see "Evidence").

---

## TL;DR

The cutover is **not** "delete dead code." It is a **promotion plus a production
correctness fix**, landed as a **merge-commit**:

1. The production build (`print.html → src/print/main.tsx → PrintApp`) currently
   compiles a **15-line placeholder stub**. The real app lives in `src/App.tsx`,
   reachable only from the *un-built* legacy `index.html → src/main.tsx` entry.
2. Therefore `npm run build` ships a blank page, and `server.jl` already serves a
   `dist/index.html` the build never emits. **Merging as-is would ship a stub to
   production.** Everything works today only because dev/e2e run Vite, which
   serves `index.html` (the real app) by path.
3. The fix is to **move the real app body into `PrintApp`**, point the backend at
   `print.html`, then delete the now-truly-dead legacy entry — with the build +
   e2e gate proving the bundle contains the real app, not the stub.

## Evidence (verified this session)

| Claim | Status | Evidence |
|---|---|---|
| `PrintApp` (`src/print/App.tsx`) is a 15-line stub | **VERIFIED** | renders only `<main data-testid="print-shell"><h1>The Print · greenfield shell</h1></main>` |
| `src/App.tsx` is the real app, legacy-entry-only | **VERIFIED** | imports `AppRoutes` + `EventSource('/api/events')` + queue `attachPersistence`/`rehydrate` + shell siblings; reachable only via `index.html → src/main.tsx` |
| `server.jl:70` serves `dist/index.html`; build emits only `dist/print.html` | **VERIFIED** | `dist/` contains `print.html` + `assets/` only; `server.jl:70` reads `joinpath(dist_dir, "index.html")` |
| `lib/queue/**` is LIVE (mislabeled "dead" by reachability-from-stub) | **VERIFIED** | imported by `src/App.tsx` (the real app) and ~15 `test/queue/**` specs |
| Backend `index_groups`/`migrate_assignments!` are migration-load-bearing | per inventory (high-confidence) | not personally re-verified; cutover does **not** touch them — verify before any future removal |

## Decisions

*(All confirmed by Jonathan, 2026-06-14.)*

- **Merge strategy: merge-commit (`--no-ff`).** The branch evolved both the
  frontend and the Julia backend/data-model; the backend migration history
  (`confirmed_index` semantic shift, rebuild-not-log-derivable pre-Plan-A
  confirmations) is disaster-recovery-load-bearing and must stay bisectable.
  Squash would destroy it; a fresh-from-main branch would discard real history.
- **Entry strategy: single entry (`print.html` canonical).** Promote `PrintApp`,
  delete the legacy `index.html`/`src/main.tsx`/`src/App.tsx`, and resolve the
  Vite dev-root via config so the dual-entry trap can't recur. (If the e2e
  dev-root gate fails on the missing `index.html`, fall back to a thin
  `index.html` that re-exports `src/print/main.tsx` — but single-entry is the goal.)

## Manifest

### PROMOTE (the load-bearing edit)
- `src/print/App.tsx` — replace the stub `PrintApp` body with `src/App.tsx`'s body:
  `AppRoutes` + `OnboardingFlow` + `NavModal` + `ToastContainer`/`LiveRegion` +
  `InfrastructureBanner`, plus the 4 effects (DEV-gated `exposeTestHelpers`;
  `EventSource('/api/events')` `curation` → `handleRemoteEvent`;
  `attachPersistence`; `rehydrate` with dropped/failed toasts). Keep the export
  name `PrintApp`. Fix relative import paths for the new location
  (`./print/shell` → `./shell`, `./lib/queue` → `../lib/queue`, etc.). Drop the
  "Foundation placeholder" docstring.
- `src/print/main.tsx` — carry the boneyard `import './bones/registry'` +
  `configureBoneyard({...})` block from `src/main.tsx:4,25-30` (QueryClient +
  ErrorBoundary + `#print-root` already match). **⚠ verified build-gate coupling:**
  that block authors raw `oklch()` literals, and `check-design.mjs:42` allowlists
  bare `"main.tsx"` while `isExcluded()` does NOT exempt `print/main.tsx` — so the
  carried block fails `lint:design` (which gates `npm run build`) unless the
  allowlist entry is renamed (see EDIT → `check-design.mjs`).

### EDIT
- `packages/HimalayaUI/src/server.jl:70` — `index.html` → `print.html` (**P0 prod
  fix**; the SPA fallback throws today on the absent file).
- `packages/HimalayaUI/frontend/scripts/check-design.mjs:42` — rename the
  `COLOR_AUTHORING_ALLOWLIST` entry `"main.tsx"` → `"print/main.tsx"` so the carried
  `configureBoneyard` oklch literals don't fail `lint:design`. **Required** for the
  build gate to pass after the boneyard carry (survey-found, verified).
- `test/smoke.test.tsx` — `import { App } from '../src/App'` →
  `import { PrintApp as App } from '../src/print/App'`. Becomes the real coverage
  of the production root.
- `test/print-shell.test.tsx` — currently asserts the *stub* `print-shell`
  landmark; re-point to the real shell or fold into `smoke.test.tsx`.
- `vite.config.ts` — delete the now-false "index.html left for dev reference"
  comment; note `print.html` is the sole entry. `rollupOptions.input` unchanged.
- `playwright.config.ts` — confirm Vite serves `print.html` at `/` after the
  delete (dev-root gate).
- `src/AGENTS.md`, `docs/superpowers/plans/2026-06-09-frontend-full-migration.md`,
  `CLAUDE.md` — reflect `PrintApp` as the sole production root.

### DELETE (only AFTER the promotion + gate)
- `index.html`, `src/main.tsx`, `src/App.tsx`.
- `test/ErrorBoundary.test.tsx` — **KEEP** (ErrorBoundary survives; the map's
  "delete" was on the false legacy-only premise).

### KEEP (unconditional)
- All of `src/print/**` (the greenfield UI, incl. `shell/AppRoutes.tsx`).
- `src/api.ts`, `src/queries.ts`, `src/state.ts`, `src/phases.ts`, `src/styles.css`,
  `src/ErrorBoundary.tsx`, all of `src/hooks/**` and `src/lib/**`
  (**incl. `lib/queue/**` — LIVE**).
- The entire Julia backend, incl. `index_groups` schema, `migrate_assignments!`,
  `seed_assignment_if_absent!`, `persist_analysis!` re-attachment. The only
  backend edit is `server.jl:70`. Add clarifying comments to
  `migrate_assignments!` (disaster-recovery warning) but DO NOT delete any
  backend file or table — removing the shadowed `persist_analysis!` auto-group
  write is a multi-release deprecation, not this cutover.

## Risks

- **P0 — shipping a stub:** merging/building without the promotion serves a blank
  page; e2e won't catch it (it runs Vite dev). Gate: confirm the built bundle
  contains real-app strings / its size jumps.
- **P0 — backend static serve:** `server.jl:70` reads a non-existent
  `index.html`; deep-link/refresh on any route 500s in prod today. The
  `print.html` edit is mandatory.
- **P1 — Vite dev-root after `index.html` delete:** `/` may not resolve a document.
  Gate on `npm run e2e`; mitigate with a thin `index.html` re-export if needed.
- **P1 — losing SSE/queue:** do NOT delete `lib/queue/**` or `test/queue/**`
  (LIVE; mislabeled by reachability-from-stub).
- **P2 — boneyard drift:** reconcile the `bones/registry` import + defaults when
  deleting `src/main.tsx`, or skeleton loading regresses.
- **P2 — staging hygiene:** working tree has `bones/registry.ts`,
  `contact-sheet.bones.json` modified and `graphify-out/` untracked. Stage only
  cutover files by name. No `git add -A` (standing mandate).

## Open questions

1. ~~**Confirm the promotion is intended.**~~ **RESOLVED (Jonathan, 2026-06-14):**
   the promotion was simply never finished — no special reason. The Print IS the
   app; proceed.
2. ~~**Boneyard in `print/main.tsx`:** carry or drop?~~ **RESOLVED:** carry the
   registry import + `configureBoneyard` defaults (skeleton loading is live).
3. ~~**Backend dead routes**~~ **RESOLVED (survey `w70n20q1o`):** NO backend route
   is removable — all 12 `register_*_routes!()` are wired; the retired `/groups`
   and the entire `/api/comparisons/*` surface are **genuinely absent** (no handler,
   no `routes_comparisons.jl`), not orphaned. Only doc residue (stale comment at
   `routes_analysis.jl:206`). The removable residue is **frontend-side** (`api.ts`
   `/api/comparisons/*` client functions, zero `print/` consumers) — a connected
   cluster (api + queries hooks + queue mutators + `mutatorRegistry` arms),
   **out of cutover scope → follow-up cleanup.**
4. ~~**`active_group_kind` consumption**~~ **RESOLVED:** NOT consumed by any
   frontend code (`routes_export.jl:28,50,91` has zero frontend callers; the export
   route itself is unused by the UI). KEEP the field/route — future-deprecation
   candidate, not this cutover.

## Survey findings (run `w70n20q1o`, 2026-06-14) — `promotion_ready: true`

Full evidence-cited study: `docs/superpowers/notes/2026-06-14-greenfield-cutover-survey.md`.
The promotion is mechanically specified end-to-end (13 imports + 4 effects + the
6-element JSX fragment; exact import-path deltas; `#print-root` already aligned;
do NOT re-add a conflict modal/bridge or a `useAppState` import to App.tsx — that
claim broke under scrutiny). New items folded in above: the `check-design.mjs`
allowlist rename (verified build-gate fix). Standing follow-ups / pre-flight:

- **Follow-up (post-cutover):** excise the frontend `/api/comparisons/*` orphan
  cluster (api + hooks + mutators + registry arms). Not this PR.
- **Pre-flight, human/ops (not codebase-answerable):** does any external CI/deploy
  script expect `dist/index.html` by name? Does any external tool parse the
  `/api/experiments/{id}/export` JSON (the `active_group_kind` field)? Confirm
  before merge.
- **Gate-decided:** the Vite dev-root after `index.html` deletion (thin-`index.html`
  re-export fallback only if `npm run e2e` fails on `/`); `print/App.tsx` styles.css
  import redundancy (the `main.tsx`-level import already covers it).

## Folded-in fix (from the session code review, `fix-then-ship`)

- **P1 — `addArmed`/`xDomain` leak across same-route sample navigation**
  (`FocusPage.tsx:198-199`). React Router does not remount FocusPage on `[`/`]`
  steps, so the "+ Peak" arm and the zoom window survive a sample switch — the
  first click on the next sample silently mutates *its* peaks. Fix: one effect
  keyed on `activeSampleId` resetting both (`setAddArmed(false); setXDomain(null)`),
  plus a navigate-then-assert regression test. (Same "page state survives
  same-route nav" family as the boneyard fix.) Land on the cutover branch.
- Non-blocking from the same review (defer or fold opportunistically): `dodgeX`
  right-edge clamp (P3, do NOT re-add the global recenter — minimal dodging was a
  deliberate request); `routes_samples.jl` rollup `ORDER BY` tie-break to match
  `comparisons.jl` (P3, determinism).

## Ordered execution steps

0. Pre-flight: working tree clean except intended cutover changes; no `git add -A`.
1. **Promote** `src/App.tsx` body → `PrintApp`; fix import paths; keep `PrintApp` name.
2. **Reconcile** `src/print/main.tsx` (carry boneyard registry import +
   `configureBoneyard` block) **and rename `check-design.mjs:42` allowlist
   `"main.tsx"` → `"print/main.tsx"`** in the same commit (else `lint:design` fails).
3. **Tests:** repoint `smoke.test.tsx` + `print-shell.test.tsx`; `npm test` green
   (esp. `test/queue/**`).
4. **Fold in** the `addArmed`/`xDomain` reset fix + its regression test.
5. **Delete** `src/App.tsx`, `src/main.tsx`, `index.html`. Clean any stray `.js`
   shadows (`tsc --noEmit -p tsconfig.build.json`, not an emitting build).
6. **Build gate:** `npm run build` — confirm `dist/print.html` + a bundle that now
   contains the **real app** (size jump / AppRoutes-derived string), not the stub.
7. **Backend:** `server.jl:70` `index.html` → `print.html`; add the
   `migrate_assignments!` disaster-recovery docstring note. No table/file deletes.
8. **Docs:** `vite.config.ts` comment, `src/AGENTS.md`, the frontend-migration
   plan, `CLAUDE.md`.
9. **e2e gate:** `npm run e2e` (mocked) — the load-bearing check that `/` still
   resolves a working app after the `index.html` delete AND that routing/SSE/queue
   work through `PrintApp`.
10. **Backend suite:** `julia --project=packages/HimalayaUI` HimalayaUI tests (slow;
    capture + grep) — confirm the static-serve change is sound.
11. **Commit** promotion + deletes + backend edit as focused commits (Co-Authored-By
    trailer). Push.
12. **Open the PR** into `main` via `gh`: call out (a) frontend+backend cutover,
    (b) the `PrintApp` promotion, (c) the `server.jl` prod-serve fix, (d) merge via
    **merge-commit** to preserve backend migration history.
13. **STOP.** Do not merge — branch stays unmerged pending Jonathan's go-ahead;
    when given, merge `--no-ff`.
