# Greenfield UI Rebuild — "The Print, from the foundation"

**Date:** 2026-05-31
**Status:** Design — pending review
**Scope:** HimalayaUI frontend (`packages/HimalayaUI/frontend`)

## Problem

The current frontend renders close to the spec and the canonical "The Print"
mockups, but the **implementation** is built on top of an earlier UI iteration
that has been patch-piled. The defect is not the target (the spec and mockups
are right) — it is that **two UI generations are tangled in the same tree**, and
neither a human nor an agent can reliably tell new-iteration code from
old-iteration code. The repeated failure mode: an agent, asked to change or add a
surface, copies a neighbouring component that was never migrated, and the old
iteration propagates forward. Discipline ("just don't reuse the old stuff") has
failed multiple times because the old code remains *reachable* — importable and
visually plausible — so it keeps getting reused.

The fix must be **structural and mechanical**, not a matter of judgment: the old
UI must become impossible to import from new code.

## Goal

Re-implement the entire frontend **visual layer** from the foundation up —
primitives first, then components, then pages — against the existing, trusted
mockups, while **preserving the application logic core unchanged**. The new tree
is a clean room that cannot reach the old one.

Non-goal: redesign. The mockups and `DESIGN.md` are the target as-is. This is a
clean re-implementation of a target we already trust, not a new design pass.

Non-goal: preserving usability of the current app during the rebuild. The
redesign is not in production and will not be deployed until complete. We invest
**zero effort** keeping the old app runnable or its UI tests green during the
effort.

## The seam: keep vs. rebuild

### Shared core — KEPT, imported by the new tree, never rewritten

- `src/api.ts`, `src/queries.ts` (TanStack Query server state), `src/state.ts`
  (Zustand client state), `src/phases.ts` (phase→color palette)
- `src/lib/queue/**` — the optimistic mutation queue (idempotency,
  replay-as-rerun, persistence, conflict bridge) and its six-layer contract tests
- SSE wiring (currently inline in `App.tsx`; the new `App` re-establishes the same
  `EventSource('/api/events')` → `handleRemoteEvent` plumbing)
- `src/lib/**` geometry & logic: `plot/`, `qRing`, `detectorOrient`,
  `assignment`, `scoping/`, `series/`, `comparison/`, `units`,
  `seriesRatio`, `customIndex`, `color-distance`, `clientId`, `clientOpId`,
  `authOpts`, `toast`
- Logic hooks: `useStateFromUrl`, `useGlobalShortcuts`, `useFocusTrap`,
  `useSyncActiveSampleFromRoute`, `useAutoPickExposure`, `useCurrentUserId`

**Known impurity — `src/lib/figure-export/**` (verified 2026-05-31):** unlike the
rest of `lib/**`, `figure-export` is **coupled to the renderer layer** — it
imports old renderer components (`components/ui/peakMark`,
`components/MemberHeatmapLayer`, `components/CrossTraceTrackingLayer`,
`components/RepresentationToggle`) to reuse their Observable-Plot *mark* logic. It
is therefore not pure shared core. It is **not imported by the foundation
`src/print/` tree**, so the clean room holds initially; it is repointed at the
`print/` renderer components when those are rebuilt (`peakMark` is seeded into
`src/print/ui/`, so a target already exists for it).

Phase 0 verifies the rest of the shared core is pure logic (no import of
`src/components/**` or `src/pages/**`); the `figure-export` coupling is the one
known, documented exception. Any *other* accidental entanglement is extracted into
the core during Phase 0, not carried forward.

### Greenfield — REBUILT clean from the mockups

- All consumer components (~84 files in `src/components/*.tsx`)
- All pages (`src/pages/**`)
- The skeleton / `bones` definitions (rebuilt for the new component shapes)
- A **fresh copy** of the `ui/` primitives — seeded from today's clean
  `src/components/ui/` extraction (the 2026-05-29 work is genuinely clean and is
  the right starting point), then honed and **expanded** to cover every reusable
  pattern extracted from the mockups

### Renderers — REBUILT visually, math REUSED

`TraceViewer`, `MillerPlot`, `MultiTracePlot`, `DetectorImage` are visual but
contain real geometry/canvas logic. They are rebuilt visually in the new tree but
**import the existing geometry helpers** (`lib/plot`, `qRing`, `detectorOrient`)
rather than re-deriving any SAXS math. The recent Plotting-redesign (#280) is
trusted; we do not touch the plotting math.

## Architecture

### Greenfield tree location + mount: separate Vite entry

New app lives under **`src/print/`**:

```
src/print/
  main.tsx        # mirrors src/main.tsx: StrictMode > ErrorBoundary >
                  # QueryClientProvider > BrowserRouter > App, boneyard config
  App.tsx         # root: routes + SSE wiring + queue persistence/rehydrate
  ui/             # honed + expanded design-system primitives (closed-look)
  components/     # rebuilt consumer components
  pages/          # rebuilt pages (one per mockup surface)
  bones/          # new skeleton captures for the new component shapes
```

Mounted via a **separate Vite entry** — a new `print.html` (→ `/src/print/main.tsx`)
added to `build.rollupOptions.input` alongside the existing `index.html`. The old
`index.html` / `src/main.tsx` is left untouched and unmaintained; it may rot and
is deleted at cutover. Dev runs the new app at `print.html`; old app is not a
concern.

The new tree imports the shared core by its existing paths
(`../api`, `../lib/...`, `../state`, etc.). It never imports `src/components/**`
or `src/pages/**`.

### Enforcement spine — the import-boundary guard

The structural guarantee that fixes the root cause. Extend the existing
`scripts/check-design.mjs` (already run as the `lint:design` build step) with a
new rule:

- **Rule: no-legacy-import.** Any file under `src/print/**` that imports
  (relatively) from `src/components/**` or `src/pages/**` is a hard error
  (exit 2). Old UI is mechanically unreachable from new UI by **direct** import;
  the build fails if anyone reaches back.
- **Known limitation — transitive reach.** The guard checks direct imports only.
  A `print/` file importing a shared-core module that *itself* imports old UI
  (today: `lib/figure-export` → `components/**`) would reach old code
  transitively without tripping the guard. Mitigation: `print/` does not import
  `figure-export` until the renderers are rebuilt, at which point `figure-export`
  is repointed at the `print/` renderer components — closing the leak. Tracked in
  the Phase 0 core-purity ledger.

Additionally:

- The design guard's `ui/`-authoring exclusion is extended to cover
  `src/print/ui/**` (appearance authored there; `src/print/**` consumers are
  placement-only, same closed-look contract as today).
- The color-authoring allowlist in `check-design.mjs` (currently naming
  `components/DetectorImage.tsx`, `components/FocusDetectorPanel.tsx`,
  `components/MemberHeatmapLayer.tsx`) gains `src/print/` equivalents as those
  renderers are rebuilt.

### Old code and the build

The production build emits **only** the new app (`vite.config.ts`
`rollupOptions.input` → `print.html`); the old `index.html` entry is not built.

Old code is **not** excluded from the `tsc --noEmit` typecheck in early phases,
for two reasons (verified 2026-05-31): (a) it is unnecessary while old still
typechecks green and we are only *adding* `src/print/`; (b) a `tsconfig`
`exclude` is *ineffective* here anyway — TS still typechecks an excluded file if
an *included* file imports it, and `lib/figure-export` (included shared core)
imports old `components/**`, dragging them back in. A real typecheck-time
exclusion of old becomes relevant only once old rots, and must be paired with
severing/repointing `figure-export`'s component imports — a later-phase task, not
foundation work. Old code remains git-tracked as reference until cutover deletes
it.

### Visual verification: Storybook + Playwright

Two complementary harnesses, both built fresh (the stranded `e2e/fidelity`
harness in the `audit-followups` worktree is **not** reused):

- **Storybook** (latest, `@storybook/react-vite` builder; exact version pinned
  during planning) — the primary visual-verification surface for **component**
  fidelity and behavior in isolation. Stories colocate with components
  (`src/print/ui/Button.stories.tsx`, etc.). `styles.css` (`@theme` tokens) is
  imported into the Storybook preview so primitives render with the real Print
  palette/typography. Each story can sit beside its mockup specimen for
  side-by-side fidelity checking. A component that only renders in Storybook
  against mock props cannot reach for old app code or hidden global state — this
  reinforces the clean room.
- **Playwright page-fidelity** (fresh, thin) — **integration** fidelity: compares
  assembled `print.html` pages against mockup renders in Phase 3.

The defended green-bar during the rebuild is the **logic core** (six-layer
contract tests, queue, idempotency, logic/unit tests), not the old UI tests.
Old component Vitest tests and old Playwright e2e are not maintained; they are
deleted with their components at cutover. New components get fresh tests + stories.

## Mockup surfaces → new pages

Mockups live in `docs/redesign-mockups/`:

| Mockup | New surface |
|---|---|
| `sample-table.html` | corpus / contact-sheet page |
| `sample-table.html` (the loupe/detail view is mocked *within* this file) | Loupe page — remains its **own page/route** (own file under `src/print/pages/`); its visual design is sourced from the loupe section of the `sample-table` mockup, not a dedicated mockup file |
| `focus-workspace.html` + `2026-05-29-focus-plot.html` | Focus workspace page + Focus plot renderer |
| `series-folio.html` | Series folio page |
| `series-scoping.html` | Series scoping page |
| `series-builder.html` + `2026-05-29-series-plot.html` | Series builder page + Series plot renderer |

The Loupe stays a standalone page. The only clarification is that its design
reference lives *inside* `sample-table.html` rather than in a separate mockup
file — there is no `loupe.html`.

## Phases

Bottom-up: foundation first, then outward. Each phase becomes its own
implementation plan; Phase 0's output feeds all later phases.

- **Phase 0 — Inventory + mockup extraction.** Two ledgers:
  - *Core-purity ledger:* confirm every shared-core file is pure logic
    (no component/page imports); list the exact geometry functions the rebuilt
    renderers will import; flag and extract any logic accidentally entangled with
    old components.
  - *Design-system catalog + surface map:* walk all 7 mockups and `DESIGN.md`;
    extract **every** reusable primitive/pattern into a complete catalog (the
    `src/print/ui/` target set, seeded from today's primitives and expanded); map
    each mockup to its new page + the components it needs.
- **Phase 1 — Foundation.** Stand up: separate Vite entry (`print.html`, build
  output only), the import-boundary guard, and Storybook. Seed `src/print/ui/`
  from today's primitives, then hone + expand to the full catalog, each primitive
  verified in Storybook against its mockup specimen. (No `tsconfig` exclusion of
  old — see "Old code and the build".)
- **Phase 2 — Components.** Rebuild consumer components surface by surface from
  the primitives, each verified in Storybook. Renderers rebuilt here, importing
  the trusted geometry math.
- **Phase 3 — Pages + integration.** Assemble pages in `src/print/pages/`, wire
  routing + SSE + the mutation queue in `src/print/App.tsx`, run against real
  data via `print.html`. Playwright page-fidelity checks against mockups.
- **Phase 4 — Cutover.** Repoint `index.html`/entry to the new app; delete old
  `src/components/**` (consumer), `src/pages/**`, old `src/components/ui/**`, old
  `App.tsx`/`AppRoutes`, old bones; ensure `lib/figure-export` now imports only
  `print/` renderer components (the last component coupling); remove the obsolete
  `print.html`/`index.html` split; flatten `src/print/**` → `src/**`; remove the
  no-legacy-import rule (or repoint it); delete dead old UI tests; full suite
  green.

## Verification spine

- **Throughout:** logic-core tests stay green (`lib/queue`, idempotency,
  six-layer contract tests, logic unit tests). `npm run build` for the new entry
  passes (`lint:design` incl. the import-boundary rule + design guard, `tsc`,
  `vite build`).
- **Per component (Phases 1–2):** Storybook story renders all variants/states;
  side-by-side mockup-specimen check.
- **Per page (Phase 3):** Playwright page-fidelity against the mockup.
- **At cutover (Phase 4):** full `npm run build` + `npm test` + `e2e` green on the
  flattened tree; backend Julia suite unaffected (no backend changes).

## Risks & mitigations

- **Big-bang cutover risk.** Mitigated: the new app is continuously runnable via
  `print.html` and per-component/per-page fidelity is verified incrementally in
  Storybook + Playwright. Only the *prod swap* is deferred to one cutover — and
  production was never serving the redesign anyway.
- **Duplicated primitives during transition.** Acceptable and temporary; the old
  `ui/` is deleted at cutover.
- **Renderer behavior regressions.** Renderers import unchanged geometry math; the
  risk is confined to the visual shell and is caught by Storybook (component) +
  Playwright (integration) checks.
- **Shared-core entanglement discovered in Phase 0.** If a "logic" file secretly
  imports an old component, it is extracted into the pure core during Phase 0
  before any rebuild depends on it.

## Open items for planning

- Exact Storybook version + Tailwind v4 preview wiring (verify via context7).
- Whether Playwright page-fidelity uses pixel snapshots vs. structural assertions
  (decide per surface; the mockups are static HTML, so rendered-mockup vs.
  rendered-page screenshot diffing is feasible).
- The flatten-at-cutover mechanics (alias vs. physical move) — a Phase 4 detail.
