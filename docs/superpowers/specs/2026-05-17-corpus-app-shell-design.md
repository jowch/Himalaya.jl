# Corpus app shell — hoisted routing + nav-bridge (issue #155 / I1.1)

**Status:** design — approved 2026-05-17.
**Issue:** [#155](https://github.com/) — "Add corpus app shell with hoisted routing and nav-bridge (I1.1)".
**Parent plan:** [`2026-05-17-himalaya-ui-redesign.md`](../plans/2026-05-17-himalaya-ui-redesign.md) §2.3, §4.1 · [issue map](../plans/2026-05-17-himalaya-ui-redesign-issue-map.md) I1.1, §3.

## 1. Purpose

This is the structural spine of the HimalayaUI redesign Phases 1–4. It introduces the
new corpus-scoped application shell (topbar with wordmark, three workflow stage-tabs, a
Beamtime facet chip) and the routing model the new surfaces mount into. No user-facing
surface ships here — the contact sheet (#160), loupe (#161), and focus workspace (#179)
mount *inside* this shell later. #155 delivers the empty shell and the routing/nav
machinery, nothing more.

The redesign makes every new surface URL-routed from the start. Today's navigation is
**dual**: Compare is URL-routed; Index/Inspect are driven by the Zustand `activePage`
field through `parseLocation`'s fragile three-page union. Rather than extend that union,
the migration runs both models side by side under one router until Phase 5 retires the
dual model.

## 2. Architecture — one `<Routes>`, two layout routes

A single hoisted top-level `<Routes>` lives in a dedicated `AppRoutes.tsx` (rendered by
`App.tsx`) — splitting the route table into its own file keeps it independently testable
without `App.tsx`'s queue/SSE/`EventSource` effects. It uses **react-router v6 layout
routes** (a pathless parent `<Route>` whose element renders `<Outlet/>`) so each body
keeps its own chrome:

```
App.tsx     — root concerns: queue/SSE/persistence effects; renders <AppRoutes/>
└── AppRoutes.tsx — theme effect + useGlobalShortcuts (above both bodies)
    └── <Routes>                                ← the single hoisted route table
    ├── <Route element={<CorpusShell/>}>         new shell: topbar + <Outlet/>
    │     └── <Route path="/samples" element={<SamplesPage/>}/>   ← #160 fills the stub
    └── <Route element={<AppShell/>}>            legacy shell: AppHeader + TabRocker + <Outlet/>
          ├── /  ·  /index/*  ·  /inspect/*                  → <PageBody/>
          ├── /experiments/:eid/compare/*  ·  /compare/all/* → <Compare/>  (+ /edit redirects)
          └── *                                              → <PageBody/>  (stale fallback)
```

Two page bodies under one router, exactly as the master plan §4.1 mandates.

### 2.1 `AppShell` is refactored, not replaced

`AppShell` currently *owns* an internal `<Routes>`. That route table is **hoisted out** —
its child routes move up into the one top-level table — and `AppShell` becomes a layout
component that renders `<Outlet/>`. It keeps its name and its grain-background +
`AppHeader` + `TabRocker` + `NavModal` chrome. Phase 5 (#182) deletes it; until then it
is the rollback path and must keep working unchanged for legacy URLs.

Route specificity ranking (not source order) governs matches: `/samples` exact-matches
the corpus child and always outranks the legacy `*` splat.

### 2.2 The route table is the extension point

The `<Route element={<CorpusShell/>}>` parent is where later issues register their
surface: #161 adds `/samples/loupe/:sampleId`, #179 adds `/sample/:sampleId`. #155 seeds
it with a single child — `/samples` → a `SamplesPage` stub that #160 fleshes into the
contact sheet.

## 3. The nav-bridge — structural isolation

`/samples` is **not** a member of `parseLocation`'s three-page union, and the master plan
forbids extending that union. Left alone, the legacy URL-sync hooks
(`useStateFromUrl`/`useUrlFromState`) would parse `/samples` as a stale path, call
`setStaleUnknownPath`, and strand the user on an empty page — the issue-#77 failure class.

**Fix — relocation, not a guard.** The two URL-sync hooks, plus `AppShell`'s existing
`activePage ↔ compare-URL` sync effects, move out of `App.tsx` and into `AppShell` (the
legacy layout component). They then mount **only when a legacy route matches**. On
`/samples` only `CorpusShell` mounts; the legacy machinery never runs, so it cannot fight
the new route. This is zero-internal-change to the #114/#118-scarred hooks — only their
call site moves.

*Alternative considered and rejected:* keep the hooks at root and add an
`if (isCorpusRoute(pathname)) return` early-return. Rejected — it edits the internals of
fragile hooks; the layout-route relocation achieves the same isolation structurally.

`useGlobalShortcuts` and the `<html>` theme effect move **up** into `AppRoutes.tsx`
(above the `<Routes>`, so both bodies share them) — out of `AppShell`, per master plan
§4.1.

## 4. Stale-`activePage` coercion — `state.ts`

`state.ts` gains an exported pure function:

```ts
export const VALID_PAGE_IDS = new Set<PageId>(["index", "inspect", "compare"]);
export function coerceActivePage(raw: unknown): PageId {
  return typeof raw === "string" && VALID_PAGE_IDS.has(raw as PageId)
    ? (raw as PageId)
    : "index";
}
```

It is wired through a `persist` **`merge`** callback that replicates zustand's default
shallow merge (`{ ...current, ...persisted }`) and then coerces `activePage`:

```ts
merge: (persisted, current) => {
  const merged = { ...current, ...(persisted as Partial<AppState>) };
  return { ...merged, activePage: coerceActivePage(merged.activePage) };
}
```

`merge` is a separate `persist` config field — **not a version bump**. `version` stays
`3`; no `migrate` callback is introduced (a bump without `migrate` would discard every
user's persisted state — master plan §2.3). A stale persisted `activePage` is coerced at
rehydration and never enters the store.

At #155 all three `PageId` values are still valid, so coercion is a no-op on real data —
the machinery and its test exist so that #1.7 / #3.6 / #4.4 can shrink `VALID_PAGE_IDS`
as each legacy surface is retired and have coercion fire automatically.

## 5. New components and files

| File | Role |
|---|---|
| `components/CorpusShell.tsx` | New layout route — renders `CorpusTopbar` + `<Outlet/>`. |
| `components/CorpusTopbar.tsx` | Topbar: `Himalaya · SAXS` wordmark; three stage-tabs; Beamtime facet chip. |
| `pages/SamplesPage.tsx` | Placeholder stub at `/samples`; #160 builds the contact sheet into it. |
| `components/AppRoutes.tsx` | **New.** The single hoisted `<Routes>`; theme effect + `useGlobalShortcuts` above both bodies; hosts `PageBody` + `EditToBareRedirect`. |
| `App.tsx` | Renders `<AppRoutes/>`; loses the two URL-sync hooks (moved to `AppShell`). |
| `components/AppShell.tsx` | Refactored to a layout route: internal `<Routes>` removed, renders `<Outlet/>`, hosts the relocated URL-sync hooks + compare-sync effects. |
| `state.ts` | `coerceActivePage` + `VALID_PAGE_IDS` exports; `mergePersistedState` wired as the `persist` `merge` callback. |
| `styles.css` | "The Print" palette tokens added to `@theme`. |

### 5.1 Stage-tabs during the migration window

The topbar renders all three stage-tabs (Samples / Index / Series) per the mockup. Only
**Samples is active** (→ `/samples`); **Index and Series are rendered inert** until their
phases build those surfaces. The new shell is self-contained in #155 — `/samples` is
reached by URL only (the `/inspect*` → `/samples` redirect is #1.7's, not #155's).

### 5.2 Beamtime facet chip

Rendered as a **presentational placeholder**. The issue map assigns `?beamtime=` URL
query state to #160's `/samples` route; #155 only places the chip in the topbar.

### 5.3 "The Print" palette

New palette tokens are added to the existing `@theme` block as **single-valued** tokens
(no `:root.theme-light` override). The redesign is light-first and "The Print" is
inherently a paper aesthetic, so the corpus topbar renders light in both theme modes; a
dark treatment of the new surfaces is a later, post-#155 concern.

## 6. Testing (Vitest + build)

- **`coerceActivePage` unit test** — valid ids pass through; `"loupe"`, `undefined`,
  `null`, non-strings → `"index"`.
- **Nav-bridge integration (MemoryRouter)** — at `/samples`, `CorpusShell` renders
  (wordmark, three stage-tabs, Beamtime chip present) and the legacy URL-sync hooks do
  **not** redirect away; at `/index`, `AppShell` renders; navigating by URL updates
  `activePage` and changing `activePage` updates the URL with no desync.
- **Stale-`activePage` coercion path** — rehydrating with a bogus persisted `activePage`
  yields a valid one and no empty `PageBody`.
- `npm run build` (`tsc --noEmit` + vite) passes.

## 7. Out of scope

- The `/samples` contact sheet and loupe bodies (#160, #161).
- The `/sample/:sampleId` focus-workspace route slot (#179).
- Deleting `AppShell`, `WorkspaceGrid`, `TabRocker`, or the `activePage` model — Phase 5
  (#182); they coexist until then.
- Any Zustand `persist` version bump — adding a `merge` callback and additive fields
  needs none.
- Extending `parseLocation`'s three-page union — new surfaces own their URLs via
  `useParams`/`useNavigate`.

## 8. Acceptance criteria (from #155)

- [ ] A single top-level `<Routes>` renders the legacy `AppShell` body for old paths and
      the new shell body for new paths.
- [ ] The new shell shows the corpus wordmark, three stage-tabs, and the Beamtime facet
      chip.
- [ ] `useGlobalShortcuts` and the theme effect work under both bodies.
- [ ] Navigating by URL updates `activePage`, and changing `activePage` updates the URL,
      with no desync.
- [ ] A stale persisted `activePage` naming a non-existent surface is coerced to a valid
      one and does not render an empty `PageBody`.
- [ ] No Zustand `persist` version bump is introduced.
- [ ] Vitest covers the nav-bridge and the stale-`activePage` coercion.
- [ ] `npm run build` passes.
