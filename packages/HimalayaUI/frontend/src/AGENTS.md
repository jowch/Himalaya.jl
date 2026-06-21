# packages/HimalayaUI/frontend/src — Frontend

React 18 + Vite + TypeScript strict + TailwindCSS 4. TanStack Query for server state, Zustand for client state.

## Where to look

| Task | Location | Notes |
|------|----------|-------|
| App entry | `print/main.tsx` | `index.html → print/main.tsx`; `StrictMode > ErrorBoundary > QueryClientProvider > BrowserRouter > PrintApp`; mounts `#app` |
| App shell | `print/App.tsx` (`PrintApp`) | Composition root: `AppRoutes` + SSE + mutation-queue effects + shell siblings |
| Server state | `queries.ts` | TanStack Query hooks; `authOpts(username, clientId, clientOpId)` helper lives in `lib/authOpts.ts` |
| API layer | `api.ts` | Fetch wrappers; AuthOpts per-call for mutations |
| Client state | `state.ts` | Zustand store — **use named actions only** |
| Phase palette | `phases.ts` | phase → color mapping |
| Tailwind theme | `styles.css` | `@theme { --color-* … }` |
| App shell + routing | `print/shell/` | `CorpusShell`, `CorpusTopbar`, `AppRoutes`, `IndexSlugRedirect`, `StaleUrlPage`, `ResolvingFallback`, `OnboardingFlow`, `NavModal`, `InfrastructureBanner`. See [print/shell/AGENTS.md](print/shell/AGENTS.md) |
| Composites | `print/components/` | Page-composing components (rails, plates, panels, rows, modals) built from the `ui/` primitives |
| UI primitives | `print/ui/` | Closed-look design-system primitives (Button, Card, SegmentedControl, PhaseChip, PhaseStrip, ModalShell, Kicker, IconButton, ScoreBar, Dot, ToastContainer, HintText, …). Appearance lives here; consumer `className` is placement-only. See "Design system" below. |
| Render layers | `print/{plot,detector,comb,export}/` (and `print/waterfall/`) | Appearance-authoring render layers (trace-plot engine, detector image, comb/residual, waterfall, the `cleanFigureSvg` figure builder). The `lint:design` appearance guard excludes `print/{plot,detector,comb,export}/` only — `print/waterfall/` is NOT exempt |
| Pages | `print/pages/` | `SamplesPage`, `LoupePage`, `FocusPage`, `SeriesFolioPage`, `SeriesScopingPage`, `SeriesBuilderPage` (all under the single `CorpusShell`; legacy Index/Inspect/Compare pages + `AppShell` retired) |
| Hooks | `hooks/` | `useFocusTrap`, `useGlobalShortcuts`, `useStateFromUrl`, … |
| Library | `lib/` | URL helpers, plot helpers, comparison helpers, figure export |
| Mutation queue | `lib/queue/` | See [lib/queue/AGENTS.md](lib/queue/AGENTS.md) |
| Skeletons | `bones/` | Committed `*.bones.json` + auto-generated `registry.ts` |

## TypeScript strict + `exactOptionalPropertyTypes`

`set({ username: undefined })` fails — optional fields declared as `string | undefined` (rather than `username?: string`) keep this ergonomic. For passing optional values through (e.g. `AuthOpts`), use the `authOpts(username, clientId, clientOpId)` helper in `lib/authOpts.ts`, which omits any undefined key (never `{ username: undefined }`).

## Zustand — named actions

Use the store's named actions (`clearUsername`, `setTutorialSeen`, `openNavModal`, …). Avoid `useAppState.setState({ ... })` — direct setState bypasses encapsulation and triggers lint warnings. New state transitions go in `state.ts` as named actions.

## Zustand persist — version-bump wipe-guard

Two invisible invariants in `state.ts` (`persist` config, `queries.ts`'s sibling — see `state.ts:494-530`). Adding a persisted field + naively bumping `version` silently wipes **every existing user's** state on first load:

- **Bumping `version` WITHOUT a `migrate` discards the entire persisted blob.** Zustand v4 drops a blob whose version it can't migrate — all prefs (username/tutorialSeen/active*/…) gone on upgrade. Every version bump must ship a `migrate`. The current migrate's only job across all prior versions is unconditional dead-key strips (`activePage` I5.2, `theme` R0a), so there's no `switch (version)`.
- **`migrate` must return the ORIGINAL object on malformed/non-object input — never `{}`.** Handing `{}` is itself a partial wipe; returning the original lets `merge` fold whatever survived (or fall back to defaults). The wipe-guard is `if (persisted && typeof persisted === "object") { …strip… } return persisted` — the `else` returns untouched.

## State split (load-bearing)

- **Zustand owns *client* state**: active sample/exposure, hoveredIndexId, username.
- **TanStack Query owns *server* state**: experiments, samples, exposures, peaks, indices, assignment.
- Mutations invalidate scoped query keys (`queryKeys.peaks(id)`, `queryKeys.indices(id)`). Don't mix the two concerns in the same hook.

## Route → store → auto-pick cascade (the exposure-clobber hazard)

`setActiveSample` cascades `activeExposureId: undefined`, so a careless call resets the user's deliberate exposure pick. Two cooperating hooks keep this in check — touching either requires understanding the whole route-hook → Zustand → auto-pick chain:

- **`useSyncActiveSampleFromRoute` calls `setActiveSample` ONLY when the parsed id differs from the stored one** (`useSyncActiveSampleFromRoute.ts:67-73`). Drop the `parsed === activeSampleId` no-op guard and every render re-seeds the sample → clears the exposure → `useAutoPickExposure` re-adopts the representative, discarding an in-session deliberate switch (the FO-NAV-SKELETON / R5 representative-switcher fix).
- **`useAutoPickExposure` writes only when the rule picks a DIFFERENT, defined exposure** (`useAutoPickExposure.ts:85-90`). Its rule (`resolveActiveExposure`) KEEPS a still-valid `currentId` and is pure, so the page resolves the same value during render — that's what paints a cached sample's trace in the same frame instead of flashing a skeleton through the store→effect seam.

## Per-tab SSE identity

SSE self-echo filtering uses a per-tab `client_id` minted into `sessionStorage` on first load (see `lib/clientId.ts`). Audit identity (`actor` / `X-Username`) is unchanged. **Two tabs of the same user are distinct subscribers** — edits in one refresh the other. `client_id` lives for the tab session: survives reload, scoped to one tab. See `docs/event-log.md` §"Client side".

## Imperative renderers in effects

Wrap any function that is both defined inside a component AND used as a `useEffect` dependency in `useCallback` with its true deps. The effect then depends on `[theCallback]` alone — no redundant dep list, no eslint-disable. The trace plot's overlay renderer follows this pattern.

**Async re-entry guards need a ref, not state.** `useFigureExport` blocks double-invocation with a synchronous `inFlight` ref (`if (inFlight.current) return` before `setPending(true)`), *not* the `pending` state alone — state flips asynchronously, so two clicks in one tick both pass a state-only guard. The ref blocks re-entry; the state drives the UI. Same shape for any async-action hook.

## Tailwind v4 theming

"The Print" palette (light warm paper, terracotta accent hue 38, Newsreader serif) is defined once in `styles.css` via `@theme { --color-* ... }` — a single identity, no theme toggle (R0a, #221). The legacy dark-era neutral-ramp shim (the duplicated `bg-`/`text-`/`border-` neutral utilities that mirrored the canonical Print names) was excised in R3-F (#259): use the canonical Print utilities directly — `bg-paper`/`bg-paper-sunk`/`bg-plate`, `text-ink`/`text-ink-soft`/`text-ink-faint`, `border-hair`/`border-hair-strong`, `text-print-accent` (or `accent` for terracotta). Reintroducing an old name is self-revealing: Tailwind won't generate the utility and its `--color-*` custom property no longer resolves. Serif titles use the `text-display`/`text-headline` roles. To add a new color, add it to `@theme` first.

## Design system — closed-look primitives + the design guard (ENFORCED)

The Print's recurring patterns live as **closed-look** primitives in `print/ui/` (Button, Card, SegmentedControl, PhaseChip, PhaseStrip, ModalShell, Kicker, IconButton, ScoreBar, Dot, ToastContainer, HintText). They own their appearance via semantic props (`variant` / `size` / `tone` / domain props); a consumer's `className` is **placement-only** (margins, position, grid). To change how a primitive looks, build a variant *into the primitive* — the idiom is a `Record<Variant,string>` map + a tiny local `cx()` join helper (no cva/clsx/tailwind-merge). Don't restyle from the outside.

This is **mechanically enforced** (2026-05-29 extraction). `scripts/check-design.mjs` runs as a pure-absolute `lint:design` step prepended to `npm run build` (plus a warn-only PostToolUse hook), and **fails the build** on any inline appearance utility *outside* `print/ui/**`:

- arbitrary type size `text-[…]` → use a named scale role (`text-xs/sm/base/lg/xl/headline-lg/display`)
- arbitrary radius `rounded-[…]` → `rounded-sm` / `rounded-md` (both 5px) / `rounded-full`
- raw colour literal (`oklch(` / `rgba(` / quoted `#hex`) → a `--color-*` token utility
- side-stripe `border-l/r` > 1px → a full border + a leading icon/word instead

Only the colour-AUTHORING files are exempt (rules #3/#5 share an allowlist: `phases.ts`, `lib/comparison/coloring.ts`, `lib/figure-export/**`, the `print/{plot,detector,comb,export}/` render-layer prefixes, `print/main.tsx`). Note `print/waterfall/` is NOT among the exempt prefixes — its appearance (line colour, bead glyphs, axis strokes) lives inside the already-exempt `print/plot/` layer it composes, so a raw SVG colour literal surfacing in `print/waterfall/` must move into `print/plot/` or become a `--color-*` token. Need a colour anywhere else → add a `--color-*` token to `@theme`, then use the utility. Visual reference: `docs/design-system.html`; full system: root `DESIGN.md`.

`check-design.mjs` also enforces a **`no-legacy-import` rule** (`scanLegacyImports`) alongside the appearance guard: it fails the build if any `src/print/**` file relatively imports from top-level `src/components/**` or `src/pages/**`. Those legacy dirs are gone post-cutover, so the rule now functions as a regression tripwire against re-introducing the retired tree.

## Skeleton loading via boneyard-js

Each load-gated card wraps content in `<Skeleton>` from `boneyard-js/react`. Full reference: `packages/HimalayaUI/docs/boneyard.md`. Two rules that bite hardest:

- **Gate on `query.isLoading`, not `isPending`.** `isLoading = isPending && isFetching` — disabled queries and background refetches stay skeleton-free; only true cold fetches animate. Wrong gating ⇒ flicker on every refetch.
- **`className` on `<Skeleton>` is load-bearing.** Boneyard adds two wrapper divs that break parent flex chains (e.g. a scrolling list inside a card collapsing to a fixed small height). Pass `flex-1 min-h-0 flex flex-col` (or `h-full w-full`) to inherit the original child's layout role. Companion CSS in `styles.css`: `[data-boneyard-content] { display: contents }`.

## Mutation queue — load-bearing one-liners

Full architecture in `docs/mutation-queue.md`; queue internals in `lib/queue/AGENTS.md`. The invariants that bite UI code outside the queue:

- **Optimistic placeholder ids are NEGATIVE.** `Peak.id < 0` means "not yet confirmed by server"; SSE confirmation overwrites with the positive server id. UI code that filters or compares peak ids must handle negatives.
- **`useExposureHasPendingPeakOps` gates any UI that reads `peaks(id)` derivatively** while a peak op is in flight (e.g. useSpeculativeSnap). Without it: flicker as optimistic / HTTP / SSE land out of order.

## Query layer gotchas

- **`queryKeys` namespace isolation — single-entity keys must NOT share the collection prefix** (`queries.ts:73-78`). Collection keys (`["exposure", id, "peaks"]`, `…, "indices"`, `…, "assignment"`) and single-entity keys (`["peak-entity", id]`, `["index-entity", id]`, `["exposure-entity", id]`, …) are disjoint on purpose. Mutators fire targeted `invalidateQueries(["exposure", id])`, which prefix-matches and blasts everything under that root — so the `-entity` suffix is the **only mechanical barrier** keeping a per-entity query from being over-invalidated. A new per-entity query placed under `["exposure", …]` silently over-invalidates. (Same discipline elsewhere: `["series", id]` detail vs `["series-list"]` listing are distinct roots.)

- **`useStableQueryMap`'s WeakMap-nonce — do not "simplify" it** (`queries.ts:241-281`). It's the only thing stopping every parent re-render from minting a fresh `Map<id, T>`, which tears down + replots the multi-trace plot on every keystroke — wheel/brush smoothness in `MultiTracePlot` is the casualty. It can't be replaced by a `JSON.stringify` signature (a variable-length deps array breaks `useMemo` memoisation when the previous length was 0 — React's elementwise loop terminates early and returns the stale value) nor by deep equality (it tracks REF identity, since TanStack reuses `q.data` refs when nothing changed; a per-hook-instance WeakMap assigns each fresh object a stable nonce). The regression is **invisible to Vitest** — verify by hand. The WeakMap write is idempotent, so it's StrictMode- and Concurrent-safe.

## Multi-layer contract testing

Every reconciliation contract has six layers (route emit → SSE payload → `applyRemoteToCache` merge → cache row → `onMutate` → `onSuccess`). When fixing a bug at one layer, add a regression row at every other layer where the same class can manifest. See `docs/contract-testing.md` for canonical paired test files (`cache-shape.test.ts`, `sseEventPayload.contract.test.ts`, `rollbackSymmetry.test.ts`, `authHeaders.test.ts`, `test_route_response_shapes.jl`, `test_idempotency_replay_invariant.jl`).

## Samples contact sheet (`SheetTable`)

`SheetTable` renders a plain `<table>` (`role="table"`, body cells `role="cell"`). The data-grid keyboard system (`src/lib/grid/`) was removed in the ingestion-redesign branch (T2.2); keyboard navigation is now owned at the page level.

## Anti-patterns

- Mint `client_op_id` inside `mutationFn`, not at hook creation time.
- Don't read `peaks(id)` derivatively during in-flight ops without `useExposureHasPendingPeakOps` gating.
- Don't assert on Tailwind class strings in tests — use `data-testid` / `data-*` attributes.
- Don't inline appearance utilities (`text-[…]`, `rounded-[…]`, raw colours, side-stripes) in a consumer — `lint:design` fails the build. Put appearance in a `print/ui/` primitive; `className` is placement-only (see "Design system").
