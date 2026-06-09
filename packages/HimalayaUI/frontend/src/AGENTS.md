# packages/HimalayaUI/frontend/src — Frontend

React 18 + Vite + TypeScript strict + TailwindCSS 4. TanStack Query for server state, Zustand for client state.

## Where to look

| Task | Location | Notes |
|------|----------|-------|
| App entry | `main.tsx` | `StrictMode > ErrorBoundary > QueryClientProvider > App` |
| App shell | `App.tsx` | Composition root; Zustand selectors + TanStack Query |
| Server state | `queries.ts` | TanStack Query hooks; `authOpts(username)` helper |
| API layer | `api.ts` | Fetch wrappers; AuthOpts per-call for mutations |
| Client state | `state.ts` | Zustand store — **use named actions only** |
| Phase palette | `phases.ts` | phase → color mapping |
| Tailwind theme | `styles.css` | `@theme { --color-* … }` |
| App shell + routing | `print/shell/` | `CorpusShell`, `CorpusTopbar`, `AppRoutes`, `IndexSlugRedirect`, `StaleUrlPage`, `ResolvingFallback`, `OnboardingFlow`, `NavModal`, `InfrastructureBanner`. See [print/shell/AGENTS.md](print/shell/AGENTS.md) |
| Composites | `print/components/` | Page-composing components (rails, plates, panels, rows, modals) built from the `ui/` primitives |
| UI primitives | `print/ui/` | Closed-look design-system primitives (Button, Card, SegmentedControl, PhaseChip, PhaseStrip, ModalShell, Kicker, IconButton, ScoreBar, Dot, ToastContainer, HintText, …). Appearance lives here; consumer `className` is placement-only. See "Design system" below. |
| Render layers | `print/{plot,detector,comb,waterfall,export}/` | Appearance-authoring render layers (trace-plot engine, detector image, comb/residual, waterfall, figure-export marks) — excluded from the `lint:design` appearance guard |
| Pages | `print/pages/` | `SamplesPage`, `LoupePage`, `FocusPage`, `SeriesFolioPage`, `SeriesScopingPage`, `SeriesBuilderPage` (all under the single `CorpusShell`; legacy Index/Inspect/Compare pages + `AppShell` retired) |
| Hooks | `hooks/` | `useFocusTrap`, `useGlobalShortcuts`, `useStateFromUrl`, … |
| Library | `lib/` | URL helpers, plot helpers, comparison helpers, figure export |
| Mutation queue | `lib/queue/` | See [lib/queue/AGENTS.md](lib/queue/AGENTS.md) |
| Skeletons | `bones/` | Committed `*.bones.json` + auto-generated `registry.ts` |

## TypeScript strict + `exactOptionalPropertyTypes`

`set({ username: undefined })` fails — optional fields declared as `string | undefined` (rather than `username?: string`) keep this ergonomic. For passing optional values through (e.g. `AuthOpts`), use the `authOpts(username)` helper in `queries.ts` which returns `{}` or `{ username }` — never `{ username: undefined }`.

## Zustand — named actions

Use the store's named actions (`clearUsername`, `setTutorialSeen`, `openNavModal`, …). Avoid `useAppState.setState({ ... })` — direct setState bypasses encapsulation and triggers lint warnings. New state transitions go in `state.ts` as named actions.

## State split (load-bearing)

- **Zustand owns *client* state**: active sample/exposure, hoveredIndexId, username.
- **TanStack Query owns *server* state**: experiments, samples, exposures, peaks, indices, assignment.
- Mutations invalidate scoped query keys (`queryKeys.peaks(id)`, `queryKeys.indices(id)`). Don't mix the two concerns in the same hook.

## Per-tab SSE identity

SSE self-echo filtering uses a per-tab `client_id` minted into `sessionStorage` on first load (see `lib/clientId.ts`). Audit identity (`actor` / `X-Username`) is unchanged. **Two tabs of the same user are distinct subscribers** — edits in one refresh the other. `client_id` lives for the tab session: survives reload, scoped to one tab. See `docs/event-log.md` §"Client side".

## Imperative renderers in effects

Wrap any function that is both defined inside a component AND used as a `useEffect` dependency in `useCallback` with its true deps. The effect then depends on `[theCallback]` alone — no redundant dep list, no eslint-disable. `TraceViewer`'s overlay renderer follows this pattern.

## Tailwind v4 theming

"The Print" palette (light warm paper, terracotta accent hue 38, Newsreader serif) is defined once in `styles.css` via `@theme { --color-* ... }` — a single identity, no theme toggle (R0a, #221). The legacy dark-era neutral-ramp shim (the duplicated `bg-`/`text-`/`border-` neutral utilities that mirrored the canonical Print names) was excised in R3-F (#259): use the canonical Print utilities directly — `bg-paper`/`bg-paper-sunk`/`bg-plate`, `text-ink`/`text-ink-soft`/`text-ink-faint`, `border-hair`/`border-hair-strong`, `text-print-accent` (or `accent` for terracotta). Reintroducing an old name is self-revealing: Tailwind won't generate the utility and its `--color-*` custom property no longer resolves. Serif titles use the `text-display`/`text-headline` roles. To add a new color, add it to `@theme` first.

## Design system — closed-look primitives + the design guard (ENFORCED)

The Print's recurring patterns live as **closed-look** primitives in `components/ui/` (Button, Card, SegmentedControl, PhaseChip, PhaseStrip, ModalShell, Kicker, IconButton, ScoreBar, Dot, ToastContainer, HintText). They own their appearance via semantic props (`variant` / `size` / `tone` / domain props); a consumer's `className` is **placement-only** (margins, position, grid). To change how a primitive looks, build a variant *into the primitive* — the idiom is a `Record<Variant,string>` map + a tiny local `cx()` join helper (no cva/clsx/tailwind-merge). Don't restyle from the outside.

This is **mechanically enforced** (2026-05-29 extraction). `scripts/check-design.mjs` runs as a pure-absolute `lint:design` step prepended to `npm run build` (plus a warn-only PostToolUse hook), and **fails the build** on any inline appearance utility *outside* `components/ui/**`:

- arbitrary type size `text-[…]` → use a named scale role (`text-xs/sm/base/lg/xl/headline-lg/display`)
- arbitrary radius `rounded-[…]` → `rounded-sm` / `rounded-md` (both 5px) / `rounded-full`
- raw colour literal (`oklch(` / `rgba(` / quoted `#hex`) → a `--color-*` token utility
- side-stripe `border-l/r` > 1px → a full border + a leading icon/word instead

Only the colour-AUTHORING files are exempt (rules #3/#5 share an allowlist: `phases.ts`, `lib/comparison/coloring.ts`, `lib/figure-export/**`, `MemberHeatmapLayer.tsx`, `DetectorImage.tsx`, `FocusDetectorPanel.tsx`, `main.tsx`). Need a colour anywhere else → add a `--color-*` token to `@theme`, then use the utility. Visual reference: `docs/design-system.html`; full system: root `DESIGN.md`.

## Skeleton loading via boneyard-js

Each load-gated card wraps content in `<Skeleton>` from `boneyard-js/react`. Full reference: `packages/HimalayaUI/docs/boneyard.md`. Two rules that bite hardest:

- **Gate on `query.isLoading`, not `isPending`.** `isLoading = isPending && isFetching` — disabled queries and background refetches stay skeleton-free; only true cold fetches animate. Wrong gating ⇒ flicker on every refetch.
- **`className` on `<Skeleton>` is load-bearing.** Boneyard adds two wrapper divs that break parent flex chains (e.g. a scrolling list inside a card collapsing to a fixed small height). Pass `flex-1 min-h-0 flex flex-col` (or `h-full w-full`) to inherit the original child's layout role. Companion CSS in `styles.css`: `[data-boneyard-content] { display: contents }`.

## Mutation queue — load-bearing one-liners

Full architecture in `docs/mutation-queue.md`; queue internals in `lib/queue/AGENTS.md`. The invariants that bite UI code outside the queue:

- **Optimistic placeholder ids are NEGATIVE.** `Peak.id < 0` means "not yet confirmed by server"; SSE confirmation overwrites with the positive server id. UI code that filters or compares peak ids must handle negatives.
- **`useExposureHasPendingPeakOps` gates any UI that reads `peaks(id)` derivatively** while a peak op is in flight (StaleIndicesBanner, useSpeculativeSnap). Without it: flicker as optimistic / HTTP / SSE land out of order.

## Multi-layer contract testing

Every reconciliation contract has six layers (route emit → SSE payload → `applyRemoteToCache` merge → cache row → `onMutate` → `onSuccess`). When fixing a bug at one layer, add a regression row at every other layer where the same class can manifest. See `docs/contract-testing.md` for canonical paired test files (`cache-shape.test.ts`, `sseEventPayload.contract.test.ts`, `rollbackSymmetry.test.ts`, `authHeaders.test.ts`, `test_route_response_shapes.jl`, `test_idempotency_replay_invariant.jl`).

## Anti-patterns

- Mint `client_op_id` inside `mutationFn`, not at hook creation time.
- Don't read `peaks(id)` derivatively during in-flight ops without `useExposureHasPendingPeakOps` gating.
- Don't assert on Tailwind class strings in tests — use `data-testid` / `data-*` attributes.
- Don't inline appearance utilities (`text-[…]`, `rounded-[…]`, raw colours, side-stripes) in a consumer — `lint:design` fails the build. Put appearance in a `components/ui/` primitive; `className` is placement-only (see "Design system").
