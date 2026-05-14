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
| Components | `components/` | See [components/AGENTS.md](components/AGENTS.md) |
| Pages | `pages/` | `ComparePage`, `ComparePageEdit` (Index/Inspect inline into AppShell) |
| Hooks | `hooks/` | `useFocusTrap`, `useGlobalShortcuts`, `useStateFromUrl`, … |
| Library | `lib/` | URL helpers, plot helpers, comparison helpers, figure export |
| Mutation queue | `lib/queue/` | See [lib/queue/AGENTS.md](lib/queue/AGENTS.md) |
| Skeletons | `bones/` | Committed `*.bones.json` + auto-generated `registry.ts` |

## TypeScript strict + `exactOptionalPropertyTypes`

`set({ username: undefined })` fails — optional fields declared as `string | undefined` (rather than `username?: string`) keep this ergonomic. For passing optional values through (e.g. `AuthOpts`), use the `authOpts(username)` helper in `queries.ts` which returns `{}` or `{ username }` — never `{ username: undefined }`.

## Zustand — named actions

Use the store's named actions (`clearUsername`, `setTheme`, `openNavModal`, …). Avoid `useAppState.setState({ ... })` — direct setState bypasses encapsulation and triggers lint warnings. New state transitions go in `state.ts` as named actions.

## State split (load-bearing)

- **Zustand owns *client* state**: active sample/exposure, hoveredIndexId, username.
- **TanStack Query owns *server* state**: experiments, samples, exposures, peaks, indices, groups.
- Mutations invalidate scoped query keys (`queryKeys.peaks(id)`, `queryKeys.groups(id)`). Don't mix the two concerns in the same hook.

## Per-tab SSE identity

SSE self-echo filtering uses a per-tab `client_id` minted into `sessionStorage` on first load (see `lib/clientId.ts`). Audit identity (`actor` / `X-Username`) is unchanged. **Two tabs of the same user are distinct subscribers** — edits in one refresh the other. `client_id` lives for the tab session: survives reload, scoped to one tab. See `docs/event-log.md` §"Client side".

## Imperative renderers in effects

Wrap any function that is both defined inside a component AND used as a `useEffect` dependency in `useCallback` with its true deps. The effect then depends on `[theCallback]` alone — no redundant dep list, no eslint-disable. `TraceViewer`'s overlay renderer follows this pattern.

## Tailwind v4 theming

The dark palette is defined once in `styles.css` via `@theme { --color-* ... }`. Component files use utility classes (`bg-bg`, `text-fg-muted`, `border-accent`). To add a new color, add it to `@theme` first.

## Skeleton loading via boneyard-js

Each load-gated card wraps content in `<Skeleton>` from `boneyard-js/react`. Full reference: `packages/HimalayaUI/docs/boneyard.md`. Two rules that bite hardest:

- **Gate on `query.isLoading`, not `isPending`.** `isLoading = isPending && isFetching` — disabled queries and background refetches stay skeleton-free; only true cold fetches animate. Wrong gating ⇒ flicker on every refetch.
- **`className` on `<Skeleton>` is load-bearing.** Boneyard adds two wrapper divs that break parent flex chains (e.g. ChatCard's message list collapsing to 60px). Pass `flex-1 min-h-0 flex flex-col` (or `h-full w-full`) to inherit the original child's layout role. Companion CSS in `styles.css`: `[data-boneyard-content] { display: contents }`.

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
