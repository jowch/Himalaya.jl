# packages/HimalayaUI/frontend/src — Frontend

## OVERVIEW
React 18 + Vite + TypeScript strict + TailwindCSS 4. TanStack Query for server state, Zustand for client state.

## WHERE TO LOOK
| Task | Location | Notes |
|------|----------|-------|
| App entry | `main.tsx` | Root render, QueryClientProvider, BrowserRouter |
| App shell | `App.tsx` | AppShell, OnboardingFlow, SSE listener, queue persistence |
| Server state | `queries.ts` | TanStack Query hooks; `authOpts(username)` helper |
| API layer | `api.ts` | Fetch wrappers, request/response types |
| Client state | `state.ts` | Zustand store — **use named actions only** |
| Components | `components/` | UI components; `components/ui/` for shared primitives |
| Pages | `pages/` | Route-level page components |
| Hooks | `hooks/` | Custom React hooks |
| Queue system | `lib/queue/` | Mutation queue framework (see its own AGENTS.md) |
| Comparison | `lib/comparison/` | Compare-page domain logic |

## CONVENTIONS
- **Zustand**: named actions only, never `useAppState.setState({ ... })`
- **State split**: Zustand = client state; TanStack Query = server state — don't mix
- **Optimistic peak ids are negative** — UI code must handle `peak.id < 0`
- **Gate skeletons on `query.isLoading`**, not `isPending`
- **`className` on `<Skeleton>` is load-bearing** — pass layout classes to preserve flex chains
- **`exactOptionalPropertyTypes: true`** — use `authOpts(username)` helper, never `{ username: undefined }`
- SSE self-echo filtering uses per-tab `client_id` from `sessionStorage`

## ANTI-PATTERNS
- Mint `client_op_id` inside `mutationFn`, not at hook creation time
- Don't read `peaks(id)` derivatively during in-flight ops without `useExposureHasPendingPeakOps` gating
