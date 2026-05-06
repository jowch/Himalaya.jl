# AGENTS.md — Himalaya.jl quick reference

## First session

Read these in order before touching code:
1. `docs/peak-finding.md` — peak-finding design (load-bearing)
2. `docs/experiment-config.md` — if touching config/manifest/cli
3. `docs/scoring.md` — if touching score/auto_group/remove_subsets
4. `docs/event-log.md` — if touching events.jl, hash.jl, SSE, or StaleIndicesBanner
5. `docs/mutation-queue.md` — if touching lib/queue/, idempotency.jl, or applyRemoteToCache.ts
6. `docs/contract-testing.md` — if fixing queue/SSE/cache reconciliation bugs

## Project structure

- **Root** = core `Himalaya` package (peak finding, phase indexing, scoring)
- **`packages/HimalayaUI/`** = full-stack web app
  - `src/` = Julia/Oxygen.jl REST backend + CLI
  - `frontend/` = React 18 + Vite + TypeScript strict

## Running tests

```bash
# Core Himalaya
julia --project=. -e 'using Pkg; Pkg.test()'

# HimalayaUI backend
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'

# Frontend unit tests — from packages/HimalayaUI/frontend/
npm test

# Frontend E2E (mocked)
npm run e2e

# Frontend build (must pass before PR)
npm run build
```

## Conventions

- **TDD**: failing test → minimal impl → verify → commit (each step is a commit)
- **One responsibility per file** — split when a file accumulates concepts
- **Regression floors**, not hard-coded counts: `recall ≥ floor`, `spurious ≤ ceiling`
- **Worktrees** for feature branches: `git worktree add .claude/worktrees/<topic> -b <topic>`

## PR Review Workflow

When reviewing a PR (as reviewer, not author):

1. **Check out the PR branch into a worktree** — never review from diff alone:
   ```bash
   git worktree add .claude/worktrees/pr-<N> -b pr-<N> origin/pull/<N>/head
   ```
   This gives you the exact commit that will be merged, isolates from your current work, and enables LSP navigation + test runs against the PR code.

2. **Read the changed files from the worktree** — full file context, not just the delta. Pay attention to docstrings, type signatures, and module-level changes (new `using` imports, etc.).

3. **Run the tests from the worktree**:
   ```bash
   cd .claude/worktrees/pr-<N>
   julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
   ```
   A diff can't tell you whether a tab-separated fixture still parses correctly after switching to CSV.jl.

4. **Review and post the comment** — use `gh pr comment <N> --body "..."`. Structure the review with severity levels (🔴 blocker, 🟡 medium, 🟢 minor) and specific line-level issues.

5. **If the author rebuts**, re-read the updated diff and verify each fix point-by-point. Post a rebuttal review confirming whether each issue is adequate, inadequate, or over-fixed.

6. **Clean up after merge**:
   ```bash
   git worktree remove .claude/worktrees/pr-<N>
   ```

## Critical Julia/SQLite patterns

- `DBInterface.lastrowid(result)` — takes the **result**, not the db
- `Tables.rowtable(DBInterface.execute(...))` — materialize rows before query closes
- SQL NULL → `missing` (not `nothing`): use `ismissing(row.field)` and normalize
- FK enforcement is ON — use `ON DELETE SET NULL` in DDL for FKs that must survive deletion
- PKs use `AUTOINCREMENT` on mention-target tables (experiments, samples, exposures, peaks, indices)
- `apply_event!(InTransaction(), ...)` inside `with_idempotency` — never the public `apply_event!(db, req; ...)` from within a queue mutation
- `persist_analysis!` and `reingest!` are transactional — add new writes inside `_persist_analysis_inner!` / `_reingest_inner!`

## Critical frontend patterns

- **Zustand**: use named actions, not `useAppState.setState({ ... })`
- **State split**: Zustand = client state; TanStack Query = server state — don't mix
- **Optimistic peak ids are negative** — UI code must handle `peak.id < 0`
- **Mint `client_op_id` inside `mutationFn`**, not at hook creation time
- **Gate skeletons on `query.isLoading`**, not `isPending`
- **`className` on `<Skeleton>` is load-bearing** — pass layout classes to preserve flex chains
- **`exactOptionalPropertyTypes: true`** — use `authOpts(username)` helper, never `{ username: undefined }`
- SSE self-echo filtering uses per-tab `client_id` from `sessionStorage`

## Load-bearing one-liners

- `analyze_exposure!` synthesizes effective peaks as `auto_peaks − excludes ∪ adds` — touch only with curation-lifecycle tests green
- `score(index)` = `coverage × consistency` — R² is NOT part of score (UI hard gate at 0.98)
- `StaleIndicesBanner` is gated on hash mismatch + `useExposureHasPendingPeakOps`
- `useExposureHasPendingPeakOps` gates any UI reading `peaks(id)` derivatively during in-flight ops
- `MutationCache.getAll()` insertion order is load-bearing for replay-as-rerun
- Experiment dirs are read-only at runtime (except `himalaya config new`)
- `parse_manifest(cfg::ExperimentConfig, source)` is the config-driven method — use in new code
- `Fm3m` is missing from `indexpeaks` dispatch — known gap, don't fix opportunistically

## Contract testing (six layers)

When fixing a queue/SSE/cache bug, add regression tests at **all** relevant layers:
1. Route emit
2. SSE payload
3. `applyRemoteToCache` merge
4. Cache row
5. `onMutate`
6. `onSuccess`

See `docs/contract-testing.md` for canonical paired test files.

## Common file locations

| What | Where |
|---|---|
| Peak finding | `src/peakfinding.jl`, `src/persistence.jl`, `src/sharpness.jl`, `src/threshold.jl` |
| Phase types & indexing | `src/phase.jl`, `src/index.jl` |
| DB schema + CRUD | `packages/HimalayaUI/src/db.jl` |
| Config/manifest | `packages/HimalayaUI/src/config.jl`, `packages/HimalayaUI/src/manifest.jl` |
| Analysis pipeline | `packages/HimalayaUI/src/pipeline.jl` |
| CLI | `packages/HimalayaUI/src/cli.jl` |
| Events + SSE | `packages/HimalayaUI/src/events.jl` |
| Idempotency | `packages/HimalayaUI/src/idempotency.jl` |
| REST routes | `packages/HimalayaUI/src/routes_*.jl` |
| Frontend state | `packages/HimalayaUI/frontend/src/state.ts`, `queries.ts` |
| Mutation queue | `packages/HimalayaUI/frontend/src/lib/queue/` |
| Skeleton loading | `packages/HimalayaUI/docs/boneyard.md` |
