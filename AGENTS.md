# AGENTS.md — Himalaya.jl quick reference

**Generated:** 2026-05-07
**Commit:** a907bea
**Branch:** main

## OVERVIEW

Julia monorepo for indexing SAXS diffraction patterns. Core `Himalaya` package (root) finds Bragg peaks and identifies liquid-crystalline phases. `HimalayaUI` (under `packages/`) is a full-stack web app — Julia/Oxygen.jl REST backend + React/Vite frontend — for running and curating analyses.

## STRUCTURE

```
.
├── src/                         # Core Himalaya package (peak finding, phase indexing)
├── test/                        # Core package tests
├── packages/
│   └── HimalayaUI/              # Full-stack web app
│       ├── src/                 # Julia backend (Oxygen.jl, SQLite, SSE)
│       ├── test/                # Backend tests
│       ├── frontend/            # React 18 + Vite frontend
│       │   ├── src/             # Frontend source
│       │   ├── test/            # Vitest unit tests
│       │   └── e2e/             # Playwright E2E specs
├── docs/                        # Architecture docs (peak-finding, queue, events, scoring)
├── scripts/                     # Build scripts (PackageCompiler sysimage)
└── bin/                         # CLI wrapper script
```

## First session

Read these in order before touching code:
0. `docs/decisions.md` — Jonathan's product decisions, one row each, with the code or test that implements it or "planned". Authors and reviewers check a change against it first; a change that contradicts a row updates the row.
1. `docs/peak-finding.md` — peak-finding design (load-bearing)
2. `docs/experiment-config.md` — if touching the New-experiment funnel, experiment create/scan, or config.jl
3. `docs/scoring.md` — if touching score/auto_group/remove_subsets
4. `docs/event-log.md` — if touching events.jl, hash.jl, SSE, or the stale-indices reanalyze affordance
5. `docs/mutation-queue.md` — if touching lib/queue/, idempotency.jl, or applyRemoteToCache.ts
6. `docs/contract-testing.md` — if fixing queue/SSE/cache reconciliation bugs

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

# Frontend build and type-check of tests and e2e specs (both must pass before PR)
npm run build
npm run typecheck:test
```

## Conventions

- **TDD**: failing test → minimal impl → verify → commit (each step is a commit)
- **One responsibility per file** — split when a file accumulates concepts
- **Regression floors**, not hard-coded counts: `recall ≥ floor`, `spurious ≤ ceiling`
- **Worktrees** for feature branches: `git worktree add .claude/worktrees/<topic> -b <topic>`

## Pre-submit verification checklist (author)

Before declaring any bugfix or behavioral change complete, do these in order:

1. **Revert the fix, confirm the test fails.**
   - Temporarily revert your code change.
   - Run the new regression test.
   - If it still passes, the test is a false positive — rewrite the assertion.
   - Re-apply the fix before continuing.

2. **Read wrapper contracts before assuming behavior.**
   - When using a wrapper (e.g., `useQueueMutation`), read its return type.
   - `mutate` = fire-and-forget; errors surface via `.error`, not by throwing.
   - `try/catch` around `mutate()` only catches synchronous exceptions, not HTTP failures.

3. **Assert on observable side-effects, not internal instances.**
   - Bad: `expect(mock.instances.length).toBe(1)` — too broad, counts anything.
   - Good: `expect(element.style.height).not.toBe("0px")` — asserts the state change.
   - Ask: "What user-visible behavior proves this works?"

4. **Use the reviewer agent.**
   - Spawn `.claude/agents/frontend-reviewer.md` (or the relevant reviewer) before sending.
   - The reviewer catches patterns you are blind to after writing the code.
   - Do not self-review complex PRs — delegation is cheaper than review cycles.

## PR Review Workflow

When reviewing a PR (as reviewer, not author):

1. **Check out the PR branch into a worktree**:
   ```bash
   git worktree add .claude/worktrees/pr-<N> -b pr-<N> origin/pull/<N>/head
   ```
2. **Read changed files from the worktree** — full context, not just diff.
3. **Run tests from the worktree**:
   ```bash
   cd .claude/worktrees/pr-<N>
   julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
   ```
4. **Review and post comment** — `gh pr comment <N> --body "..."`. Use severity levels (🔴 blocker, 🟡 medium, 🟢 minor).
5. **If author rebuts**, verify fixes point-by-point.
6. **Clean up after merge**:
   ```bash
   git worktree remove .claude/worktrees/pr-<N>
   ```

## Critical Julia/SQLite patterns

- `DBInterface.lastrowid(result)` — takes the **result**, not the db
- `Tables.rowtable(DBInterface.execute(...))` — materialize rows before query closes
- SQL NULL → `missing` (not `nothing`): use `ismissing(row.field)` and normalize
- FK enforcement is ON — use `ON DELETE SET NULL` in DDL for FKs that must survive deletion
- PKs use `AUTOINCREMENT` on entity tables (experiments, samples, exposures, peaks, indices) so a deleted id is never reused
- `apply_event!(InTransaction(), ...)` inside `with_idempotency` — never the public `apply_event!(db, req; ...)` from within a queue mutation
- `persist_analysis!` is transactional, and `scan_and_group!` persists loads/samples/exposures in one write transaction (its geometry fill is a separate, earlier one) — add new writes inside `_persist_analysis_inner!` / that persist transaction

## Critical frontend patterns

- **Zustand**: use named actions, not `useAppState.setState({ ... })`
- **State split**: Zustand = client state; TanStack Query = server state — don't mix
- **Optimistic peak ids are negative** — UI code must handle `peak.id < 0`
- **Mint `client_op_id` inside `mutationFn`**, not at hook creation time
- **`exactOptionalPropertyTypes: true`** — use `authOpts(username)` helper, never `{ username: undefined }`
- SSE self-echo filtering uses per-tab `client_id` from `sessionStorage`

## Load-bearing one-liners

- `analyze_exposure!` synthesizes effective peaks as `auto_peaks − excludes ∪ adds` — touch only with curation-lifecycle tests green
- `score(index)` = `coverage × consistency` — R² is NOT part of score and gates nothing
- `MutationCache.getAll()` insertion order is load-bearing for replay-as-rerun
- Experiment dirs are read-only at runtime
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
| Ingest funnel + scan | `packages/HimalayaUI/src/routes_fs.jl`, `packages/HimalayaUI/src/routes_experiments.jl`, `packages/HimalayaUI/src/ingest.jl`, `packages/HimalayaUI/src/grouping.jl` |
| Analysis pipeline | `packages/HimalayaUI/src/pipeline.jl` |
| CLI | `packages/HimalayaUI/src/cli.jl` |
| Events + SSE | `packages/HimalayaUI/src/events.jl` |
| Idempotency | `packages/HimalayaUI/src/idempotency.jl` |
| REST routes | `packages/HimalayaUI/src/routes_*.jl` |
| Frontend state | `packages/HimalayaUI/frontend/src/state.ts`, `queries.ts` |
| Mutation queue | `packages/HimalayaUI/frontend/src/lib/queue/` |

See `CLAUDE.md` for more info and gotchas.
