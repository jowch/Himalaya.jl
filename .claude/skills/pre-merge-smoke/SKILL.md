---
name: pre-merge-smoke
description: Run the manual smoke checklist for queue-touching changes before merging. Spins up dev server, exercises the queue end-to-end through the browser preview, reports findings. The smoke steps in §3 are the canonical manual checklist.
disable-model-invocation: true
---

# pre-merge-smoke

Walks the manual smoke checklist for queue-touching changes. The mutation queue framework's correctness depends on user-observable behaviors that unit + integration tests can't fully cover (cache flicker, banner timing, two-tab sync, sessionStorage rehydrate). This skill exercises each scenario in the actual browser and reports findings.

## When to use

Before merging any PR that touches:
- `lib/queue/*` (framework changes)
- `lib/queue/mutators/*` (mutator additions or refactors)
- `routes_peaks.jl`, `routes_analysis.jl`, or any route emitting `apply_event!`
- `applyRemoteToCache.ts` (SSE-driven cache merge)
- `StaleIndicesBanner.tsx`, `SpeculativeBuilder.tsx`, `Toast.tsx`, `InfrastructureBanner.tsx`
- `App.tsx` SSE wiring or `attachPersistence`/`rehydrate` calls

For non-queue PRs, a clean unit + integration test run is enough; this skill is overkill.

## Procedure

### 1. Scope

The browser-driving smoke steps in §3 below are the canonical manual checklist (the queue's design plan was consolidated away — `docs/event-log.md` §3a and `docs/mutation-queue.md` are the living architecture references). Automated coverage (`Pkg.test("HimalayaUI")`, `npm test`, `npm run build`, the two-context Playwright multiplayer replay test) runs separately; this skill is only the manual browser sweep that those can't cover.

### 2. Build + serve

```bash
cd packages/HimalayaUI/frontend
npm run build
```

Then in the worktree root:

```bash
# Use a dev experiment dir; defaults to ~/.himalaya/himalaya.db
bin/himalaya serve /path/to/dev/experiment --port 8080 &
```

(If the user hasn't set up a dev experiment, ask before scaffolding one — `himalaya config new --dir` writes to disk.)

Open the preview at `http://127.0.0.1:8080` (backend serves frontend from `frontend/dist/`).

### 3. Walk each smoke step

Use the `preview_*` browser tools to exercise each scenario. Report pass/fail per step with a concise note.

**Step 1 — Peak exclude is instant + no flicker**
- Navigate to an exposure with auto peaks.
- Alt-click a peak (or whatever exclude UX is current — verify with `preview_snapshot`).
- Assert: peak's `data-excluded` flips within one frame; no banner appears; indices update without observable round-trip.

**Step 2 — Two-tab cross-tab sync**
- Open the same exposure in a second tab (use `tabs_create_mcp`).
- In tab A, exclude a different peak.
- In tab B (within ~1s), assert the same peak shows excluded without manual refresh. Use `preview_snapshot` or poll `data-excluded`.

**Step 3 — Validation toast on network failure**
- Use `preview_eval` to inject `window._origFetch = fetch; window.fetch = () => Promise.reject(new TypeError('offline'))` (or use DevTools network throttling if available).
- Trigger a curation. Assert: validation toast appears OR infrastructure banner appears (depending on classification — 4xx vs network).
- Restore fetch: `preview_eval` `window.fetch = window._origFetch`.

**Step 4 — sessionStorage rehydrate on reload**
- Trigger a curation that takes long enough to inspect (use a slow network throttle).
- Reload the tab mid-flight via `preview_eval` `location.reload()`.
- Assert: after reload, the queue rehydrates and the op replays via HTTP. Cache settles to the correct state. No duplicate row in `peak_curations` (verify via SQL if SQLite MCP is available, or via the peaks list count).

**Step 5 — Speculative builder during pending peak op**
- Open the speculative builder modal.
- Trigger a peak op while the modal is open (e.g. via the trace).
- Assert: modal shows "updating to latest…" subtext; snap suggestions are gated; after the peak op settles, snap suggestions populate.

### 4. Telemetry checks (optional, query-based)

If a SQLite MCP server is configured against `~/.himalaya/himalaya.db`:

```sql
-- Fast-skip latency sanity (post-recent-curations)
SELECT json_extract(payload, '$.duration_ms') AS dur,
       json_extract(payload, '$.findpeaks_skipped') AS fp_skip,
       json_extract(payload, '$.indexpeaks_skipped') AS ip_skip
FROM user_actions WHERE action = 'analyze_run'
ORDER BY id DESC LIMIT 20;

-- post_state size distribution
SELECT json_extract(payload, '$.post_state_size_bytes') AS bytes
FROM user_actions WHERE action = 'analyze_run'
  AND payload IS NOT NULL ORDER BY id DESC LIMIT 100;
-- Expect P50 < 3KB, P99 < 8KB.

-- Idempotent retries
SELECT COUNT(*) FROM idempotent_responses;
SELECT COUNT(DISTINCT client_op_id) FROM user_actions WHERE client_op_id IS NOT NULL;
```

Without SQLite MCP, skip this section or query via a one-off `julia` REPL.

### 5. Report

Report format:

```
## pre-merge-smoke findings

**Scope:** <what was tested>

### Pass
- Step 1: <one-line result>
- ...

### Fail / regression
- Step N: <what went wrong, with screenshot path if visual>

### Skipped
- Step N: <why — e.g., no SQLite MCP for telemetry>
```

If all steps pass: "All N smoke steps pass; OK to merge."

### 6. Cleanup

- Stop the backend: `kill %1` (or whichever job number).
- Close any browser tabs the skill opened.
- Reset any `window.fetch` overrides via `preview_eval`.

## Fallback triggers

If a smoke step exposes a user-visible regression the test suite missed, don't merge to ship faster — escalate to a design conversation. The queue's three pre-committed fallback triggers (the conditions under which its optimistic-merge approach should be backed out) are:

1. **A mutator needs to inspect *other* pending mutations to compute its effect.** That breaks the framework's `(op, baseState) → (cacheEffect, request)` shape — the abstraction has leaked. Fall back to per-mutation `onMutate` optimistic updates without queue infrastructure; lose replay-on-remote-event, keep the autoReanalyze elimination.
2. **The `analyze_exposure!` fast-skip no longer reduces the no-change path to microseconds** (steady state ~150µs; ceiling ~500µs). The synchronous-reanalyze pattern then relocates latency rather than eliminating it — revert to client-side autoReanalyze and treat the queue as cache-merging-only.
3. **The response-body cache produces user-visible drift** between rehydrated retries and current server state more than a few times per long session (canonical case: a deploy that changes a response shape mid-session). Schema-version the cached responses or scope retries by deploy-version.

## What this skill is NOT for

- Unit/integration test runs — use the regular `npm test` / `Pkg.test` flow.
- Ad-hoc UI exploration — just open the browser; this skill is for the deliberate pre-merge sweep.
- Performance benchmarking beyond the fast-skip sanity check — use a separate benchmark harness if needed.
