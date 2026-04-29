---
name: review-pr
description: Review a Himalaya.jl PR and post findings as a GitHub comment in one step. Checks code quality, project conventions, and the full Himalaya gotcha checklist. Usage: /review-pr <number>
---

# review-pr

Reviews a pull request and posts the findings directly as a PR comment without prompting.

## Steps

1. **Get the PR number.** Read it from the skill args. If none was given, run `gh pr list` and ask the user which PR to review.

2. **Fetch PR context:**
   ```bash
   gh pr view <number>
   gh pr diff <number>
   ```

3. **Analyze the diff** against all of the following:

   ### General quality
   - Correctness and logic errors
   - Security issues (injection, auth bypass, secret exposure)
   - Test coverage for new behaviour
   - Performance implications (N+1, unnecessary recomputation)
   - Documentation gaps for non-obvious decisions

   ### Julia backend — SQLite / Oxygen.jl
   - **`DBInterface.lastrowid`** takes the query *result*, not the db. `DBInterface.lastrowid(db)` is wrong.
   - **Materialize rows** via `Tables.rowtable(DBInterface.execute(...))` before accessing fields. Raw cursor values silently vanish after the query closes.
   - **FK columns referencing `users(id)`** that must survive user deletion need `ON DELETE SET NULL` in schema DDL — not at call sites.
   - **Transactional writes in `pipeline.jl`** go inside `_persist_analysis_inner!`. Same pattern for `_reingest_inner!` in `cli.jl`. New multi-write steps anywhere in these paths must stay inside those inner functions.
   - **Stdlib deps must be explicit** in `Project.toml [deps]`. New `using Sockets`/`Printf`/`SparseArrays`/`DelimitedFiles`/`TOML` without a matching entry will break on fresh installs.
   - **Oxygen.jl singleton API:** `@get "/path/{id}" function(req, id::Int) ... end`. JSON body with unqualified `json(req)`, not `Oxygen.json(req)`.
   - **Phase serialization:** `string(nameof(P))` not `string(P)` — the latter returns `"Himalaya.Pn3m"` which breaks SQLite roundtrips.
   - **Detector TIFFs are Q0f31 fixed-point.** `Float32.(channelview(raw))` divides by 2³¹. Use `reinterpret.(Int32, channelview(raw))` to recover photon counts.
   - **Image route `Cache-Control: no-store`** in `routes_exposures.jl` — don't weaken without exposure-id-tied invalidation.
   - **Index scoring:** `score = coverage × consistency`. R² is NOT in the score — it's a UI gate in PhasePanel at 0.98. Guard `cv` against zero mean.
   - **`experiment.toml` is read-only at runtime.** Himalaya never writes inside experiment directories except `himalaya config new`. Any new code path that writes there is a bug.
   - **`config_from_db` / `load_config` share `_build_config`.** Changes to config parsing must go through that shared helper.

   ### Frontend — TypeScript strict / React / Zustand / TanStack Query
   - **`exactOptionalPropertyTypes: true`** — `set({ username: undefined })` fails. Fields must be `string | undefined`, not `username?: string`. Use `authOpts(username)` helper in `queries.ts` for optional auth.
   - **Zustand named actions only.** No direct `useAppState.setState({ ... })` outside `state.ts`. New state transitions need a named action.
   - **State split:** Zustand = client state (active ids, hover, username). TanStack Query = server state. Mutations invalidate scoped keys (`queryKeys.peaks(id)` etc.).
   - **`ImageBitmap.close()` neuters `width`/`height` to 0.** Capture dims before closing: `const { width, height } = bitmap; bitmap.close();`.
   - **`DetectorImage` auto-rotate** is ResizeObserver-driven. JS must set `maxWidth`/`maxHeight` on the canvas element — pure CSS `transform: rotate` doesn't change layout box.
   - **TraceViewer y-fit is floor-only.** `yDomain = [p05·0.7, fullTraceMax·1.2]`. Upper bound uses full-trace max, not windowed — preserves peaks-vs-beam scale.
   - **Observable Plot `.scale().invert()`** is a runtime method not in DOM types — cast with `(el as unknown as { scale: ... })`.
   - **`useCallback` for imperative render deps.** Functions used as `useEffect` deps must be wrapped in `useCallback` with their true deps. Effect then depends on `[theCallback]` only.
   - **`QNumInput` focus-gate pattern.** Any numeric input that can be updated by external events (wheel zoom, etc.) must only sync external value to draft state when not focused.
   - **E2E selectors:** `data-testid`, `role`, `data-*` only. Never assert on Tailwind class strings.
   - **Playwright port:** `127.0.0.1:5173` not `localhost`. Live Julia backend on :8080 leaks past route mocks for query-string URLs.

   ### Build / infra
   - `Manifest.toml` is gitignored. Any change to `Project.toml` deps in a worktree needs a re-resolve.
   - `scripts/Manifest.toml` is also gitignored — `make sysimage` in a fresh worktree needs it copied from main or instantiated.
   - Sysimage path is `build/himalaya.so`. `make check-sysimage` validates Julia version stamp.

4. **Format the review** with these sections (omit any that don't apply to this diff):

   ```
   ## Overview
   <one paragraph — what the PR does>

   ## Issues
   <numbered list of blocking problems — each with file:line and the specific rule violated>
   _(omit section if none)_

   ## Suggestions
   <non-blocking improvements>
   _(omit section if none)_

   ## Strengths
   <what's done well — always at least one point>

   ## Summary
   <one sentence verdict — e.g. "Ready to merge after fixing X.">
   ```

5. **Post the comment:**
   ```bash
   gh pr comment <number> --body "<formatted review>"
   ```

6. Print the returned comment URL.
