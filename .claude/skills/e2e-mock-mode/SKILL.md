---
name: e2e-mock-mode
description: Run the Playwright E2E suite cleanly without the live Julia backend on :8080 leaking past route mocks. Detects port 8080, prompts the user before stopping, runs E2E with the documented Vite host binding, restores port state. Use when smoke tests appear flaky/failing — usually it's the live backend leaking through.
disable-model-invocation: true
---

# e2e-mock-mode

Playwright route mocks like `**/api/samples/10/exposures` don't match URLs with query strings (e.g. `?exclude_rejected=true`), so the request leaks through Vite's `/api` proxy to the live Julia backend on `:8080`. When the backend is up, mocked tests can silently consume real data and fail in confusing ways. This skill makes that interaction explicit.

## Procedure

1. **Detect the conflict.**
   ```bash
   julia_pid=$(lsof -ti:8080 2>/dev/null | head -1)
   vite_pid=$(lsof -ti:5173 2>/dev/null | head -1)
   ```
   If `julia_pid` is set, the backend is running and will leak past route mocks.

2. **Surface the problem and ask before stopping.** Do NOT stop the backend without user confirmation — they may have it up intentionally for live preview / curl / Inspect-page work.
   ```
   Live Julia backend running on :8080 (PID <julia_pid>). Playwright route mocks don't match URLs with query strings, so this will leak past them and 3 smoke tests will fail with confusing data-active assertion errors.

   Stop the Julia backend for this E2E run? [y/N]
   ```

3. **If user agrees, stop the backend.**
   ```bash
   kill "$julia_pid" 2>/dev/null
   sleep 1
   lsof -ti:8080 2>/dev/null && echo "still running" || echo "stopped"
   ```
   Remember the PID so you can offer to restart it after.

4. **Make sure Vite is bound to 127.0.0.1, not localhost.** `playwright.config.ts` expects `127.0.0.1:5173`. If Vite is already running on `:5173`, leave it (Playwright reuses). If not, start it bound correctly:
   ```bash
   nohup npm run dev -- --host 127.0.0.1 > /tmp/vite-e2e.log 2>&1 &
   ```
   from `packages/HimalayaUI/frontend/`.

5. **Run the E2E suite.**
   ```bash
   cd packages/HimalayaUI/frontend
   node_modules/.bin/playwright test
   ```
   Or a subset if the user named one: `--grep "<name>"` or a specific file.

6. **Report results.** Pass/fail counts, plus any specific failure reasons.

7. **Offer to restart the backend** if it was stopped:
   ```
   E2E run complete (X passed, Y failed). Restart Julia backend on :8080? [y/N]
   ```
   If yes, run the documented serve command from CLAUDE.md:
   ```bash
   julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
     serve <experiment-dir> --port 8080
   ```
   The user supplies the experiment dir (it's not safe to guess).

## Args

- `--grep "<name>"` — pass through to Playwright to run a single test
- `<file>` — pass through to run a single spec file (e.g. `e2e/inspect.spec.ts`)
- `--keep-backend` — skip steps 2/3, run anyway (user takes responsibility for the leak)

## Why user-only (`disable-model-invocation: true`)

Stopping the Julia backend can disrupt other work the user is doing — manual preview testing, curl-based debugging, an open browser session. Always confirm before killing it.

## Gotchas

- **Port reuse:** Playwright's `webServer` config has `reuseExistingServer: !process.env.CI`. If Vite is already running, Playwright reuses it. If your Vite was started with `npm run dev` (default `--host`), it binds to `localhost` (potentially `::1` IPv6) and Playwright's `127.0.0.1` checks fail with a 60-second timeout. Always use `--host 127.0.0.1` when starting Vite for E2E.
- **Test count:** 14 tests total (5 inspect.spec.ts + 9 smoke.spec.ts). All should pass when the backend is stopped. If smoke tests still fail with the backend down, that's a real bug, not the leak.
- **The leak is asymmetric:** read-only routes leak through and return real data. Mutations (POST/PATCH) might 404 or silently succeed against the wrong sample. Don't trust E2E mutation assertions when the backend is up.
