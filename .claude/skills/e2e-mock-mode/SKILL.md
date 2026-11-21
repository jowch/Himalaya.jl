---
name: e2e-mock-mode
description: Run the mocked Playwright E2E suite (npm run e2e) locally or in a cloud session, and read its failures. Use when a mocked spec fails or you need to run a subset.
disable-model-invocation: true
---

# e2e-mock-mode

The mocked suite needs no backend and no server started by hand. `npm run e2e` builds this checkout's bundle into `dist-e2e/` and serves it with `vite preview` on `127.0.0.1:5191`. It never reuses a running server, so a demo or dev server on :5173 can't stand in for your code, and it doesn't touch the demo's port or `dist/`. Every spec takes `test` from `e2e/fixtures.ts`, which answers any `/api` request the spec did not mock with a 501 and fails the test. A backend running on :8080 can't leak into the results, so there's no need to stop it.

## Procedure

1. From `packages/HimalayaUI/frontend/`, run `npx playwright test --trace off`. For a subset, add `--grep "<name>"` or a spec file.
2. In a cloud session the project's Playwright build has no matching browser. Run with a throwaway config that sets `launchOptions.executablePath: "/opt/pw-browsers/chromium"` (see `e2e/AGENTS.md`).
3. Leave the demo and any dev server running. If :5191 is taken, the run fails at startup; free it with `lsof -ti:5191 | xargs kill -9`.
4. Report pass/fail counts and the reason for each failure.

## Reading failures

- `every /api request the page makes must be answered by a mock` lists the unmocked requests. Mock each with data the page would really get, not an empty body that leaves the page in an error state.
- A handler that should let a request through calls `route.fallback()`, which hands it to the catch-all. Never `route.continue()`: that sends the request to the real proxy.
