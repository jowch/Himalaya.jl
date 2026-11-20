# frontend/e2e — Playwright E2E

Two suites:

- **`e2e/`** — default mocked suite. `page.route` intercepts `/api/...`. Run via `npm run e2e`; CI runs it in the `Frontend E2E (mocked)` job. It builds the bundle and serves it with `vite preview` on 127.0.0.1:5173 (or reuses a server already there, outside CI).
- **`e2e/live/`** — live-integration. Real backend + dev DB. Run via `npm run e2e:live`. Operator brings up backend (port 8090) + Vite (5180) manually. Full runbook: [live/README.md](live/README.md).

## Every request is mocked

Specs import `test` and `expect` from `./fixtures`, not `@playwright/test` (`test/e2e-fixtures.test.ts` enforces it). The fixture answers any `/api/` request the spec did not mock with a 501 and fails the test, listing the requests. Without it, unmocked requests went through the Vite proxy to whatever backend ran on :8080, or to nothing, and specs passed on error states. When a spec trips it, mock the request with data the page would really get.

In a cloud session the project's Playwright build differs from the preinstalled Chromium; point a throwaway config at it:

```ts
import base from "./playwright.config";
export default { ...base, use: { ...base.use, launchOptions: { executablePath: "/opt/pw-browsers/chromium" } } };
```

## Selectors

Use stable `data-*` attributes, `data-testid`, or `role`. **Never assert on Tailwind class strings** — they change when styling evolves. Stable testid families: `data-index-id`, `data-active`, `data-sample-id` (cold-assign panel rows only).

Single test by name:

```bash
node_modules/.bin/playwright test --grep "batch-reject"
```

## Seeding Zustand state

There is no `data-exposure-id`, and `data-sample-id` only exists on the cold-assign panel rows (`ColdAssignPanel`), not on the main sample navigation surfaces. To navigate to a specific sample/exposure, seed Zustand state in `localStorage` before `page.goto("/")`:

```ts
await page.addInitScript((state) => {
  localStorage.setItem("himalaya-ui:state", JSON.stringify({ state, version: 5 }));
}, { activeSampleId: 1, activeExposureId: 2 });
await page.goto("/");
```

The persisted version is currently `5` (see `src/state.ts`); any older version number also works because the store's `migrate` function upgrades old seeds, so the exact number is not load-bearing.

See the live specs for the full pattern.

## Port binding

`npm run e2e` builds the bundle and serves it on `http://127.0.0.1:5173` (not `localhost`) with `vite preview`. There's no need to start a server first. Outside CI it reuses whatever already answers on :5173 (a dev server or the demo). Because every `/api` request is answered inside the browser, the backend behind that server is never reached. If something else holds :5173 and doesn't serve the app, free it: `lsof -ti:5173 | xargs kill -9`.

## Live-integration timing

`page.goto("/")` in live mode races SSE subscription. **Wait ~800ms before any mutation** that expects an SSE echo, otherwise the test browser misses the broadcast. Full guidance in [live/README.md](live/README.md).

## Anti-patterns

- Don't use `localhost` in any URL — Playwright config and CORS expectations all assume `127.0.0.1`.
- Don't `page.route` and `npm run e2e:live` together — live tests must hit the real backend.
- Don't assume CSS-class changes are E2E-stable — only `data-*` and `role` are.
