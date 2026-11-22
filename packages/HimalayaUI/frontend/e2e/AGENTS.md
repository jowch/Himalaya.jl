# frontend/e2e — Playwright E2E

Two suites:

- **`e2e/`** — default mocked suite. `page.route` intercepts `/api/...`. Run via `npm run e2e`; CI runs it in the `Frontend E2E (mocked)` job. It builds this checkout's bundle into `dist-e2e/` and serves it with `vite preview` on 127.0.0.1:5191. It has its own port and build folder, and never reuses a running server.
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

The workspace sample list rows carry `data-sample-id` and the detector overlay carries `data-exposure-id`. To land on a specific sample and frame, go to its URL: `page.goto("/experiments/1/workspace/10?exposure=2")`. The persisted store holds only the user and tutorial fields; seed those in `localStorage` (`himalaya-ui:state`) before `page.goto`. Any version number works because the store's `migrate` drops retired keys.

See the live specs for the full pattern.

## Port binding

`npm run e2e` builds this checkout's bundle into `dist-e2e/` and serves it on `http://127.0.0.1:5191` (not `localhost`) with `vite preview`. There's no need to start a server first, and a running demo can stay up: the suite touches neither its port nor its `dist/`. It never reuses a running server, because a reused server serves another checkout's bundle, so the suite would pass or fail on code that isn't yours (this hid a missing mock on #347). If :5191 is taken, the run fails at startup; free it with `lsof -ti:5191 | xargs kill -9`.

## Live-integration timing

`page.goto("/")` in live mode races SSE subscription. **Wait ~800ms before any mutation** that expects an SSE echo, otherwise the test browser misses the broadcast. Full guidance in [live/README.md](live/README.md).

## Anti-patterns

- Don't use `localhost` in any URL — Playwright config and CORS expectations all assume `127.0.0.1`.
- Don't `page.route` and `npm run e2e:live` together — live tests must hit the real backend.
- Don't assume CSS-class changes are E2E-stable — only `data-*` and `role` are.
- A `page.route` glob matches the whole URL, query string included, so `"**/api/samples"` never matches `/api/samples?experiment_id=4` (the workspace's sample list). For an endpoint that takes query parameters, match on the path: `page.route((url) => url.pathname === "/api/samples", …)`.
