# frontend/e2e — Playwright E2E

Two suites:

- **`e2e/`** — default mocked suite. `page.route` intercepts `/api/...`. Run via `npm run e2e`.
- **`e2e/live/`** — live-integration. Real backend + dev DB. Run via `npm run e2e:live`. Operator brings up backend (port 8090) + Vite (5180) manually. Full runbook: [live/README.md](live/README.md).

## Selectors

Use stable `data-*` attributes, `data-testid`, or `role`. **Never assert on Tailwind class strings** — they change when styling evolves. Stable testid families: `data-index-id`, `data-alternative-id`, `data-active`, `data-low-r2`, `data-speculative`.

Single test by name:

```bash
node_modules/.bin/playwright test --grep "Reject → Other"
```

## Seeding Zustand state

There is no `data-sample-id` / `data-exposure-id`. To navigate to a specific sample/exposure, seed Zustand state in `localStorage` before `page.goto("/")`:

```ts
await page.addInitScript((state) => {
  localStorage.setItem("himalaya-ui:state", JSON.stringify({ state, version: 3 }));
}, { activeSampleId: 1, activeExposureId: 2 });
await page.goto("/");
```

See the live specs for the full pattern.

## Port binding

`playwright.config.ts` expects the dev server on `http://127.0.0.1:5173`, **not `localhost`**. If another process holds that port, tests hang for 60s then fail:

```bash
lsof -ti:5173 | xargs kill -9
```

If starting Vite separately before `npm run e2e`, bind it explicitly:

```bash
npm run dev -- --host 127.0.0.1
```

## Live-integration timing

`page.goto("/")` in live mode races SSE subscription. **Wait ~800ms before any mutation** that expects an SSE echo, otherwise the test browser misses the broadcast. Full guidance in [live/README.md](live/README.md).

## Anti-patterns

- Don't use `localhost` in any URL — Playwright config and CORS expectations all assume `127.0.0.1`.
- Don't `page.route` and `npm run e2e:live` together — live tests must hit the real backend.
- Don't assume CSS-class changes are E2E-stable — only `data-*` and `role` are.
