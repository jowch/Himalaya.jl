# Live integration tests

Playwright E2E tests under this directory hit a **real backend + dev DB**.
They're distinct from the default mocked `npm run e2e` suite (which
intercepts `/api/*` via `page.route` and never starts a Julia process).

Use this category for any check that needs SSE, real DB state transitions,
cross-process atomicity, or any other behaviour `page.route` can't simulate.

```
npm run e2e:live    # uses playwright.live.config.ts
```

---

## Operator workflow

`playwright.live.config.ts` deliberately omits `webServer` so it can't
accidentally proxy to the wrong backend. You bring up the backend and Vite
manually before running the tests:

```bash
# Terminal 1 — backend on port 8090, test DB
HIMALAYA_DB_PATH=/tmp/himalaya-test/himalaya.db \
  julia --project=packages/HimalayaUI -e 'using HimalayaUI; HimalayaUI.serve(open_db("/tmp/himalaya-test/himalaya.db"); port=8090)'

# Terminal 2 — Vite on port 5180, proxying /api → :8090
(cd packages/HimalayaUI/frontend && \
   VITE_BACKEND_PORT=8090 npm run dev -- --host 127.0.0.1 --port 5180)

# Terminal 3 — run the live tests
(cd packages/HimalayaUI/frontend && npm run e2e:live)
```

The test DB lives at `/tmp/himalaya-test/himalaya.db` (a copy of prod
data) so prod's 8080 / 5173 / `data/himalaya.db` stays untouched.

Override either side via env vars:

- `PLAYWRIGHT_BASE_URL` — where Playwright drives the browser (default `http://127.0.0.1:5180`)
- `BACKEND_BASE` — where the backend is reachable (default `http://127.0.0.1:8090`)

---

## Reset between flows

```bash
bash packages/HimalayaUI/scripts/reset-test-db.sh
```

This re-copies `data/himalaya.db` over the test DB and bounces the backend
(takes ~5s end-to-end). Required for tests that need a clean slate on the
same exposure — e.g. `post-pr37.spec.ts`'s "Bug 1c" run needs no existing
custom group and would false-pass on a re-run otherwise.

The reset script also **warms Vite's proxy connection pool**. Without that
warmup, the next page load hangs ~30s while Vite retries dead keepAlive
connections to the old backend process. The warmup is built into the
script, so you don't have to think about it — just run it before any test
that needs a fresh DB.

---

## Why `workers: 1`

`playwright.live.config.ts` pins `workers: 1`. Live tests share one
backend; parallel mutations on the same row collide on `with_idempotency`
locks or on FK enforcement during teardown. Running serially trades
wall-clock for determinism, which is the right trade for an integration
suite.

The mocked `npm run e2e` suite has no such constraint and runs parallel.

---

## EventSource warmup

After `page.goto("/")`, **wait ~800 ms before triggering a mutation that
expects an SSE echo.**

Why: the server-side SSE subscriber registration completes ~50–200 ms
after the GET (the backend's HTTP handler returns first; the streaming
response keeps writing). A click that fires before the subscriber lands
in `SSE_SUBSCRIBERS[]` will broadcast `post_state` to no one, the test
browser misses the indices update, and `StaleIndicesBanner` sticks until
the next polling refetch.

Pattern:

```ts
await page.goto("/");
await page.waitForTimeout(800);   // SSE warmup
await page.click('[data-testid="add-peak"]');
await expect(page.locator(".stale-banner")).toHaveCount(0);
```

The mocked suite has no such warmup because there's no real subscriber to
register.

---

## Vite ignores `e2e/` in its watcher

`vite.config.ts` adds `e2e/**` to the watcher's `ignored` list. Editing a
spec during a live Playwright run would otherwise trigger an HMR
page-reload that crashes the dev server mid-test.

If you add a *new* test directory at the same level (e.g. an
`integration/` sibling), add it to the ignore list too.

---

## What to test here vs. in the mocked suite

**Use live integration when the test depends on:**

- SSE multiplayer (cross-tab / cross-user fanout)
- Real DB state transitions (autoincrement id, FK cascades, AUTOINCREMENT
  migration paths)
- `with_idempotency` cache hit/miss behaviour
- `analyze_exposure!` running for real (peak detection, index scoring)
- Cross-process atomicity (e.g. the route emit ↔ broadcast ordering)

**Use the mocked suite when:**

- The test only needs UI behaviour given a known API response
- The test focuses on form validation, modal flow, focus management
- You want fast, parallel, deterministic regression coverage in CI

---

## Existing specs in this directory

| Spec | What it covers |
|---|---|
| `peak-add-no-stale-banner.spec.ts` | Synchronous reanalyze closes the stale-banner window |
| `post-pr37.spec.ts` | Multiple bug fixes from PR #37 review (queue ordering, custom group reset) |
| `sample-rename-preserves-fields.spec.ts` | Sample rename doesn't blank other metadata fields |
| `sample-tag-add-and-delete.spec.ts` | Tag CRUD round-trip via SSE |
| `speculative-create.spec.ts` | Speculative index create + delete |

Add new specs as `<feature>.spec.ts`. Keep one feature per file — the
serial worker means each spec adds linearly to wall-clock.

---

## Files

- `packages/HimalayaUI/frontend/playwright.live.config.ts` — config (no webServer, workers: 1)
- `packages/HimalayaUI/scripts/reset-test-db.sh` — DB reset + Vite proxy warmup
- `packages/HimalayaUI/frontend/vite.config.ts` — `e2e/` ignore rule
