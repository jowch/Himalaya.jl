import { defineConfig } from "@playwright/test";

/**
 * Live-server integration tests. These hit a real backend and dev DB —
 * NO page.route mocks. Run with `npm run e2e:live` after bringing up the
 * backend and Vite manually.
 *
 * Default port pair: Vite on 5180, backend on 8090, DB at
 * /tmp/himalaya-test/himalaya.db (a copy of prod data, NOT prod itself).
 * This pair is deliberately offset from prod's 8080/5173 so a stale prod
 * serve can keep running on the same machine without colliding.
 *
 *   HIMALAYA_DB_PATH=/tmp/himalaya-test/himalaya.db \
 *     julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' \
 *     -- serve --port 8090 --host 127.0.0.1
 *   (cd packages/HimalayaUI/frontend && \
 *      VITE_API_PORT=8090 npm run dev -- --host 127.0.0.1 --port 5180 --strictPort)
 *
 * Override either side via env (useful when a worktree wants its own pair):
 *   PLAYWRIGHT_BASE_URL=http://127.0.0.1:5181 npm run e2e:live
 *   BACKEND_BASE=http://127.0.0.1:8091 npm run e2e:live
 *
 * Distinct from the default `playwright.config.ts` (mocked-route Playwright
 * suite) because:
 *   - Live tests need state in a real DB; not portable to CI's empty
 *     fixture directory.
 *   - `webServer` is intentionally omitted — the operator brings up Vite +
 *     backend with a known port mapping (VITE_API_PORT). Letting
 *     Playwright auto-start `npm run dev` would silently proxy to the
 *     wrong backend.
 *   - Timeout is longer because real network + DB queries dominate.
 */
const BASE_URL = process.env["PLAYWRIGHT_BASE_URL"] ?? "http://127.0.0.1:5180";

export default defineConfig({
  testDir: "./e2e/live",
  timeout: 30_000,
  // Workers=1: live tests share a single backend + DB. Parallel mutations
  // collide (one test's tag-add interferes with another's sample patch on
  // the same row). The deterministic mocked suite still runs in parallel.
  workers: 1,
  use: {
    baseURL: BASE_URL,
    browserName: "chromium",
    headless: true,
  },
});
