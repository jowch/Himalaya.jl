import { defineConfig } from "@playwright/test";

/**
 * Live-server integration tests. These hit a real backend and dev DB —
 * NO page.route mocks. Run with `npm run e2e:live` after bringing up the
 * backend and Vite manually:
 *
 *   bin/himalaya serve <experiment-dir> --port 8080
 *   (cd packages/HimalayaUI/frontend && \
 *      VITE_API_PORT=8080 npm run dev -- --host 127.0.0.1)
 *
 * Distinct from the default `playwright.config.ts` (mocked-route Playwright
 * suite) because:
 *   - Live tests need state in a real DB; not portable to CI's empty
 *     fixture directory.
 *   - `webServer` is intentionally omitted — the receiving operator brings
 *     up Vite + backend with a known port mapping (VITE_API_PORT). Letting
 *     Playwright auto-start `npm run dev` would silently proxy to the
 *     wrong backend (prod's :8080).
 *   - Timeout is longer because real network + DB queries dominate.
 */
export default defineConfig({
  testDir: "./e2e/live",
  timeout: 30_000,
  use: {
    baseURL: "http://127.0.0.1:5173",
    browserName: "chromium",
    headless: true,
  },
});
