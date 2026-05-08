import { defineConfig } from "@playwright/test";

export default defineConfig({
  testDir: "./e2e",
  // Live-server integration tests live under `e2e/live/` and need a real
  // backend + dev DB. Excluded from the default `npm run e2e` (which expects
  // page.route mocks); run with `npm run e2e:live` instead.
  testIgnore: ["**/live/**"],
  timeout: 15_000,
  use: {
    baseURL: "http://127.0.0.1:5173",
    browserName: "chromium",
    headless: true,
    permissions: ["clipboard-read", "clipboard-write"],
  },
  webServer: {
    // `--host 127.0.0.1` is load-bearing: default `npm run dev` binds Vite to
    // `localhost` (sometimes IPv6 ::1) and Playwright's 127.0.0.1 probe times
    // out for 60s. Always bind explicitly here.
    command: "npm run dev -- --host 127.0.0.1",
    url: "http://127.0.0.1:5173",
    reuseExistingServer: !process.env["CI"],
    stdout: "pipe",
    stderr: "pipe",
  },
});
