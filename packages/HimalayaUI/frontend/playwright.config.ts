import { defineConfig } from "@playwright/test";

export default defineConfig({
  testDir: "./e2e",
  // Live-server integration tests live under `e2e/live/` and need a real
  // backend + dev DB. Excluded from the default `npm run e2e` (which expects
  // page.route mocks); run with `npm run e2e:live` instead.
  testIgnore: ["**/live/**"],
  timeout: 15_000,
  forbidOnly: !!process.env["CI"],
  retries: 0,
  reporter: process.env["CI"] ? [["list"], ["html", { open: "never" }]] : "list",
  use: {
    baseURL: "http://127.0.0.1:5173",
    browserName: "chromium",
    headless: true,
  },
  webServer: {
    command: "vite build && vite preview --host 127.0.0.1 --port 5173 --strictPort",
    url: "http://127.0.0.1:5173",
    reuseExistingServer: !process.env["CI"],
    stdout: "pipe",
    stderr: "pipe",
  },
});
