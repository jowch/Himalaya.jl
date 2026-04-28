import { defineConfig } from "@playwright/test";

export default defineConfig({
  testDir: "./e2e",
  timeout: 15_000,
  use: {
    baseURL: "http://127.0.0.1:5173",
    browserName: "chromium",
    headless: true,
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
