import { defineConfig } from "@playwright/test";

const E2E_PORT = 5191;

export default defineConfig({
  testDir: "./e2e",
  timeout: 15_000,
  forbidOnly: !!process.env["CI"],
  retries: 0,
  reporter: process.env["CI"] ? [["list"], ["html", { open: "never" }]] : "list",
  use: {
    baseURL: `http://127.0.0.1:${E2E_PORT}`,
    browserName: "chromium",
    headless: true,
  },
  webServer: {
    command: `vite build --outDir dist-e2e && vite preview --outDir dist-e2e --host 127.0.0.1 --port ${E2E_PORT} --strictPort`,
    url: `http://127.0.0.1:${E2E_PORT}`,
    reuseExistingServer: false,
    stdout: "pipe",
    stderr: "pipe",
  },
});
