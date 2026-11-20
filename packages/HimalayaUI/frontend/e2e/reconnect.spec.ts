import { test, expect, type Page } from "./fixtures";

const experiment = (name: string) => ({
  id: 3, name, path: "/data", data_dir: "/data", analysis_dir: "/analysis",
  manifest_path: null, created_at: "2026-01-01", last_scanned_at: "2026-01-01T10:00:00",
  ingest_status: "complete", stats: { loads: 1, samples: 1, exposures: 1, sessions: 1 },
});

type Answer = { status: number; body: string; contentType: string };
const stream = (body: string): Answer => ({ status: 200, contentType: "text/event-stream", body });
const KEEP_QUIET = stream("retry: 600000\n\n");

async function mockApp(page: Page, answers: Answer[]): Promise<{ connections: () => number }> {
  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({ state: { username: "alice", tutorialSeen: true }, version: 5 }));
  });
  let connections = 0;
  let listed = 0;
  await page.route("**/api/users", (r) => r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/samples", (r) => r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) => {
    listed += 1;
    const afterRestart = listed > 1 && connections > 1;
    return r.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify([experiment(afterRestart ? "Renamed while the stream was down" : "Before the restart")]),
    });
  });
  await page.route("**/api/events", (r) => {
    connections += 1;
    return r.fulfill(answers[connections - 1] ?? KEEP_QUIET);
  });
  return { connections: () => connections };
}

test("after the stream ends and the browser reconnects, the page shows what changed while it was down", async ({ page }) => {
  await mockApp(page, [stream("retry: 1\n\n")]);
  await page.goto("/experiments");
  await expect(page.getByText("Renamed while the stream was down")).toBeVisible();
});

test("after the proxy answers 500 while the backend is down, the app opens a new stream and refreshes", async ({ page }) => {
  const app = await mockApp(page, [
    stream("retry: 1\n\n"),
    { status: 500, contentType: "text/plain", body: "upstream down" },
  ]);
  await page.goto("/experiments");
  await expect.poll(app.connections).toBeGreaterThanOrEqual(3);
  await expect(page.getByText("Renamed while the stream was down")).toBeVisible();
});
