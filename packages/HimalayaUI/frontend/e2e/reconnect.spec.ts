import { test, expect } from "./fixtures";

const experiment = (name: string) => ({
  id: 3, name, path: "/data", data_dir: "/data", analysis_dir: "/analysis",
  manifest_path: null, created_at: "2026-01-01", last_scanned_at: "2026-01-01T10:00:00",
  ingest_status: "complete", stats: { loads: 1, samples: 1, exposures: 1, sessions: 1 },
});

test("after the live-update stream drops and reconnects, the page shows what changed while it was down", async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({ state: { username: "alice", tutorialSeen: true }, version: 5 }));
  });
  let restarted = false;
  await page.route("**/api/users", (r) => r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/samples", (r) => r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify([experiment(restarted ? "Renamed while the stream was down" : "Before the restart")]),
    }));
  let connections = 0;
  await page.route("**/api/events", (r) => {
    connections += 1;
    return r.fulfill({
      status: 200, contentType: "text/event-stream",
      body: connections === 1 ? "retry: 300\n\n" : "retry: 600000\n\n",
    });
  });

  await page.goto("/experiments");
  await expect(page.getByText("Before the restart")).toBeVisible();
  restarted = true;
  await expect.poll(() => connections).toBeGreaterThanOrEqual(2);
  await expect(page.getByText("Renamed while the stream was down")).toBeVisible();
});
