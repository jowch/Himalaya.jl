import { test, expect, type Page } from "@playwright/test";

// I3.4 (#174): the Phase-3 scoping mocked spec. Exercises /series/new —
// propose an ordering variable from corpus tags, confirm-and-build, batch
// write with source='scoping', then land on the folio.

const EXPERIMENT = {
  id: 1, name: "SSRL Test", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};

function pickerRow(id: number, name: string, ratio: string) {
  return {
    sample: {
      id, experiment_id: 1, name, display_name: null, notes: null,
      tags: [{ id, key: "ratio", value: ratio, source: "manual" }],
    },
    indexing_exposure_id: null,
    all_exposures: [],
  };
}

async function mockScoping(page: Page, captured: unknown[]): Promise<void> {
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPERIMENT]) }));
  await page.route("**/api/sample-tags", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([{ key: "ratio", value: "1:1" }, { key: "ratio", value: "2:1" }]) }));
  await page.route("**/api/picker-samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([pickerRow(10, "A", "1:1"), pickerRow(11, "B", "2:1")]) }));
  await page.route("**/api/samples/tags/batch", async (route) => {
    captured.push(route.request().postDataJSON());
    await route.fulfill({ status: 201, contentType: "application/json",
      body: JSON.stringify([{ id: 1, sample_id: 10, key: "ratio", value: "1:1", source: "scoping" }]) });
  });
  // The folio (post-confirm landing) fetches its listing on mount.
  await page.route("**/api/series", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
}

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => localStorage.setItem("himalaya-ui:state",
    JSON.stringify({ state: { username: "alice", theme: "dark", tutorialSeen: true }, version: 3 })));
});

test("scoping: propose → confirm → batch write with source=scoping → folio", async ({ page }) => {
  const captured: unknown[] = [];
  await mockScoping(page, captured);
  await page.goto("/series/new");
  await expect(page.getByTestId("scoping-page")).toBeVisible();
  await expect(page.getByTestId("ordering-key")).toHaveText("ratio");
  // Build gate opens once both rows seed (both non-flagged).
  await expect(page.getByTestId("scoping-open-confirm")).toBeEnabled();
  await page.getByTestId("scoping-open-confirm").click();
  await expect(page.getByTestId("scoping-confirm-modal")).toBeVisible();
  await page.getByTestId("scoping-confirm-build").click();
  await expect.poll(() => captured[0]).toMatchObject({ key: "ratio", source: "scoping" });
  // D1: lands on the folio after the write (not /series/:id).
  await expect(page).toHaveURL(/\/series$/);
});
