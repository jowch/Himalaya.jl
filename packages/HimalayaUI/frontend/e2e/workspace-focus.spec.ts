import { test, expect } from "./fixtures";

const SAMPLES = ["JC001", "JC002", "JC003", "BUF001"].map((name, i) => ({
  id: 40 + i, experiment_id: 17, display_name: name, name,
  notes: null, q_units: "A-1",
}));

function exposureOf(sampleId: number, frame: number) {
  return {
    id: sampleId * 10 + frame, sample_id: sampleId, filename: `sample_S${sampleId}_0_00${frame + 1}`,
    kind: "file", status: null, image_path: null, image_version: "", sources: [], trace_hash: null,
    analysis_inputs_hash: null, exposure_time: null, counting_time: null, transmitted_counts: null,
    incident_monitor: null, timestamp: null, beam_energy_ev: null, pipe_length_m: null, detector: null,
  };
}

const focusVisible = (page: import("@playwright/test").Page) =>
  page.evaluate(() => document.querySelector(":focus-visible")?.textContent?.trim() ?? null);

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.clear();
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: { username: "alice", tutorialSeen: true },
      version: 6,
    }));
  });
  const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });
  const experiment = { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "",
    manifest_path: null, created_at: "", q_units: null };
  await page.route("**/api/experiments", (r) => r.fulfill(json([experiment])));
  await page.route("**/api/experiments/17", (r) => r.fulfill(json(experiment)));
  await page.route((url) => url.pathname === "/api/samples", (r) => r.fulfill(json(SAMPLES)));
  await page.route("**/api/users", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/bridges", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/compares", (r) => r.fulfill(json({ compares: [], subtractions: [] })));
  await page.route("**/api/label-colors", (r) => r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments/17/scan", (r) =>
    r.fulfill(json({ status: "complete", experiment_id: 17, changed: false })));
  await page.route("**/api/experiments/17/labels", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/problems", (r) => r.fulfill(json({ ingest_error: null, exposures: [] })));
  await page.route(/\/api\/samples\/\d+\/call$/, (r) => {
    const id = Number(/samples\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill(json({ sample_id: id, state: "indexed", phase: null, call_phases: [], members: [] }));
  });
  await page.route("**/api/samples/*/exposures*", (r) => {
    const id = Number(/samples\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill(json([0, 1, 2].map((f) => exposureOf(id, f))));
  });
  await page.route(/\/api\/exposures\/\d+\/trace/, (r) => r.fulfill(json({ q: [0.02, 0.2], I: [10, 1], sigma: [1, 1] })));
  await page.route(/\/api\/exposures\/\d+\/(peaks|indices)/, (r) => r.fulfill(json([])));
  await page.route(/\/api\/exposures\/\d+\/assignment/, (r) => {
    const eid = Number(/exposures\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill(json({ exposure_id: eid, state: "indexed", members: [] }));
  });
  await page.route("**/api/events", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
});

test("a clicked sample row shows no focus ring after the arrow keys move the selection", async ({ page }) => {
  await page.goto("/experiments/17/workspace/40");
  const list = page.getByTestId("workspace-samples");
  await list.getByRole("button", { name: /JC002/ }).click();
  await expect(page).toHaveURL(/\/workspace\/41$/);
  await page.keyboard.press("ArrowDown");
  await expect(page).toHaveURL(/\/workspace\/42$/);
  await page.keyboard.press("ArrowUp");
  await page.keyboard.press("ArrowUp");
  await expect(page).toHaveURL(/\/workspace\/40$/);
  expect(await focusVisible(page)).toBeNull();
});

test("keyboard focus in the sample list follows the selection", async ({ page }) => {
  await page.goto("/experiments/17/workspace/40");
  const list = page.getByTestId("workspace-samples");
  await list.getByRole("button", { name: /JC001/ }).focus();
  await page.keyboard.press("ArrowDown");
  await expect(page).toHaveURL(/\/workspace\/41$/);
  await expect(list.getByRole("button", { name: /JC002/ })).toBeFocused();
  expect(await focusVisible(page)).toMatch(/JC002/);
});

test("a clicked image shows no focus ring after the arrow keys move to another image", async ({ page }) => {
  await page.goto("/experiments/17/workspace/40");
  const strip = page.getByTestId("workspace-exposures");
  await strip.getByRole("button", { name: "Image 2" }).click();
  await page.keyboard.press("ArrowRight");
  await expect(strip.getByRole("button", { name: "Image 3" })).toHaveAttribute("aria-current", "true");
  expect(await page.evaluate(() => document.querySelector(":focus-visible")?.getAttribute("aria-label") ?? null)).toBeNull();
});
