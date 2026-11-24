import { test, expect } from "./fixtures";

const SAMPLES = [{ id: 41, name: "JC002" }].map((s) => ({ ...s, experiment_id: 17, display_name: s.name, notes: null, q_units: "A-1" }));

const q = Array.from({ length: 200 }, (_, i) => 0.02 + i * 0.0015);
const bump = (x: number, c: number) => 50 * Math.exp(-(((x - c) / 0.002) ** 2));
const trace = { q, I: q.map((x) => 1 / (x * x) + bump(x, 0.1) + bump(x, 0.1414)), sigma: q.map(() => 1) };
const peak = (id: number, at: number) => ({ id, exposure_id: 410, q: at, intensity: 50, prominence: null, sharpness: null, source: "auto", excluded: false });

const exposure = {
  id: 410, sample_id: 41, filename: "sample_S41_0_001", kind: "file", status: null,
  image_path: null, image_version: "", sources: [], trace_hash: null, analysis_inputs_hash: null,
  exposure_time: null, counting_time: null, transmitted_counts: null, incident_monitor: null,
  timestamp: null, beam_energy_ev: null, pipe_length_m: null, detector: null,
};

test.use({ viewport: { width: 1440, height: 900 } });

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.clear();
    localStorage.setItem("himalaya-ui:state", JSON.stringify({ state: { username: "alice", tutorialSeen: true }, version: 6 }));
  });
  const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });
  const experiment = { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "", manifest_path: null, created_at: "", q_units: null };
  await page.route("**/api/experiments", (r) => r.fulfill(json([experiment])));
  await page.route("**/api/experiments/17", (r) => r.fulfill(json(experiment)));
  await page.route((url) => url.pathname === "/api/samples", (r) => r.fulfill(json(SAMPLES)));
  await page.route("**/api/users", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/bridges", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/compares", (r) => r.fulfill(json({ compares: [], subtractions: [] })));
  await page.route("**/api/experiments/17/problems", (r) => r.fulfill(json({ ingest_error: null, exposures: [] })));
  await page.route("**/api/label-colors", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/labels", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/scan", (r) => r.fulfill(json({ status: "complete", experiment_id: 17, changed: false })));
  await page.route(/\/api\/samples\/\d+\/call$/, (r) => r.fulfill(json({ sample_id: 41, state: "indexed", phase: null, call_phases: [], members: [] })));
  await page.route(/\/api\/samples\/\d+\/exposures/, (r) => r.fulfill(json([exposure])));
  await page.route(/\/api\/exposures\/\d+\/trace/, (r) => r.fulfill(json(trace)));
  await page.route(/\/api\/exposures\/\d+\/peaks/, (r) => r.fulfill(json([peak(1, 0.1), peak(2, 0.1414)])));
  await page.route(/\/api\/exposures\/\d+\/indices/, (r) => r.fulfill(json([])));
  await page.route(/\/api\/exposures\/\d+\/assignment/, (r) => r.fulfill(json({ exposure_id: 410, state: "indexed", members: [] })));
  await page.route("**/api/events", (r) => r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
});

async function recordCustomIndex(page: import("@playwright/test").Page): Promise<unknown[]> {
  const posted: unknown[] = [];
  await page.route("**/api/exposures/410/custom-index", (r) => {
    posted.push(r.request().postDataJSON());
    return r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify({
      id: 9001, exposure_id: 410, phase: "Pn3m", basis: 0.0707, score: null, r_squared: null, lattice_d: 88.9,
      status: "candidate", kind: "speculative", inputs_hash: null, peaks: [], predicted_q: [], event_id: 1, view_row_id: null,
    }) });
  });
  return posted;
}

async function openCustomIndex(page: import("@playwright/test").Page): Promise<void> {
  await page.goto("/experiments/17/workspace/41");
  await expect(page.locator('[data-peak-key="auto:1"]')).toBeVisible();
  await page.keyboard.press("i");
  await expect(page.getByTestId("workspace-custom-index")).toBeVisible();
}

async function snapAndAdd(page: import("@playwright/test").Page, posted: unknown[]): Promise<void> {
  await page.locator('[data-peak-key="auto:1"]').click();
  await expect(page.locator('[data-peak-key="auto:1"]')).not.toBeFocused();
  await page.keyboard.press("Enter");
  await expect(page.getByTestId("workspace-custom-index")).toBeHidden();
  expect(posted).toHaveLength(1);
}

test("clicking a peak to snap the custom index, then Enter, adds it", async ({ page }) => {
  const posted = await recordCustomIndex(page);
  await openCustomIndex(page);
  await snapAndAdd(page, posted);
});

test("picking a symmetry, clicking a peak to snap, then Enter, adds it", async ({ page }) => {
  const posted = await recordCustomIndex(page);
  await openCustomIndex(page);
  await page.getByRole("group", { name: "Symmetry" }).getByRole("button").nth(1).click();
  await snapAndAdd(page, posted);
});

test("with the filter focused, clicking a peak to snap, then Enter, adds it", async ({ page }) => {
  const posted = await recordCustomIndex(page);
  await openCustomIndex(page);
  await page.getByLabel("Filter samples").focus();
  await snapAndAdd(page, posted);
});

test("a value typed in the lattice field gives way to the peak the click snaps to", async ({ page }) => {
  const posted = await recordCustomIndex(page);
  await openCustomIndex(page);
  const field = page.getByRole("spinbutton", { name: /Lattice . value/ });
  await field.fill("150");
  await page.locator('[data-peak-key="auto:1"]').click();
  await expect(field).not.toHaveValue("150");
  const snapped = await field.inputValue();
  await page.keyboard.press("Enter");
  await expect(page.getByTestId("workspace-custom-index")).toBeHidden();
  expect(posted).toHaveLength(1);
  expect(snapped).not.toBe("150.0");
});
