import { test, expect } from "./fixtures";

const SAMPLES = ["JC001", "JC002", "JC003", "BUF001"].map((name, i) => ({
  id: 40 + i, experiment_id: 17, display_name: name, name,
  notes: null, q_units: "A-1",
}));

function exposureOf(sampleId: number, frame: number) {
  return {
    id: sampleId * 10 + frame, sample_id: sampleId, filename: `JC_S${sampleId}_0_00${frame + 1}`,
    kind: "file", status: null, image_path: null, image_version: "", sources: [], trace_hash: null,
    analysis_inputs_hash: null, exposure_time: null, counting_time: null, transmitted_counts: null,
    incident_monitor: null, timestamp: null, beam_energy_ev: null, pipe_length_m: null, detector: null,
  };
}

const ring = (page: import("@playwright/test").Page) =>
  page.evaluate(() => {
    const el = document.activeElement;
    if (!(el instanceof HTMLElement) || el === document.body) return null;
    return { text: el.textContent?.trim() ?? "", outline: getComputedStyle(el).outlineStyle };
  });

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
  await page.route(/\/api\/exposures\/\d+\/peaks/, (r) => r.fulfill(json([])));
  await page.route(/\/api\/exposures\/\d+\/indices/, (r) => {
    const eid = Number(/exposures\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill(json(["Pn3m", "Im3m", "Ia3d"].map((phase, i) => ({
      id: i + 1, exposure_id: eid, phase, basis: 0.1, score: 0.9 - i / 10, r_squared: 0.99, lattice_d: 100 + i,
      status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [0.1, 0.141],
    }))));
  });
  await page.route(/\/api\/exposures\/\d+\/assignment/, (r) => {
    const eid = Number(/exposures\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill(json({ exposure_id: eid, state: "indexed", members: [] }));
  });
  await page.route("**/api/events", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
});

test("a phase choice clicked with the mouse shows no ring after →, and Tab brings it back", async ({ page }) => {
  await page.goto("/experiments/17/workspace/40");
  const choices = page.getByRole("group", { name: "Phase choices" });
  const im3m = choices.getByRole("button", { name: /Im3m/ });
  await im3m.click();
  await expect(im3m).toBeFocused();
  await expect(page.getByTestId("workspace-exposures").getByRole("button", { name: "Image 3" })).toHaveAttribute("aria-current", "true");
  await page.keyboard.press("ArrowRight");
  expect(await ring(page)).toMatchObject({ text: expect.stringMatching(/Im3m/), outline: "none" });
  await page.keyboard.press("Tab");
  expect(await ring(page)).toMatchObject({ text: expect.stringMatching(/Ia3d/), outline: "solid" });
});

test("a phase choice reached with Tab shows the ring", async ({ page }) => {
  await page.goto("/experiments/17/workspace/40");
  await expect(page.getByRole("group", { name: "Phase choices" }).getByRole("button", { name: /Pn3m/ })).toBeVisible();
  await page.keyboard.press("Tab");
  expect(await ring(page)).toMatchObject({ text: expect.stringMatching(/Pn3m/), outline: "solid" });
});
