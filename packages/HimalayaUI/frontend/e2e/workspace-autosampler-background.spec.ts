import { test, expect } from "./fixtures";
import { expectButtonRowsAligned } from "./buttonRows";

const FIRST_BUFFER = "wat 01A (B0101)";
const SAMPLE = "wat 01B (S0102)";
const CLOSING_BUFFER = "wat 01H (B0103)";

const SAMPLES = [
  { id: 50, name: FIRST_BUFFER, slot: 101 },
  { id: 51, name: SAMPLE, slot: 102 },
  { id: 52, name: CLOSING_BUFFER, slot: 103 },
].map((s) => ({
  id: s.id, experiment_id: 17, name: s.name, display_name: s.name, notes: null, q_units: "A-1", setup_id: 1,
  analyzed_exposure_id: s.id * 10, assignment_state: "indexed", phase: null, call_phases: [], slot_index: s.slot,
}));

const q = Array.from({ length: 200 }, (_, i) => 0.02 + i * 0.0015);
const trace = (level: number) => ({ q, I: q.map((x) => level + 40 * Math.exp(-(((x - 0.1) / 0.03) ** 2))), sigma: q.map(() => 0.5) });

const PNG = Buffer.from("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==", "base64");

interface Background { sample_id: number; ref_sample_id: number | null; source: "default" | "user"; default_ref_sample_id: number | null; is_buffer: boolean; run_buffer_ids: number[] }

function backgrounds(chosen: Map<number, number | null>): Background[] {
  return SAMPLES.map((s) => {
    const buffer = s.id !== 51;
    const fallback = buffer ? null : 50;
    const user = chosen.has(s.id);
    return { sample_id: s.id, ref_sample_id: user ? chosen.get(s.id)! : fallback, source: user ? "user" : "default", default_ref_sample_id: fallback, is_buffer: buffer, run_buffer_ids: [50, 52] };
  });
}

test.use({ viewport: { width: 1440, height: 900 } });

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.clear();
    localStorage.setItem("himalaya-ui:state", JSON.stringify({ state: { username: "alice", tutorialSeen: true }, version: 6 }));
  });
  const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });
  const setup = {
    id: 1, experiment_id: 17, folder: "autosampler", display_name: null, name: "autosampler", position: 0, included: true,
    mode: "autosampler", mode_source: "detected", root_dir: "/demo/autosampler", data_dir: "/demo/autosampler/data",
    analysis_dir: "/demo/autosampler/analysis", image_pattern: "{name}_0_001.tif", metadata_pattern: "{name}_0_001.prp",
    integration_pattern: "{name}_tot.dat", q_units: "A-1", beam_center_x: null, beam_center_y: null, pixel_size_um: null,
    energy_kev: null, flight_path_m: null, energy_kev_source: "default", flight_path_m_source: "default",
    beam_center_x_source: "default", beam_center_y_source: "default", pixel_size_um_source: "default",
    trace_normalisation: "i2", normalisation_base: 2.5e7,
  };
  const experiment = { id: 17, name: "solution", path: "/demo", data_dir: "/demo/autosampler/data", analysis_dir: "", manifest_path: null, created_at: "", q_units: null, ingest_status: "complete", setups: [setup] };
  const chosen = new Map<number, number | null>();
  const state = () => ({ compares: [], subtractions: [], backgrounds: backgrounds(chosen) });

  await page.route("**/api/experiments", (r) => r.fulfill(json([experiment])));
  await page.route("**/api/experiments/17", (r) => r.fulfill(json(experiment)));
  await page.route((url) => url.pathname === "/api/samples", (r) => r.fulfill(json(SAMPLES)));
  await page.route("**/api/users", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/bridges", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/compares", (r) => r.fulfill(json(state())));
  await page.route(/\/api\/samples\/\d+\/compares$/, (r) => {
    const id = Number(new URL(r.request().url()).pathname.split("/")[3]);
    const change = r.request().postDataJSON() as { subtraction?: number | null | "default" };
    if (change.subtraction === "default") chosen.delete(id);
    else if (change.subtraction !== undefined) chosen.set(id, change.subtraction);
    return r.fulfill(json(state()));
  });
  await page.route("**/api/experiments/17/problems", (r) => r.fulfill(json({ ingest_error: null, exposures: [] })));
  await page.route("**/api/experiments/17/labels", (r) => r.fulfill(json([])));
  await page.route("**/api/label-colors", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/scan", (r) => r.fulfill(json({ status: "complete", experiment_id: 17, changed: false })));
  await page.route(/\/api\/samples\/\d+\/call$/, (r) => {
    const id = Number(/samples\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill(json({ sample_id: id, state: "indexed", phase: null, call_phases: [], members: [] }));
  });
  await page.route(/\/api\/samples\/\d+\/exposures/, (r) => {
    const id = Number(/samples\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill(json([{
      id: id * 10, sample_id: id, filename: `wat_S${id}`, kind: "file", status: null, selected: true, image_path: null, image_version: "v1",
      sources: [], trace_hash: null, analysis_inputs_hash: null, exposure_time: 1, counting_time: 1.8, transmitted_counts: 3.0e6,
      incident_monitor: 1.5, timestamp: "2026-04-27T09:31:00", beam_energy_ev: 9000, pipe_length_m: 1.7, detector: "Pilatus 1M",
      frames_total: 10, frames_rejected: 1, setup_id: 1,
    }]));
  });
  await page.route(/\/api\/exposures\/\d+\/trace/, (r) => {
    const id = Number(/exposures\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill(json({ ...trace(id === 510 ? 140 : 100), monitors: { transmitted_counts: 3.0e6, exposure_time: 1, counting_time: 1.8, normalised_by: "i2", beam_off: false, i2_collapsed: false, mode: "autosampler" } }));
  });
  await page.route(/\/api\/exposures\/\d+\/image/, (r) => r.fulfill({ status: 200, contentType: "image/png", body: PNG, headers: { "X-Image-Width": "1", "X-Image-Height": "1" } }));
  await page.route(/\/api\/exposures\/\d+\/(peaks|indices)/, (r) => r.fulfill(json([])));
  await page.route(/\/api\/exposures\/\d+\/assignment/, (r) => {
    const id = Number(/exposures\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill(json({ exposure_id: id, state: "indexed", members: [] }));
  });
  await page.route("**/api/events", (r) => r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
});

test("an autosampler sample opens in Minus on its run's buffer, and the Background module changes, clears and restores it", async ({ page }) => {
  await page.goto("/experiments/17/workspace/51");
  const name = page.getByTestId("background-name");
  const difference = page.locator("[data-role='difference-line']");
  const minus = (buffer: string) => page.getByRole("button", { name: `Minus ${buffer}` });

  await expect(name).toContainText(FIRST_BUFFER);
  await expect(name).toContainText("run's buffer");
  await expect(difference).toHaveCount(1);
  await expect(minus(FIRST_BUFFER)).toHaveAttribute("aria-pressed", "true");
  await expect(page.getByTestId("workspace-candidates")).toHaveCount(0);
  await expect(page.getByTestId("workspace-custom-index")).toHaveCount(0);
  await expectButtonRowsAligned(page, "Change and None");

  await page.getByRole("button", { name: /^Change/ }).click();
  await page.getByRole("dialog", { name: "Choose the background" }).locator(`[data-sample-id='52']`).click();
  await expect(name).toContainText(CLOSING_BUFFER);
  await expect(name).toContainText("chosen");
  await expect(name).not.toContainText("run's buffer");
  await expect(minus(CLOSING_BUFFER)).toHaveAttribute("aria-pressed", "true");
  await expect(difference).toHaveCount(1);

  await page.getByRole("button", { name: "None", exact: true }).click();
  await expect(name).toHaveText("None");
  await expect(difference).toHaveCount(0);
  await expect(page.getByRole("button", { name: "Run's buffer" })).toBeVisible();
  await expectButtonRowsAligned(page, "Change and Run's buffer");

  await page.getByRole("button", { name: "Run's buffer" }).click();
  await expect(name).toContainText(FIRST_BUFFER);
  await expect(name).toContainText("run's buffer");
  await expect(minus(FIRST_BUFFER)).toHaveAttribute("aria-pressed", "true");
  await expect(difference).toHaveCount(1);
});

test("a buffer says so and offers no background", async ({ page }) => {
  await page.goto("/experiments/17/workspace/50");
  await expect(page.getByTestId("background-name")).toHaveText("This is a buffer");
  await expect(page.getByRole("button", { name: /^Change/ })).toHaveCount(0);
  await expect(page.locator("[data-role='difference-line']")).toHaveCount(0);
});
