import { test, expect, type Page } from "./fixtures";

const LONG = "a_sample_with_a_long_descriptive_name_S24";
const SAMPLES = [
  { id: 40, name: "JC001" },
  { id: 41, name: "JC002" },
  { id: 42, name: `${LONG}01` },
  { id: 43, name: `${LONG}02_${LONG}_${LONG}` },
  { id: 44, name: `${LONG}03` },
].map((s) => ({ ...s, experiment_id: 17, display_name: s.name, notes: null, tags: [], q_units: "A-1" }));

const q = Array.from({ length: 200 }, (_, i) => 0.02 + i * 0.0015);
const trace = { q, I: q.map((x) => 1 / (x * x) + 50 * Math.exp(-(((x - 0.1) / 0.002) ** 2))), sigma: q.map(() => 1) };

function exposureOf(sampleId: number, monitors: boolean, frame = 0, averaged = false) {
  return {
    id: sampleId * 10 + frame, sample_id: sampleId, filename: `JC_S${sampleId}_0_00${frame + 1}`,
    kind: averaged ? "averaged" : "file", selected: true, status: null,
    image_path: null, image_version: "", tags: [], sources: averaged ? [sampleId * 10] : [], trace_hash: null, analysis_inputs_hash: null,
    exposure_time: monitors ? 15 : null, counting_time: monitors ? 15 : null,
    transmitted_counts: monitors ? 7000 : null, incident_monitor: monitors ? 1.5 : null,
    timestamp: monitors ? "2026-04-26T18:20:27" : null, beam_energy_ev: monitors ? 11000 : null,
    pipe_length_m: monitors ? 1.7 : null, detector: monitors ? "Pilatus 1M" : null,
  };
}

const PHASES = ["Pn3m", "Im3m", "Ia3d", "HII", "Lamellar", "Fd3m"];
function indicesOf(exposureId: number, n: number) {
  return PHASES.slice(0, n).map((phase, i) => ({
    id: exposureId * 10 + i, exposure_id: exposureId, phase, basis: 0.1, score: 0.9 - i * 0.1, r_squared: 0.99,
    lattice_d: 100 + i, ngc: null, status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [0.1, 0.141],
  }));
}

const CANDIDATES: Record<number, number> = { 40: 1, 41: 6, 42: 3, 43: 0, 44: 2 };
const FRAMES: Record<number, number> = { 44: 7 };
const TAGS = ["buffer", "37C", "heated twice", "retake", "dilute 1 to 10", "cloudy", "after spin", "second capillary", "check flare", "low counts"];
const LABELS = [
  ...TAGS.slice(0, 2).map((label, i) => ({ id: i + 1, sample_id: 42, label, state: "confirmed", source: "user" })),
  ...TAGS.map((label, i) => ({ id: i + 11, sample_id: 44, label, state: i === 0 ? "suggested" : "confirmed", source: "user" })),
];
const VIEWPORTS = [
  { width: 1440, height: 900 },
  { width: 1366, height: 650 },
];

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.clear();
    localStorage.setItem("himalaya-ui:state", JSON.stringify({ state: { username: "alice", tutorialSeen: true }, version: 5 }));
  });
  const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });
  const experiment = { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "", manifest_path: null, created_at: "", q_units: null };
  await page.route("**/api/experiments", (r) => r.fulfill(json([experiment])));
  await page.route("**/api/experiments/17", (r) => r.fulfill(json(experiment)));
  await page.route((url) => url.pathname === "/api/samples", (r) => r.fulfill(json(SAMPLES)));
  await page.route("**/api/users", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/bridges", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/problems", (r) => r.fulfill(json({ ingest_error: null, exposures: [] })));
  await page.route("**/api/experiments/17/labels", (r) => r.fulfill(json(LABELS)));
  await page.route("**/api/experiments/17/scan", (r) => r.fulfill(json({ status: "complete", experiment_id: 17, changed: false })));
  await page.route(/\/api\/samples\/\d+\/call$/, (r) => {
    const id = Number(/samples\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill(json({ sample_id: id, state: "indexed", phase: null, call_phases: [], members: [] }));
  });
  await page.route(/\/api\/samples\/\d+\/exposures/, (r) => {
    const id = Number(/samples\/(\d+)\//.exec(r.request().url())![1]);
    const n = FRAMES[id] ?? 1;
    return r.fulfill(json(Array.from({ length: n }, (_, f) => exposureOf(id, id % 2 === 0, f, n > 1 && f === n - 1))));
  });
  await page.route(/\/api\/exposures\/\d+\/trace/, (r) => r.fulfill(json(trace)));
  await page.route(/\/api\/exposures\/\d+\/peaks/, (r) => r.fulfill(json([])));
  await page.route(/\/api\/exposures\/\d+\/indices/, (r) => {
    const eid = Number(/exposures\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill(json(indicesOf(eid, CANDIDATES[Math.floor(eid / 10)] ?? 0)));
  });
  await page.route(/\/api\/exposures\/\d+\/assignment/, (r) => {
    const eid = Number(/exposures\/(\d+)\//.exec(r.request().url())![1]);
    const members = eid === 410 ? indicesOf(eid, 3).map((ix) => ix.id) : [];
    return r.fulfill(json({ exposure_id: eid, state: "indexed", members }));
  });
  await page.route("**/api/events", (r) => r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
});

async function layout(page: Page): Promise<Record<string, number[]>> {
  return page.evaluate(() => {
    const box = (el: Element | null | undefined) => {
      const r = el?.getBoundingClientRect();
      return r ? [Math.round(r.top), Math.round(r.height)] : [];
    };
    const sampleTitle = document.querySelector("main h2");
    const section = (title: string) => [...document.querySelectorAll("aside section")].find((s) => s.querySelector("h2")?.textContent === title);
    return {
      title: box(sampleTitle),
      tags: box(document.querySelector("[data-testid=workspace-tags]")),
      compareBar: box(document.querySelector("[data-testid=workspace-compare-bar]")),
      plot: box(document.querySelector("[data-testid=workspace-trace]")),
      exposures: box(document.querySelector("[data-testid=workspace-exposures]")),
      indexing: box(section("Indexing")),
      indexingPlot: box(document.querySelector("[data-testid=workspace-indexing] svg")),
      exposure: box(section("Exposure")),
    };
  });
}

async function measured(page: Page): Promise<Record<string, number[]>> {
  const boxes = await layout(page);
  for (const [name, b] of Object.entries(boxes)) expect(b.length, `${name} is on the page`).toBe(2);
  return boxes;
}

async function openSample(page: Page, id: number): Promise<void> {
  await page.locator(`[data-sample-id="${id}"]`).click();
  await expect(page.locator("[data-testid=workspace-samples] [aria-current='true']")).toHaveAttribute("data-sample-id", String(id));
  await expect(page.locator("[data-testid=workspace-trace]")).toBeVisible();
  await expect(page.locator("[data-testid=workspace-indexing] svg")).toBeVisible();
  await expect(page.locator("[data-testid=workspace-candidates] li")).toHaveCount(CANDIDATES[id] || 0);
  if (!CANDIDATES[id]) await expect(page.getByText("No candidate indices for this exposure.")).toBeVisible();
  await expect(page.locator("[data-testid=workspace-exposures] > li")).toHaveCount(FRAMES[id] ?? 1);
  await expect(page.locator("aside").getByText(`JC_S${id}_0_00${FRAMES[id] ?? 1}`)).toBeVisible();
}

for (const viewport of VIEWPORTS) {
  test(`the plot, the indexing plot and the exposure details stay put from sample to sample at ${viewport.width}×${viewport.height}`, async ({ page }) => {
    await page.setViewportSize(viewport);
    await page.goto("/experiments/17/workspace/40");
    await openSample(page, 40);
    const above = async () => {
      const { exposures, ...rest } = await measured(page);
      return { ...rest, exposuresTop: exposures[0] };
    };
    const first = await above();
    for (const id of [41, 42, 43, 44, 40]) {
      await openSample(page, id);
      expect(await above(), `layout on sample ${id}`).toEqual(first);
    }
  });

  test(`opening a custom index keeps the indexing plot and the exposure details in place at ${viewport.width}×${viewport.height}`, async ({ page }) => {
    await page.setViewportSize(viewport);
    await page.goto("/experiments/17/workspace/41");
    await openSample(page, 41);
    const before = await measured(page);
    await page.getByRole("button", { name: /Custom index/ }).click();
    await expect(page.getByTestId("workspace-custom-index")).toBeVisible();
    const after = await measured(page);
    expect({ indexing: after.indexing, exposure: after.exposure }).toEqual({ indexing: before.indexing, exposure: before.exposure });
  });
}

test("stepping between frames and hovering frames in another row keeps everything in place", async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.goto("/experiments/17/workspace/44");
  await openSample(page, 44);
  const first = await measured(page);
  const cards = page.locator("[data-testid=workspace-exposures] > li");
  await expect(page.locator("aside").getByText("Averaged")).toBeVisible();
  await cards.nth(0).hover();
  expect(await measured(page), "hovering a frame in the other row").toEqual(first);
  await page.mouse.move(0, 0);
  await page.keyboard.press("ArrowLeft");
  await expect(page.locator("aside").getByText("JC_S44_0_006")).toBeVisible();
  await expect(page.locator("aside").getByText("Single")).toBeVisible();
  await expect(page.locator("[data-testid=workspace-trace]")).toBeVisible();
  expect(await measured(page), "on the single frame").toEqual(first);
  await cards.nth(6).hover();
  expect(await measured(page), "hovering the averaged frame alone in its row").toEqual(first);
});

test("hovering a candidate does not move the candidate list, even when three phases are ticked", async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.goto("/experiments/17/workspace/41");
  await openSample(page, 41);
  const rows = page.locator("[data-testid=workspace-candidates] button");
  await expect(page.locator("[data-testid=workspace-indexing-fit]")).toHaveCount(3);
  const before = await rows.nth(3).boundingBox();
  await rows.nth(0).hover();
  await expect(page.locator("[data-testid=workspace-indexing-fit]")).toHaveCount(1);
  expect(await rows.nth(3).boundingBox()).toEqual(before);
});

test("three comparisons with long names keep the compare bar on one line", async ({ page }) => {
  await page.setViewportSize({ width: 1280, height: 900 });
  await page.goto("/experiments/17/workspace/40");
  await openSample(page, 40);
  const before = await measured(page);
  for (const name of [`${LONG}01`, `${LONG}02`, `${LONG}03`]) {
    await page.keyboard.press("c");
    await page.keyboard.type(name);
    await page.keyboard.press("Enter");
  }
  await expect(page.locator("[data-compare-chip]")).toHaveCount(3);
  const after = await measured(page);
  expect(after.compareBar).toEqual(before.compareBar);
  expect(after.plot[0]).toEqual(before.plot[0]);
});
