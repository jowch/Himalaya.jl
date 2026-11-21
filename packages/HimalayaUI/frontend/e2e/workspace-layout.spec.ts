import { test, expect, type Page } from "./fixtures";

const LONG = "a_sample_with_a_long_descriptive_name_S24";
const SAMPLES = [
  { id: 40, name: "JC001" },
  { id: 41, name: "JC002" },
  { id: 42, name: `${LONG}01` },
  { id: 43, name: `${LONG}02_${LONG}_${LONG}` },
  { id: 44, name: `${LONG}03` },
].map((s) => ({ ...s, experiment_id: 17, display_name: s.name, notes: null, q_units: "A-1" }));

const q = Array.from({ length: 200 }, (_, i) => 0.02 + i * 0.0015);
const trace = { q, I: q.map((x) => 1 / (x * x) + 50 * Math.exp(-(((x - 0.1) / 0.002) ** 2))), sigma: q.map(() => 1) };

function exposureOf(sampleId: number, monitors: boolean, frame = 0, averaged = false) {
  return {
    id: sampleId * 10 + frame, sample_id: sampleId, filename: `JC_S${sampleId}_0_00${frame + 1}`,
    kind: averaged ? "averaged" : "file", selected: true, status: null,
    image_path: null, image_version: "", sources: averaged ? [sampleId * 10] : [], trace_hash: null, analysis_inputs_hash: null,
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
const NAMES = ["buffer", "37C", "heated twice", "retake", "dilute 1 to 10", "cloudy", "after spin", "second capillary", "check flare", "low counts"];
const LABELS = [
  ...NAMES.slice(0, 2).map((label, i) => ({ id: i + 1, sample_id: 42, label, state: "confirmed", source: "user" })),
  ...NAMES.map((label, i) => ({ id: i + 11, sample_id: 44, label, state: i === 0 ? "suggested" : "confirmed", source: "user" })),
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
  await page.route("**/api/label-colors", (r) => r.fulfill(json([])));
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
    if (r.request().method() === "POST") {
      const { index_id } = r.request().postDataJSON() as { index_id: number };
      return r.fulfill(json({ exposure_id: eid, state: "indexed", members: [index_id], event_id: 1, view_row_id: null }));
    }
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
    const section = (title: string) => document.querySelector(`aside section[data-section="${title}"]`);
    return {
      title: box(sampleTitle),
      labels: box(document.querySelector("[data-testid=workspace-labels]")),
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

function exposureToggle(page: Page) {
  return page.locator("aside").getByRole("button", { name: "Exposure" });
}

function callToggle(page: Page) {
  return page.locator("aside").getByRole("button", { name: "Call" });
}

async function exposureOpen(page: Page): Promise<boolean> {
  const toggle = exposureToggle(page);
  return (await toggle.count()) === 0 || (await toggle.getAttribute("aria-expanded")) === "true";
}

async function openSample(page: Page, id: number): Promise<void> {
  await page.locator(`[data-sample-id="${id}"]`).click();
  await expect(page.locator("[data-testid=workspace-samples] [aria-current='true']")).toHaveAttribute("data-sample-id", String(id));
  await expect(page.locator("[data-testid=workspace-trace]")).toBeVisible();
  await expect(page.locator("[data-testid=workspace-indexing] svg")).toBeVisible();
  if ((await callToggle(page).count()) === 0 || (await callToggle(page).getAttribute("aria-expanded")) === "true") {
    await expect(page.locator("[data-testid=workspace-candidates] li")).toHaveCount(CANDIDATES[id] || 0);
    if (!CANDIDATES[id]) await expect(page.getByText("No candidate indices for this exposure.")).toBeVisible();
  }
  await expect(page.locator("[data-testid=workspace-exposures] > li")).toHaveCount(FRAMES[id] ?? 1);
  const file = `JC_S${id}_0_00${FRAMES[id] ?? 1}`;
  if (await exposureOpen(page)) await expect(page.locator("aside").getByText(file)).toBeVisible();
  else await expect(page.locator("[data-testid=workspace-exposures] [aria-current='true']")).toBeVisible();
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

test("on a short window Call and Exposure fold as an accordion, and stay put from sample to sample", async ({ page }) => {
  await page.setViewportSize({ width: 1366, height: 650 });
  await page.goto("/experiments/17/workspace/41");
  await openSample(page, 41);
  await expect(callToggle(page)).toHaveAttribute("aria-expanded", "true");
  await expect(exposureToggle(page)).toHaveAttribute("aria-expanded", "false");
  await expect(page.locator("aside").getByText("JC_S41_0_001")).toHaveCount(0);
  const aside = page.locator("aside");
  const fits = () => aside.evaluate((a) => a.scrollHeight <= a.clientHeight);
  expect(await fits(), "the tools column itself never scrolls").toBe(true);

  await exposureToggle(page).click();
  await expect(exposureToggle(page)).toHaveAttribute("aria-expanded", "true");
  await expect(callToggle(page)).toHaveAttribute("aria-expanded", "false");
  await expect(page.locator("[data-testid=workspace-candidates]")).toHaveCount(0);
  await expect(page.locator("aside").getByText("JC_S41_0_001")).toBeVisible();
  expect(await fits()).toBe(true);
  const lastRow = () => page.locator("aside section[data-section=Exposure] dl > :last-child > dd").evaluate((el) => Math.round(el.getBoundingClientRect().top));
  const first = await measured(page);
  const firstLastRow = await lastRow();
  for (const id of [44, 43, 40, 41]) {
    await openSample(page, id);
    await expect(exposureToggle(page)).toHaveAttribute("aria-expanded", "true");
    const { exposures, ...rest } = await measured(page);
    const { exposures: firstExposures, ...firstRest } = first;
    expect({ ...rest, exposuresTop: exposures[0] }, `layout on sample ${id}`).toEqual({ ...firstRest, exposuresTop: firstExposures[0] });
    expect(await lastRow(), `last exposure row on sample ${id}`).toBe(firstLastRow);
  }

  await page.keyboard.press("i");
  await expect(page.getByTestId("workspace-custom-index")).toBeVisible();
  await expect(exposureToggle(page), "Exposure is not a dead toggle while a custom index is open").toHaveCount(0);
  await expect(page.locator("aside section[data-section=Exposure]")).toHaveAttribute("data-collapsed", "true");
  await page.keyboard.press("Escape");
  await expect(exposureToggle(page)).toHaveAttribute("aria-expanded", "true");
});

test("on a short window a number key opens the Call list before making the call", async ({ page }) => {
  await page.setViewportSize({ width: 1366, height: 650 });
  await page.goto("/experiments/17/workspace/42");
  await openSample(page, 42);
  const first = page.locator("[data-testid=workspace-candidates] button[aria-pressed]").first();
  await expect(first).toHaveAttribute("aria-pressed", "false");
  await exposureToggle(page).click();
  await expect(callToggle(page)).toHaveAttribute("aria-expanded", "false");
  await page.keyboard.press("1");
  await expect(callToggle(page)).toHaveAttribute("aria-expanded", "true");
  await expect(first).toHaveAttribute("aria-pressed", "true");
});

test("on a medium window the Call list scrolls inside a column that never scrolls", async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 880 });
  await page.goto("/experiments/17/workspace/41");
  await openSample(page, 41);
  await expect(exposureToggle(page)).toHaveCount(0);
  const state = await page.evaluate(() => {
    const aside = document.querySelector("aside")!;
    const middle = document.querySelector("[data-testid=workspace-tools-middle]") as HTMLElement;
    return { asideFits: aside.scrollHeight <= aside.clientHeight, middleScrolls: middle.scrollHeight > middle.clientHeight };
  });
  expect(state).toEqual({ asideFits: true, middleScrolls: true });
});

test("on a very short window the open Exposure section scrolls inside a column that never scrolls", async ({ page }) => {
  await page.setViewportSize({ width: 1366, height: 520 });
  await page.goto("/experiments/17/workspace/44");
  await openSample(page, 44);
  await exposureToggle(page).click();
  await expect(page.locator("aside").getByText("JC_S44_0_007")).toBeVisible();
  const state = await page.evaluate(() => {
    const aside = document.querySelector("aside")!;
    const exposure = aside.querySelector("section[data-section=Exposure]") as HTMLElement;
    return { asideFits: aside.scrollHeight <= aside.clientHeight, exposureScrolls: exposure.scrollHeight > exposure.clientHeight };
  });
  expect(state).toEqual({ asideFits: true, exposureScrolls: true });
});

for (const viewport of [...VIEWPORTS, { width: 1024, height: 700 }]) {
  test(`every frame's Drop or Restore button stays inside its frame card at ${viewport.width}×${viewport.height}`, async ({ page }) => {
    await page.setViewportSize(viewport);
    await page.route(/\/api\/samples\/44\/exposures/, (r) => r.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify(Array.from({ length: 7 }, (_, f) => ({ ...exposureOf(44, true, f, f === 6), status: f % 2 === 0 ? "rejected" : null }))),
    }));
    await page.goto("/experiments/17/workspace/44");
    const frames = page.locator("[data-testid=workspace-exposures] > li");
    await expect(frames).toHaveCount(7);
    for (let i = 0; i < 7; i++) {
      await frames.nth(i).hover();
      const button = frames.nth(i).getByRole("button", { name: /^(Drop|Restore)/ });
      await expect(button).toBeVisible();
      const [card, b] = await Promise.all([frames.nth(i).boundingBox(), button.boundingBox()]);
      expect(b!.x + b!.width, `frame ${i + 1} button right edge`).toBeLessThanOrEqual(card!.x + card!.width + 0.5);
      expect(b!.x, `frame ${i + 1} button left edge`).toBeGreaterThanOrEqual(card!.x - 0.5);
    }
  });
}
