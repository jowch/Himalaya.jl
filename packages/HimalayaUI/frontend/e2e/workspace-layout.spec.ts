import { test, expect, type Page } from "./fixtures";

const LONG = "a_sample_with_a_long_descriptive_name_S24";
const SAMPLES = [
  { id: 40, name: "JC001" },
  { id: 41, name: "JC002" },
  { id: 42, name: `${LONG}01` },
  { id: 43, name: `${LONG}02_${LONG}_${LONG}` },
  { id: 44, name: `${LONG}03` },
  { id: 45, name: "JC045_many_peaks" },
].map((s) => ({ ...s, experiment_id: 17, display_name: s.name, notes: null, q_units: "A-1" }));

const q = Array.from({ length: 200 }, (_, i) => 0.02 + i * 0.0015);
const trace = { q, I: q.map((x) => 1 / (x * x) + 50 * Math.exp(-(((x - 0.1) / 0.002) ** 2))), sigma: q.map(() => 1) };

function exposureOf(sampleId: number, monitors: boolean, frame = 0, averaged = false) {
  return {
    id: sampleId * 10 + frame, sample_id: sampleId, filename: `JC_S${sampleId}_0_00${frame + 1}`,
    kind: averaged ? "averaged" : "file", status: null,
    image_path: null, image_version: "", sources: averaged ? [sampleId * 10] : [], trace_hash: null, analysis_inputs_hash: null,
    exposure_time: monitors ? 15 : null, counting_time: monitors ? 15 : null,
    transmitted_counts: monitors ? 7000 : null, incident_monitor: monitors ? 1.5 : null,
    timestamp: monitors ? "2026-04-26T18:20:27" : null, beam_energy_ev: monitors ? 11000 : null,
    pipe_length_m: monitors ? 1.7 : null, detector: monitors ? "Pilatus 1M" : null,
  };
}

const PHASES = ["Pn3m", "Im3m", "Ia3d", "HII", "Lamellar", "Fd3m"];
function indicesOf(exposureId: number, n: number) {
  return Array.from({ length: n }, (_, i) => PHASES[i % PHASES.length]!).map((phase, i) => ({
    id: exposureId * 100 + i, exposure_id: exposureId, phase, basis: 0.1, score: 0.9 - i * 0.01, r_squared: 0.99,
    lattice_d: 100 + i, status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [0.1, 0.141],
  }));
}

const CANDIDATES: Record<number, number> = { 40: 1, 41: 6, 42: 3, 43: 0, 44: 2, 45: 10 };
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
    localStorage.setItem("himalaya-ui:state", JSON.stringify({ state: { username: "alice", tutorialSeen: true }, version: 6 }));
  });
  const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });
  const experiment = { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "", manifest_path: null, created_at: "", q_units: null };
  await page.route("**/api/experiments", (r) => r.fulfill(json([experiment])));
  await page.route("**/api/experiments/17", (r) => r.fulfill(json(experiment)));
  await page.route((url) => url.pathname === "/api/samples", (r) => r.fulfill(json(SAMPLES)));
  await page.route("**/api/users", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/bridges", (r) => r.fulfill(json([])));
  const compares = { compares: [] as { sample_id: number; ref_sample_id: number; added_by: string | null }[], subtractions: [] as { sample_id: number; ref_sample_id: number }[] };
  await page.route("**/api/experiments/17/compares", (r) => r.fulfill(json(compares)));
  await page.route(/\/api\/samples\/\d+\/compares$/, (r) => {
    const sampleId = Number(new URL(r.request().url()).pathname.split("/")[3]);
    const change = r.request().postDataJSON() as { add?: number[]; remove?: number[]; subtraction?: number | null };
    for (const ref of change.add ?? []) compares.compares.push({ sample_id: sampleId, ref_sample_id: ref, added_by: "alice" });
    compares.compares = compares.compares.filter((c) => !(c.sample_id === sampleId && (change.remove ?? []).includes(c.ref_sample_id)));
    if ("subtraction" in change) {
      compares.subtractions = compares.subtractions.filter((t) => t.sample_id !== sampleId);
      if (change.subtraction !== null && change.subtraction !== undefined) compares.subtractions.push({ sample_id: sampleId, ref_sample_id: change.subtraction });
    }
    return r.fulfill(json(compares));
  });
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
    return r.fulfill(json(Array.from({ length: n }, (_, f) => ({ ...exposureOf(id, id % 2 === 0, f, n > 1 && f === n - 1), selected: f === n - 1 }))));
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
      exposure: box(section("Image")),
    };
  });
}

async function measured(page: Page): Promise<Record<string, number[]>> {
  const boxes = await layout(page);
  for (const [name, b] of Object.entries(boxes)) expect(b.length, `${name} is on the page`).toBe(2);
  return boxes;
}

function exposureToggle(page: Page) {
  return page.locator("aside").getByRole("button", { name: "Image" });
}

function callToggle(page: Page) {
  return page.locator("aside").getByRole("button", { name: "Assignment" });
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
    if (!CANDIDATES[id]) await expect(page.getByText("No candidate indices for this image.")).toBeVisible();
  }
  await expect(page.locator("[data-testid=workspace-exposures] li")).toHaveCount(FRAMES[id] ?? 1);
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
    for (const id of [41, 42, 45, 43, 44, 40]) {
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
  const cards = page.locator("[data-testid=workspace-exposures] li");
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

test("a Custom cubic row keeps its curvature and unit whole on one line", async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.route(/\/api\/exposures\/410\/indices/, (r) =>
    r.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify([["Im3m", 699.9], ["Ia3d", 50], ["Pn3m", 499.9], ["Hexagonal", 499.9]].map(([phase, lattice_d], i) => ({
        id: 41000 + i, exposure_id: 410, phase, basis: 0.1, score: null, r_squared: 0.99, lattice_d, 
        status: "candidate", kind: "speculative", inputs_hash: null, peaks: [], predicted_q: [0.1, 0.141],
      }))),
    }),
  );
  await page.goto("/experiments/17/workspace/41");
  const rows = page.locator("[data-testid=workspace-candidates] button");
  await expect(rows).toHaveCount(4);
  const plain = await rows.filter({ hasText: "Hexagonal" }).boundingBox();
  for (const row of await rows.filter({ hasNotText: "Hexagonal" }).all()) {
    await expect(row).toContainText("Custom");
    await expect(row).toContainText("nm⁻²");
    const inside = await row.evaluate((b) => {
      const cell = b.children[2]!.getBoundingClientRect();
      return b.querySelector("[data-testid=workspace-candidate-curvature]")!.getBoundingClientRect().right <= cell.right;
    });
    expect(inside).toBe(true);
    expect((await row.boundingBox())!.height).toBe(plain!.height);
  }
});

test("an indented Custom cubic in a Bonnet group keeps its curvature and unit whole on one line", async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  const peaksAt = (qs: readonly number[]) => qs.map((q_observed, i) => ({ peak_id: 1 + i, ratio_position: 1 + i, residual: 0, q_observed }));
  const group: [string, number, number[]][] = [
    ["Im3m", 699.9, [0.01, 0.02, 0.03]], ["Pn3m", 699.9 / 1.2793, [0.04, 0.05, 0.06]], ["Ia3d", 699.9 * 1.2319, [0.07, 0.08, 0.09]],
  ];
  await page.route(/\/api\/exposures\/410\/indices/, (r) =>
    r.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify([...group, ["Hexagonal", 499.9, [0.11, 0.12, 0.13]] as const].map(([phase, lattice_d, qs], i) => ({
        id: 41000 + i, exposure_id: 410, phase, basis: 0.1, score: null, r_squared: 0.99, lattice_d,
        status: "candidate", kind: "speculative", inputs_hash: null, peaks: peaksAt(qs), predicted_q: [0.1, 0.141],
      }))),
    }),
  );
  await page.goto("/experiments/17/workspace/41");
  const header = page.locator("[data-testid=workspace-bonnet-group]");
  await expect(header).toContainText("nm⁻²");
  const rows = page.locator("[data-testid=workspace-candidates] li li button");
  await expect(rows).toHaveCount(3);
  const plain = await page.locator("[data-testid=workspace-candidates] button").filter({ hasText: "Hexagonal" }).boundingBox();
  expect((await header.boundingBox())!.height).toBe(plain!.height);
  const headerInside = await header.evaluate((b) =>
    b.querySelector("[data-testid=workspace-bonnet-curvature]")!.getBoundingClientRect().right <= b.children[1]!.getBoundingClientRect().right);
  expect(headerInside).toBe(true);
  for (const row of await rows.all()) {
    await expect(row).toContainText("Custom");
    await expect(row).toContainText("nm⁻²");
    const inside = await row.evaluate((b) => {
      const cell = b.children[2]!.getBoundingClientRect();
      return b.querySelector("[data-testid=workspace-candidate-curvature]")!.getBoundingClientRect().right <= cell.right;
    });
    expect(inside).toBe(true);
    expect((await row.boundingBox())!.height).toBe(plain!.height);
  }
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

for (const viewport of VIEWPORTS) {
  const short = viewport.height <= 860;
  const order = short ? "Assignment, then Image, with Indexing pinned at the bottom" : "Assignment, then Indexing, then Image";
  test(`the tools column reads ${order} at ${viewport.width}×${viewport.height}`, async ({ page }) => {
    await page.setViewportSize(viewport);
    await page.goto("/experiments/17/workspace/41");
    await openSample(page, 41);
    const tops = await page.evaluate(() => {
      const top = (sel: string) => Math.round(document.querySelector(sel)!.getBoundingClientRect().top);
      return { call: top("aside section[data-section=Assignment]"), indexing: top("aside section[data-section=Indexing]"), exposure: top("aside section[data-section=Image]") };
    });
    const sequence = Object.entries(tops).sort((a, b) => a[1] - b[1]).map(([name]) => name);
    expect(sequence).toEqual(short ? ["call", "exposure", "indexing"] : ["call", "indexing", "exposure"]);
  });
}

test("a long Call list (ten candidates, more than the number keys reach) keeps Indexing and Frame in place, and the column never scrolls", async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.goto("/experiments/17/workspace/45");
  await openSample(page, 45);
  const before = await measured(page);
  const middle = page.getByTestId("workspace-tools-middle");
  expect(await middle.evaluate((m) => m.scrollHeight > m.clientHeight), "the full list overflows its slot").toBe(true);
  const last = page.locator("[data-testid=workspace-candidates] li").last();
  await last.scrollIntoViewIfNeeded();
  await expect(last).toBeInViewport();
  expect(await middle.evaluate((m) => m.scrollTop), "the list scrolled, not the column").toBeGreaterThan(0);
  const after = await measured(page);
  expect({ indexing: after.indexing, exposure: after.exposure }).toEqual({ indexing: before.indexing, exposure: before.exposure });
  expect(await page.locator("aside").evaluate((a) => a.scrollHeight <= a.clientHeight), "the tools column itself never scrolls").toBe(true);
});

test("subtracting a comparison, switching between Raw and Minus and stopping move nothing", async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.goto("/experiments/17/workspace/40");
  await openSample(page, 40);
  await page.keyboard.press("c");
  await page.keyboard.type(`${LONG}01`);
  await page.keyboard.press("Enter");
  const toggle = page.getByTestId("workspace-subtraction-toggle");
  const chip = page.locator("[data-compare-chip] button[aria-pressed]");
  await expect(chip).toBeVisible();
  await expect(toggle).toHaveCount(0);
  const before = await measured(page);
  await chip.click({ modifiers: ["ControlOrMeta"] });
  await expect(toggle.getByRole("button", { name: /^Minus / })).toHaveAttribute("aria-pressed", "true");
  await expect(page.locator("[data-role=difference-line]")).toHaveCount(1);
  expect(await measured(page)).toEqual(before);
  await toggle.getByRole("button", { name: "Raw" }).click();
  await expect(page.locator("[data-role=trace-line]")).toHaveCount(1);
  expect(await measured(page)).toEqual(before);
  await chip.click({ modifiers: ["ControlOrMeta"] });
  await expect(toggle).toHaveCount(0);
  expect(await measured(page)).toEqual(before);
});

test("on a short window Assignment and Frame fold as an accordion, and stay put from sample to sample", async ({ page }) => {
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
  const lastRow = () => page.locator("aside section[data-section=Image] dl > :last-child > dd").evaluate((el) => Math.round(el.getBoundingClientRect().top));
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
  await expect(page.locator("aside section[data-section=Image]")).toHaveAttribute("data-collapsed", "true");
  await page.keyboard.press("Escape");
  await expect(exposureToggle(page)).toHaveAttribute("aria-expanded", "true");
});

test("on a short window, folding Assignment or Frame and opening a custom index never moves the indexing chart", async ({ page }) => {
  await page.setViewportSize({ width: 1366, height: 650 });
  await page.goto("/experiments/17/workspace/41");
  await openSample(page, 41);
  const indexing = async () => (await measured(page)).indexing;
  const start = await indexing();
  await exposureToggle(page).click();
  await expect(exposureToggle(page)).toHaveAttribute("aria-expanded", "true");
  expect(await indexing(), "after opening Image").toEqual(start);
  await page.getByRole("button", { name: /Custom index/ }).click();
  await expect(page.getByTestId("workspace-custom-index")).toBeVisible();
  expect(await indexing(), "after opening a custom index from Image").toEqual(start);
  await page.getByRole("button", { name: /Custom index/ }).click();
  await expect(page.getByTestId("workspace-custom-index")).toHaveCount(0);
  expect(await indexing(), "after closing it with the same button").toEqual(start);
  await callToggle(page).click();
  expect(await indexing(), "after opening Assignment again").toEqual(start);
});

test("the indexing chart and the exposure details hold their place while a sample's exposures load", async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.goto("/experiments/17/workspace/40");
  await openSample(page, 40);
  const loaded = await measured(page);
  let release: () => void = () => {};
  const gate = new Promise<void>((r) => { release = r; });
  await page.route(/\/api\/samples\/41\/exposures/, async (r) => { await gate; await r.fallback(); });
  await page.locator('[data-sample-id="41"]').click();
  await expect(page.locator("[data-testid=workspace-samples] [aria-current='true']")).toHaveAttribute("data-sample-id", "41");
  const loading = await layout(page);
  release();
  expect({ indexing: loading.indexing, exposure: loading.exposure }).toEqual({ indexing: loaded.indexing, exposure: loaded.exposure });
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

test("on a very short window the open Frame section scrolls inside a column that never scrolls", async ({ page }) => {
  await page.setViewportSize({ width: 1366, height: 520 });
  await page.goto("/experiments/17/workspace/44");
  await openSample(page, 44);
  await exposureToggle(page).click();
  await expect(page.locator("aside").getByText("JC_S44_0_007")).toBeVisible();
  const state = await page.evaluate(() => {
    const aside = document.querySelector("aside")!;
    const exposure = aside.querySelector("section[data-section=Image]") as HTMLElement;
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
    const frames = page.locator("[data-testid=workspace-exposures] li");
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

for (const width of [1024, 1280]) {
  for (const name of ["DOPC_DOPE_1to1_S2417 (S02P02)", "buffer_blank_S2424 (S02P05)"]) {
    for (const id of [42, 43]) {
      test(`at ${width} px "${name}" stays whole on the title row ${id === 42 ? "with two labels" : "with no labels"}`, async ({ page }) => {
        const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });
        await page.setViewportSize({ width, height: 760 });
        await page.route((url) => url.pathname === "/api/samples", (r) =>
          r.fulfill(json([{ id, experiment_id: 17, display_name: name, name, notes: null, q_units: "A-1" }])));
        await page.goto(`/experiments/17/workspace/${id}`);
        const title = page.getByRole("heading", { level: 2, name });
        await expect(page.locator("[data-testid=workspace-exposures] li").first()).toBeVisible();
        const meta = page.getByText("Image 1 of 1", { exact: true });
        if (width === 1280) await expect(meta).toBeVisible();
        expect(await title.evaluate((el) => el.scrollWidth <= el.clientWidth)).toBe(true);
        const add = page.getByRole("button", { name: "+ Label L" });
        await expect(add).toBeVisible();
        const heading = (await title.evaluate((el) => el.parentElement!.parentElement!.getBoundingClientRect().toJSON())) as DOMRect;
        const btn = (await add.boundingBox())!;
        expect(btn.x + btn.width).toBeLessThanOrEqual(heading.x + heading.width);
        if (id === 42) await expect(page.getByTestId("workspace-label-chips").locator("[data-label], [data-testid=workspace-labels-more]").first()).toBeVisible();
      });
    }
  }
}

const CHIP_PEG = "MO_PEG5_S2406 (S01P03)";

async function openChips(page: Page, width: number, refs: readonly number[], subtracting: boolean): Promise<void> {
  const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });
  await page.setViewportSize({ width, height: 800 });
  const setup = (id: number, name: string) => ({ id, experiment_id: 17, folder: name, display_name: null, name, position: id, included: true, root_dir: "", data_dir: "", analysis_dir: "" });
  const experiment = { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "", manifest_path: null, created_at: "", q_units: null, setups: [setup(1, "setup1"), setup(2, "setup2")] };
  await page.route("**/api/experiments/17", (r) => r.fulfill(json(experiment)));
  const samples = [
    { id: 40, name: "MO_H2O_S2401 (S01P01)", setup_id: 1 },
    { id: 42, name: CHIP_PEG, setup_id: 1 },
    { id: 44, name: CHIP_PEG, setup_id: 2 },
    { id: 46, name: "MO_PEG10_S2409 (S01P04)", setup_id: 1 },
    { id: 48, name: "MO_PEG10_S2409 (S01P04)", setup_id: 2 },
    { id: 50, name: "vesicles_ULV_S2425 (S02P06)", setup_id: 1 },
  ].map((s) => ({ ...s, experiment_id: 17, display_name: s.name, notes: null, q_units: "A-1" }));
  await page.route((url) => url.pathname === "/api/samples", (r) => r.fulfill(json(samples)));
  let subtracted: number | null = subtracting && refs.length > 0 ? refs[0]! : null;
  const state = () => ({
    compares: refs.map((ref) => ({ sample_id: 40, ref_sample_id: ref, added_by: "alice" })),
    subtractions: subtracted === null ? [] : [{ sample_id: 40, ref_sample_id: subtracted }],
  });
  await page.route("**/api/experiments/17/compares", (r) => r.fulfill(json(state())));
  await page.route("**/api/samples/40/compares", (r) => {
    const change = r.request().postDataJSON() as { subtraction?: number | null };
    if (change.subtraction !== undefined) subtracted = change.subtraction;
    return r.fulfill(json(state()));
  });
  await page.goto("/experiments/17/workspace/40");
  await expect(page.locator("[data-compare-chip]")).toHaveCount(refs.length);
  await expect(page.getByTestId("workspace-trace")).toBeVisible();
  if (refs.length > 0) await expect(page.locator("[data-scale-value]")).toHaveCount(refs.length);
}

function barRects(page: Page) {
  return page.getByTestId("workspace-compare-bar").locator("[data-compare-chip], [data-compare-chips] > button").evaluateAll((ns) =>
    ns.map((n) => { const r = n.getBoundingClientRect(); return [Math.round(r.left), Math.round(r.top), Math.round(r.width)]; }));
}

for (const width of [1024, 1280]) {
  for (const refs of [[42], [42, 44, 46]]) {
    test(`at ${width} px with ${refs.length} compare chip${refs.length === 1 ? "" : "s"}, turning subtraction on and off moves no chip`, async ({ page }) => {
      await openChips(page, width, refs, false);
      const before = await barRects(page);
      const first = page.locator(`[data-compare-chip="${refs[0]}"] button[aria-pressed]`);
      await first.click({ modifiers: ["ControlOrMeta"] });
      await expect(page.getByTestId("workspace-subtraction-toggle")).toBeVisible();
      expect(await barRects(page)).toEqual(before);
      await first.click({ modifiers: ["ControlOrMeta"] });
      await expect(page.getByTestId("workspace-subtraction-toggle")).toHaveCount(0);
      expect(await barRects(page)).toEqual(before);
    });
  }
  test(`at ${width} px the plot starts at the same height with no compare chips and with five`, async ({ page }) => {
    await openChips(page, width, [], false);
    const none = await layout(page);
    await openChips(page, width, [42, 44, 46, 48, 50], true);
    const five = await layout(page);
    expect(five.plot![0]).toBe(none.plot![0]);
    expect(five.compareBar).toEqual(none.compareBar);
  });
}

for (const width of [1024, 1280, 1440]) {
  for (const subtracting of [false, true]) {
    test(`at ${width} px five compare chips${subtracting ? " and the subtraction toggle" : ""} fit in the two reserved rows, each whole`, async ({ page }) => {
      await openChips(page, width, [42, 44, 46, 48, 50], subtracting);
      const peg = CHIP_PEG;
      const bar = page.getByTestId("workspace-compare-bar");
      await expect(bar.getByRole("button", { name: peg, exact: true })).toBeVisible();
      await expect(bar.getByRole("button", { name: `${peg}, setup2`, exact: true })).toBeVisible();
      await expect(bar.locator("[data-scale-value]")).toHaveCount(5);
      if (subtracting) await expect(bar.getByTestId("workspace-subtraction-toggle")).toBeVisible();
      expect(await bar.locator("[data-compare-chip] [data-chip-setup]").count()).toBe(2);
      const area = await bar.locator("[data-compare-chips]").evaluate((n) => ({ scroll: n.scrollHeight, client: n.clientHeight, right: n.getBoundingClientRect().right }));
      expect(area.scroll).toBeLessThanOrEqual(area.client);
      const chips = await bar.locator("[data-compare-chip]").evaluateAll((ns) => ns.map((n) => {
        const box = (sel: string) => n.querySelector(sel)!.getBoundingClientRect();
        const r = n.getBoundingClientRect();
        return {
          left: r.left, right: r.right, top: r.top,
          inner: Math.max(...Array.from(n.children).map((c) => c.getBoundingClientRect().right)),
          value: box("[data-scale-value]"), wrapper: box("[data-scale-source]"), remove: box("button[aria-label^='Remove']"),
          name: box("button[aria-pressed] > span:first-child").width,
        };
      }));
      for (const [i, c] of chips.entries()) {
        expect(c.inner).toBeLessThanOrEqual(c.right + 1);
        expect(c.right).toBeLessThanOrEqual(area.right + 1);
        expect(c.wrapper.width + 0.5).toBeGreaterThanOrEqual(c.value.width);
        expect(c.remove.left + 0.5).toBeGreaterThanOrEqual(c.value.right);
        expect(c.name).toBeGreaterThan(width >= 1280 && !subtracting ? 60 : 20);
        const prev = chips[i - 1];
        if (prev && prev.top === c.top) expect(prev.right).toBeLessThanOrEqual(c.left);
      }
      if (subtracting) {
        const toggle = await bar.getByTestId("workspace-subtraction-toggle").boundingBox();
        const last = chips[chips.length - 1]!;
        expect(toggle!.x + toggle!.width).toBeLessThanOrEqual(area.right + 1);
        if (Math.abs(toggle!.y - last.top) < 4) expect(last.right).toBeLessThanOrEqual(toggle!.x);
      }
    });
  }
}

for (const width of [1024, 1440]) {
  for (const subtracting of [false, true]) {
  test(`at ${width} px the plot view switcher sits above the ${subtracting ? "difference" : "trace"} and below the peak markers in every view`, async ({ page }) => {
    const peakAt = (id: number, q: number) => ({ id, exposure_id: 1, q, intensity: null, prominence: null, sharpness: null, source: "auto", excluded: false });
    await page.route(/\/api\/exposures\/\d+\/peaks/, (r) => r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([peakAt(1, 0.1), peakAt(2, q[q.length - 2]!)]) }));
    const plain = { q, I: q.map((x) => 1 / (x * x)), sigma: q.map(() => 1) };
    await page.route(/\/api\/exposures\/\d+\/trace/, (r) => {
      const eid = Number(/exposures\/(\d+)\//.exec(r.request().url())![1]);
      return r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(Math.floor(eid / 10) === 40 ? trace : plain) });
    });
    await openChips(page, width, [42, 44, 46, 48, 50], subtracting);
    if (subtracting) {
      await page.getByTestId("workspace-subtraction-toggle").getByRole("button", { name: /^Minus / }).click();
      await expect(page.locator('[data-role="difference-line"] path')).toBeVisible();
    }
    const control = page.getByTestId("workspace-plot-view");
    const before = await control.boundingBox();
    for (const name of ["Log", "Linear", "Log-log", "q²I"]) {
      await control.getByRole("button", { name: `${name === "q²I" ? "q²I (Kratky)" : name} plot`, exact: true }).click();
      const box = (await control.boundingBox())!;
      expect(box).toEqual(before);
      const line = (await page.locator(`[data-role="${subtracting ? "difference-line" : "trace-line"}"] path`).boundingBox())!;
      expect(box.y + box.height).toBeLessThanOrEqual(line.y);
      const marks = await page.locator('[data-role="peak-glyph"] polygon').evaluateAll((ns) => ns.map((n) => n.getBoundingClientRect().bottom));
      expect(marks.length).toBe(2);
      for (const bottom of marks) expect(box.y).toBeGreaterThanOrEqual(bottom);
    }
  });
  }
}

test("the plot view switcher stays faint until the pointer is on it", async ({ page }) => {
  await openChips(page, 1440, [], false);
  const corner = page.getByTestId("plot-card-corner");
  const opacity = () => corner.evaluate((n) => Number(getComputedStyle(n).opacity));
  await page.getByTestId("workspace-trace").hover();
  await expect.poll(opacity).toBeLessThan(0.5);
  await page.getByTestId("workspace-plot-view").getByRole("button", { name: "Linear plot", exact: true }).hover();
  await expect.poll(opacity).toBe(1);
  await page.mouse.move(0, 0);
  await expect.poll(opacity).toBeLessThan(0.5);
});
