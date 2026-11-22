import { test, expect, type Page } from "./fixtures";

const EXPERIMENT_ID = 7;

function makeExperiment(ingest_status = "complete"): object {
  return {
    id: EXPERIMENT_ID,
    name: "SSRL Beamtime",
    path: "/data",
    data_dir: "/data",
    analysis_dir: "/analysis",
    manifest_path: null,
    created_at: "2026-01-01",
    last_scanned_at: "2026-01-01T10:00:00",
    ingest_status,
    stats: { loads: 1, samples: 3, exposures: 3, sessions: 1 },
  };
}

function makeSamples(n: number): object[] {
  return Array.from({ length: n }, (_, i) => ({
    id: 100 + i,
    experiment_id: EXPERIMENT_ID,
    name: `Sample ${i + 1}`,
    notes: null,
  }));
}

function makeLoads(sampleCount: number): object[] {
  const samples = Array.from({ length: sampleCount }, (_, i) => ({
    sample_id: 100 + i,
    name: `Sample ${i + 1}`,
    slot_index: i + 1,
    grouping_source: "auto_position",
    name_source: "auto",
    merged_into_id: null,
    exposures: [
      {
        id: 200 + i,
        filename: `sample_${i + 1}.dat`,
        horizontal_position: 10.0 + i,
        timestamp: "2026-01-01T10:00:00",
      },
      ...(i === 0
        ? [{ id: 300, filename: "sample_1_flare.dat", horizontal_position: 10.0, timestamp: "2026-01-01T10:01:00" }]
        : []),
    ],
  }));
  return [
    {
      load_id: 1,
      load_index: 1,
      session_id: null,
      start_time: null,
      end_time: null,
      frame_count: sampleCount,
      note: null,
      samples,
    },
  ];
}

async function seedState(page: Page): Promise<void> {
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state: { username: "alice", tutorialSeen: true, theme: "dark" }, version: 5 }),
    );
  });
}

async function mockExperiment(
  page: Page,
  opts: { id: number; ingest_status: string; samples: number; problems?: object },
): Promise<void> {
  const exp = makeExperiment(opts.ingest_status);
  const samples = makeSamples(opts.samples);
  const loads = makeLoads(opts.samples);
  const problems = opts.problems ?? { ingest_error: null, exposures: [] };

  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([exp]) }));
  await page.route(`**/api/experiments/${opts.id}`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(exp) }));
  await page.route(`**/api/experiments/${opts.id}/samples`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(samples) }));
  await page.route(`**/api/experiments/${opts.id}/problems`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(problems) }));
  await page.route((url) => url.pathname === "/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(samples) }));
  for (const sample of loads.flatMap((l) => (l as { samples: { sample_id: number; exposures: { id: number; filename: string }[] }[] }).samples)) {
    const exposures = sample.exposures.map((e) => ({
      id: e.id, sample_id: sample.sample_id, filename: e.filename, kind: "file",
      selected: false, status: e.id === 300 ? "rejected" : null, image_path: null,
    }));
    await page.route(`**/api/samples/${sample.sample_id}/exposures`, (r) =>
      r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(exposures) }));
  }
  await page.route("**/api/events*", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
  await page.route(`**/api/experiments/${opts.id}/bridges`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/label-colors", (r) => r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route(`**/api/experiments/${opts.id}/scan`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify({ status: "complete", experiment_id: opts.id, changed: false }) }));
  await page.route(`**/api/experiments/${opts.id}/labels`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  const q = Array.from({ length: 50 }, (_, i) => 0.01 + i * 0.006);
  await page.route(/\/api\/exposures\/(\d+)\/trace$/, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify({ q, I: q.map((x) => 1 / x), sigma: q.map(() => 0.1) }) }));
  await page.route(/\/api\/exposures\/(\d+)\/image/, (r) =>
    r.fulfill({ status: 200, contentType: "image/png", body: Buffer.from("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR4nGNgYGD4DwABBAEAwS2OUAAAAABJRU5ErkJggg==", "base64") }));
  await page.route(/\/api\/exposures\/(\d+)\/(peaks|indices)$/, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route(/\/api\/samples\/(\d+)\/call$/, (r) => {
    const id = Number(/samples\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify({ sample_id: id, state: "null", phase: null, call_phases: [], members: [] }) });
  });
  await page.route(/\/api\/exposures\/(\d+)\/assignment$/, (r) => {
    const id = Number(/exposures\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify({ exposure_id: id, state: "null", members: [] }) });
  });
}

test.beforeEach(async ({ page }) => {
  await seedState(page);
});

test("opening a scanned experiment lands on its first sample in the workspace", async ({ page }) => {
  await mockExperiment(page, { id: 7, ingest_status: "complete", samples: 3 });
  await page.goto("/experiments/7");
  await expect(page).toHaveURL(/\/experiments\/7\/workspace/);
  await expect(page.getByRole("heading", { level: 2, name: "Sample 1" })).toBeVisible();
  await expect(page.getByTestId("workspace-notice")).toHaveCount(0);
});

test("an old corpus link opens the workspace", async ({ page }) => {
  await mockExperiment(page, { id: 7, ingest_status: "complete", samples: 2 });
  await page.goto("/experiments/7/corpus");
  await expect(page).toHaveURL(/\/experiments\/7\/workspace/);
});

test("problem exposures show in the sample pane and open their sample", async ({ page }) => {
  await mockExperiment(page, {
    id: 7, ingest_status: "complete", samples: 3,
    problems: {
      ingest_error: null,
      exposures: [
        { exposure_id: 201, filename: "sample_2", sample_id: 101, sample_name: "Sample 2", problem: "no_trace", reason: "no trace" },
        { exposure_id: 202, filename: "sample_3", sample_id: 102, sample_name: "Sample 3", problem: "failed", reason: "could not read the trace header" },
      ],
    },
  });
  await page.goto("/experiments/7");
  const notice = page.getByTestId("workspace-notice");
  await expect(notice).toContainText("1 frame waiting for a trace · 1 analysis failed");
  await notice.getByRole("button", { name: "Details" }).click();
  const problems = page.getByTestId("workspace-problems");
  await expect(problems).toContainText("could not read the trace header");
  await problems.getByRole("button", { name: /sample_3/ }).click();
  await expect(page).toHaveURL(/\/experiments\/7\/workspace\/102/);
  await expect(page.getByRole("heading", { level: 2, name: "Sample 3" })).toBeVisible();
});

test("a failed ingest opens the workspace and says why", async ({ page }) => {
  await mockExperiment(page, {
    id: 7, ingest_status: "failed", samples: 0,
    problems: { ingest_error: "SystemError: opening the data folder", exposures: [] },
  });
  await page.goto("/experiments/7");
  await expect(page).toHaveURL(/\/experiments\/7\/workspace/);
  const notice = page.getByTestId("workspace-notice");
  await expect(notice).toContainText("The last ingest failed.");
  await notice.getByRole("button", { name: "Details" }).click();
  await expect(page.getByTestId("problems-ingest-error")).toContainText("SystemError: opening the data folder");
});

test("a live rescan SSE frame shows its progress in the workspace", async ({ page }) => {
  await mockExperiment(page, { id: 7, ingest_status: "scanning", samples: 3 });
  await page.route("**/api/events*", (r) =>
    r.fulfill({
      status: 200,
      contentType: "text/event-stream",
      body:
        'event: curation\n' +
        'data: {"kind":"ingest_progress","payload":{"experiment_id":7,"processed":5,"total":10,"stage":"analyzing","phase":"rescan"}}\n\n',
    }));
  await page.goto("/experiments/7");
  await expect(page.getByTestId("workspace-ingest-status")).toHaveText("Analyzing 5 of 10");
});
