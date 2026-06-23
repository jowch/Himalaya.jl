/**
 * E2E: ExperimentCorpusPage — mocked Playwright spec (Task 4.3).
 *
 * Exercises the `/experiments/:id/corpus` route and its §6.2 state machine:
 *   1. `clean` — ingest_status=complete + no grouping flags → SheetTable visible
 *      with one row per sample (header + N rows).
 *   2. `scanning` — inFlight.status=scanning → GroupingReviewPage rendered.
 *   3. `has-flags` — loads with flagged samples → banner + SheetTable.
 *   4. `failed` — ingest_status=failed → ScanFailedPage rendered.
 */

import { test, expect, type Page } from "@playwright/test";

// ---------------------------------------------------------------------------
// Fixture helpers
// ---------------------------------------------------------------------------

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
    // Non-null so a `scanning` row renders the inline rescan ProgressBar
    // (rescanning), not the initial-scan grouping-takeover redirect —
    // last_scanned_at is null only before the first scan completes.
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
    tags: [],
  }));
}

function makeLoads(sampleCount: number, withFlag = false): object[] {
  const samples = Array.from({ length: sampleCount }, (_, i) => ({
    sample_id: 100 + i,
    name: `Sample ${i + 1}`,
    slot_index: i + 1,
    grouping_source: "auto_position",
    name_source: "auto",
    merged_into_id: null,
    flag: withFlag && i === 0
      ? { kind: "merge", merge_with_sample_id: 101, merge_with_label: "Sample 2" }
      : null,
    exposures: [
      {
        id: 200 + i,
        filename: `sample_${i + 1}.dat`,
        horizontal_position: 10.0 + i,
        timestamp: "2026-01-01T10:00:00",
      },
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

/** Seed localStorage so the tutorial and username prompts don't gate navigation. */
async function seedState(page: Page): Promise<void> {
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state: { username: "alice", tutorialSeen: true, theme: "dark" }, version: 5 }),
    );
  });
}

/**
 * Register base API mocks for the experiment corpus route.
 * `mockExperiment(page, { id: 7, ingest_status: 'complete', samples: 3 })`
 */
async function mockExperiment(
  page: Page,
  opts: { id: number; ingest_status: string; samples: number; withFlag?: boolean },
): Promise<void> {
  const exp = makeExperiment(opts.ingest_status);
  const samples = makeSamples(opts.samples);
  const loads = makeLoads(opts.samples, opts.withFlag ?? false);

  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([exp]) }));
  await page.route(`**/api/experiments/${opts.id}`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(exp) }));
  await page.route(`**/api/experiments/${opts.id}/samples`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(samples) }));
  await page.route(`**/api/experiments/${opts.id}/loads`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(loads) }));
  // Corpus-wide sample list (NavModal / SheetTable / SamplesPage may request it).
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(samples) }));
  await page.route("**/api/sample-tags", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/series", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  // SSE endpoint — return an empty stream so the SSE hook doesn't hang.
  await page.route("**/api/events*", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
}

// ---------------------------------------------------------------------------
// Setup
// ---------------------------------------------------------------------------

test.beforeEach(async ({ page }) => {
  await seedState(page);
});

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test("experiment corpus shows the scoped sheet when clean", async ({ page }) => {
  await mockExperiment(page, { id: 7, ingest_status: "complete", samples: 3 });
  await page.goto("/experiments/7/corpus");
  await expect(page.getByRole("table")).toBeVisible();
  await expect(page.getByRole("row")).toHaveCount(4); // header + 3
});

test("experiment corpus shows the scoped sheet when complete with no flags", async ({ page }) => {
  await mockExperiment(page, { id: 7, ingest_status: "complete", samples: 2 });
  await page.goto("/experiments/7/corpus");
  await expect(page.getByRole("table")).toBeVisible();
  await expect(page.getByRole("row")).toHaveCount(3); // header + 2
  // No banner when there are no grouping flags.
  await expect(page.getByTestId("grouping-review-banner")).toHaveCount(0);
});

test("experiment corpus shows review banner when loads have flagged samples", async ({ page }) => {
  await mockExperiment(page, { id: 7, ingest_status: "complete", samples: 3, withFlag: true });
  await page.goto("/experiments/7/corpus");
  // Sheet still visible.
  await expect(page.getByRole("table")).toBeVisible();
  // Banner appears above the sheet.
  await expect(page.getByTestId("grouping-review-banner")).toBeVisible();
  await expect(page.getByTestId("grouping-review-link")).toBeVisible();
});

test("experiment corpus shows ScanFailedPage when ingest_status=failed", async ({ page }) => {
  await mockExperiment(page, { id: 7, ingest_status: "failed", samples: 0 });
  await page.goto("/experiments/7/corpus");
  // The failed branch renders ScanFailedPage — the primary action is "Open Configuration".
  await expect(page.getByRole("button", { name: /Open Configuration/i })).toBeVisible();
  // The table must NOT be visible on the failed branch.
  await expect(page.getByRole("table")).toHaveCount(0);
});

test("a live rescan SSE frame drives the corpus ProgressBar (6f)", async ({ page }) => {
  // End-to-end seam: backend broadcast_progress! → /api/events frame →
  // App.tsx ingest listener → ingestInFlight → corpus rescanning ProgressBar.
  // The unit halves are covered separately; this ties the whole chain together.
  await mockExperiment(page, { id: 7, ingest_status: "scanning", samples: 3 });
  // Override /api/events with a single rescan progress frame. EventSource parses
  // it on connect and App.tsx's "curation" listener maps phase:"rescan" →
  // analyzing. Registered AFTER mockExperiment so this handler wins (Playwright
  // invokes the most-recently-added matching route first).
  await page.route("**/api/events*", (r) =>
    r.fulfill({
      status: 200,
      contentType: "text/event-stream",
      body:
        'event: curation\n' +
        'data: {"kind":"ingest_progress","payload":{"experiment_id":7,"processed":5,"total":10,"phase":"rescan"}}\n\n',
    }));
  await page.goto("/experiments/7/corpus");

  const bar = page.getByTestId("live-ingest-slot").getByRole("progressbar");
  await expect(bar).toBeVisible();
  // The bar reflects the FRAME's processed/total, proving the SSE→UI wiring.
  await expect(bar).toHaveAttribute("aria-valuenow", "5");
  await expect(bar).toHaveAttribute("aria-valuemax", "10");
});
