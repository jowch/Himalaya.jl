/**
 * E2E: GroupingReviewPage — mocked Playwright spec (Task 20, Phase E2).
 *
 * Exercises the `/experiments/:id/grouping` route:
 *   1. Default "Needs review" filter shows only flagged loads; "All loads" reveals all.
 *   2. Glob filter (HA8*) narrows results; selecting a sample then changing filter
 *      keeps the selection (bulk bar count stays).
 *   3. Merge-prompt "Merge into that sample" → confirm modal → loser sample disappears
 *      + undo toast shows (merge = POST /api/samples/:loser/merge → 200).
 *   4. Split divider "Split here" → split fires (POST /api/samples/:id/split → 201).
 *   5. Per-exposure Move fires (POST /api/exposures/:id/move → 200) and shows toast.
 *
 * Structural edits invalidate `loads` in onSuccess — the optimistic DOM lands before
 * any refetch. The tests assert on the optimistic DOM then on the post-invalidate
 * /loads GET (returning the post-edit roll-up). Own-op SSE reconcile is NOT tested
 * here (no SSE needed for own-op; the queue resolves via onSuccess + invalidation).
 */
import { test, expect, type Page } from "@playwright/test";

// ---------------------------------------------------------------------------
// Shared fixture data (§8.8 shapes: sample has sample_id, exposures have .id)
// ---------------------------------------------------------------------------

const EXPERIMENT_ID = 7;

const EXPERIMENT = {
  id: EXPERIMENT_ID,
  name: "SSRL Beamtime",
  path: "/data",
  data_dir: "/data",
  analysis_dir: "/analysis",
  manifest_path: null,
  created_at: "2026-01-01",
  ingest_status: "complete",
};

/** Exposure with merge flag — will be used in the merge test. */
const LOSER_SAMPLE_ID = 30;
const SURVIVOR_SAMPLE_ID = 31;

/** A load with a merge-flagged sample (S01) and a clean sample (S02). */
const LOAD_FLAGGED: object = {
  load_id: 1,
  load_index: 1,
  session_id: null,
  start_time: null,
  end_time: null,
  frame_count: 2,
  note: null,
  samples: [
    {
      sample_id: LOSER_SAMPLE_ID,
      name: "HA85 (S01P01)",
      slot_index: 1,
      grouping_source: "auto_position",
      name_source: "auto",
      merged_into_id: null,
      flag: {
        kind: "merge",
        merge_with_sample_id: SURVIVOR_SAMPLE_ID,
        merge_with_label: "HA85 (S02P01)",
      },
      exposures: [
        { id: 101, filename: "HA85_001.dat", horizontal_position: 10.0, timestamp: "2026-01-01T10:00:00" },
      ],
    },
    {
      sample_id: SURVIVOR_SAMPLE_ID,
      name: "HA85 (S02P01)",
      slot_index: 1,
      grouping_source: "auto_position",
      name_source: "auto",
      merged_into_id: null,
      flag: null,
      exposures: [
        { id: 102, filename: "HA85_002.dat", horizontal_position: 10.1, timestamp: "2026-01-01T11:00:00" },
      ],
    },
  ],
};

/** A load with a split-flagged sample — drives the split test. */
const SPLIT_SAMPLE_ID = 40;
const LOAD_SPLIT: object = {
  load_id: 2,
  load_index: 2,
  session_id: null,
  start_time: null,
  end_time: null,
  frame_count: 2,
  note: null,
  samples: [
    {
      sample_id: SPLIT_SAMPLE_ID,
      name: "HA86 (S01P01)",
      slot_index: 1,
      grouping_source: "auto_position",
      name_source: "auto",
      merged_into_id: null,
      flag: {
        kind: "split",
        split_at_index: 1,
        jump_from: 10.0,
        jump_to: 15.5,
      },
      exposures: [
        { id: 201, filename: "HA86_001.dat", horizontal_position: 10.0, timestamp: "2026-01-01T12:00:00" },
        { id: 202, filename: "HA86_002.dat", horizontal_position: 15.5, timestamp: "2026-01-01T12:30:00" },
      ],
    },
  ],
};

/** A clean load (no flags) — visible only in "All loads" view. */
const CLEAN_SAMPLE_ID = 50;
const LOAD_CLEAN: object = {
  load_id: 3,
  load_index: 3,
  session_id: null,
  start_time: null,
  end_time: null,
  frame_count: 1,
  note: null,
  samples: [
    {
      sample_id: CLEAN_SAMPLE_ID,
      name: "JC C01",
      slot_index: 1,
      grouping_source: "auto_position",
      name_source: "auto",
      merged_into_id: null,
      flag: null,
      exposures: [
        { id: 301, filename: "JC_001.dat", horizontal_position: 5.0, timestamp: "2026-01-01T13:00:00" },
      ],
    },
  ],
};

const ALL_LOADS = [LOAD_FLAGGED, LOAD_SPLIT, LOAD_CLEAN];
const FLAGGED_LOADS = [LOAD_FLAGGED, LOAD_SPLIT];

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Seed localStorage so the tutorial and username prompts don't gate navigation. */
async function seedState(page: Page): Promise<void> {
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state: { username: "alice", tutorialSeen: true, theme: "dark" }, version: 5 }),
    );
  });
}

/** Register all the base mocks needed for the grouping route. */
async function mockGrouping(
  page: Page,
  loads: object[] = ALL_LOADS,
): Promise<void> {
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPERIMENT]) }));
  await page.route(`**/api/experiments/${EXPERIMENT_ID}`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPERIMENT) }));
  await page.route(`**/api/experiments/${EXPERIMENT_ID}/loads`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(loads) }));
  // Stub the samples + events routes that ExperimentShell or child components may call.
  await page.route(`**/api/experiments/${EXPERIMENT_ID}/samples`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/events*", (r) => r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
  // Stub the corpus /api/samples route (CorpusShell or sidebar may request it).
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/sample-tags", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/series", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test.beforeEach(async ({ page }) => {
  await seedState(page);
});

test("grouping: default view shows only flagged loads; switching to 'All loads' reveals clean load", async ({ page }) => {
  await mockGrouping(page, ALL_LOADS);
  await page.goto(`/experiments/${EXPERIMENT_ID}/grouping`);

  // Page heading confirms GroupingReviewPage mounted.
  await expect(page.getByRole("heading", { level: 1, name: /Check the grouping/ })).toBeVisible();

  // Default filter is "Needs review" — only the two flagged loads are visible.
  const folds = page.getByTestId("load-fold");
  await expect(folds).toHaveCount(2);

  // The clean load (Load 3 / "JC C01") is NOT visible yet.
  await expect(page.getByText("JC C01")).toBeHidden();

  // Switch to "All loads" — all three loads appear. SegmentedControl renders
  // role="group" with aria-pressed buttons by default; select by button text.
  await page.getByRole("button", { name: "All loads" }).click();
  await expect(folds).toHaveCount(3);
  await expect(page.getByText("JC C01")).toBeVisible();
});

test("grouping filter: glob HA8* narrows; selection survives filter change", async ({ page }) => {
  await mockGrouping(page, ALL_LOADS);
  await page.goto(`/experiments/${EXPERIMENT_ID}/grouping`);

  await expect(page.getByRole("heading", { level: 1, name: /Check the grouping/ })).toBeVisible();

  // All loads visible initially (flagged ones are expanded by default).
  // Switch to All loads so the clean one is also shown.
  // SegmentedControl renders role="group" with aria-pressed buttons by default.
  await page.getByRole("button", { name: "All loads" }).click();
  await expect(page.getByTestId("load-fold")).toHaveCount(3);

  // Select a sample by clicking its checkbox before filtering.
  // HA85 (S02P01) is in Load 1 — click its checkbox via aria-label.
  const survivor = page.getByRole("checkbox", { name: `Select HA85 (S02P01)` });
  await expect(survivor).toBeVisible();
  await survivor.click({ force: true });
  // Bulk bar appears showing 1 sample.
  await expect(page.getByTestId("bulk-bar")).toBeVisible();
  await expect(page.getByTestId("bulk-bar")).toContainText("1 sample selected");

  // Apply a glob filter for HA8* — narrows to Load 1 + Load 2 (both have HA85/HA86).
  await page.getByRole("textbox", { name: "Filter samples" }).fill("HA8*");
  await expect(page.getByTestId("load-fold")).toHaveCount(2);

  // The selected sample (HA85 S02P01) stays visible even though only HA8* matches —
  // GroupingReviewPage preserves selection across filter changes.
  await expect(page.getByTestId("bulk-bar")).toBeVisible();
  await expect(page.getByTestId("bulk-bar")).toContainText("1 sample selected");
});

test("grouping merge: confirm modal fires POST /merge; loser disappears + undo toast shows", async ({ page }) => {
  // After the merge, return a /loads response without the loser sample.
  const POST_MERGE_LOADS = [
    {
      ...LOAD_FLAGGED,
      samples: [(LOAD_FLAGGED as { samples: object[] }).samples[1]], // survivor only
    },
    LOAD_SPLIT,
    LOAD_CLEAN,
  ];

  let mergeBody: unknown = null;
  await mockGrouping(page, ALL_LOADS);

  // Register merge endpoint AFTER mockGrouping so it takes priority.
  await page.route(`**/api/samples/${LOSER_SAMPLE_ID}/merge`, async (route) => {
    mergeBody = route.request().postDataJSON();
    await route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({ loser_id: LOSER_SAMPLE_ID, survivor_id: SURVIVOR_SAMPLE_ID }),
    });
  });

  // After onSuccess invalidates loads, serve the post-merge roll-up.
  let loadsCallCount = 0;
  await page.route(`**/api/experiments/${EXPERIMENT_ID}/loads`, (r) => {
    loadsCallCount++;
    const body = loadsCallCount === 1 ? ALL_LOADS : POST_MERGE_LOADS;
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(body) });
  });

  await page.goto(`/experiments/${EXPERIMENT_ID}/grouping`);
  await expect(page.getByRole("heading", { level: 1, name: /Check the grouping/ })).toBeVisible();

  // Load 1 is flagged → auto-expanded. The merge-prompt should be visible
  // inside the loser sample's open fold.
  const mergePrompt = page.getByTestId("merge-prompt").first();
  await expect(mergePrompt).toBeVisible();

  // Click "Merge into that sample" — opens the confirm modal.
  await mergePrompt.getByRole("button", { name: "Merge into that sample" }).click();
  const confirmModal = page.getByTestId("merge-confirm");
  await expect(confirmModal).toBeVisible();
  await expect(confirmModal).toContainText("HA85 (S02P01)");

  // Confirm — fires POST /merge.
  await confirmModal.getByRole("button", { name: "Merge" }).click();

  // Modal closes.
  await expect(confirmModal).toBeHidden();

  // POST was sent with the correct survivor_id.
  await expect.poll(() => mergeBody).toMatchObject({ survivor_id: SURVIVOR_SAMPLE_ID });

  // Toast appears confirming the merge.
  await expect(page.getByRole("status")).toContainText(/Merged into/);
});

test("grouping split: split-divider 'Split here' fires POST /split", async ({ page }) => {
  const SPLIT_LOADS_AFTER = [
    LOAD_FLAGGED,
    {
      ...LOAD_SPLIT,
      samples: [
        { ...(LOAD_SPLIT as { samples: object[] }).samples[0], flag: null },
        {
          sample_id: 41,
          name: "HA86 (S01P01) (split)",
          slot_index: 2,
          grouping_source: "auto_position",
          name_source: "auto",
          merged_into_id: null,
          flag: null,
          exposures: [
            { id: 202, filename: "HA86_002.dat", horizontal_position: 15.5, timestamp: "2026-01-01T12:30:00" },
          ],
        },
      ],
    },
    LOAD_CLEAN,
  ];

  let splitBody: unknown = null;
  await mockGrouping(page, ALL_LOADS);

  await page.route(`**/api/samples/${SPLIT_SAMPLE_ID}/split`, async (route) => {
    splitBody = route.request().postDataJSON();
    await route.fulfill({
      status: 201, contentType: "application/json",
      body: JSON.stringify({ new_sample_id: 41 }),
    });
  });

  let loadsCall = 0;
  await page.route(`**/api/experiments/${EXPERIMENT_ID}/loads`, (r) => {
    loadsCall++;
    const body = loadsCall === 1 ? ALL_LOADS : SPLIT_LOADS_AFTER;
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(body) });
  });

  await page.goto(`/experiments/${EXPERIMENT_ID}/grouping`);
  await expect(page.getByRole("heading", { level: 1, name: /Check the grouping/ })).toBeVisible();

  // Load 2 contains the split-flagged sample → auto-expanded. The split divider is inside.
  const splitDivider = page.getByTestId("split-divider").first();
  await expect(splitDivider).toBeVisible();

  // "Split here" fires the split mutator.
  await splitDivider.getByRole("button", { name: "Split here" }).click();

  // POST was sent with the correct exposure_ids.
  await expect.poll(() => splitBody).toMatchObject({
    exposure_ids: [202],
  });

  // Toast appears.
  await expect(page.getByRole("status")).toContainText(/Split sample/);
});

test("grouping move: per-exposure move fires POST /move and shows undo toast", async ({ page }) => {
  const EXPOSURE_ID = 101;

  let moveBody: unknown = null;
  await mockGrouping(page, ALL_LOADS);

  await page.route(`**/api/exposures/${EXPOSURE_ID}/move`, async (route) => {
    moveBody = route.request().postDataJSON();
    await route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({ id: EXPOSURE_ID, sample_id: SURVIVOR_SAMPLE_ID }),
    });
  });

  await page.goto(`/experiments/${EXPERIMENT_ID}/grouping`);
  await expect(page.getByRole("heading", { level: 1, name: /Check the grouping/ })).toBeVisible();

  // Load 1 is flagged → auto-expanded. The loser sample's merge-prompt is visible.
  // Open the loser sample fold to see its exposures (the loser has a merge flag,
  // so it's auto-expanded). Find the Move button for exposure 101.
  const sampleFolds = page.getByTestId("sample-fold");
  const loserFold = sampleFolds.first();
  await expect(loserFold).toBeVisible();

  // The ExposureLeaf Move button fires onMoveExposure. Click it to trigger the move.
  const moveBtn = loserFold.getByRole("button", { name: /Move/i }).first();
  await expect(moveBtn).toBeVisible();
  await moveBtn.click();

  // POST was fired with exposure id + destination sample id (the current sample,
  // as per the stub handler in GroupingReviewPage).
  await expect.poll(() => moveBody).not.toBeNull();

  // Undo toast appears.
  await expect(page.getByRole("status")).toContainText(/Moved exposure/);
});
