import { test, expect, type Page } from "@playwright/test";

// M-A Task 8 (new-series creation flow): the full mocked B-path. Drives the
// contact-sheet sample-grain checkbox column → ComposeBar → "+ New series" →
// the pre-seeded SeriesScopingPage → "Confirm & build" → the two sequential
// writes (POST /api/samples/tags/batch then POST /api/series) → navigation to
// /series/:id.
//
// Greenfield DOM only (src/print/**): SamplesPage + SheetTable + SampleTableRow
// + Checkbox + ComposeBar on the picker side; SeriesScopingPage + ScopePlate +
// ScopeSampleRow on the scoping side. All assertions are on stable data-* /
// role selectors per e2e/AGENTS.md.
//
// Critical pin (commit 0ded204): the two POSTs carry DISTINCT X-Client-Op-Id
// headers — the tags batch under the op's bare client_op_id, the series create
// under a `${op}:series`-suffixed id. Without that, the backend's
// client_op_id-keyed idempotency cache would replay the tags 201 in place of
// creating the series. This spec captures both headers and asserts they differ.

const EXPERIMENT = {
  id: 1, name: "SSRL Test", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};

/** Two samples that both carry a "ratio" tag → proposeOrdering proposes
 *  "ratio" as the ordering variable (the warm ScopePlate path). */
const SAMPLES = [
  {
    id: 10, experiment_id: 1, display_name: "DOPC-1:1", name: "run10",
    notes: null, q_units: "A-1",
    tags: [{ id: 1, key: "ratio", value: "1:1", source: "manual" }],
  },
  {
    id: 11, experiment_id: 1, display_name: "DOPC-2:1", name: "run11",
    notes: null, q_units: "A-1",
    tags: [{ id: 2, key: "ratio", value: "2:1", source: "manual" }],
  },
];
const EXPOSURES_10 = [
  {
    id: 1, sample_id: 10, filename: "f1.dat", kind: "file", selected: true,
    status: "accepted", image_path: null, image_version: "", tags: [],
    sources: [], trace_hash: null, analysis_inputs_hash: null,
  },
];
const EXPOSURES_11 = [
  {
    id: 2, sample_id: 11, filename: "f2.dat", kind: "file", selected: true,
    status: "accepted", image_path: null, image_version: "", tags: [],
    sources: [], trace_hash: null, analysis_inputs_hash: null,
  },
];
const NEW_SERIES = {
  id: 42, title: "Series by ratio", description: null, state: "draft",
  ordering_variable: "ratio", order_rule: "asc", forked_from_id: null,
  forked_at_hash: null, content_hash: "xyz",
  view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
  members: [], samples: [],
};

interface Captured {
  tagsBatch: Record<string, unknown> | null;
  tagsOpId: string | null;
  seriesCreate: Record<string, unknown> | null;
  seriesOpId: string | null;
}

function freshCaptured(): Captured {
  return { tagsBatch: null, tagsOpId: null, seriesCreate: null, seriesOpId: null };
}

/** Seed identity so the "Who are you?" onboarding modal never blocks the page.
 *  Mirrors corpus-culling / series-scoping specs: version 3, theme dark. */
async function seedIdentity(page: Page): Promise<void> {
  await page.addInitScript(() =>
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({
        state: { username: "alice", theme: "dark", tutorialSeen: true },
        version: 3,
      }),
    ),
  );
}

/** Mock the common reads shared by every test (identity, corpus, SSE drain). */
async function mockCommon(page: Page): Promise<void> {
  await seedIdentity(page);
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPERIMENT]) }));
  await page.route("**/api/experiments/1", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPERIMENT) }));
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(SAMPLES) }));
  await page.route("**/api/samples/10/exposures*", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPOSURES_10) }));
  await page.route("**/api/samples/11/exposures*", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPOSURES_11) }));
  // SSE drain — empty event-stream so the App-level EventSource doesn't hang.
  await page.route("**/api/events", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
}

/** Mock the scoping reads + the two write endpoints + the builder reads.
 *  `warm` chooses whether the corpus carries a shared "ratio" tag (warm =
 *  ScopePlate proposes an ordering variable) or none (cold = ColdAssignPanel). */
async function mockScopingAndWrites(
  page: Page,
  captured: Captured,
  warm: boolean,
): Promise<void> {
  // Corpus-wide tag pairs — drive proposeOrdering. Warm: "ratio" appears twice.
  await page.route("**/api/sample-tags", (r) =>
    r.fulfill({
      status: 200, contentType: "application/json",
      body: warm
        ? JSON.stringify([{ key: "ratio", value: "1:1" }, { key: "ratio", value: "2:1" }])
        : "[]",
    }));
  // Picker rows — both seeded samples. Warm rows carry the "ratio" tag; cold
  // rows carry none (so no ordering variable is proposable).
  await page.route("**/api/picker-samples", (r) =>
    r.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify([
        { sample: warm ? SAMPLES[0] : { ...SAMPLES[0], tags: [] }, indexing_exposure_id: null, all_exposures: [] },
        { sample: warm ? SAMPLES[1] : { ...SAMPLES[1], tags: [] }, indexing_exposure_id: null, all_exposures: [] },
      ]),
    }));

  // Write #1 — the scoping tag batch. Capture body + X-Client-Op-Id header.
  await page.route("**/api/samples/tags/batch", async (route) => {
    captured.tagsBatch = route.request().postDataJSON();
    captured.tagsOpId = route.request().headers()["x-client-op-id"] ?? null;
    await route.fulfill({
      status: 201, contentType: "application/json",
      body: JSON.stringify([
        { id: 1, sample_id: 10, key: "ratio", value: "1:1", source: "scoping" },
        { id: 2, sample_id: 11, key: "ratio", value: "2:1", source: "scoping" },
      ]),
    });
  });

  // Write #2 — the series create (POST) / folio listing (GET). Capture body +
  // X-Client-Op-Id header on the POST.
  await page.route("**/api/series", async (route) => {
    if (route.request().method() === "POST") {
      captured.seriesCreate = route.request().postDataJSON();
      captured.seriesOpId = route.request().headers()["x-client-op-id"] ?? null;
      await route.fulfill({
        status: 201, contentType: "application/json", body: JSON.stringify(NEW_SERIES),
      });
    } else {
      await route.fulfill({ status: 200, contentType: "application/json", body: "[]" });
    }
  });

  // Builder reads after navigation to /series/42.
  await page.route("**/api/series/42", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(NEW_SERIES) }));
  await page.route("**/api/series/42/traces", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "{}" }));
  await page.route("**/api/series/42/forks", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
}

test.describe("new-series creation flow (M-A)", () => {
  test("full B-path: check samples → compose bar → scoping → confirm → /series/:id", async ({ page }) => {
    const captured = freshCaptured();
    await mockCommon(page);
    await mockScopingAndWrites(page, captured, /* warm */ true);
    await page.goto("/samples");

    // The contact sheet loads with both sample rows.
    await expect(page.getByTestId("samples-page")).toBeVisible();
    await expect(page.getByTestId("sample-table-row")).toHaveCount(2);

    // Compose bar starts hidden (no sample checked).
    await expect(page.getByTestId("compose-bar")).toHaveAttribute("data-show", "false");

    // Check the first sample → compose bar appears with count 1.
    await page.getByTestId("sample-table-row").first().getByRole("checkbox").click();
    await expect(page.getByTestId("compose-bar")).toHaveAttribute("data-show", "true");
    await expect(page.getByTestId("compose-bar")).toContainText("1 sample");

    // Check the second sample → count 2.
    await page.getByTestId("sample-table-row").nth(1).getByRole("checkbox").click();
    await expect(page.getByTestId("compose-bar")).toContainText("2 samples");

    // The frame-grain CullBar stays hidden — sample-grain selection is distinct.
    await expect(page.getByTestId("cull-bar")).toHaveAttribute("data-show", "false");

    // "+ New series" → land on /series/new, pre-seeded to the two checked samples.
    await page.getByRole("button", { name: /new series/i }).click();
    await expect(page).toHaveURL(/\/series\/new/);
    await expect(page.getByTestId("scoping-page")).toBeVisible();
    // Seeded: exactly the 2 selected samples become scope rows (warm path).
    await expect(page.getByTestId("scope-sample-row")).toHaveCount(2);

    // Confirm & build is enabled (both reads kept).
    const buildBtn = page.getByRole("button", { name: /confirm & build/i });
    await expect(buildBtn).toBeEnabled();
    await buildBtn.click();

    // Both writes fire: tags batch with key="ratio" + source="scoping".
    await expect.poll(() => captured.tagsBatch?.["key"], { timeout: 4000 }).toBe("ratio");
    await expect.poll(() => captured.tagsBatch?.["source"]).toBe("scoping");
    // …then the series create with a title mentioning the ordering variable.
    await expect.poll(() => captured.seriesCreate?.["title"], { timeout: 4000 }).toMatch(/ratio/i);

    // CRITICAL PIN: the two POSTs carry DISTINCT X-Client-Op-Id headers — the
    // series create is the tags op-id suffixed with ":series". This is the
    // idempotency fix (without it, the create would replay the tags 201).
    await expect.poll(() => captured.seriesOpId).not.toBeNull();
    expect(captured.tagsOpId).not.toBeNull();
    expect(captured.seriesOpId).not.toBe(captured.tagsOpId);
    expect(captured.seriesOpId).toBe(`${captured.tagsOpId}:series`);

    // Navigated to the new series builder at /series/42.
    await expect(page).toHaveURL(/\/series\/42/);
  });

  test("Clear on the compose bar resets the sample selection", async ({ page }) => {
    const captured = freshCaptured();
    await mockCommon(page);
    await mockScopingAndWrites(page, captured, /* warm */ true);
    await page.goto("/samples");
    await expect(page.getByTestId("sample-table-row")).toHaveCount(2);

    await page.getByTestId("sample-table-row").first().getByRole("checkbox").click();
    await expect(page.getByTestId("compose-bar")).toHaveAttribute("data-show", "true");

    await page.getByTestId("compose-bar").getByRole("button", { name: /clear/i }).click();
    await expect(page.getByTestId("compose-bar")).toHaveAttribute("data-show", "false");
  });

  test("direct visit to /series/new shows the full corpus (no seed)", async ({ page }) => {
    const captured = freshCaptured();
    await mockCommon(page);
    await mockScopingAndWrites(page, captured, /* warm */ true);
    // No checker round-trip → location.state has no seed → full-corpus path.
    await page.goto("/series/new");
    await expect(page.getByTestId("scoping-page")).toBeVisible();
    // Both corpus samples appear as scope rows (unfiltered).
    await expect(page.getByTestId("scope-sample-row")).toHaveCount(2);
  });

  test("cold-corpus path: no proposable variable → ColdAssignPanel, build disabled", async ({ page }) => {
    const captured = freshCaptured();
    await mockCommon(page);
    await mockScopingAndWrites(page, captured, /* warm */ false);

    // Arrive via the picker so a seed is carried (the cold panel only renders
    // when the user deliberately selected samples but none share a tag key).
    await page.goto("/samples");
    await expect(page.getByTestId("sample-table-row")).toHaveCount(2);
    await page.getByTestId("sample-table-row").first().getByRole("checkbox").click();
    await page.getByTestId("sample-table-row").nth(1).getByRole("checkbox").click();
    await page.getByRole("button", { name: /new series/i }).click();

    // No ordering variable proposable → the cold assign panel appears.
    await expect(page.getByTestId("cold-assign-panel")).toBeVisible();
    // Confirm & build is disabled until the key + every value are filled.
    await expect(page.getByRole("button", { name: /confirm & build/i })).toBeDisabled();
  });
});
