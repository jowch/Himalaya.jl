import { test, expect, type Page } from "@playwright/test";

// I3.4 (#174): the Phase-3 scoping mocked spec. Exercises /series/new —
// propose an ordering variable from corpus tags, confirm-and-build, batch
// write with source='scoping', then land on the folio.
//
// Migrated to the greenfield DOM (SeriesScopingPage + ScopePlate +
// ScopeSampleRow + FlagButton): no inline value editing, no ScopingConfirmModal,
// no "+ Add to series" button. The only durable effect is POSTing
// /api/samples/tags/batch and navigating to /series on success.

const EXPERIMENT = {
  id: 1, name: "SSRL Test", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};

/** A picker row whose sample has a known value for the ordering key ("ratio"). */
function memberRow(id: number, name: string, ratio: string) {
  return {
    sample: {
      id, experiment_id: 1, name, display_name: null, notes: null,
      tags: [{ id, key: "ratio", value: ratio, source: "manual" }],
    },
    indexing_exposure_id: null,
    all_exposures: [],
  };
}

/** A picker row whose sample lacks the ordering key entirely — becomes a candidate. */
function candidateRow(id: number, name: string) {
  return {
    sample: {
      id, experiment_id: 1, name, display_name: null, notes: null,
      tags: [],
    },
    indexing_exposure_id: null,
    all_exposures: [],
  };
}

async function mockScoping(page: Page, captured: { body: unknown }): Promise<void> {
  // Identity-seed: suppress the "Who are you?" onboarding modal (mirrors folio spec).
  await page.addInitScript(() =>
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state: { username: "alice", theme: "dark", tutorialSeen: true }, version: 3 }),
    ),
  );

  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([EXPERIMENT]) }));

  // Corpus-wide tag pairs. "ratio" appears twice — it will be proposed as the
  // ordering variable (highest frequency, deterministic per proposeOrdering).
  await page.route("**/api/sample-tags", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([
        { key: "ratio", value: "1:1" },
        { key: "ratio", value: "2:1" },
      ]) }));

  // Two member samples (have "ratio") + one candidate sample (lacks "ratio").
  await page.route("**/api/picker-samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([
        memberRow(10, "A_1to1", "1:1"),
        memberRow(11, "B_2to1", "2:1"),
        candidateRow(12, "C_unlabelled"),
      ]) }));

  // No indexing exposures in this mock — traces and indices never fetched.
  // (exposureIds will be empty because all indexing_exposure_id are null.)

  // The batch write endpoint — capture the POST body for assertions.
  await page.route("**/api/samples/tags/batch", async (route) => {
    captured.body = route.request().postDataJSON();
    await route.fulfill({ status: 201, contentType: "application/json",
      body: JSON.stringify([
        { id: 1, sample_id: 10, key: "ratio", value: "1:1", source: "scoping" },
        { id: 2, sample_id: 11, key: "ratio", value: "2:1", source: "scoping" },
      ]) });
  });

  // The folio listing that /series renders after the successful write.
  await page.route("**/api/series", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));

  // SSE endpoint — drain immediately so the App-level EventSource doesn't hang
  // (see e2e/AGENTS.md + series-folio.spec.ts).
  await page.route("**/api/events", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
}

test.describe("series scoping — greenfield DOM", () => {
  test("renders member rows and the candidate list for unlabelled samples", async ({ page }) => {
    const captured = { body: null as unknown };
    await mockScoping(page, captured);
    await page.goto("/series/new");

    // The root scoping-page wrapper must be visible.
    await expect(page.getByTestId("scoping-page")).toBeVisible();

    // Two member rows (samples with "ratio" value) are rendered.
    await expect(page.getByTestId("scope-sample-row")).toHaveCount(2);

    // The informational candidate list appears for the unlabelled sample (C_unlabelled).
    await expect(page.getByTestId("scope-candidates")).toBeVisible();
    await expect(page.getByTestId("scope-candidate")).toHaveCount(1);

    // Candidate list carries NO add button (Option A: candidates are discovery-only).
    await expect(page.getByTestId("scope-candidates").getByRole("button")).toHaveCount(0);
  });

  test("Confirm & build is enabled when members exist and navigates to /series on success", async ({ page }) => {
    const captured = { body: null as unknown };
    await mockScoping(page, captured);
    await page.goto("/series/new");

    await expect(page.getByTestId("scoping-page")).toBeVisible();

    // The build button is enabled (both rows kept — none flagged by default).
    const buildBtn = page.getByRole("button", { name: /confirm & build/i });
    await expect(buildBtn).toBeEnabled();

    // Click confirm — triggers POST /api/samples/tags/batch.
    await buildBtn.click();

    // The batch must have been posted with key="ratio" and source="scoping".
    await expect.poll(() => (captured.body as Record<string, unknown> | null)?.["key"], { timeout: 3000 }).toBe("ratio");
    await expect.poll(() => (captured.body as Record<string, unknown> | null)?.["source"]).toBe("scoping");

    // On success the page navigates to the series folio.
    await expect(page).toHaveURL(/\/series$/);
    // Greenfield folio: folio-header is the mount signal.
    await expect(page.getByTestId("folio-header")).toBeVisible();
  });

  test("skipping a member via its value control excludes it from the batch payload", async ({ page }) => {
    const captured = { body: null as unknown };
    await mockScoping(page, captured);
    await page.goto("/series/new");

    await expect(page.getByTestId("scoping-page")).toBeVisible();

    // Click the first row's FlagButton to skip it (toggles data-state ok→flagged,
    // meaning "skip this read from the write").
    const firstFlagBtn = page.getByTestId("scope-sample-row").first().getByTestId("flag-button");
    await expect(firstFlagBtn).toHaveAttribute("data-state", "ok");
    await firstFlagBtn.click();
    await expect(firstFlagBtn).toHaveAttribute("data-state", "flagged");

    // Build is still enabled (one kept member remains).
    const buildBtn = page.getByRole("button", { name: /confirm & build/i });
    await expect(buildBtn).toBeEnabled();
    await buildBtn.click();

    // The batch payload must exclude the skipped sample (only 1 tag written).
    await expect.poll(
      () => {
        const b = captured.body as Record<string, unknown> | null;
        const tags = b?.["tags"];
        return Array.isArray(tags) ? tags.length : null;
      },
      { timeout: 3000 },
    ).toBe(1);

    await expect(page).toHaveURL(/\/series$/);
  });

  test("Confirm & build is disabled when all members are skipped", async ({ page }) => {
    const captured = { body: null as unknown };
    await mockScoping(page, captured);
    await page.goto("/series/new");

    await expect(page.getByTestId("scoping-page")).toBeVisible();

    // Skip both member rows.
    const flagBtns = page.getByTestId("scope-sample-row").getByTestId("flag-button");
    await flagBtns.nth(0).click();
    await flagBtns.nth(1).click();

    // Build gate must close — nothing kept.
    await expect(page.getByRole("button", { name: /confirm & build/i })).toBeDisabled();
  });

  test("scoping-discard navigates to /series without writing", async ({ page }) => {
    const captured = { body: null as unknown };
    await mockScoping(page, captured);
    await page.goto("/series/new");

    await expect(page.getByTestId("scoping-page")).toBeVisible();

    // Click the discard button — must navigate to /series without any batch POST.
    await page.getByTestId("scoping-discard").click();
    await expect(page).toHaveURL(/\/series$/);

    // No batch was written.
    expect(captured.body).toBeNull();
  });
});
