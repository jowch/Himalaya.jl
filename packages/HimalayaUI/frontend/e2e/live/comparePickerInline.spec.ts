/**
 * Live-mode spec for the inline picker panel on the Compare edit page.
 *
 * Exercises the roundtrip: open /experiments/:eid/compare/new → see the
 * picker panel mounted in the right slot → click a sample row's checkbox →
 * see the draft member land in the plot gutter → expand the override caret →
 * pick a second exposure → observe the row's resolved exposure id update.
 *
 * Lives in `e2e/live/` because the picker needs the real
 * `GET /api/experiments/:eid/picker-samples` response from the backend.
 * The mocked suite can cover this shape with `page.route`, but the live
 * spec confirms the actual endpoint payload drives the panel correctly.
 *
 * Operator setup (see README.md):
 *   - Backend on port 8090 with the test DB.
 *   - Vite on port 5180 proxying /api → 8090.
 *   Run: npm run e2e:live
 */
import { test, expect, type Page } from "@playwright/test";

const BACKEND_BASE = process.env["BACKEND_BASE"] ?? "http://127.0.0.1:8090";

type Sample = { id: number; experiment_id: number };
type PickerSampleExposure = { id: number; sample_id: number; filename: string | null; selected: boolean };
type PickerSampleRow = {
  sample: Sample;
  indexing_exposure_id: number | null;
  all_exposures: PickerSampleExposure[];
};

interface Fixture {
  experimentId: number;
  /** A sample row with a non-null indexing_exposure_id (required for checkbox). */
  pickerRow: PickerSampleRow;
}

/**
 * Find the first experiment that has at least one picker-sample row with
 * a non-null `indexing_exposure_id`. Rows with null are rendered as
 * disabled (no checkbox) and cannot be checked.
 */
async function findFixture(): Promise<Fixture> {
  const exps = await fetch(`${BACKEND_BASE}/api/experiments`).then(r => r.json());
  for (const exp of exps) {
    const rows: PickerSampleRow[] = await fetch(
      `${BACKEND_BASE}/api/experiments/${exp.id}/picker-samples`,
    ).then(r => r.json());
    const usable = rows.find(r => r.indexing_exposure_id !== null);
    if (usable) {
      return { experimentId: exp.id, pickerRow: usable };
    }
  }
  throw new Error(
    `No experiment with a picker-sample row with indexing_exposure_id set at ${BACKEND_BASE}. ` +
    `Run \`himalaya analyze\` against the test DB first.`,
  );
}

/**
 * Seed Zustand-persisted state so:
 *   - OnboardingFlow is skipped (username + tutorialSeen).
 *   - The active experiment is pre-set (compare picker scopes by experimentId).
 * Then navigate to "/" first for the SSE subscriber to register, then go to
 * the compare/new route.
 */
async function seedAndNavigate(page: Page, fx: Fixture): Promise<void> {
  await page.addInitScript((args) => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: {
        username: args.username,
        firstName: args.username,
        lastName: "tester",
        activeExperimentId: args.expId,
        activePage: "compare",
        tutorialSeen: true,
        theme: "dark",
      },
      version: 3,
    }));
  }, { username: "compare-picker-tester", expId: fx.experimentId });

  // "/" → SSE warmup per live runbook (subscribers register ~50–200ms after GET).
  await page.goto("/");
  await page.waitForTimeout(800);

  // Navigate to the create shell for this experiment.
  await page.goto(`/experiments/${fx.experimentId}/compare/new`);
}

test.describe("Compare picker inline panel — live mode", () => {
  let fx: Fixture;

  test.beforeAll(async () => {
    fx = await findFixture();
    console.log(
      `[compare-picker] using experiment=${fx.experimentId} ` +
      `sample=${fx.pickerRow.sample.id} ` +
      `indexing_exposure=${fx.pickerRow.indexing_exposure_id}`,
    );
  });

  test("picker panel is visible and checking a row adds a draft member", async ({ page }) => {
    await seedAndNavigate(page, fx);

    // Right slot: picker panel must be present.
    await expect(page.getByTestId("comparison-picker-panel")).toBeVisible();

    // Find the row for our fixture sample. Wait for picker rows to populate
    // (the picker query runs after mount).
    // `data-sample-id` is on the same element as `data-testid="sample-picker-row"`,
    // so `has:` doesn't apply (it expects descendants). Match both attrs in one selector.
    const sampleRow = page.locator(
      `[data-testid="sample-picker-row"][data-sample-id="${fx.pickerRow.sample.id}"]`,
    );
    await expect(sampleRow).toBeVisible({ timeout: 8000 });

    // Before checking: the plot empty-state placeholder should be visible.
    await expect(page.getByTestId("compare-edit-plot-empty")).toBeVisible();

    // Click the checkbox.
    await sampleRow.getByTestId("sample-picker-row-checkbox").click();

    // Row should now be checked.
    await expect(
      sampleRow.getByTestId("sample-picker-row-checkbox"),
    ).toBeChecked();

    // Adding a member replaces the empty-state with the plot + gutter.
    // The plot-empty placeholder should hide immediately on the next render.
    await expect(page.getByTestId("compare-edit-plot-empty")).not.toBeVisible();

    // The gutter renders inside the Skeleton wrapper that gates on `tracesLoading`,
    // so we wait long enough for the cold trace fetch to land. Live mode against
    // real exposures can take several seconds per trace on a cold cache.
    await expect(page.getByTestId("compare-edit-gutter")).toBeVisible({ timeout: 30_000 });
  });

  test("override caret expands exposure list and switching radio updates data-exposure-id", async ({ page }) => {
    await seedAndNavigate(page, fx);

    // Skip test if the fixture sample has only one exposure — nothing to switch to.
    if (fx.pickerRow.all_exposures.length < 2) {
      console.log(
        `[compare-picker] skipping override-caret test: ` +
        `sample ${fx.pickerRow.sample.id} has only 1 exposure`,
      );
      test.skip();
      return;
    }

    // Wait for the picker row to load.
    // `data-sample-id` is on the same element as `data-testid="sample-picker-row"`,
    // so `has:` doesn't apply (it expects descendants). Match both attrs in one selector.
    const sampleRow = page.locator(
      `[data-testid="sample-picker-row"][data-sample-id="${fx.pickerRow.sample.id}"]`,
    );
    await expect(sampleRow).toBeVisible({ timeout: 8000 });

    // Read the initial resolved exposure id (set by the component to indexing_exposure_id).
    const initialEid = await sampleRow.getAttribute("data-exposure-id");
    expect(initialEid).not.toBeNull();

    // Check the row (not required for the caret, but matches realistic usage).
    await sampleRow.getByTestId("sample-picker-row-checkbox").click();
    await expect(sampleRow.getByTestId("sample-picker-row-checkbox")).toBeChecked();

    // Expand the override list.
    await sampleRow.getByTestId("sample-picker-row-caret").click();
    await expect(sampleRow.locator('[role="radiogroup"]')).toBeVisible();

    // The radiogroup should contain at least two radios.
    const radios = sampleRow.locator('input[type="radio"]');
    await expect(radios).toHaveCount(fx.pickerRow.all_exposures.length);

    // Click the second radio (index 1). The first radio corresponds to the
    // default `indexing_exposure_id`; the second is a different exposure.
    await radios.nth(1).click();
    await expect(radios.nth(1)).toBeChecked();

    // The row's data-exposure-id should now reflect the new selection.
    // We compare against the second exposure's id from the fixture.
    const secondExposureId = String(fx.pickerRow.all_exposures[1]!.id);
    await expect(sampleRow).toHaveAttribute("data-exposure-id", secondExposureId);
    expect(secondExposureId).not.toBe(initialEid);
  });
});
