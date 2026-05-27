import { test, expect } from "@playwright/test";

/**
 * figure-export.spec.ts — Playwright coverage for the Copy/Save figure
 * export buttons (spec: 2026-05-08-figure-export-design.md, issue #90).
 *
 * Two scenarios:
 *   1. Copy → clipboard contains image/png   (Chromium-only)
 *   2. Download → file lands; filename matches the `himalaya-trace-…-YYYY-MM-DD.{png,svg}` pattern.
 */

const SAMPLE_FIXTURE = {
  exposure: { id: 1, sample_id: 10, name: "exp1.dat", selected: true },
  trace: {
    q: Array.from({ length: 50 }, (_, i) => 0.05 + i * 0.02),
    I: Array.from({ length: 50 }, (_, i) => 1000 / (1 + i)),
    sigma: Array.from({ length: 50 }, () => 5),
  },
  peaks: [{
    id: 1, exposure_id: 1, q: 0.18, intensity: 420, prominence: 180,
    sharpness: 2.1, source: "auto", excluded: false,
  }],
};

async function mockApi(page: import("@playwright/test").Page): Promise<void> {
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ json: [{ id: 7, name: "JC23", path: "/x", data_dir: "/x", analysis_dir: "/x", manifest_path: null, created_at: "2026", q_units: "A-1" }] }));
  await page.route("**/api/experiments/7", (r) =>
    r.fulfill({ json: { id: 7, name: "JC23", path: "/x", data_dir: "/x", analysis_dir: "/x", manifest_path: null, created_at: "2026", q_units: "A-1" } }));
  await page.route("**/api/experiments/7/samples", (r) =>
    r.fulfill({ json: [{ id: 10, experiment_id: 7, display_name: "S4", name: "Sample 4", notes: null, tags: [] }] }));
  // Corpus list — the focus workspace (/sample/:id, I4.4) learns sample 10's
  // experiment_id from here.
  await page.route("**/api/samples", (r) =>
    r.fulfill({ json: [{ id: 10, experiment_id: 7, display_name: "S4", name: "Sample 4", notes: null, tags: [] }] }));
  await page.route("**/api/samples/10/messages", (r) =>
    r.fulfill({ json: [] }));
  await page.route("**/api/samples/10/exposures*", (r) =>
    r.fulfill({ json: [SAMPLE_FIXTURE.exposure] }));
  await page.route("**/api/exposures/1/trace", (r) =>
    r.fulfill({ json: SAMPLE_FIXTURE.trace }));
  await page.route("**/api/exposures/1/peaks", (r) =>
    r.fulfill({ json: SAMPLE_FIXTURE.peaks }));
  await page.route("**/api/exposures/1/indices", (r) =>
    r.fulfill({ json: [] }));
  await page.route("**/api/exposures/1/groups", (r) =>
    r.fulfill({ json: [] }));
}

async function seedState(page: import("@playwright/test").Page): Promise<void> {
  await page.addInitScript(() => {
    const state = {
      version: 3,
      state: {
        username: "alice",
        tutorialSeen: true,
        activeExperimentId: 7,
        activeSampleId: 10,
        activeExposureId: 1,
        activePage: "compare",
        theme: "dark",
        // Other Zustand fields fall back to defaults from state.ts.
      },
    };
    localStorage.setItem("himalaya-ui:state", JSON.stringify(state));
  });
}

// Clipboard permissions are scoped to this spec only (rather than globally
// in playwright.config.ts) so prompts don't leak into other specs that
// didn't opt in. The Copy test reads the clipboard via page.evaluate; the
// Download tests don't need permissions but inheriting them is harmless.
test.use({ permissions: ["clipboard-read", "clipboard-write"] });

test.describe("Figure export — focus workspace TraceViewer (#181)", () => {
  test("Copy puts a PNG on the clipboard (Chromium)", async ({ page, browserName }) => {
    test.skip(browserName !== "chromium", "Clipboard read requires Chromium permissions");
    await mockApi(page);
    await seedState(page);
    await page.goto("/sample/10");

    const copyBtn = page.getByRole("button", { name: /copy trace plot to clipboard/i });
    await expect(copyBtn).toBeEnabled();
    await copyBtn.click();

    // Wait for the success toast to fire (eats any in-flight render time).
    // Asserting on either the toast text or the clipboard content directly —
    // we go with the clipboard since that's the actual contract.
    await page.waitForTimeout(500);
    const types = await page.evaluate(async () => {
      const items = await navigator.clipboard.read();
      const out: string[] = [];
      for (const item of items) for (const t of item.types) out.push(t);
      return out;
    });
    expect(types).toContain("image/png");
  });

  test("Download → PNG file lands with himalaya-trace-… filename", async ({ page }) => {
    await mockApi(page);
    await seedState(page);
    await page.goto("/sample/10");

    const dlBtn = page.getByRole("button", { name: /download trace plot as png/i });
    await expect(dlBtn).toBeEnabled();

    const downloadPromise = page.waitForEvent("download");
    await dlBtn.click();
    const download = await downloadPromise;

    expect(download.suggestedFilename()).toMatch(/^himalaya-trace-.*\d{4}-\d{2}-\d{2}\.png$/);
  });

  test("Download → SVG via chevron menu", async ({ page }) => {
    await mockApi(page);
    await seedState(page);
    await page.goto("/sample/10");

    const chevron = page.getByRole("button", { name: /other download formats/i });
    await chevron.click();
    const svgRow = page.getByText(/download as svg/i);

    const downloadPromise = page.waitForEvent("download");
    await svgRow.click();
    const download = await downloadPromise;

    expect(download.suggestedFilename()).toMatch(/^himalaya-trace-.*\d{4}-\d{2}-\d{2}\.svg$/);
  });
});
