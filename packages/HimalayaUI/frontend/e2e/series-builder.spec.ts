import { test, expect, type Page } from "@playwright/test";

const MEMBER = {
  id: 1, series_id: 5, exposure_id: 101, display_order: 0,
  band_height: 1, y_offset: 0, normalization: "max",
  color_override: null, label_override: null,
  q_window_min: null, q_window_max: null, peak_display: null,
  snapshot: null, is_stale: false, created_by: 1, created_at: null,
};

const SERIES = {
  id: 5, title: "LL37 titration", description: null, content_hash: "h",
  created_by: 1, created_at: "2026-05-01", updated_at: "2026-05-02",
  forked_from_id: null, forked_at_hash: null, forked_from_title: null,
  view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
  ordering_variable: "LL37 : lipid ratio", order_rule: "ascending",
  state: "committed", members: [MEMBER], samples: [],
};

const TRACE = {
  q: Array.from({ length: 50 }, (_, i) => 0.05 + i * 0.02),
  I: Array.from({ length: 50 }, (_, i) => 1000 / (1 + i)),
  sigma: Array.from({ length: 50 }, () => 5),
};

const EXPOSURE = {
  id: 101, sample_id: 10, filename: "exp1.dat", kind: "file",
  selected: true, status: "accepted", image_path: null, image_version: "",
  tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
};

const SAMPLE = {
  id: 10, experiment_id: 1, name: "JC042", display_name: "DOPE 80%",
  notes: null, tags: [], q_units: "A-1",
};

async function mockCore(page: Page): Promise<void> {
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([{ id: 1, username: "jc" }]) }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/series/5", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(SERIES) }));
  await page.route("**/api/exposures/101/trace", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(TRACE) }));
  await page.route("**/api/exposures/101", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPOSURE) }));
  await page.route("**/api/samples/10", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(SAMPLE) }));
  // SSE endpoint — drain immediately so the multiplayer EventSource doesn't
  // hang the page (see e2e/AGENTS.md + compare.spec.ts).
  await page.route("**/api/events", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
}

/** Seed a username + tutorialSeen so the first-run onboarding overlay (which
 *  intercepts pointer events) does not block clicks. Mirrors figure-export. */
async function seedState(page: Page): Promise<void> {
  await page.addInitScript(() => {
    const state = {
      version: 3,
      state: { username: "jc", tutorialSeen: true, theme: "dark" },
    };
    localStorage.setItem("himalaya-ui:state", JSON.stringify(state));
  });
}

test.describe("series builder", () => {
  test.use({ permissions: ["clipboard-read", "clipboard-write"] });

  test("renders the builder waterfall under the corpus shell", async ({ page }) => {
    await mockCore(page);
    await seedState(page);
    await page.goto("/series/5");
    await expect(page.getByTestId("series-builder-page")).toBeVisible();
    await expect(page.getByTestId("series-builder-plot")).toBeVisible();
    await expect(page.getByTestId("series-builder-rail")).toBeVisible();
    // heatmap is the deferred (#208) option — present but disabled
    await expect(page.getByTestId("repr-heatmap")).toBeDisabled();
  });

  test("Copy puts a PNG on the clipboard (Chromium)", async ({ page, browserName }) => {
    test.skip(browserName !== "chromium", "Clipboard read requires Chromium permissions");
    await mockCore(page);
    await seedState(page);
    await page.goto("/series/5");

    const copyBtn = page.getByRole("button", { name: /copy series figure to clipboard/i });
    await expect(copyBtn).toBeEnabled();
    await copyBtn.click();

    await page.waitForTimeout(500);
    const types = await page.evaluate(async () => {
      const items = await navigator.clipboard.read();
      const out: string[] = [];
      for (const item of items) for (const t of item.types) out.push(t);
      return out;
    });
    expect(types).toContain("image/png");
  });

  test("collapses the rail to full-bleed and restores it", async ({ page }) => {
    await mockCore(page);
    await seedState(page);
    await page.goto("/series/5");
    await page.getByTestId("rail-collapse-toggle").click();
    await expect(page.getByTestId("rail-restore")).toBeVisible();
    await page.getByTestId("rail-restore").click();
    await expect(page.getByTestId("series-builder-rail")).toBeVisible();
  });
});
