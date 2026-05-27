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
  // Corpus sample list — the recipe editor's "Add sample" dropdown source.
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([
      { ...SAMPLE },
      { id: 20, experiment_id: 1, name: "JC050", display_name: "DOPE 90%",
        notes: null, tags: [], q_units: "A-1" },
    ]) }));
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

  test("build → edit → commit → export round-trips", async ({ page }) => {
    await mockCore(page);
    await seedState(page);

    // Capture the recipe-save PATCH + the commit POST bodies.
    let patchBody: unknown = null;
    let commitBody: unknown = null;
    await page.route("**/api/series/5", async (route) => {
      const req = route.request();
      if (req.method() === "PATCH") {
        patchBody = req.postDataJSON();
        await route.fulfill({ status: 200, contentType: "application/json",
          body: JSON.stringify(SERIES) });
        return;
      }
      await route.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify(SERIES) });
    });
    await page.route("**/api/series/5/commit", async (route) => {
      commitBody = route.request().postDataJSON();
      await route.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify({ ...SERIES, content_hash: "h2" }) });
    });

    await page.goto("/series/5");
    await expect(page.getByTestId("series-builder-page")).toBeVisible();

    // Enter edit mode → the recipe editor appears.
    await page.getByTestId("series-builder-edit").click();
    await expect(page.getByTestId("series-recipe-editor")).toBeVisible();

    // Add a sample → an optimistic placeholder row appears immediately.
    // SERIES.samples is empty, so the added row is the recipe's first row.
    await page.getByTestId("recipe-add-sample").selectOption("20");
    await expect(page.getByTestId("recipe-row")).toHaveCount(1);

    // Save the recipe (optimistic → PATCH; NO expected_content_hash).
    await page.getByTestId("recipe-save").click();
    await expect.poll(() => patchBody).not.toBeNull();
    expect(JSON.stringify(patchBody)).not.toContain("expected_content_hash");

    // Commit the plate (spinner → POST commit, carries expected_content_hash).
    await page.getByTestId("recipe-commit").click();
    await expect.poll(() => commitBody).not.toBeNull();
    // Positively pin the commit-only-hash invariant end-to-end (round-1 nit).
    expect(JSON.stringify(commitBody)).toContain("expected_content_hash");

    // After commit success the editor tears down (back to read mode).
    await expect(page.getByTestId("series-builder-edit")).toBeVisible();

    // Export still works (figure-export round-trip) — the control is enabled.
    await expect(
      page.getByRole("button", { name: /copy series figure to clipboard/i }),
    ).toBeEnabled();
  });

  test("a zero-member series is editable — the recipe editor mounts and the first sample can be added", async ({ page }) => {
    await mockCore(page);
    await seedState(page);
    // A zero-member series — reachable from the folio (the backend listing
    // includes zero-member series). The empty-plate placeholder must not lock
    // out edit mode: the rail + recipe editor still mount so the first sample
    // can be added (round-1 blocking fix).
    const EMPTY = { ...SERIES, id: 7, members: [], samples: [] };
    await page.route("**/api/series/7", (r) =>
      r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EMPTY) }));

    await page.goto("/series/7");
    await expect(page.getByTestId("series-builder-page")).toBeVisible();
    await expect(page.getByTestId("series-builder-empty")).toBeVisible();

    // Edit mode is reachable for the empty series.
    await page.getByTestId("series-builder-edit").click();
    await expect(page.getByTestId("series-recipe-editor")).toBeVisible();

    // The empty placeholder persists in the plot area, but the editor lets the
    // user add the FIRST sample — the flow the gate used to lock out.
    await expect(page.getByTestId("series-builder-empty")).toBeVisible();
    await page.getByTestId("recipe-add-sample").selectOption("20");
    await expect(page.getByTestId("recipe-row")).toHaveCount(1);
  });

  test("commit 409 opens the series conflict modal", async ({ page }) => {
    await mockCore(page);
    await seedState(page);
    await page.route("**/api/series/5/commit", (route) =>
      route.fulfill({ status: 409, contentType: "application/json",
        body: JSON.stringify({ current_hash: "h2", current_state: { ...SERIES, content_hash: "h2" } }) }));

    await page.goto("/series/5");
    await page.getByTestId("series-builder-edit").click();
    await page.getByTestId("recipe-commit").click();

    // The series-commit conflict modal opens (shared chrome, series semantics).
    await expect(page.getByTestId("conflict-modal")).toBeVisible();
    await expect(page.getByText(/Series changed while you were editing/)).toBeVisible();
    // No Fork action on the series modal.
    await expect(page.getByTestId("conflict-fork")).toHaveCount(0);
  });

  test("permalink: /series/5 is a URL-owned deep link that re-renders on reload", async ({ page }) => {
    await mockCore(page);
    await seedState(page);
    await page.goto("/series/5");
    await expect(page.getByTestId("series-builder-page")).toBeVisible();
    // Reload the same URL — the page re-fetches via useSeries (react-router
    // owns the id; no slug round-trip, no stale flag).
    await page.reload();
    await expect(page.getByTestId("series-builder-page")).toBeVisible();
    await expect(page).toHaveURL(/\/series\/5$/);
  });
});
