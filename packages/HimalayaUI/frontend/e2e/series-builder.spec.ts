import { test, expect, type Page } from "@playwright/test";

// Greenfield series-builder spec (SeriesBuilderPage + SeriesPlate + BuilderRail +
// MemberList). The legacy DOM (recipe-title / recipe-add-sample / recipe-save /
// recipe-commit / recipe-cancel + the commit-409 conflict modal) is GONE.
//
// The greenfield page is a LAZY-DRAFT surface on a single always-"Compose"
// screen:
//   • Opening /series/:id renders the committed series READ-ONLY. The rail's
//     "Confirm series" button is present but DISABLED (no live draft → no
//     onConfirm), and "Adjust" is the entry point.
//   • The FIRST recipe edit (title edit OR add-sample) silently STARTS a draft
//     → "Confirm series" becomes ENABLED + a Cancel button appears.
//   • "Confirm series" = a Save→Commit CHAIN: PATCH /api/series/:id (save,
//     persists the RECIPE) THEN, on its success, POST /api/series/:id/commit
//     with the plate RESOLVED FROM THE SAVED RECIPE (recipe samples → picker
//     indexing exposures; the PATCH does not rebuild members — BU-RECIPENOOP).
//     After commit success the page returns to read state (no navigation).
//   • There is NO separate Save button and NO conflict modal (409 relaxed to
//     last-write-wins — Plan 6a).
//   • Cancel discards the draft with no request.
//
// Selectors are taken verbatim from the component test
// (test/print-pages/SeriesBuilderPage.test.tsx): series-plate, series-member-row,
// builder-add-sample-select, aria-label="Series title", role=button names
// "Confirm series" / "Adjust" / "Cancel".

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
  state: "committed",
  members: [MEMBER],
  // The recipe (samples) backs the add-sample "addable" projection. One member
  // already in the recipe → only the OTHER corpus sample is addable.
  samples: [{ id: 1, series_id: 5, sample_id: 10, position: 0, pinned: false, excluded: false }],
};

// The series returned by a successful save: the RECIPE now carries BOTH
// samples (10 and the added 20) while the cached members still hold only one
// row — the dev-DB BU-RECIPENOOP shape. The commit body must be resolved from
// the recipe (samples 10, 20 → exposures 101, 202), not echo these members.
const SAVED_MEMBER = { ...MEMBER, id: 2, exposure_id: 202, display_order: 0 };
const SAVED_SERIES = {
  ...SERIES, content_hash: "h2", title: "LL37 titration v2",
  members: [SAVED_MEMBER],
  samples: [
    { id: 1, series_id: 5, sample_id: 10, position: 0, pinned: false, excluded: false },
    { id: 2, series_id: 5, sample_id: 20, position: 1, pinned: false, excluded: false },
  ],
};

const TRACE = {
  q: Array.from({ length: 50 }, (_, i) => 0.05 + i * 0.02),
  I: Array.from({ length: 50 }, (_, i) => 1000 / (1 + i)),
  sigma: Array.from({ length: 50 }, () => 5),
};

const CORPUS = [
  { id: 10, experiment_id: 1, name: "JC042", display_name: "DOPE 80%",
    notes: null, tags: [], q_units: "A-1" },
  { id: 20, experiment_id: 1, name: "JC050", display_name: "DOPE 90%",
    notes: null, tags: [], q_units: "A-1" },
];

async function mockCore(page: Page): Promise<void> {
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([{ id: 1, username: "jc" }]) }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  // Full nested series detail (GET). PATCH is overridden per-test where the
  // save chain is exercised.
  await page.route("**/api/series/5", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(SERIES) }));
  // Batch member traces (exposure_id → Trace). Both the committed (101) and the
  // post-save (202) plate exposures resolve so neither waterfall is empty.
  await page.route("**/api/series/5/traces", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify({ 101: TRACE, 202: TRACE }) }));
  // Corpus sample list — the add-sample select's option source.
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(CORPUS) }));
  // Corpus picker projection — the Confirm chain's sample→exposure resolution
  // source (BU-RECIPENOOP). Without it Confirm stays gated (disabled).
  await page.route("**/api/picker-samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([
        { sample: CORPUS[0], indexing_exposure_id: 101, all_exposures: [] },
        { sample: CORPUS[1], indexing_exposure_id: 202, all_exposures: [] },
      ]) }));
  // SSE endpoint — drain immediately so the multiplayer EventSource doesn't
  // hang the page (see e2e/AGENTS.md + the sibling series specs).
  await page.route("**/api/events", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
}

/** Seed a username + tutorialSeen so the first-run onboarding overlay (which
 *  intercepts pointer events) does not block clicks. Mirrors the sibling specs. */
async function seedState(page: Page): Promise<void> {
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state: { username: "jc", tutorialSeen: true, theme: "dark" }, version: 3 }),
    );
  });
}

test.describe("series builder — greenfield DOM", () => {
  test("opening /series/:id renders the committed series read-only; Confirm series is DISABLED", async ({ page }) => {
    await mockCore(page);
    await seedState(page);
    await page.goto("/series/5");

    // The committed plate renders.
    await expect(page.getByTestId("series-plate")).toBeVisible();
    // The committed title surfaces in the editable plate-title field.
    await expect(page.getByLabel(/series title/i)).toHaveValue("LL37 titration");
    // Read state: "Adjust" is the live entry point; no Cancel (no draft).
    await expect(page.getByRole("button", { name: /adjust/i })).toBeVisible();
    await expect(page.getByRole("button", { name: /^cancel$/i })).toHaveCount(0);
    // controls-don't-lie: "Confirm series" is present but DISABLED in read state.
    await expect(page.getByRole("button", { name: /confirm series/i })).toBeDisabled();
    // No separate Save button and no conflict modal in the read-state DOM.
    await expect(page.getByRole("button", { name: /^save$/i })).toHaveCount(0);
    await expect(page.getByTestId("conflict-modal")).toHaveCount(0);
  });

  test("adding a sample STARTS a draft → Confirm series becomes ENABLED + Cancel appears", async ({ page }) => {
    await mockCore(page);
    await seedState(page);
    await page.goto("/series/5");

    await expect(page.getByRole("button", { name: /confirm series/i })).toBeDisabled();

    // Add the only addable corpus sample (id 20 — sample 10 is already in the recipe).
    await page.getByTestId("builder-add-sample-select").selectOption("20");

    // The draft is now live: Confirm enables, Cancel appears, an editable recipe row renders.
    await expect(page.getByRole("button", { name: /confirm series/i })).toBeEnabled();
    await expect(page.getByRole("button", { name: /^cancel$/i })).toBeVisible();
    await expect(page.getByTestId("builder-recipe-row")).toHaveCount(2);
  });

  test("Confirm series fires PATCH then POST /commit (in that order) and returns to read state", async ({ page }) => {
    await mockCore(page);
    await seedState(page);

    // Record the request order. The commit (POST) mock asserts the PATCH was
    // seen first — gating, not timing, makes the ordering assertion robust.
    const order: string[] = [];
    let patchSeen = false;

    // Override GET/PATCH on /api/series/5. PATCH returns the saved full Series:
    // its RECIPE holds samples 10 + 20 while its members hold only exposure 202
    // — the commit body must be resolved from the recipe, not echo the members.
    await page.route("**/api/series/5", async (route) => {
      const req = route.request();
      if (req.method() === "PATCH") {
        order.push("PATCH");
        patchSeen = true;
        await route.fulfill({ status: 200, contentType: "application/json",
          body: JSON.stringify(SAVED_SERIES) });
        return;
      }
      await route.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify(SERIES) });
    });

    let commitBody: { members?: Array<{ exposure_id: number }> } | null = null;
    await page.route("**/api/series/5/commit", async (route) => {
      order.push("POST");
      commitBody = route.request().postDataJSON();
      await route.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify({ ...SAVED_SERIES, content_hash: "h3" }) });
    });

    await page.goto("/series/5");

    // Start a draft via a title edit, then Confirm.
    await page.getByLabel(/series title/i).fill("LL37 titration v2");
    const confirm = page.getByRole("button", { name: /confirm series/i });
    await expect(confirm).toBeEnabled();
    await confirm.click();

    // The PATCH lands first.
    await expect.poll(() => patchSeen, { timeout: 3000 }).toBe(true);
    // Then the POST /commit lands.
    await expect.poll(() => commitBody, { timeout: 3000 }).not.toBeNull();
    // ORDER assertion (load-bearing): PATCH strictly precedes POST.
    expect(order).toEqual(["PATCH", "POST"]);
    // Provenance (BU-RECIPENOOP): the commit body is the plate RESOLVED FROM
    // THE SAVED RECIPE — sample 10 → exposure 101 AND the added sample 20 →
    // exposure 202 — not an echo of the save response's stale single-member
    // plate ([202]) and not the old committed plate ([101]).
    expect(commitBody!.members!.map((m) => m.exposure_id)).toEqual([101, 202]);

    // After commit success the page returns to read state: Confirm DISABLED, no Cancel.
    await expect(page.getByRole("button", { name: /confirm series/i })).toBeDisabled();
    await expect(page.getByRole("button", { name: /^cancel$/i })).toHaveCount(0);
  });

  test("Cancel after an edit discards the draft — no PATCH, no POST", async ({ page }) => {
    await mockCore(page);
    await seedState(page);

    let patched = false;
    let committed = false;
    await page.route("**/api/series/5", async (route) => {
      const req = route.request();
      if (req.method() === "PATCH") {
        patched = true;
        await route.fulfill({ status: 200, contentType: "application/json",
          body: JSON.stringify(SAVED_SERIES) });
        return;
      }
      await route.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify(SERIES) });
    });
    await page.route("**/api/series/5/commit", async (route) => {
      committed = true;
      await route.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(SERIES) });
    });

    await page.goto("/series/5");

    // Start a draft, then cancel it.
    await page.getByLabel(/series title/i).fill("scratch edit");
    await expect(page.getByRole("button", { name: /confirm series/i })).toBeEnabled();
    await page.getByRole("button", { name: /^cancel$/i }).click();

    // Back to read state, with no request fired.
    await expect(page.getByRole("button", { name: /confirm series/i })).toBeDisabled();
    await expect(page.getByRole("button", { name: /^cancel$/i })).toHaveCount(0);
    expect(patched).toBe(false);
    expect(committed).toBe(false);
  });

  test("there is no separate Save button and no conflict modal anywhere in the flow", async ({ page }) => {
    await mockCore(page);
    await seedState(page);
    await page.goto("/series/5");

    // Read state.
    await expect(page.getByRole("button", { name: /^save$/i })).toHaveCount(0);
    await expect(page.getByTestId("conflict-modal")).toHaveCount(0);

    // Draft state — still no Save button, still no conflict modal.
    await page.getByLabel(/series title/i).fill("x");
    await expect(page.getByRole("button", { name: /confirm series/i })).toBeEnabled();
    await expect(page.getByRole("button", { name: /^save$/i })).toHaveCount(0);
    await expect(page.getByTestId("conflict-modal")).toHaveCount(0);
  });
});
