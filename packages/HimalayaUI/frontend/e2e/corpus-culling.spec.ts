import { test, expect, type Page } from "@playwright/test";

// I1.7 (#163): the Phase-1 mocked spec for the corpus culling flow, ceded here
// by I1.6. Exercises the contact sheet (cull / batch-reject / representative)
// and the loupe (flip), all against mocked /api routes. Replaces the deleted
// Inspect E2E (inspect.spec.ts).

const EXPERIMENT = {
  id: 1, name: "SSRL Test", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};
const SAMPLES = [
  {
    id: 10, experiment_id: 1, display_name: "D1", name: "run03",
    notes: null, tags: [], q_units: "A-1",
  },
];
const EXPOSURES = [
  {
    id: 1, sample_id: 10, filename: "pos1.dat", kind: "file", selected: true,
    status: "accepted", image_path: null, image_version: "", tags: [],
    sources: [], trace_hash: null, analysis_inputs_hash: null,
  },
  {
    id: 2, sample_id: 10, filename: "pos2.dat", kind: "file", selected: false,
    status: "accepted", image_path: null, image_version: "", tags: [],
    sources: [], trace_hash: null, analysis_inputs_hash: null,
  },
  {
    id: 3, sample_id: 10, filename: "pos3.dat", kind: "file", selected: false,
    status: "accepted", image_path: null, image_version: "", tags: [],
    sources: [], trace_hash: null, analysis_inputs_hash: null,
  },
];

async function mockCorpus(page: Page): Promise<void> {
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPERIMENT]) }));
  await page.route("**/api/experiments/1", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPERIMENT) }));
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(SAMPLES) }));
  await page.route("**/api/samples/10/exposures*", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPOSURES) }));
  await page.route("**/api/samples/10/messages", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
}

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() =>
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({
        state: { username: "alice", theme: "dark", tutorialSeen: true },
        version: 3,
      }),
    ),
  );
});

test("cull: rejecting a selected exposure via X dims its thumbnail and PATCHes status", async ({ page }) => {
  // Greenfield contact sheet: single-thumb reject is "click-to-select, then X".
  // The thumb gains data-state="selected"; the page-global CullBar appears; the
  // X key drops the selection (the same path the CullBar "Drop" button drives).
  // After the PATCH resolves + the exposures query updates, the optimistic
  // status flips to "rejected" so the thumb re-renders dimmed.
  let patchedBody: unknown = null;
  await mockCorpus(page);
  await page.route("**/api/exposures/1/status", async (route) => {
    patchedBody = route.request().postDataJSON();
    await route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({ id: 1, status: "rejected" }),
    });
  });

  await page.goto("/samples");
  const row = page.getByTestId("sample-table-row").first();
  await expect(row).toBeVisible();

  // Mock EXPOSURES order is ids 1,2,3 → nth(0) is exposure 1.
  const thumb1 = row.getByTestId("thumbnail").nth(0);
  await thumb1.click();
  await expect(thumb1).toHaveAttribute("data-state", /selected/);
  await expect(page.getByTestId("cull-bar")).toHaveAttribute("data-show", "true");

  await page.keyboard.press("x");

  // The dropped thumb re-renders rejected (state token + dimmed image).
  await expect(thumb1).toHaveAttribute("data-state", /rejected/);
  await expect(thumb1.locator('[data-dimmed="true"]')).toBeVisible();
  await expect.poll(() => patchedBody).toMatchObject({ status: "rejected" });
  // Drop clears the selection → the cull bar hides.
  await expect(page.getByTestId("cull-bar")).toHaveAttribute("data-show", "false");
});

test("batch-reject: multi-select then reject PATCHes each selected exposure", async ({ page }) => {
  // Click two thumb bodies to add them to the page-global selection, then hit
  // the CullBar "Drop" button. Each selected exposure fires a status PATCH.
  const patched: string[] = [];
  await mockCorpus(page);
  await page.route(/\/api\/exposures\/\d+\/status$/, async (route) => {
    patched.push(new URL(route.request().url()).pathname);
    await route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({ status: "rejected" }),
    });
  });

  await page.goto("/samples");
  const row = page.getByTestId("sample-table-row").first();
  await expect(row).toBeVisible();

  // EXPOSURES order ids 1,2,3 → nth 0,2 are exposures 1 and 3.
  await row.getByTestId("thumbnail").nth(0).click();
  await row.getByTestId("thumbnail").nth(2).click();

  const cullBar = page.getByTestId("cull-bar");
  await expect(cullBar).toHaveAttribute("data-show", "true");
  await cullBar.getByRole("button", { name: /Drop/ }).click();

  await expect.poll(() => patched.length).toBe(2);
  expect(patched.some((p) => p.endsWith("/exposures/1/status"))).toBe(true);
  expect(patched.some((p) => p.endsWith("/exposures/3/status"))).toBe(true);
  // Exposure 2 was not selected → no op.
  expect(patched.some((p) => p.endsWith("/exposures/2/status"))).toBe(false);
});

test("representative: picking a representative in the loupe PATCHes select", async ({ page }) => {
  // R2-M11 (#207): the per-thumb rep ⊙ button is gone — representative pick
  // lives in the loupe sidebar's "Set as representative" affordance (and the
  // `R` keystroke). Cover the loupe path here so the contact-sheet spec keeps
  // exercising the underlying selectExposure mutator.
  let selected = false;
  await mockCorpus(page);
  await page.route("**/api/exposures/2/select", async (route) => {
    selected = true;
    await route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({ id: 2, selected: true }),
    });
  });

  await page.goto("/samples/loupe/10");
  await expect(page.getByTestId("loupe-page")).toBeVisible();
  // Open the loupe on exposure 2 by clicking its strip thumbnail (frameNo 2 →
  // nth(1)). Every greenfield strip thumb shares data-testid="thumbnail".
  await page.getByTestId("thumbnail").nth(1).click();
  await page.getByRole("button", { name: "Set as representative" }).click();

  await expect.poll(() => selected).toBe(true);
});

test("loupe-flip: arrow keys move between exposures in the loupe", async ({ page }) => {
  await mockCorpus(page);
  await page.goto("/samples/loupe/10");
  await expect(page.getByTestId("loupe-page")).toBeVisible();
  // Loupe opens on the representative (exposure 1, frame 1). The greenfield
  // page has no filename row; assert the active frame via the BigFrame caption.
  await expect(page.locator('[data-role="frame-caption"]')).toContainText("frame 1 of");
  await page.keyboard.press("ArrowRight");
  await expect(page.locator('[data-role="frame-caption"]')).toContainText("frame 2 of");
});

// Regression: a many-exposure loupe must NOT let the filmstrip balloon the grid's
// `1fr` column and shove the side panel off-screen. (Found in the Phase-4 loupe
// walkthrough: a CSS `1fr` track has an implicit `min-width:auto`, so the 20-thumb
// strip forced the column to ~1552px and pushed the side panel past the viewport.
// Fixed via `minmax(0,1fr)` + `min-w-0` on the column + `justify-center-safe` on
// the strip.) The mockup only ever drew a handful of thumbs, so it never exposed this.
test("loupe layout: a many-exposure filmstrip keeps the side panel on-screen", async ({ page }) => {
  await page.setViewportSize({ width: 1200, height: 800 });
  const MANY = Array.from({ length: 16 }, (_, i) => ({
    id: 500 + i, sample_id: 11, filename: `f${i + 1}.dat`, kind: "file",
    selected: i === 0, status: i === 0 ? "accepted" : null,
    image_path: null, image_version: "", tags: [],
    sources: [], trace_hash: null, analysis_inputs_hash: null,
  }));
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPERIMENT]) }));
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([
      { id: 11, experiment_id: 1, display_name: "D-many", name: "run11", notes: null, tags: [], q_units: "A-1" },
    ]) }));
  await page.route("**/api/samples/11/exposures*", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(MANY) }));
  await page.route("**/api/samples/11/messages", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));

  await page.goto("/samples/loupe/11");
  await expect(page.getByTestId("loupe-page")).toBeVisible();

  // The filmstrip genuinely overflows its column (so the guard is meaningful)…
  const strip = page.getByTestId("thumbnail-gallery");
  expect(await strip.evaluate((g) => g.scrollWidth > g.clientWidth + 1)).toBe(true);

  // …yet the side panel's right edge stays within the viewport.
  const box = await page.getByTestId("loupe-side-panel").boundingBox();
  expect(box).not.toBeNull();
  expect(box!.x + box!.width).toBeLessThanOrEqual(1200 + 1);
});
