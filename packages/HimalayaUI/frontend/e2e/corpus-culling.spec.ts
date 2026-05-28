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
  // R2-M11 (#207): the per-thumb reject ✕ button is gone — single-thumb
  // reject is now "click-to-select, then X" (the same path the cull-bar's
  // Drop button exercises, but driven from the keyboard).
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
  await expect(page.getByTestId("sample-row-10")).toBeVisible();
  await page.getByTestId("exposure-thumb-1").click();
  await page.keyboard.press("x");

  await expect(page.getByTestId("exposure-thumb-1")).toHaveAttribute("data-rejected", "true");
  await expect.poll(() => patchedBody).toMatchObject({ status: "rejected" });
});

test("batch-reject: multi-select then reject PATCHes each selected exposure", async ({ page }) => {
  // R2-M11: no per-thumb checkbox; click the thumb body to select.
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
  await expect(page.getByTestId("sample-row-10")).toBeVisible();
  await page.getByTestId("exposure-thumb-1").click();
  await page.getByTestId("exposure-thumb-3").click();
  await page.getByTestId("batch-reject").click();

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
  // Open the loupe on exposure 2 by clicking its strip thumbnail.
  await page.getByTestId("thumb-cell-2").click();
  await page.getByTestId("loupe-set-representative").click();

  await expect.poll(() => selected).toBe(true);
});

test("loupe-flip: arrow keys move between exposures in the loupe", async ({ page }) => {
  await mockCorpus(page);
  await page.goto("/samples/loupe/10");
  await expect(page.getByTestId("loupe-page")).toBeVisible();
  // Loupe opens on the representative (exposure 1 → pos1.dat). Assert via the
  // sidebar's filename meta-row (data-testid, per e2e/AGENTS.md), not text.
  await expect(page.getByTestId("loupe-meta-filename")).toHaveText("pos1.dat");
  await page.keyboard.press("ArrowRight");
  await expect(page.getByTestId("loupe-meta-filename")).toHaveText("pos2.dat");
});
