import { test, expect, type Page } from "@playwright/test";

// Spec §8.2 — Playwright mocked. /api/* is intercepted so this runs without a
// backend.
//
// I4.4 (#181): the Index surface is retired. Old `/index/<exp>/<sample>`
// permalink deep-links are PRESERVED but now redirect to the focus workspace
// (`/sample/:id`) via IndexSlugRedirect, which resolves the slug pair through
// `/api/resolve`. Sampleless `/index*` and bare `/` redirect to the corpus
// contact sheet (`/samples`). A failed resolve also lands on `/samples`.

const RESOLVE_RE = /\/api\/resolve\?/;

const SAMPLE = {
  id: 42, experiment_id: 17, display_name: "D1", name: "JC001",
  notes: null, tags: [], q_units: "A-1",
};

async function seedSession(page: Page): Promise<void> {
  await page.addInitScript(() => {
    localStorage.clear();
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: {
        username: "alice", firstName: undefined, lastName: undefined,
        tutorialSeen: true, theme: "dark",
        activeExperimentId: undefined,
        activeSampleId: undefined,
        activeExposureId: undefined,
      },
      version: 3,
    }));
  });
}

test.beforeEach(async ({ page }) => {
  await page.route("**/api/experiments", (route) => {
    route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify([
        { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "",
          manifest_path: null, created_at: "", q_units: null },
      ]),
    });
  });
  await page.route("**/api/experiments/17", (route) =>
    route.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify({ id: 17, name: "lipid", path: "", data_dir: "",
        analysis_dir: "", manifest_path: null, created_at: "", q_units: null }) }));
  await page.route("**/api/experiments/17/samples", (route) => {
    route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify([SAMPLE]),
    });
  });
  // Corpus list — the focus workspace learns the sample's experiment from here.
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([SAMPLE]) }));
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/samples/42/exposures*", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/samples/42/messages", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
});

test("legacy /index/<exp>/<sample> deep-link resolves the slug and redirects to the focus workspace", async ({ page }) => {
  await seedSession(page);
  await page.route(RESOLVE_RE, (route) => {
    route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    });
  });
  await page.goto("/index/lipid/JC001");
  // IndexSlugRedirect resolves the slug → /sample/42 (the focus workspace).
  await expect(page).toHaveURL(/\/sample\/42$/);
  await expect(page.getByTestId("focus-workspace")).toBeVisible();
});

test("legacy /index deep-link whose slug 404s falls back to the corpus contact sheet", async ({ page }) => {
  await seedSession(page);
  await page.route(RESOLVE_RE, (route) => {
    route.fulfill({
      status: 404, contentType: "application/json",
      body: JSON.stringify({
        error: "not_found", missing: "experiment", missing_value: "lipid-typo",
        experiment_resolved: undefined, sample_resolved: undefined,
      }),
    });
  });
  await page.goto("/index/lipid-typo/JC001");
  await expect(page).toHaveURL(/\/samples$/);
  await expect(page.getByTestId("samples-page")).toBeVisible();
});

test("sampleless /index redirects to the corpus contact sheet", async ({ page }) => {
  await seedSession(page);
  await page.goto("/index");
  await expect(page).toHaveURL(/\/samples$/);
  await expect(page.getByTestId("samples-page")).toBeVisible();
});

test("bare / cold-mount redirects to the corpus contact sheet (§4.1)", async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: {
        activeExperimentId: 17,
        activeSampleId: 42,
        username: "test", firstName: undefined, lastName: undefined,
        activeExposureId: undefined,
        tutorialSeen: true, theme: "dark",
      },
      version: 3,
    }));
  });
  await page.goto("/");
  await expect(page).toHaveURL(/\/samples$/);
  await expect(page.getByTestId("samples-page")).toBeVisible();
});
