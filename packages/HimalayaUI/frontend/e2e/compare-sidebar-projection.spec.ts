/**
 * Compare-sidebar listing-projection e2e (Compare UX Phase F, F-4).
 *
 * Mocks `/api/experiments/1/comparisons` to return a row in the Phase A
 * listing-projection shape and asserts the redesigned sidebar row renders
 * the phase summary, author byline, and relative age.
 *
 * All `/api/*` endpoints are mocked via `page.route`; a catch-all returns
 * `[]` so the surrounding ComparePage shell mounts without console noise.
 */
import { test, expect, type Page, type Route } from "@playwright/test";

const EXPERIMENT = {
  id: 1, name: "SSRL Test", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};

const USERS = [{ id: 1, username: "alice", first_name: null, last_name: null }];

// One comparison in the new listing-projection shape (#136/#137).
const COMPARISON_ROW = {
  id: 7,
  title: "Cubic vs Hex",
  description: null,
  content_hash: "h7",
  created_by: 1, // alice's id → current user → "by you"
  created_at: "2026-05-01T00:00:00Z",
  updated_at: "2026-05-10T00:00:00Z",
  forked_from_id: null,
  forked_at_hash: null,
  view_grouping_mode: null,
  view_show_peak_ticks: null,
  view_show_peak_labels: null,
  last_event_at: "2026-05-10T00:00:00Z",
  author_username: "alice",
  member_count: 3,
  member_phases: ["Pn3m", "Hex"],
  member_phase_count: 2,
  has_stale_members: false,
};

async function jsonOK(route: Route, body: unknown): Promise<void> {
  await route.fulfill({
    status: 200,
    contentType: "application/json",
    body: JSON.stringify(body),
  });
}

async function mockApi(page: Page): Promise<void> {
  // Catch-all first → specific routes registered after it take precedence.
  await page.route("**/api/**", (r) => jsonOK(r, []));
  await page.route("**/api/users", (r) => jsonOK(r, USERS));
  await page.route("**/api/experiments", (r) => jsonOK(r, [EXPERIMENT]));
  await page.route("**/api/experiments/1", (r) => jsonOK(r, EXPERIMENT));
  await page.route("**/api/experiments/1/comparisons", (r) =>
    jsonOK(r, [COMPARISON_ROW]));
  await page.route("**/api/users/me/comparison-pins", (r) => jsonOK(r, []));
}

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.clear();
    sessionStorage.clear();
  });
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({
        state: {
          username: "alice",
          activePage: "compare",
          tutorialSeen: true,
          theme: "dark",
          activeExperimentId: 1,
        },
        version: 3,
      }),
    );
  });
});

test("sidebar shows phase summary and author + age", async ({ page }) => {
  await mockApi(page);
  await page.goto("/experiments/1/compare");

  const sidebar = page.getByTestId("comparison-sidebar");
  await expect(sidebar).toBeVisible();

  // Row title from the projection.
  await expect(sidebar.getByText("Cubic vs Hex")).toBeVisible();
  // Phase summary: member_phases joined · member_count traces.
  await expect(sidebar.getByText("Pn3m · Hex · 3 traces")).toBeVisible();
  // Author byline ("by you" — current user authored it) + relative age.
  await expect(sidebar.getByText(/by you · edited/)).toBeVisible();
});
