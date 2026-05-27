import { test, expect, type Page } from "@playwright/test";

const SERIES_LISTING = [
  {
    id: 1, title: "Lipid dose response", description: null, content_hash: "h1",
    created_by: 1, created_at: "2026-05-01", updated_at: "2026-05-02",
    forked_from_id: null, forked_at_hash: null, view_grouping_mode: null,
    view_show_peak_ticks: null, view_show_peak_labels: null,
    last_event_at: "2026-05-02 10:00:00", author_username: "jc",
    member_count: 3, member_phases: ["Pn3m", "Lamellar"], member_phase_count: 2,
    has_stale_members: false,
  },
  {
    id: 2, title: "Cubic temperature sweep", description: null, content_hash: "",
    created_by: 1, created_at: "2026-05-03", updated_at: "2026-05-03",
    forked_from_id: null, forked_at_hash: null, view_grouping_mode: null,
    view_show_peak_ticks: null, view_show_peak_labels: null,
    last_event_at: "2026-05-03 09:00:00", author_username: "jc",
    member_count: 6, member_phases: ["Im3m"], member_phase_count: 1,
    has_stale_members: true,
  },
];

async function mockCore(page: Page): Promise<void> {
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([{ id: 1, username: "jc" }]) }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  // The folio's only data dependency:
  await page.route("**/api/series", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(SERIES_LISTING) }));
  // SSE endpoint — drain immediately with an empty stream so the multiplayer
  // EventSource doesn't hang the page (see e2e/AGENTS.md + compare.spec.ts).
  await page.route("**/api/events", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
}

test.describe("series folio", () => {
  test("renders the corpus masonry of saved series", async ({ page }) => {
    await mockCore(page);
    await page.goto("/series");
    await expect(page.getByTestId("series-folio-page")).toBeVisible();
    await expect(page.getByTestId("series-card-1")).toBeVisible();
    await expect(page.getByTestId("series-card-2")).toBeVisible();
    await expect(page.getByTestId("series-folio-count")).toHaveText("2");
    // draft card marked
    await expect(page.getByTestId("series-card-2")).toHaveAttribute("data-draft", "true");
    // stale card marked
    await expect(page.getByTestId("series-card-2")).toHaveAttribute("data-stale", "true");
  });

  test("filters by the search box", async ({ page }) => {
    await mockCore(page);
    await page.goto("/series");
    await page.getByTestId("series-folio-search").fill("cubic");
    await expect(page.getByTestId("series-card-1")).toHaveCount(0);
    await expect(page.getByTestId("series-card-2")).toBeVisible();
  });

  test("the Series stage-tab is active on /series", async ({ page }) => {
    await mockCore(page);
    await page.goto("/series");
    await expect(page.getByTestId("stage-tab-series")).toHaveAttribute("aria-current", "page");
  });
});
