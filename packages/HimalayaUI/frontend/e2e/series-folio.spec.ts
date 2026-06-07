import { test, expect, type Page } from "@playwright/test";

const SERIES_LISTING = [
  {
    id: 1, title: "Lipid dose response", description: null, content_hash: "h1",
    created_by: 1, created_at: "2026-05-01", updated_at: "2026-05-02",
    forked_from_id: null, forked_at_hash: null, view_grouping_mode: null,
    view_show_peak_ticks: null, view_show_peak_labels: null,
    last_event_at: "2026-05-02 10:00:00", author_username: "jc",
    member_count: 3, member_phases: ["Pn3m", "Lamellar"], member_phase_count: 2,
    has_stale_members: false, ordering_variable: null, spans_experiments: false,
  },
  {
    id: 2, title: "Cubic temperature sweep", description: null, content_hash: "",
    created_by: 1, created_at: "2026-05-03", updated_at: "2026-05-03",
    forked_from_id: null, forked_at_hash: null, view_grouping_mode: null,
    view_show_peak_ticks: null, view_show_peak_labels: null,
    last_event_at: "2026-05-03 09:00:00", author_username: "jc",
    member_count: 6, member_phases: ["Im3m"], member_phase_count: 1,
    has_stale_members: true, ordering_variable: null, spans_experiments: false,
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
  // Per-series detail + traces fetches (FolioCard mounts a useSeries + useSeriesTraces
  // per card; mock both so the cards can render).
  for (const s of SERIES_LISTING) {
    await page.route(`**/api/series/${s.id}`, (r) =>
      r.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify({ ...s, members: [] }) }));
    await page.route(`**/api/series/${s.id}/traces`, (r) =>
      r.fulfill({ status: 200, contentType: "application/json", body: "{}" }));
  }
  // SSE endpoint — drain immediately with an empty stream so the multiplayer
  // EventSource doesn't hang the page (see e2e/AGENTS.md + compare.spec.ts).
  await page.route("**/api/events", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
}

test.describe("series folio", () => {
  // Greenfield folio testid mapping (replaces legacy series-folio-page DOM):
  //   series-folio-page  → folio-header     (FolioHeader root, data-testid="folio-header")
  //   series-card-N      → series-card      (all SeriesCards share one testid; use .nth())
  //   series-folio-count → folio-count      (count div inside FolioHeader)
  //   series-folio-search → search-input    (SearchInput wrapper; <input> reached via getByRole)
  //   data-stale         → DROPPED          (greenfield SeriesCard has no data-stale attribute;
  //                                          has_stale_members is not surfaced as a card attribute
  //                                          in the greenfield design — the stale-indices path
  //                                          belongs to the Focus page, not the folio card)

  test("renders the corpus masonry of saved series", async ({ page }) => {
    await mockCore(page);
    await page.goto("/series");
    // Folio header is the greenfield folio-mounted signal (replaces series-folio-page).
    await expect(page.getByTestId("folio-header")).toBeVisible();
    // Both cards render (two series in the mock listing).
    await expect(page.getByTestId("series-card")).toHaveCount(2);
    // FolioHeader count shows the total number of series.
    await expect(page.getByTestId("folio-count")).toHaveText("2");
    // Series 2 has content_hash="" → draft=true → SeriesCard carries data-draft="true".
    await expect(page.getByTestId("series-card").nth(1)).toHaveAttribute("data-draft", "true");
    // Series 1 is a saved series (non-empty content_hash) → not a draft.
    await expect(page.getByTestId("series-card").nth(0)).toHaveAttribute("data-draft", "false");
  });

  test("filters by the search box", async ({ page }) => {
    await mockCore(page);
    await page.goto("/series");
    // SearchInput wrapper has data-testid="search-input"; the actual <input> is inside.
    await page.getByTestId("search-input").getByRole("textbox").fill("cubic");
    // "Lipid dose response" (series 1) should be filtered out.
    await expect(page.getByTestId("series-card")).toHaveCount(1);
    // The surviving card is "Cubic temperature sweep".
    await expect(page.getByTestId("series-card").first()).toHaveAttribute("data-draft", "true");
  });

  test("the Series stage-tab is active on /series", async ({ page }) => {
    await mockCore(page);
    await page.goto("/series");
    await expect(page.getByTestId("stage-tab-series")).toHaveAttribute("aria-current", "page");
  });
});
