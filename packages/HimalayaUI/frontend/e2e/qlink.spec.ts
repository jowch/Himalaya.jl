import { test, expect, type Page } from "@playwright/test";

// q-link (#180): on the focus workspace (/sample/:id), the trace peaks and the
// detector rings share an ephemeral `hoveredQ`. Hovering a detector ring lights
// it (and cross-highlights the matching trace peak). The ring is the reliable
// Playwright target — trace-peak hover needs Observable Plot pixel math, which
// is exercised by the Vitest contract test; here we drive the ring source/sink
// through a real layout + the rotation-aware overlay.

const EXPERIMENT = {
  id: 1, name: "SSRL May 2026", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};
const SAMPLE = {
  id: 10, experiment_id: 1, display_name: "D1", name: "cubic_run03",
  notes: null, tags: [], q_units: "A-1",
};
const EXPOSURES = [
  { id: 50, sample_id: 10, filename: "cubic_run03-001.dat", kind: "file",
    selected: true, status: "accepted", image_path: "/img/50.png",
    image_version: "v1", tags: [], sources: [],
    trace_hash: null, analysis_inputs_hash: null },
];
const TRACE = { q: [0.04, 0.1, 0.2, 0.3], I: [10, 80, 40, 20], sigma: [1, 1, 1, 1] };
const PEAKS = [
  { id: 1, exposure_id: 50, q: 0.045, intensity: 100, prominence: 50,
    sharpness: 1.2, source: "auto", excluded: false },
  { id: 2, exposure_id: 50, q: 0.103, intensity: 80, prominence: 40,
    sharpness: 1.1, source: "auto", excluded: false },
  { id: 3, exposure_id: 50, q: 0.206, intensity: 50, prominence: 30,
    sharpness: 1.0, source: "auto", excluded: false },
];

// 1x1 transparent PNG.
const PNG_1x1 = Buffer.from(
  "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg==",
  "base64",
);

async function mockFocus(page: Page): Promise<void> {
  const json = (body: unknown) => ({
    status: 200, contentType: "application/json", body: JSON.stringify(body),
  });
  await page.route("**/api/users", (r) =>
    r.fulfill(json([{ id: 1, username: "alice" }])));
  await page.route("**/api/experiments", (r) => r.fulfill(json([EXPERIMENT])));
  await page.route("**/api/experiments/1", (r) => r.fulfill(json(EXPERIMENT)));
  await page.route("**/api/experiments/1/samples", (r) => r.fulfill(json([SAMPLE])));
  // Corpus list — keyed by sample id alone; how the layout learns experiment_id.
  await page.route("**/api/samples", (r) => r.fulfill(json([SAMPLE])));
  await page.route("**/api/samples/10/exposures*", (r) => r.fulfill(json(EXPOSURES)));
  await page.route("**/api/samples/10/messages", (r) => r.fulfill(json([])));
  await page.route("**/api/exposures/50/trace", (r) => r.fulfill(json(TRACE)));
  await page.route("**/api/exposures/50/peaks", (r) => r.fulfill(json(PEAKS)));
  await page.route("**/api/exposures/50/indices", (r) => r.fulfill(json([])));
  await page.route("**/api/exposures/50/groups", (r) => r.fulfill(json([])));
  await page.route("**/api/exposures/50/image*", (r) =>
    r.fulfill({ status: 200, contentType: "image/png", body: PNG_1x1 }));
}

async function seedUser(page: Page): Promise<void> {
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({
        state: { username: "alice", tutorialSeen: true, theme: "dark" },
        version: 3,
      }),
    );
  });
}

test.beforeEach(async ({ page }) => {
  await seedUser(page);
  await mockFocus(page);
});

test("hovering a detector ring lights it (q-link source + sink)", async ({ page }) => {
  await page.goto("/sample/10");
  await expect(page.getByTestId("focus-workspace-page")).toBeVisible();
  // Rings render once peaks load.
  const overlay = page.getByTestId("detector-ring-overlay");
  await expect(overlay).toBeVisible();
  // Drive the q-link from the hit-ring. SVG-in-overlay hover via the pointer
  // ray is brittle under the rotation transform / pointer-events layering, so
  // dispatch the same `mouseenter` the React handler binds — this exercises
  // the real source→hoveredQ→sink chain (the visible ring lights) without
  // depending on Playwright's SVG hit-testing.
  const hit = page.locator('[data-testid^="detector-ring-hit-"]').first();
  const ring = page.locator('[data-testid^="detector-ring-q-"]').first();
  await expect(ring).toHaveAttribute("data-hot", "false");
  // React 18 delegates at the root and synthesizes onMouseEnter from a
  // bubbling `mouseover` (a direct `mouseenter` dispatch is NOT seen). Use
  // mouseover/mouseout so the React handler actually runs.
  await hit.dispatchEvent("mouseover");
  await expect(ring).toHaveAttribute("data-hot", "true");
  await hit.dispatchEvent("mouseout");
  await expect(ring).toHaveAttribute("data-hot", "false");
});
