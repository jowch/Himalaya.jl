import { test, expect, type Page } from "./fixtures";

const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });

const experiment = (id: number, name: string) => ({
  id, name, description: null, path: `/lab/${id}`, data_dir: `/lab/${id}/data`, analysis_dir: `/lab/${id}/an`,
  manifest_path: null, created_at: "2026-04-12T00:00:00", q_units: "A^-1", beam_center_x: 1, beam_center_y: 1,
  pixel_size_um: 172, energy_kev: 9, flight_path_m: 1.8, energy_kev_source: "prp", flight_path_m_source: "setup",
  beam_center_x_source: "setup", beam_center_y_source: "setup", pixel_size_um_source: "prp", q_units_source: "prp",
  last_scanned_at: null, ingest_status: "complete", image_pattern: null, metadata_pattern: null, integration_pattern: null,
  stats: { loads: 2, samples: 13, exposures: 28, sessions: 1, span_hours: 2, started_at: "2026-04-01T10:00:00" },
});

async function mockShell(page: Page): Promise<void> {
  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({ state: { username: "alice", tutorialSeen: true }, version: 6 }));
  });
  await page.route("**/api/users", (r) => r.fulfill(json([])));
  await page.route("**/api/events*", (r) => r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
  await page.route("**/api/label-colors", (r) => r.fulfill(json([])));
  await page.route(/\/api\/fs\/suggest\?/, (r) => r.fulfill(json({ suggestions: [] })));
  await page.route(/\/api\/fs\/preview\?/, (r) => r.fulfill(json({ path: "", setups: [] })));
  await page.route((url) => url.pathname === "/api/experiments", (r) =>
    r.fulfill(json([experiment(1, "First"), experiment(2, "Second"), experiment(3, "Third")])));
}

const focusedCard = (page: Page) => page.evaluate(() => document.activeElement?.getAttribute("data-testid"));

test("arrow keys move between experiment cards in the grid, and Enter opens one", async ({ page }) => {
  await mockShell(page);
  await page.goto("/experiments");
  await expect(page.getByTestId("experiment-card-3")).toBeVisible();
  await page.keyboard.press("ArrowDown");
  expect(await focusedCard(page)).toBe("experiment-card-1");
  await page.keyboard.press("ArrowRight");
  expect(await focusedCard(page)).toBe("experiment-card-2");
  await page.keyboard.press("ArrowDown");
  expect(await focusedCard(page)).toBe("experiment-card-3");
  await page.keyboard.press("ArrowUp");
  expect(await focusedCard(page)).toBe("experiment-card-1");
  await page.route("**/api/experiments/1", (r) => r.fulfill(json(experiment(1, "First"))));
  await page.route("**/api/experiments/1/**", (r) => r.fulfill(json([])));
  await page.route((url) => url.pathname === "/api/samples", (r) => r.fulfill(json([])));
  await page.keyboard.press("Enter");
  await expect(page).toHaveURL(/\/experiments\/1/);
});

test("? lists the keys for the page it is pressed on", async ({ page }) => {
  await mockShell(page);
  await page.goto("/experiments");
  await expect(page.getByTestId("experiment-card-1")).toBeVisible();
  await page.keyboard.press("?");
  const help = page.getByRole("dialog", { name: "Keyboard shortcuts" });
  await expect(help).toContainText("Move between experiments");
  await page.keyboard.press("Escape");
  await expect(help).toHaveCount(0);
});

test("the new-experiment page takes typing straight into the folder field", async ({ page }) => {
  await mockShell(page);
  await page.goto("/experiments/new");
  const field = page.getByTestId("dirpicker-input").locator("input");
  await expect(field).toBeFocused();
  await page.keyboard.type("/lab");
  await expect(field).toHaveValue("/lab");
});

test("? on the new-experiment page opens the shortcuts until the folder field has text", async ({ page }) => {
  await mockShell(page);
  await page.goto("/experiments/new");
  const field = page.getByTestId("dirpicker-input").locator("input");
  await expect(field).toBeFocused();
  await page.keyboard.press("?");
  const help = page.getByRole("dialog", { name: "Keyboard shortcuts" });
  await expect(help).toBeVisible();
  await expect(field).toHaveValue("");
  await page.keyboard.press("Escape");
  await expect(help).toHaveCount(0);
  await field.focus();
  await page.keyboard.type("/lab?");
  await expect(field).toHaveValue("/lab?");
  await expect(help).toHaveCount(0);
});
