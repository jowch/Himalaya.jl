import { test, expect } from "./fixtures";

const SAMPLES = ["JC001", "JC002", "JC003", "BUF001"].map((name, i) => ({
  id: 40 + i, experiment_id: 17, display_name: name, name,
  notes: null, tags: [], q_units: "A-1",
}));

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.clear();
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: { username: "alice", tutorialSeen: true },
      version: 5,
    }));
  });
  const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });
  const experiment = { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "",
    manifest_path: null, created_at: "", q_units: null };
  await page.route("**/api/experiments", (r) => r.fulfill(json([experiment])));
  await page.route("**/api/experiments/17", (r) => r.fulfill(json(experiment)));
  await page.route("**/api/samples", (r) => r.fulfill(json(SAMPLES)));
  await page.route("**/api/users", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/bridges", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/labels", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/loads", (r) => r.fulfill(json([])));
  await page.route(/\/api\/samples\/\d+\/call$/, (r) => {
    const id = Number(/samples\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill(json({ sample_id: id, state: "indexed", phase: null, call_phases: [], members: [] }));
  });
  await page.route("**/api/samples/*/exposures*", (r) => r.fulfill(json([])));
  await page.route("**/api/events", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
});

test("Escape leaves the sample filter without clearing it", async ({ page }) => {
  await page.goto("/experiments/17/workspace/40");
  const rows = page.locator("[data-testid='workspace-samples'] li");
  await expect(rows).toHaveCount(4);
  await page.keyboard.press("/");
  await page.keyboard.type("BUF");
  await expect(rows).toHaveCount(1);
  await page.keyboard.press("Escape");
  await expect(page.getByRole("searchbox", { name: "Filter samples" })).toHaveValue("BUF");
  await expect(page.getByRole("searchbox", { name: "Filter samples" })).not.toBeFocused();
  await expect(rows).toHaveCount(1);
});

test("Enter in the sample filter opens the first match", async ({ page }) => {
  await page.goto("/experiments/17/workspace/40");
  await expect(page.locator("[data-testid='workspace-samples'] li")).toHaveCount(4);
  await page.keyboard.press("/");
  await page.keyboard.type("JC003");
  await page.keyboard.press("Enter");
  await expect(page).toHaveURL(/\/experiments\/17\/workspace\/42$/);
  await expect(page.locator("[data-testid='workspace-samples'] [aria-current='true']")).toHaveText(/JC003/);
});
