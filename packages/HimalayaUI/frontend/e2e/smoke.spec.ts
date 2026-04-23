import { test, expect, type Page } from "@playwright/test";

const EXP = {
  id: 1, name: "TestRun", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-01-01",
};
const SAMPLES = [
  { id: 1, experiment_id: 1, label: "D1", name: "UX1", notes: null,
    tags: [{ id: 1, key: "lipid", value: "DOPC", source: "manual" }] },
  { id: 2, experiment_id: 1, label: "D2", name: "UX2", notes: null, tags: [] },
];

async function mockApi(page: Page, users: { id: number; username: string }[] = []): Promise<void> {
  await page.route("**/api/users", async (route) => {
    const req = route.request();
    if (req.method() === "GET") {
      await route.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify(users) });
    } else if (req.method() === "POST") {
      const data = JSON.parse(req.postData() ?? "{}");
      await route.fulfill({ status: 201, contentType: "application/json",
        body: JSON.stringify({ id: users.length + 1, username: data.username }) });
    } else {
      await route.continue();
    }
  });
  await page.route("**/api/experiments/1", async (route) =>
    route.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXP) }));
  await page.route("**/api/experiments/1/samples", async (route) =>
    route.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(SAMPLES) }));
}

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => { localStorage.clear(); });
});

test("shell loads with navbar, layout placeholders, and sample list", async ({ page }) => {
  await mockApi(page);
  await page.goto("/");

  await expect(page.locator("[data-testid='nav-logo']")).toHaveText("Himalaya");
  await expect(page.locator("[data-testid='placeholder']")).toHaveCount(4);
  await expect(page.locator("[data-sample-id]")).toHaveCount(2);
  await expect(page.locator('[data-sample-id="1"] [data-testid="sample-label"]')).toHaveText("D1");
});

test("first-visit prompts user modal and submits new user", async ({ page }) => {
  await mockApi(page);
  await page.goto("/");

  await expect(page.locator("[role='dialog']")).toBeVisible();
  await expect(page.locator("[role='dialog'] h2")).toHaveText("Who are you?");

  await page.locator("input[placeholder='Enter username']").fill("alice");
  await page.locator("[role='dialog'] button", { hasText: "Continue" }).click();

  await expect(page.locator("[role='dialog']")).not.toBeVisible();
  await expect(page.locator("[data-testid='nav-user']")).toHaveText("alice");
});

test("clicking a sample updates the active state and breadcrumb", async ({ page }) => {
  await mockApi(page, [{ id: 1, username: "alice" }]);
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state: { username: "alice" }, version: 1 }),
    );
  });
  await page.goto("/");

  await expect(page.locator("[role='dialog']")).not.toBeVisible();
  await page.locator('[data-sample-id="2"]').click();
  await expect(page.locator('[data-sample-id="2"]')).toHaveAttribute("data-active", "true");
  await expect(page.locator("[data-testid='nav-breadcrumb']")).toContainText("D2");
});

test("filter narrows the sample list", async ({ page }) => {
  await mockApi(page, [{ id: 1, username: "alice" }]);
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state: { username: "alice" }, version: 1 }),
    );
  });
  await page.goto("/");

  await page.locator("input[placeholder='Filter samples\u2026']").fill("DOPC");
  await expect(page.locator("[data-sample-id]")).toHaveCount(1);
  await expect(page.locator("[data-sample-id]")).toHaveAttribute("data-sample-id", "1");
});
