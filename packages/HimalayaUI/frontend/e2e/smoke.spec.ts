import { test, expect, type Page } from "@playwright/test";

const EXPERIMENT = {
  id: 1, name: "SSRL May 2026", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};
const SAMPLES = [
  { id: 10, experiment_id: 1, label: "D1", name: "cubic_run03", notes: null, tags: [] },
  { id: 11, experiment_id: 1, label: "D2", name: "hex_run01",   notes: null, tags: [] },
];

async function mockCore(page: Page, users: { id: number; username: string }[] = []): Promise<void> {
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
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([EXPERIMENT]) }));
  await page.route("**/api/experiments/1", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(EXPERIMENT) }));
  await page.route("**/api/experiments/1/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(SAMPLES) }));
  for (const s of SAMPLES) {
    await page.route(`**/api/samples/${s.id}/exposures`, (r) =>
      r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
    await page.route(`**/api/samples/${s.id}/messages`, (r) =>
      r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  }
}

async function seedState(page: Page, extra: Record<string, unknown>): Promise<void> {
  await page.addInitScript((state) => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state, version: 2 }),
    );
  }, { username: "alice", activePage: "index", tutorialSeen: true, theme: "dark", ...extra });
}

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => { localStorage.clear(); });
});

test("first-run onboarding overlay is shown when no username", async ({ page }) => {
  await mockCore(page);
  await page.goto("/");
  await expect(page.getByTestId("onboarding-overlay")).toBeVisible();
  await expect(page.getByTestId("onboarding-name")).toBeVisible();
});

test("picking a new user triggers the tutorial and dismisses with 'Got it'", async ({ page }) => {
  await mockCore(page);
  await page.goto("/");

  await expect(page.getByTestId("onboarding-name")).toBeVisible();
  await page.locator('input[placeholder="Enter username"]').fill("alice");
  await page.getByTestId("onboarding-continue").click();

  await expect(page.getByTestId("onboarding-tutorial")).toBeVisible();
  await page.getByTestId("tutorial-next").click();
  await page.getByTestId("tutorial-next").click();
  await page.getByTestId("tutorial-next").click();
  await page.getByTestId("tutorial-done").click();

  await expect(page.getByTestId("onboarding-overlay")).not.toBeVisible();
  await expect(page.getByTestId("plot-title")).toBeVisible();
});

test("plot-title opens the nav modal at sample step when experiment is set", async ({ page }) => {
  await seedState(page, { activeExperimentId: 1, activeSampleId: 10 });
  await mockCore(page, [{ id: 1, username: "alice" }]);
  await page.goto("/");

  await expect(page.getByTestId("plot-title")).toBeVisible();
  await page.getByTestId("plot-title").click();
  await expect(page.getByTestId("nav-modal")).toBeVisible();
  await expect(page.getByTestId("nav-chip-experiment")).toHaveText(/SSRL May 2026/);
});

test("`.` key advances to the next sample", async ({ page }) => {
  await seedState(page, { activeExperimentId: 1, activeSampleId: 10 });
  await mockCore(page, [{ id: 1, username: "alice" }]);
  await page.goto("/");

  // Wait until the first sample name is visible
  await expect(page.getByTestId("plot-title")).toContainText("cubic_run03");
  await page.keyboard.press("."); // next
  await expect(page.getByTestId("plot-title")).toContainText("hex_run01");
  await page.keyboard.press(","); // back
  await expect(page.getByTestId("plot-title")).toContainText("cubic_run03");
});

test("chat posts a message via /api/samples/:id/messages", async ({ page }) => {
  let posted: { body: string } | null = null;
  await seedState(page, { activeExperimentId: 1, activeSampleId: 10 });
  await mockCore(page, [{ id: 1, username: "alice" }]);
  await page.route("**/api/samples/10/messages", async (route) => {
    const req = route.request();
    if (req.method() === "POST") {
      posted = req.postDataJSON() as { body: string };
      return route.fulfill({
        status: 201, contentType: "application/json",
        body: JSON.stringify({
          id: 1, sample_id: 10, author_id: 1, author: "alice", body: posted.body,
          created_at: "2026-04-24 10:00:00",
        }),
      });
    }
    return route.fulfill({ status: 200, contentType: "application/json", body: "[]" });
  });

  await page.goto("/");
  const compose = page.getByTestId("chat-compose");
  await compose.fill("looks cubic to me");
  await compose.press("Enter");
  await expect.poll(() => posted?.body, { timeout: 2000 }).toBe("looks cubic to me");
});

test("tab rocker switches to the Compare page placeholder", async ({ page }) => {
  await seedState(page, { activeExperimentId: 1, activeSampleId: 10 });
  await mockCore(page, [{ id: 1, username: "alice" }]);
  await page.goto("/");

  await expect(page.getByTestId("index-page")).toBeVisible();
  await page.getByTestId("tab-compare").click();
  await expect(page.getByTestId("compare-page")).toBeVisible();
  await expect(page.getByText(/Coming soon/i)).toBeVisible();
});
