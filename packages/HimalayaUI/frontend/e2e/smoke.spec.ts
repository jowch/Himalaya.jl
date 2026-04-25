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

test("curate: clicking + adds a candidate to the active set", async ({ page }) => {
  const EXPOSURE = {
    id: 5, sample_id: 10, filename: "scan1.dat", kind: "file",
    selected: true, tags: [], sources: [],
  };
  const CANDIDATES = [
    {
      id: 1, exposure_id: 5, phase: "Pn3m", basis: 0.1, score: 0.95,
      r_squared: 0.99, lattice_d: 12.5, status: "candidate",
      predicted_q: [0.1, 0.14], peaks: [],
    },
    {
      id: 2, exposure_id: 5, phase: "Im3m", basis: 0.15, score: 0.7,
      r_squared: 0.85, lattice_d: 9.0, status: "candidate",
      predicted_q: [0.12, 0.17], peaks: [],
    },
  ];
  const BASE_GROUP = { id: 1, exposure_id: 5, kind: "auto", active: true };
  let groupMembers: number[] = [];
  let addedIndexId: number | null = null;

  await seedState(page, { activeExperimentId: 1, activeSampleId: 10, activeExposureId: 5 });
  await mockCore(page, [{ id: 1, username: "alice" }]);

  // Override exposures for sample 10 to include our test exposure.
  await page.route("**/api/samples/10/exposures", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPOSURE]) }));

  await page.route("**/api/exposures/5/trace", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify({ q: [0.1, 0.2], I: [10, 20], sigma: [1, 1] }) }));
  await page.route("**/api/exposures/5/peaks", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/exposures/5/indices", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(CANDIDATES) }));

  // Groups reflect the current state of `groupMembers` on every fetch.
  await page.route("**/api/exposures/5/groups", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([{ ...BASE_GROUP, members: [...groupMembers] }]) }));

  await page.route("**/api/groups/1/members", async (route) => {
    const data = route.request().postDataJSON() as { index_id: number };
    addedIndexId = data.index_id;
    groupMembers = [data.index_id];
    await route.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify({ ...BASE_GROUP, members: groupMembers }) });
  });

  await page.goto("/");

  // Wait for PhasePanel to show the Pn3m candidate.
  await expect(page.getByText("Pn3m")).toBeVisible();

  // Click the "+" button on the first candidate (index id 1).
  await page.getByRole("button", { name: "Add index 1" }).click();

  // POST should have been sent with the correct index_id.
  await expect.poll(() => addedIndexId, { timeout: 2000 }).toBe(1);

  // After groups refetch, index 1 should appear in the Active set section.
  await expect(page.locator('[data-index-id="1"][data-active]')).toBeVisible();
});

test("reanalyze: stale-indices banner fires POST /analyze when clicked", async ({ page }) => {
  const EXPOSURE = {
    id: 5, sample_id: 10, filename: "scan1.dat", kind: "file",
    selected: true, tags: [], sources: [],
  };
  const STALE_INDEX = {
    id: 3, exposure_id: 5, phase: "Pn3m", basis: 0.1, score: 0.9,
    r_squared: 0.99, lattice_d: 12.5, status: "stale",
    predicted_q: [0.1, 0.14], peaks: [],
  };
  let analyzeCalled = false;

  await seedState(page, { activeExperimentId: 1, activeSampleId: 10, activeExposureId: 5 });
  await mockCore(page, [{ id: 1, username: "alice" }]);

  await page.route("**/api/samples/10/exposures", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPOSURE]) }));
  await page.route("**/api/exposures/5/trace", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify({ q: [0.1, 0.2], I: [10, 20], sigma: [1, 1] }) }));
  await page.route("**/api/exposures/5/peaks", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/exposures/5/indices", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([STALE_INDEX]) }));
  await page.route("**/api/exposures/5/groups", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([{ id: 1, exposure_id: 5, kind: "auto", active: true, members: [3] }]) }));

  await page.route("**/api/exposures/5/analyze", async (route) => {
    analyzeCalled = true;
    await route.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify({ id: 5, analyzed: true }) });
  });

  await page.goto("/");

  // StaleIndicesBanner should appear because the index has status "stale".
  await expect(page.getByRole("alert")).toBeVisible();
  await expect(page.getByRole("alert")).toContainText("stale");

  await page.getByRole("button", { name: /Re-analyze/ }).click();

  await expect.poll(() => analyzeCalled, { timeout: 2000 }).toBe(true);
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
