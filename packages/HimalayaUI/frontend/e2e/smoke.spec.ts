import { test, expect, type Page } from "./fixtures";

const EXPERIMENT = {
  id: 1, name: "SSRL May 2026", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};
const SAMPLES = [
  { id: 10, experiment_id: 1, display_name: "D1", name: "cubic_run03", notes: null },
  { id: 11, experiment_id: 1, display_name: "D2", name: "hex_run01",   notes: null },
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
      await route.fallback();
    }
  });
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([EXPERIMENT]) }));
  await page.route("**/api/experiments/1", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(EXPERIMENT) }));
  await page.route((url) => url.pathname === "/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(SAMPLES) }));
  for (const s of SAMPLES) {
    await page.route(`**/api/samples/${s.id}/exposures*`, (r) =>
      r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  }
  await page.route("**/api/events", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
}

async function seedState(page: Page, extra: Record<string, unknown>): Promise<void> {
  await page.addInitScript((state) => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state, version: 6 }),
    );
  }, { username: "alice", tutorialSeen: true, ...extra });
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
  await page.goto("/experiments");

  await expect(page.getByTestId("onboarding-name")).toBeVisible();
  await page.getByTestId("onboarding-new-handle").fill("alice");
  await page.getByTestId("onboarding-continue").click();

  await expect(page.getByTestId("onboarding-tutorial")).toBeVisible();
  await page.getByTestId("tutorial-next").click();
  await page.getByTestId("tutorial-next").click();
  await page.getByTestId("tutorial-done").click();

  await expect(page.getByTestId("onboarding-overlay")).not.toBeVisible();
  await expect(page.getByText("All experiments")).toBeVisible();
});

for (const path of ["/experiments/1/compare", "/index/lipid/JC001", "/sample/10"]) {
  test(`an old link ${path} goes home and says the page is gone`, async ({ page }) => {
    await seedState(page, {});
    await mockCore(page, [{ id: 1, username: "alice" }]);
    await page.goto(path);
    await expect(page).toHaveURL(/\/experiments$/);
    await expect(page.getByText("That page no longer exists.")).toBeVisible();
    await expect(page.getByTestId("new-experiment-cta")).toBeVisible();
  });
}

test("the experiments home follows a dark pick from the workspace, and its status pill has no extra letter-spacing", async ({ page }) => {
  await mockCore(page);
  await seedState(page, {});
  await page.addInitScript(() => { localStorage.setItem("himalaya-ui:theme", "dark"); });
  await page.goto("/experiments");
  const pill = page.getByTestId("notice-pill").filter({ hasText: "Up to date" });
  await expect(pill).toBeVisible();
  const read = () => page.evaluate(() => {
    const probe = document.createElement("div");
    probe.style.color = "oklch(0.165 0.005 280)";
    document.body.appendChild(probe);
    const darkPaper = getComputedStyle(probe).color;
    probe.remove();
    const pillEl = document.querySelector('[data-testid="notice-pill"]')!;
    return { onDarkPaper: getComputedStyle(document.body).backgroundColor === darkPaper, spacing: getComputedStyle(pillEl).letterSpacing };
  });
  await expect.poll(async () => (await read()).onDarkPaper).toBe(true);
  expect((await read()).spacing).toBe("normal");
});

test("the collection ribbon takes the same height whether or not a beamtime has image times", async ({ page }) => {
  await mockCore(page);
  const stats = { loads: 1, samples: 2, exposures: 4, sessions: 1, span_hours: 1, started_at: "2026-05-01T10:00:00", review: { reviewed: 1, total: 2 } };
  const timeline = { start: "2026-05-01T10:00:00", samples: [[0, 300, 1], [360, 700, 0]], beam_off: [700] };
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([
      { ...EXPERIMENT, id: 1, name: "Timed", ingest_status: "complete", stats: { ...stats, timeline } },
      { ...EXPERIMENT, id: 2, name: "Untimed", ingest_status: "complete", stats: { ...stats, timeline: null } },
    ]) }));
  await seedState(page, {});
  await page.goto("/experiments");
  const timed = page.getByTestId("experiment-card-1").getByTestId("collection-ribbon");
  const untimed = page.getByTestId("experiment-card-2").getByTestId("collection-ribbon");
  await expect(timed.getByText("1 to review")).toBeVisible();
  await expect(untimed.getByText("No image times yet")).toBeVisible();
  const [a, b] = await Promise.all([timed.boundingBox(), untimed.boundingBox()]);
  expect(a!.height).toBeGreaterThan(0);
  expect(b!.height).toBe(a!.height);
});

test("a dark pick is on the page before the app script runs, so a dark user never sees light paper first", async ({ page }) => {
  await page.route(/\/assets\/.*\.js$/, (r) => r.abort());
  await page.addInitScript(() => { localStorage.setItem("himalaya-ui:theme", "dark"); });
  await page.goto("/experiments", { waitUntil: "domcontentloaded" });
  const root = await page.evaluate(() => ({
    theme: document.documentElement.getAttribute("data-theme"),
    appMounted: document.getElementById("app")!.childElementCount > 0,
  }));
  expect(root).toEqual({ theme: "dark", appMounted: false });
});
