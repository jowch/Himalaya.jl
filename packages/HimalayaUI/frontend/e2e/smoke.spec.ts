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
  // Permalinks: cold mount with seeded activeExperimentId/activeSampleId
  // triggers useStateFromUrl's resolve-by-id fallback when the TanStack
  // cache hasn't hydrated yet. Same path also fires on every navigate
  // (e.g. `.` → next sample). The mock inspects the query so multi-sample
  // tests get the right (experiment, sample) pair back.
  await page.route(/\/api\/resolve\?/, (route) => {
    const url = new URL(route.request().url());
    const params = url.searchParams;
    // Match by name first (slug URL navigation), then by id (cold-mount).
    const sampleName = params.get("sample") ?? undefined;
    const sampleIdStr = params.get("sample_id") ?? undefined;
    let sample = SAMPLES[0];
    if (sampleName !== undefined) {
      sample = SAMPLES.find((s) => s.name === sampleName) ?? SAMPLES[0];
    } else if (sampleIdStr !== undefined) {
      const sId = Number(sampleIdStr);
      sample = SAMPLES.find((s) => s.id === sId) ?? SAMPLES[0];
    }
    route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({
        experiment_id: 1, experiment_name: "SSRL May 2026",
        sample_id: sample.id, sample_name: sample.name,
        exposure_id: undefined, exposure_filename: undefined,
      }),
    });
  });
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
    // Trailing `*` covers any future query string; without it requests fall through Vite proxy to whatever's on :8080.
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
      JSON.stringify({ state, version: 3 }),
    );
  }, { username: "alice", tutorialSeen: true, theme: "dark", ...extra });
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

test("a deleted /compare URL shows the stale-URL page", async ({ page }) => {
  await seedState(page, { activeExperimentId: 1, activeSampleId: 10 });
  await mockCore(page, [{ id: 1, username: "alice" }]);
  await page.goto("/experiments/1/compare");
  await expect(page.getByTestId("stale-url-page")).toBeVisible();
});

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
