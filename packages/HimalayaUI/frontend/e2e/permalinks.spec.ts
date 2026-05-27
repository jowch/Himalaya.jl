import { test, expect, type Page } from "@playwright/test";

// Spec §8.2 — Playwright mocked. /api/* is intercepted so this runs
// without a backend.
//
// Forced-enabler additions vs. plan verbatim (lines 2806–2905):
//   - Seed username/tutorialSeen in localStorage for tests 1–3 so the
//     onboarding overlay (which gates on `username === undefined`) doesn't
//     intercept clicks on the rocker / stale-url CTA. The plan code didn't
//     seed it; the overlay covers the page until dismissed.
//   - Use a regex URL matcher for `/api/resolve` rather than the glob
//     `**/api/resolve?**`. Playwright globs are micromatch-style and `?`
//     is "match any single char" — the literal query-string `?` either
//     under- or over-matches depending on encoding. Regex sidesteps it.
//   - Add `/api/resolve` mock to test 4. Cold-mount with seeded ids has
//     an empty TanStack cache, so `useStateFromUrl`'s root branch falls
//     through to `api.resolve({experiment_id, sample_id})` to recover
//     the slug names. Without the mock, the catch redirects to /index.
//   - Add a small delay to test 1's resolve mock so the ResolvingFallback
//     (`data-testid='resolving'`) is observable. With instant fulfill the
//     transition window is microseconds and `toBeAttached()` polls miss it.

const RESOLVE_RE = /\/api\/resolve\?/;

async function seedSession(page: Page): Promise<void> {
  // Match the smoke spec pattern: clear localStorage, then seed minimal
  // identity + tutorial-seen so OnboardingFlow stays mounted-but-hidden.
  await page.addInitScript(() => {
    localStorage.clear();
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: {
        username: "alice", firstName: undefined, lastName: undefined,
        tutorialSeen: true, theme: "dark",
        activePage: "index",
        activeExperimentId: undefined,
        activeSampleId: undefined,
        activeExposureId: undefined,
      },
      version: 3,
    }));
  });
}

test.beforeEach(async ({ page }) => {
  // Common stubs: list experiments, list samples for the relevant experiment.
  await page.route("**/api/experiments", (route) => {
    route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify([
        { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "",
          manifest_path: null, created_at: "", q_units: null },
      ]),
    });
  });
  await page.route("**/api/experiments/17/samples", (route) => {
    route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify([
        { id: 42, experiment_id: 17, name: "JC001", display_name: null,
          notes: null, tags: [] },
      ]),
    });
  });
  // Catch-all guards so unmocked API calls don't fall through to a real
  // backend (Vite proxy → :8080) or hang.
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/samples/42/exposures*", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/samples/42/messages", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
});

test("paste deep URL: lands on right page, no flash of wrong content", async ({ page }) => {
  await seedSession(page);
  await page.route(RESOLVE_RE, async (route) => {
    // Tiny delay so ResolvingFallback is observable. Without it the route
    // fulfills synchronously and React commits resolving:true→false in one
    // microtask window, missing Playwright's polling.
    await new Promise((r) => setTimeout(r, 150));
    route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    });
  });
  await page.goto("/index/lipid/JC001");
  // ResolvingFallback should appear briefly (or already past); page must not
  // show a different sample's content.
  await expect(page.locator("[data-testid='resolving']")).toBeAttached();
  // After resolve, page should be at the right sample.
  await expect(page).toHaveURL(/\/index\/lipid\/JC001$/);
});

test("paste stale URL: 404 page → CTA opens NavModal at right step", async ({ page }) => {
  await seedSession(page);
  await page.route(RESOLVE_RE, (route) => {
    route.fulfill({
      status: 404, contentType: "application/json",
      body: JSON.stringify({
        error: "not_found", missing: "experiment", missing_value: "lipid-typo",
        experiment_resolved: undefined, sample_resolved: undefined,
      }),
    });
  });
  await page.goto("/index/lipid-typo/JC001");
  await expect(page.locator("[data-testid='stale-url-page']")).toBeVisible();
  await expect(page.locator("[data-testid='stale-url-page']")).toHaveAttribute("data-missing", "experiment");
  await page.locator("[data-testid='stale-url-cta']").click();
  await expect(page.locator("[data-testid='nav-modal']")).toBeVisible();
});

test("TabRocker: the Index tab holds the /<page>/<exp>/<sample> slug URL", async ({ page }) => {
  // I1.7 (#163): Inspect is retired, so the old index→inspect→back continuity
  // leg is gone (Compare emits no slug URL — it returns `current`). This now
  // just pins that the Index tab keeps the resolved slug URL.
  await seedSession(page);
  await page.route(RESOLVE_RE, (route) => {
    route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    });
  });
  await page.goto("/index/lipid/JC001");
  await expect(page.locator("[data-testid='tab-index']")).toBeVisible();
  await expect(page).toHaveURL(/\/index\/lipid\/JC001$/);
  // The retired Inspect tab must no longer render.
  await expect(page.locator("[data-testid='tab-inspect']")).toHaveCount(0);
});

test("/ cold-mount: replaces to last-active slug URL", async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: {
        activePage: "index",
        activeExperimentId: 17,
        activeSampleId: 42,
        username: "test", firstName: undefined, lastName: undefined,
        activeExposureId: undefined,
        tutorialSeen: true, theme: "dark",
      },
      version: 3,
    }));
  });
  // Cold-mount root-redirect path: empty TanStack cache → resolve-by-id.
  await page.route(RESOLVE_RE, (route) => {
    route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    });
  });
  await page.goto("/");
  await expect(page).toHaveURL(/\/index\/lipid\/JC001$/);
});
