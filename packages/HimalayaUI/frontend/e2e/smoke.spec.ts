import { test, expect, type Page } from "@playwright/test";

const EXPERIMENT = {
  id: 1, name: "SSRL May 2026", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};
const SAMPLES = [
  { id: 10, experiment_id: 1, display_name: "D1", name: "cubic_run03", notes: null, tags: [] },
  { id: 11, experiment_id: 1, display_name: "D2", name: "hex_run01",   notes: null, tags: [] },
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
  // Corpus sample list — the focus workspace (/sample/:id, I4.4) learns each
  // sample's experiment_id from here. Keyed by sample id alone.
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(SAMPLES) }));
  for (const s of SAMPLES) {
    // Trailing `*` covers any future query string; without it requests fall through Vite proxy to whatever's on :8080.
    await page.route(`**/api/samples/${s.id}/exposures*`, (r) =>
      r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
    await page.route(`**/api/samples/${s.id}/messages`, (r) =>
      r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  }
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
  // I4.4 (#181): the three-card Index at `/` is retired; the focus workspace
  // (/sample/:id) is the surface that renders the trace plot. Onboarding shows
  // over it when no username is set.
  await page.goto("/sample/10");

  await expect(page.getByTestId("onboarding-name")).toBeVisible();
  await page.getByTestId("onboarding-new-handle").fill("alice");
  await page.getByTestId("onboarding-continue").click();

  await expect(page.getByTestId("onboarding-tutorial")).toBeVisible();
  await page.getByTestId("tutorial-next").click();
  await page.getByTestId("tutorial-next").click();
  await page.getByTestId("tutorial-next").click();
  await page.getByTestId("tutorial-done").click();

  await expect(page.getByTestId("onboarding-overlay")).not.toBeVisible();
  // Greenfield focus workspace: the retired three-card Index `plot-title` is
  // replaced by the TracePlate trace plot on the focus surface.
  await expect(page.getByTestId("trace-plate")).toBeVisible();
});

// I4.4 (#181): the legacy three-card Index's shell affordances — opening the
// NavModal from the plot title, the `,`/`.` Zustand-driven sample-step, and
// the in-Index ChatCard — are retired with the Index surface. The focus
// workspace (/sample/:id) is URL-routed (sample stepping happens via the
// corpus picker / URL, not a Zustand keyboard step) and does not host the
// NavModal or the chat card. Those interactions are covered on the focus
// surface by FocusWorkspacePage.*.test.tsx + qlink.spec.ts, so the three
// legacy-Index E2E cases that exercised them are removed here rather than
// repointed (they have no equivalent on the replacement surface).

// I4.4 + plotting redesign (#280): the auto-peaks/curation model — the `+`
// "Add index N" candidate buttons, the curated "Active set" with
// `data-index-id`/`data-active` rows, the `/groups/:id/members` POST, and the
// stale-indices alert (`role="alert"` + POST `/analyze` "Re-analyze") — were
// RETIRED by the greenfield rebuild. The focus surface now uses the assignment
// model (AssignmentRail/AssignmentCart driven by `confirmed_index`); there is no
// curation `+`, no active-set membership row, and no stale-indices banner /
// reanalyze affordance anywhere in `src/print` or `src/components`. The three
// E2E cases that exercised those surfaces ("curate: clicking + adds a
// candidate", "reanalyze: stale-indices banner fires POST /analyze", and
// "curate → reanalyze … round-trip") are removed here rather than repointed —
// they have no equivalent on the replacement surface (matching the legacy-Index
// removal precedent above). Replacement coverage for the assignment model lives
// in FocusWorkspacePage.*.test.tsx + the assignment Vitest suites.

test("a /compare* URL redirects to the series folio (Compare retired, #177)", async ({ page }) => {
  // I3.6 (#177): Compare is retired and replaced by the series stage. Every
  // /compare* deep-link now redirects to the series folio (/series); the old
  // Compare page + sidebar are gone.
  await seedState(page, { activeExperimentId: 1, activeSampleId: 10 });
  await mockCore(page, [{ id: 1, username: "alice" }]);
  // The series folio fetches its corpus-wide listing on mount.
  await page.route("**/api/series", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.goto("/experiments/1/compare");

  await expect(page).toHaveURL(/\/series$/);
  // Greenfield folio: folio-header replaces the retired series-folio-page testid.
  await expect(page.getByTestId("folio-header")).toBeVisible();
  await expect(page.getByTestId("compare-page")).toHaveCount(0);
});
