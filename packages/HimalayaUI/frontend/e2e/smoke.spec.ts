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
  await expect(page.getByTestId("plot-title")).toBeVisible();
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
  await page.route("**/api/samples/10/exposures*", (r) =>
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

  await page.goto("/sample/10");

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
  // Stale index: inputs_hash differs from the exposure's analysis_inputs_hash.
  // The banner derives staleness from hash mismatch, not the legacy
  // status='stale' enum (which R3 retired).
  const EXPOSURE = {
    id: 5, sample_id: 10, filename: "scan1.dat", kind: "file",
    selected: true, tags: [], sources: [],
    trace_hash: "newhash", analysis_inputs_hash: "newhash",
  };
  const STALE_INDEX = {
    id: 3, exposure_id: 5, phase: "Pn3m", basis: 0.1, score: 0.9,
    r_squared: 0.99, lattice_d: 12.5, status: "candidate",
    inputs_hash: "oldhash",
    predicted_q: [0.1, 0.14], peaks: [],
  };
  let analyzeCalled = false;

  await seedState(page, { activeExperimentId: 1, activeSampleId: 10, activeExposureId: 5 });
  await mockCore(page, [{ id: 1, username: "alice" }]);

  await page.route("**/api/exposures/5", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPOSURE) }));
  await page.route("**/api/samples/10/exposures*", (r) =>
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

  await page.goto("/sample/10");

  // StaleIndicesBanner should appear because the index has status "stale".
  await expect(page.getByRole("alert")).toBeVisible();
  await expect(page.getByRole("alert")).toContainText("stale");

  await page.getByRole("button", { name: /Re-analyze/ }).click();

  await expect.poll(() => analyzeCalled, { timeout: 2000 }).toBe(true);
});

test("curate → reanalyze: active-set membership survives a reanalysis round-trip", async ({ page }) => {
  // Index is stale because its inputs_hash differs from the exposure's
  // current analysis_inputs_hash — the post-R3 derivation. The legacy
  // status='stale' enum was retired.
  const EXPOSURE = {
    id: 5, sample_id: 10, filename: "scan1.dat", kind: "file",
    selected: true, tags: [], sources: [],
    trace_hash: "newhash", analysis_inputs_hash: "newhash",
  };
  const BASE_INDEX = {
    id: 1, exposure_id: 5, phase: "Pn3m", basis: 0.1, score: 0.95,
    r_squared: 0.99, lattice_d: 12.5, predicted_q: [0.1, 0.14], peaks: [],
    status: "candidate",
  };
  let indexInputsHash = "oldhash";
  let analyzeCalled = false;

  await seedState(page, { activeExperimentId: 1, activeSampleId: 10, activeExposureId: 5 });
  await mockCore(page, [{ id: 1, username: "alice" }]);

  await page.route("**/api/exposures/5", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPOSURE) }));
  await page.route("**/api/samples/10/exposures*", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPOSURE]) }));
  await page.route("**/api/exposures/5/trace", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify({ q: [0.1, 0.2], I: [10, 20], sigma: [1, 1] }) }));
  await page.route("**/api/exposures/5/peaks", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/exposures/5/indices", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([{ ...BASE_INDEX, inputs_hash: indexInputsHash }]) }));
  await page.route("**/api/exposures/5/groups", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([{ id: 1, exposure_id: 5, kind: "auto", active: true, members: [1] }]) }));

  await page.route("**/api/exposures/5/analyze", async (route) => {
    analyzeCalled = true;
    indexInputsHash = "newhash";  // post-reanalyze, hashes match → banner clears
    await route.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify({ id: 5, analyzed: true }) });
  });

  await page.goto("/sample/10");

  // Index 1 is active (curated) but stale — banner must appear.
  await expect(page.locator('[data-index-id="1"][data-active]')).toBeVisible();
  await expect(page.getByRole("alert")).toContainText("stale");

  // Trigger reanalysis.
  await page.getByRole("button", { name: /Re-analyze/ }).click();
  await expect.poll(() => analyzeCalled, { timeout: 2000 }).toBe(true);

  // Group membership must survive — index 1 stays active after the refetch.
  await expect(page.locator('[data-index-id="1"][data-active]')).toBeVisible();
  // Stale banner must be gone once indices are fresh.
  await expect(page.getByRole("alert")).not.toBeVisible();
});

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
