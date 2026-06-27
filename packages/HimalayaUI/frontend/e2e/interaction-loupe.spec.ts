/**
 * E2E: Loupe interaction migration — mocked Playwright spec (Task 3.2).
 *
 * Exercises the shell-owned InteractionDock (two steppers: sample + frame axes),
 * compound cursor (useListCursor for frames + useStepperOnly for samples), and
 * the keyboard layer's cold-load scope-focus after the Loupe page was migrated
 * to the shell-level interaction registry (commits 41f4be3b + f57ab387).
 *
 * Dock testid contract (InteractionDock.tsx / DockStepper / DockUpLink):
 *   dock-prev-sample / dock-next-sample / dock-sample-count  — extra stepper (URL-driven)
 *   dock-prev-frame  / dock-next-frame  / dock-frame-count   — cursor stepper (in-page)
 *   dock-primary     — openFocus (Focus button; Enter target)
 *   dock-up-link     — DockUpLink "‹ Corpus"
 *   dock-action-cull / dock-action-representative  (Drop is the sole cull verb, a toggle; no Keep)
 *
 * Covered cases:
 *   IL-1  Both dock-sample-count AND dock-frame-count render (two steppers)
 *   IL-2  ArrowRight and dock-next-frame both advance the frame caption
 *   IL-3  ArrowDown and dock-next-sample both navigate to the sibling sample URL
 *   IL-4  Cold-load `x` (no prior click) fires the Drop PATCH — the scope div
 *          auto-focuses via callback ref when the detector column mounts; the
 *          keyboard layer sees the target inside [data-interaction-scope] and fires.
 *          jsdom unit tests FALSE-GREEN this (keyDown on window bypasses the scope
 *          guard); real Chromium + this spec is the only effective gate.
 *   IL-5  Enter navigates to /sample/<id> (Focus route, no /loupe suffix)
 *   IL-6  dock-primary click also navigates to Focus
 */

import { test, expect, type Page } from "@playwright/test";

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

const EXPERIMENT_ID = 1;

const EXPERIMENT = {
  id: EXPERIMENT_ID,
  name: "SSRL Test",
  path: "/data",
  data_dir: "/data/data",
  analysis_dir: "/data/analysis",
  manifest_path: null,
  created_at: "2026-01-01",
  last_scanned_at: "2026-01-01T10:00:00",
  ingest_status: "complete",
  stats: { loads: 1, samples: 2, exposures: 4, sessions: 1 },
};

// Two samples in the corpus so the sample stepper has a sibling.
const SAMPLE_A = {
  id: 10, experiment_id: EXPERIMENT_ID,
  name: "JC010", display_name: "JC010",
  notes: null, q_units: "A-1", tags: [],
};
const SAMPLE_B = {
  id: 11, experiment_id: EXPERIMENT_ID,
  name: "JC011", display_name: "JC011",
  notes: null, q_units: "A-1", tags: [],
};
const SAMPLES = [SAMPLE_A, SAMPLE_B];

// Three exposures for SAMPLE_A — gives a 3-frame axis to exercise.
// Exposure 101 is the representative (selected:true) → becomes defaultExposureId
// → seeded as frame 1 of 3 at cold load.
const EXPOSURE_101 = {
  id: 101, sample_id: 10, filename: "pos1.dat", kind: "file" as const,
  selected: true, status: null as null, image_path: null, image_version: "",
  tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
};
const EXPOSURE_102 = {
  id: 102, sample_id: 10, filename: "pos2.dat", kind: "file" as const,
  selected: false, status: null as null, image_path: null, image_version: "",
  tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
};
const EXPOSURE_103 = {
  id: 103, sample_id: 10, filename: "pos3.dat", kind: "file" as const,
  selected: false, status: "accepted" as const, image_path: null, image_version: "",
  tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
};
const EXPOSURES_A = [EXPOSURE_101, EXPOSURE_102, EXPOSURE_103];

// One exposure for SAMPLE_B (needed if the next page renders after navigation).
const EXPOSURE_201 = {
  id: 201, sample_id: 11, filename: "pos1.dat", kind: "file" as const,
  selected: true, status: null as null, image_path: null, image_version: "",
  tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
};
const EXPOSURES_B = [EXPOSURE_201];

/** Seed localStorage so onboarding/tutorial prompts do not gate navigation. */
async function seedState(page: Page): Promise<void> {
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({
        state: { username: "alice", tutorialSeen: true, theme: "dark" },
        version: 5,
      }),
    );
  });
}

/**
 * Register mocked API routes for the Loupe page.
 * Mirrors the route set from loupe-tags.spec.ts plus experiment/loads stubs
 * used by the app shell.
 */
async function mockLoupe(page: Page): Promise<void> {
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPERIMENT]) }));
  await page.route(`**/api/experiments/${EXPERIMENT_ID}`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPERIMENT) }));
  await page.route(`**/api/experiments/${EXPERIMENT_ID}/loads`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(SAMPLES) }));
  await page.route(`**/api/samples/${SAMPLE_A.id}/exposures*`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPOSURES_A) }));
  await page.route(`**/api/samples/${SAMPLE_B.id}/exposures*`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPOSURES_B) }));
  await page.route(`**/api/samples/${SAMPLE_A.id}/messages`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route(`**/api/samples/${SAMPLE_B.id}/messages`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/sample-tags", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  // SSE: empty stream so the hook settles without hanging the test.
  await page.route("**/api/events*", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
}

// ---------------------------------------------------------------------------
// Setup
// ---------------------------------------------------------------------------

test.beforeEach(async ({ page }) => {
  await seedState(page);
});

// ---------------------------------------------------------------------------
// IL-1: Both dock steppers (sample + frame) render
// ---------------------------------------------------------------------------

test("IL-1: dock-sample-count AND dock-frame-count both visible (two steppers in dock)", async ({ page }) => {
  await mockLoupe(page);
  await page.goto(`/sample/${SAMPLE_A.id}/loupe?experiment=${EXPERIMENT_ID}`);

  // Wait for the frame to load (Skeleton resolves → BigFrame caption visible).
  await expect(page.locator('[data-role="frame-caption"]')).toBeVisible();

  // Both steppers must appear in the shell dock.
  await expect(page.getByTestId("dock-sample-count")).toBeVisible();
  await expect(page.getByTestId("dock-frame-count")).toBeVisible();
});

// ---------------------------------------------------------------------------
// IL-2: Frame ArrowRight / dock-next-frame both advance the active frame
// ---------------------------------------------------------------------------

test("IL-2: ArrowRight and dock-next-frame both advance the frame caption (frame axis)", async ({ page }) => {
  await mockLoupe(page);
  await page.goto(`/sample/${SAMPLE_A.id}/loupe?experiment=${EXPERIMENT_ID}`);

  const caption = page.locator('[data-role="frame-caption"]');
  await expect(caption).toBeVisible();

  // EXPOSURE_101 is the representative → default → frame 1 of 3.
  await expect(caption).toHaveText(/frame 1 of 3/);

  // ArrowRight: the [data-interaction-scope] div is auto-focused on mount
  // (scopeRef callback), so the key fires its onKeyDown handler directly.
  await page.keyboard.press("ArrowRight");
  await expect(caption).toHaveText(/frame 2 of 3/);

  // ArrowLeft retreats back to frame 1.
  await page.keyboard.press("ArrowLeft");
  await expect(caption).toHaveText(/frame 1 of 3/);

  // dock-next-frame button advances the frame via the cursor stepper.
  await page.getByTestId("dock-next-frame").click();
  await expect(caption).toHaveText(/frame 2 of 3/);
});

// ---------------------------------------------------------------------------
// IL-3: Sample ArrowDown / dock-next-sample navigate to the sibling URL
// ---------------------------------------------------------------------------

test("IL-3: ArrowDown and dock-next-sample navigate to the sibling sample URL", async ({ page }) => {
  await mockLoupe(page);
  await page.goto(`/sample/${SAMPLE_A.id}/loupe?experiment=${EXPERIMENT_ID}`);
  await expect(page.locator('[data-role="frame-caption"]')).toBeVisible();

  // ArrowDown → sampleStepper.onNext() → gotoSample(SAMPLE_B.id).
  // [data-interaction-scope] onKeyDown handles ArrowDown.
  await page.keyboard.press("ArrowDown");
  await expect(page).toHaveURL(new RegExp(`/sample/${SAMPLE_B.id}/`));

  // Navigate back to SAMPLE_A so we can test the stepper button path.
  await mockLoupe(page);
  await page.goto(`/sample/${SAMPLE_A.id}/loupe?experiment=${EXPERIMENT_ID}`);
  await expect(page.locator('[data-role="frame-caption"]')).toBeVisible();

  // dock-next-sample button also calls sampleStepper.onNext().
  await page.getByTestId("dock-next-sample").click();
  await expect(page).toHaveURL(new RegExp(`/sample/${SAMPLE_B.id}/`));
});

// ---------------------------------------------------------------------------
// IL-4: Cold-load `x` (Drop) fires PATCH without any prior click — KEY CASE
// ---------------------------------------------------------------------------

test("IL-4: cold-load `x` fires Drop PATCH without prior click (scope auto-focus gate)", async ({ page }) => {
  // INTENT: the [data-interaction-scope] div receives focus via a callback ref
  // (scopeRef) when it mounts after data loads. Bare key `x` must fire the Drop
  // action through the keyboard layer because the auto-focus places `e.target`
  // inside [data-interaction-scope], satisfying inPageScope(). If auto-focus
  // were broken, `x` would be swallowed and the PATCH would never fire.
  // jsdom unit tests FALSE-GREEN this (keyDown on window bypasses the scope
  // guard); only real Chromium + this assertion catches the regression.
  await mockLoupe(page);

  let patchedStatus: string | null | undefined;
  // Register the PATCH handler AFTER mockLoupe so Playwright's most-recently-
  // added route wins for **/api/exposures/*/status.
  await page.route("**/api/exposures/*/status", async (route) => {
    const body = route.request().postDataJSON() as { status: string | null };
    patchedStatus = body.status;
    await route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({ id: EXPOSURE_101.id, status: body.status }),
    });
  });

  await page.goto(`/sample/${SAMPLE_A.id}/loupe?experiment=${EXPERIMENT_ID}`);

  // Wait for the frame caption — confirms the Skeleton resolved, the scope div
  // mounted, and scopeRef has focused it. No click before this point.
  await expect(page.locator('[data-role="frame-caption"]')).toHaveText(/frame 1 of 3/);

  // Press `x` — must fire via the keyboard layer (bare key inside scope).
  // DO NOT click first: clicking would trivially satisfy any focus guard and
  // would mask the exact auto-focus bug this case verifies.
  await page.keyboard.press("x");

  await expect.poll(() => patchedStatus).toBe("rejected");
});

// ---------------------------------------------------------------------------
// IL-5: Enter → Focus (no /loupe suffix)
// ---------------------------------------------------------------------------

test("IL-5: Enter navigates to /sample/<id> (Focus route)", async ({ page }) => {
  await mockLoupe(page);
  await page.goto(`/sample/${SAMPLE_A.id}/loupe?experiment=${EXPERIMENT_ID}`);
  await expect(page.locator('[data-role="frame-caption"]')).toBeVisible();

  // Enter fires the openFocus core action → navigate(`/sample/${sampleId}`).
  // Enter is NOT a bare key (e.key.length > 1), so inPageScope is not required —
  // it fires regardless of focus position.
  await page.keyboard.press("Enter");
  await expect(page).toHaveURL(new RegExp(`/sample/${SAMPLE_A.id}(?!/loupe)`));
});

// ---------------------------------------------------------------------------
// IL-6: dock-primary click → Focus
// ---------------------------------------------------------------------------

test("IL-6: dock-primary click navigates to Focus (/sample/<id> without /loupe)", async ({ page }) => {
  await mockLoupe(page);
  await page.goto(`/sample/${SAMPLE_A.id}/loupe?experiment=${EXPERIMENT_ID}`);
  await expect(page.locator('[data-role="frame-caption"]')).toBeVisible();

  // dock-primary renders the Focus button (variant="accent") bound to openFocus.
  await page.getByTestId("dock-primary").click();
  await expect(page).toHaveURL(new RegExp(`/sample/${SAMPLE_A.id}(?!/loupe)`));
});
