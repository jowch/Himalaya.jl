/**
 * E2E: Focus page interaction migration — mocked Playwright spec (Task 4.2).
 *
 * Exercises the shell-owned InteractionDock, candidate cursor (useListCursor
 * in headless / scope-focus mode), Enter-toggle assignment, the p/Escape mode
 * ladder, and the sample stepper after FocusPage was migrated to the
 * shell-level interaction registry (commits 288dbbce + 13a4796e).
 *
 * Dock testid contract (InteractionDock / DockStepper / DockUpLink):
 *   dock-prev-sample / dock-next-sample / dock-sample-count    — sample stepper (extra)
 *   dock-prev-candidate / dock-next-candidate / dock-candidate-count — cursor stepper
 *   dock-primary     — "Apply" (Enter target; gated on previewWasExplicit)
 *   dock-up-link     — DockUpLink "‹ Corpus"
 *   dock-action-addPeak  — "+ Peak" page action bound to key "p"
 *   dock-action-openLoupe — "Loupe" page action bound to key "l"
 *
 * Covered cases:
 *   IF-1  dock-sample-count AND dock-candidate-count both render; primary text "Apply"
 *   IF-2  ArrowRight/ArrowLeft and dock-next-candidate advance the candidate cursor;
 *          the previewed CandidateRow button gains data-previewed="true"
 *   IF-3  Enter fires assignment mutation (POST add); URL stays /sample/:id.
 *          dock-primary click does the same toggle. Both require previewWasExplicit.
 *   IF-4  Enter on a native <button> (dock-action-addPeak) does NOT fire assignment —
 *          isNativeInteractiveTarget guard intercepts it and lets the button activate
 *   IF-5  `p` cold-load (NO prior click) arms addPeak via scope auto-focus; Escape
 *          ladder: disarm addPeak → clear explicit preview → navigate back
 *   IF-6  ArrowDown / dock-next-sample navigate to the sibling sample URL
 */

import { test, expect, type Page } from "@playwright/test";

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

const EXPERIMENT_ID = 1;
const EXPOSURE_A_ID = 101;
const EXPOSURE_B_ID = 201;

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
  stats: { loads: 1, samples: 2, exposures: 2, sessions: 1 },
};

// Two samples in the corpus so the sample stepper has a sibling.
const SAMPLE_A = {
  id: 10, experiment_id: EXPERIMENT_ID,
  name: "JC010",
  notes: null, q_units: "A-1", tags: [],
};
const SAMPLE_B = {
  id: 11, experiment_id: EXPERIMENT_ID,
  name: "JC011",
  notes: null, q_units: "A-1", tags: [],
};
const SAMPLES = [SAMPLE_A, SAMPLE_B];

// One exposure for SAMPLE_A; selected=true so useAutoPickExposure picks it.
const EXPOSURE_A = {
  id: EXPOSURE_A_ID, sample_id: 10,
  filename: "pos1.dat", kind: "file" as const,
  selected: true, status: null as null,
  image_path: null, image_version: "",
  tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
};
// One exposure for SAMPLE_B (needed when the page re-renders after navigation).
const EXPOSURE_B = {
  id: EXPOSURE_B_ID, sample_id: 11,
  filename: "pos1.dat", kind: "file" as const,
  selected: true, status: null as null,
  image_path: null, image_version: "",
  tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
};

// Minimal trace with real q values.
const TRACE_A = {
  q: [0.05, 0.07, 0.10, 0.14, 0.17],
  I: [5000, 3000, 1800, 900, 450],
  sigma: [50, 35, 22, 12, 7],
};

// One auto peak on EXPOSURE_A (gives the index candidates something to claim).
const PEAK_A = {
  id: 501, exposure_id: EXPOSURE_A_ID,
  q: 0.07, intensity: 3000, prominence: 500, sharpness: 0.8,
  source: "auto" as const, excluded: false,
};

// TWO non-speculative (kind: "auto") candidates — required so the cursor has a
// pool to step AND Enter can toggle assignment.
const CANDIDATE_A = {
  id: 1, exposure_id: EXPOSURE_A_ID,
  phase: "Pn3m", basis: 14.2, score: 0.87,
  r_squared: null, lattice_d: 142, ngc: null,
  status: "candidate" as const, kind: "auto" as const,
  inputs_hash: null,
  peaks: [{ peak_id: 501, ratio_position: 1.0, residual: 0.001, q_observed: 0.07 }],
  predicted_q: [0.07],
};
const CANDIDATE_B = {
  id: 2, exposure_id: EXPOSURE_A_ID,
  phase: "Im3m", basis: 16.0, score: 0.65,
  r_squared: null, lattice_d: 160, ngc: null,
  status: "candidate" as const, kind: "auto" as const,
  inputs_hash: null,
  peaks: [{ peak_id: 501, ratio_position: 1.0, residual: 0.003, q_observed: 0.07 }],
  predicted_q: [0.07],
};

// Empty assignment (no phases in the call yet).
const ASSIGNMENT_A = {
  exposure_id: EXPOSURE_A_ID,
  state: "indexed" as const,
  members: [],
};

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
 * Register mocked API routes for the Focus page.
 * Mirrors interaction-loupe.spec.ts plus Focus-specific exposure routes
 * (trace, peaks, indices, assignment).
 */
async function mockFocus(page: Page): Promise<void> {
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPERIMENT]) }));
  await page.route(`**/api/experiments/${EXPERIMENT_ID}`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPERIMENT) }));
  await page.route(`**/api/experiments/${EXPERIMENT_ID}/loads`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  // GET /api/samples → corpus-wide listing (useCorpusSamples).
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(SAMPLES) }));
  await page.route(`**/api/samples/${SAMPLE_A.id}/exposures*`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPOSURE_A]) }));
  await page.route(`**/api/samples/${SAMPLE_B.id}/exposures*`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPOSURE_B]) }));
  await page.route(`**/api/samples/${SAMPLE_A.id}/messages`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route(`**/api/samples/${SAMPLE_B.id}/messages`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/sample-tags", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  // Exposure-scoped data for SAMPLE_A (representative exposure 101).
  await page.route(`**/api/exposures/${EXPOSURE_A_ID}/trace`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(TRACE_A) }));
  await page.route(`**/api/exposures/${EXPOSURE_A_ID}/peaks`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([PEAK_A]) }));
  await page.route(`**/api/exposures/${EXPOSURE_A_ID}/indices`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([CANDIDATE_A, CANDIDATE_B]) }));
  await page.route(`**/api/exposures/${EXPOSURE_A_ID}/assignment`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(ASSIGNMENT_A) }));
  // Exposure-scoped data for SAMPLE_B (minimal — navigation assertions only).
  await page.route(`**/api/exposures/${EXPOSURE_B_ID}/trace`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify({ q: [], I: [], sigma: [] }) }));
  await page.route(`**/api/exposures/${EXPOSURE_B_ID}/peaks`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route(`**/api/exposures/${EXPOSURE_B_ID}/indices`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route(`**/api/exposures/${EXPOSURE_B_ID}/assignment`, (r) =>
    r.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({ exposure_id: EXPOSURE_B_ID, state: "indexed", members: [] }),
    }));
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
// IF-1: Dock shape — both steppers visible; primary button text = "Apply"
// ---------------------------------------------------------------------------

test("IF-1: dock renders both sample+candidate steppers and primary button text is Apply", async ({ page }) => {
  await mockFocus(page);
  await page.goto(`/sample/${SAMPLE_A.id}`);

  // Wait for indices to load so candidatePool is non-empty and the cursor stepper
  // shows the 2-candidate count ("1 / 2"). This also confirms the scope div
  // mounted (isLoading=false) and usePageActions registered.
  await expect(page.getByTestId("dock-candidate-count")).toHaveText(/1 \/ 2/);

  // Both steppers must appear in the shell dock.
  await expect(page.getByTestId("dock-sample-count")).toBeVisible();
  await expect(page.getByTestId("dock-candidate-count")).toBeVisible();

  // Primary button label is "Apply" — core openFocus is overridden with
  // label: "Apply" in FocusPage's usePageActions declaration.
  await expect(page.getByTestId("dock-primary")).toContainText("Apply");
});

// ---------------------------------------------------------------------------
// IF-2: Candidate cursor — ArrowRight/ArrowLeft + dock-next-candidate
// ---------------------------------------------------------------------------

test("IF-2: ArrowRight/ArrowLeft and dock-next-candidate advance the candidate cursor", async ({ page }) => {
  await mockFocus(page);
  await page.goto(`/sample/${SAMPLE_A.id}`);

  // useListCursor initialises cursorId at ids[0] (CANDIDATE_A, id=1).
  // Stepper count = "1 / 2".
  await expect(page.getByTestId("dock-candidate-count")).toHaveText(/1 \/ 2/);

  // --- ArrowRight: moveBy(1) → cursor to CANDIDATE_B (id=2); previewWasExplicit=true ---
  // FocusPage's scope is auto-focused once isLoading→false (skeleton reveals
  // content). ArrowRight is handled by the scope's own onKeyDown; no prior click
  // or explicit focus is needed. This path is the cold-load regression guard.
  await page.keyboard.press("ArrowRight");
  await expect(page.getByTestId("dock-candidate-count")).toHaveText(/2 \/ 2/);
  // Im3m (CANDIDATE_B) is now previewed: CandidateRow → Card → data-previewed="true".
  await expect(page.getByRole("button", { name: /Im3m/ })).toHaveAttribute("data-previewed", "true");

  // --- ArrowLeft: moveBy(-1) → cursor back to CANDIDATE_A (id=1) ---
  await page.keyboard.press("ArrowLeft");
  await expect(page.getByTestId("dock-candidate-count")).toHaveText(/1 \/ 2/);
  // Pn3m (CANDIDATE_A) is now previewed.
  await expect(page.getByRole("button", { name: /Pn3m/ })).toHaveAttribute("data-previewed", "true");

  // --- dock-next-candidate: also calls moveBy(1) → CANDIDATE_B again ---
  await page.getByTestId("dock-next-candidate").click();
  await expect(page.getByTestId("dock-candidate-count")).toHaveText(/2 \/ 2/);
  await expect(page.getByRole("button", { name: /Im3m/ })).toHaveAttribute("data-previewed", "true");
});

// ---------------------------------------------------------------------------
// IF-3: Enter fires assignment mutation; dock-primary click does the same toggle
// ---------------------------------------------------------------------------

test("IF-3: Enter toggles assignment (not navigation); dock-primary click does the same toggle", async ({ page }) => {
  await mockFocus(page);

  // Track the two assignment mutation endpoints (registered AFTER mockFocus so
  // Playwright's most-recently-added route wins for these specific paths).
  // addAssignmentPhase:    POST   /api/exposures/:id/assignment/members
  // removeAssignmentPhase: DELETE /api/exposures/:id/assignment/members/:indexId
  let addFired = false;
  let removeFired = false;
  await page.route(`**/api/exposures/${EXPOSURE_A_ID}/assignment/members`, async (r) => {
    if (r.request().method() === "POST") {
      addFired = true;
      await r.fulfill({
        status: 200, contentType: "application/json",
        body: JSON.stringify({
          exposure_id: EXPOSURE_A_ID, state: "indexed",
          members: [CANDIDATE_B.id], event_id: 1, view_row_id: 1,
        }),
      });
    } else {
      await r.continue();
    }
  });
  await page.route(`**/api/exposures/${EXPOSURE_A_ID}/assignment/members/*`, async (r) => {
    if (r.request().method() === "DELETE") {
      removeFired = true;
      await r.fulfill({
        status: 200, contentType: "application/json",
        body: JSON.stringify({
          exposure_id: EXPOSURE_A_ID, state: "indexed",
          members: [], event_id: 2, view_row_id: null,
        }),
      });
    } else {
      await r.continue();
    }
  });

  await page.goto(`/sample/${SAMPLE_A.id}`);
  await expect(page.getByTestId("dock-candidate-count")).toHaveText(/1 \/ 2/);

  // Establish an explicit preview: ArrowRight moves cursor to CANDIDATE_B (id=2)
  // and sets previewWasExplicit=true. The scope is auto-focused after load, so
  // the ArrowRight reaches the scope's onKeyDown without any prior click.
  await page.keyboard.press("ArrowRight");
  await expect(page.getByTestId("dock-primary")).toBeEnabled();

  // Enter fires openFocus (key: ["Enter"]) through the keyboard layer.
  // Target is the scope div (focus-workspace) — isNativeInteractiveTarget(e)
  // returns false (not a button/link/input) so the guard does NOT block it.
  // candidateCursor.activate() → toggleAssignmentForId(2) → addAssignmentPhase.mutate(2)
  // → POST /api/exposures/101/assignment/members.
  await page.keyboard.press("Enter");
  // (a) Assignment mutation fired:
  await expect.poll(() => addFired).toBeTruthy();
  // (b) URL did NOT navigate — still on Focus route, no /loupe suffix:
  await expect(page).toHaveURL(new RegExp(`/sample/${SAMPLE_A.id}(?!/loupe)`));

  // dock-primary click calls primary.run() → same candidateCursor.activate().
  // CANDIDATE_B is now in the assignment (optimistic members=[2]) so activate()
  // calls removeAssignmentPhase.mutate(2) → DELETE.
  await page.getByTestId("dock-primary").click();
  await expect.poll(() => removeFired).toBeTruthy();
  await expect(page).toHaveURL(new RegExp(`/sample/${SAMPLE_A.id}(?!/loupe)`));
});

// ---------------------------------------------------------------------------
// IF-4: Enter on a native button does NOT toggle assignment (foundation guard)
// ---------------------------------------------------------------------------

test("IF-4: Enter on dock-action-addPeak does not fire assignment — isNativeInteractiveTarget guard", async ({ page }) => {
  await mockFocus(page);

  // Count any assignment mutation calls (POST or DELETE). Should stay 0.
  let assignmentCalls = 0;
  await page.route(`**/api/exposures/${EXPOSURE_A_ID}/assignment/members*`, async (r) => {
    const method = r.request().method();
    if (method === "POST" || method === "DELETE") assignmentCalls++;
    await r.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({ exposure_id: EXPOSURE_A_ID, state: "indexed", members: [], event_id: 1, view_row_id: null }),
    });
  });

  await page.goto(`/sample/${SAMPLE_A.id}`);
  await expect(page.getByTestId("dock-candidate-count")).toHaveText(/1 \/ 2/);

  // Establish explicit preview (previewWasExplicit=true) so openFocus is enabled
  // and WOULD fire if focus were on the scope div. The scope is auto-focused
  // after load, so ArrowRight reaches the scope's onKeyDown without a prior click.
  await page.keyboard.press("ArrowRight");
  await expect(page.getByTestId("dock-primary")).toBeEnabled();

  // Focus the dock-action-addPeak button — a native <button>.
  // When Enter is pressed, the keyboard layer checks isNativeInteractiveTarget:
  //   e.target.closest("button, a, [role=button], [aria-sort]") → the button itself
  //   → returns true → the layer skips openFocus (no assignment mutation).
  // The browser's native Enter behaviour activates the button instead (fires its
  // click handler → setAddArmed(v => !v) → addArmed=true).
  await page.getByTestId("dock-action-addPeak").focus();
  await page.keyboard.press("Enter");

  // Positive evidence: the button click DID fire — addArmed toggled to true →
  // TracePlate "+ Peak" button gets data-armed="true" (Button armed prop).
  // This confirms the Enter reached the button natively, not via the page action.
  await expect(page.locator('[data-armed="true"]')).toBeVisible();

  // The critical negative assertion: no assignment mutation was called.
  expect(assignmentCalls).toBe(0);
});

// ---------------------------------------------------------------------------
// IF-5: p cold-load (no prior click); Escape ladder: disarm → clear → navigate
// ---------------------------------------------------------------------------

test("IF-5: cold-load p arms addPeak via scope auto-focus; Escape ladder disarms then clears preview then navigates", async ({ page }) => {
  // INTENT: the [data-interaction-scope] div auto-focuses via its scopeRef
  // callback when data loads and the div mounts. Bare key `p` must arm addPeak
  // through the keyboard layer without any prior click into the page. If
  // auto-focus were broken, `p` would be blocked by inPageScope() and addPeak
  // would not arm. jsdom unit tests FALSE-GREEN this (they dispatch keyDown on
  // window directly, bypassing the scope guard); only real Chromium catches it.
  await mockFocus(page);
  await page.goto(`/sample/${SAMPLE_A.id}`);

  // Wait for the scope to auto-focus (isLoading→false; the useEffect fires and
  // el.focus() is NOT dropped because the skeleton has revealed real content).
  // "1 / 2" confirms indices loaded and the candidate cursor is seeded.
  // DO NOT click before this point — the scope must auto-focus for bare keys.
  await expect(page.getByTestId("dock-candidate-count")).toHaveText(/1 \/ 2/);

  // Press p — bare key; isBareKey=true so the keyboard layer checks inPageScope.
  // The scope div is now auto-focused (cold-load fix), so the layer fires the
  // addPeak action. This cold-load test is valid: no prior click anywhere on
  // the page.
  await page.keyboard.press("p");

  // addPeak armed: TracePlate "+ Peak" button receives armed=true →
  // Button renders data-armed="true" and aria-pressed="true".
  await expect(page.locator('[data-armed="true"]')).toBeVisible();

  // Arm the Escape ladder's middle rung by setting an explicit preview.
  // ArrowRight is handled by the scope's own onKeyDown. The scope is auto-
  // focused, so ArrowRight reaches it without an explicit .focus() call.
  await page.keyboard.press("ArrowRight");
  // dock-primary enables when previewWasExplicit=true && cursorId!==null.
  await expect(page.getByTestId("dock-primary")).toBeEnabled();

  // ----- Escape 1: highest-priority rung — disarm addPeak -----
  await page.keyboard.press("Escape");
  await expect(page.locator('[data-armed="true"]')).toHaveCount(0); // disarmed
  await expect(page).toHaveURL(new RegExp(`/sample/${SAMPLE_A.id}`));  // still here

  // ----- Escape 2: clear explicit preview (previewWasExplicit → false) -----
  // dock-primary goes disabled: enabled() = cursorId!==null && previewWasExplicit
  //   = true && false = false.
  await page.keyboard.press("Escape");
  await expect(page.getByTestId("dock-primary")).toBeDisabled();
  await expect(page).toHaveURL(new RegExp(`/sample/${SAMPLE_A.id}`));  // still here

  // ----- Escape 3: navigate back -----
  // goBack() → navigate(`/experiments/${experimentId}/corpus`)
  // = navigate("/experiments/1/corpus") since SAMPLE_A.experiment_id = 1.
  await page.keyboard.press("Escape");
  await expect(page).not.toHaveURL(new RegExp(`/sample/${SAMPLE_A.id}`));
});

// ---------------------------------------------------------------------------
// IF-6: Sample ↑/↓ + dock-next-sample navigate to the sibling URL
// ---------------------------------------------------------------------------

test("IF-6: ArrowDown and dock-next-sample navigate to the sibling sample URL", async ({ page }) => {
  await mockFocus(page);
  await page.goto(`/sample/${SAMPLE_A.id}`);
  await expect(page.getByTestId("dock-candidate-count")).toHaveText(/1 \/ 2/);

  // ArrowDown → scope onKeyDown fires sampleStepper.onNext()
  // → navigate(`/sample/${SAMPLE_B.id}`). The scope is auto-focused after load.
  await page.keyboard.press("ArrowDown");
  await expect(page).toHaveURL(new RegExp(`/sample/${SAMPLE_B.id}`));

  // Navigate back to SAMPLE_A to test the stepper button path.
  await mockFocus(page);
  await page.goto(`/sample/${SAMPLE_A.id}`);
  await expect(page.getByTestId("dock-candidate-count")).toHaveText(/1 \/ 2/);

  // dock-next-sample button calls sampleStepper.onNext() → same navigation.
  await page.getByTestId("dock-next-sample").click();
  await expect(page).toHaveURL(new RegExp(`/sample/${SAMPLE_B.id}`));
});
