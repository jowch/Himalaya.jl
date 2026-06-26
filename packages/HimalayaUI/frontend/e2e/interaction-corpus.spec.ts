/**
 * E2E: Corpus interaction migration — mocked Playwright spec (Task 2.3).
 *
 * Exercises the roving-tabindex cursor, the shell-level InteractionDock
 * integration, and the cull-from-selection flow after the hand-built page
 * Dock was migrated to the shell-level InteractionDock (commit 4dd53b07).
 *
 * Dock testid contract (InteractionDock.tsx / DockStepper / DockUpLink):
 *   dock-prev-sample / dock-next-sample / dock-sample-count  — DockStepper
 *   dock-primary      — openFocus (Focus button, always; disabled when no cursor)
 *   dock-up-link      — DockUpLink "‹ Experiments"
 *   dock-action-<id>  — page actions: openLoupe, cull, keep, restore
 *
 * Covered cases:
 *   IC-1  Click a row → data-cursored="true"; Enter navigates to /sample/<id>
 *   IC-2  dock-next-sample and ArrowDown land on the SAME cursored row (parity)
 *   IC-3  Exactly one [role="row"][tabindex="0"] at a time (roving invariant)
 *   IC-4  Checking a row reveals dock-action-cull enabled; clicking it fires
 *          the status PATCH (mocked mutation route)
 *   IC-5  Pressing x when selection is active also fires the cull via the
 *          keyboard layer (the Checkbox is a <span>, not an <input>, so
 *          isTyping() returns false even while the checkbox has focus)
 */

import { test, expect, type Page } from "@playwright/test";

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

const EXPERIMENT_ID = 7;

const EXPERIMENT = {
  id: EXPERIMENT_ID,
  name: "SSRL Beamtime",
  path: "/data",
  data_dir: "/data",
  analysis_dir: "/analysis",
  manifest_path: null,
  created_at: "2026-01-01",
  last_scanned_at: "2026-01-01T10:00:00",
  ingest_status: "complete",
  stats: { loads: 1, samples: 3, exposures: 3, sessions: 1 },
};

function makeSamples(n: number): object[] {
  return Array.from({ length: n }, (_, i) => ({
    id: 100 + i,
    experiment_id: EXPERIMENT_ID,
    name: `Sample ${i + 1}`,
    display_name: `Sample ${i + 1}`,
    notes: null,
    tags: [],
    q_units: "A-1",
  }));
}

/** One accepted exposure per sample, id = 200 + (sampleId - 100). */
function makeExposures(sampleId: number): object[] {
  return [
    {
      id: 200 + (sampleId - 100),
      sample_id: sampleId,
      filename: `frame_${sampleId}.dat`,
      kind: "file",
      selected: true,
      status: "accepted",
      image_path: null,
      image_version: "",
      tags: [],
      sources: [],
      trace_hash: null,
      analysis_inputs_hash: null,
    },
  ];
}

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
 * Register mocked API routes for the corpus page.
 * Per-sample exposure routes are included so useCorpusExposures populates
 * byId, which cullSelected reads to build the batch mutation.
 */
async function mockCorpus(page: Page, sampleCount = 3): Promise<void> {
  const samples = makeSamples(sampleCount);

  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPERIMENT]) }));
  await page.route(`**/api/experiments/${EXPERIMENT_ID}`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPERIMENT) }));
  // Loads: empty → grouping-review-banner renders with data-state="clear".
  await page.route(`**/api/experiments/${EXPERIMENT_ID}/loads`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(samples) }));
  await page.route("**/api/sample-tags", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/series", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  // SSE: empty stream so the hook doesn't stay in-flight.
  await page.route("**/api/events*", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));

  // Per-sample exposures — needed by cullSelected (reads corpusExposures.byId).
  for (const s of samples as Array<{ id: number }>) {
    await page.route(`**/api/samples/${s.id}/exposures*`, (r) =>
      r.fulfill({
        status: 200, contentType: "application/json",
        body: JSON.stringify(makeExposures(s.id)),
      }));
  }
}

// ---------------------------------------------------------------------------
// Setup
// ---------------------------------------------------------------------------

test.beforeEach(async ({ page }) => {
  await seedState(page);
});

// ---------------------------------------------------------------------------
// IC-1: Click a row → cursored; Enter → navigate
// ---------------------------------------------------------------------------

test("IC-1: clicking a row parks the cursor there; Enter navigates to /sample/<id>", async ({ page }) => {
  await mockCorpus(page, 3);
  await page.goto(`/experiments/${EXPERIMENT_ID}/corpus`);

  const rows = page.getByTestId("sample-table-row");
  await expect(rows.first()).toBeVisible();

  // Initially the first row is cursored (useListCursor initializes to ids[0]).
  await expect(rows.first()).toHaveAttribute("data-cursored", "true");

  // Click the second row (sample id=101) to move the cursor.
  await rows.nth(1).click();
  await expect(rows.nth(1)).toHaveAttribute("data-cursored", "true");
  await expect(rows.nth(1)).toHaveAttribute("tabindex", "0");
  // First row loses the cursor.
  await expect(rows.first()).toHaveAttribute("data-cursored", "false");

  // Enter fires the openFocus action through the keyboard layer → navigate.
  await page.keyboard.press("Enter");
  await expect(page).toHaveURL(/\/sample\/101/);
});

// ---------------------------------------------------------------------------
// IC-2: Stepper ↔ ArrowDown parity
// ---------------------------------------------------------------------------

test("IC-2: dock-next-sample and ArrowDown move to the same cursored row (parity)", async ({ page }) => {
  await mockCorpus(page, 3);
  await page.goto(`/experiments/${EXPERIMENT_ID}/corpus`);

  const rows = page.getByTestId("sample-table-row");
  await expect(rows.first()).toBeVisible();

  // Stepper path: dock-next-sample → second row.
  await page.getByTestId("dock-next-sample").click();
  await expect(rows.nth(1)).toHaveAttribute("data-cursored", "true");

  // Reset: click the first row → cursor back to row 0.
  await rows.first().click();
  await expect(rows.first()).toHaveAttribute("data-cursored", "true");

  // Arrow path: ArrowDown is handled by the data-interaction-scope onKeyDown;
  // the cursored row has tabIndex=0 and focus (useListCursor effect), so the
  // event bubbles to the scope div and sampleCursor.moveBy(1) fires.
  await page.keyboard.press("ArrowDown");
  await expect(rows.nth(1)).toHaveAttribute("data-cursored", "true");
});

// ---------------------------------------------------------------------------
// IC-3: Roving tabindex invariant
// ---------------------------------------------------------------------------

test("IC-3: exactly one [role='row'][tabindex='0'] exists at a time (roving invariant)", async ({ page }) => {
  await mockCorpus(page, 3);
  await page.goto(`/experiments/${EXPERIMENT_ID}/corpus`);

  await expect(page.getByTestId("sample-table-row").first()).toBeVisible();

  // Initial state.
  await expect(page.locator('[role="row"][tabindex="0"]')).toHaveCount(1);

  // After first stepper advance.
  await page.getByTestId("dock-next-sample").click();
  await expect(page.locator('[role="row"][tabindex="0"]')).toHaveCount(1);

  // After second stepper advance (now on last row).
  await page.getByTestId("dock-next-sample").click();
  await expect(page.locator('[role="row"][tabindex="0"]')).toHaveCount(1);
});

// ---------------------------------------------------------------------------
// IC-4: Checking a row → dock-action-cull visible → clicking fires PATCH
// ---------------------------------------------------------------------------

test("IC-4: checking a row shows dock-action-cull enabled; clicking fires status PATCH", async ({ page }) => {
  await mockCorpus(page, 1);

  // Capture the PATCH request to /api/exposures/{id}/status (registered AFTER
  // mockCorpus so Playwright's most-recently-added route wins).
  let patchedStatus: string | null | undefined;
  await page.route("**/api/exposures/*/status", async (route) => {
    const body = route.request().postDataJSON() as { status: string | null };
    patchedStatus = body.status;
    await route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({ id: 200, status: body.status }),
    });
  });

  await page.goto(`/experiments/${EXPERIMENT_ID}/corpus`);
  await expect(page.getByTestId("sample-table-row").first()).toBeVisible();

  // Browse mode: mode-gated dock-action-cull is hidden (InteractionDock filters
  // out enabled()=false mode actions entirely — not just disabled, absent).
  await expect(page.getByTestId("dock-action-cull")).toHaveCount(0);

  // Check the row checkbox (role="checkbox" aria-label="Select sample") →
  // sampleCursor.selected gains the id → mode flips to "selection".
  await page.getByRole("checkbox", { name: "Select sample" }).first().click();

  // Selection mode: cull button appears and is enabled.
  await expect(page.getByTestId("dock-action-cull")).toBeVisible();
  await expect(page.getByTestId("dock-action-cull")).toBeEnabled();

  // Click the cull button → cullSelected("rejected") → PATCH /api/exposures/200/status.
  await page.getByTestId("dock-action-cull").click();

  await expect.poll(() => patchedStatus).toBe("rejected");
});

// ---------------------------------------------------------------------------
// IC-5: `x` keyboard shortcut fires the cull when a sample is selected
// ---------------------------------------------------------------------------

test("IC-5: pressing x with a sample selected fires the cull via keyboard layer", async ({ page }) => {
  await mockCorpus(page, 1);

  let patchedStatus: string | null | undefined;
  await page.route("**/api/exposures/*/status", async (route) => {
    const body = route.request().postDataJSON() as { status: string | null };
    patchedStatus = body.status;
    await route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({ id: 200, status: body.status }),
    });
  });

  await page.goto(`/experiments/${EXPERIMENT_ID}/corpus`);
  await expect(page.getByTestId("sample-table-row").first()).toBeVisible();

  // Check the row checkbox → selection active.
  // The Checkbox primitive is a <span role="checkbox">, not an <input>, so
  // matchKey.isTyping() returns false even while the span has DOM focus.
  // The keyboard layer therefore fires bare-key `x` correctly.
  await page.getByRole("checkbox", { name: "Select sample" }).first().click();
  await expect(page.getByTestId("dock-action-cull")).toBeVisible();

  // Press x — fires cull action via useKeyboardLayer → dropSelected().
  await page.keyboard.press("x");

  await expect.poll(() => patchedStatus).toBe("rejected");
});
