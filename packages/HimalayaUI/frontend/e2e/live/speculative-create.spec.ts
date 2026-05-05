/**
 * Live-server spec for issue #35 Bugs 3 & 6 — speculative-index create.
 *
 * Two bugs, one UI flow: clicking "+ Add speculative" → picking phase + anchor
 * → "Save". Both bugs land on the resulting IndexEntry.
 *
 *   Bug 3 — `createSpeculative` SSE-wins phantom IndexEntry. Pre-fix, when the
 *     SSE frame for `speculative_created` beat the HTTP response back to the
 *     browser, the mutator's `onSuccess` ran with the synth payload (which
 *     spreads only the SSE event fields), leaving `phase: undefined` and
 *     `score: NaN` on the spliced cache row. The card rendered "NaN" and a
 *     blank phase chip. Post-fix: `onSuccess` guards on `response.phase` and
 *     falls back to `invalidateQueries({groups, indices})`.
 *
 *   Bug 6 — `insert_speculative_index!` missing `inputs_hash`. Pre-fix, the
 *     INSERT into `indices` left `inputs_hash` NULL. StaleIndicesBanner gates
 *     on (index.inputs_hash !== exposure.analysis_inputs_hash) → NULL ≠ hash
 *     → banner spuriously fired immediately after every speculative create.
 *     Post-fix: the inserted row inherits the exposure's current
 *     analysis_inputs_hash (correct by construction; same effective peak set).
 *
 * Both fixes are observable end-to-end: after Save, the new card MUST show a
 * real numeric score (not "NaN") and a non-empty phase chip; the banner MUST
 * stay hidden.
 */
import { test, expect, type Page } from "@playwright/test";

const BACKEND_BASE = process.env["BACKEND_BASE"] ?? "http://127.0.0.1:8090";

type Sample = { id: number; experiment_id: number };
type Exposure = { id: number; sample_id: number; selected: boolean };
type Peak = { id: number; q: number; source: string; excluded: boolean };
type IndexEntry = { id: number; kind: string };

interface Fixture {
  experimentId: number;
  sampleId: number;
  exposureId: number;
}

async function findFixture(): Promise<Fixture> {
  // Need a sample whose `selected:true` exposure has ≥2 non-excluded peaks
  // (the SpeculativeBuilder filters peaks by `!excluded`, and Lamellar needs
  // anchor + ratio 2). Page redirects to the selected exposure; using any
  // other one means the test runs against the wrong fixture.
  const exps = await fetch(`${BACKEND_BASE}/api/experiments`).then(r => r.json());
  for (const exp of exps) {
    const samples: Sample[] = await fetch(
      `${BACKEND_BASE}/api/experiments/${exp.id}/samples`).then(r => r.json());
    for (const s of samples) {
      const exposures: Exposure[] = await fetch(
        `${BACKEND_BASE}/api/samples/${s.id}/exposures`).then(r => r.json());
      const selected = exposures.find(e => e.selected);
      if (!selected) continue;
      const peaks: Peak[] = await fetch(
        `${BACKEND_BASE}/api/exposures/${selected.id}/peaks`).then(r => r.json());
      const usable = peaks.filter(p => !p.excluded);
      if (usable.length < 2) continue;
      return { experimentId: exp.id, sampleId: s.id, exposureId: selected.id };
    }
  }
  throw new Error(
    `No selected exposure with ≥2 peaks at ${BACKEND_BASE}.`);
}

async function seedAndOpen(page: Page, fx: Fixture): Promise<void> {
  // Seed Zustand-persisted state to land directly on the Index page for this
  // exposure — bypasses OnboardingFlow, NavModal, and any tutorial overlays.
  await page.addInitScript((args) => {
    const { username, expId, sampId, expoId } = args;
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: {
        username, firstName: username, lastName: "tester",
        activeExperimentId: expId,
        activeSampleId: sampId,
        activeExposureId: expoId,
        activePage: "index",
        tutorialSeen: true,
        theme: "dark",
      },
      version: 3,
    }));
  }, { username: "speculative-tester", expId: fx.experimentId,
       sampId: fx.sampleId, expoId: fx.exposureId });
  await page.goto("/");
  // PhasePanel renders the speculative entry button once indices load.
  await expect(page.getByTestId("add-speculative-button")).toBeVisible();
}

async function existingSpeculativeIds(exposureId: number): Promise<Set<number>> {
  const indices: IndexEntry[] = await fetch(
    `${BACKEND_BASE}/api/exposures/${exposureId}/indices`).then(r => r.json());
  return new Set(indices.filter(ix => ix.kind === "speculative").map(ix => ix.id));
}

test.describe("issue #35 speculative-create reconciliation", () => {
  let fx: Fixture;

  test.beforeAll(async () => {
    fx = await findFixture();
    console.log(
      `[speculative-create] using experiment=${fx.experimentId} ` +
      `sample=${fx.sampleId} exposure=${fx.exposureId}`);
  });

  async function buildSpeculative(page: Page): Promise<void> {
    await page.getByTestId("add-speculative-button").click();
    await expect(page.getByTestId("speculative-builder")).toBeVisible();
    await page.getByTestId("spec-phase-select").selectOption("Lamellar");
    // Wait for the peaks query to resolve and populate options. Index 0 is
    // the disabled "Choose…", so we need ≥2 options before we can pick a peak.
    const anchor = page.getByTestId("spec-anchor-select");
    await expect(anchor.locator("option")).not.toHaveCount(1, { timeout: 8000 });
    await anchor.selectOption({ index: 1 });
    await page.getByTestId("spec-save-button").click();
    await expect(page.getByTestId("speculative-builder"))
      .not.toBeVisible({ timeout: 8000 });
  }

  test("speculative renders with real phase and score (Bug 3)", async ({ page }) => {
    const before = await existingSpeculativeIds(fx.exposureId);
    await seedAndOpen(page, fx);
    await buildSpeculative(page);

    // Identify the new speculative by polling the indices list — the SSE
    // round-trip after the HTTP response can lag a few hundred ms.
    let newId: number | undefined;
    await expect.poll(async () => {
      const after = await existingSpeculativeIds(fx.exposureId);
      const fresh = [...after].filter(id => !before.has(id));
      if (fresh.length === 1) newId = fresh[0];
      return fresh.length;
    }, { timeout: 8000 }).toBe(1);

    const card = page.locator(`[data-index-id="${newId!}"]`);
    await expect(card).toBeVisible();

    // Bug 3 symptom: rendered card had `phase: undefined` (blank chip) and
    // `score: NaN` (literal "NaN" string from formatScore). Post-fix, both
    // come from the canonical refetch.
    const cardText = await card.innerText();
    expect(cardText, "score must not be NaN").not.toContain("NaN");
    expect(cardText, "phase chip must show real phase").toContain("Lamellar");
  });

  test("StaleIndicesBanner stays hidden after speculative create (Bug 6)", async ({ page }) => {
    const before = await existingSpeculativeIds(fx.exposureId);
    await seedAndOpen(page, fx);
    // Sanity: no banner at rest.
    await expect(page.getByRole("alert")).toHaveCount(0);

    await buildSpeculative(page);

    // Confirm the create actually landed — otherwise the no-banner assertion
    // would be a false positive (no spec ⇒ no inherited hash ⇒ no
    // mismatch ⇒ no banner, regardless of the fix).
    await expect.poll(async () => {
      const after = await existingSpeculativeIds(fx.exposureId);
      return [...after].filter(id => !before.has(id)).length;
    }, { timeout: 8000 }).toBe(1);

    // Pre-fix, the banner appeared within the StaleIndicesBanner debounce
    // (150ms) because the new index's inputs_hash was NULL. Wait 1 s to give
    // it a fair chance to surface, then assert it didn't.
    await page.waitForTimeout(1000);
    const banner = page.getByRole("alert").filter({ hasText: /stale/i });
    await expect(banner, "banner must stay hidden — new spec inherits hash").toHaveCount(0);
  });
});
