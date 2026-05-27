/**
 * Live-server spec for issue #35 Bug 1 — own-op `post_state` propagation gap.
 *
 * Pre-fix, when the SSE frame carrying `post_state` for an own peak op beat
 * the HTTP body back to the browser, the self-echo path in
 * `replayCoordinator.handleRemoteEvent` ran but never applied `post_state`
 * to the cache. That left `exposure.analysis_inputs_hash` stale, and
 * StaleIndicesBanner — gated on (index.inputs_hash !== exposure.hash) —
 * incorrectly fired and persisted until the next refresh.
 *
 * Post-fix: `applyPostStateOnly(remote, qc)` is extracted from
 * `applyRemoteToCache` and called from BOTH Case 1 (own-op SSE) and the
 * self-echo path. The HTTP-wins path was always fine (mutator's `onSuccess`
 * writes `analysis_inputs_hash` directly from the route response).
 *
 * The user-visible promise is the same regardless of which path won: after
 * a UI-driven peak op, the StaleIndicesBanner does NOT appear. Asserting
 * that end-to-end exercises both paths in production conditions; the
 * deterministic SSE-wins coverage lives in `replayCoordinator.test.ts`.
 *
 * Driving the peak op via UI: clicks on the trace-viewer interior route to
 * `addPeak` (empty area), `togglePeakExclusion` (auto peak), or `removePeak`
 * (manual peak) depending on pixel proximity. All three paths run through
 * `useQueueMutation` and return `{analysis_inputs_hash}`. Test asserts the
 * peak-set changed and the banner stays hidden — bug 1 fails on either.
 */
import { test, expect, type Page } from "@playwright/test";

const BACKEND_BASE = process.env["BACKEND_BASE"] ?? "http://127.0.0.1:8090";

type Sample = { id: number; experiment_id: number };
type Exposure = { id: number; sample_id: number; selected: boolean };
type Peak = { id: number; q: number; source: string; excluded: boolean };

interface Fixture {
  experimentId: number;
  sampleId: number;
  exposureId: number;
}

async function findFixture(): Promise<Fixture> {
  // Need an analyzed exposure with ≥1 peak that's marked `selected:true` on
  // its sample. The Index page's PlotCard auto-switches activeExposureId to
  // the selected exposure (PlotCard.tsx:106), overriding any seeded value —
  // so picking the unselected one means the test would mutate the wrong
  // exposure.
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
      if (peaks.length < 1) continue;
      return { experimentId: exp.id, sampleId: s.id, exposureId: selected.id };
    }
  }
  throw new Error(
    `No analyzed exposure with peaks at ${BACKEND_BASE}. ` +
    `Need a sample with a selected:true exposure that has ≥1 peak.`);
}

async function seedAndOpen(page: Page, fx: Fixture): Promise<void> {
  await page.addInitScript((args) => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: {
        username: args.username, firstName: args.username, lastName: "tester",
        activeExperimentId: args.expId,
        activeSampleId: args.sampId,
        activeExposureId: args.expoId,
        tutorialSeen: true,
        theme: "dark",
      },
      version: 3,
    }));
  }, { username: "peak-tester", expId: fx.experimentId,
       sampId: fx.sampleId, expoId: fx.exposureId });
  await page.goto("/");
  // Trace must finish loading before we can click it.
  const trace = page.getByTestId("trace-viewer");
  await expect(trace).toBeVisible();
  // Trace SVG also needs to lay out; PlotCard is the broader shell.
  await expect(page.getByTestId("plot-card")).toBeVisible();
  // EventSource('/api/events') is opened in App.tsx's mount effect. The
  // server-side subscriber registration completes ~50–200 ms after the GET
  // — but if the click fires before that, the post_state SSE frame is
  // emitted to a not-yet-registered subscriber and lost. 800 ms buffers
  // for slow first-load JS bundle + SSE handshake without making the test
  // noticeably slower.
  await page.waitForTimeout(800);
}

async function fetchPeaks(exposureId: number): Promise<Peak[]> {
  return fetch(`${BACKEND_BASE}/api/exposures/${exposureId}/peaks`).then(r => r.json());
}

function peakSetSignature(peaks: Peak[]): string {
  // Insensitive to ordering; sensitive to id, source, excluded — enough to
  // detect any of {add, remove, exclude-toggle}.
  return peaks.map(p => `${p.id}:${p.source}:${p.excluded}`).sort().join("|");
}

test.describe("issue #35 own-op peak op leaves banner hidden (Bug 1)", () => {
  let fx: Fixture;

  test.beforeAll(async () => {
    fx = await findFixture();
    console.log(
      `[peak-add] using experiment=${fx.experimentId} ` +
      `sample=${fx.sampleId} exposure=${fx.exposureId}`);
  });

  test("clicking the trace triggers a peak op without surfacing the stale banner",
    async ({ page }) => {
      const beforeSig = peakSetSignature(await fetchPeaks(fx.exposureId));

      await seedAndOpen(page, fx);

      // Click on the trace at a position safely inside the plot interior.
      // TraceViewer's plot has marginLeft=50 marginTop=36 marginBottom=40
      // marginRight=14, and `insideInterior` rejects clicks in axis areas.
      // We click 110 px from the wrapper's left edge (well past the 50 px
      // margin) and vertically centered. Use locator-based click with
      // explicit position so Playwright dispatches the mouse events on the
      // element that owns the click handler (Observable Plot's SVG is the
      // child of trace-viewer; locator clicks bubble correctly).
      // If the click lands within PEAK_HIT_PX (10) of an existing peak, it
      // toggles exclusion / removes a manual — all three paths exercise
      // Bug 1's post_state propagation.
      const trace = page.getByTestId("trace-viewer");
      const box = await trace.boundingBox();
      if (!box) throw new Error("trace-viewer has no bounding box");
      await trace.click({
        position: { x: 110, y: Math.floor(box.height * 0.50) },
      });

      // Wait for the peak set to change on the server. Polls /peaks until
      // the signature shifts — covers add, remove, exclude-toggle. 15 s
      // because the cold-start mutation round-trip + post-commit broadcast
      // queue can take several seconds the first time the worker hits a
      // route on a freshly-bounced backend.
      await expect.poll(
        async () => peakSetSignature(await fetchPeaks(fx.exposureId)),
        { timeout: 15_000, intervals: [200, 400, 800, 1500] },
      ).not.toBe(beforeSig);

      // Wait for the post-mutation hash on indices and exposure to align on
      // the server. Banner gates on local-cache mismatch, but the test browser
      // still needs to have observed the new hashes — either via SSE post_state
      // (the fix's primary path) or via background refetch (ultimate fallback).
      // Polling the API confirms server-side reanalysis has settled before we
      // assert on UI.
      await expect.poll(async () => {
        const ix = await fetch(`${BACKEND_BASE}/api/exposures/${fx.exposureId}/indices`)
          .then(r => r.json());
        const exp = await fetch(`${BACKEND_BASE}/api/exposures/${fx.exposureId}`)
          .then(r => r.json());
        if (ix.length === 0) return true;
        return ix.every((i: { inputs_hash?: string }) =>
          i.inputs_hash === exp.analysis_inputs_hash);
      }, { timeout: 10_000 }).toBe(true);

      // Bug 1's symptom is "StaleIndicesBanner stays VISIBLE after own peak
      // op until manual refresh". Post-fix the banner clears as soon as
      // applyPostStateOnly propagates (own-op SSE) OR the indices query
      // refetches in the background. Pre-fix, neither happens — banner
      // sticks indefinitely. 12 s is generous enough to cover slow refetch
      // settling on a freshly-bounced backend, while still failing the
      // pre-fix behavior (which is "indefinite").
      const banner = page.getByRole("alert").filter({ hasText: /stale/i });
      await expect(
        banner,
        "StaleIndicesBanner must clear after own peak op (don't stick)",
      ).toHaveCount(0, { timeout: 12_000 });
    });
});
