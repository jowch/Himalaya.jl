/**
 * Multi-tab cross-tab sync via SSE replay-as-rerun.
 *
 * This is the E2E half of the deep-scan response. The unit-test layer
 * encodes types (which can drift from routes); only an E2E that runs the
 * mutator end-to-end through real HTTP + real onSuccess catches contract
 * drift like the `PeakAddResponse.peak: Peak` bug we found by hand.
 *
 * Test affordances on `window.__himalayaTest` (DEV-only; tree-shaken in
 * prod) let Playwright drive mutators without the heavy Plot canvas UI.
 *
 * Scenarios:
 *   1. Tab A adds a peak; Tab A's cache reflects the SERVER row (real id,
 *      not a placeholder) — catches Bug #1 / shape mismatch.
 *   2. SSE delivers Tab A's add to Tab B; Tab B's cache picks up the same
 *      row using `payload.peak_curation_id` (issue #2).
 *   3. SSE self-echo to Tab A (same client_id, no matching deferred) is
 *      DROPPED — no duplicate row (issue #8).
 *   4. Tab B's reanalyze writes the new analysis_inputs_hash inline on
 *      HTTP success (deep-scan Bug #3).
 */
import { test, expect, type Page, type BrowserContext } from "@playwright/test";

// ---------------------------------------------------------------------------
// Mock-network setup shared between tabs.
// ---------------------------------------------------------------------------

interface MockState {
  peaks: Array<{
    id: number;
    exposure_id: number;
    q: number;
    intensity: number | null;
    prominence: number | null;
    sharpness: number | null;
    source: "auto" | "manual";
    excluded: boolean;
  }>;
  nextPeakId: number;
  nextEventId: number;
  inputsHash: string;
  hashCounter: number;
}

const SAMPLE_ID = 10;
const EXPOSURE_ID = 5;

function makeMockState(): MockState {
  return {
    peaks: [],
    nextPeakId: 100,
    nextEventId: 1,
    inputsHash: "h_initial",
    hashCounter: 0,
  };
}

async function injectFakeEventSource(page: Page): Promise<void> {
  // Replace EventSource with a FakeEventSource so the test can dispatch
  // SSE frames directly via window.__sseDispatch — no streaming-mock
  // gymnastics needed.
  await page.addInitScript(() => {
    const subscribers: { __emit(type: string, data: string): void }[] = [];
    class FakeES {
      private listeners: Record<string, ((e: unknown) => void)[]> = {};
      constructor(public url: string) {
        subscribers.push(this);
      }
      addEventListener(type: string, fn: (e: unknown) => void) {
        (this.listeners[type] ||= []).push(fn);
      }
      close() {
        const i = subscribers.indexOf(this);
        if (i >= 0) subscribers.splice(i, 1);
      }
      __emit(type: string, data: string) {
        for (const fn of (this.listeners[type] ?? [])) {
          fn({ data } as unknown);
        }
      }
    }
    (window as unknown as { EventSource: unknown }).EventSource = FakeES;
    (window as unknown as { __sseDispatch: (frame: unknown) => void }).__sseDispatch =
      (frame: unknown) => {
        const data = JSON.stringify(frame);
        for (const sub of subscribers) sub.__emit("curation", data);
      };
  });
}

async function mockBackend(page: Page, state: MockState, username: string): Promise<void> {
  // Users (just one)
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([{ id: 1, username }]) }));
  await page.route("**/api/health", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "{}" }));

  // Experiment / sample / exposure shells
  const exp = {
    id: 1, name: "Experiment", path: "/p", data_dir: "/p/data",
    analysis_dir: "/p/a", manifest_path: null, created_at: "2026-05-03",
  };
  const sample = {
    id: SAMPLE_ID, experiment_id: 1, display_name: "D1", name: "s1", notes: null, tags: [],
  };
  const exposure = {
    id: EXPOSURE_ID, sample_id: SAMPLE_ID, filename: "scan", kind: "file",
    selected: true, status: null, image_path: null,
    analysis_inputs_hash: state.inputsHash,
    tags: [], sources: [], image_version: 1,
  };
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([exp]) }));
  await page.route("**/api/experiments/1", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(exp) }));
  await page.route("**/api/experiments/1/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([sample]) }));
  // Corpus list — the focus workspace (/sample/:id, I4.4) learns the sample's
  // experiment_id from here.
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([sample]) }));
  await page.route(`**/api/samples/${SAMPLE_ID}/exposures*`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([exposure]) }));
  await page.route(`**/api/samples/${SAMPLE_ID}/messages`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route(`**/api/exposures/${EXPOSURE_ID}`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify({ ...exposure, analysis_inputs_hash: state.inputsHash }) }));
  await page.route(`**/api/exposures/${EXPOSURE_ID}/indices`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route(`**/api/exposures/${EXPOSURE_ID}/groups`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route(`**/api/exposures/${EXPOSURE_ID}/trace`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify({ q: [0.1], I: [10], sigma: [1] }) }));

  // Peaks: GET serves shared list; POST inserts and returns FLAT shape
  // (mirrors routes_peaks.jl — important: this is the shape under test,
  // not the type's claim). DELETE removes.
  await page.route(`**/api/exposures/${EXPOSURE_ID}/peaks`, async (route) => {
    const req = route.request();
    if (req.method() === "POST") {
      const data = req.postDataJSON() as { q: number };
      const peak = {
        id: state.nextPeakId++, exposure_id: EXPOSURE_ID, q: data.q,
        intensity: null, prominence: null, sharpness: null,
        source: "manual" as const, excluded: false,
      };
      state.peaks.push(peak);
      state.hashCounter++;
      state.inputsHash = `h_after_${state.hashCounter}`;
      // FLAT response — Peak fields + metadata. Matches the actual route.
      // The unit test fixtures used to nest under {peak: ...}; that bug
      // was caught only by reading the route source. This mock encodes
      // the route, so it'd fail if the mutator regressed.
      await route.fulfill({
        status: 201, contentType: "application/json",
        body: JSON.stringify({
          ...peak,
          event_id: state.nextEventId++,
          view_row_id: peak.id,
          analysis_inputs_hash: state.inputsHash,
        }),
      });
    } else {
      await route.fulfill({
        status: 200, contentType: "application/json",
        body: JSON.stringify(state.peaks),
      });
    }
  });

  await page.route(`**/api/exposures/${EXPOSURE_ID}/analyze`, async (route) => {
    state.hashCounter++;
    state.inputsHash = `h_after_${state.hashCounter}`;
    await route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({ id: EXPOSURE_ID, analyzed: true,
                             analysis_inputs_hash: state.inputsHash }),
    });
  });
}

async function seedSession(page: Page, username: string): Promise<void> {
  await page.addInitScript((u) => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: {
        username: u,
        activePage: "compare",
        tutorialSeen: true,
        theme: "dark",
        activeExperimentId: 1,
        activeSampleId: SAMPLE_ID,
        activeExposureId: EXPOSURE_ID,
      },
      version: 3,
    }));
  }, username);
}

async function setupTab(
  ctx: BrowserContext, username: string, state: MockState,
): Promise<Page> {
  const page = await ctx.newPage();
  await injectFakeEventSource(page);
  await mockBackend(page, state, username);
  await seedSession(page, username);
  await page.goto(`/sample/${SAMPLE_ID}`);
  // Wait for the test helpers to attach (App's useEffect runs after mount).
  await page.waitForFunction(() => Boolean((window as unknown as {
    __himalayaTest?: unknown;
  }).__himalayaTest));
  // Seed the per-id exposure cache so mutators that conditionally update it
  // (`old ? {...old, ...} : old`) have something to update against. The
  // production app populates this via mention resolution / inspect page;
  // the index-only test path doesn't, so seed deterministically.
  await page.evaluate((args) => {
    (window as unknown as { __himalayaTest: {
      seedExposure(id: number, e: unknown): void;
    }}).__himalayaTest.seedExposure(args.id, {
      id: args.id, sample_id: args.sampleId, filename: "scan", kind: "file",
      selected: true, status: null, image_path: null,
      analysis_inputs_hash: args.hash,
      tags: [], sources: [], image_version: 1,
    });
  }, { id: EXPOSURE_ID, sampleId: SAMPLE_ID, hash: state.inputsHash });
  return page;
}

// ---------------------------------------------------------------------------
// Scenarios
// ---------------------------------------------------------------------------

test("two-tab: peak_added flows from Tab A onSuccess into Tab B's cache via SSE", async ({ browser }) => {
  const state = makeMockState();
  const ctxA = await browser.newContext();
  const ctxB = await browser.newContext();
  try {
    const pageA = await setupTab(ctxA, "alice", state);
    const pageB = await setupTab(ctxB, "bob", state);

    // Read each tab's client_id (they should differ — sessionStorage is
    // per-context).
    const clientA = await pageA.evaluate(() =>
      (window as unknown as { __himalayaTest: { clientId(): string } })
        .__himalayaTest.clientId());
    const clientB = await pageB.evaluate(() =>
      (window as unknown as { __himalayaTest: { clientId(): string } })
        .__himalayaTest.clientId());
    expect(clientA).not.toBe(clientB);

    // Run peakAdd through the real mutator on Tab A. The HTTP mock returns
    // the FLAT shape; the mutator's onSuccess must destructure metadata
    // and write the Peak into queryKeys.peaks(EXPOSURE_ID). If the type
    // were stale (response.peak.id), this would crash.
    const opIdA = await pageA.evaluate(() =>
      (window as unknown as { __himalayaTest: { newClientOpId(): string } })
        .__himalayaTest.newClientOpId());
    await pageA.evaluate(async (args) => {
      const t = (window as unknown as { __himalayaTest: {
        runMutator(name: string, flat: unknown): Promise<unknown>;
      }}).__himalayaTest;
      await t.runMutator("peak_added", {
        kind: "peak_added",
        clientOpId: args.opId,
        exposureId: 5,
        username: "alice",
        clientId: args.clientId,
        q: 0.42,
        payload: { q: 0.42 },
      });
    }, { opId: opIdA, clientId: clientA });

    // Tab A: cache contains the SERVER peak with positive id.
    const peaksA = await pageA.evaluate(() =>
      (window as unknown as { __himalayaTest: { getPeaks(id: number): unknown[] } })
        .__himalayaTest.getPeaks(5));
    expect(peaksA).toHaveLength(1);
    expect((peaksA as { id: number; q: number }[])[0].q).toBe(0.42);
    expect((peaksA as { id: number }[])[0].id).toBeGreaterThan(0);
    const serverPeakId = (peaksA as { id: number }[])[0].id;

    // Tab A also: exposure cache should have the new analysis_inputs_hash.
    const expA = await pageA.evaluate(() =>
      (window as unknown as { __himalayaTest: { getExposure(id: number): unknown } })
        .__himalayaTest.getExposure(5));
    expect((expA as { analysis_inputs_hash: string }).analysis_inputs_hash)
      .toMatch(/^h_after_/);

    // Now dispatch the SSE frame for that add to Tab B (foreign tab).
    await pageB.evaluate((args) => {
      (window as unknown as { __sseDispatch: (f: unknown) => void })
        .__sseDispatch({
          id: 1, kind: "peak_added", entity_type: "exposure", entity_id: 5,
          actor: "alice",
          client_id: args.fromClient, // different from B's tab id
          client_op_id: args.opId,
          payload: { q: 0.42, peak_curation_id: args.peakId },
          post_state: { analysis_inputs_hash: args.hash, indices: [] },
        });
    }, { fromClient: clientA, opId: opIdA, peakId: serverPeakId, hash: "h_after_1" });

    const peaksB = await pageB.evaluate(() =>
      (window as unknown as { __himalayaTest: { getPeaks(id: number): unknown[] } })
        .__himalayaTest.getPeaks(5));
    expect(peaksB).toHaveLength(1);
    expect((peaksB as { id: number; q: number }[])[0].id).toBe(serverPeakId);
    expect((peaksB as { q: number }[])[0].q).toBe(0.42);
  } finally {
    await ctxA.close();
    await ctxB.close();
  }
});

test("self-echo: own SSE frame with no matching deferred is dropped (no duplicate row)", async ({ browser }) => {
  const state = makeMockState();
  const ctx = await browser.newContext();
  try {
    const page = await setupTab(ctx, "alice", state);
    const clientId = await page.evaluate(() =>
      (window as unknown as { __himalayaTest: { clientId(): string } })
        .__himalayaTest.clientId());

    // Add a peak via the real mutator. HTTP-first wins; the deferred is
    // cleared in mutationFn finally (the runMutator helper clears too).
    const opId = await page.evaluate(() =>
      (window as unknown as { __himalayaTest: { newClientOpId(): string } })
        .__himalayaTest.newClientOpId());
    await page.evaluate(async (args) => {
      const t = (window as unknown as { __himalayaTest: {
        runMutator(name: string, flat: unknown): Promise<unknown>;
      }}).__himalayaTest;
      await t.runMutator("peak_added", {
        kind: "peak_added",
        clientOpId: args.opId,
        exposureId: 5,
        username: "alice",
        clientId: args.clientId,
        q: 0.5,
        payload: { q: 0.5 },
      });
    }, { opId, clientId });

    // One peak in the cache.
    let peaks = await page.evaluate(() =>
      (window as unknown as { __himalayaTest: { getPeaks(id: number): unknown[] } })
        .__himalayaTest.getPeaks(5));
    expect(peaks).toHaveLength(1);
    const id1 = (peaks as { id: number }[])[0].id;

    // Now the late SSE for our own op arrives. client_id matches our tab.
    // No deferred is registered (mutator already finished). Without the
    // self-echo guard, applyRemoteToCache.peak_added would re-insert the
    // row — duplicate.
    await page.evaluate((args) => {
      (window as unknown as { __sseDispatch: (f: unknown) => void })
        .__sseDispatch({
          id: 1, kind: "peak_added", entity_type: "exposure", entity_id: 5,
          actor: "alice",
          client_id: args.clientId, // OWN client_id
          client_op_id: args.opId,
          payload: { q: 0.5, peak_curation_id: args.peakId },
          post_state: { analysis_inputs_hash: "h_after_1", indices: [] },
        });
    }, { clientId, opId, peakId: id1 });

    peaks = await page.evaluate(() =>
      (window as unknown as { __himalayaTest: { getPeaks(id: number): unknown[] } })
        .__himalayaTest.getPeaks(5));
    expect(peaks).toHaveLength(1); // still one — self-echo dropped
    expect((peaks as { id: number }[])[0].id).toBe(id1);
  } finally {
    await ctx.close();
  }
});

test("reanalyze response writes analysis_inputs_hash to exposure cache (deep-scan Bug #3)", async ({ browser }) => {
  const state = makeMockState();
  state.inputsHash = "h_initial";
  const ctx = await browser.newContext();
  try {
    const page = await setupTab(ctx, "alice", state);
    const clientId = await page.evaluate(() =>
      (window as unknown as { __himalayaTest: { clientId(): string } })
        .__himalayaTest.clientId());

    // Pre-condition: exposure cache has h_initial.
    const before = await page.evaluate(() =>
      (window as unknown as { __himalayaTest: { getExposure(id: number): unknown } })
        .__himalayaTest.getExposure(5));
    expect((before as { analysis_inputs_hash: string }).analysis_inputs_hash)
      .toBe("h_initial");

    // Run reanalyze through the real mutator. The HTTP mock advances the
    // hash; the mutator's onSuccess MUST write it onto the exposure cache.
    // Pre-fix the onSuccess was a no-op — banner would flicker between
    // HTTP success and SSE arrival.
    const opId = await page.evaluate(() =>
      (window as unknown as { __himalayaTest: { newClientOpId(): string } })
        .__himalayaTest.newClientOpId());
    await page.evaluate(async (args) => {
      const t = (window as unknown as { __himalayaTest: {
        runMutator(name: string, flat: unknown): Promise<unknown>;
      }}).__himalayaTest;
      await t.runMutator("reanalyze_exposure", {
        kind: "reanalyze_exposure",
        clientOpId: args.opId,
        exposureId: 5,
        username: "alice",
        clientId: args.clientId,
        payload: {},
      });
    }, { opId, clientId });

    const after = await page.evaluate(() =>
      (window as unknown as { __himalayaTest: { getExposure(id: number): unknown } })
        .__himalayaTest.getExposure(5));
    expect((after as { analysis_inputs_hash: string }).analysis_inputs_hash)
      .toMatch(/^h_after_/);
  } finally {
    await ctx.close();
  }
});
