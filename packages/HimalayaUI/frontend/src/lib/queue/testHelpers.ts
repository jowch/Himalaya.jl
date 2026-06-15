import type { QueryClient, MutationCache } from "@tanstack/react-query";
import type { Mutator, SseEvent, FlatPayload } from "./types";
import { handleRemoteEvent } from "./replayCoordinator";
import { peakAddMutator } from "./mutators/peakAdd";
import { reanalyzeExposureMutator } from "./mutators/reanalyzeExposure";
import { newClientOpId } from "../clientOpId";
import { getClientId } from "../clientId";
import { queryKeys } from "../../queries";

/**
 * E2E test affordances exposed on `window.__himalayaTest` in DEV builds only.
 * Gated by `import.meta.env.DEV`; production bundles do NOT carry this surface
 * (Vite tree-shakes the conditional `if (!import.meta.env.DEV) return;` block).
 *
 * What it exposes:
 *
 * - `runMutator(mutatorName, flat)`: runs a mutator end-to-end (onMutate →
 *   request → onSuccess) WITHOUT React, so Playwright tests can exercise the
 *   real HTTP path and real cache writes against mock-network fixtures. The
 *   E2E layer is the only place that catches contract drift between the
 *   Julia route response shape and the TypeScript mutator's read pattern
 *   (the unit test fixtures encode the type, not the route).
 *
 * - `injectSse(frame)`: pushes an SSE frame through `handleRemoteEvent` to
 *   exercise the foreign-event replay-as-rerun path. Multi-tab tests use
 *   this to deliver one tab's mutation to another tab's cache.
 *
 * - `getQueryData(...)`: read out cache state for assertions.
 *
 * - `clientId()`: per-tab identity, so tests can spoof self-echo frames.
 *
 * The set of mutators exposed is intentionally small (peak_added,
 * reanalyze_exposure) — extend as new test scenarios need them.
 */
export interface HimalayaTestHelpers {
  runMutator(name: "peak_added",
             flat: FlatPayload<{ q: number },
                               { exposureId: number; username: string | undefined; clientId: string }>): Promise<unknown>;
  runMutator(name: "reanalyze_exposure",
             flat: FlatPayload<Record<string, never>,
                               { exposureId: number; username: string | undefined; clientId: string }>): Promise<unknown>;
  injectSse(frame: SseEvent): void;
  getPeaks(exposureId: number): unknown;
  getExposure(exposureId: number): unknown;
  /** Seed the exposure cache so mutators with `old ? {...old, ...} : old`
   *  guards have something to update against. Pre-fetching via api.* in
   *  tests is racy; this lets the test stage the cache deterministically. */
  seedExposure(exposureId: number, exposure: unknown): void;
  clientId(): string;
  newClientOpId(): string;
}

declare global {
  interface Window {
    __himalayaTest?: HimalayaTestHelpers;
  }
}

const MUTATORS: Record<string, Mutator<any, any, any>> = {
  peak_added: peakAddMutator,
  reanalyze_exposure: reanalyzeExposureMutator,
};

async function runMutatorImpl(
  qc: QueryClient,
  m: Mutator<any, any, any>,
  flat: any,
): Promise<unknown> {
  const ctx = m.onMutate(flat, qc);
  try {
    const response = await m.request(flat, new AbortController().signal);
    m.onSuccess(flat, response, qc);
    return response;
  } catch (err) {
    ctx.restore?.();
    throw err;
  }
}

export function exposeTestHelpers(qc: QueryClient, mc: MutationCache): void {
  // Gate on DEV — vite tree-shakes the helpers out of production bundles.
  if (!import.meta.env.DEV) return;
  window.__himalayaTest = {
    runMutator: ((name: string, flat: any) => {
      const m = MUTATORS[name];
      if (!m) return Promise.reject(new Error(`unknown mutator: ${name}`));
      return runMutatorImpl(qc, m, flat);
    }) as HimalayaTestHelpers["runMutator"],
    injectSse: (frame) => handleRemoteEvent(frame, qc, mc),
    getPeaks: (exposureId) => qc.getQueryData(queryKeys.peaks(exposureId)),
    getExposure: (exposureId) => qc.getQueryData(queryKeys.exposure(exposureId)),
    seedExposure: (exposureId, exposure) =>
      qc.setQueryData(queryKeys.exposure(exposureId), exposure),
    clientId: () => getClientId(),
    newClientOpId: () => newClientOpId(),
  };
}
