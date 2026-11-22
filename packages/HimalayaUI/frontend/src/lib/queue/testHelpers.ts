import type { QueryClient, MutationCache } from "@tanstack/react-query";
import type { Mutator, SseEvent, FlatPayload } from "./types";
import { handleRemoteEvent } from "./replayCoordinator";
import { peakAddMutator } from "./mutators/peakAdd";
import { reanalyzeExposureMutator } from "./mutators/reanalyzeExposure";
import { newClientOpId } from "../clientOpId";
import { getClientId } from "../clientId";
import { queryKeys } from "../../queries";

export interface HimalayaTestHelpers {
  runMutator(name: "peak_added",
             flat: FlatPayload<{ q: number },
                               { exposureId: number; username: string | undefined; clientId: string }>): Promise<unknown>;
  runMutator(name: "reanalyze_exposure",
             flat: FlatPayload<Record<string, never>,
                               { exposureId: number; username: string | undefined; clientId: string }>): Promise<unknown>;
  injectSse(frame: SseEvent): void;
  getPeaks(exposureId: number): unknown;
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
  if (!import.meta.env.DEV) return;
  window.__himalayaTest = {
    runMutator: ((name: string, flat: any) => {
      const m = MUTATORS[name];
      if (!m) return Promise.reject(new Error(`unknown mutator: ${name}`));
      return runMutatorImpl(qc, m, flat);
    }) as HimalayaTestHelpers["runMutator"],
    injectSse: (frame) => handleRemoteEvent(frame, qc, mc),
    getPeaks: (exposureId) => qc.getQueryData(queryKeys.peaks(exposureId)),
    clientId: () => getClientId(),
    newClientOpId: () => newClientOpId(),
  };
}