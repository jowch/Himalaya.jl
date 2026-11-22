/**
 * The server's peak_listing decides which auto peaks are excluded or replaced
 * by a user peak (grid-aware, nearest-first). The client must not re-derive
 * that with its own tolerance, so every confirmed peak mutation, own or
 * remote, invalidates the exposure's peak list and lets the refetch settle it.
 */
import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import { peakAddMutator } from "../../src/lib/queue/mutators/peakAdd";
import { peakRemoveMutator } from "../../src/lib/queue/mutators/peakRemove";
import {
  peakExcludeMutator, peakUnexcludeMutator,
} from "../../src/lib/queue/mutators/peakSetExcluded";
import { queryKeys } from "../../src/queries";
import type { Peak } from "../../src/api";

const EXPOSURE = 5;
const peak = (over: Partial<Peak>): Peak => ({
  id: 1, exposure_id: EXPOSURE, q: 0.2, intensity: 1, prominence: 1,
  sharpness: 1, source: "auto", excluded: false, ...over,
});

function setup() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  qc.setQueryData(queryKeys.peaks(EXPOSURE), [peak({}), peak({ id: 9, source: "manual", q: 0.5 })]);
  const spy = vi.spyOn(qc, "invalidateQueries");
  const peaksInvalidated = () =>
    spy.mock.calls.some(([f]) =>
      JSON.stringify((f as { queryKey?: unknown })?.queryKey) === JSON.stringify(queryKeys.peaks(EXPOSURE)));
  return { qc, peaksInvalidated };
}

const remote = (kind: string, payload: Record<string, unknown>) => ({
  id: 1, kind, entity_type: "exposure", entity_id: EXPOSURE, payload,
});

describe("remote peak events refetch the peak list", () => {
  const cases: [string, Record<string, unknown>][] = [
    ["peak_added", { q: 0.6, peak_curation_id: 12 }],
    ["peak_excluded", { q: 0.2, auto_peak_id: 1 }],
    ["peak_unexcluded", { q: 0.2, auto_peak_id: 1 }],
    ["peak_removed", { q: 0.5, peak_curation_id: 9 }],
  ];
  for (const [kind, payload] of cases) {
    it(`${kind}`, () => {
      const { qc, peaksInvalidated } = setup();
      applyRemoteToCache(remote(kind, payload), qc);
      expect(peaksInvalidated()).toBe(true);
    });
  }
});

describe("own peak mutations refetch the peak list after confirmation", () => {
  const scope = { exposureId: EXPOSURE, username: "a", clientId: "c", clientOpId: "o" };
  const meta = { event_id: 1, view_row_id: 1, analysis_inputs_hash: "h" };

  it("peak_added", () => {
    const { qc, peaksInvalidated } = setup();
    peakAddMutator.onSuccess({ ...scope, q: 0.6 } as never,
      { ...peak({ id: 12, source: "manual", q: 0.6 }), ...meta } as never, qc);
    expect(peaksInvalidated()).toBe(true);
  });
  for (const [name, m] of [["peak_excluded", peakExcludeMutator], ["peak_unexcluded", peakUnexcludeMutator]] as const) {
    it(name, () => {
      const { qc, peaksInvalidated } = setup();
      m.onSuccess({ ...scope, peakId: 1 } as never, { ...peak({}), ...meta } as never, qc);
      expect(peaksInvalidated()).toBe(true);
    });
  }
  it("peak_removed", () => {
    const { qc, peaksInvalidated } = setup();
    peakRemoveMutator.onSuccess({ ...scope, peakId: 9 } as never, {} as never, qc);
    expect(peaksInvalidated()).toBe(true);
  });
});
