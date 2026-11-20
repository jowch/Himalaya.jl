import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { usePeaks, useRemovePeak, queryKeys } from "../src/queries";
import { pendingDeferreds } from "../src/lib/queue/deferred";
import { createReconnectRefresh } from "../src/lib/queue/reconnectRefresh";
import { useAppState } from "../src/state";
import * as api from "../src/api";
import type { Peak, PeakRemoveResponse } from "../src/api";

const EXPOSURE_ID = 42;
const peak = (id: number): Peak => ({
  id, exposure_id: EXPOSURE_ID, q: 0.1 * id, intensity: 1, prominence: 1, sharpness: 1, source: "manual", excluded: false,
});

describe("refresh after the live-update stream reconnects", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
    useAppState.setState({ username: "alice" });
  });

  it("marks every cached query stale except raw traces", async () => {
    const qc = makeClient();
    qc.setQueryData(queryKeys.peaks(EXPOSURE_ID), [peak(1)]);
    qc.setQueryData(["samples", "corpus"], []);
    qc.setQueryData(queryKeys.trace(EXPOSURE_ID), { q: [], I: [], sigma: [] });
    createReconnectRefresh(qc, qc.getMutationCache())();
    await waitFor(() => expect(qc.getQueryState(queryKeys.peaks(EXPOSURE_ID))?.isInvalidated).toBe(true));
    expect(qc.getQueryState(["samples", "corpus"])?.isInvalidated).toBe(true);
    expect(qc.getQueryState(queryKeys.trace(EXPOSURE_ID))?.isInvalidated).toBe(false);
  });

  it("a peak removed while the backend was down stays removed once the refetch lands", async () => {
    let server = [peak(7), peak(8)];
    vi.spyOn(api, "listPeaks").mockImplementation(async () => server.map((p) => ({ ...p })));
    let finishRemove!: (r: PeakRemoveResponse) => void;
    vi.spyOn(api, "removePeak").mockImplementation(() => new Promise((resolve) => { finishRemove = resolve; }));

    const qc = makeClient();
    const wrapper = ({ children }: { children: ReactNode }) => (
      <QueryClientProvider client={qc}>{children}</QueryClientProvider>
    );
    const { result } = renderHook(() => ({ peaks: usePeaks(EXPOSURE_ID), remove: useRemovePeak(EXPOSURE_ID) }), { wrapper });
    await waitFor(() => expect(result.current.peaks.data?.map((p) => p.id)).toEqual([7, 8]));

    act(() => result.current.remove.mutate(7));
    await waitFor(() => expect(result.current.peaks.data?.map((p) => p.id)).toEqual([8]));

    createReconnectRefresh(qc, qc.getMutationCache())();
    await new Promise((r) => setTimeout(r, 20));
    expect(result.current.peaks.data?.map((p) => p.id)).toEqual([8]);

    server = [peak(8)];
    await act(async () => {
      finishRemove({ event_id: 1, view_row_id: 7, analysis_inputs_hash: "h" } as PeakRemoveResponse);
    });
    await waitFor(() => expect(api.listPeaks).toHaveBeenCalledTimes(2));
    await waitFor(() => expect(result.current.peaks.isFetching).toBe(false));
    expect(result.current.peaks.data?.map((p) => p.id)).toEqual([8]);
  });
});
