import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "../test-utils";
import { usePeaks, useRemovePeak } from "../../src/queries";
import { pendingDeferreds } from "../../src/lib/queue/deferred";
import { useAppState } from "../../src/state";
import * as api from "../../src/api";
import type { Peak, PeakRemoveResponse } from "../../src/api";

const E = 42;
const peak = (id: number): Peak => ({
  id, exposure_id: E, q: 0.1 * id, intensity: 1, prominence: 1, sharpness: 1,
  source: "manual", excluded: false,
});
const confirmed = (id: number) =>
  ({ event_id: id, view_row_id: id, analysis_inputs_hash: "h" }) as PeakRemoveResponse;

function setup() {
  let server = [peak(7), peak(8)];
  const listPeaks = vi.spyOn(api, "listPeaks").mockImplementation(async () => server.map((p) => ({ ...p })));
  const finish: ((r: PeakRemoveResponse) => void)[] = [];
  vi.spyOn(api, "removeManualPeak").mockImplementation(
    () => new Promise((resolve) => { finish.push(resolve); }));
  const qc = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={qc}>{children}</QueryClientProvider>
  );
  const hook = renderHook(() => ({ peaks: usePeaks(E), remove: useRemovePeak(E) }), { wrapper });
  return { hook, finish, listPeaks, setServer: (s: Peak[]) => { server = s; } };
}

describe("peak list refetch waits for pending peak mutations", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
    useAppState.setState({ username: "alice" });
  });

  it("a pending remove stays hidden while an earlier op's confirmation lands", async () => {
    const { hook, finish, listPeaks, setServer } = setup();
    await waitFor(() => expect(hook.result.current.peaks.data).toBeDefined());
    act(() => hook.result.current.remove.mutate(7));
    act(() => hook.result.current.remove.mutate(8));
    await waitFor(() => expect(hook.result.current.peaks.data?.map((p) => p.id)).toEqual([]));
    await waitFor(() => expect(finish).toHaveLength(1));
    const callsBefore = listPeaks.mock.calls.length;
    setServer([peak(8)]);
    await act(async () => { finish[0]!(confirmed(1)); });
    await new Promise((r) => setTimeout(r, 50));
    expect(hook.result.current.peaks.data?.map((p) => p.id)).toEqual([]);
    expect(listPeaks.mock.calls.length).toBe(callsBefore);
  });

  it("refetches once the last pending op confirms", async () => {
    const { hook, finish, listPeaks, setServer } = setup();
    await waitFor(() => expect(hook.result.current.peaks.data).toBeDefined());
    act(() => hook.result.current.remove.mutate(7));
    act(() => hook.result.current.remove.mutate(8));
    await waitFor(() => expect(finish).toHaveLength(1));
    setServer([peak(8)]);
    await act(async () => { finish[0]!(confirmed(1)); });
    await waitFor(() => expect(finish).toHaveLength(2));
    setServer([]);
    const callsBefore = listPeaks.mock.calls.length;
    await act(async () => { finish[1]!(confirmed(2)); });
    await waitFor(() => expect(listPeaks.mock.calls.length).toBe(callsBefore + 1));
    await waitFor(() => expect(hook.result.current.peaks.isFetching).toBe(false));
    expect(hook.result.current.peaks.data).toEqual([]);
  });
});
