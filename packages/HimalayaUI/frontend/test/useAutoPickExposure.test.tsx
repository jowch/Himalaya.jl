import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { useAutoPickExposure } from "../src/hooks/useAutoPickExposure";
import { useAppState } from "../src/state";

const SAMPLE_ID = 10;

// rep = the `selected` exposure; the other two are switchable non-reps.
const EXPOSURES = [
  { id: 100, sample_id: SAMPLE_ID, filename: "e100.tif", kind: "file",
    selected: false, status: "accepted", image_path: null, image_version: "",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null },
  { id: 101, sample_id: SAMPLE_ID, filename: "e101.tif", kind: "file",
    selected: true, status: "accepted", image_path: null, image_version: "",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null },
  { id: 102, sample_id: SAMPLE_ID, filename: "e102.tif", kind: "file",
    selected: false, status: "accepted", image_path: null, image_version: "",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null },
];

function wrapper({ children }: { children: ReactNode }): JSX.Element {
  return <QueryClientProvider client={makeClient()}>{children}</QueryClientProvider>;
}

beforeEach(() => {
  vi.restoreAllMocks();
  useAppState.setState({ activeExposureId: undefined });
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify(EXPOSURES), {
      status: 200, headers: { "Content-Type": "application/json" },
    }),
  );
});

describe("useAutoPickExposure", () => {
  it("cold arrival (undefined) adopts the flagged representative", async () => {
    renderHook(() => useAutoPickExposure(SAMPLE_ID), { wrapper });
    await waitFor(() =>
      expect(useAppState.getState().activeExposureId).toBe(101),
    );
  });

  it("falls back to the first acceptable exposure when none is flagged", async () => {
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify(EXPOSURES.map((e) => ({ ...e, selected: false }))), {
        status: 200, headers: { "Content-Type": "application/json" },
      }),
    );
    renderHook(() => useAutoPickExposure(SAMPLE_ID), { wrapper });
    await waitFor(() =>
      expect(useAppState.getState().activeExposureId).toBe(100),
    );
  });

  it("does NOT clobber a deliberate switch to a valid non-rep exposure (F-11)", async () => {
    // The switcher has set a valid, non-representative exposure as active.
    useAppState.setState({ activeExposureId: 102 });
    renderHook(() => useAutoPickExposure(SAMPLE_ID), { wrapper });
    // Give the effect time to (wrongly) yank it back to the rep, then assert
    // it held the user's pick.
    await new Promise((r) => setTimeout(r, 50));
    expect(useAppState.getState().activeExposureId).toBe(102);
  });

  it("re-picks the flagged exposure when the current id is stale for the sample", async () => {
    // A leftover id from another sample is not valid here → adopt the flag.
    useAppState.setState({ activeExposureId: 999 });
    renderHook(() => useAutoPickExposure(SAMPLE_ID), { wrapper });
    await waitFor(() =>
      expect(useAppState.getState().activeExposureId).toBe(101),
    );
  });
});
