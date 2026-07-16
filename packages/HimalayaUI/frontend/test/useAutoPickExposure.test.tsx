import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useAutoPickExposure, acceptableExposures, noUsableExposureState,
  resolveActiveExposure,
} from "../src/hooks/useAutoPickExposure";
import { useAppState } from "../src/state";
import type { Exposure } from "../src/api";

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

describe("resolveActiveExposure (render-time exposure resolution)", () => {
  const accepted = EXPOSURES as Exposure[]; // id 101 is the flagged rep

  it("returns undefined while exposures are unknown (still fetching) → skeleton holds", () => {
    expect(resolveActiveExposure(undefined, undefined)).toBeUndefined();
    // a stale id from another sample cannot resolve against unknown exposures
    expect(resolveActiveExposure(999, undefined)).toBeUndefined();
  });

  it("keeps a valid current pick (never clobbers a deliberate switch)", () => {
    expect(resolveActiveExposure(102, accepted)).toBe(102);
  });

  it("adopts the flagged representative for a cold/undefined or stale pick", () => {
    expect(resolveActiveExposure(undefined, accepted)).toBe(101);
    expect(resolveActiveExposure(999, accepted)).toBe(101); // stale → rep
  });

  it("falls back to the first acceptable when none is flagged", () => {
    const noFlag = accepted.map((e) => ({ ...e, selected: false }));
    expect(resolveActiveExposure(undefined, noFlag)).toBe(100);
  });

  it("skips rejected exposures (a pick on a rejected frame re-resolves to the rep)", () => {
    const repRejected = accepted.map((e) =>
      e.id === 102 ? { ...e, status: "rejected" as const } : e,
    );
    expect(resolveActiveExposure(102, repRejected)).toBe(101); // 102 rejected → rep
  });

  it("returns undefined when no exposure is usable (loaded, all rejected/none)", () => {
    expect(resolveActiveExposure(undefined, [])).toBeUndefined();
    const allRejected = accepted.map((e) => ({ ...e, status: "rejected" as const }));
    expect(resolveActiveExposure(101, allRejected)).toBeUndefined();
  });
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

// Pure helpers backing the auto-pick effect AND PlotCard's zero-exposure empty
// state (#24) — one source of truth for "which exposures count".
function exp(id: number, status: "accepted" | "rejected" | null): Exposure {
  return {
    id, sample_id: SAMPLE_ID, filename: `e${id}.tif`, kind: "file",
    selected: false, status, image_path: null, image_version: "",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
  };
}

describe("acceptableExposures", () => {
  it("drops rejected exposures, keeps accepted and null-status", () => {
    const out = acceptableExposures([exp(1, "accepted"), exp(2, "rejected"), exp(3, null)]);
    expect(out.map((e) => e.id)).toEqual([1, 3]);
  });
  it("treats undefined (not-yet-loaded) as empty", () => {
    expect(acceptableExposures(undefined)).toEqual([]);
  });
});

describe("noUsableExposureState", () => {
  it("is not usable-empty while exposures are still loading (undefined)", () => {
    // Guards against flashing the empty state during the fetch.
    expect(noUsableExposureState(undefined)).toEqual({ noUsable: false, allRejected: false });
  });
  it("is usable when at least one exposure is acceptable", () => {
    expect(noUsableExposureState([exp(1, "rejected"), exp(2, "accepted")]))
      .toEqual({ noUsable: false, allRejected: false });
  });
  it("flags all-rejected when the sample has exposures but none acceptable", () => {
    expect(noUsableExposureState([exp(1, "rejected"), exp(2, "rejected")]))
      .toEqual({ noUsable: true, allRejected: true });
  });
  it("flags no-exposures (not all-rejected) when the sample has none at all", () => {
    expect(noUsableExposureState([])).toEqual({ noUsable: true, allRejected: false });
  });
});
