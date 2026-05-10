// Regression test for issue #118 — `,`/`.` sample-step keypress.
//
// Two bugs interacted to produce a doubled history op (replace + push to the
// same URL) on backward `,` presses:
//
// 1. `useLocation()` lags `window.location` by one render when Zustand's
//    `useSyncExternalStore`-driven re-render fires before BrowserRouter's
//    `React.useState` commit. The equality guard in `useUrlFromState`
//    compared against the stale react-router `location`, so a second
//    navigate fired even though the URL had already been written.
//
// 2. PlotCard's auto-pick effect calls `setActiveExposure(first.id)` after
//    a sample switch. `setActiveExposure` armed `emitReplaceNext()`
//    unconditionally, clobbering the intended PUSH mode for the sample
//    switch — even on the Index page where the URL doesn't include the
//    exposure at all. Result: rapid `,` presses didn't accumulate history
//    entries, breaking back-button.
//
// Reproduction needs the auto-pick effect to fire *in the same commit
// cycle* as the keypress's URL emit. That happens when the destination
// sample's exposures are already cached — typical for backward `,`
// (revisiting a previous sample) but not forward `.` (visiting a new one).

import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, act } from "@testing-library/react";
import { BrowserRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useEffect, type ReactNode } from "react";
import { useUrlFromState } from "../src/hooks/useUrlFromState";
import { useStateFromUrl } from "../src/hooks/useStateFromUrl";
import { useAppState } from "../src/state";
import { useExposures, queryKeys } from "../src/queries";
import { _resetEmitMode } from "../src/lib/url/emitMode";
import type { Experiment, Sample, Exposure } from "../src/api";

function makeWrapper() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  qc.setQueryData<Experiment[]>(queryKeys.experiments, [
    { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "",
      manifest_path: null, created_at: "", q_units: null },
  ]);
  qc.setQueryData<Sample[]>(queryKeys.samples(17), [
    { id: 41, experiment_id: 17, name: "C1", display_name: null, notes: null, tags: [] },
    { id: 42, experiment_id: 17, name: "C2", display_name: null, notes: null, tags: [] },
    { id: 43, experiment_id: 17, name: "C3", display_name: null, notes: null, tags: [] },
  ]);
  // Pre-cache exposures for both C2 and C3 — backward-step path simulates
  // returning to a previously-visited sample whose exposures are still in
  // TanStack's cache.
  qc.setQueryData<Exposure[]>(queryKeys.exposures(42), [
    { id: 200, sample_id: 42, filename: "C2-001", kind: "simple", selected: 0,
      status: "accepted", reject_reason: null },
  ] as Exposure[]);
  qc.setQueryData<Exposure[]>(queryKeys.exposures(43), [
    { id: 300, sample_id: 43, filename: "C3-001", kind: "simple", selected: 0,
      status: "accepted", reject_reason: null },
  ] as Exposure[]);

  const Wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={qc}>
      <BrowserRouter>{children}</BrowserRouter>
    </QueryClientProvider>
  );
  return { qc, Wrapper };
}

// Mirrors PlotCard's auto-pick effect (PlotCard.tsx:106-119). Runs as a
// child useEffect so it fires within the same commit phase as
// useUrlFromState's effect — that interleaving is what reproduces #118.
function FakePlotCard(): null {
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const activeExposureId = useAppState((s) => s.activeExposureId);
  const setActiveExposure = useAppState((s) => s.setActiveExposure);
  const exposuresQ = useExposures(activeSampleId);
  useEffect(() => {
    const exposures = (exposuresQ.data ?? []).filter((e) => e.status !== "rejected");
    if (exposures.length === 0) return;
    const flagged = exposures.find((e) => e.selected);
    if (flagged) {
      if (activeExposureId !== flagged.id) setActiveExposure(flagged.id);
      return;
    }
    const stillValid = exposures.some((e) => e.id === activeExposureId);
    if (!stillValid) setActiveExposure(exposures[0]!.id);
  }, [exposuresQ.data, activeExposureId, setActiveExposure]);
  return null;
}

function App(): JSX.Element {
  useStateFromUrl();
  useUrlFromState();
  return <FakePlotCard />;
}

beforeEach(() => {
  history.replaceState(null, "", "/");
  _resetEmitMode();
  useAppState.setState({
    activePage: "index",
    activeExperimentId: undefined, activeSampleId: undefined, activeExposureId: undefined,
    staleUrlContext: null, resolving: false,
  });
});

function hookOps(pushSpy: ReturnType<typeof vi.spyOn>,
                 replaceSpy: ReturnType<typeof vi.spyOn>) {
  // Filter out BrowserRouter's internal `{ idx: N }` markers (state-only,
  // no URL string in the third arg).
  const pushes = pushSpy.mock.calls
    .filter((c) => typeof c[2] === "string")
    .map((c) => ({ kind: "push", url: c[2] as string }));
  const replaces = replaceSpy.mock.calls
    .filter((c) => typeof c[2] === "string")
    .map((c) => ({ kind: "replace", url: c[2] as string }));
  return [...replaces, ...pushes];
}

describe("issue #118 — , sample-step keypress on Index page", () => {
  it("backward step with cached exposures: emits exactly one PUSH (no double-emit)", () => {
    const { Wrapper } = makeWrapper();
    history.replaceState(null, "", "/index/lipid/C3");
    useAppState.setState({
      activePage: "index", activeExperimentId: 17, activeSampleId: 43,
      activeExposureId: 300,
    });

    render(<App />, { wrapper: Wrapper });

    const pushSpy = vi.spyOn(history, "pushState");
    const replaceSpy = vi.spyOn(history, "replaceState");
    pushSpy.mockClear(); replaceSpy.mockClear();

    // Simulate the `,` keypress (useGlobalShortcuts.ts:77).
    act(() => { useAppState.getState().setActiveSample(42); });

    const ops = hookOps(pushSpy, replaceSpy);
    expect(ops).toHaveLength(1);
    // Must be PUSH so back-button can undo each `,` press.
    expect(ops[0]).toMatchObject({ kind: "push", url: "/index/lipid/C2" });
    expect(location.pathname).toBe("/index/lipid/C2");
  });

  it("forward step on uncached sample: still emits exactly one PUSH (regression baseline)", () => {
    const { qc, Wrapper } = makeWrapper();
    // Drop C2's exposure cache to simulate a forward `.` step into a
    // not-yet-visited sample.
    qc.removeQueries({ queryKey: queryKeys.exposures(42) });
    history.replaceState(null, "", "/index/lipid/C1");
    useAppState.setState({
      activePage: "index", activeExperimentId: 17, activeSampleId: 41,
      activeExposureId: undefined,
    });

    render(<App />, { wrapper: Wrapper });

    const pushSpy = vi.spyOn(history, "pushState");
    const replaceSpy = vi.spyOn(history, "replaceState");
    pushSpy.mockClear(); replaceSpy.mockClear();

    act(() => { useAppState.getState().setActiveSample(42); });

    const ops = hookOps(pushSpy, replaceSpy);
    expect(ops).toHaveLength(1);
    expect(ops[0]).toMatchObject({ kind: "push", url: "/index/lipid/C2" });
  });
});
