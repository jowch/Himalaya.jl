import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook } from "@testing-library/react";
import { BrowserRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { useUrlFromState } from "../src/hooks/useUrlFromState";
import { useAppState } from "../src/state";
import { queryKeys } from "../src/queries";
import { _resetEmitMode } from "../src/lib/url/emitMode";
import type { Experiment, Sample } from "../src/api";

// Spec §4.3

function makeWrapper() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  qc.setQueryData<Experiment[]>(queryKeys.experiments, [
    { id: 17, name: "lipid-screen", path: "", data_dir: "", analysis_dir: "",
      manifest_path: null, created_at: "", q_units: null },
  ]);
  qc.setQueryData<Sample[]>(queryKeys.samples(17), [
    { id: 42, experiment_id: 17, name: "JC001", display_name: null, notes: null, tags: [] },
  ]);
  // BrowserRouter reads from `window.history` / `window.location`, which is
  // what the tests manipulate via `history.replaceState`. Mirrors the
  // wrapper used in `useStateFromUrl.test.tsx`.
  const Wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={qc}>
      <BrowserRouter>{children}</BrowserRouter>
    </QueryClientProvider>
  );
  return { qc, Wrapper };
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

describe("useUrlFromState", () => {
  it("active sample → /index/<exp>/<sample>", () => {
    const { Wrapper } = makeWrapper();
    useAppState.setState({ activePage: "index", activeExperimentId: 17, activeSampleId: 42 });
    renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    expect(location.pathname).toBe("/index/lipid-screen/JC001");
  });

  it("page change → push", () => {
    const { Wrapper } = makeWrapper();
    history.replaceState(null, "", "/index/lipid-screen/JC001");
    useAppState.setState({ activePage: "index", activeExperimentId: 17, activeSampleId: 42 });
    const pushSpy = vi.spyOn(history, "pushState");
    const replaceSpy = vi.spyOn(history, "replaceState");
    const { rerender } = renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    pushSpy.mockClear(); replaceSpy.mockClear();
    useAppState.setState({ activePage: "inspect" });
    rerender();
    // react-router's BrowserRouter calls pushState with a state object (not
    // raw null) and the resolved URL as the third arg. Assert push (not
    // replace) was invoked with the expected URL.
    expect(pushSpy).toHaveBeenCalled();
    expect(replaceSpy).not.toHaveBeenCalled();
    expect(pushSpy.mock.calls[0]?.[2]).toBe("/inspect/lipid-screen/JC001");
  });

  it("exposure-only change → replace", () => {
    const { qc, Wrapper } = makeWrapper();
    qc.setQueryData(queryKeys.exposures(42), [
      { id: 100, sample_id: 42, filename: "JC001-007", kind: "simple", selected: 1 },
    ]);
    history.replaceState(null, "", "/inspect/lipid-screen/JC001");
    useAppState.setState({ activePage: "inspect", activeExperimentId: 17, activeSampleId: 42 });
    const pushSpy = vi.spyOn(history, "pushState");
    const replaceSpy = vi.spyOn(history, "replaceState");
    const { rerender } = renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    pushSpy.mockClear(); replaceSpy.mockClear();
    // Use the action (not raw setState) so emitReplaceNext fires; that's
    // the contract — exposure-only changes that route through
    // setActiveExposure must replace, not push.
    useAppState.getState().setActiveExposure(100);
    rerender();
    expect(replaceSpy).toHaveBeenCalled();
    expect(pushSpy).not.toHaveBeenCalled();
    expect(location.search).toBe("?exposure=JC001-007");
  });

  it("equality guard: identical URL does not emit", () => {
    const { Wrapper } = makeWrapper();
    history.replaceState(null, "", "/index/lipid-screen/JC001");
    useAppState.setState({ activePage: "index", activeExperimentId: 17, activeSampleId: 42 });
    const pushSpy = vi.spyOn(history, "pushState");
    const replaceSpy = vi.spyOn(history, "replaceState");
    const { rerender } = renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    pushSpy.mockClear(); replaceSpy.mockClear();
    rerender();
    expect(pushSpy).not.toHaveBeenCalled();
    expect(replaceSpy).not.toHaveBeenCalled();
  });

  it("replay-as-rerun: identical optimistic + confirmed slug → no spurious emit", () => {
    // Simulate the trivial replay case: cache row gets replaced (foreign event)
    // but the same id-name mapping holds. URL recompute should see the same
    // slug and not emit.
    const { qc, Wrapper } = makeWrapper();
    history.replaceState(null, "", "/index/lipid-screen/JC001");
    useAppState.setState({ activePage: "index", activeExperimentId: 17, activeSampleId: 42 });
    const replaceSpy = vi.spyOn(history, "replaceState");
    const { rerender } = renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    replaceSpy.mockClear();
    // Simulate applyRemoteToCache rewriting samples in place.
    qc.setQueryData<Sample[]>(queryKeys.samples(17), [
      { id: 42, experiment_id: 17, name: "JC001", display_name: "JC001 (touched)",
        notes: null, tags: [] },
    ]);
    rerender();
    expect(replaceSpy).not.toHaveBeenCalled();
  });
});
