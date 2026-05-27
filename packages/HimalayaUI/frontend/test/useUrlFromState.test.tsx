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
//
// I4.4 (#181) + I1.7 (#163): Index and Inspect are retired. `activePage` can
// only be "compare", for which `buildUrl` returns `current` — so this hook no
// longer emits any Zustand-derived URL. Compare owns its own URL via
// useNavigate/ComparePage. The remaining contract is purely "never tug the
// URL away from a Compare surface, and bail while resolving / stale".

function makeWrapper() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  qc.setQueryData<Experiment[]>(queryKeys.experiments, [
    { id: 17, name: "lipid-screen", path: "", data_dir: "", analysis_dir: "",
      manifest_path: null, created_at: "", q_units: null },
  ]);
  qc.setQueryData<Sample[]>(queryKeys.samples(17), [
    { id: 42, experiment_id: 17, name: "JC001", display_name: null, notes: null, tags: [] },
  ]);
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
    activePage: "compare",
    activeExperimentId: undefined, activeSampleId: undefined, activeExposureId: undefined,
    staleUrlContext: null, resolving: false,
  });
});

describe("useUrlFromState", () => {
  it("does not tug a Compare URL away (buildUrl returns current → no emit)", () => {
    const { Wrapper } = makeWrapper();
    history.replaceState(null, "", "/experiments/17/compare");
    useAppState.setState({ activePage: "compare", activeExperimentId: 17, activeSampleId: 42 });
    const pushSpy = vi.spyOn(history, "pushState");
    const replaceSpy = vi.spyOn(history, "replaceState");
    renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    const hookPushes = pushSpy.mock.calls.filter((c) => typeof c[2] === "string");
    const hookReplaces = replaceSpy.mock.calls.filter((c) => typeof c[2] === "string");
    expect(hookPushes).toHaveLength(0);
    expect(hookReplaces).toHaveLength(0);
    expect(location.pathname).toBe("/experiments/17/compare");
  });

  it("equality guard: identical URL does not emit", () => {
    const { Wrapper } = makeWrapper();
    history.replaceState(null, "", "/compare/all");
    useAppState.setState({ activePage: "compare", activeExperimentId: 17, activeSampleId: 42 });
    const pushSpy = vi.spyOn(history, "pushState");
    const replaceSpy = vi.spyOn(history, "replaceState");
    const { rerender } = renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    pushSpy.mockClear(); replaceSpy.mockClear();
    rerender();
    expect(pushSpy).not.toHaveBeenCalled();
    expect(replaceSpy).not.toHaveBeenCalled();
  });

  it("bails (no emit) while resolving", () => {
    const { Wrapper } = makeWrapper();
    history.replaceState(null, "", "/compare/all");
    useAppState.setState({
      activePage: "compare", activeExperimentId: 17, activeSampleId: 42, resolving: true,
    });
    const pushSpy = vi.spyOn(history, "pushState");
    const replaceSpy = vi.spyOn(history, "replaceState");
    renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    const hookPushes = pushSpy.mock.calls.filter((c) => typeof c[2] === "string");
    const hookReplaces = replaceSpy.mock.calls.filter((c) => typeof c[2] === "string");
    expect(hookPushes).toHaveLength(0);
    expect(hookReplaces).toHaveLength(0);
  });

  it("bails (no emit) while a stale URL context is parked", () => {
    const { Wrapper } = makeWrapper();
    history.replaceState(null, "", "/foo/bar");
    useAppState.setState({
      activePage: "compare", activeExperimentId: 17, activeSampleId: 42,
      staleUrlContext: { kind: "unknown_path", raw: "/foo/bar" },
    });
    const pushSpy = vi.spyOn(history, "pushState");
    const replaceSpy = vi.spyOn(history, "replaceState");
    renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    const hookPushes = pushSpy.mock.calls.filter((c) => typeof c[2] === "string");
    const hookReplaces = replaceSpy.mock.calls.filter((c) => typeof c[2] === "string");
    expect(hookPushes).toHaveLength(0);
    expect(hookReplaces).toHaveLength(0);
  });
});
