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

  // I1.7 (#163): the "page change → push (/inspect)" and "exposure-only
  // change → replace (?exposure=)" tests are retired with Inspect — it was
  // the only legacy tab that emitted a navigable slug distinct from Index and
  // the only surface whose URL carried ?exposure=.

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

  it("does not emit while experiments cache is unhydrated (deep-link race)", () => {
    // Cold-mount race: applySuccess populated activeExperimentId / activeSampleId
    // before useExperiments() finished. If the hook emitted now, buildUrl
    // would resolve to /index (no slugs) and useStateFromUrl would wipe the
    // just-populated active ids. The cache-hydration gate prevents that.
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    // Deliberately do NOT setQueryData — simulate cache loading.
    const Wrapper = ({ children }: { children: ReactNode }) => (
      <QueryClientProvider client={qc}>
        <BrowserRouter>{children}</BrowserRouter>
      </QueryClientProvider>
    );
    history.replaceState(null, "", "/index/lipid-screen/JC001");
    useAppState.setState({
      activePage: "index", activeExperimentId: 17, activeSampleId: 42,
    });
    const pushSpy = vi.spyOn(history, "pushState");
    const replaceSpy = vi.spyOn(history, "replaceState");
    renderHook(() => useUrlFromState(), { wrapper: Wrapper });
    // Effect runs but should bail because experiments cache is empty.
    // BrowserRouter's mount calls history.replaceState(state, "") with two
    // args to inject its own internal `{ idx: 0 }` marker — filter for
    // hook-driven navigations (which pass a URL string as the 3rd arg).
    const hookPushes = pushSpy.mock.calls.filter((c) => typeof c[2] === "string");
    const hookReplaces = replaceSpy.mock.calls.filter((c) => typeof c[2] === "string");
    expect(hookPushes).toHaveLength(0);
    expect(hookReplaces).toHaveLength(0);
    // URL must remain the deep link — the bug would have rewritten it to /index.
    expect(location.pathname).toBe("/index/lipid-screen/JC001");
  });

  // I1.7 (#163): the exposures-cache-unhydrated deep-link race test is retired
  // with Inspect — no URL carries ?exposure= and useUrlFromState no longer
  // subscribes to the exposures cache.

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
