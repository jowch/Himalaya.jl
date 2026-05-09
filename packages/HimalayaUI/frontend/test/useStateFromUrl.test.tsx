import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { BrowserRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { useStateFromUrl } from "../src/hooks/useStateFromUrl";
import { useAppState } from "../src/state";
import type { ResolveSuccess, ResolveError404 } from "../src/api";
import { queryKeys } from "../src/queries";
import { _resetEmitMode } from "../src/lib/url/emitMode";

// BrowserRouter reads from `window.history` / `window.location`, which is
// what the tests manipulate via `history.replaceState`. A QueryClient
// wrapper is required because the hook calls `useQueryClient()` for the
// root-redirect cache-hit fast path.
function wrapper({ children }: { children: ReactNode }) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return (
    <QueryClientProvider client={qc}>
      <BrowserRouter>{children}</BrowserRouter>
    </QueryClientProvider>
  );
}

// Spec §4.2

const ok = (body: ResolveSuccess) => Promise.resolve({
  ok: true, status: 200, json: () => Promise.resolve(body),
} as Response);

const notFound = (body: ResolveError404) => Promise.resolve({
  ok: false, status: 404, json: () => Promise.resolve(body),
} as Response);

beforeEach(() => {
  useAppState.setState({
    staleUrlContext: null, resolving: false,
    activeExperimentId: undefined, activeSampleId: undefined, activeExposureId: undefined,
    activePage: "index",
  });
  history.replaceState(null, "", "/");
  _resetEmitMode();
});

afterEach(() => {
  vi.restoreAllMocks();
});

describe("useStateFromUrl", () => {
  it("on mount with /index/<exp>/<sample>: dispatches setActive*, clears stale", async () => {
    history.replaceState(null, "", "/index/lipid/JC001");
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(() =>
      ok({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    );
    renderHook(() => useStateFromUrl(), { wrapper });
    await waitFor(() => {
      expect(useAppState.getState().activeExperimentId).toBe(17);
    });
    expect(useAppState.getState().activeSampleId).toBe(42);
    expect(useAppState.getState().activePage).toBe("index");
    expect(useAppState.getState().resolving).toBe(false);
    expect(fetchSpy).toHaveBeenCalledOnce();
  });

  it("on mount with /<unknown>/path: sets staleUrlContext kind=unknown_path, no fetch", () => {
    history.replaceState(null, "", "/foo/bar");
    const fetchSpy = vi.spyOn(global, "fetch");
    renderHook(() => useStateFromUrl(), { wrapper });
    const ctx = useAppState.getState().staleUrlContext;
    expect(ctx?.kind).toBe("unknown_path");
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it("404: sets staleUrlContext from response body", async () => {
    history.replaceState(null, "", "/index/lipid-typo");
    vi.spyOn(global, "fetch").mockImplementation(() =>
      notFound({
        error: "not_found", missing: "experiment", missing_value: "lipid-typo",
        experiment_resolved: undefined, sample_resolved: undefined,
      }),
    );
    renderHook(() => useStateFromUrl(), { wrapper });
    await waitFor(() => {
      const ctx = useAppState.getState().staleUrlContext;
      expect(ctx?.kind).toBe("not_found");
    });
  });

  it("origin-tag race: navigation mid-resolve discards the response", async () => {
    history.replaceState(null, "", "/index/lipid/JC001");
    let resolveFetch: ((v: Response) => void) | undefined;
    vi.spyOn(global, "fetch").mockImplementation(() =>
      new Promise<Response>((res) => { resolveFetch = res; }),
    );
    renderHook(() => useStateFromUrl(), { wrapper });
    expect(useAppState.getState().resolving).toBe(true);

    // Simulate user navigating mid-flight via TabRocker (Zustand mutation, no popstate).
    history.replaceState(null, "", "/inspect/other/SAMPLE");
    // Now satisfy the original fetch — its origin tag should mismatch.
    resolveFetch?.({
      ok: true, status: 200,
      json: () => Promise.resolve({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    } as Response);
    await new Promise((r) => setTimeout(r, 0));
    // Active state should NOT have been updated by the discarded resolve.
    expect(useAppState.getState().activeExperimentId).toBeUndefined();
  });

  it("pre-fetch clear: staleUrlContext is cleared at start of recognized fetch", async () => {
    useAppState.setState({
      staleUrlContext: { kind: "unknown_path", raw: "/old" },
    });
    history.replaceState(null, "", "/index/lipid/JC001");
    vi.spyOn(global, "fetch").mockImplementation(() =>
      ok({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    );
    renderHook(() => useStateFromUrl(), { wrapper });
    // The clear happens synchronously before fetch; assert without waiting.
    expect(useAppState.getState().staleUrlContext).toBeNull();
  });
});

// Issue #114 — slug-equality fast path. When the parsed slugs already
// match what the current activeIds + TanStack cache would resolve to,
// the URL is just a stable representation of state we already have.
// Skip the /api/resolve round-trip (and the resolving:true that swaps
// PageBody to ResolvingFallback, which is what visually flashes).
describe("useStateFromUrl — slug-equality fast path", () => {
  function makeFastPathWrapper(qc: QueryClient) {
    return function FastPathWrapper({ children }: { children: ReactNode }) {
      return (
        <QueryClientProvider client={qc}>
          <BrowserRouter>{children}</BrowserRouter>
        </QueryClientProvider>
      );
    };
  }

  function seedExperimentSample(qc: QueryClient) {
    qc.setQueryData(queryKeys.experiments, [
      { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "",
        manifest_path: null, created_at: "", q_units: null },
    ]);
    qc.setQueryData(queryKeys.samples(17), [
      { id: 42, experiment_id: 17, name: "JC001", display_name: null, notes: null, tags: [] },
    ]);
  }

  function seedExposure(qc: QueryClient) {
    qc.setQueryData(queryKeys.exposures(42), [
      { id: 99, sample_id: 42, filename: "img001.tif", kind: "file",
        selected: true, status: null, image_path: null, image_version: "",
        tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null },
    ]);
  }

  it("Index URL: slugs match cache+activeIds → no fetch, resolving stays false", () => {
    history.replaceState(null, "", "/index/lipid/JC001");
    useAppState.setState({
      activePage: "index", activeExperimentId: 17, activeSampleId: 42, activeExposureId: undefined,
    });
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    seedExperimentSample(qc);
    const fetchSpy = vi.spyOn(global, "fetch");
    renderHook(() => useStateFromUrl(), { wrapper: makeFastPathWrapper(qc) });
    expect(fetchSpy).not.toHaveBeenCalled();
    expect(useAppState.getState().resolving).toBe(false);
    expect(useAppState.getState().activeExperimentId).toBe(17);
    expect(useAppState.getState().activeSampleId).toBe(42);
  });

  it("Inspect URL: slugs+exposure match cache+activeIds → no fetch", () => {
    history.replaceState(null, "", "/inspect/lipid/JC001?exposure=img001.tif");
    useAppState.setState({
      activePage: "inspect", activeExperimentId: 17, activeSampleId: 42, activeExposureId: 99,
    });
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    seedExperimentSample(qc);
    seedExposure(qc);
    const fetchSpy = vi.spyOn(global, "fetch");
    renderHook(() => useStateFromUrl(), { wrapper: makeFastPathWrapper(qc) });
    expect(fetchSpy).not.toHaveBeenCalled();
    expect(useAppState.getState().resolving).toBe(false);
    expect(useAppState.getState().activeExposureId).toBe(99);
  });

  it("URL kind differs from activePage but slugs match → no fetch, activePage flips", () => {
    history.replaceState(null, "", "/inspect/lipid/JC001");
    useAppState.setState({
      activePage: "index", activeExperimentId: 17, activeSampleId: 42, activeExposureId: undefined,
    });
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    seedExperimentSample(qc);
    const fetchSpy = vi.spyOn(global, "fetch");
    renderHook(() => useStateFromUrl(), { wrapper: makeFastPathWrapper(qc) });
    expect(fetchSpy).not.toHaveBeenCalled();
    expect(useAppState.getState().activePage).toBe("inspect");
  });

  it("cold mount with empty cache → still fetches /api/resolve", async () => {
    history.replaceState(null, "", "/index/lipid/JC001");
    useAppState.setState({
      activePage: "index", activeExperimentId: 17, activeSampleId: 42, activeExposureId: undefined,
    });
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    // Intentionally NO seedExperimentSample(qc) — cache is empty (cold mount).
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(() =>
      ok({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    );
    renderHook(() => useStateFromUrl(), { wrapper: makeFastPathWrapper(qc) });
    await waitFor(() => {
      expect(fetchSpy).toHaveBeenCalledOnce();
    });
  });

  it("popstate to a new sample slug (cache lookup wrong name) → still fetches", async () => {
    // Zustand still on old sample id 42 (name "JC001"), URL points at "JC002".
    // The cache row for 42 disagrees with the URL slug → can't fast-path.
    history.replaceState(null, "", "/index/lipid/JC002");
    useAppState.setState({
      activePage: "index", activeExperimentId: 17, activeSampleId: 42, activeExposureId: undefined,
    });
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    seedExperimentSample(qc);
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(() =>
      ok({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 43, sample_name: "JC002",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    );
    renderHook(() => useStateFromUrl(), { wrapper: makeFastPathWrapper(qc) });
    await waitFor(() => {
      expect(fetchSpy).toHaveBeenCalledOnce();
    });
  });

  it("activeIds undefined despite URL slugs (eg. fresh tab on deep link) → still fetches", async () => {
    history.replaceState(null, "", "/index/lipid/JC001");
    useAppState.setState({
      activePage: "index", activeExperimentId: undefined, activeSampleId: undefined, activeExposureId: undefined,
    });
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    seedExperimentSample(qc);  // cache present but activeIds aren't
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(() =>
      ok({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    );
    renderHook(() => useStateFromUrl(), { wrapper: makeFastPathWrapper(qc) });
    await waitFor(() => {
      expect(fetchSpy).toHaveBeenCalledOnce();
    });
  });
});

describe("useStateFromUrl — / redirect", () => {
  // The redirect path uses `useNavigate()` which under BrowserRouter calls
  // `window.history.replaceState`. We assert on the URL arg only (`mock.calls[i][2]`)
  // because BrowserRouter passes a state object (not raw `null`) as the first arg.
  // Wrapper must include <BrowserRouter> (router context) AND <QueryClientProvider>
  // (the redirect's synchronous getQueryData fast path).
  function makeRedirectWrapper(qc: QueryClient) {
    return function RedirectWrapper({ children }: { children: ReactNode }) {
      return (
        <QueryClientProvider client={qc}>
          <BrowserRouter>{children}</BrowserRouter>
        </QueryClientProvider>
      );
    };
  }

  it("with persisted experiment+sample: synchronous getQueryData → replaceState to slug URL", () => {
    history.replaceState(null, "", "/");
    useAppState.setState({
      activePage: "index", activeExperimentId: 17, activeSampleId: 42,
    });
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    qc.setQueryData(queryKeys.experiments, [
      { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "",
        manifest_path: null, created_at: "", q_units: null },
    ]);
    qc.setQueryData(queryKeys.samples(17), [
      { id: 42, experiment_id: 17, name: "JC001", display_name: null, notes: null, tags: [] },
    ]);
    const Wrapper = makeRedirectWrapper(qc);
    const replaceSpy = vi.spyOn(history, "replaceState");
    renderHook(() => useStateFromUrl(), { wrapper: Wrapper });
    // Expect a synchronous replaceState to /index/lipid/JC001 (no fetch).
    expect(
      replaceSpy.mock.calls.some((c) => c[2] === "/index/lipid/JC001"),
    ).toBe(true);
  });

  it("with persisted ids missing from cache: falls back to /api/resolve?experiment_id=…", async () => {
    history.replaceState(null, "", "/");
    useAppState.setState({
      activePage: "index", activeExperimentId: 17, activeSampleId: 42,
    });
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    const Wrapper = makeRedirectWrapper(qc);
    vi.spyOn(global, "fetch").mockImplementation(() =>
      ok({
        experiment_id: 17, experiment_name: "lipid",
        sample_id: 42, sample_name: "JC001",
        exposure_id: undefined, exposure_filename: undefined,
      }),
    );
    const replaceSpy = vi.spyOn(history, "replaceState");
    renderHook(() => useStateFromUrl(), { wrapper: Wrapper });
    await waitFor(() => {
      expect(
        replaceSpy.mock.calls.some((c) => c[2] === "/index/lipid/JC001"),
      ).toBe(true);
    });
  });

  it("with no persisted state: replaceState to /index", () => {
    history.replaceState(null, "", "/");
    useAppState.setState({
      activePage: "index", activeExperimentId: undefined, activeSampleId: undefined,
    });
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    const Wrapper = makeRedirectWrapper(qc);
    const replaceSpy = vi.spyOn(history, "replaceState");
    renderHook(() => useStateFromUrl(), { wrapper: Wrapper });
    expect(
      replaceSpy.mock.calls.some((c) => c[2] === "/index"),
    ).toBe(true);
  });
});
