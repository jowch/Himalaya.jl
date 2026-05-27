import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { BrowserRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { useStateFromUrl } from "../src/hooks/useStateFromUrl";
import { useAppState } from "../src/state";
import { _resetEmitMode } from "../src/lib/url/emitMode";

// BrowserRouter reads from `window.history` / `window.location`, which is
// what the tests manipulate via `history.replaceState`.
//
// I4.4 (#181): the Index surface is retired. `parseLocation` no longer returns
// a slug-bearing `index` kind, so useStateFromUrl no longer resolves
// `/index/<exp>/<sample>` URLs into Zustand active ids (the router redirects
// them via IndexSlugRedirect). The surviving behaviors:
//   - `compare`  → activePage="compare", resolving cleared.
//   - `root` (/) → redirect to /samples (§4.1 corpus home).
//   - `stale`    → setStaleUnknownPath.
function wrapper({ children }: { children: ReactNode }) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return (
    <QueryClientProvider client={qc}>
      <BrowserRouter>{children}</BrowserRouter>
    </QueryClientProvider>
  );
}

// Spec §4.2

beforeEach(() => {
  useAppState.setState({
    staleUrlContext: null, resolving: false,
    activeExperimentId: undefined, activeSampleId: undefined, activeExposureId: undefined,
    activePage: "compare",
  });
  history.replaceState(null, "", "/");
  _resetEmitMode();
});

afterEach(() => {
  vi.restoreAllMocks();
});

describe("useStateFromUrl", () => {
  it("on mount with /<unknown>/path: sets staleUrlContext kind=unknown_path, no fetch", () => {
    history.replaceState(null, "", "/foo/bar");
    const fetchSpy = vi.spyOn(global, "fetch");
    renderHook(() => useStateFromUrl(), { wrapper });
    const ctx = useAppState.getState().staleUrlContext;
    expect(ctx?.kind).toBe("unknown_path");
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it("I4.4: a stray /index path is now a stale unknown_path (no resolve)", () => {
    // /index* is redirected at the router; if the parser ever sees it, it is
    // classified stale rather than triggering a resolve fetch.
    history.replaceState(null, "", "/index/lipid/JC001");
    const fetchSpy = vi.spyOn(global, "fetch");
    renderHook(() => useStateFromUrl(), { wrapper });
    expect(useAppState.getState().staleUrlContext?.kind).toBe("unknown_path");
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it("compare URL: sets activePage='compare', clears resolving, no fetch", () => {
    history.replaceState(null, "", "/compare/all");
    useAppState.setState({ resolving: true });
    const fetchSpy = vi.spyOn(global, "fetch");
    renderHook(() => useStateFromUrl(), { wrapper });
    expect(useAppState.getState().activePage).toBe("compare");
    expect(useAppState.getState().resolving).toBe(false);
    expect(fetchSpy).not.toHaveBeenCalled();
  });
});

describe("useStateFromUrl — / redirect (§4.1 corpus home)", () => {
  // The redirect path uses `useNavigate()` which under BrowserRouter calls
  // `window.history.replaceState`. We assert on the URL arg (`mock.calls[i][2]`)
  // because BrowserRouter passes a state object (not raw `null`) as the first arg.
  function makeRedirectWrapper(qc: QueryClient) {
    return function RedirectWrapper({ children }: { children: ReactNode }) {
      return (
        <QueryClientProvider client={qc}>
          <BrowserRouter>{children}</BrowserRouter>
        </QueryClientProvider>
      );
    };
  }

  it("bare / redirects to /samples regardless of persisted state, no fetch", async () => {
    history.replaceState(null, "", "/");
    useAppState.setState({
      activePage: "compare", activeExperimentId: 17, activeSampleId: 42,
    });
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    const Wrapper = makeRedirectWrapper(qc);
    const fetchSpy = vi.spyOn(global, "fetch");
    const replaceSpy = vi.spyOn(history, "replaceState");
    renderHook(() => useStateFromUrl(), { wrapper: Wrapper });
    await waitFor(() => {
      expect(replaceSpy.mock.calls.some((c) => c[2] === "/samples")).toBe(true);
    });
    expect(fetchSpy).not.toHaveBeenCalled();
  });
});
