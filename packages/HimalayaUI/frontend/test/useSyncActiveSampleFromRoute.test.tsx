import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import {
  useSyncActiveSampleFromRoute,
  resolveRouteSampleStatus,
} from "../src/hooks/useSyncActiveSampleFromRoute";
import { useAppState } from "../src/state";
import type { CorpusSample } from "../src/api";

// Corpus samples fixture served by the fetch interceptor for GET /api/samples.
const CORPUS = [
  { id: 11, experiment_id: 1, name: "JC011", display_name: "Sample 11",
    notes: null, tags: [], q_units: "A-1" },
  { id: 22, experiment_id: 1, name: "JC022", display_name: "Sample 22",
    notes: null, tags: [], q_units: "A-1" },
];

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    activeSampleId: undefined,
    activeExposureId: undefined,
    username: "tester",
  });
  vi.stubGlobal(
    "fetch",
    vi.fn(async (url: string) => {
      if (String(url).includes("/api/samples")) {
        return new Response(JSON.stringify(CORPUS), {
          status: 200,
          headers: { "content-type": "application/json" },
        });
      }
      return new Response("[]", {
        status: 200,
        headers: { "content-type": "application/json" },
      });
    }),
  );
});

// A throwaway host component that runs the hook and surfaces the returned
// route-resolution status as text so tests can assert it.
function Host(): JSX.Element {
  const status = useSyncActiveSampleFromRoute();
  return <div data-testid="host">{status}</div>;
}

function renderAt(path: string) {
  const qc = new QueryClient({
    defaultOptions: { queries: { retry: false } },
  });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}>
        <Routes>
          <Route path="/sample/:sampleId" element={<Host />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("useSyncActiveSampleFromRoute", () => {
  it("drives activeSampleId from a valid route param", async () => {
    renderAt("/sample/22");
    await waitFor(() =>
      expect(useAppState.getState().activeSampleId).toBe(22),
    );
  });

  it("ignores a route param that names no known sample", async () => {
    renderAt("/sample/999999");
    // host is mounted (RTL retry-aware query) => the effect has run; the
    // corpus fetch is awaited via findBy so the unknown-id branch was reached.
    expect(await screen.findByTestId("host")).toBeInTheDocument();
    expect(useAppState.getState().activeSampleId).toBeUndefined();
  });

  it("ignores a non-numeric route param", async () => {
    renderAt("/sample/not-a-number");
    expect(await screen.findByTestId("host")).toBeInTheDocument();
    expect(useAppState.getState().activeSampleId).toBeUndefined();
  });

  it("does not re-call setActiveSample when the id is unchanged (no exposure clobber)", async () => {
    useAppState.setState({ activeSampleId: 22, activeExposureId: 7 });
    const spy = vi.spyOn(useAppState.getState(), "setActiveSample");
    renderAt("/sample/22");
    expect(await screen.findByTestId("host")).toBeInTheDocument();
    expect(spy).not.toHaveBeenCalled();
    // exposure preserved because setActiveSample (which clears it) never fired
    expect(useAppState.getState().activeExposureId).toBe(7);
  });
});

describe("useSyncActiveSampleFromRoute: route-resolution status (F-STALEURL)", () => {
  it("returns 'pending' while the corpus cache is unresolved", () => {
    // A never-resolving fetch keeps the cache undefined: the hook cannot judge
    // the param yet, so it must not claim found or unknown.
    vi.stubGlobal("fetch", vi.fn(() => new Promise<Response>(() => {})));
    renderAt("/sample/11");
    expect(screen.getByTestId("host")).toHaveTextContent("pending");
  });

  it("returns 'found' once the cache resolves and the param names a known sample", async () => {
    renderAt("/sample/22");
    await waitFor(() =>
      expect(screen.getByTestId("host")).toHaveTextContent("found"),
    );
  });

  it("returns 'unknown' for a numeric param matching no sample", async () => {
    renderAt("/sample/999999");
    await waitFor(() =>
      expect(screen.getByTestId("host")).toHaveTextContent("unknown"),
    );
  });

  it("returns 'unknown' for a non-numeric param once the cache is ready", async () => {
    renderAt("/sample/not-a-number");
    await waitFor(() =>
      expect(screen.getByTestId("host")).toHaveTextContent("unknown"),
    );
  });

  it("'unknown' leaves the previous activeSampleId untouched (mid-session bogus URL)", async () => {
    useAppState.setState({ activeSampleId: 22 });
    renderAt("/sample/999999");
    await waitFor(() =>
      expect(screen.getByTestId("host")).toHaveTextContent("unknown"),
    );
    expect(useAppState.getState().activeSampleId).toBe(22);
  });
});

describe("resolveRouteSampleStatus (shared pure predicate)", () => {
  const samples = CORPUS as CorpusSample[];

  it("returns 'pending' while the samples cache is undefined", () => {
    expect(resolveRouteSampleStatus("22", undefined)).toBe("pending");
    expect(resolveRouteSampleStatus("not-a-number", undefined)).toBe("pending");
  });

  it("returns 'found' for a param naming a known sample", () => {
    expect(resolveRouteSampleStatus("22", samples)).toBe("found");
  });

  it("returns 'unknown' for a numeric param matching no sample", () => {
    expect(resolveRouteSampleStatus("999999", samples)).toBe("unknown");
  });

  it("returns 'unknown' for a non-numeric or missing param once the cache is ready", () => {
    expect(resolveRouteSampleStatus("not-a-number", samples)).toBe("unknown");
    expect(resolveRouteSampleStatus(undefined, samples)).toBe("unknown");
  });
});
