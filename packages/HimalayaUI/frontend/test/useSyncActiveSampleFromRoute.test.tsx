import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useSyncActiveSampleFromRoute } from "../src/hooks/useSyncActiveSampleFromRoute";
import { useAppState } from "../src/state";

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
    theme: "dark",
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

// A throwaway host component that just runs the hook.
function Host(): JSX.Element {
  useSyncActiveSampleFromRoute();
  return <div data-testid="host" />;
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
