import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { FocusDetectorPanel } from "../src/components/FocusDetectorPanel";
import { useAppState } from "../src/state";

const EXPOSURES = [
  { id: 5, sample_id: 1, filename: "JC001-005.dat", kind: "file" as const,
    selected: true, status: "accepted" as const,
    image_path: "/img/5.png", image_version: "v1", tags: [], sources: [],
    trace_hash: null, analysis_inputs_hash: null },
];

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    activeSampleId: 1, activeExposureId: 5,
    username: "tester", theme: "dark",
  });
  vi.stubGlobal("ResizeObserver", class {
    observe() {} unobserve() {} disconnect() {}
  });
  vi.stubGlobal("fetch", vi.fn(async (url: string) => {
    if (String(url).includes("/exposures")) {
      return new Response(JSON.stringify(EXPOSURES), {
        status: 200, headers: { "content-type": "application/json" },
      });
    }
    return new Response("[]", { status: 200,
      headers: { "content-type": "application/json" } });
  }));
});

function renderPanel() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <FocusDetectorPanel />
    </QueryClientProvider>,
  );
}

describe("FocusDetectorPanel", () => {
  it("renders the detector panel region", async () => {
    renderPanel();
    expect(await screen.findByTestId("focus-detector-panel")).toBeInTheDocument();
  });

  it("shows a pick-a-sample hint when no sample is active", async () => {
    useAppState.setState({ activeSampleId: undefined, activeExposureId: undefined });
    renderPanel();
    expect(await screen.findByTestId("focus-detector-empty")).toBeInTheDocument();
  });
});
