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
const PEAKS = [
  { id: 1, exposure_id: 5, q: 0.045, intensity: 100, prominence: 50,
    sharpness: 1.2, source: "auto" as const, excluded: false },
  { id: 2, exposure_id: 5, q: 0.103, intensity: 80, prominence: 40,
    sharpness: 1.1, source: "auto" as const, excluded: false },
];

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    activeSampleId: 1, activeExposureId: 5,
    username: "tester",
  });
  vi.stubGlobal("ResizeObserver", class {
    observe() {} unobserve() {} disconnect() {}
  });
  vi.stubGlobal("fetch", vi.fn(async (url: string) => {
    const u = String(url);
    const json = (b: unknown) => new Response(JSON.stringify(b), {
      status: 200, headers: { "content-type": "application/json" } });
    // The detector PNG route: return non-ok so DetectorImage.renderImage bails
    // at `if (!res.ok) return;` BEFORE calling createImageBitmap (absent in
    // JSDOM). We assert the ring overlay, not the rendered canvas pixels.
    if (u.includes("/image")) return new Response(null, { status: 404 });
    // /peaks (/api/exposures/:id/peaks) before the looser /exposures check.
    if (u.includes("/peaks")) return json(PEAKS);
    if (u.includes("/exposures")) return json(EXPOSURES);
    return json([]);
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

  it("renders the q-link ring overlay over the detector image", async () => {
    renderPanel();
    expect(await screen.findByTestId("detector-ring-overlay")).toBeInTheDocument();
  });
});
