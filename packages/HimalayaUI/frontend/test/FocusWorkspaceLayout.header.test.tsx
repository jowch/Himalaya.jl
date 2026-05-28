import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../src/components/AppRoutes";
import { useAppState } from "../src/state";

// Corpus list — how the layout learns experiment_id from the route sampleId.
const CORPUS = [
  { id: 1, experiment_id: 7, name: "smp_09", display_name: "Lipid 1-2 + LL37 1:0.5",
    notes: null, tags: [], q_units: "A-1" },
];
const EXP_SAMPLES = [
  { id: 1, experiment_id: 7, name: "smp_09", display_name: "Lipid 1-2 + LL37 1:0.5",
    notes: null, tags: [] },
];
const EXPERIMENT = {
  id: 7, name: "SSRL Apr 2026", path: "/x", data_dir: "/x/d",
  analysis_dir: "/x/a", manifest_path: null, created_at: "2026-04-01",
  q_units: "A-1",
};
const EXPOSURES = [
  { id: 5, sample_id: 1, filename: "smp_09_e03.dat", kind: "file",
    selected: true, status: "accepted", image_path: "/img/5.png",
    image_version: "v1", tags: [], sources: [],
    trace_hash: null, analysis_inputs_hash: null },
];
const TRACE = { q: [0.04, 0.1, 0.3], I: [10, 80, 20], sigma: [1, 1, 1] };

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    activeSampleId: undefined, activeExposureId: undefined,
    activeExperimentId: undefined, username: "tester",
  });
  vi.stubGlobal("ResizeObserver", class {
    observe() {} unobserve() {} disconnect() {}
  });
  vi.stubGlobal("fetch", vi.fn(async (url: string) => {
    const u = String(url);
    const json = (b: unknown) => new Response(JSON.stringify(b), {
      status: 200, headers: { "content-type": "application/json" } });
    // DetectorImage (co-resident in FocusDetectorPanel) bails at `if (!res.ok)`
    // before createImageBitmap/OffscreenCanvas (absent in JSDOM) — return 404.
    if (/\/api\/exposures\/\d+\/image/.test(u)) {
      return new Response(null, { status: 404 });
    }
    if (/\/api\/experiments\/\d+\/samples/.test(u)) return json(EXP_SAMPLES);
    if (/\/api\/experiments\/\d+$/.test(u)) return json(EXPERIMENT);
    if (u.includes("/trace")) return json(TRACE);
    if (u.includes("/peaks")) return json([]);
    if (u.includes("/indices")) return json([]);
    if (u.includes("/groups")) return json([]);
    // The exposures sub-resource (/api/samples/:id/exposures) must be matched
    // BEFORE the bare /api/samples corpus route, which would otherwise swallow
    // it and hand back the corpus list.
    if (/\/exposures(\?|$)/.test(u)) return json(EXPOSURES);
    if (u.includes("/api/samples")) return json(CORPUS);
    return json([]);
  }));
});

function renderAt(path: string) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}>
        <AppRoutes />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("focus workspace header", () => {
  it("shows the serif sample name + provenance, never the picker", async () => {
    renderAt("/sample/1");

    // Serif title is the sample's display name (L-7).
    const title = await screen.findByTestId("focus-plot-title");
    expect(title).toHaveTextContent("Lipid 1-2 + LL37 1:0.5");

    // Terracotta "Integration" kicker (L-8).
    expect(screen.getByTestId("focus-plot-kicker")).toHaveTextContent("Integration");

    // Mono subline: code · beamtime · representative exposure (L-6).
    const sub = (await screen.findByTestId("focus-plot-sub")).textContent ?? "";
    expect(sub).toContain("smp_09");
    expect(sub).toContain("SSRL Apr 2026");
    expect(sub).toContain("representative exposure smp_09_e03");

    // The leftover experiment-picker must be gone.
    expect(screen.queryByTestId("plot-title")).toBeNull();
    expect(screen.queryByText(/pick an experiment/i)).toBeNull();
  });
});
