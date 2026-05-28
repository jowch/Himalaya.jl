import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { FocusDetectorPanel } from "../src/components/FocusDetectorPanel";
import { useAppState } from "../src/state";

const EXPOSURES = [
  { id: 5, sample_id: 1, filename: "JC001-005.dat", kind: "file" as const,
    selected: true, status: "accepted" as const,
    image_path: "/img/5.png", image_version: "v1", tags: [], sources: [],
    trace_hash: null, analysis_inputs_hash: null },
  { id: 6, sample_id: 1, filename: "JC001-006.dat", kind: "file" as const,
    selected: false, status: "accepted" as const,
    image_path: "/img/6.png", image_version: "v1", tags: [], sources: [],
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
  vi.stubGlobal("fetch", vi.fn(async (url: string, init?: RequestInit) => {
    const u = String(url);
    const json = (b: unknown) => new Response(JSON.stringify(b), {
      status: 200, headers: { "content-type": "application/json" } });
    // The detector PNG route: return non-ok so DetectorImage.renderImage bails
    // at `if (!res.ok) return;` BEFORE calling createImageBitmap (absent in
    // JSDOM). We assert the ring overlay, not the rendered canvas pixels.
    if (u.includes("/image")) return new Response(null, { status: 404 });
    // PATCH /api/exposures/:id/select — set-representative echoes the row.
    if (/\/exposures\/\d+\/select/.test(u) && init?.method === "PATCH") {
      const id = Number(u.match(/\/exposures\/(\d+)\/select/)![1]);
      return json({ id, selected: true });
    }
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

  it("R3-N3 (#209): the panel <section> wears the Plate Lift `.card` class", async () => {
    // Round-2 finding R3-N3: the detector panel shipped with the bare
    // `rounded border bg-plate` triple — no Plate Lift, so it read as
    // flat against the warm paper. DESIGN.md §Elevation: plate-like
    // cards belong to the lifted family.
    renderPanel();
    const panel = await screen.findByTestId("focus-detector-panel");
    expect(panel.className).toContain("card");
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

  // ── F-11: representative-exposure switcher ───────────────────────────────
  it("renders a thumbnail per exposure in the switcher strip", async () => {
    renderPanel();
    await screen.findByTestId("focus-detector-panel");
    expect(await screen.findByTestId("exposure-switcher")).toBeInTheDocument();
    expect(await screen.findByTestId("exposure-thumb-5")).toBeInTheDocument();
    expect(await screen.findByTestId("exposure-thumb-6")).toBeInTheDocument();
  });

  it("marks the selected exposure as the representative", async () => {
    renderPanel();
    // exposure 5 has selected:true in the fixture
    const thumb5 = await screen.findByTestId("exposure-thumb-5");
    expect(thumb5).toHaveAttribute("data-rep", "true");
    expect((await screen.findByTestId("exposure-thumb-6"))).not.toHaveAttribute("data-rep");
  });

  it("clicking a thumbnail switches the active (viewed) exposure", async () => {
    renderPanel();
    const thumb6 = await screen.findByTestId("exposure-thumb-6");
    fireEvent.click(thumb6);
    await waitFor(() =>
      expect(useAppState.getState().activeExposureId).toBe(6),
    );
  });

  it("set-representative persists via the select route", async () => {
    renderPanel();
    // switch to exposure 6 first, then set it representative
    fireEvent.click(await screen.findByTestId("exposure-thumb-6"));
    await waitFor(() => expect(useAppState.getState().activeExposureId).toBe(6));
    fireEvent.click(await screen.findByTestId("exposure-set-rep"));
    await waitFor(() => {
      const calls = (global.fetch as unknown as { mock: { calls: unknown[][] } }).mock.calls;
      const hit = calls.some(([u, init]) =>
        /\/exposures\/6\/select/.test(String(u)) &&
        (init as RequestInit | undefined)?.method === "PATCH");
      expect(hit).toBe(true);
    });
  });

  it("does not render the switcher when there is only one exposure", async () => {
    vi.stubGlobal("fetch", vi.fn(async (url: string) => {
      const u = String(url);
      const json = (b: unknown) => new Response(JSON.stringify(b), {
        status: 200, headers: { "content-type": "application/json" } });
      if (u.includes("/image")) return new Response(null, { status: 404 });
      if (u.includes("/peaks")) return json(PEAKS);
      if (u.includes("/exposures")) return json([EXPOSURES[0]]);
      return json([]);
    }));
    renderPanel();
    await screen.findByTestId("focus-detector-panel");
    await waitFor(() =>
      expect(screen.queryByTestId("exposure-switcher")).not.toBeInTheDocument(),
    );
  });
});
