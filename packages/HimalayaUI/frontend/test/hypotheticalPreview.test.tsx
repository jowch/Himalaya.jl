import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, fireEvent } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { PhasePanel } from "../src/components/PhasePanel";
import { CombPanel } from "../src/components/CombPanel";
import { useAppState } from "../src/state";

const INDICES = [
  { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.0318, score: 0.9,
    r_squared: 0.99, lattice_d: 197, ngc: -1.5, status: "candidate",
    kind: "auto", inputs_hash: null,
    peaks: [{ peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.045 }],
    predicted_q: [0.045, 0.055] },
  { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.02, score: 0.4,
    r_squared: 0.7, lattice_d: 252, ngc: null, status: "candidate",
    kind: "auto", inputs_hash: null,
    peaks: [{ peak_id: 2, ratio_position: 1, residual: 0, q_observed: 0.103 }],
    predicted_q: [0.103, 0.146] },
];
const PEAKS = [
  { id: 1, exposure_id: 42, q: 0.045, intensity: 1, prominence: 1, sharpness: 1, source: "auto", excluded: false },
  { id: 2, exposure_id: 42, q: 0.103, intensity: 1, prominence: 1, sharpness: 1, source: "auto", excluded: false },
];
// Pn3m (10) is in the cart; Im3m (11) is a candidate to preview.
const ASSIGNMENT = { exposure_id: 42, state: "indexed", members: [10] };

function mountFetch() {
  vi.stubGlobal("fetch", vi.fn(async (input: RequestInfo | URL, init?: RequestInit) => {
    const u = typeof input === "string" ? input : (input as Request).url;
    const json = (b: unknown) => new Response(JSON.stringify(b),
      { status: 200, headers: { "content-type": "application/json" } });
    if (u.endsWith("/assignment")) return json(ASSIGNMENT);
    if (u.endsWith("/indices")) return json(INDICES);
    if (u.endsWith("/peaks")) return json(PEAKS);
    if (u.match(/\/experiments\/\d+$/)) return json({ id: 0, name: "E", q_units: "A-1" });
    return json([]);
  }) as typeof fetch);
}

describe("hypothetical-candidate preview (Plan D-7)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    localStorage.clear();
    useAppState.setState({
      activeExposureId: 42, activeExperimentId: 0,
      previewIndexId: undefined, hoveredIndexId: undefined,
    });
    mountFetch();
  });

  it("hovering a candidate sets previewIndexId; leaving clears it", async () => {
    renderWithProviders(<PhasePanel exposureId={42} />);
    const cand = await screen.findByRole("checkbox", { name: /Im3m/i });
    fireEvent.mouseEnter(cand);
    expect(useAppState.getState().previewIndexId).toBe(11);
    fireEvent.mouseLeave(cand);
    expect(useAppState.getState().previewIndexId).toBeUndefined();
  });

  it("blur clears the preview (no stale ghost masking the cart)", async () => {
    renderWithProviders(<PhasePanel exposureId={42} />);
    const cand = await screen.findByRole("checkbox", { name: /Im3m/i });
    fireEvent.focus(cand);
    expect(useAppState.getState().previewIndexId).toBe(11);
    fireEvent.blur(cand);
    expect(useAppState.getState().previewIndexId).toBeUndefined();
  });

  it("preview never enqueues a mutation (no POST/DELETE while previewing)", async () => {
    renderWithProviders(<PhasePanel exposureId={42} />);
    const cand = await screen.findByRole("checkbox", { name: /Im3m/i });
    const spy = global.fetch as unknown as { mock: { calls: unknown[][] } };
    const before = spy.mock.calls.length;
    fireEvent.mouseEnter(cand);
    fireEvent.mouseLeave(cand);
    const writes = spy.mock.calls.slice(before).filter((c) => {
      const init = c[1] as RequestInit | undefined;
      return init?.method === "POST" || init?.method === "DELETE";
    });
    expect(writes).toHaveLength(0);
  });

  it("CombPanel renders a dashed ghost comb row for the previewed candidate", async () => {
    useAppState.setState({ previewIndexId: 11 });
    renderWithProviders(<CombPanel />);
    // Im3m (11) is NOT in the cart → it shows as a ghost preview row.
    expect(await screen.findByTestId("comb-tooth-11-1")).toBeInTheDocument();
  });
});
