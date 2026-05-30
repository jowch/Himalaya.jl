import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, fireEvent, within } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { CombPanel } from "../src/components/CombPanel";
import { useAppState } from "../src/state";

// Pn3m claims peaks 10 & 11 at √2,√3; predicts a √4 order that is absent.
const INDICES = [
  {
    id: 1, exposure_id: 50, phase: "Pn3m", basis: 0.0318, score: 0.9,
    r_squared: 0.99, lattice_d: 197, ngc: -1.5, status: "candidate",
    kind: "auto", inputs_hash: null,
    peaks: [
      { peak_id: 10, ratio_position: 1, residual: 0.0001, q_observed: 0.045 },
      { peak_id: 11, ratio_position: 2, residual: 0.0002, q_observed: 0.055 },
    ],
    // √2, √3, √4 → the third order (0.0636) is predicted but unobserved
    predicted_q: [0.045, 0.055, 0.0636],
  },
];
const PEAKS = [
  { id: 10, exposure_id: 50, q: 0.045, intensity: 1, prominence: 1, sharpness: 1, source: "auto", excluded: false },
  { id: 11, exposure_id: 50, q: 0.055, intensity: 1, prominence: 1, sharpness: 1, source: "auto", excluded: false },
  // a leftover observed peak no assigned phase claims
  { id: 12, exposure_id: 50, q: 0.103, intensity: 1, prominence: 1, sharpness: 1, source: "auto", excluded: false },
];
const ASSIGNMENT = { exposure_id: 50, state: "indexed", members: [1] };

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({ activeExposureId: 50, hoveredQ: undefined });
  vi.stubGlobal("fetch", vi.fn(async (url: string) => {
    const u = String(url);
    const json = (b: unknown) => new Response(JSON.stringify(b),
      { status: 200, headers: { "content-type": "application/json" } });
    if (u.endsWith("/assignment")) return json(ASSIGNMENT);
    if (u.endsWith("/indices")) return json(INDICES);
    if (u.endsWith("/peaks")) return json(PEAKS);
    if (u.endsWith("/trace")) return json({ q: [0.04, 0.06, 0.11], intensity: [1, 2, 1], q_units: "A-1" });
    return json([]);
  }));
});

describe("<CombPanel> (Plan D-5)", () => {
  it("renders a tooth per predicted order for an assigned phase", async () => {
    renderWithProviders(<CombPanel />);
    await screen.findByTestId("comb-tooth-1-1");
    // 2 observed orders (√2, √3) render as solid teeth; √4 is absent.
    const teeth = screen.getAllByTestId(/^comb-tooth-/);
    expect(teeth.length).toBeGreaterThanOrEqual(2);
  });

  it("renders a hollow caret for a predicted-but-absent order", async () => {
    renderWithProviders(<CombPanel />);
    // the √4 order (0.0636) is predicted but unobserved → absent marker
    expect(await screen.findByTestId("comb-absent-1-3")).toBeInTheDocument();
  });

  it("renders the leftover row with the unclaimed observed peak", async () => {
    renderWithProviders(<CombPanel />);
    expect(await screen.findByTestId("comb-leftover-12")).toBeInTheDocument();
  });

  it("the residual toggle switches to the residual view", async () => {
    renderWithProviders(<CombPanel />);
    await screen.findByTestId("comb-svg");
    expect(screen.queryByTestId("comb-residual-zero-line")).not.toBeInTheDocument();
    fireEvent.click(screen.getByRole("radio", { name: /residual/i }));
    expect(await screen.findByTestId("comb-residual-zero-line")).toBeInTheDocument();
  });

  it("hovering hoveredQ lights the matching tooth (hot = grow + ring, not recolour)", async () => {
    renderWithProviders(<CombPanel />);
    await screen.findByTestId("comb-svg");
    // set hoveredQ to the first Pn3m peak q
    useAppState.setState({ hoveredQ: 0.045 });
    const tooth = await screen.findByTestId("comb-tooth-1-1");
    expect(tooth).toHaveAttribute("data-hot", "true");
  });

  it("hovering a tooth sets hoveredQ (q-link source)", async () => {
    renderWithProviders(<CombPanel />);
    const tooth = await screen.findByTestId("comb-tooth-1-1");
    fireEvent.mouseEnter(tooth);
    expect(useAppState.getState().hoveredQ).toBeCloseTo(0.045, 4);
  });
});
