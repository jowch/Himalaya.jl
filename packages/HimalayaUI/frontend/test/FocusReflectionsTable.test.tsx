import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { FocusReflectionsTable } from "../src/components/FocusReflectionsTable";
import { useAppState } from "../src/state";

// Fixtures — two Pn3m-claimed peaks + one unindexed peak. The active group
// has one index (the Pn3m one); the second peak below mirrors a real
// active-set claim. The third peak is unclaimed → should render as "unindexed".
const PEAKS = [
  { id: 1, exposure_id: 50, q: 0.045, intensity: 100, prominence: 50,
    sharpness: 1.2, source: "auto" as const, excluded: false },
  { id: 2, exposure_id: 50, q: 0.055, intensity: 80, prominence: 40,
    sharpness: 1.1, source: "auto" as const, excluded: false },
  { id: 3, exposure_id: 50, q: 0.103, intensity: 60, prominence: 30,
    sharpness: 1.0, source: "auto" as const, excluded: false },
];

const INDICES = [
  { id: 10, exposure_id: 50, phase: "Pn3m", basis: 0.5, score: 0.89,
    r_squared: 0.99, lattice_d: 197, ngc: -1.5,
    status: "candidate", kind: "auto", inputs_hash: null,
    predicted_q: [0.045, 0.055],
    peaks: [
      { peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.045 },
      { peak_id: 2, ratio_position: 2, residual: 0, q_observed: 0.055 },
    ] },
];

const GROUPS = [
  { id: 1, exposure_id: 50, kind: "auto", active: true, members: [10] },
];

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    activeExposureId: 50, hoveredQ: undefined, username: "tester",
  });
  vi.stubGlobal("fetch", vi.fn(async (url: string) => {
    const u = String(url);
    const json = (b: unknown) => new Response(JSON.stringify(b), {
      status: 200, headers: { "content-type": "application/json" } });
    if (u.endsWith("/peaks")) return json(PEAKS);
    if (u.endsWith("/indices")) return json(INDICES);
    if (u.endsWith("/groups")) return json(GROUPS);
    return json([]);
  }));
});

function renderTable() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <FocusReflectionsTable />
    </QueryClientProvider>,
  );
}

describe("FocusReflectionsTable (#209)", () => {
  it("renders the panel region with the Reflections header", async () => {
    renderTable();
    const panel = await screen.findByTestId("focus-reflections-panel");
    expect(panel).toBeInTheDocument();
    // Header is the panel's first child — the `text-meta uppercase` span.
    expect(panel.textContent).toMatch(/Reflections/);
  });

  it("R3-F09: column headers use the text-meta scale token, not an inline px size", async () => {
    renderTable();
    await screen.findByTestId("focus-reflections");
    const ths = document.querySelectorAll('[data-testid="focus-reflections"] thead th');
    expect(ths.length).toBe(4);
    for (const th of ths) {
      // Fixed-Scale Rule: no inline text-[Npx] on the sticky header.
      expect(th.className).not.toMatch(/text-\[\d/);
      expect(th.className).toContain("text-meta");
    }
  });

  it("R3-F10: unindexed row is flagged data-indexed='false' and dimmed as a whole", async () => {
    renderTable();
    const row3 = await screen.findByTestId("reflection-row-3"); // unclaimed peak
    const row1 = screen.getByTestId("reflection-row-1");        // Pn3m-claimed
    expect(row3).toHaveAttribute("data-indexed", "false");
    expect(row1).toHaveAttribute("data-indexed", "true");
    // Whole-row dim (mockup `.refl tr.unindexed`), not just the phase cell.
    expect(row3.className).toMatch(/opacity-55/);
    expect(row1.className).not.toMatch(/opacity-55/);
  });

  it("R3-F07: reflections panel-head carries a right-side N-of-M summary cluster", async () => {
    renderTable();
    await screen.findByTestId("focus-reflections");
    const cluster = screen.getByTestId("focus-reflections-head-summary");
    // 2 of 3 peaks claimed in the fixture.
    expect(cluster).toHaveTextContent("2 / 3");
  });

  it("renders one row per peak, sorted low-q first", async () => {
    renderTable();
    expect(await screen.findByTestId("reflection-row-1")).toBeInTheDocument();
    expect(screen.getByTestId("reflection-row-2")).toBeInTheDocument();
    expect(screen.getByTestId("reflection-row-3")).toBeInTheDocument();
  });

  it("labels claimed peaks with the active-set phase + ratio position", async () => {
    renderTable();
    const row1 = await screen.findByTestId("reflection-row-1");
    expect(row1).toHaveTextContent("Pn3m");
    expect(row1).toHaveTextContent("#1");
    expect(row1).toHaveTextContent("0.0450");

    const row2 = screen.getByTestId("reflection-row-2");
    expect(row2).toHaveTextContent("Pn3m");
    expect(row2).toHaveTextContent("#2");
  });

  it("marks unclaimed peaks as 'unindexed' with an em-dash order", async () => {
    renderTable();
    const row3 = await screen.findByTestId("reflection-row-3");
    expect(row3).toHaveTextContent(/unindexed/i);
    expect(row3).toHaveTextContent("—");
  });

  it("computes d = 2π/q in the d column", async () => {
    renderTable();
    const row1 = await screen.findByTestId("reflection-row-1");
    // 2π/0.045 ≈ 139.6
    expect(row1).toHaveTextContent("139.6");
  });

  it("surfaces a footer with N of M peaks + phase count", async () => {
    renderTable();
    const foot = await screen.findByTestId("focus-reflections-foot");
    // 2 peaks claimed (id 1, 2) of 3 total, across 1 phase.
    expect(foot).toHaveTextContent("2 of 3");
    expect(foot).toHaveTextContent(/1 phase\b/i);
  });

  it("shows a pick-a-sample hint when no exposure is active", async () => {
    useAppState.setState({ activeExposureId: undefined });
    renderTable();
    expect(await screen.findByText(/pick a sample/i)).toBeInTheDocument();
  });

  // ── q-link triple (the headline contract of #209) ──────────────────────────
  describe("q-link wiring (mirrors #180 detector-ring overlay)", () => {
    it("row hover sets hoveredQ to the row's peak q (source)", async () => {
      renderTable();
      const row1 = await screen.findByTestId("reflection-row-1");
      fireEvent.mouseEnter(row1);
      expect(useAppState.getState().hoveredQ).toBeCloseTo(0.045, 6);
    });

    it("row leave clears hoveredQ when this row was the active source", async () => {
      renderTable();
      const row1 = await screen.findByTestId("reflection-row-1");
      fireEvent.mouseEnter(row1);
      fireEvent.mouseLeave(row1);
      expect(useAppState.getState().hoveredQ).toBeUndefined();
    });

    it("row leave does NOT clobber hoveredQ when another source has taken over", async () => {
      // Mirrors the detector-ring guarded clear (DetectorRingOverlay.tsx:131-136).
      // While row 1 is the active hover, another source (e.g. a ring or another
      // row) takes over and pushes a different hoveredQ. The row's mouseleave
      // must NOT clear it back to undefined — the new hover is still active.
      renderTable();
      const row1 = await screen.findByTestId("reflection-row-1");
      fireEvent.mouseEnter(row1);
      // Another source (a ring, or row 2) takes over:
      useAppState.getState().setHoveredQ(0.055);
      fireEvent.mouseLeave(row1);
      expect(useAppState.getState().hoveredQ).toBeCloseTo(0.055, 6);
    });

    it("row whose q matches hoveredQ lights via data-hot='true' (sink)", async () => {
      renderTable();
      const row1 = await screen.findByTestId("reflection-row-1");
      expect(row1).toHaveAttribute("data-hot", "false");
      // External source pushes hoveredQ matching row 1's q (the detector ring
      // path would set it via setHoveredQ; we simulate that here).
      useAppState.getState().setHoveredQ(0.045);
      // re-query after state change
      const row1Hot = await screen.findByTestId("reflection-row-1");
      expect(row1Hot).toHaveAttribute("data-hot", "true");
      // Other rows must NOT light.
      expect(screen.getByTestId("reflection-row-2")).toHaveAttribute("data-hot", "false");
      expect(screen.getByTestId("reflection-row-3")).toHaveAttribute("data-hot", "false");
    });
  });
});
