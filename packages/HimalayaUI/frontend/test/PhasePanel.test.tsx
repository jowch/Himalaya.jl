import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, waitFor, fireEvent } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { PhasePanel } from "../src/components/PhasePanel";
import { useAppState } from "../src/state";

beforeEach(() => { vi.restoreAllMocks(); });

function mockAll(indices: unknown[], groups: unknown[]): void {
  vi.spyOn(global, "fetch").mockImplementation(async (input) => {
    const u = typeof input === "string" ? input : (input as Request).url;
    if (u.endsWith("/indices")) {
      return new Response(JSON.stringify(indices),
        { status: 200, headers: { "Content-Type": "application/json" } });
    }
    if (u.endsWith("/groups")) {
      return new Response(JSON.stringify(groups),
        { status: 200, headers: { "Content-Type": "application/json" } });
    }
    if (u.match(/\/api\/groups\/\d+\/members\/\d+$/)) {
      return new Response(JSON.stringify({
        id: 2, exposure_id: 42, kind: "custom", active: true, members: [],
      }), { status: 200, headers: { "Content-Type": "application/json" } });
    }
    return new Response("not found", { status: 404 });
  });
}

describe("<PhasePanel> — base", () => {
  it("renders a hint when no exposure is active", () => {
    renderWithProviders(<PhasePanel exposureId={undefined} />);
    expect(screen.getByText(/no exposure selected/i)).toBeInTheDocument();
  });
});

describe("<PhasePanel> — phase-call output block (R4 L-9/L-10)", () => {
  it("renders a Phase call block per active-set phase with serif name, score and series ratio", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 0.89,
          r_squared: 0.998, lattice_d: 197, ngc: -1.5, status: "candidate", kind: "auto",
          predicted_q: [0.045, 0.055, 0.064],
          peaks: [
            { peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.045 },
            { peak_id: 2, ratio_position: 2, residual: 0, q_observed: 0.055 },
            { peak_id: 3, ratio_position: 3, residual: 0, q_observed: 0.064 },
          ] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const block = await screen.findByTestId("phase-call-block-10");
    expect(block).toHaveTextContent("Pn3m");
    expect(block).toHaveTextContent("0.89");
    expect(block).toHaveTextContent("197");
    expect(block).toHaveTextContent("√2 : √3 : 2");
  });

  it("shows a Coexistence header when two phases are in the call", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 0.89,
          r_squared: 0.99, lattice_d: 197, ngc: -1.5, status: "candidate", kind: "auto",
          predicted_q: [0.045], peaks: [{ peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.045 }] },
        { id: 11, exposure_id: 42, phase: "Lamellar", basis: 0.3, score: 0.96,
          r_squared: 0.99, lattice_d: 61, ngc: null, status: "candidate", kind: "auto",
          predicted_q: [0.103], peaks: [{ peak_id: 4, ratio_position: 1, residual: 0, q_observed: 0.103 }] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10, 11] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    await screen.findByTestId("phase-call-block-10");
    expect(screen.getByTestId("coexistence-tag")).toHaveTextContent(/Coexistence.*2 phases/i);
  });

  it("shows an empty phase-call message when no index is in the call", async () => {
    mockAll(
      [
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
          r_squared: 0.71, lattice_d: 9.1, ngc: null, status: "candidate", kind: "auto",
          predicted_q: [0.4], peaks: [{ peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.4 }] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    expect(await screen.findByTestId("phase-call-empty")).toBeInTheDocument();
  });
});

describe("<PhasePanel> — candidate multi-select (R4 L-10)", () => {
  it("renders candidates as checkboxes reflecting active-set membership", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 0.89,
          r_squared: 0.99, lattice_d: 197, ngc: -1.5, status: "candidate", kind: "auto",
          predicted_q: [0.045], peaks: [{ peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.045 }] },
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
          r_squared: 0.71, lattice_d: 9.1, ngc: null, status: "candidate", kind: "auto",
          predicted_q: [0.4], peaks: [{ peak_id: 2, ratio_position: 1, residual: 0, q_observed: 0.4 }] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const inCall = await screen.findByRole("checkbox", { name: /Pn3m/i });
    expect(inCall).toBeChecked();
    const candidate = screen.getByRole("checkbox", { name: /Im3m/i });
    expect(candidate).not.toBeChecked();
  });

  it("surfaces the series ratio on a candidate row", async () => {
    mockAll(
      [
        { id: 11, exposure_id: 42, phase: "Lamellar", basis: 0.3, score: 0.6,
          r_squared: 0.99, lattice_d: 61, ngc: null, status: "candidate", kind: "auto",
          predicted_q: [0.103, 0.206],
          peaks: [
            { peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.103 },
            { peak_id: 2, ratio_position: 2, residual: 0, q_observed: 0.206 },
          ] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const row = await screen.findByTestId("candidate-row-11");
    expect(row).toHaveTextContent("1 : 2");
    // The row is keyboard-operable (role=checkbox, tabIndex); it must show a
    // focus-visible ring so keyboard users can see where they are.
    expect(row.className).toContain("focus-visible:outline");
  });

  it("toggling an unchecked candidate posts add-to-group", async () => {
    vi.spyOn(global, "fetch").mockImplementation(async (input) => {
      const u = typeof input === "string" ? input : (input as Request).url;
      if (u.endsWith("/indices")) return new Response(JSON.stringify([
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 0.89, r_squared: 0.99,
          lattice_d: 197, ngc: -1.5, status: "candidate", kind: "auto", predicted_q: [0.045],
          peaks: [{ peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.045 }] },
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6, r_squared: 0.71,
          lattice_d: 9.1, ngc: null, status: "candidate", kind: "auto", predicted_q: [0.4],
          peaks: [{ peak_id: 2, ratio_position: 1, residual: 0, q_observed: 0.4 }] },
      ]), { status: 200, headers: { "Content-Type": "application/json" } });
      if (u.endsWith("/groups")) return new Response(JSON.stringify([
        { id: 2, exposure_id: 42, kind: "custom", active: true, members: [10] },
      ]), { status: 200, headers: { "Content-Type": "application/json" } });
      if (u.endsWith("/api/groups/2/members")) return new Response(JSON.stringify({
        id: 2, exposure_id: 42, kind: "custom", active: true, members: [10, 11],
      }), { status: 200, headers: { "Content-Type": "application/json" } });
      return new Response("not found", { status: 404 });
    });
    renderWithProviders(<PhasePanel exposureId={42} />);
    const candidate = await screen.findByRole("checkbox", { name: /Im3m/i });
    fireEvent.click(candidate);
    await waitFor(() => {
      const spy = global.fetch as unknown as { mock: { calls: unknown[][] } };
      const urls = spy.mock.calls.map((c) => typeof c[0] === "string" ? c[0] : (c[0] as Request).url);
      expect(urls).toContain("/api/groups/2/members");
    });
  });

  it("toggling a checked candidate posts remove-from-group", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 0.89, r_squared: 0.99,
          lattice_d: 197, ngc: -1.5, status: "candidate", kind: "auto", predicted_q: [0.045],
          peaks: [{ peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.045 }] },
      ],
      [{ id: 2, exposure_id: 42, kind: "custom", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const checked = await screen.findByRole("checkbox", { name: /Pn3m/i });
    fireEvent.click(checked);
    await waitFor(() => {
      const spy = global.fetch as unknown as { mock: { calls: unknown[][] } };
      const urls = spy.mock.calls.map((c) => typeof c[0] === "string" ? c[0] : (c[0] as Request).url);
      expect(urls.some((u) => u === "/api/groups/2/members/10")).toBe(true);
    });
  });

  it("hovering a candidate sets hoveredIndexId; leaving clears it", async () => {
    useAppState.setState({ hoveredIndexId: undefined });
    mockAll(
      [
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6, r_squared: 0.71,
          lattice_d: 9.1, ngc: null, status: "candidate", kind: "auto", predicted_q: [0.4],
          peaks: [{ peak_id: 2, ratio_position: 1, residual: 0, q_observed: 0.4 }] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const row = await screen.findByTestId("candidate-row-11");
    fireEvent.mouseEnter(row);
    expect(useAppState.getState().hoveredIndexId).toBe(11);
    fireEvent.mouseLeave(row);
    expect(useAppState.getState().hoveredIndexId).toBeUndefined();
  });
});

describe("<PhasePanel> — R4-N1 copy dedup (#209)", () => {
  it("renders the explanatory sentence in exactly one rail slot, not both", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 0.89,
          r_squared: 0.99, lattice_d: 197, ngc: -1.5, status: "candidate", kind: "auto",
          predicted_q: [0.045],
          peaks: [{ peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.045 }] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    await screen.findByTestId("phase-call-block-10");
    // The "Check every phase…" framing sentence appeared verbatim in BOTH
    // the card-header subtitle AND the rail-note paragraph (PhasePanel.tsx
    // :273-275 and :343-346, round-2 finding R4-N1). It must appear exactly
    // once now — the rail-note copy, which the mockup has below candidates.
    const occurrences = screen.queryAllByText(/check every phase that is present/i);
    expect(occurrences).toHaveLength(1);
  });
});

describe("<PhasePanel> — speculative", () => {
  it("renders speculative indices as candidate rows with a delete affordance + add button", async () => {
    mockAll(
      [
        { id: 20, exposure_id: 42, phase: "Lamellar", basis: 0.3, score: 0.55, r_squared: 1.0,
          lattice_d: 52, ngc: null, status: "candidate", kind: "speculative", predicted_q: [0.21, 0.42],
          peaks: [
            { peak_id: 5, ratio_position: 1, residual: 0, q_observed: 0.21 },
            { peak_id: 6, ratio_position: 2, residual: 0, q_observed: 0.42 },
          ] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    expect(await screen.findByTestId("candidate-row-20")).toBeInTheDocument();
    const del = screen.getByTestId("spec-delete-20");
    expect(del).toBeInTheDocument();
    // The delete control is an ink-stroke SVG, not the off-brand 🗑 emoji
    // (the only emoji that used to live in the rail).
    expect(del.querySelector("svg")).not.toBeNull();
    expect(del.textContent).not.toContain("🗑");
    expect(screen.getByTestId("add-speculative-button")).toBeInTheDocument();
  });

  it("R3-F02: Speculative lives behind a collapsed disclosure when none exist", async () => {
    // No speculative index — only an auto candidate. The rail should read as
    // the calm two-section output (Phase call + Candidates), with Speculative
    // tucked below the fold behind a collapsed <details>.
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 0.89,
          r_squared: 0.99, lattice_d: 197, ngc: -1.5, status: "candidate", kind: "auto",
          predicted_q: [0.045], peaks: [{ peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.045 }] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const disclosure = await screen.findByTestId("speculative-disclosure");
    // Contract: a NATIVE <details> collapsed by default (the browser hides the
    // body, incl. the CTA, until opened).
    expect(disclosure.tagName.toLowerCase()).toBe("details");
    expect(disclosure).not.toHaveAttribute("open");
    // The CTA is a DESCENDANT of the (closed) disclosure — behind-the-fold,
    // not an unconditional sibling. NOTE: JSDOM does not honor the native
    // display:none of closed <details>, so queryByTestId would still FIND the
    // CTA — assert containment, not absence.
    const cta = screen.getByTestId("add-speculative-button");
    expect(disclosure).toContainElement(cta);
  });

  it("R3-F02: disclosure is open when speculatives exist (user builds stay visible)", async () => {
    mockAll(
      [
        { id: 20, exposure_id: 42, phase: "Lamellar", basis: 0.3, score: 0.55, r_squared: 1.0,
          lattice_d: 52, ngc: null, status: "candidate", kind: "speculative", predicted_q: [0.21],
          peaks: [{ peak_id: 5, ratio_position: 1, residual: 0, q_observed: 0.21 }] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const disclosure = await screen.findByTestId("speculative-disclosure");
    expect(disclosure).toHaveAttribute("open");
  });
});
