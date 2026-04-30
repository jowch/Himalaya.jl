import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, waitFor, fireEvent } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { SpeculativeBuilder } from "../src/components/SpeculativeBuilder";

interface MockHandlers {
  peaks?: unknown[];
  snap?: unknown[];
  /** Captures POST body for assertions. */
  onPost?: (body: unknown) => void;
}

function mockApi(opts: MockHandlers): void {
  const peaks = opts.peaks ?? [
    { id: 1, exposure_id: 42, q: 0.10, intensity: 100, prominence: 50, sharpness: 0.8, source: "auto", excluded: false },
    { id: 2, exposure_id: 42, q: 0.20, intensity: 80,  prominence: 40, sharpness: 0.5, source: "auto", excluded: false },
    { id: 3, exposure_id: 42, q: 0.30, intensity: 60,  prominence: 30, sharpness: 0.4, source: "auto", excluded: false },
  ];
  const snap = opts.snap ?? [
    { ratio_position: 1, predicted_q: 0.10, suggested_peak_id: null, suggested_q: null, suggested_residual: null, is_anchor: true },
    { ratio_position: 2, predicted_q: 0.20, suggested_peak_id: 2,    suggested_q: 0.20, suggested_residual: 0.0001, is_anchor: false },
    { ratio_position: 3, predicted_q: 0.30, suggested_peak_id: 3,    suggested_q: 0.30, suggested_residual: 0.0001, is_anchor: false },
  ];
  vi.spyOn(global, "fetch").mockImplementation(async (input, init) => {
    const u = typeof input === "string" ? input : (input as Request).url;
    if (u.endsWith("/peaks")) {
      return new Response(JSON.stringify(peaks),
        { status: 200, headers: { "Content-Type": "application/json" } });
    }
    if (u.includes("/speculative-snap")) {
      return new Response(JSON.stringify(snap),
        { status: 200, headers: { "Content-Type": "application/json" } });
    }
    if (u.endsWith("/speculative") && (init?.method === "POST")) {
      const body = init.body ? JSON.parse(String(init.body)) : null;
      opts.onPost?.(body);
      return new Response(JSON.stringify({
        id: 99, exposure_id: 42, phase: body?.phase ?? "Lamellar", basis: 0.1,
        score: 0.5, r_squared: 1.0, lattice_d: 60, status: "candidate",
        kind: "speculative", peaks: [], predicted_q: [],
      }), { status: 200, headers: { "Content-Type": "application/json" } });
    }
    return new Response("not found", { status: 404 });
  });
}

describe("<SpeculativeBuilder>", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("renders phase and anchor selectors", async () => {
    mockApi({});
    renderWithProviders(<SpeculativeBuilder exposureId={42} onClose={() => {}} />);
    expect(await screen.findByTestId("spec-phase-select")).toBeInTheDocument();
    expect(screen.getByTestId("spec-anchor-select")).toBeInTheDocument();
  });

  it("auto-picks the strongest peak as anchor and fetches snap suggestions", async () => {
    mockApi({});
    renderWithProviders(<SpeculativeBuilder exposureId={42} onClose={() => {}} />);
    // Anchor select should default to peak 1 (highest sharpness in fixtures).
    await waitFor(() => {
      const sel = screen.getByTestId("spec-anchor-select") as HTMLSelectElement;
      expect(sel.value).toBe("1");
    });
    // Snap row for ratio 2 appears.
    expect(await screen.findByTestId("spec-snap-row-2")).toBeInTheDocument();
  });

  it("save POSTs anchor + checked additional peaks and calls onClose", async () => {
    let captured: unknown = null;
    mockApi({ onPost: (body) => { captured = body; } });
    const onClose = vi.fn();
    renderWithProviders(<SpeculativeBuilder exposureId={42} onClose={onClose} />);

    // Wait for the snap rows to render
    await screen.findByTestId("spec-snap-row-2");
    await screen.findByTestId("spec-snap-row-3");

    // Toggle ratio 2 on (ratio 3 stays off to verify selective inclusion).
    const r2 = screen.getByTestId("spec-snap-row-2").querySelector("input[type=checkbox]") as HTMLInputElement;
    fireEvent.click(r2);
    expect(r2.checked).toBe(true);

    fireEvent.click(screen.getByTestId("spec-save-button"));

    await waitFor(() => { expect(onClose).toHaveBeenCalled(); });
    expect(captured).toMatchObject({
      phase: "Lamellar",
      anchor_peak_id: 1,
      anchor_ratio: 1,
      additional: [{ ratio_position: 2, peak_id: 2 }],
      active: false,
    });
  });

  it("Cancel button closes without POST", async () => {
    let captured: unknown = null;
    mockApi({ onPost: (body) => { captured = body; } });
    const onClose = vi.fn();
    renderWithProviders(<SpeculativeBuilder exposureId={42} onClose={onClose} />);
    await screen.findByTestId("spec-phase-select");
    fireEvent.click(screen.getByRole("button", { name: /cancel/i }));
    expect(onClose).toHaveBeenCalled();
    expect(captured).toBeNull();
  });
});
