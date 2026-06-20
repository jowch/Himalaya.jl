import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { ConfigurationBody } from "../src/print/components/ConfigurationBody";

const updateMutate = vi.fn();

vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  return {
    ...actual,
    // E1's Experiment carries PER-FIELD *_source columns (NOT a combined
    // beam_center_source): energy_kev_source, flight_path_m_source,
    // beam_center_x_source, beam_center_y_source, pixel_size_um_source,
    // q_units_source -- plus ingest_status/scan_signature/last_scanned_at.
    useExperiment: () => ({
      data: {
        id: 7, name: "SSRL · 1p7m",
        description: "April 2026 beamtime at SSRL 1p7m",
        path: "/d", data_dir: "/d", analysis_dir: "/d/analysis",
        image_pattern: "{name}.tiff", metadata_pattern: "{name}.prp", integration_pattern: "{name}.dat",
        manifest_path: null, created_at: "2026-04-12", q_units: "A^-1",
        energy_kev: 9.0, energy_kev_source: "prp",
        flight_path_m: 1.81, flight_path_m_source: "setup",
        beam_center_x: 421.4, beam_center_x_source: "setup",
        beam_center_y: 836.9, beam_center_y_source: "setup",
        pixel_size_um: 172, pixel_size_um_source: "prp", q_units_source: "prp",
        last_scanned_at: "2026-04-12T10:00:00", scan_signature: "sig", ingest_status: "idle",
      },
      isLoading: false,
    }),
    useLoads: () => ({ data: [], isLoading: false }),
    useUpdateExperiment: () => ({ mutate: updateMutate, isPending: false }),
  };
});

const wrap = (n: ReactNode) => render(<QueryClientProvider client={new QueryClient()}>{n}</QueryClientProvider>);

describe("ConfigurationBody", () => {
  beforeEach(() => { updateMutate.mockClear(); });

  it("renders the description, Geometry, Acquisition, and Sources cards", () => {
    wrap(<ConfigurationBody experimentId={7} />);
    // description (DECISION: E1-owned column; E2 renders it here)
    expect(screen.getByText(/April 2026 beamtime/)).toBeInTheDocument();
    expect(screen.getByText("Geometry")).toBeInTheDocument();
    expect(screen.getByText("Acquisition")).toBeInTheDocument();
    expect(screen.getByText("Sources")).toBeInTheDocument();
  });
  it("derives a geometry row per typed field with its *_source provenance", () => {
    wrap(<ConfigurationBody experimentId={7} />);
    expect(screen.getByText("Beam energy")).toBeInTheDocument();
    // flight_path_m, beam_center_x, beam_center_y all have source="setup" -- use getAllByText
    expect(screen.getAllByText("setup files").length).toBeGreaterThan(0);
  });
  it("renders the 3 pattern rows as editable (E1 columns: image/metadata/integration_pattern)", () => {
    wrap(<ConfigurationBody experimentId={7} />);
    expect(screen.getByText("{name}.tiff")).toBeInTheDocument();  // image_pattern
  });

  // --- inline override: real-change triggers PATCH ---
  it("Override activates inline input seeded with raw numeric value; typing a new value + Enter fires PATCH", () => {
    wrap(<ConfigurationBody experimentId={7} />);

    // Click Override for "Beam energy" (first override button)
    fireEvent.click(screen.getAllByRole("button", { name: /override beam energy/i })[0]!);

    // An inline input should appear seeded with the raw numeric string (no units)
    const inp = screen.getByRole("textbox", { name: /override beam energy/i });
    expect((inp as HTMLInputElement).value).toBe("9");

    // Change the value and commit via Enter
    fireEvent.change(inp, { target: { value: "10" } });
    fireEvent.keyDown(inp, { key: "Enter" });

    // updateMutate must have been called with the NEW parsed value
    expect(updateMutate).toHaveBeenCalledWith(
      expect.objectContaining({ energy_kev: 10 })
    );
  });

  // --- inline override: no-change does NOT fire PATCH ---
  it("Override + commit with the SAME value does NOT call updateMutate (no mis-stamp)", () => {
    wrap(<ConfigurationBody experimentId={7} />);

    fireEvent.click(screen.getAllByRole("button", { name: /override beam energy/i })[0]!);

    const inp = screen.getByRole("textbox", { name: /override beam energy/i });
    // The draft is seeded to "9" (raw value). Do NOT change it -- commit with Enter.
    expect((inp as HTMLInputElement).value).toBe("9");
    fireEvent.keyDown(inp, { key: "Enter" });

    expect(updateMutate).not.toHaveBeenCalled();
  });

  // --- Escape cancels without PATCH ---
  it("Override + Escape cancels without calling updateMutate", () => {
    wrap(<ConfigurationBody experimentId={7} />);

    fireEvent.click(screen.getAllByRole("button", { name: /override beam energy/i })[0]!);
    const inp = screen.getByRole("textbox", { name: /override beam energy/i });
    fireEvent.keyDown(inp, { key: "Escape" });

    expect(updateMutate).not.toHaveBeenCalled();
    // The value display should be back (no input visible)
    expect(screen.queryByRole("textbox", { name: /override beam energy/i })).toBeNull();
  });
});
