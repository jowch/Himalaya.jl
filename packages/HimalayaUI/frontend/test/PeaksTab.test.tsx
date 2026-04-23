import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, waitFor } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { PeaksTab } from "../src/components/PeaksTab";
import { useAppState } from "../src/state";

beforeEach(() => {
  vi.restoreAllMocks();
  useAppState.setState({ activeExposureId: undefined });
});

function mockPeaks(items: unknown[]): void {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify(items), {
      status: 200, headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("<PeaksTab>", () => {
  it("shows a hint when no exposure is active", () => {
    renderWithProviders(<PeaksTab />);
    expect(screen.getByText(/no exposure selected/i)).toBeInTheDocument();
  });

  it("shows an empty-state message when the exposure has no peaks", async () => {
    useAppState.setState({ activeExposureId: 42 });
    mockPeaks([]);
    renderWithProviders(<PeaksTab />);
    await waitFor(() => expect(screen.getByText(/no peaks/i)).toBeInTheDocument());
  });

  it("renders one table row per peak with q, prominence, sharpness, source", async () => {
    useAppState.setState({ activeExposureId: 42 });
    mockPeaks([
      { id: 1, exposure_id: 42, q: 0.123, intensity: 10,
        prominence: 0.8, sharpness: 2.5, source: "auto" },
      { id: 2, exposure_id: 42, q: 0.456, intensity: 20,
        prominence: null, sharpness: null, source: "manual" },
    ]);
    renderWithProviders(<PeaksTab />);
    await waitFor(() => expect(screen.getByText("0.1230")).toBeInTheDocument());
    expect(screen.getByText("0.4560")).toBeInTheDocument();
    expect(screen.getByText("0.80")).toBeInTheDocument();
    expect(screen.getByText("2.50")).toBeInTheDocument();
    expect(screen.getByText("auto")).toBeInTheDocument();
    expect(screen.getByText("manual")).toBeInTheDocument();
  });
});
