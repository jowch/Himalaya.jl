import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import type { Exposure } from "../src/api";
import { LoupeFrame } from "../src/components/LoupeFrame";

// DetectorImage touches fetch / createImageBitmap / OffscreenCanvas (absent in
// JSDOM); mock it — LoupeFrame's own behaviour does not depend on its render.
vi.mock("../src/components/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));

function exposure(over: Partial<Exposure> = {}): Exposure {
  return {
    id: 100, sample_id: 7, filename: "JC042-001.dat", kind: "file",
    selected: false, status: null, image_path: null, image_version: "",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
    ...over,
  };
}

describe("LoupeFrame", () => {
  it("renders one strip thumbnail per exposure", () => {
    const exposures = [
      exposure({ id: 1 }), exposure({ id: 2 }), exposure({ id: 3 }),
    ];
    render(
      <LoupeFrame
        exposure={exposures[0]}
        exposures={exposures}
        onSelectExposure={() => {}}
      />,
    );
    expect(screen.getByTestId("thumb-cell-1")).toBeInTheDocument();
    expect(screen.getByTestId("thumb-cell-2")).toBeInTheDocument();
    expect(screen.getByTestId("thumb-cell-3")).toBeInTheDocument();
  });

  it("calls onSelectExposure with the clicked exposure id", () => {
    const onSelect = vi.fn();
    const exposures = [exposure({ id: 1 }), exposure({ id: 2 })];
    render(
      <LoupeFrame
        exposure={exposures[0]}
        exposures={exposures}
        onSelectExposure={onSelect}
      />,
    );
    fireEvent.click(screen.getByTestId("thumb-cell-2"));
    expect(onSelect).toHaveBeenCalledWith(2);
  });

  it("shows the Dropped badge when the active exposure is rejected", () => {
    const e = exposure({ id: 1, status: "rejected" });
    const { rerender } = render(
      <LoupeFrame exposure={e} exposures={[e]} onSelectExposure={() => {}} />,
    );
    expect(screen.getByTestId("loupe-dropped-badge")).toBeInTheDocument();

    const kept = exposure({ id: 1, status: "accepted" });
    rerender(
      <LoupeFrame exposure={kept} exposures={[kept]} onSelectExposure={() => {}} />,
    );
    expect(screen.queryByTestId("loupe-dropped-badge")).not.toBeInTheDocument();
  });

  it("draws the grease-pencil reject ✕ on the big frame when rejected (M-10)", () => {
    const e = exposure({ id: 1, status: "rejected" });
    const { rerender } = render(
      <LoupeFrame exposure={e} exposures={[e]} onSelectExposure={() => {}} />,
    );
    expect(screen.getByTestId("loupe-reject-xmark")).toBeInTheDocument();

    const kept = exposure({ id: 1, status: "accepted" });
    rerender(
      <LoupeFrame exposure={kept} exposures={[kept]} onSelectExposure={() => {}} />,
    );
    expect(screen.queryByTestId("loupe-reject-xmark")).not.toBeInTheDocument();
  });
});
