import { describe, it, expect, vi } from "vitest";
import { render } from "@testing-library/react";
import { TraceViewer } from "../src/components/TraceViewer";

vi.mock("@observablehq/plot", () => ({
  plot: vi.fn(() => {
    const el = document.createElement("div");
    el.setAttribute("data-testid", "plot-svg");
    return el;
  }),
  areaY: vi.fn(() => ({ _kind: "areaY" })),
  line:  vi.fn(() => ({ _kind: "line" })),
  dot:   vi.fn(() => ({ _kind: "dot" })),
}));

describe("<TraceViewer>", () => {
  it("renders a container with testid 'trace-viewer'", () => {
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 30], sigma: [1, 1, 1] };
    const { container } = render(
      <TraceViewer
        trace={trace}
        peaks={[]}
        onAddPeak={() => {}}
        onRemovePeak={() => {}}
      />,
    );
    expect(container.querySelector('[data-testid="trace-viewer"]')).not.toBeNull();
  });

  it("invokes Plot.plot with marks including line and areaY for trace", async () => {
    const Plot = await import("@observablehq/plot");
    const trace = { q: [0.1, 0.2], I: [10, 20], sigma: [1, 1] };
    render(
      <TraceViewer trace={trace} peaks={[]} onAddPeak={() => {}} onRemovePeak={() => {}} />,
    );
    expect(Plot.plot).toHaveBeenCalled();
    expect(Plot.areaY).toHaveBeenCalled();
    expect(Plot.line).toHaveBeenCalled();
  });

  it("passes auto peaks and manual peaks to separate dot marks", async () => {
    const Plot = await import("@observablehq/plot");
    (Plot.dot as unknown as { mockClear: () => void }).mockClear();
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 30], sigma: [1, 1, 1] };
    const peaks = [
      { id: 1, exposure_id: 1, q: 0.1, intensity: null, prominence: null,
        sharpness: null, source: "auto"   as const },
      { id: 2, exposure_id: 1, q: 0.2, intensity: null, prominence: null,
        sharpness: null, source: "manual" as const },
    ];
    render(
      <TraceViewer trace={trace} peaks={peaks}
        onAddPeak={() => {}} onRemovePeak={() => {}} />,
    );
    expect((Plot.dot as unknown as { mock: { calls: unknown[] } }).mock.calls.length).toBe(2);
  });
});
