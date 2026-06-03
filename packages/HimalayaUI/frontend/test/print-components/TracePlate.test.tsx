import { render, screen, fireEvent } from "@testing-library/react";
import { TracePlate } from "../../src/print/components/TracePlate";
import type { TraceModel } from "../../src/print/plot/TracePlot";

const model: TraceModel = {
  trace: { q: [0.02, 0.05, 0.1, 0.2], I: [10, 40, 20, 5], sigma: [1, 1, 1, 1] },
  peaks: [{ id: 0, q: 0.05, intensity: 40, source: "auto" }],
  phase: "Pn3m",
};

const base = {
  title: "Lipid 1-2 + LL37",
  trace: model,
  scale: "log" as const,
  onScaleChange: () => {},
};

describe("TracePlate", () => {
  it("renders the title and the trace plot region", () => {
    render(<TracePlate {...base} kicker="Integration" subtitle="smp_09" />);
    expect(screen.getByTestId("trace-plate")).toBeInTheDocument();
    expect(screen.getByText("Lipid 1-2 + LL37")).toBeInTheDocument();
    expect(screen.getByText("Integration")).toBeInTheDocument();
    expect(screen.getByText("smp_09")).toBeInTheDocument();
  });
  it("calls onScaleChange when the scale toggle is used", () => {
    const onScaleChange = vi.fn();
    render(<TracePlate {...base} onScaleChange={onScaleChange} />);
    fireEvent.click(screen.getByText("linear q"));
    expect(onScaleChange).toHaveBeenCalledWith("lin");
  });
  it("shows Auto-fit and fires onAutoFit", () => {
    const onAutoFit = vi.fn();
    render(<TracePlate {...base} onAutoFit={onAutoFit} />);
    fireEvent.click(screen.getByText("Auto-fit"));
    expect(onAutoFit).toHaveBeenCalled();
  });
  it("reflects the armed '+ Peak' toggle and fires onToggleAddPeak", () => {
    const onToggleAddPeak = vi.fn();
    render(<TracePlate {...base} addPeakArmed onToggleAddPeak={onToggleAddPeak} />);
    const peak = screen.getByText("+ Peak");
    expect(peak).toHaveAttribute("aria-pressed", "true");
    fireEvent.click(peak);
    expect(onToggleAddPeak).toHaveBeenCalled();
  });
  it("forwards a placement-only className", () => {
    render(<TracePlate {...base} className="mt-6" />);
    expect(screen.getByTestId("trace-plate").className).toContain("mt-6");
  });
});
