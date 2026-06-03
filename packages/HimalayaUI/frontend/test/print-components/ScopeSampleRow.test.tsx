import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { ScopeSampleRow } from "../../src/print/components/ScopeSampleRow";

const TRACE = { q: [0.02, 0.04, 0.08, 0.16], I: [10, 80, 30, 5] };

describe("<ScopeSampleRow>", () => {
  it("renders name, id, sparkline, grip, and the value", () => {
    render(<ScopeSampleRow name="Lipid 1-2 + LL37 1:0.25" sampleId="smp_07" trace={TRACE} phase="Pn3m" value="1 : 0.25" />);
    const row = screen.getByTestId("scope-sample-row");
    expect(row).toHaveTextContent("Lipid 1-2 + LL37 1:0.25");
    expect(row).toHaveTextContent("smp_07");
    expect(screen.getByTestId("sparkline")).toBeInTheDocument();
    expect(screen.getByTestId("grip-handle")).toBeInTheDocument();
    expect(screen.getByTestId("flag-button")).toHaveTextContent("1 : 0.25");
  });
  it("marks the flagged state on the row and the value", () => {
    render(<ScopeSampleRow name="Lipid 1-2, no LL37" sampleId="smp_04" trace={TRACE} phase="Pn3m" value="1 : 0" flagged />);
    expect(screen.getByTestId("scope-sample-row").getAttribute("data-flagged")).toBe("true");
    expect(screen.getByTestId("flag-button").getAttribute("data-state")).toBe("flagged");
  });
  it("forwards the flag toggle", () => {
    const onToggleFlag = vi.fn();
    render(<ScopeSampleRow name="x" sampleId="smp_04" trace={TRACE} value="1 : 0" flagged onToggleFlag={onToggleFlag} />);
    fireEvent.click(screen.getByTestId("flag-button"));
    expect(onToggleFlag).toHaveBeenCalledOnce();
  });
});
