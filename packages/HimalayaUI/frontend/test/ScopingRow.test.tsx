import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ScopingRow } from "../src/components/ScopingRow";

const baseRow = { sampleId: 10, sampleName: "Lipid 1-2 + LL37 1:1", value: "1:1", flagged: false, include: true };

describe("ScopingRow", () => {
  it("renders the sample name and a sparkline frame, with no fake drag grip", () => {
    render(<ScopingRow row={baseRow} trace={undefined} phase={null}
      onChangeValue={() => {}} onToggleFlag={() => {}} />);
    const row = screen.getByTestId("scoping-row-10");
    expect(row).toHaveTextContent("Lipid 1-2 + LL37 1:1");
    expect(screen.getByTestId("scoping-sparkline-empty")).toBeInTheDocument();
    // M2: the decorative drag grip is gone (rows are value-ordered, not dragged).
    expect(screen.queryByTestId("scoping-grip-10")).toBeNull();
  });
  it("marks the row data-flagged when flagged", () => {
    render(<ScopingRow row={{ ...baseRow, flagged: true }} trace={undefined} phase={null}
      onChangeValue={() => {}} onToggleFlag={() => {}} />);
    expect(screen.getByTestId("scoping-row-10")).toHaveAttribute("data-flagged", "true");
  });
  it("threads value edits up via onChangeValue", () => {
    const onChangeValue = vi.fn();
    render(<ScopingRow row={baseRow} trace={undefined} phase={null}
      onChangeValue={onChangeValue} onToggleFlag={() => {}} />);
    fireEvent.click(screen.getByTestId("scoping-value-10"));
    const input = screen.getByTestId("scoping-value-input-10");
    fireEvent.change(input, { target: { value: "3:1" } });
    fireEvent.blur(input);
    expect(onChangeValue).toHaveBeenCalledWith("3:1");
  });
});
