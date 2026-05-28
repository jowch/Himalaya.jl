import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ScopingValueCell } from "../src/components/ScopingValueCell";

describe("ScopingValueCell", () => {
  it("a confident value renders as ink text, not an input", () => {
    render(<ScopingValueCell sampleId={10} sampleName="A" value="1:1" flagged={false}
      onChange={() => {}} onToggleFlag={() => {}} />);
    expect(screen.getByTestId("scoping-value-10")).toHaveTextContent("1:1");
    expect(screen.queryByRole("textbox")).toBeNull();
  });
  it("a flagged value shows the amber 'check the read' affordance", () => {
    render(<ScopingValueCell sampleId={11} sampleName="B" value="1 : 0" flagged
      onChange={() => {}} onToggleFlag={() => {}} />);
    const cell = screen.getByTestId("scoping-value-11");
    expect(cell).toHaveAttribute("data-flagged", "true");
    expect(cell).toHaveTextContent(/check the read/i);
  });
  it("clicking a confident value re-opens it as an input", () => {
    render(<ScopingValueCell sampleId={10} sampleName="A" value="1:1" flagged={false}
      onChange={() => {}} onToggleFlag={() => {}} />);
    fireEvent.click(screen.getByTestId("scoping-value-10"));
    expect(screen.getByTestId("scoping-value-input-10")).toHaveValue("1:1");
  });
  it("committing an edit calls onChange with the new value", () => {
    const onChange = vi.fn();
    render(<ScopingValueCell sampleId={10} sampleName="A" value="1:1" flagged={false}
      onChange={onChange} onToggleFlag={() => {}} />);
    fireEvent.click(screen.getByTestId("scoping-value-10"));
    const input = screen.getByTestId("scoping-value-input-10");
    fireEvent.change(input, { target: { value: "2:1" } });
    fireEvent.blur(input);
    expect(onChange).toHaveBeenCalledWith("2:1");
  });
  it("clicking a flagged value accepts the read (onToggleFlag)", () => {
    const onToggleFlag = vi.fn();
    render(<ScopingValueCell sampleId={11} sampleName="B" value="1 : 0" flagged
      onChange={() => {}} onToggleFlag={onToggleFlag} />);
    fireEvent.click(screen.getByTestId("scoping-value-11"));
    expect(onToggleFlag).toHaveBeenCalledTimes(1);
  });
});
