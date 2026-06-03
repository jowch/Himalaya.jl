import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Field } from "../../src/print/ui/Field";

describe("<Field>", () => {
  it("renders the value and is a button", () => {
    render(<Field value="LL37 : lipid ratio" />);
    const f = screen.getByTestId("field");
    expect(f.tagName).toBe("BUTTON");
    expect(f).toHaveTextContent("LL37 : lipid ratio");
  });
  it("fires onClick", () => {
    const onClick = vi.fn();
    render(<Field value="x" onClick={onClick} />);
    fireEvent.click(screen.getByTestId("field"));
    expect(onClick).toHaveBeenCalledOnce();
  });
  it("shows the placeholder when value is empty", () => {
    render(<Field value="" placeholder="Choose a variable" />);
    expect(screen.getByTestId("field")).toHaveTextContent("Choose a variable");
  });
  it("opens a menu and fires onSelect in dropdown mode", () => {
    const onSelect = vi.fn();
    render(<Field value="Time" options={[{value:"Time",label:"Time"},{value:"Dose",label:"Dose"}]} onSelect={onSelect} />);
    expect(screen.queryByTestId("menu")).toBeNull();          // closed initially
    fireEvent.click(screen.getByTestId("field"));
    expect(screen.getByTestId("menu")).toBeInTheDocument();    // opens
    fireEvent.click(screen.getByRole("menuitem", { name: "Dose" }));
    expect(onSelect).toHaveBeenCalledWith("Dose");
    expect(screen.queryByTestId("menu")).toBeNull();          // closes on select
  });
});
