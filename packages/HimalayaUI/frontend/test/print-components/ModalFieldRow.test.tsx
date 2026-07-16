import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { ModalFieldRow } from "../../src/print/components/ModalFieldRow";

describe("<ModalFieldRow>", () => {
  it("renders with data-testid=modal-field-row", () => {
    render(
      <ModalFieldRow label="Symmetry">
        <input placeholder="control" />
      </ModalFieldRow>
    );
    expect(screen.getByTestId("modal-field-row")).toBeInTheDocument();
  });

  it("renders the label text", () => {
    render(
      <ModalFieldRow label="Symmetry">
        <input placeholder="control" />
      </ModalFieldRow>
    );
    expect(screen.getByTestId("modal-field-label")).toHaveTextContent("Symmetry");
  });

  it("label element has style.width of 100px", () => {
    render(
      <ModalFieldRow label="Symmetry">
        <input placeholder="control" />
      </ModalFieldRow>
    );
    const label = screen.getByTestId("modal-field-label");
    expect(label).toHaveStyle({ width: "100px" });
  });

  it("renders labelSuffix text when given", () => {
    render(
      <ModalFieldRow label="Lattice" labelSuffix="a">
        <input placeholder="control" />
      </ModalFieldRow>
    );
    expect(screen.getByTestId("modal-field-label")).toHaveTextContent("a");
  });

  it("does not render a labelSuffix span when omitted", () => {
    render(
      <ModalFieldRow label="Symmetry">
        <input placeholder="control" />
      </ModalFieldRow>
    );
    // labelSuffix absent — only the label text present
    expect(screen.getByTestId("modal-field-label").querySelector("span")).toBeNull();
  });

  it("renders the control-slot child", () => {
    render(
      <ModalFieldRow label="Symmetry">
        <input placeholder="my-control" />
      </ModalFieldRow>
    );
    expect(screen.getByPlaceholderText("my-control")).toBeInTheDocument();
  });
});
