// test/DirectoryPickerField.test.tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { DirectoryPickerField } from "../src/print/components/DirectoryPickerField";

describe("DirectoryPickerField (Phase E1)", () => {
  it("shows suggestions and Tab completes the top one", () => {
    const onChange = vi.fn();
    render(
      <DirectoryPickerField
        value="/Volumes/data/ssrl/2026_04/1"
        onChange={onChange}
        suggestions={["/Volumes/data/ssrl/2026_04/1p7m", "/Volumes/data/ssrl/2026_04/10x"]}
        validation={null}
      />,
    );
    // Input puts data-testid on its WRAPPER; onKeyDown (via ...rest) lands on
    // the inner <input>. Fire on the inner input so the handler runs.
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.keyDown(input, { key: "Tab" });
    expect(onChange).toHaveBeenCalledWith("/Volumes/data/ssrl/2026_04/1p7m");
  });

  it("arrow-down moves the active suggestion, Enter completes it", () => {
    const onChange = vi.fn();
    render(
      <DirectoryPickerField
        value="/Volumes/data/ssrl/2026_04/1"
        onChange={onChange}
        suggestions={["/a/one", "/a/two"]}
        validation={null}
      />,
    );
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.keyDown(input, { key: "ArrowDown" }); // moves to index 1
    fireEvent.keyDown(input, { key: "Enter" });
    expect(onChange).toHaveBeenCalledWith("/a/two");
  });

  it("renders a positive validation line", () => {
    render(
      <DirectoryPickerField
        value="/d" onChange={() => {}} suggestions={[]}
        validation={{ ok: true, matched: 682, scanned: 700, message: null }}
      />,
    );
    expect(screen.getByTestId("dirpicker-validation")).toHaveTextContent("682");
  });

  it("renders a failure validation line with the message", () => {
    render(
      <DirectoryPickerField
        value="/d" onChange={() => {}} suggestions={[]}
        validation={{ ok: false, matched: 0, scanned: 0, message: "No exposures found" }}
      />,
    );
    const line = screen.getByTestId("dirpicker-validation");
    expect(line).toHaveAttribute("data-ok", "false");
    expect(line).toHaveTextContent("No exposures found");
  });
});
