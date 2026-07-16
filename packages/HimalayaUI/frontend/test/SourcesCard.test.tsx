// test/SourcesCard.test.tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { SourcesCard, type SourceRow } from "../src/print/components/SourcesCard";

const ROWS: SourceRow[] = [
  { key: "data_dir", label: "Data directory", value: "/Volumes/data/ssrl/2026_04/1p7m", editable: false },
  { key: "analysis_dir", label: "Analysis directory", value: "/Volumes/analysis/ssrl/2026_04/1p7m", editable: false },
  { key: "image_pattern", label: "Image pattern", value: "{name}.tiff", editable: true },
  { key: "metadata_pattern", label: "Metadata pattern", value: "{name}.dat", editable: true },
  { key: "integration_pattern", label: "Integration pattern", value: "{name}.chi", editable: true },
];

describe("SourcesCard", () => {
  it("renders each source row label", () => {
    render(<SourcesCard rows={ROWS} onEdit={() => {}} onRescan={() => {}} />);
    expect(screen.getByText("Data directory")).toBeInTheDocument();
    expect(screen.getByText("Image pattern")).toBeInTheDocument();
    expect(screen.getByText("{name}.tiff")).toBeInTheDocument();
  });

  it("editable pattern row: clicking value opens an input and Enter fires onEdit", () => {
    const onEdit = vi.fn();
    render(<SourcesCard rows={ROWS} onEdit={onEdit} onRescan={() => {}} />);
    // The editable value is rendered as a button that becomes an input on click
    fireEvent.click(screen.getByText("{name}.tiff"));
    const input = screen.getByDisplayValue("{name}.tiff");
    fireEvent.change(input, { target: { value: "{name}.tif" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(onEdit).toHaveBeenCalledWith("image_pattern", "{name}.tif");
  });

  it("opening an editable row focuses AND selects the input (immediate retype)", () => {
    render(<SourcesCard rows={ROWS} onEdit={() => {}} onRescan={() => {}} />);
    fireEvent.click(screen.getByText("{name}.tiff"));
    const input = screen.getByDisplayValue("{name}.tiff") as HTMLInputElement;
    expect(document.activeElement).toBe(input);
    expect(input.selectionStart).toBe(0);
    expect(input.selectionEnd).toBe("{name}.tiff".length);
  });

  it("read-only directory rows have no edit control", () => {
    render(<SourcesCard rows={ROWS} onEdit={() => {}} onRescan={() => {}} />);
    // The directory value should be plain text, not a button or input
    const dirValue = screen.getByText("/Volumes/data/ssrl/2026_04/1p7m");
    expect(dirValue.tagName).not.toBe("BUTTON");
    expect(dirValue.tagName).not.toBe("INPUT");
    // Clicking the read-only value should NOT open an input
    fireEvent.click(dirValue);
    expect(screen.queryByRole("textbox")).toBeNull();
  });

  it("Rescan now calls onRescan", () => {
    const onRescan = vi.fn();
    render(<SourcesCard rows={ROWS} onEdit={() => {}} onRescan={onRescan} />);
    fireEvent.click(screen.getByRole("button", { name: /rescan now/i }));
    expect(onRescan).toHaveBeenCalled();
  });

  it("all 3 pattern rows are editable (have click-to-edit buttons)", () => {
    const onEdit = vi.fn();
    render(<SourcesCard rows={ROWS} onEdit={onEdit} onRescan={() => {}} />);

    // image_pattern
    fireEvent.click(screen.getByText("{name}.tiff"));
    const imgInput = screen.getByDisplayValue("{name}.tiff");
    fireEvent.change(imgInput, { target: { value: "{name}.tif" } });
    fireEvent.keyDown(imgInput, { key: "Escape" });

    // metadata_pattern
    fireEvent.click(screen.getByText("{name}.dat"));
    const metaInput = screen.getByDisplayValue("{name}.dat");
    fireEvent.change(metaInput, { target: { value: "{name}.data" } });
    fireEvent.keyDown(metaInput, { key: "Enter" });
    expect(onEdit).toHaveBeenCalledWith("metadata_pattern", "{name}.data");

    // integration_pattern
    fireEvent.click(screen.getByText("{name}.chi"));
    const intInput = screen.getByDisplayValue("{name}.chi");
    fireEvent.change(intInput, { target: { value: "{name}.xy" } });
    fireEvent.keyDown(intInput, { key: "Enter" });
    expect(onEdit).toHaveBeenCalledWith("integration_pattern", "{name}.xy");
  });
});
