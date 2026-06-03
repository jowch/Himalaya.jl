import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ExportButton } from "../../src/print/components/ExportButton";

function setup(overrides: Partial<React.ComponentProps<typeof ExportButton>> = {}) {
  const props = {
    onCopy: vi.fn(),
    onDownloadPng: vi.fn(),
    onDownloadSvg: vi.fn(),
    ariaContext: "the trace",
    ...overrides,
  };
  render(<ExportButton {...props} />);
  return props;
}

describe("ExportButton", () => {
  it("renders the Copy action with an aria-label built from ariaContext", () => {
    setup();
    expect(screen.getByRole("button", { name: /copy the trace to clipboard/i }))
      .toBeInTheDocument();
  });

  it("does not show the download menu until the chevron is clicked", () => {
    setup();
    expect(screen.queryByRole("menu", { name: /download formats/i })).toBeNull();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    expect(screen.getByRole("menu", { name: /download formats/i })).toBeInTheDocument();
    expect(screen.getByRole("menuitem", { name: /download as png/i })).toBeInTheDocument();
    expect(screen.getByRole("menuitem", { name: /download as svg/i })).toBeInTheDocument();
  });

  it("primary click calls onCopy", () => {
    const props = setup();
    fireEvent.click(screen.getByTestId("export-copy"));
    expect(props.onCopy).toHaveBeenCalledTimes(1);
  });

  it("selecting PNG calls onDownloadPng and closes the menu", () => {
    const props = setup();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    fireEvent.click(screen.getByRole("menuitem", { name: /download as png/i }));
    expect(props.onDownloadPng).toHaveBeenCalledTimes(1);
    expect(screen.queryByRole("menu", { name: /download formats/i })).toBeNull();
  });

  it("selecting SVG calls onDownloadSvg", () => {
    const props = setup();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    fireEvent.click(screen.getByRole("menuitem", { name: /download as svg/i }));
    expect(props.onDownloadSvg).toHaveBeenCalledTimes(1);
  });

  it("disables Copy when copyDisabled", () => {
    setup({ copyDisabled: true });
    expect(screen.getByTestId("export-copy")).toBeDisabled();
  });

  it("disables the PNG menuitem when pngDisabled (SVG stays enabled)", () => {
    setup({ pngDisabled: true });
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    expect(screen.getByRole("menuitem", { name: /download as png/i })).toBeDisabled();
    expect(screen.getByRole("menuitem", { name: /download as svg/i })).not.toBeDisabled();
  });

  it("the page-level disabled prop ORs into Copy and the trigger", () => {
    setup({ disabled: true });
    expect(screen.getByTestId("export-copy")).toBeDisabled();
    expect(screen.getByRole("button", { name: /download formats/i })).toBeDisabled();
  });

  it("Escape closes the open menu", () => {
    setup();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    const menu = screen.getByRole("menu", { name: /download formats/i });
    fireEvent.keyDown(menu, { key: "Escape" });
    expect(screen.queryByRole("menu", { name: /download formats/i })).toBeNull();
  });

  it("an outside pointerdown closes the open menu", () => {
    setup();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    expect(screen.getByRole("menu", { name: /download formats/i })).toBeInTheDocument();
    fireEvent.pointerDown(document.body);
    expect(screen.queryByRole("menu", { name: /download formats/i })).toBeNull();
  });
});
