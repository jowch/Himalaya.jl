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

  it("a disabled Export states its reason and ties both actions to it (BU-EXPORT-EMPTYALLOWED)", () => {
    setup({ disabled: true, disabledReason: "Traces are still loading." });
    const reason = screen.getByTestId("export-disabled-reason");
    expect(reason).toHaveTextContent("Traces are still loading.");
    expect(reason).toHaveAttribute("role", "note");
    // Both actions point at the reason for assistive tech.
    expect(screen.getByTestId("export-copy")).toHaveAttribute("aria-describedby", reason.id);
    expect(screen.getByTestId("export-menu-trigger")).toHaveAttribute("aria-describedby", reason.id);
  });

  it("does NOT render a reason while enabled, even if disabledReason is passed", () => {
    setup({ disabled: false, disabledReason: "should not show" });
    expect(screen.queryByTestId("export-disabled-reason")).toBeNull();
    expect(screen.getByTestId("export-copy")).not.toHaveAttribute("aria-describedby");
  });

  it("a disabled Export with no reason renders no reason node (back-compat)", () => {
    setup({ disabled: true });
    expect(screen.queryByTestId("export-disabled-reason")).toBeNull();
  });

  it("Escape closes the open menu", () => {
    setup();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    const menu = screen.getByRole("menu", { name: /download formats/i });
    fireEvent.keyDown(menu, { key: "Escape" });
    expect(screen.queryByRole("menu", { name: /download formats/i })).toBeNull();
  });

  // APG menu-button keyboard pattern.
  it("ArrowDown on the closed trigger opens and focuses the first enabled item", () => {
    setup();
    const trigger = screen.getByTestId("export-menu-trigger");
    trigger.focus();
    fireEvent.keyDown(trigger, { key: "ArrowDown" });
    expect(document.activeElement).toBe(
      screen.getByRole("menuitem", { name: /download as png/i }),
    );
  });

  it("ArrowUp on the closed trigger opens and focuses the LAST enabled item", () => {
    setup();
    const trigger = screen.getByTestId("export-menu-trigger");
    trigger.focus();
    fireEvent.keyDown(trigger, { key: "ArrowUp" });
    expect(document.activeElement).toBe(
      screen.getByRole("menuitem", { name: /download as svg/i }),
    );
  });

  it("Escape inside the menu closes and RETURNS focus to the trigger", () => {
    setup();
    const trigger = screen.getByTestId("export-menu-trigger");
    fireEvent.click(trigger);
    fireEvent.keyDown(screen.getByRole("menu", { name: /download formats/i }), {
      key: "Escape",
    });
    expect(screen.queryByRole("menu", { name: /download formats/i })).toBeNull();
    expect(document.activeElement).toBe(trigger);
  });

  it("select-refocus survives the pending window (parked focus across the disabled trigger)", () => {
    // Real sequence: selecting a download flips `pending` synchronously, which
    // disables the trigger (`disabled || pending`). Real browsers blur a
    // focused element the moment it disables — JSDOM doesn't, so model the
    // blur — and the parked-focus effect must restore the trigger when
    // pending transitions back to false (RepresentativeBox precedent).
    const props = {
      onCopy: vi.fn(),
      onDownloadPng: vi.fn(),
      onDownloadSvg: vi.fn(),
      ariaContext: "the trace",
    };
    const { rerender } = render(<ExportButton {...props} />);
    const trigger = screen.getByTestId("export-menu-trigger");
    fireEvent.click(trigger);
    fireEvent.click(screen.getByRole("menuitem", { name: /download as png/i }));
    rerender(<ExportButton {...props} pending />);
    expect(trigger).toBeDisabled();
    // Model the blur-on-disable JSDOM omits (its blur() doesn't even unfocus —
    // move focus off the trigger explicitly).
    screen.getByTestId("export-copy").focus();
    expect(document.activeElement).not.toBe(trigger);
    rerender(<ExportButton {...props} />); // render lands → pending ends
    expect(document.activeElement).toBe(trigger);
  });

  it("a pending cycle WITHOUT a local menu activation never yanks focus (foreign flip)", () => {
    const props = {
      onCopy: vi.fn(),
      onDownloadPng: vi.fn(),
      onDownloadSvg: vi.fn(),
      ariaContext: "the trace",
    };
    const { rerender } = render(<ExportButton {...props} />);
    rerender(<ExportButton {...props} pending />);
    rerender(<ExportButton {...props} />);
    expect(document.activeElement).not.toBe(screen.getByTestId("export-menu-trigger"));
  });

  it("an outside pointerdown closes the open menu", () => {
    setup();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    expect(screen.getByRole("menu", { name: /download formats/i })).toBeInTheDocument();
    fireEvent.pointerDown(document.body);
    expect(screen.queryByRole("menu", { name: /download formats/i })).toBeNull();
  });

  it("focus leaving the widget (Tab away) closes the open menu (UI-MENUTAB)", () => {
    setup();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    expect(screen.getByRole("menu", { name: /download formats/i })).toBeInTheDocument();
    // APG: Tab closes the menu. jsdom does not move focus on Tab, so model the
    // focusout whose relatedTarget is an element OUTSIDE the export widget
    // (real Tab is covered by e2e). focusOut bubbles to the wrapper's onBlur.
    fireEvent.focusOut(screen.getByRole("menuitem", { name: /download as png/i }), {
      relatedTarget: document.body,
    });
    expect(screen.queryByRole("menu", { name: /download formats/i })).toBeNull();
  });

  it("focus moving to the trigger (inside the widget) keeps the menu open (no toggle-reopen)", () => {
    // The regression close-on-blur would introduce if scoped to the menu alone:
    // focus → trigger would close, then the trigger's toggle onClick reopens.
    // Scoping the close to the whole widget (wrapRef) keeps focus-to-trigger in.
    setup();
    fireEvent.click(screen.getByRole("button", { name: /download formats/i }));
    fireEvent.focusOut(screen.getByRole("menuitem", { name: /download as png/i }), {
      relatedTarget: screen.getByTestId("export-menu-trigger"),
    });
    expect(screen.getByRole("menu", { name: /download formats/i })).toBeInTheDocument();
  });
});
