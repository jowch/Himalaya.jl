import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ModalShell } from "../../src/print/ui/ModalShell";

describe("ModalShell", () => {
  it("renders nothing when open=false", () => {
    render(
      <ModalShell open={false} onClose={() => {}} aria-label="X" testId="m">
        <button>hi</button>
      </ModalShell>,
    );
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(screen.queryByTestId("m")).toBeNull();
  });

  it("renders a role=dialog frame with aria-modal and the passed aria-label when open", () => {
    render(
      <ModalShell open onClose={() => {}} aria-label="Build thing" testId="m">
        <button>hi</button>
      </ModalShell>,
    );
    const dialog = screen.getByRole("dialog");
    expect(dialog).toHaveAttribute("aria-modal", "true");
    expect(dialog).toHaveAttribute("aria-label", "Build thing");
    expect(dialog).toHaveAttribute("data-testid", "m");
  });

  it("forwards aria-labelledby / aria-describedby to the frame", () => {
    render(
      <ModalShell open onClose={() => {}} aria-labelledby="t" aria-describedby="s" testId="m">
        <button>hi</button>
      </ModalShell>,
    );
    const dialog = screen.getByRole("dialog");
    expect(dialog).toHaveAttribute("aria-labelledby", "t");
    expect(dialog).toHaveAttribute("aria-describedby", "s");
  });

  it("role=dialog is on the frame, NOT on the scrim (regression guard for the scrim-role bug)", () => {
    render(
      <ModalShell open onClose={() => {}} aria-label="X" testId="m">
        <button>hi</button>
      </ModalShell>,
    );
    const scrim = screen.getByTestId("m-scrim");
    expect(scrim).not.toHaveAttribute("role", "dialog");
    expect(scrim).toHaveAttribute("role", "presentation");
  });

  it("exposes data-variant / data-size / data-align on the frame", () => {
    render(
      <ModalShell open onClose={() => {}} size="lg" align="top" aria-label="X" testId="m">
        <button>hi</button>
      </ModalShell>,
    );
    const dialog = screen.getByRole("dialog");
    expect(dialog).toHaveAttribute("data-variant", "dialog");
    expect(dialog).toHaveAttribute("data-size", "lg");
    expect(dialog).toHaveAttribute("data-align", "top");
  });

  it("Esc fires onClose by default; closeOnEsc=false suppresses it", () => {
    const onClose = vi.fn();
    const { rerender } = render(
      <ModalShell open onClose={onClose} aria-label="X" testId="m"><button>hi</button></ModalShell>,
    );
    fireEvent.keyDown(document, { key: "Escape" });
    expect(onClose).toHaveBeenCalledTimes(1);

    onClose.mockClear();
    rerender(
      <ModalShell open onClose={onClose} closeOnEsc={false} aria-label="X" testId="m"><button>hi</button></ModalShell>,
    );
    fireEvent.keyDown(document, { key: "Escape" });
    expect(onClose).not.toHaveBeenCalled();
  });

  it("clicking the scrim fires onClose; clicking inside the frame does not", () => {
    const onClose = vi.fn();
    render(
      <ModalShell open onClose={onClose} aria-label="X" testId="m">
        <button data-testid="inner">hi</button>
      </ModalShell>,
    );
    fireEvent.click(screen.getByTestId("inner"));
    expect(onClose).not.toHaveBeenCalled();
    fireEvent.click(screen.getByTestId("m-scrim"));
    expect(onClose).toHaveBeenCalledTimes(1);
  });

  it("closeOnOutsideClick=false suppresses scrim-click dismiss", () => {
    const onClose = vi.fn();
    render(
      <ModalShell open onClose={onClose} closeOnOutsideClick={false} aria-label="X" testId="m">
        <button>hi</button>
      </ModalShell>,
    );
    fireEvent.click(screen.getByTestId("m-scrim"));
    expect(onClose).not.toHaveBeenCalled();
  });

  it("does NOT steal initial focus to the frame (children own focus)", () => {
    const trigger = document.createElement("button");
    document.body.appendChild(trigger);
    trigger.focus();
    render(
      <ModalShell open onClose={() => {}} aria-label="X" testId="m"><button>hi</button></ModalShell>,
    );
    // ModalShell imposes no autofocus → focus stays where the caller left it.
    expect(document.activeElement).toBe(trigger);
    document.body.removeChild(trigger);
  });

  it("appends placement className to the frame", () => {
    render(
      <ModalShell open onClose={() => {}} aria-label="X" testId="m" className="max-h-[80vh]">
        <button>hi</button>
      </ModalShell>,
    );
    expect(screen.getByRole("dialog").className).toContain("max-h-[80vh]");
  });

  it("drawer variant: data-variant=drawer on the frame, scrim still present + dismissable", () => {
    const onClose = vi.fn();
    render(
      <ModalShell open onClose={onClose} variant="drawer" aria-label="Notes" testId="d">
        <textarea />
      </ModalShell>,
    );
    expect(screen.getByRole("dialog")).toHaveAttribute("data-variant", "drawer");
    fireEvent.click(screen.getByTestId("d-scrim"));
    expect(onClose).toHaveBeenCalledTimes(1);
  });

  it("traps focus inside a textarea-only frame (the Notes-drawer case)", () => {
    render(
      <ModalShell open onClose={() => {}} variant="drawer" aria-label="Notes" testId="d">
        <textarea data-testid="ta" />
      </ModalShell>,
    );
    const ta = screen.getByTestId("ta");
    ta.focus();
    const evt = fireEvent.keyDown(ta, { key: "Tab", bubbles: true });
    expect(evt).toBe(false); // trap acted (textarea is the only focusable, wraps to itself)
    expect(document.activeElement).toBe(ta);
  });
});
