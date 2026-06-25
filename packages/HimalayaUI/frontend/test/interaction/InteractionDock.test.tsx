import { describe, it, expect, vi, afterEach } from "vitest";
import { render, screen, fireEvent, cleanup } from "@testing-library/react";
import { InteractionDock } from "../../src/print/interaction/InteractionDock";
import { useInteraction } from "../../src/print/interaction/registry";
import { core, page } from "../../src/print/interaction/core";
import type { ListCursor } from "../../src/print/interaction/types";

const fakeCursor = (over: Partial<ListCursor> = {}): ListCursor => ({
  cursorId: 1, selected: new Set(), setCursor: () => {}, moveBy: vi.fn(),
  activate: () => {}, toggleSelect: () => {}, rowProps: () => ({} as never),
  stepperProps: () => ({ label: "Sample", axis: "vertical", testIdBase: "sample", count: "1 / 3", onPrev: () => {}, onNext: () => {}, prevDisabled: true, nextDisabled: false }),
  ...over,
});

afterEach(() => { useInteraction.getState().clearPage(); cleanup(); });

describe("InteractionDock", () => {
  it("renders nothing when the store is empty", () => {
    const { container } = render(<InteractionDock />);
    expect(container.querySelector('[data-testid="dock"]')).toBeNull();
  });

  it("renders the stepper from cursor.stepperProps()", () => {
    useInteraction.getState().setPage(fakeCursor(), []);
    render(<InteractionDock />);
    expect(screen.getByTestId("dock-sample-count")).toHaveTextContent("1 / 3");
    expect(screen.getByTestId("dock-prev-sample")).toBeDisabled();
  });

  it("stepper next button calls cursor.moveBy(1)", () => {
    const moveBy = vi.fn();
    useInteraction.getState().setPage(fakeCursor({ moveBy }), []);
    render(<InteractionDock />);
    fireEvent.click(screen.getByTestId("dock-next-sample"));
    expect(moveBy).toHaveBeenCalledWith(1);
  });

  it("renders dock:true action buttons and greys disabled ones", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [
      page("cull", { label: "Cull", keys: ["x"], group: "Act", enabled: () => false, dock: true, run }),
    ]);
    render(<InteractionDock />);
    const btn = screen.getByRole("button", { name: /cull/i });
    expect(btn).toBeDisabled();
  });

  it("renders the primary action prominently and runs it on click", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [core("openFocus", { run, dock: "primary" })]);
    render(<InteractionDock />);
    fireEvent.click(screen.getByTestId("dock-primary"));
    expect(run).toHaveBeenCalledTimes(1);
  });

  it("renders the up-link when core('back') is declared", () => {
    useInteraction.getState().setPage(null, [core("back", { run: () => {}, label: "Corpus", dock: true })]);
    render(<InteractionDock />);
    expect(screen.getByTestId("dock-up-link")).toHaveTextContent("Corpus");
  });
});
