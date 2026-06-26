import { describe, it, expect, vi, afterEach } from "vitest";
import { render, screen, fireEvent, cleanup } from "@testing-library/react";
import { InteractionDock } from "../../src/print/interaction/InteractionDock";
import { useInteraction } from "../../src/print/interaction/registry";
import { core, page } from "../../src/print/interaction/core";
import type { CursorStepperProps, ListCursor } from "../../src/print/interaction/types";

const fakeCursor = (over: Partial<ListCursor> = {}): ListCursor => {
  const moveBy = (over.moveBy ?? vi.fn()) as ListCursor["moveBy"];
  return {
    cursorId: 1,
    selected: new Set(),
    setCursor: () => {},
    activate: () => {},
    toggleSelect: () => {},
    rowProps: () => ({} as never),
    stepperProps: () => ({
      label: "Sample", axis: "vertical", testIdBase: "sample", count: "1 / 3",
      onPrev: () => moveBy(-1), onNext: () => moveBy(1),
      prevDisabled: true, nextDisabled: false,
    }),
    ...over,
    moveBy, // ensure the exposed moveBy is the SAME ref stepperProps closes over
  };
};

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

  it("renders extraSteppers BEFORE the cursor stepper when both are present", () => {
    const extraStepper: CursorStepperProps = {
      label: "Sample", axis: "vertical", testIdBase: "sample",
      count: "2 / 5", onPrev: vi.fn(), onNext: vi.fn(),
      prevDisabled: false, nextDisabled: false,
    };
    // fakeCursor uses testIdBase "sample" → reuse but supply a frame-axis cursor
    const frameCursor: ListCursor = {
      ...fakeCursor(),
      stepperProps: () => ({
        label: "Frame", axis: "horizontal", testIdBase: "frame", count: "1 / 3",
        onPrev: vi.fn(), onNext: vi.fn(), prevDisabled: true, nextDisabled: false,
      }),
    };
    useInteraction.getState().setPage(frameCursor, [], [extraStepper]);
    render(<InteractionDock />);
    // Both steppers are present
    expect(screen.getByTestId("dock-sample-count")).toHaveTextContent("2 / 5");
    expect(screen.getByTestId("dock-frame-count")).toHaveTextContent("1 / 3");
    // Extra stepper (sample) appears BEFORE cursor stepper (frame) in DOM
    const sampleCount = screen.getByTestId("dock-sample-count");
    const frameCount = screen.getByTestId("dock-frame-count");
    expect(sampleCount.compareDocumentPosition(frameCount) & Node.DOCUMENT_POSITION_FOLLOWING).toBeTruthy();
  });

  it("renders dockExtra content when provided (dock shows even without cursor/actions)", () => {
    useInteraction.getState().setPage(null, [], [], <div data-testid="dock-extra-slot">identity</div>);
    render(<InteractionDock />);
    // The dock renders because dockExtra is truthy
    expect(screen.getByTestId("dock-extra-slot")).toHaveTextContent("identity");
    expect(screen.getByTestId("dock")).toBeInTheDocument();
  });
});
