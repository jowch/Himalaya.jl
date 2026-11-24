import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";
import { BackgroundModule, type BackgroundView } from "../../src/ui/BackgroundModule";
import { SCALE_NOTES } from "../../src/ui/CompareBar";
import type { CompareScale } from "../../src/workspace/compare";

const set = (over: Partial<Extract<BackgroundView, { kind: "set" }>> = {}): BackgroundView => ({
  kind: "set", name: "lipidA 01A (B0100)", source: "default", trace: "ready", scale: { scale: 1, source: "normalised" }, canUseRunBuffer: false, ...over,
});

function show(view: BackgroundView, over: { editable?: boolean; picker?: JSX.Element } = {}) {
  const on = { onChange: vi.fn(), onNone: vi.fn(), onDefault: vi.fn() };
  render(<BackgroundModule view={view} editable={over.editable ?? true} picker={over.picker} {...on} />);
  return on;
}

describe("BackgroundModule", () => {
  it("names the background and says it is the run's buffer", () => {
    show(set());
    expect(screen.getByRole("heading", { name: "Background" })).toBeInTheDocument();
    expect(screen.getByText("lipidA 01A (B0100)")).toBeInTheDocument();
    expect(screen.getByText("The run's first buffer")).toBeInTheDocument();
    expect(screen.queryByRole("button", { name: "Use the run's buffer" })).toBeNull();
  });

  it("says when the background is the user's choice and offers the run's buffer back", () => {
    const on = show(set({ source: "user", canUseRunBuffer: true }));
    expect(screen.getByText("Your choice")).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: "Use the run's buffer" }));
    expect(on.onDefault).toHaveBeenCalledTimes(1);
  });

  it("changes, or turns the background off, from its buttons", () => {
    const on = show(set());
    fireEvent.click(screen.getByRole("button", { name: /Change/ }));
    fireEvent.click(screen.getByRole("button", { name: "None" }));
    expect(on.onChange).toHaveBeenCalledTimes(1);
    expect(on.onNone).toHaveBeenCalledTimes(1);
  });

  it("shows how the background is scaled, in the chip's words", () => {
    show(set({ scale: { scale: 1.234, source: "matched" } }));
    expect(screen.getByText("No transmitted-beam reading, so matched at high q")).toBeInTheDocument();
    expect(screen.getByText("×1.23")).toBeInTheDocument();
  });

  it("shows a collapsed transmitted beam", () => {
    show(set({ scale: { scale: 2, source: "i2_collapsed" } }));
    expect(screen.getByText("Transmitted beam dropped on this image, so matched at high q")).toBeInTheDocument();
  });

  it("says when the background has no trace", () => {
    show(set({ trace: "missing", scale: undefined }));
    expect(screen.getByText("No trace for this background")).toBeInTheDocument();
  });

  it("says while the background's trace loads", () => {
    show(set({ trace: "loading", scale: undefined }));
    expect(screen.getByText("Loading the trace")).toBeInTheDocument();
  });

  it("says when the setup of the background is hidden", () => {
    show(set({ trace: "hidden", scale: undefined }));
    expect(screen.getByText("The setup of this background is hidden")).toBeInTheDocument();
  });

  it("reads no background when turned off, and does not offer None again", () => {
    const on = show({ kind: "off", reason: "chosen", canUseRunBuffer: true });
    expect(screen.getByText("No background. Nothing is subtracted.")).toBeInTheDocument();
    expect(screen.queryByRole("button", { name: "None" })).toBeNull();
    fireEvent.click(screen.getByRole("button", { name: "Use the run's buffer" }));
    expect(on.onDefault).toHaveBeenCalled();
  });

  it("says when the run has no buffer to use, and does not offer None, which would stop a later buffer being the default", () => {
    show({ kind: "off", reason: "no-buffer", canUseRunBuffer: false });
    expect(screen.getByText("No buffer in this run. Choose one to subtract.")).toBeInTheDocument();
    expect(screen.queryByRole("button", { name: "None" })).toBeNull();
    expect(screen.getByRole("button", { name: /Change/ })).toBeInTheDocument();
  });

  it.each(Object.keys(SCALE_NOTES) as CompareScale["source"][])("shows the whole %s note on the page, not only in a tooltip", (source) => {
    show(set({ scale: { scale: 1.5, source } }));
    expect(screen.getByText(SCALE_NOTES[source])).toBeVisible();
    expect(screen.getByText("×1.50")).toBeVisible();
    expect(screen.queryByRole("tooltip")).toBeNull();
  });

  it("says a buffer has nothing subtracted and offers no controls", () => {
    show({ kind: "buffer" });
    expect(screen.getByText("This is a buffer, so nothing is subtracted from it.")).toBeInTheDocument();
    expect(screen.queryByRole("button")).toBeNull();
  });

  it("keeps its heading and controls in place while the data loads, and cannot be used yet", () => {
    show({ kind: "loading" });
    expect(screen.getByRole("heading", { name: "Background" })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /Change/ })).toBeDisabled();
    expect(screen.getByRole("button", { name: "None" })).toBeDisabled();
  });

  it("disables the controls for someone who is not signed in", () => {
    show(set(), { editable: false });
    expect(screen.getByRole("button", { name: /Change/ })).toBeDisabled();
    expect(screen.getByRole("button", { name: "None" })).toBeDisabled();
  });

  it("holds the picker the page passes in", () => {
    show(set(), { picker: <div role="dialog" aria-label="Subtract a sample" /> });
    expect(screen.getByRole("dialog", { name: "Subtract a sample" })).toBeInTheDocument();
  });

  it("renders one name row and one controls row in every state, so the module's height does not depend on the sample", () => {
    const states: BackgroundView[] = [
      { kind: "loading" },
      { kind: "buffer" },
      { kind: "off", reason: "chosen", canUseRunBuffer: true },
      { kind: "off", reason: "no-buffer", canUseRunBuffer: false },
      set(),
      set({ source: "user", canUseRunBuffer: true, trace: "missing", scale: undefined }),
    ];
    for (const view of states) {
      const { unmount } = render(<BackgroundModule view={view} editable picker={undefined} onChange={vi.fn()} onNone={vi.fn()} onDefault={vi.fn()} />);
      expect(screen.getAllByTestId("background-controls")).toHaveLength(1);
      expect(screen.getAllByTestId("background-name")).toHaveLength(1);
      unmount();
    }
  });
});
