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
  it("names the background and tags it as the run buffer", () => {
    show(set());
    expect(screen.getByRole("heading", { name: "Background" })).toBeInTheDocument();
    expect(screen.getByText("lipidA 01A (B0100)")).toBeInTheDocument();
    expect(screen.getByText("run buffer")).toBeInTheDocument();
    expect(screen.queryByRole("button", { name: "Run buffer" })).toBeNull();
  });

  it("tags the background as chosen and offers the run buffer back", () => {
    const on = show(set({ source: "user", canUseRunBuffer: true }));
    expect(screen.getByText("chosen")).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: "Run buffer" }));
    expect(on.onDefault).toHaveBeenCalledTimes(1);
  });

  it("changes, or turns the background off, from its buttons", () => {
    const on = show(set());
    fireEvent.click(screen.getByRole("button", { name: /Change/ }));
    fireEvent.click(screen.getByRole("button", { name: "None" }));
    expect(on.onChange).toHaveBeenCalledTimes(1);
    expect(on.onNone).toHaveBeenCalledTimes(1);
  });

  it("says nothing on the note row when both traces are already normalised", () => {
    show(set());
    expect(screen.getByTestId("background-note")).toBeEmptyDOMElement();
  });

  it("shows a matched scale as a short label with its value", () => {
    show(set({ scale: { scale: 1.234, source: "matched" } }));
    expect(screen.getByTestId("background-note")).toHaveTextContent("matched at high q ×1.23");
  });

  it("shows a collapsed transmitted beam", () => {
    show(set({ scale: { scale: 2, source: "i2_collapsed" } }));
    expect(screen.getByTestId("background-note")).toHaveTextContent("I₂ dropped, matched ×2.00");
  });

  it("says when the background has no trace", () => {
    show(set({ trace: "missing", scale: undefined }));
    expect(screen.getByTestId("background-note")).toHaveTextContent("no trace");
  });

  it("says while the background's trace loads", () => {
    show(set({ trace: "loading", scale: undefined }));
    expect(screen.getByTestId("background-note")).toHaveTextContent("loading");
  });

  it("says when the setup of the background is hidden", () => {
    show(set({ trace: "hidden", scale: undefined }));
    expect(screen.getByTestId("background-note")).toHaveTextContent("setup hidden");
  });

  it("reads no background when turned off, and does not offer None again", () => {
    const on = show({ kind: "off", reason: "chosen", canUseRunBuffer: true });
    expect(screen.getByTestId("background-name")).toHaveTextContent("None");
    expect(screen.getByTestId("background-note")).toBeEmptyDOMElement();
    expect(screen.queryByRole("button", { name: "None" })).toBeNull();
    fireEvent.click(screen.getByRole("button", { name: "Run buffer" }));
    expect(on.onDefault).toHaveBeenCalled();
  });

  it("says when the run has no buffer to use, and does not offer None, which would stop a later buffer being the default", () => {
    show({ kind: "off", reason: "no-buffer", canUseRunBuffer: false });
    expect(screen.getByTestId("background-name")).toHaveTextContent("No buffer in this run");
    expect(screen.getByTestId("background-note")).toBeEmptyDOMElement();
    expect(screen.queryByRole("button", { name: "None" })).toBeNull();
    expect(screen.getByRole("button", { name: /Change/ })).toBeInTheDocument();
  });

  const SHORT: Record<CompareScale["source"], string> = {
    normalised: "",
    monitor: "scaled ×1.50",
    matched: "matched at high q ×1.50",
    unknown_normalisation: "matched at high q ×1.50",
    unknown_mode: "matched at high q ×1.50",
    other_setup: "other setup, matched ×1.50",
    beam_off: "beam off, matched ×1.50",
    i2_collapsed: "I₂ dropped, matched ×1.50",
    capillary: "matched at high q ×1.50",
  };

  it.each(Object.keys(SCALE_NOTES) as CompareScale["source"][])("shows %s as a short label and gives the whole sentence as its description", (source) => {
    show(set({ scale: { scale: source === "normalised" ? 1 : 1.5, source } }));
    const note = screen.getByTestId("background-note");
    if (source === "normalised") {
      expect(note).toBeEmptyDOMElement();
      expect(note).not.toHaveAttribute("aria-description");
    } else {
      expect(note).toHaveTextContent(SHORT[source]);
      expect(note).toHaveAttribute("aria-description", SCALE_NOTES[source]);
    }
    expect(screen.queryByText(SCALE_NOTES[source])).toBeNull();
    expect(screen.queryByRole("tooltip")).toBeNull();
  });

  it("tags a buffer and offers no controls", () => {
    show({ kind: "buffer" });
    expect(screen.getByTestId("background-name")).toHaveTextContent("Buffer");
    expect(screen.getByTestId("background-note")).toBeEmptyDOMElement();
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

  it("renders one name row, one note row and one controls row in every state, so the module's height does not depend on the sample", () => {
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
      expect(screen.getAllByTestId("background-note")).toHaveLength(1);
      unmount();
    }
  });
});
