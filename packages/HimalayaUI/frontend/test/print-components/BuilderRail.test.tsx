import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { BuilderRail } from "../../src/print/components/BuilderRail";

const base = {
  orderedBy: "LL37 : lipid ratio",
  offset: 1.2,
  onOffsetChange: () => {},
  traces: <div data-testid="traces-slot" />,
};

describe("<BuilderRail>", () => {
  it("renders the compose head, field, display controls, and traces slot", () => {
    render(<BuilderRail {...base} />);
    expect(screen.getByTestId("builder-rail")).toBeInTheDocument();
    expect(screen.getByText("Compose")).toBeInTheDocument();
    expect(screen.getByTestId("field")).toHaveTextContent("LL37 : lipid ratio");
    expect(screen.getByTestId("slider")).toBeInTheDocument();
    expect(screen.getByTestId("traces-slot")).toBeInTheDocument();
  });
  it("ditches the grouping-summary card (BU-EDIT-BUTTON): no auto-group, no 'N samples' summary", () => {
    // The "N samples, ordered by …" compose card was redundant with the
    // Ordering-variable field + the Traces list and confused the edit flow.
    // It is gone; the rail leads with the action buttons instead.
    render(<BuilderRail {...base} onAdjust={() => {}} />);
    expect(screen.queryByTestId("auto-group")).toBeNull();
    expect(screen.queryByText(/samples, ordered by/i)).toBeNull();
  });
  it("does NOT render 'Copy as PNG' (figure export lives on the plate head)", () => {
    // Same single-contextual-control reasoning as the q-scale toggle: the
    // export affordance is the ExportButton in SeriesPlate's actions slot.
    render(<BuilderRail {...base} />);
    expect(screen.queryByText(/copy as png/i)).toBeNull();
  });
  it("does NOT render 'Copy as PNG' even when the foot row renders (onAddSample passed)", () => {
    render(<BuilderRail {...base} onAddSample={() => {}} />);
    expect(screen.queryByText(/copy as png/i)).toBeNull();
  });
  it("OMITS the collapse and add-sample controls when their handlers are withheld", () => {
    // controls-don't-lie: the live page passes neither onCollapse (no collapse
    // behavior exists) nor onAddSample (its add path is the traces-slot select),
    // so neither control may render inert.
    render(<BuilderRail {...base} />);
    expect(screen.queryByRole("button", { name: /collapse rail/i })).toBeNull();
    expect(screen.queryByRole("button", { name: /add sample/i })).toBeNull();
  });
  it("renders and fires the collapse and add-sample controls when wired", () => {
    const onAddSample = vi.fn(), onCollapse = vi.fn();
    render(<BuilderRail {...base} onAddSample={onAddSample} onCollapse={onCollapse} />);
    fireEvent.click(screen.getByRole("button", { name: /add sample/i }));
    expect(onAddSample).toHaveBeenCalledOnce();
    fireEvent.click(screen.getByRole("button", { name: /collapse rail/i }));
    expect(onCollapse).toHaveBeenCalledOnce();
  });
  it("OMITS the heatmap representation toggle and the track-reflections switch", () => {
    render(<BuilderRail {...base} />);
    expect(screen.queryByText(/heatmap/i)).toBeNull();
    expect(screen.queryByText(/track reflections/i)).toBeNull();
  });
  it("does NOT render a q-scale toggle (it lives on the plate, not the rail)", () => {
    render(<BuilderRail {...base} />);
    // The single q-scale control is contextual to the figure (the plate); the
    // rail's redundant copy was removed.
    expect(screen.queryByRole("button", { name: /log q/i })).toBeNull();
    expect(screen.queryByRole("button", { name: /linear q/i })).toBeNull();
  });
  it("names the static ordering field with its label so it is not a bare value (WCAG 4.1.2)", () => {
    // The builder page passes neither orderOptions nor onChangeOrder, so the
    // Field renders in static mode; its readable text must still compose
    // "Ordering variable {value}" and never drop the value.
    render(<BuilderRail {...base} />);
    expect(screen.getByTestId("field")).toHaveTextContent(
      /Ordering variable\s+LL37 : lipid ratio/,
    );
  });
  it("names the dropdown ordering field with label + value (WCAG 4.1.2)", () => {
    render(
      <BuilderRail
        {...base}
        orderOptions={[
          { value: "LL37 : lipid ratio", label: "LL37 : lipid ratio" },
          { value: "Time", label: "Time" },
        ]}
        onOrderSelect={() => {}}
      />,
    );
    // With options the field is a real dropdown trigger button; its accessible
    // name must read label + value, not just the value.
    expect(
      screen.getByRole("button", { name: /ordering variable\s+LL37 : lipid ratio/i }),
    ).toBeInTheDocument();
  });
  it("renders the Edit action when onAdjust is provided (read state)", () => {
    render(<BuilderRail {...base} onAdjust={() => {}} />);
    expect(screen.getByRole("button", { name: /^edit$/i })).toBeInTheDocument();
    // BU-AUTOGROUP-STALE: real verbs, not "Adjust"/"Confirm series".
    expect(screen.queryByRole("button", { name: /adjust/i })).toBeNull();
    expect(screen.queryByRole("button", { name: /confirm series/i })).toBeNull();
  });
  it("OMITS the Edit action when onAdjust is withheld (draft live), shows Save changes + Cancel as buttons", () => {
    // controls-don't-lie: a live draft withholds onAdjust → no redundant Edit.
    // BU-DRAFT-ACTIONS: a live draft (reorderable) presents Save changes + Cancel
    // as a proper button row, not a card link + a buried ghost Cancel.
    render(<BuilderRail {...base} reorderable onConfirm={() => {}} onCancel={() => {}} />);
    expect(screen.queryByRole("button", { name: /^edit$/i })).toBeNull();
    expect(screen.getByTestId("builder-save")).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /save changes/i })).toBeEnabled();
    expect(screen.getByRole("button", { name: /^cancel$/i })).toBeInTheDocument();
  });
  it("renders the default WYSIWYG foot caption when the caption prop is omitted", () => {
    render(<BuilderRail {...base} />);
    expect(
      screen.getByText(
        "The plate above is the figure as it will export. What you compose is what you publish.",
      ),
    ).toBeInTheDocument();
  });
  it("a provided caption REPLACES the default WYSIWYG sentence (copy-doesn't-lie override)", () => {
    render(<BuilderRail {...base} caption="The plate still shows the last confirmed figure." />);
    expect(
      screen.getByText("The plate still shows the last confirmed figure."),
    ).toBeInTheDocument();
    expect(
      screen.queryByText(/what you compose is what you publish/i),
    ).toBeNull();
  });
  it("fires the offset-slider handler", () => {
    const onOffsetChange = vi.fn();
    render(<BuilderRail {...base} onOffsetChange={onOffsetChange} />);
    fireEvent.change(screen.getByTestId("slider"), { target: { value: "0.8" } });
    expect(onOffsetChange).toHaveBeenCalled();
  });

  // ── Traces section label tells the truth about reorderability ─────────────
  it("by default the Traces section makes NO drag-to-reorder promise (read mode)", () => {
    // The read-mode traces slot is a non-reorderable MemberList; the section
    // label must not instruct a drag that does nothing (copy-doesn't-lie).
    render(<BuilderRail {...base} />);
    expect(screen.getByText("Traces")).toBeInTheDocument();
    expect(screen.queryByText(/drag to reorder/i)).toBeNull();
  });
  it("when reorderable, the Traces section label states 'drag to reorder' (draft mode)", () => {
    render(<BuilderRail {...base} reorderable />);
    expect(screen.getByText("Traces · drag to reorder")).toBeInTheDocument();
  });
});
