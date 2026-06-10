import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { BuilderRail } from "../../src/print/components/BuilderRail";

const base = {
  grouping: <>Himalaya read <strong>6 samples</strong> as one series.</>,
  orderedBy: "LL37 : lipid ratio",
  offset: 1.2,
  onOffsetChange: () => {},
  traces: <div data-testid="traces-slot" />,
};

describe("<BuilderRail>", () => {
  it("renders the compose head, autogroup, field, display controls, traces slot, and foot", () => {
    render(<BuilderRail {...base} />);
    expect(screen.getByTestId("builder-rail")).toBeInTheDocument();
    expect(screen.getByText("Compose")).toBeInTheDocument();
    expect(screen.getByText("Auto-grouped")).toBeInTheDocument();
    expect(screen.getByTestId("field")).toHaveTextContent("LL37 : lipid ratio");
    expect(screen.getByTestId("slider")).toBeInTheDocument();
    expect(screen.getByTestId("traces-slot")).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /add sample/i })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /copy as png/i })).toBeInTheDocument();
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
  it("renders the Adjust action when onAdjust is provided (read state)", () => {
    render(<BuilderRail {...base} onAdjust={() => {}} />);
    expect(screen.getByRole("button", { name: /adjust/i })).toBeInTheDocument();
  });
  it("OMITS the Adjust action when onAdjust is withheld (draft live)", () => {
    // controls-don't-lie: a live draft withholds onAdjust → no redundant Adjust.
    render(<BuilderRail {...base} onConfirm={() => {}} />);
    expect(screen.queryByRole("button", { name: /adjust/i })).toBeNull();
    expect(screen.getByRole("button", { name: /confirm series/i })).toBeInTheDocument();
  });
  it("fires the display + foot + collapse handlers", () => {
    const onOffsetChange = vi.fn(), onAddSample = vi.fn(), onCollapse = vi.fn();
    render(<BuilderRail {...base} onOffsetChange={onOffsetChange}
      onAddSample={onAddSample} onCollapse={onCollapse} />);
    fireEvent.change(screen.getByTestId("slider"), { target: { value: "0.8" } });
    expect(onOffsetChange).toHaveBeenCalled();
    fireEvent.click(screen.getByRole("button", { name: /add sample/i }));
    expect(onAddSample).toHaveBeenCalledOnce();
    fireEvent.click(screen.getByRole("button", { name: /collapse rail/i }));
    expect(onCollapse).toHaveBeenCalledOnce();
  });
});
