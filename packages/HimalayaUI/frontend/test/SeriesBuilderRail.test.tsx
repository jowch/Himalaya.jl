import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { SeriesBuilderRail } from "../src/components/SeriesBuilderRail";

function setup(over: Partial<React.ComponentProps<typeof SeriesBuilderRail>> = {}) {
  const props = {
    collapsed: false,
    onToggleCollapsed: vi.fn(),
    representation: "waterfall" as const,
    onRepresentationChange: vi.fn(),
    orderingVariable: "LL37 : lipid ratio",
    offset: 1.2,
    onOffsetChange: vi.fn(),
    scaleMode: "log" as const,
    onScaleModeChange: vi.fn(),
    trackOn: false,
    onTrackOnChange: vi.fn(),
    sampleCount: 6,
    onConfirmSeries: vi.fn(),
    onAdjustSeries: vi.fn(),
    exportControls: <div data-testid="mock-export" />,
    ...over,
  };
  render(<SeriesBuilderRail {...props} />);
  return props;
}

describe("SeriesBuilderRail", () => {
  it("renders the ordering variable and the representation toggle", () => {
    setup();
    // The ordering variable appears in both the autogroup card body and the
    // dedicated ordering-variable line — assert on the dedicated line.
    expect(screen.getAllByText("LL37 : lipid ratio").length).toBeGreaterThan(0);
    expect(screen.getByTestId("representation-toggle")).toBeInTheDocument();
  });

  it("renders the injected export controls", () => {
    setup();
    expect(screen.getByTestId("mock-export")).toBeInTheDocument();
  });

  it("calls onToggleCollapsed when the collapse button is clicked", () => {
    const props = setup();
    fireEvent.click(screen.getByTestId("rail-collapse-toggle"));
    expect(props.onToggleCollapsed).toHaveBeenCalled();
  });

  it("hides the rail body but keeps the restore tab when collapsed", () => {
    setup({ collapsed: true });
    expect(screen.queryByText("LL37 : lipid ratio")).not.toBeInTheDocument();
    expect(screen.getByTestId("rail-restore")).toBeInTheDocument();
  });

  it("renders the offset slider, scale toggle, and autogroup card", () => {
    setup();
    expect(screen.getByTestId("offset-slider")).toBeInTheDocument();
    expect(screen.getByTestId("scale-toggle")).toBeInTheDocument();
    expect(screen.getByTestId("autogroup-card")).toBeInTheDocument();
  });

  it("forwards offset slider input to onOffsetChange", () => {
    const props = setup();
    fireEvent.change(screen.getByTestId("offset-slider"), { target: { value: "0.8" } });
    expect(props.onOffsetChange).toHaveBeenCalledWith(0.8);
  });

  it("forwards scale toggle to onScaleModeChange", () => {
    const props = setup();
    fireEvent.click(screen.getByTestId("scale-linear"));
    expect(props.onScaleModeChange).toHaveBeenCalledWith("linear");
  });

  it("hides the autogroup card and shows the recipe editor in edit mode", () => {
    setup({ editControls: <div data-testid="mock-edit" /> });
    expect(screen.queryByTestId("autogroup-card")).not.toBeInTheDocument();
    expect(screen.getByTestId("mock-edit")).toBeInTheDocument();
    // compose controls stay available in edit mode (they shape the figure too)
    expect(screen.getByTestId("offset-slider")).toBeInTheDocument();
  });

  it("renders the Track-reflections checkbox in the Display section and forwards toggles", () => {
    const props = setup();
    const toggle = screen.getByTestId("track-toggle-input") as HTMLInputElement;
    expect(toggle).toBeInTheDocument();
    expect(toggle.checked).toBe(false);
    fireEvent.click(toggle);
    expect(props.onTrackOnChange).toHaveBeenCalledWith(true);
  });

  it("reflects trackOn=true as a checked checkbox", () => {
    setup({ trackOn: true });
    expect((screen.getByTestId("track-toggle-input") as HTMLInputElement).checked).toBe(true);
  });

  it("recedes the COMPOSE header behind the plate (R3-Y03)", () => {
    setup();
    // The header is a recessed label, not a competing title: it carries the
    // receded data marker the component sets when it drops to ink-faint.
    expect(screen.getByTestId("rail-compose-header")).toHaveAttribute(
      "data-recede",
      "true",
    );
  });

  it("renders the Track-reflections checkbox in the terracotta accent (R3-Y04)", () => {
    setup();
    expect(screen.getByTestId("track-toggle-input")).toHaveAttribute(
      "data-accent",
      "print",
    );
  });

  it("scopes edit-input styling so it does not target slider thumbs (R3-Y09)", () => {
    setup({ editControls: <input data-testid="ec" /> });
    expect(screen.getByTestId("rail-edit")).toHaveAttribute("data-rail-edit-inputs", "");
  });
});
