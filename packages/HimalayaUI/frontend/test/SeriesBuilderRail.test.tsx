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
    exportControls: <div data-testid="mock-export" />,
    ...over,
  };
  render(<SeriesBuilderRail {...props} />);
  return props;
}

describe("SeriesBuilderRail", () => {
  it("renders the ordering variable and the representation toggle", () => {
    setup();
    expect(screen.getByText("LL37 : lipid ratio")).toBeInTheDocument();
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
});
