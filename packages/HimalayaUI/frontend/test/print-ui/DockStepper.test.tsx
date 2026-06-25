import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { DockStepper } from "../../src/print/ui/DockStepper";

const base = {
  label: "Sample",
  axis: "vertical",
  testIdBase: "sample",
  onPrev: () => {},
  onNext: () => {},
} as const;

describe("<DockStepper>", () => {
  it("renders the label, vertical glyphs, and stable testids", () => {
    render(<DockStepper {...base} count="3 / 12" />);
    expect(screen.getByText("Sample")).toBeInTheDocument();
    const prev = screen.getByTestId("dock-prev-sample");
    const next = screen.getByTestId("dock-next-sample");
    expect(prev).toHaveTextContent("↑");
    expect(next).toHaveTextContent("↓");
    expect(prev).toHaveAttribute("aria-label", "Previous sample");
    expect(next).toHaveAttribute("aria-label", "Next sample");
    expect(screen.getByTestId("dock-sample-count")).toHaveTextContent("3 / 12");
  });

  it("uses horizontal glyphs for the frame axis", () => {
    render(
      <DockStepper label="Frame" axis="horizontal" testIdBase="frame" onPrev={() => {}} onNext={() => {}} count="2 / 5" />,
    );
    expect(screen.getByTestId("dock-prev-frame")).toHaveTextContent("←");
    expect(screen.getByTestId("dock-next-frame")).toHaveTextContent("→");
    expect(screen.getByTestId("dock-frame-count")).toBeInTheDocument();
  });

  it("hides the readout when count is omitted", () => {
    render(<DockStepper {...base} />);
    expect(screen.queryByTestId("dock-sample-count")).not.toBeInTheDocument();
  });

  it("disables prev/next per props and fires the right handler", () => {
    const onPrev = vi.fn();
    const onNext = vi.fn();
    render(<DockStepper {...base} onPrev={onPrev} onNext={onNext} prevDisabled count="1 / 9" />);
    expect(screen.getByTestId("dock-prev-sample")).toBeDisabled();
    fireEvent.click(screen.getByTestId("dock-next-sample"));
    expect(onNext).toHaveBeenCalledTimes(1);
    expect(onPrev).not.toHaveBeenCalled();
  });
});
