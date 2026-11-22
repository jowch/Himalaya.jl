import { act, fireEvent, render, screen } from "@testing-library/react";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { CompareBar, type CompareChip } from "../../src/ui/CompareBar";
import { resetTooltipWarmth, TOOLTIP_OPEN_DELAY_MS } from "../../src/ui/Tooltip";

const chip = (over: Partial<CompareChip> = {}): CompareChip => ({ id: 2, name: "MO_PEG5", color: "#888", scale: undefined, ...over });

function bar(chips: CompareChip[], subtractId?: number) {
  const noop = vi.fn();
  render(<CompareBar chips={chips} subtractId={subtractId} hidden={false} modKey="⌘" onToggleShown={noop} onSubtract={noop} onRemove={noop} onAdd={noop} />);
}

function hover(name: string): HTMLElement {
  fireEvent.mouseEnter(screen.getByRole("button", { name }));
  act(() => {
    vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS);
  });
  return screen.getByRole("tooltip");
}

describe("CompareBar chip tooltip", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    resetTooltipWarmth();
  });
  afterEach(() => {
    vi.useRealTimers();
  });

  it("names the chip and its setup and says how to hide and subtract", () => {
    bar([chip({ setup: "setup2" })]);
    expect(hover("MO_PEG5, setup2")).toHaveTextContent("MO_PEG5 · setup2Click to hide. ⌘-click to subtract it from the current trace.");
  });

  it("says how to stop subtracting on the subtracted chip, with no setup when it is the current one", () => {
    bar([chip()], 2);
    expect(hover("MO_PEG5")).toHaveTextContent("MO_PEG5Click to hide. ⌘-click to stop subtracting.");
  });

  it("repeats the scale's note on the chip, so keyboard focus reaches it", () => {
    bar([chip({ scale: { scale: 1, source: "monitor" } as CompareChip["scale"] })]);
    expect(hover("MO_PEG5")).toHaveTextContent("Scaled by the transmitted beam");
  });

  it("carries no native title on the chip button", () => {
    bar([chip({ setup: "setup2" })]);
    expect(screen.getByRole("button", { name: "MO_PEG5, setup2" })).not.toHaveAttribute("title");
  });
});

describe("CompareBar scale tooltip", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    resetTooltipWarmth();
  });
  afterEach(() => {
    vi.useRealTimers();
  });

  it("explains a high-q match on the scale, with no native title", () => {
    bar([chip({ scale: { scale: 1, source: "matched" } as CompareChip["scale"] })]);
    const scale = document.querySelector("[data-scale-source]") as HTMLElement;
    expect(scale).not.toHaveAttribute("title");
    fireEvent.mouseEnter(scale);
    act(() => {
      vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS);
    });
    expect(screen.getByRole("tooltip")).toHaveTextContent("No beam monitors, so matched at high q");
  });
});
