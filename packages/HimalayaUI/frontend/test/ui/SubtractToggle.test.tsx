import { act, fireEvent, render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { SubtractToggle } from "../../src/ui/SubtractToggle";
import { minusExplanation, RAW } from "../../src/ui/explanations";
import { resetTooltipWarmth, TOOLTIP_OPEN_DELAY_MS } from "../../src/ui/Tooltip";

describe("SubtractToggle", () => {
  beforeEach(() => {
    resetTooltipWarmth();
  });

  it("explains Raw and Minus on their buttons", () => {
    render(<SubtractToggle name="buffer_blank" on={false} onChange={() => {}} />);
    expect(screen.getByRole("button", { name: "Raw" })).toHaveAccessibleDescription(RAW);
    expect(screen.getByRole("button", { name: "Minus buffer_blank" })).toHaveAccessibleDescription(minusExplanation("buffer_blank"));
  });

  it("adds no Tab stops of its own: Tab goes Raw, then Minus", async () => {
    render(<><SubtractToggle name="buffer_blank" on={false} onChange={() => {}} /><button type="button">after</button></>);
    await userEvent.tab();
    expect(screen.getByRole("button", { name: "Raw" })).toHaveFocus();
    await userEvent.tab();
    expect(screen.getByRole("button", { name: "Minus buffer_blank" })).toHaveFocus();
    await userEvent.tab();
    expect(screen.getByRole("button", { name: "after" })).toHaveFocus();
  });

  it("shows Minus's explanation when the keyboard reaches it", async () => {
    render(<SubtractToggle name="buffer_blank" on={false} onChange={() => {}} />);
    await userEvent.tab();
    await userEvent.tab();
    expect(screen.getByRole("tooltip")).toHaveTextContent(minusExplanation("buffer_blank"));
  });
});

describe("SubtractToggle on hover", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    resetTooltipWarmth();
  });
  afterEach(() => {
    vi.useRealTimers();
  });

  it("explains Raw on hover", () => {
    render(<SubtractToggle name="buffer_blank" on onChange={() => {}} />);
    fireEvent.mouseEnter(screen.getByText("Raw"));
    act(() => {
      vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS);
    });
    expect(screen.getByRole("tooltip")).toHaveTextContent(RAW);
  });
});
