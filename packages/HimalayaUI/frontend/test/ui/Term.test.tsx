import { act, fireEvent, render, screen } from "@testing-library/react";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { Term } from "../../src/ui/Term";
import { resetTooltipWarmth, TOOLTIP_OPEN_DELAY_MS } from "../../src/ui/Tooltip";

const EXPLAIN = "Multiplied by 1.10 so the two traces agree at high q, leaving out peaks.";

describe("Term", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    resetTooltipWarmth();
  });
  afterEach(() => {
    vi.useRealTimers();
  });

  const term = (): HTMLElement => screen.getByText("×1.10");

  it("explains itself on hover", () => {
    render(<Term explain={EXPLAIN}>×1.10</Term>);
    fireEvent.mouseEnter(term());
    act(() => {
      vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS);
    });
    expect(screen.getByRole("tooltip")).toHaveTextContent(EXPLAIN);
  });

  it("is a Tab stop", () => {
    render(<Term explain={EXPLAIN}>×1.10</Term>);
    expect(term()).toHaveAttribute("tabindex", "0");
  });

  it("gives a screen reader the same text with the tooltip closed", () => {
    render(<Term explain={EXPLAIN}>×1.10</Term>);
    expect(screen.queryByRole("tooltip")).toBeNull();
    expect(term()).toHaveAccessibleDescription(EXPLAIN);
  });

  it("is described once while its tooltip is open", () => {
    render(<Term explain={EXPLAIN}>×1.10</Term>);
    fireEvent.mouseEnter(term());
    act(() => {
      vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS);
    });
    expect(term().getAttribute("aria-describedby")?.split(" ")).toHaveLength(1);
    expect(term()).toHaveAccessibleDescription(EXPLAIN);
  });

  it("keeps its explanation out of view", () => {
    render(<Term explain={EXPLAIN}>×1.10</Term>);
    expect(screen.getByText(EXPLAIN)).not.toBeVisible();
  });

  it("carries no native title", () => {
    render(<Term explain={EXPLAIN}>×1.10</Term>);
    expect(term()).not.toHaveAttribute("title");
  });
});
