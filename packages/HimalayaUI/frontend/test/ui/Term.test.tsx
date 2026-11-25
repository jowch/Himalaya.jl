import { act, fireEvent, render, screen } from "@testing-library/react";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { Explanation, Term } from "../../src/ui/Term";
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

describe("Term inside a control", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    resetTooltipWarmth();
  });
  afterEach(() => {
    vi.useRealTimers();
  });

  it("explains itself on hover but is no second Tab stop, and the control carries the explanation", () => {
    render(
      <button type="button" aria-describedby="why">
        Pn3m <Term explain={EXPLAIN} nested>0.87</Term>
        <Explanation id="why" texts={[EXPLAIN]} />
      </button>,
    );
    const term = screen.getByText("0.87");
    expect(term).not.toHaveAttribute("tabindex");
    expect(term).not.toHaveAttribute("aria-describedby");
    expect(screen.getByRole("button")).toHaveAccessibleDescription(EXPLAIN);
    fireEvent.mouseEnter(term);
    act(() => {
      vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS);
    });
    expect(screen.getByRole("tooltip")).toHaveTextContent(EXPLAIN);
  });

  it("joins several explanations into one description", () => {
    render(
      <button type="button" aria-describedby="why">
        x
        <Explanation id="why" texts={["One.", "Two."]} />
      </button>,
    );
    expect(screen.getByRole("button")).toHaveAccessibleDescription("One. Two.");
  });
});
