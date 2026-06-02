import { render, screen, fireEvent, within } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { PhaseBlock } from "../../src/print/components/PhaseBlock";

describe("<PhaseBlock> rendering", () => {
  it("renders the root [data-testid=phase-block]", () => {
    render(<PhaseBlock phase="Pn3m" score={0.92} meta="a = 14.2 nm · 5 reflections" />);
    expect(screen.getByTestId("phase-block")).toBeInTheDocument();
  });

  it("renders the phase name inside a [data-phase-label]", () => {
    render(<PhaseBlock phase="Pn3m" score={0.92} meta="a = 14.2 nm · 5 reflections" />);
    const name = screen.getByText("Pn3m");
    expect(name.closest("[data-phase-label]")).not.toBeNull();
  });

  it("renders the score formatted to 2 decimals inside a [data-phase-label]", () => {
    render(<PhaseBlock phase="Pn3m" score={0.92} meta="a = 14.2 nm · 5 reflections" />);
    const val = screen.getByText("0.92");
    expect(val.closest("[data-phase-label]")).not.toBeNull();
  });

  it("formats whole-ish scores to 2 decimals", () => {
    render(<PhaseBlock phase="Pn3m" score={0.4} meta="weak" />);
    expect(screen.getByText("0.40")).toBeInTheDocument();
  });

  it("renders the meta text", () => {
    render(<PhaseBlock phase="Pn3m" score={0.92} meta="a = 14.2 nm · 5 reflections" />);
    expect(screen.getByText("a = 14.2 nm · 5 reflections")).toBeInTheDocument();
  });

  it("renders a ScoreBar ([data-score-bar]) with data-phase equal to the phase", () => {
    const { container } = render(
      <PhaseBlock phase="Im3m" score={0.71} meta="a = 9.1 nm · 4 reflections" />,
    );
    const bar = container.querySelector("[data-score-bar]");
    expect(bar).not.toBeNull();
    expect(bar).toHaveAttribute("data-phase", "Im3m");
  });
});

describe("<PhaseBlock> remove control", () => {
  it("renders a remove button with accessible name 'Remove <phase>'", () => {
    render(<PhaseBlock phase="Pn3m" score={0.92} meta="meta" onRemove={vi.fn()} />);
    expect(screen.getByRole("button", { name: /Remove Pn3m/ })).toBeInTheDocument();
  });

  it("fires onRemove when the remove button is clicked", () => {
    const onRemove = vi.fn();
    render(<PhaseBlock phase="Pn3m" score={0.92} meta="meta" onRemove={onRemove} />);
    fireEvent.click(screen.getByRole("button", { name: /Remove Pn3m/ }));
    expect(onRemove).toHaveBeenCalledTimes(1);
  });

  it("always renders the remove button even without onRemove; clicking is a no-op", () => {
    render(<PhaseBlock phase="Pn3m" score={0.92} meta="meta" />);
    const btn = screen.getByRole("button", { name: /Remove Pn3m/ });
    expect(btn).toBeInTheDocument();
    expect(() => fireEvent.click(btn)).not.toThrow();
  });
});

describe("<PhaseBlock> placement className", () => {
  it("appends a placement className to the root", () => {
    render(<PhaseBlock phase="Pn3m" score={0.92} meta="meta" className="border-t" />);
    expect(screen.getByTestId("phase-block").className).toContain("border-t");
  });

  it("scopes the phase label query to this block's root", () => {
    render(<PhaseBlock phase="Lamellar" score={0.42} meta="a = 6.0 nm · 2 reflections" />);
    const root = screen.getByTestId("phase-block");
    expect(within(root).getByText("Lamellar")).toBeInTheDocument();
  });
});
