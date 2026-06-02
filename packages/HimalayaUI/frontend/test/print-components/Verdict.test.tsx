import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Verdict } from "../../src/print/components/Verdict";

describe("<Verdict> kept state (dropped=false)", () => {
  it('shows "Kept" state text', () => {
    render(<Verdict dropped={false} />);
    expect(screen.getByText("Kept")).toBeInTheDocument();
  });

  it("shows the default kept hint", () => {
    render(<Verdict dropped={false} />);
    expect(
      screen.getByText("Everything is kept until you drop it.")
    ).toBeInTheDocument();
  });

  it("action button is labeled /drop/i", () => {
    render(<Verdict dropped={false} onToggle={() => {}} />);
    expect(screen.getByRole("button", { name: /drop/i })).toBeInTheDocument();
  });

  it("clicking the action button calls onToggle", () => {
    const onToggle = vi.fn();
    render(<Verdict dropped={false} onToggle={onToggle} />);
    fireEvent.click(screen.getByRole("button", { name: /drop/i }));
    expect(onToggle).toHaveBeenCalledTimes(1);
  });

  it("status Dot has data-tone=success when kept", () => {
    const { container } = render(<Verdict dropped={false} />);
    expect(container.querySelector("[data-tone='success']")).toBeInTheDocument();
  });
});

describe("<Verdict> dropped state (dropped=true)", () => {
  it('shows "Dropped" state text', () => {
    render(<Verdict dropped={true} />);
    expect(screen.getByText("Dropped")).toBeInTheDocument();
  });

  it("shows the default dropped hint", () => {
    render(<Verdict dropped={true} />);
    expect(
      screen.getByText("Restore to keep this exposure.")
    ).toBeInTheDocument();
  });

  it("action button is labeled /restore/i", () => {
    render(<Verdict dropped={true} onToggle={() => {}} />);
    expect(screen.getByRole("button", { name: /restore/i })).toBeInTheDocument();
  });

  it("clicking the action button calls onToggle", () => {
    const onToggle = vi.fn();
    render(<Verdict dropped={true} onToggle={onToggle} />);
    fireEvent.click(screen.getByRole("button", { name: /restore/i }));
    expect(onToggle).toHaveBeenCalledTimes(1);
  });

  it("status Dot has data-tone=accent when dropped", () => {
    const { container } = render(<Verdict dropped={true} />);
    expect(container.querySelector("[data-tone='accent']")).toBeInTheDocument();
  });
});

describe("<Verdict> hint override", () => {
  it("renders a custom hint when provided", () => {
    render(<Verdict dropped={false} hint="Custom hint text." />);
    expect(screen.getByText("Custom hint text.")).toBeInTheDocument();
    expect(
      screen.queryByText("Everything is kept until you drop it.")
    ).not.toBeInTheDocument();
  });
});

describe("<Verdict> no onToggle", () => {
  it("does not throw when clicked without onToggle", () => {
    render(<Verdict dropped={false} />);
    expect(() =>
      fireEvent.click(screen.getByRole("button", { name: /drop/i }))
    ).not.toThrow();
  });
});
