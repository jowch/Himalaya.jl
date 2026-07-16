import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { FilterChip } from "../../src/print/ui/FilterChip";

describe("FilterChip", () => {
  it("renders the label text", () => {
    render(<FilterChip label="Coexistence" active={false} onClick={() => {}} />);
    expect(screen.getByTestId("filter-chip").textContent).toContain(
      "Coexistence",
    );
  });

  it('defaults to data-size="md"', () => {
    render(<FilterChip label="Coexistence" active={false} onClick={() => {}} />);
    expect(screen.getByTestId("filter-chip").getAttribute("data-size")).toBe(
      "md",
    );
  });

  it("reflects active=true via aria-pressed and data-active", () => {
    render(<FilterChip label="Coexistence" active={true} onClick={() => {}} />);
    const chip = screen.getByTestId("filter-chip");
    expect(chip.getAttribute("aria-pressed")).toBe("true");
    expect(chip.getAttribute("data-active")).toBe("true");
  });

  it("reflects active=false via aria-pressed and data-active", () => {
    render(<FilterChip label="Coexistence" active={false} onClick={() => {}} />);
    const chip = screen.getByTestId("filter-chip");
    expect(chip.getAttribute("aria-pressed")).toBe("false");
    expect(chip.getAttribute("data-active")).toBe("false");
  });

  it("calls onClick when clicked", () => {
    const onClick = vi.fn();
    render(
      <FilterChip label="Coexistence" active={false} onClick={onClick} />,
    );
    fireEvent.click(screen.getByTestId("filter-chip"));
    expect(onClick).toHaveBeenCalledTimes(1);
  });

  it("threads an optional title tooltip onto the chip (F5)", () => {
    render(
      <FilterChip
        label="Coexistence"
        active={false}
        onClick={() => {}}
        title="Series whose members span more than one phase"
      />,
    );
    expect(screen.getByTestId("filter-chip")).toHaveAttribute(
      "title",
      "Series whose members span more than one phase",
    );
  });

  it("omits the title attribute when no tooltip is given (byte-identical default)", () => {
    render(<FilterChip label="Coexistence" active={false} onClick={() => {}} />);
    expect(screen.getByTestId("filter-chip")).not.toHaveAttribute("title");
  });

  it("exposes an optional description to AT via aria-describedby + a visually-hidden span, keeping the NAME the label (FIX 2)", () => {
    render(
      <FilterChip
        label="Has transition"
        active={false}
        onClick={() => {}}
        description="Series whose members span more than one phase"
      />,
    );
    // Accessible NAME stays the label (the description must NOT bleed into it).
    const chip = screen.getByRole("button", { name: "Has transition" });
    // The description is wired and reachable.
    const descId = chip.getAttribute("aria-describedby");
    expect(descId).toBeTruthy();
    const desc = document.getElementById(descId!);
    expect(desc).not.toBeNull();
    expect(desc).toHaveTextContent("Series whose members span more than one phase");
    // The description carrier is visually hidden, not laid out as visible text.
    expect(desc).toHaveClass("sr-only");
  });

  it("omits aria-describedby when no description is given (byte-identical default)", () => {
    render(<FilterChip label="Coexistence" active={false} onClick={() => {}} />);
    expect(screen.getByTestId("filter-chip")).not.toHaveAttribute("aria-describedby");
  });
});
