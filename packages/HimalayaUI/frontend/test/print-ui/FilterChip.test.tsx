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
});
