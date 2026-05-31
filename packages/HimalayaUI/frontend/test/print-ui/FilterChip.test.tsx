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
});
