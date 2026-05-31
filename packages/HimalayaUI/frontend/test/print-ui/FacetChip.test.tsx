import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { FacetChip } from "../../src/print/ui/FacetChip";

describe("FacetChip", () => {
  it("renders the label text", () => {
    render(<FacetChip label="Beamtime" />);
    expect(screen.getByTestId("facet-chip").textContent).toContain("Beamtime");
  });

  it('has data-testid="facet-chip"', () => {
    render(<FacetChip label="Beamtime" />);
    expect(screen.getByTestId("facet-chip")).toBeTruthy();
  });

  it('defaults to data-size="md"', () => {
    render(<FacetChip label="Beamtime" />);
    expect(screen.getByTestId("facet-chip").getAttribute("data-size")).toBe(
      "md",
    );
  });

  it("calls onClick when clicked", () => {
    const onClick = vi.fn();
    render(<FacetChip label="Beamtime" onClick={onClick} />);
    fireEvent.click(screen.getByTestId("facet-chip"));
    expect(onClick).toHaveBeenCalledTimes(1);
  });

  it("renders an aria-hidden chevron with textContent ▾", () => {
    render(<FacetChip label="Beamtime" />);
    const chevron = screen
      .getByTestId("facet-chip")
      .querySelector('[aria-hidden="true"]');
    expect(chevron).toBeTruthy();
    expect(chevron?.textContent).toBe("▾");
  });
});
