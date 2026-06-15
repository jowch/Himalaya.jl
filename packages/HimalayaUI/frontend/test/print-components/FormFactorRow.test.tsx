import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { FormFactorRow } from "../../src/print/components/FormFactorRow";

describe("FormFactorRow", () => {
  it("renders the declaration label + sub-note and reflects unselected state", () => {
    render(<FormFactorRow selected={false} onToggle={() => {}} />);
    const row = screen.getByTestId("form-factor-row");
    expect(row).toHaveTextContent("Form factor");
    expect(row).toHaveTextContent("No Bragg peaks to index");
    expect(row).toHaveAttribute("aria-pressed", "false");
  });

  it("marks itself selected when the sample is declared form factor", () => {
    render(<FormFactorRow selected onToggle={() => {}} />);
    const row = screen.getByTestId("form-factor-row");
    expect(row).toHaveAttribute("aria-pressed", "true");
    // The accessible name distinguishes the declared state for SR users.
    expect(row).toHaveAccessibleName(/declared/);
  });

  it("fires onToggle on click", () => {
    const onToggle = vi.fn();
    render(<FormFactorRow selected={false} onToggle={onToggle} />);
    fireEvent.click(screen.getByTestId("form-factor-row"));
    expect(onToggle).toHaveBeenCalledTimes(1);
  });
});
