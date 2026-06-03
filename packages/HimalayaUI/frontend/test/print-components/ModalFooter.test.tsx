import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { ModalFooter } from "../../src/print/components/ModalFooter";

describe("<ModalFooter>", () => {
  it("renders with data-testid=modal-foot", () => {
    render(<ModalFooter actions={<button>OK</button>} />);
    expect(screen.getByTestId("modal-foot")).toBeInTheDocument();
  });

  it("renders note text when provided", () => {
    render(
      <ModalFooter note="This cannot be undone." actions={<button>OK</button>} />
    );
    expect(screen.getByTestId("modal-foot")).toHaveTextContent(
      "This cannot be undone."
    );
  });

  it("renders actions children", () => {
    render(
      <ModalFooter
        actions={
          <>
            <button>Cancel</button>
            <button>Add</button>
          </>
        }
      />
    );
    expect(screen.getByRole("button", { name: "Cancel" })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "Add" })).toBeInTheDocument();
  });

  it("renders both note and actions together", () => {
    render(
      <ModalFooter
        note="Some note"
        actions={<button>Confirm</button>}
      />
    );
    expect(screen.getByTestId("modal-foot")).toHaveTextContent("Some note");
    expect(screen.getByRole("button", { name: "Confirm" })).toBeInTheDocument();
  });

  it("renders without note (no crash)", () => {
    expect(() =>
      render(<ModalFooter actions={<button>OK</button>} />)
    ).not.toThrow();
  });
});
