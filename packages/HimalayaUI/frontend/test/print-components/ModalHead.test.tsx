import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { ModalHead } from "../../src/print/components/ModalHead";

describe("<ModalHead>", () => {
  it("renders the kicker text", () => {
    render(<ModalHead kicker="Speculative" title="Custom index" onClose={() => {}} />);
    expect(screen.getByText("Speculative")).toBeInTheDocument();
  });

  it("renders the title text in an h2", () => {
    render(<ModalHead kicker="Speculative" title="Custom index" onClose={() => {}} />);
    const h2 = screen.getByRole("heading", { level: 2 });
    expect(h2).toHaveTextContent("Custom index");
  });

  it("applies titleId to the h2 when given", () => {
    render(
      <ModalHead
        kicker="Speculative"
        title="Custom index"
        titleId="my-modal-title"
        onClose={() => {}}
      />
    );
    const h2 = screen.getByRole("heading", { level: 2 });
    expect(h2).toHaveAttribute("id", "my-modal-title");
  });

  it("does not set an id on the h2 when titleId is omitted", () => {
    render(<ModalHead kicker="Speculative" title="Custom index" onClose={() => {}} />);
    const h2 = screen.getByRole("heading", { level: 2 });
    expect(h2).not.toHaveAttribute("id");
  });

  it("renders a Close button with aria-label='Close'", () => {
    render(<ModalHead kicker="Speculative" title="Custom index" onClose={() => {}} />);
    expect(screen.getByRole("button", { name: "Close" })).toBeInTheDocument();
  });

  it("calls onClose when the Close button is clicked", () => {
    const onClose = vi.fn();
    render(<ModalHead kicker="Speculative" title="Custom index" onClose={onClose} />);
    fireEvent.click(screen.getByRole("button", { name: "Close" }));
    expect(onClose).toHaveBeenCalledTimes(1);
  });

  it("renders with data-testid=modal-head", () => {
    render(<ModalHead kicker="Speculative" title="Custom index" onClose={() => {}} />);
    expect(screen.getByTestId("modal-head")).toBeInTheDocument();
  });
});
