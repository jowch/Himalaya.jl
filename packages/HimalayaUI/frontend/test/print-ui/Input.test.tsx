import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Input } from "../../src/print/ui/Input";

describe("<Input>", () => {
  it("reflects value", () => {
    render(<Input value="lipid" onValueChange={() => {}} aria-label="name" />);
    const input = screen.getByLabelText("name") as HTMLInputElement;
    expect(input.value).toBe("lipid");
  });

  it("fires onValueChange with the new string on typing", () => {
    const onValueChange = vi.fn();
    render(<Input value="" onValueChange={onValueChange} aria-label="name" />);
    const input = screen.getByLabelText("name");
    fireEvent.change(input, { target: { value: "abc" } });
    expect(onValueChange).toHaveBeenCalledWith("abc");
  });

  it("sets data-invalid and aria-invalid only when invalid", () => {
    const { rerender } = render(
      <Input value="" onValueChange={() => {}} aria-label="name" />,
    );
    expect(screen.getByTestId("input").getAttribute("data-invalid")).toBeNull();
    expect(screen.getByLabelText("name").getAttribute("aria-invalid")).toBeNull();

    rerender(<Input value="" onValueChange={() => {}} aria-label="name" invalid />);
    expect(screen.getByTestId("input").getAttribute("data-invalid")).toBe("true");
    expect(screen.getByLabelText("name").getAttribute("aria-invalid")).toBe("true");
  });

  it("renders leading and trailing slot nodes", () => {
    render(
      <Input
        value=""
        onValueChange={() => {}}
        aria-label="name"
        leading={<span data-testid="lead">L</span>}
        trailing={<span data-testid="trail">T</span>}
      />,
    );
    expect(screen.getByTestId("lead")).toBeInTheDocument();
    expect(screen.getByTestId("trail")).toBeInTheDocument();
  });

  it("forwards placeholder via ...rest", () => {
    render(
      <Input value="" onValueChange={() => {}} placeholder="Search series" />,
    );
    expect(screen.getByPlaceholderText("Search series")).toBeInTheDocument();
  });

  it("overrides the wrapper data-testid via testId, default 'input'", () => {
    const { rerender } = render(
      <Input value="" onValueChange={() => {}} aria-label="name" />,
    );
    expect(screen.getByTestId("input")).toBeInTheDocument();

    rerender(
      <Input value="" onValueChange={() => {}} aria-label="name" testId="x" />,
    );
    const wrapper = screen.getByTestId("x");
    expect(wrapper).toBeInTheDocument();
    // testId must not leak onto the inner input
    expect(wrapper).toContainElement(screen.getByLabelText("name"));
    expect(screen.queryByTestId("input")).toBeNull();
  });
});
