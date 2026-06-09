// test/Checkbox.test.tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { Checkbox } from "../src/print/ui/Checkbox";

describe("Checkbox", () => {
  it("renders unchecked by default with correct aria-checked", () => {
    render(<Checkbox aria-label="select sample" />);
    const cb = screen.getByRole("checkbox", { name: "select sample" });
    expect(cb).toHaveAttribute("aria-checked", "false");
    expect(cb).toHaveAttribute("data-checked", "false");
  });

  it("renders checked state", () => {
    render(<Checkbox checked aria-label="select sample" />);
    const cb = screen.getByRole("checkbox", { name: "select sample" });
    expect(cb).toHaveAttribute("aria-checked", "true");
    expect(cb).toHaveAttribute("data-checked", "true");
  });

  it("renders indeterminate state", () => {
    render(<Checkbox indeterminate aria-label="select sample" />);
    const cb = screen.getByRole("checkbox", { name: "select sample" });
    expect(cb).toHaveAttribute("aria-checked", "mixed");
    expect(cb).toHaveAttribute("data-indeterminate", "true");
  });

  it("calls onChange when clicked", async () => {
    const user = userEvent.setup();
    const onChange = vi.fn();
    render(<Checkbox onChange={onChange} aria-label="select sample" />);
    await user.click(screen.getByRole("checkbox", { name: "select sample" }));
    expect(onChange).toHaveBeenCalledTimes(1);
  });

  it("does not call onChange when disabled", async () => {
    const user = userEvent.setup();
    const onChange = vi.fn();
    render(<Checkbox disabled onChange={onChange} aria-label="select sample" />);
    await user.click(screen.getByRole("checkbox", { name: "select sample" }));
    expect(onChange).not.toHaveBeenCalled();
  });
});
