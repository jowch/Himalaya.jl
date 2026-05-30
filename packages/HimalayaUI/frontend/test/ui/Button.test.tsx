import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Button } from "../../src/components/ui/Button";

describe("<Button>", () => {
  it("renders a button with its children and defaults to the ghost variant", () => {
    render(<Button>Save</Button>);
    const btn = screen.getByRole("button", { name: "Save" });
    expect(btn.getAttribute("data-variant")).toBe("ghost");
  });

  it("reflects each variant on data-variant", () => {
    const { rerender } = render(<Button variant="solid">x</Button>);
    expect(screen.getByRole("button").getAttribute("data-variant")).toBe("solid");
    rerender(<Button variant="accent">x</Button>);
    expect(screen.getByRole("button").getAttribute("data-variant")).toBe("accent");
    rerender(<Button variant="danger">x</Button>);
    expect(screen.getByRole("button").getAttribute("data-variant")).toBe("danger");
  });

  it("forwards disabled and fires onClick", () => {
    const onClick = vi.fn();
    const { rerender } = render(<Button onClick={onClick}>go</Button>);
    fireEvent.click(screen.getByRole("button"));
    expect(onClick).toHaveBeenCalledTimes(1);
    rerender(<Button disabled onClick={onClick}>go</Button>);
    expect(screen.getByRole("button")).toBeDisabled();
  });
});
