import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Button } from "../../src/print/ui/Button";

describe("<Button> armed state", () => {
  it("defaults to not-armed and sets no aria-pressed", () => {
    render(<Button>+ Peak</Button>);
    const b = screen.getByRole("button", { name: "+ Peak" });
    expect(b.getAttribute("data-armed")).toBe(null);
    expect(b.getAttribute("aria-pressed")).toBe(null);
  });
  it("reflects armed on data-armed and aria-pressed", () => {
    render(<Button armed>+ Peak</Button>);
    const b = screen.getByRole("button", { name: "+ Peak" });
    expect(b.getAttribute("data-armed")).toBe("true");
    expect(b.getAttribute("aria-pressed")).toBe("true");
  });
  it("keeps the variant data-attr when armed", () => {
    render(<Button variant="ghost" armed>x</Button>);
    expect(screen.getByRole("button").getAttribute("data-variant")).toBe("ghost");
  });
  it("still fires onClick when armed", () => {
    const onClick = vi.fn();
    render(<Button armed onClick={onClick}>x</Button>);
    fireEvent.click(screen.getByRole("button"));
    expect(onClick).toHaveBeenCalledTimes(1);
  });
});

describe("<Button> danger variant", () => {
  it("carries data-variant=danger and renders its label", () => {
    render(<Button variant="danger">Delete sample</Button>);
    const b = screen.getByRole("button", { name: "Delete sample" });
    expect(b.getAttribute("data-variant")).toBe("danger");
  });
  it("fires onClick", () => {
    const onClick = vi.fn();
    render(<Button variant="danger" onClick={onClick}>Delete sample</Button>);
    fireEvent.click(screen.getByRole("button", { name: "Delete sample" }));
    expect(onClick).toHaveBeenCalledTimes(1);
  });
});
