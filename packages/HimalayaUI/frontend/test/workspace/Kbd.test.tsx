import { describe, it, expect } from "vitest";
import { render } from "@testing-library/react";
import { Kbd } from "../../src/workspace/ui/Kbd";

describe("Kbd", () => {
  it("renders a key cap on its own plate", () => {
    const { container } = render(<Kbd>⌘K</Kbd>);
    const k = container.querySelector("kbd")!;
    expect(k.textContent).toBe("⌘K");
    expect(k.className).toContain("bg-plate");
    expect(k.className).not.toContain("font-mono");
  });

  it("takes its colour from a filled button with onFill", () => {
    const { container } = render(<Kbd onFill>↵</Kbd>);
    const k = container.querySelector("kbd")!;
    expect(k.className).toContain("bg-current/15");
    expect(k.className).toContain("text-current");
    expect(k.className).not.toContain("bg-plate");
  });
});
