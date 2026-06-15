import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { PhaseLabel } from "../../src/print/ui/PhaseLabel";
import { phaseColor } from "../../src/phases";

describe("<PhaseLabel> (print)", () => {
  it("renders its children text content", () => {
    const { getByText } = render(<PhaseLabel phase="Pn3m">Pn3m call</PhaseLabel>);
    expect(getByText("Pn3m call")).not.toBeNull();
  });

  it("marks the element with [data-phase-label]", () => {
    const { container } = render(<PhaseLabel phase="Pn3m">x</PhaseLabel>);
    expect(container.querySelector("[data-phase-label]")).not.toBeNull();
  });

  it("reflects the phase on data-phase", () => {
    const { container } = render(<PhaseLabel phase="Im3m">x</PhaseLabel>);
    expect(container.querySelector("[data-phase-label]")!).toHaveAttribute("data-phase", "Im3m");
  });

  it("sources the real phase color via inline style.color", () => {
    const { container } = render(<PhaseLabel phase="Pn3m">x</PhaseLabel>);
    const el = container.querySelector("[data-phase-label]") as HTMLElement;
    // The DOM canonicalizes the OKLCH string on write (e.g. trims trailing
    // zeros), so we round-trip phaseColor() through a sibling element's style
    // and compare the read-backs. This proves the component sources
    // phaseColor(phase) exactly.
    const ref = document.createElement("span");
    ref.style.color = phaseColor("Pn3m");
    expect(el.style.color).toBe(ref.style.color);
  });

  it("renders a span by default", () => {
    const { container } = render(<PhaseLabel phase="Ia3d">x</PhaseLabel>);
    expect((container.querySelector("[data-phase-label]") as HTMLElement).tagName).toBe("SPAN");
  });

  it("renders a div when as='div'", () => {
    const { container } = render(
      <PhaseLabel phase="Lamellar" as="div">
        x
      </PhaseLabel>,
    );
    expect((container.querySelector("[data-phase-label]") as HTMLElement).tagName).toBe("DIV");
  });

  it("forwards a placement className", () => {
    const { container } = render(
      <PhaseLabel phase="Pn3m" className="text-display">
        x
      </PhaseLabel>,
    );
    expect((container.querySelector("[data-phase-label]") as HTMLElement).className).toContain("text-display");
  });
});
