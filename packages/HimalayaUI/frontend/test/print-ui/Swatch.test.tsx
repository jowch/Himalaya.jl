import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { Swatch } from "../../src/print/ui/Swatch";
import { phaseColor } from "../../src/phases";

describe("<Swatch> (print)", () => {
  it("renders an element with [data-swatch]", () => {
    const { container } = render(<Swatch phase="Pn3m" />);
    expect(container.querySelector("[data-swatch]")).not.toBeNull();
  });

  it("reflects the phase on data-phase", () => {
    const { container } = render(<Swatch phase="Im3m" />);
    expect(container.querySelector("[data-swatch]")!).toHaveAttribute("data-phase", "Im3m");
  });

  it("sources the real phase color via inline style.background", () => {
    const { container } = render(<Swatch phase="Pn3m" />);
    const el = container.querySelector("[data-swatch]") as HTMLElement;
    // Compare via an identically-normalized reference: the DOM canonicalizes
    // the OKLCH string on write (e.g. trims trailing zeros), so we round-trip
    // phaseColor() through a sibling element's style and compare the read-backs.
    // This proves the component sources phaseColor(phase) exactly.
    const ref = document.createElement("span");
    ref.style.background = phaseColor("Pn3m");
    expect(el.style.background).toBe(ref.style.background);
  });

  it("is decorative (aria-hidden=true)", () => {
    const { container } = render(<Swatch phase="Ia3d" />);
    expect(container.querySelector("[data-swatch]")!).toHaveAttribute("aria-hidden", "true");
  });

  it("carries rounded-sm geometry by default (square)", () => {
    const { container } = render(<Swatch phase="Lamellar" />);
    expect((container.querySelector("[data-swatch]") as HTMLElement).className).toContain("rounded-sm");
  });

  it("carries rounded-full geometry when shape='circle'", () => {
    const { container } = render(<Swatch phase="Lamellar" shape="circle" />);
    expect((container.querySelector("[data-swatch]") as HTMLElement).className).toContain("rounded-full");
  });
});
