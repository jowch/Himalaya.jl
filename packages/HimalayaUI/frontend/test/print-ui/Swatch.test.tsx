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

  it("exposes square shape by default on data-shape", () => {
    const { container } = render(<Swatch phase="Lamellar" />);
    expect(container.querySelector("[data-swatch]")!).toHaveAttribute("data-shape", "square");
  });

  it("exposes circle shape on data-shape when shape='circle'", () => {
    const { container } = render(<Swatch phase="Lamellar" shape="circle" />);
    expect(container.querySelector("[data-swatch]")!).toHaveAttribute("data-shape", "circle");
  });

  it("renders a coexistence gradient blending both phases", () => {
    const { getByTestId } = render(<div data-testid="w"><Swatch phase="Pn3m" coexistWith="Lamellar" shape="circle" /></div>);
    const sw = getByTestId("w").firstChild as HTMLElement;
    expect(sw.getAttribute("data-coexist")).toBe("Lamellar");
    // The DOM canonicalizes the OKLCH substrings on write, so compare against an
    // identically-built reference (round-trip), as the solid-color test does.
    const ref = document.createElement("span");
    ref.style.background = `linear-gradient(135deg, ${phaseColor("Pn3m")} 48%, ${phaseColor("Lamellar")} 52%)`;
    expect(sw.style.background).toBe(ref.style.background);
  });

  it("renders an empty (form-factor) swatch with no phase fill", () => {
    const { getByTestId } = render(<div data-testid="w"><Swatch phase="Pn3m" empty /></div>);
    const sw = getByTestId("w").firstChild as HTMLElement;
    expect(sw.getAttribute("data-empty")).toBe("true");
    expect(sw.style.background).not.toContain(phaseColor("Pn3m"));
  });

  it("size md is 11px", () => {
    const { getByTestId } = render(<div data-testid="w"><Swatch phase="Pn3m" size="md" /></div>);
    const sw = getByTestId("w").firstChild as HTMLElement;
    expect(sw.className).toContain("w-[11px]");
  });
});
