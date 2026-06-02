import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { CustomPreview } from "../../src/print/comb/CustomPreview";
import { PN3M } from "../../src/print/comb/comb.fixtures";
import type { CombSeries } from "../../src/print/comb/combModel";
import { phaseColor } from "../../src/phases";

/** Read back a color through JSDOM's canonicalizing style write, so phase-color
 *  assertions compare apples to apples (never the literal phaseColor string). */
function canonicalStroke(color: string): string {
  const ref = document.createElement("div");
  ref.style.stroke = color;
  return ref.style.stroke;
}

describe("CustomPreview", () => {
  it("renders an accessible img labelled for custom index, tagged with the phase", () => {
    const { getByRole } = render(
      <CustomPreview series={PN3M} observed={[0.0712, 0.0872, 0.156]} />,
    );
    const svg = getByRole("img", { name: /custom index reflections/i });
    expect(svg).toBeTruthy();
    expect(svg.getAttribute("data-phase")).toBe("Pn3m");
  });

  it("draws exactly one observed marker per observed q", () => {
    const { container } = render(
      <CustomPreview series={PN3M} observed={[0.0712, 0.0872, 0.156]} />,
    );
    expect(container.querySelectorAll("[data-observed-q]").length).toBe(3);
  });

  it("draws one tooth per predicted reflection, with at least one observed tooth", () => {
    const { container } = render(
      <CustomPreview series={PN3M} observed={[0.0712, 0.0872, 0.156]} />,
    );
    expect(container.querySelectorAll("[data-tooth-q]").length).toBe(PN3M.teeth.length);
    expect(
      container.querySelectorAll('[data-tooth-observed="true"]').length,
    ).toBeGreaterThanOrEqual(1);
  });

  it("draws a predicted-absent tooth as a dashed caret", () => {
    // Clone PN3M and force one tooth absent so the dashed-caret branch is exercised
    // regardless of fixture contents.
    const variant: CombSeries = {
      ...PN3M,
      teeth: PN3M.teeth.map((t, i) => (i === 0 ? { ...t, observed: false } : t)),
    };
    const { container } = render(<CustomPreview series={variant} observed={[0.0712]} />);
    const absent = container.querySelector('[data-tooth-observed="false"]');
    expect(absent).toBeTruthy();
    expect(absent!.getAttribute("stroke-dasharray")).toBeTruthy();
  });

  it("colors observed teeth with the phase color (read-back round-trip)", () => {
    const { container } = render(
      <CustomPreview series={PN3M} observed={[0.0712, 0.0872, 0.156]} />,
    );
    const solid = container.querySelector(
      '[data-tooth-observed="true"]',
    ) as SVGElement;
    expect(solid).toBeTruthy();
    expect(solid.style.stroke).toBe(canonicalStroke(phaseColor("Pn3m")));
  });
});
