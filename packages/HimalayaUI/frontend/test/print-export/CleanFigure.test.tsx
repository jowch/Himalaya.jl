import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { CleanFigure } from "../../src/print/export/CleanFigure";
import { TRANSITION } from "../../src/print/waterfall/waterfall.fixtures";

describe("CleanFigure", () => {
  it("renders the card, title, footer, and one trace path per row", () => {
    const { container } = render(
      <CleanFigure
        rows={TRANSITION}
        title="LL37 Titration: Pn3m → Lamellar"
        footer="Lipid 1-2 · SSRL Apr 2026"
      />,
    );

    expect(screen.getByTestId("clean-figure")).toBeInTheDocument();
    expect(
      screen.getByText("LL37 Titration: Pn3m → Lamellar"),
    ).toBeInTheDocument();
    expect(screen.getByText(/Lipid 1-2/)).toBeInTheDocument();
    expect(
      screen.getByRole("img", { name: /LL37 Titration/ }),
    ).toBeInTheDocument();

    const paths = container.querySelectorAll("[data-phase]");
    expect(paths).toHaveLength(TRANSITION.length);
  });

  it("draws literal hex strokes per phase and renders both axis titles", () => {
    const { container } = render(
      <CleanFigure
        rows={TRANSITION}
        title="LL37 Titration: Pn3m → Lamellar"
        footer="Lipid 1-2 · SSRL Apr 2026"
      />,
    );

    // TRANSITION contains Im3m (HEX-mapped) and Ia3d (not in HEX → "#777"
    // fallback) — two distinct phases with distinct literal hex strokes.
    // Direct string compare: these are plain hex on the path, NOT OKLCH.
    const im3m = container.querySelector('[data-phase="Im3m"]');
    expect(im3m).not.toBeNull();
    expect(im3m!.getAttribute("stroke")).toBe("#4a4ba8");

    const ia3d = container.querySelector('[data-phase="Ia3d"]');
    expect(ia3d).not.toBeNull();
    expect(ia3d!.getAttribute("stroke")).toBe("#777");

    expect(screen.getByText("q (Å⁻¹)")).toBeInTheDocument();
    expect(screen.getByText("Intensity (a.u.)")).toBeInTheDocument();
  });
});
