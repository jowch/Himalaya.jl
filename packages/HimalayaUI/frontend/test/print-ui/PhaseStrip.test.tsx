import { render, screen, within } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { PhaseStrip } from "../../src/print/ui/PhaseStrip";
import type { PhaseSegment } from "../../src/print/ui/PhaseStrip";

function segs(...segments: PhaseSegment[]): PhaseSegment[] {
  return segments;
}

describe("<PhaseStrip> single phase", () => {
  it("renders a solid background (no gradient) and no data-coexist-count", () => {
    render(<PhaseStrip segments={segs({ phase: "Pn3m" })} />);
    const cell = screen.getByTestId("ps-seg");
    const bg = cell.style.background;
    expect(bg).not.toContain("linear-gradient");
    expect(bg.length).toBeGreaterThan(0);
    expect(cell.getAttribute("data-coexist-count")).toBeNull();
    expect(cell.getAttribute("aria-label")).toBe("Pn3m");
    expect(cell.getAttribute("title")).toBe("Pn3m");
  });
});

describe("<PhaseStrip> 2-phase coexistence", () => {
  it("labels both phases, sets data-coexist-count=2, and paints an N-band gradient", () => {
    render(
      <PhaseStrip
        segments={segs({ phase: "Pn3m", coexistWith: ["Lamellar"] })}
      />,
    );
    const cell = screen.getByTestId("ps-seg");
    expect(cell.getAttribute("aria-label")).toBe("Pn3m + Lamellar (coexistence)");
    expect(cell.getAttribute("title")).toBe("Pn3m + Lamellar (coexistence)");
    expect(cell.getAttribute("data-coexist-count")).toBe("2");

    const bg = cell.style.background;
    expect(bg.startsWith("linear-gradient(100deg")).toBe(true);
    // Equal hard bands: the 1/2 boundary lands at 50.00%.
    expect(bg).toContain("50.00%");
    // Two distinct color stops (the dominant + the coexisting phase).
    const colorStops = bg.match(/oklch|rgb|hsl|#|var\(/g) ?? [];
    expect(colorStops.length).toBeGreaterThanOrEqual(2);
  });
});

describe("<PhaseStrip> 3-phase coexistence", () => {
  it("labels all three phases, sets data-coexist-count=3, and paints a gradient", () => {
    render(
      <PhaseStrip
        segments={segs({
          phase: "Pn3m",
          coexistWith: ["Lamellar", "Hexagonal"],
        })}
      />,
    );
    const cell = screen.getByTestId("ps-seg");
    expect(cell.getAttribute("aria-label")).toBe(
      "Pn3m + Lamellar + Hexagonal (coexistence)",
    );
    expect(cell.getAttribute("title")).toBe(
      "Pn3m + Lamellar + Hexagonal (coexistence)",
    );
    expect(cell.getAttribute("data-coexist-count")).toBe("3");
    expect(cell.style.background.startsWith("linear-gradient(100deg")).toBe(
      true,
    );
  });
});

describe("<PhaseStrip> assistive-tech exposure", () => {
  it("exposes every segment as role=img with the right accessible names", () => {
    render(
      <PhaseStrip
        segments={segs(
          { phase: "Pn3m" },
          { phase: null },
          { phase: null, state: "form_factor" },
          { phase: null, state: "null" },
          { phase: "Pn3m", coexistWith: ["Lamellar"] },
        )}
      />,
    );
    const cells = screen.getAllByRole("img");
    expect(cells).toHaveLength(5);
    expect(screen.getByRole("img", { name: "Pn3m" })).toBeInTheDocument();
    expect(screen.getByRole("img", { name: "Unindexed" })).toBeInTheDocument();
    expect(
      screen.getByRole("img", { name: "Form factor (no Bragg peaks)" }),
    ).toBeInTheDocument();
    expect(screen.getByRole("img", { name: "No phase" })).toBeInTheDocument();
    expect(
      screen.getByRole("img", { name: "Pn3m + Lamellar (coexistence)" }),
    ).toBeInTheDocument();
  });

  it("speaks the transition relation: arrow stays decorative, 'to' is accessible", () => {
    render(
      <PhaseStrip segments={segs({ phase: "Pn3m" }, { phase: "Lamellar" })} />,
    );
    const cap = screen.getByTestId("ps-cap");
    // The visual arrow is hidden from AT...
    const arrow = cap.querySelector('[aria-hidden="true"]');
    expect(arrow).not.toBeNull();
    expect(arrow).toHaveTextContent("→");
    // ...and the relation word "to" is exposed between the phase names.
    const to = within(cap).getByText("to");
    expect(to).not.toHaveAttribute("aria-hidden");
  });
});

describe("<PhaseStrip> coexistWith empty/null", () => {
  it("treats an empty coexistWith array as single-phase", () => {
    render(
      <PhaseStrip segments={segs({ phase: "Pn3m", coexistWith: [] })} />,
    );
    const cell = screen.getByTestId("ps-seg");
    expect(cell.style.background).not.toContain("linear-gradient");
    expect(cell.getAttribute("data-coexist-count")).toBeNull();
    expect(cell.getAttribute("aria-label")).toBe("Pn3m");
  });

  it("treats a null coexistWith as single-phase", () => {
    render(
      <PhaseStrip segments={segs({ phase: "Pn3m", coexistWith: null })} />,
    );
    const cell = screen.getByTestId("ps-seg");
    expect(cell.style.background).not.toContain("linear-gradient");
    expect(cell.getAttribute("data-coexist-count")).toBeNull();
  });
});
