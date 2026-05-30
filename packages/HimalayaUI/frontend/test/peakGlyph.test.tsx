import { describe, it, expect } from "vitest";
import { render } from "@testing-library/react";
import { PeakGlyph } from "../src/components/ui/peakMark";
import { peakGlyph } from "../src/components/ui/peakMark";

const PHASE = "oklch(0.570 0.150 58)";

function glyphFor(opts: Parameters<typeof peakGlyph>[0]) {
  const { container } = render(
    <svg>
      <PeakGlyph descriptor={peakGlyph(opts)} x={40} y={20} />
    </svg>,
  );
  return container;
}

describe("<PeakGlyph> — legend + overlay geometry", () => {
  it("auto renders a downward-pointing triangle (polygon, apex below the bar)", () => {
    const c = glyphFor({ source: "auto", color: PHASE });
    const poly = c.querySelector("polygon");
    expect(poly).not.toBeNull();
    expect(poly!.getAttribute("data-shape")).toBe("triangle-down");
    // apex (lowest y) is the peak point; the two top vertices share the higher edge.
    const pts = poly!
      .getAttribute("points")!
      .trim()
      .split(/\s+/)
      .map((p) => p.split(",").map(Number));
    const ys = pts.map((p) => p[1]!);
    // downward triangle: exactly one vertex has the max y (the apex points down)
    const maxY = Math.max(...ys);
    expect(ys.filter((y) => y === maxY).length).toBe(1);
  });

  it("manual renders a diamond, NOT magenta", () => {
    const c = glyphFor({ source: "manual", color: PHASE });
    const poly = c.querySelector("polygon");
    expect(poly).not.toBeNull();
    expect(poly!.getAttribute("data-shape")).toBe("diamond");
    expect(poly!.getAttribute("fill")).toBe(PHASE);
    expect(poly!.getAttribute("fill")).not.toMatch(/340/);
  });

  it("predicted-absent renders a hollow caret (polyline, no fill)", () => {
    const c = glyphFor({ source: "auto", color: PHASE, predictedAbsent: true });
    const caret = c.querySelector('[data-shape="caret"]');
    expect(caret).not.toBeNull();
    expect(caret!.getAttribute("fill")).toBe("none");
  });

  it("optimistic renders outline-only (fill none)", () => {
    const c = glyphFor({ source: "manual", color: PHASE, optimistic: true });
    const poly = c.querySelector("polygon");
    expect(poly!.getAttribute("fill")).toBe("none");
  });

  it("hot adds a ring element", () => {
    const c = glyphFor({ source: "auto", color: PHASE, hot: true });
    expect(c.querySelector('[data-role="hot-ring"]')).not.toBeNull();
  });

  it("excluded is hollow (fill none)", () => {
    const c = glyphFor({ source: "auto", color: PHASE, excluded: true });
    const poly = c.querySelector("polygon");
    expect(poly!.getAttribute("fill")).toBe("none");
  });
});
