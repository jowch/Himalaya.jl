import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { DetectorRings } from "../../src/print/detector/DetectorRings";
import type { RingPlacement } from "../../src/print/detector/detectorGeometry";

const rings: RingPlacement[] = [
  { q: 0.10, r: 0.12, color: "var(--color-success)" },
  { q: 0.20, r: 0.22, color: "var(--color-success)" },
  { q: 0.30, r: 0.34, color: "var(--color-success)", ghost: true },
  { q: 0.40, r: 0.46 }, // leftover, neutral
];
const beamCenter = { x: 0.5, y: 0.5 };

describe("DetectorRings", () => {
  it("renders one ring group per placement, each with casing + sharp + hit", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={rings} />);
    expect(container.querySelectorAll('[data-role="det-ring"]').length).toBe(4);
    const g = container.querySelector('[data-role="det-ring"]')!;
    // Dark casing under the coloured stroke so rings read against any background.
    expect(g.querySelector('[data-role="ring-casing"]')).toBeTruthy();
    expect(g.querySelector('[data-role="ring-sharp"]')).toBeTruthy();
    expect(g.querySelector('[data-role="ring-hit"]')).toBeTruthy();
  });

  it("the casing sits under (before) the coloured stroke in paint order", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={[{ q: 0.1, r: 0.3, color: "var(--color-success)" }]} />);
    const g = container.querySelector('[data-role="det-ring"]')!;
    const roles = Array.from(g.querySelectorAll("ellipse")).map((e) => e.getAttribute("data-role"));
    expect(roles.indexOf("ring-casing")).toBeLessThan(roles.indexOf("ring-sharp"));
    // The casing is wider than the stroke (it outlines it).
    const casing = g.querySelector('[data-role="ring-casing"]')!;
    const sharp = g.querySelector('[data-role="ring-sharp"]')!;
    expect(Number(casing.getAttribute("stroke-width"))).toBeGreaterThan(Number(sharp.getAttribute("stroke-width")));
  });

  it("centers every ring at the beam center (not the viewBox middle)", () => {
    const { container } = render(
      <DetectorRings beamCenter={{ x: 0.8, y: 0.1 }} rings={[{ q: 0.1, r: 0.3 }]} />,
    );
    const sharp = container.querySelector('[data-role="ring-sharp"]')!;
    expect(Number(sharp.getAttribute("cx"))).toBeCloseTo(0.8, 5);
    expect(Number(sharp.getAttribute("cy"))).toBeCloseTo(0.1, 5);
    // Default aspect 1 → rx === ry === r (a true circle on a square frame).
    expect(Number(sharp.getAttribute("rx"))).toBeCloseTo(0.3, 5);
    expect(Number(sharp.getAttribute("ry"))).toBeCloseTo(0.3, 5);
  });

  it("aspect-corrects ry so a ring renders as a true circle on a non-square frame", () => {
    // Portrait image (981×1043): preserveAspectRatio=none would squash a circle
    // in x; we pre-scale ry by width/height so the rendered pixel radii match.
    const aspect = 981 / 1043; // ≈ 0.9405
    const { container } = render(
      <DetectorRings beamCenter={{ x: 0.5, y: 0.5 }} rings={[{ q: 0.1, r: 0.3 }]} imageAspect={aspect} />,
    );
    const sharp = container.querySelector('[data-role="ring-sharp"]')!;
    expect(Number(sharp.getAttribute("rx"))).toBeCloseTo(0.3, 5);
    expect(Number(sharp.getAttribute("ry"))).toBeCloseTo(0.3 * aspect, 5);
    expect(Number(sharp.getAttribute("ry"))).toBeLessThan(Number(sharp.getAttribute("rx")));
  });

  it("draws a ghost ring dotted + hollow (zero-dash + round cap = dots, not dashes)", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={rings} />);
    const ghost = container.querySelector('[data-ghost="true"] [data-role="ring-sharp"]')!;
    expect(ghost.getAttribute("stroke-dasharray")).toBe("0 7");
    expect(ghost.getAttribute("stroke-linecap")).toBe("round");
    expect(ghost.getAttribute("fill")).toBe("none");
  });

  it("a solid (non-ghost) ring has no dash and no round cap", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={[{ q: 0.1, r: 0.2 }]} />);
    const sharp = container.querySelector('[data-role="ring-sharp"]')!;
    expect(sharp.getAttribute("stroke-dasharray")).toBeNull();
    expect(sharp.getAttribute("stroke-linecap")).toBeNull();
  });

  it("a leftover ring (no color) strokes neutral ink-faint", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={[{ q: 0.4, r: 0.4 }]} />);
    expect(container.querySelector('[data-role="ring-sharp"]')!.getAttribute("stroke")).toBe("var(--color-ink-faint)");
  });

  it("the ring matching hoveredQ goes hot via width/opacity, keeping its own color (no terracotta recolor)", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={rings} hoveredQ={0.20} />);
    const hot = container.querySelector('[data-hot="true"]')!;
    expect(hot.getAttribute("data-ring-q")).toBe("0.2");
    const sharp = hot.querySelector('[data-role="ring-sharp"]')!;
    // Keeps the ring's own colour — NOT recoloured to the terracotta accent.
    expect(sharp.getAttribute("stroke")).toBe("var(--color-success)");
    expect(sharp.getAttribute("stroke")).not.toBe("var(--color-accent)");
    // Emphasis is a wider, more opaque stroke than a resting ring.
    expect(Number(sharp.getAttribute("stroke-width"))).toBeGreaterThan(1.0);
    expect(Number(sharp.getAttribute("opacity"))).toBeGreaterThan(0.7);
  });

  it("no ring is hot when hoveredQ is undefined", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={rings} />);
    expect(container.querySelector('[data-hot="true"]')).toBeNull();
  });

  it("calls onHoverQ(q) on hit enter and onHoverQ(undefined) on leave", () => {
    const spy = vi.fn();
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={rings} onHoverQ={spy} />);
    const hit = container.querySelectorAll('[data-role="ring-hit"]')[1]; // q=0.20
    fireEvent.mouseEnter(hit); expect(spy).toHaveBeenCalledWith(0.20);
    fireEvent.mouseLeave(hit); expect(spy).toHaveBeenCalledWith(undefined);
  });

  it("hit ring is inert (pointer-events none) when onHoverQ is absent", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={rings} />);
    const hit = container.querySelector('[data-role="ring-hit"]') as SVGElement;
    expect(hit.style.pointerEvents).toBe("none");
  });
});
