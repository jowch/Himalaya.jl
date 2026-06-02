import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { Sparkline } from "../../src/print/plot/Sparkline";

const sampleTrace = {
  q: [0.05, 0.10, 0.15, 0.20, 0.25, 0.30],
  I: [500, 1200, 800, 400, 200, 100],
};

describe("Sparkline", () => {
  it("(a) renders [data-testid='sparkline'] container", () => {
    const { container } = render(<Sparkline trace={sampleTrace} />);
    expect(container.querySelector('[data-testid="sparkline"]')).not.toBeNull();
  });

  it("(b) renders exactly one [data-role='spark-line'] with a non-empty d attribute", () => {
    const { container } = render(<Sparkline trace={sampleTrace} />);
    const lines = container.querySelectorAll('[data-role="spark-line"]');
    expect(lines.length).toBe(1);
    const d = lines[0]!.getAttribute("d");
    expect(d).toBeTruthy();
    expect(d!.length).toBeGreaterThan(0);
  });

  it("(c) renders [data-role='spark-baseline']", () => {
    const { container } = render(<Sparkline trace={sampleTrace} />);
    expect(container.querySelector('[data-role="spark-baseline"]')).not.toBeNull();
  });

  it("(d) phase drives hue: Lamellar uses oklch, null uses ink-faint token, and they differ", () => {
    const { container: lamellarContainer } = render(
      <Sparkline trace={sampleTrace} phase="Lamellar" />,
    );
    const { container: nullContainer } = render(
      <Sparkline trace={sampleTrace} phase={null} />,
    );

    const lamellarStroke = lamellarContainer
      .querySelector('[data-role="spark-line"]')!
      .getAttribute("stroke");
    const nullStroke = nullContainer
      .querySelector('[data-role="spark-line"]')!
      .getAttribute("stroke");

    expect(lamellarStroke).toBeTruthy();
    expect(nullStroke).toBeTruthy();
    expect(lamellarStroke).not.toBe(nullStroke);
    expect(lamellarStroke!).toContain("oklch");
  });

  it("(e) degenerate trace (single point) renders baseline and no spark-line, no throw", () => {
    expect(() => {
      const { container } = render(
        <Sparkline trace={{ q: [0.1], I: [5] }} />,
      );
      expect(container.querySelector('[data-role="spark-baseline"]')).not.toBeNull();
      expect(container.querySelector('[data-role="spark-line"]')).toBeNull();
    }).not.toThrow();
  });
});
