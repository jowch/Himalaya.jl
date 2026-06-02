import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { makeProjection } from "../../src/print/plot/projection";
import { TraceLine } from "../../src/print/plot/marks/TraceLine";
import type { Trace } from "../../src/api";

describe("TraceLine", () => {
  const proj = makeProjection({
    xDomain: [0.1, 0.3],
    yDomain: [1, 30],
    plotWidth: 300,
    plotHeight: 200,
    xType: "log",
    yType: "log",
  });

  it("draws neutral ink line (var(--color-ink-soft), strokeWidth 1.8) and faint-ink sigma band (opacity 0.12)", () => {
    const { container } = render(
      <svg>
        <TraceLine
          trace={{ q: [0.1, 0.2, 0.3], I: [10, 20, 15], sigma: [1, 1, 1] }}
          projection={proj}
        />
      </svg>,
    );
    const paths = container.querySelectorAll('[data-role="trace-line"] path');
    const linePath = paths[paths.length - 1]!; // line drawn last (on top of band)
    expect(linePath.getAttribute("stroke")).toBe("var(--color-ink-soft)");
    expect(linePath.getAttribute("stroke-width")).toBe("1.8");
    const band = container.querySelector(
      '[data-role="trace-line"] path[opacity="0.12"]',
    );
    expect(band?.getAttribute("fill")).toBe("var(--color-ink-soft)");
  });

  it("draws a band path and a line path (line starts with M)", () => {
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 15], sigma: [1, 2, 1.5] };
    const { container } = render(
      <svg>
        <TraceLine trace={trace} projection={proj} />
      </svg>,
    );
    const paths = container.querySelectorAll('[data-role="trace-line"] path');
    expect(paths.length).toBe(2); // band + line
    expect(paths[1]!.getAttribute("d")).toMatch(/^M/);
  });

  it("omits the band when band={false}", () => {
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 15], sigma: [1, 2, 1.5] };
    const { container } = render(
      <svg>
        <TraceLine trace={trace} projection={proj} band={false} />
      </svg>,
    );
    const paths = container.querySelectorAll('[data-role="trace-line"] path');
    expect(paths.length).toBe(1);
  });
});

const colorTrace: Trace = { q: [0.02, 0.05, 0.1], I: [100, 50, 10], sigma: [1, 1, 1] };
const colorProjection = makeProjection({
  xDomain: [0.02, 0.1], yDomain: [10, 100],
  plotWidth: 200, plotHeight: 100, xType: "log", yType: "log",
});

// The line is the second <path> when a band is present; with band={false} it is
// the only path. Use band={false} so the line path is unambiguous.
function linePath(container: HTMLElement): SVGPathElement {
  const paths = container.querySelectorAll('[data-role="trace-line"] path');
  return paths[paths.length - 1] as SVGPathElement;
}

describe("TraceLine color", () => {
  it("defaults the line stroke to ink-soft", () => {
    const { container } = render(
      <svg><TraceLine trace={colorTrace} projection={colorProjection} band={false} /></svg>,
    );
    expect(linePath(container).getAttribute("stroke")).toBe("var(--color-ink-soft)");
  });

  it("uses the provided color for the line stroke", () => {
    const { container } = render(
      <svg><TraceLine trace={colorTrace} projection={colorProjection} band={false} color="oklch(0.570 0.150 58)" /></svg>,
    );
    expect(linePath(container).getAttribute("stroke")).toBe("oklch(0.570 0.150 58)");
  });
});
