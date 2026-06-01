import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { makeProjection } from "../../src/print/plot/projection";
import { TraceLine } from "../../src/print/plot/marks/TraceLine";

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
