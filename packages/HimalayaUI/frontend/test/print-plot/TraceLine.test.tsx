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

  it("draws a band path and a line path (line starts with M)", () => {
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 15], sigma: [1, 2, 1.5] };
    const { container } = render(
      <svg>
        <TraceLine trace={trace} projection={proj} color="oklch(0.5 0.1 30)" />
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
        <TraceLine
          trace={trace}
          projection={proj}
          color="oklch(0.5 0.1 30)"
          band={false}
        />
      </svg>,
    );
    const paths = container.querySelectorAll('[data-role="trace-line"] path');
    expect(paths.length).toBe(1);
  });
});
