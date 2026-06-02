import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { makeProjection } from "../../src/print/plot/projection";
import { PlotPeaks, type PlotPeak } from "../../src/print/plot/marks/PlotPeaks";

describe("PlotPeaks", () => {
  const proj = makeProjection({
    xDomain: [0.1, 0.3],
    yDomain: [1, 30],
    plotWidth: 300,
    plotHeight: 200,
    xType: "log",
    yType: "log",
  });

  it("renders one glyph per real peak and skips id < 0", () => {
    const peaks: PlotPeak[] = [
      { id: 1, q: 0.15, source: "auto", intensity: 10 },
      { id: -2, q: 0.2, source: "auto", intensity: 5 },
    ];
    const { container } = render(
      <svg>
        <PlotPeaks peaks={peaks} projection={proj} color="oklch(0.5 0.1 30)" />
      </svg>,
    );
    expect(
      container.querySelectorAll('[data-role="peak-glyph"]').length,
    ).toBe(1);
  });

  it("draws a q-link line only for hot peaks", () => {
    const { container, rerender } = render(
      <svg>
        <PlotPeaks
          peaks={[{ id: 1, q: 0.2, intensity: 20, source: "auto" }]}
          projection={proj}
          color="var(--color-accent)"
        />
      </svg>,
    );
    expect(container.querySelector('[data-role="peak-qline"]')).toBeNull();
    rerender(
      <svg>
        <PlotPeaks
          peaks={[{ id: 1, q: 0.2, intensity: 20, source: "auto", hot: true }]}
          projection={proj}
          color="var(--color-accent)"
        />
      </svg>,
    );
    expect(container.querySelector('[data-role="peak-qline"]')).toBeTruthy();
  });

  it("uses per-peak color when present, ignoring the layer color", () => {
    const { container } = render(
      <svg>
        <PlotPeaks
          peaks={[{ id: 1, q: 0.2, intensity: 20, source: "auto", color: "var(--color-success)" }]}
          projection={proj}
          color="var(--color-accent)"
        />
      </svg>,
    );
    // The auto peak is filled → polygon fill should be the per-peak color, not the layer color.
    expect(
      container.querySelector('[data-role="peak-glyph"] polygon')?.getAttribute('fill'),
    ).toBe("var(--color-success)");
  });

  it("falls back to layer color when peak has no per-peak color", () => {
    const { container } = render(
      <svg>
        <PlotPeaks
          peaks={[{ id: 1, q: 0.2, intensity: 20, source: "auto" }]}
          projection={proj}
          color="var(--color-accent)"
        />
      </svg>,
    );
    expect(
      container.querySelector('[data-role="peak-glyph"] polygon')?.getAttribute('fill'),
    ).toBe("var(--color-accent)");
  });
});
