import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
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
    const qline = container.querySelector('[data-role="peak-qline"]')!;
    expect(qline).toBeTruthy();
    // Extends DOWN from the marker to the baseline (y1 = marker < y2 = baseline),
    // and is slightly thicker than the resting hairline.
    expect(Number(qline.getAttribute("y1"))).toBeLessThan(
      Number(qline.getAttribute("y2")),
    );
    expect(qline.getAttribute("stroke-width")).toBe("1.5");
  });

  it("a hot peak draws a concentric same-shape halo, not a terracotta recolour", () => {
    const { container } = render(
      <svg>
        <PlotPeaks
          peaks={[{ id: 1, q: 0.2, intensity: 20, source: "auto", hot: true }]}
          projection={proj}
          color="var(--color-success)"
        />
      </svg>,
    );
    // A separate, larger, same-shape outline element is drawn.
    expect(container.querySelector('[data-role="peak-halo"]')).toBeTruthy();
    // The mark itself keeps its phase colour — it is NOT recoloured terracotta.
    const mark = container.querySelector('[data-role="peak-glyph"] [data-shape]')!;
    expect(mark.getAttribute("stroke")).not.toBe("var(--color-accent)");
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

  it("peak <g> has tabIndex/role/aria-label when onPeakFocus is provided", () => {
    const { container } = render(
      <svg>
        <PlotPeaks
          peaks={[{ id: 1, q: 0.2, intensity: 20, source: "auto" }]}
          projection={proj}
          color="var(--color-accent)"
          onPeakFocus={vi.fn()}
        />
      </svg>,
    );
    const peakG = container.querySelector('[data-role="plot-peaks"] > g');
    expect(peakG).toBeTruthy();
    expect(peakG!.getAttribute("tabindex")).toBe("0");
    expect(peakG!.getAttribute("role")).toBe("button");
    expect(peakG!.getAttribute("aria-label")).toMatch(/peak at q/);
  });

  it("calls onPeakFocus(id) on focus and onPeakFocus(null) on blur", () => {
    const spy = vi.fn();
    const { container } = render(
      <svg>
        <PlotPeaks
          peaks={[{ id: 42, q: 0.2, intensity: 20, source: "auto" }]}
          projection={proj}
          color="var(--color-accent)"
          onPeakFocus={spy}
        />
      </svg>,
    );
    const peakG = container.querySelector('[data-role="plot-peaks"] > g')!;
    fireEvent.focus(peakG);
    expect(spy).toHaveBeenCalledWith(42);
    fireEvent.blur(peakG);
    expect(spy).toHaveBeenCalledWith(null);
  });

  it("does NOT add tabIndex/role when onPeakFocus is not provided", () => {
    const { container } = render(
      <svg>
        <PlotPeaks
          peaks={[{ id: 1, q: 0.2, intensity: 20, source: "auto" }]}
          projection={proj}
          color="var(--color-accent)"
        />
      </svg>,
    );
    const peakG = container.querySelector('[data-role="plot-peaks"] > g');
    // No interaction wired — tabindex should not be present
    expect(peakG?.getAttribute("tabindex")).toBeNull();
    expect(peakG?.getAttribute("role")).toBeNull();
  });

  describe("highlightPeakIds — index/phase highlight", () => {
    const peaks3: PlotPeak[] = [
      { id: 1, q: 0.13, source: "auto", intensity: 10 },
      { id: 2, q: 0.18, source: "auto", intensity: 15 },
      { id: 3, q: 0.25, source: "auto", intensity: 20 },
    ];

    it("non-highlighted peaks have data-dimmed='true' and glyph fill=ink-faint", () => {
      const { container } = render(
        <svg>
          <PlotPeaks
            peaks={peaks3}
            projection={proj}
            color="oklch(0.5 0.1 30)"
            highlightPeakIds={new Set([1])}
          />
        </svg>,
      );
      const allPeakGs = container.querySelectorAll('[data-role="plot-peaks"] > g');
      // Find the <g> for peak id=1 (highlighted) — must NOT be dimmed
      const g1 = Array.from(allPeakGs).find(
        (g) => g.querySelector('[data-peak-id="1"]'),
      );
      expect(g1).toBeTruthy();
      expect(g1!.getAttribute("data-dimmed")).toBeNull();

      // Peaks 2 and 3 must be dimmed
      const g2 = Array.from(allPeakGs).find(
        (g) => g.querySelector('[data-peak-id="2"]'),
      );
      const g3 = Array.from(allPeakGs).find(
        (g) => g.querySelector('[data-peak-id="3"]'),
      );
      expect(g2).toBeTruthy();
      expect(g3).toBeTruthy();
      expect(g2!.getAttribute("data-dimmed")).toBe("true");
      expect(g3!.getAttribute("data-dimmed")).toBe("true");

      // Dimmed glyph fill must be ink-faint
      expect(
        g2!.querySelector('[data-role="peak-glyph"] polygon')?.getAttribute("fill"),
      ).toBe("var(--color-ink-faint)");
      expect(
        g3!.querySelector('[data-role="peak-glyph"] polygon')?.getAttribute("fill"),
      ).toBe("var(--color-ink-faint)");
    });

    it("no <g> has data-dimmed when highlightPeakIds is not provided", () => {
      const { container } = render(
        <svg>
          <PlotPeaks peaks={peaks3} projection={proj} color="oklch(0.5 0.1 30)" />
        </svg>,
      );
      const dimmed = container.querySelectorAll(
        '[data-role="plot-peaks"] > g[data-dimmed]',
      );
      expect(dimmed.length).toBe(0);
    });

    it("no <g> has data-dimmed when highlightPeakIds is an empty set", () => {
      const { container } = render(
        <svg>
          <PlotPeaks
            peaks={peaks3}
            projection={proj}
            color="oklch(0.5 0.1 30)"
            highlightPeakIds={new Set()}
          />
        </svg>,
      );
      const dimmed = container.querySelectorAll(
        '[data-role="plot-peaks"] > g[data-dimmed]',
      );
      expect(dimmed.length).toBe(0);
    });

    it("hot peak is exempt from dimming even if not in highlightPeakIds", () => {
      const peaksWithHot: PlotPeak[] = [
        { id: 1, q: 0.13, source: "auto", intensity: 10 },
        { id: 2, q: 0.18, source: "auto", intensity: 15, hot: true },
        { id: 3, q: 0.25, source: "auto", intensity: 20 },
      ];
      const { container } = render(
        <svg>
          <PlotPeaks
            peaks={peaksWithHot}
            projection={proj}
            color="oklch(0.5 0.1 30)"
            highlightPeakIds={new Set([1])}
          />
        </svg>,
      );
      const allPeakGs = container.querySelectorAll('[data-role="plot-peaks"] > g');
      // Peak 2 is hot and NOT in highlightPeakIds → must NOT be dimmed
      const g2 = Array.from(allPeakGs).find(
        (g) => g.querySelector('[data-peak-id="2"]'),
      );
      expect(g2).toBeTruthy();
      expect(g2!.getAttribute("data-dimmed")).toBeNull();
    });
  });
});
