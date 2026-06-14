import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { makeProjection } from "../../src/print/plot/projection";
import {
  PlotLabels,
  type PlotLabelsProps,
} from "../../src/print/plot/marks/PlotLabels";
import { type PlotPeak } from "../../src/print/plot/marks/PlotPeaks";

// Shared projection: linear x so we can reason about pixel positions exactly.
const proj = makeProjection({
  xDomain: [0, 300],
  yDomain: [1, 30],
  plotWidth: 300,
  plotHeight: 200,
  xType: "linear",
  yType: "log",
});
// With linear xDomain=[0,300] plotWidth=300, proj.x.to(q) === q (1:1 mapping).

const color = "oklch(0.5 0.1 30)";

function wrap(props: PlotLabelsProps) {
  return render(
    <svg>
      <PlotLabels {...props} />
    </svg>,
  );
}

describe("PlotLabels", () => {
  // ── 3a-1 ── renders one <text> per labelled peak --------------------------
  it("renders 2 peak-label texts for 2 labelled peaks", () => {
    const peaks: PlotPeak[] = [
      { id: 1, q: 50, source: "auto", intensity: 10, label: "α" },
      { id: 2, q: 200, source: "auto", intensity: 10, label: "β" },
    ];
    const { container } = wrap({ peaks, projection: proj, color });
    const texts = container.querySelectorAll('[data-role="peak-label"]');
    expect(texts.length).toBe(2);
    const contents = Array.from(texts).map((t) => t.textContent);
    expect(contents).toContain("α");
    expect(contents).toContain("β");
  });

  // ── 3a-2 ── empty label, undefined label, excluded (even if labelled) ----
  it("renders no label for empty-string label", () => {
    const peaks: PlotPeak[] = [
      { id: 1, q: 100, source: "auto", intensity: 10, label: "" },
    ];
    const { container } = wrap({ peaks, projection: proj, color });
    expect(
      container.querySelectorAll('[data-role="peak-label"]').length,
    ).toBe(0);
  });

  it("renders no label when label is undefined", () => {
    const peaks: PlotPeak[] = [
      { id: 1, q: 100, source: "auto", intensity: 10 },
    ];
    const { container } = wrap({ peaks, projection: proj, color });
    expect(
      container.querySelectorAll('[data-role="peak-label"]').length,
    ).toBe(0);
  });

  it("renders no label for excluded peaks even if labelled", () => {
    const peaks: PlotPeak[] = [
      { id: 1, q: 100, source: "auto", intensity: 10, label: "X", excluded: true },
    ];
    const { container } = wrap({ peaks, projection: proj, color });
    expect(
      container.querySelectorAll('[data-role="peak-label"]').length,
    ).toBe(0);
  });

  // ── 3a-3 ── id < 0 renders no label -----------------------------------------
  it("renders no label for peaks with id < 0", () => {
    const peaks: PlotPeak[] = [
      { id: -1, q: 100, source: "auto", intensity: 10, label: "neg" },
    ];
    const { container } = wrap({ peaks, projection: proj, color });
    expect(
      container.querySelectorAll('[data-role="peak-label"]').length,
    ).toBe(0);
  });

  // ── 3a-4 ── dodge: two close peaks get distinct x positions ---------------
  it("dodges two close peaks to distinct x positions (by their measured width)", () => {
    // With 1:1 mapping, q=100 → px=100, q=105 → px=105 (5px apart). Single-char
    // labels are narrow, so they push apart by their measured width + gap — far
    // less than the old fixed 30px (the over-aggressive dodge this replaces).
    const peaks: PlotPeak[] = [
      { id: 1, q: 100, source: "auto", intensity: 10, label: "A" },
      { id: 2, q: 105, source: "auto", intensity: 10, label: "B" },
    ];
    const { container } = wrap({ peaks, projection: proj, color });
    const texts = container.querySelectorAll('[data-role="peak-label"]');
    expect(texts.length).toBe(2);
    const xs = Array.from(texts).map((t) =>
      parseFloat(t.getAttribute("x") ?? "NaN"),
    );
    expect(xs[0]).not.toBeNaN();
    expect(xs[1]).not.toBeNaN();
    const gap = Math.abs(xs[1]! - xs[0]!);
    expect(gap).toBeGreaterThan(5); // pushed apart from their natural 5px
    expect(gap).toBeLessThan(30); // but far less than the old uniform width
  });

  // ── highlightPeakIds — label dimming ---------------------------------------
  describe("highlightPeakIds — label dimming", () => {
    const peaks3: PlotPeak[] = [
      { id: 1, q: 50, source: "auto", intensity: 10, label: "α" },
      { id: 2, q: 150, source: "auto", intensity: 15, label: "β" },
      { id: 3, q: 250, source: "auto", intensity: 20, label: "γ" },
    ];

    it("non-highlighted labels have data-dimmed='true' and fill=ink-faint", () => {
      const { container } = wrap({
        peaks: peaks3,
        projection: proj,
        color,
        highlightPeakIds: new Set([1]),
      });
      const allTexts = container.querySelectorAll('[data-role="peak-label"]');
      expect(allTexts.length).toBe(3);

      // Find text elements by content
      const t1 = Array.from(allTexts).find((t) => t.textContent === "α");
      const t2 = Array.from(allTexts).find((t) => t.textContent === "β");
      const t3 = Array.from(allTexts).find((t) => t.textContent === "γ");
      expect(t1).toBeTruthy();
      expect(t2).toBeTruthy();
      expect(t3).toBeTruthy();

      // Highlighted: NOT dimmed
      expect(t1!.getAttribute("data-dimmed")).toBeNull();

      // Non-highlighted: dimmed
      expect(t2!.getAttribute("data-dimmed")).toBe("true");
      expect(t3!.getAttribute("data-dimmed")).toBe("true");

      // Dim colour lands on CSS `color` (ink-faint) with fill=currentColor so
      // the dim/release eases through the global 120ms ease-out `color`
      // transition instead of snapping.
      expect((t2 as HTMLElement).style.color).toBe("var(--color-ink-faint)");
      expect((t3 as HTMLElement).style.color).toBe("var(--color-ink-faint)");
      expect((t2 as HTMLElement).style.fill).toBe("currentColor");
      expect((t3 as HTMLElement).style.fill).toBe("currentColor");
    });

    it("no label has data-dimmed when highlightPeakIds is not provided", () => {
      const { container } = wrap({ peaks: peaks3, projection: proj, color });
      const dimmed = container.querySelectorAll(
        '[data-role="peak-label"][data-dimmed]',
      );
      expect(dimmed.length).toBe(0);
    });

    it("no label has data-dimmed when highlightPeakIds is an empty set", () => {
      const { container } = wrap({
        peaks: peaks3,
        projection: proj,
        color,
        highlightPeakIds: new Set(),
      });
      const dimmed = container.querySelectorAll(
        '[data-role="peak-label"][data-dimmed]',
      );
      expect(dimmed.length).toBe(0);
    });
  });

  // dodgeX's pure unit tests live in print-plot/labelDodge.test.ts (the
  // width-aware signature).
});
