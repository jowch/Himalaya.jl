import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { makeProjection } from "../../src/print/plot/projection";
import {
  PlotLabels,
  dodgeX,
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
  it("dodges two close peaks to distinct x positions", () => {
    // With 1:1 mapping, q=100 → px=100, q=105 → px=105 (5px apart).
    // labelWidthPx=30 → they must be pushed apart.
    const peaks: PlotPeak[] = [
      { id: 1, q: 100, source: "auto", intensity: 10, label: "A" },
      { id: 2, q: 105, source: "auto", intensity: 10, label: "B" },
    ];
    const { container } = wrap({
      peaks,
      projection: proj,
      color,
      labelWidthPx: 30,
    });
    const texts = container.querySelectorAll('[data-role="peak-label"]');
    expect(texts.length).toBe(2);
    const xs = Array.from(texts).map((t) =>
      parseFloat(t.getAttribute("x") ?? "NaN"),
    );
    expect(xs[0]).not.toBeNaN();
    expect(xs[1]).not.toBeNaN();
    // The two x-positions must differ by at least labelWidthPx.
    expect(Math.abs(xs[1]! - xs[0]!)).toBeGreaterThanOrEqual(29.9);
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

      // Dimmed fill must be ink-faint
      expect((t2 as HTMLElement).style.fill).toBe("var(--color-ink-faint)");
      expect((t3 as HTMLElement).style.fill).toBe("var(--color-ink-faint)");
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

  // ── dodgeX unit tests -------------------------------------------------------
  describe("dodgeX (pure helper)", () => {
    it("returns identity when positions are far apart", () => {
      const xs = [0, 100, 200];
      const result = dodgeX(xs, 30);
      expect(result).toEqual(xs);
    });

    it("pushes apart and recenters two overlapping positions", () => {
      // Two positions at px 100 and 105, labelWidthPx=30
      const xs = [100, 105];
      const result = dodgeX(xs, 30);
      expect(result.length).toBe(2);
      // After dodge the gap must be >= 30
      expect(result[1]! - result[0]!).toBeGreaterThanOrEqual(30 - 1e-9);
      // Recentered: mean of result should equal mean of input (102.5)
      const mean = (result[0]! + result[1]!) / 2;
      expect(mean).toBeCloseTo(102.5, 5);
    });

    it("handles a single element (no-op)", () => {
      expect(dodgeX([50], 30)).toEqual([50]);
    });

    it("handles empty array", () => {
      expect(dodgeX([], 30)).toEqual([]);
    });
  });
});
