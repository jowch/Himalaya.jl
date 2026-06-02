import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { TracePlot, type TraceModel } from "../../src/print/plot/TracePlot";

const model: TraceModel = {
  trace: { q: [0.1, 0.2, 0.3], I: [10, 20, 15], sigma: [1, 1, 1] },
  peaks: [{ id: 1, q: 0.2, source: "auto", intensity: 20 }],
  phase: "Ia3d",
};

// Helper: the single peak <g> is the first child of [data-role="plot-peaks"]
function getPeakG(container: HTMLElement): Element | null {
  return container.querySelector('[data-role="plot-peaks"] > g');
}

describe("TracePlot", () => {
  it("renders axes, a trace line, and peak glyphs for a single trace", () => {
    const { container } = render(
      <TracePlot traces={[model]} width={500} height={300} data-testid="tp" />,
    );
    expect(container.querySelector('svg[data-testid="tp"]')).toBeTruthy();
    expect(
      container.querySelector('[data-role="trace-line"] path'),
    ).toBeTruthy();
    expect(
      container.querySelectorAll('[data-role="peak-glyph"]').length,
    ).toBe(1);
    expect(container.querySelector('[data-role="axis-bottom"]')).toBeTruthy();
    expect(container.querySelector('[data-role="axis-left"]')).toBeTruthy();
  });

  it("renders no axes when axes={false} (mini level)", () => {
    const { container } = render(
      <TracePlot traces={[model]} width={120} height={48} axes={false} />,
    );
    expect(container.querySelector('[data-role="axis-bottom"]')).toBeNull();
    expect(
      container.querySelector('[data-role="trace-line"] path'),
    ).toBeTruthy();
  });

  it("calls onAddPeak when clicking empty space", () => {
    const onAddPeak = vi.fn();
    const { container } = render(
      <TracePlot
        traces={[model]}
        width={500}
        height={300}
        interaction={{ onXDomain: () => {}, onAddPeak }}
      />,
    );
    // Click in plot body: plotPx = 200-60 = 140 ∈ [0,416], plotPy = 100-20 = 80 ∈ [0,232].
    fireEvent.click(container.querySelector("svg")!, {
      clientX: 200,
      clientY: 100,
    });
    expect(onAddPeak).toHaveBeenCalledTimes(1);
  });

  it("does not add a peak when clicking the axis margins (interior guard)", () => {
    const onAddPeak = vi.fn();
    const { container } = render(
      <TracePlot
        traces={[model]}
        width={500}
        height={300}
        interaction={{ onXDomain: () => {}, onAddPeak }}
      />,
    );
    const svg = container.querySelector("svg")!;
    // Left axis gutter (clientX=10 < margins.left=60): plotPx = 10-60 = -50 < 0.
    fireEvent.click(svg, { clientX: 10, clientY: 100 });
    // Bottom label strip (clientY=285 > top+plotHeight = 20+232 = 252): plotPy = 285-20 = 265 > 232.
    fireEvent.click(svg, { clientX: 200, clientY: 285 });
    expect(onAddPeak).not.toHaveBeenCalled();
  });

  it("renders σ band and peaks by default (no show prop)", () => {
    const { container } = render(
      <TracePlot traces={[model]} width={500} height={300} />,
    );
    // σ band path: the opacity=0.12 path inside trace-line
    expect(
      container.querySelector('[data-role="trace-line"] path[opacity="0.12"]'),
    ).toBeTruthy();
    // peaks layer
    expect(container.querySelector('[data-role="plot-peaks"]')).toBeTruthy();
  });

  it("hides the σ band when show={{ band: false }}", () => {
    const { container } = render(
      <TracePlot traces={[model]} width={500} height={300} show={{ band: false }} />,
    );
    // band path should be gone
    expect(
      container.querySelector('[data-role="trace-line"] path[opacity="0.12"]'),
    ).toBeNull();
    // but the line path (no opacity attr) should still be present
    expect(
      container.querySelector('[data-role="trace-line"] path:not([opacity])'),
    ).toBeTruthy();
  });

  it("hides the peaks layer when show={{ peaks: false }}", () => {
    const { container } = render(
      <TracePlot traces={[model]} width={500} height={300} show={{ peaks: false }} />,
    );
    expect(container.querySelector('[data-role="plot-peaks"]')).toBeNull();
  });

  describe("highlightPeakIds threading", () => {
    const modelWith2Peaks: TraceModel = {
      trace: { q: [0.1, 0.2, 0.3], I: [10, 20, 15], sigma: [1, 1, 1] },
      peaks: [
        { id: 1, q: 0.15, source: "auto", intensity: 10 },
        { id: 2, q: 0.25, source: "auto", intensity: 15 },
      ],
      phase: "Ia3d",
    };

    it("non-highlighted peaks get data-dimmed when highlightPeakIds is set", () => {
      const { container } = render(
        <TracePlot
          traces={[modelWith2Peaks]}
          width={500}
          height={300}
          highlightPeakIds={new Set([1])}
        />,
      );
      const allGs = container.querySelectorAll('[data-role="plot-peaks"] > g');
      const g2 = Array.from(allGs).find(
        (g) => g.querySelector('[data-peak-id="2"]'),
      );
      expect(g2).toBeTruthy();
      expect(g2!.getAttribute("data-dimmed")).toBe("true");
    });

    it("no peaks dimmed when highlightPeakIds is not provided", () => {
      const { container } = render(
        <TracePlot traces={[modelWith2Peaks]} width={500} height={300} />,
      );
      const dimmed = container.querySelectorAll(
        '[data-role="plot-peaks"] > g[data-dimmed]',
      );
      expect(dimmed.length).toBe(0);
    });
  });

  describe("hover / focus (engine-owned)", () => {
    it("D1: focusing the peak <g> shows the q-readout chip and blur hides it", () => {
      const { container } = render(
        <TracePlot
          traces={[model]}
          width={500}
          height={300}
          interaction={{ onXDomain: () => {} }}
        />,
      );
      const peakG = getPeakG(container);
      expect(peakG).toBeTruthy();
      // Before focus: no readout
      expect(container.querySelector('[data-role="q-readout"]')).toBeNull();
      // Focus → readout appears
      fireEvent.focus(peakG!);
      expect(container.querySelector('[data-role="q-readout"]')).toBeTruthy();
      // Blur → readout disappears
      fireEvent.blur(peakG!);
      expect(container.querySelector('[data-role="q-readout"]')).toBeNull();
    });

    it("D2: with interaction=false, peaks are not focusable (no role=button) and no readout appears", () => {
      const { container } = render(
        <TracePlot
          traces={[model]}
          width={500}
          height={300}
          interaction={false}
        />,
      );
      // No role=button on peak <g>
      expect(container.querySelector('[role="button"]')).toBeNull();
      // Attempt focus anyway — no readout
      const peakG = getPeakG(container);
      if (peakG) fireEvent.focus(peakG);
      expect(container.querySelector('[data-role="q-readout"]')).toBeNull();
    });

    it("D3: pointerLeave on the container clears the readout", () => {
      const { container } = render(
        <TracePlot
          traces={[model]}
          width={500}
          height={300}
          interaction={{ onXDomain: () => {} }}
        />,
      );
      const peakG = getPeakG(container)!;
      // Show readout via focus (deterministic in jsdom)
      fireEvent.focus(peakG);
      expect(container.querySelector('[data-role="q-readout"]')).toBeTruthy();
      // pointerLeave on the outer container div clears it
      const div = container.querySelector("div")!;
      fireEvent.pointerLeave(div);
      expect(container.querySelector('[data-role="q-readout"]')).toBeNull();
    });
  });
});
