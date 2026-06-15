import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { TracePlot, type TraceModel } from "../../src/print/plot/TracePlot";
import { phaseColor } from "../../src/phases";

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
      <TracePlot trace={model} width={500} height={300} data-testid="tp" />,
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

  it("clips the trace to the plot rect so it never overdraws the axes", () => {
    const { container } = render(
      <TracePlot trace={model} width={500} height={300} />,
    );
    // The trace-line group is wrapped in a clip-path'd <g>…
    const wrapped = container
      .querySelector('[data-role="trace-line"]')!
      .closest("g[clip-path]");
    expect(wrapped).toBeTruthy();
    // …referencing a <clipPath> whose rect spans the plot area.
    const ref = wrapped!.getAttribute("clip-path")!; // "url(#trace-clip-…)"
    const id = ref.slice(5, -1);
    const clip = container.querySelector(`clipPath#${id}`);
    expect(clip).toBeTruthy();
    const rect = clip!.querySelector("rect")!;
    expect(Number(rect.getAttribute("width"))).toBeGreaterThan(0);
    expect(Number(rect.getAttribute("height"))).toBeGreaterThan(0);
  });

  it("renders no axes when axes={false} (mini level)", () => {
    const { container } = render(
      <TracePlot trace={model} width={120} height={48} axes={false} />,
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
        trace={model}
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
        trace={model}
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
      <TracePlot trace={model} width={500} height={300} />,
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
      <TracePlot trace={model} width={500} height={300} show={{ band: false }} />,
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
      <TracePlot trace={model} width={500} height={300} show={{ peaks: false }} />,
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
          trace={modelWith2Peaks}
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
        <TracePlot trace={modelWith2Peaks} width={500} height={300} />,
      );
      const dimmed = container.querySelectorAll(
        '[data-role="plot-peaks"] > g[data-dimmed]',
      );
      expect(dimmed.length).toBe(0);
    });
  });

  describe("TracePlot single-trace API", () => {
    const singleModel: TraceModel = {
      trace: { q: [0.02, 0.05, 0.1], I: [100, 50, 10], sigma: [1, 1, 1] },
      peaks: [{ id: 1, q: 0.05, intensity: 50, source: "auto" }],
      phase: "Ia3d",
    };

    it("colours the trace line by phase (not ink-soft)", () => {
      const { container } = render(<TracePlot trace={singleModel} height={200} />);
      const paths = container.querySelectorAll('[data-role="trace-line"] path');
      const line = paths[paths.length - 1] as SVGPathElement;
      expect(line.getAttribute("stroke")).toBe(phaseColor("Ia3d"));
    });

    it("falls back to the present ink-soft neutral for a null-phase trace (BU-TRACE-DIM)", () => {
      // Unassigned traces read at full strength at rest, not the washed-out
      // ink-faint that looked like a stuck hover-dim next to phase colors.
      const { container } = render(
        <TracePlot trace={{ ...singleModel, phase: null }} height={200} />,
      );
      const paths = container.querySelectorAll('[data-role="trace-line"] path');
      const line = paths[paths.length - 1] as SVGPathElement;
      expect(line.getAttribute("stroke")).toBe("var(--color-ink-soft)");
    });
  });

  describe("hover / focus (engine-owned)", () => {
    it("D1: focusing the peak <g> shows the q-readout chip and blur hides it", () => {
      const { container } = render(
        <TracePlot
          trace={model}
          width={500}
          height={300}
          interaction={{ onXDomain: () => {} }}
        />,
      );
      const peakG = getPeakG(container);
      expect(peakG).toBeTruthy();
      // Before focus: no readout
      expect(container.querySelector('[data-role="q-readout"]')).toBeNull();
      // Focus → readout appears, showing the q-value at 3 decimals (q=0.2).
      fireEvent.focus(peakG!);
      const readout = container.querySelector('[data-role="q-readout"]');
      expect(readout).toBeTruthy();
      expect(readout!.querySelector("text")!.textContent).toBe("0.200");
      // Blur → readout disappears
      fireEvent.blur(peakG!);
      expect(container.querySelector('[data-role="q-readout"]')).toBeNull();
    });

    it("D2: with interaction=false, peaks are not focusable (no role=button) and no readout appears", () => {
      const { container } = render(
        <TracePlot
          trace={model}
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
          trace={model}
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
    // Note: the pointer-hover path and the "q-readout gated when peaks hidden"
    // edge are not unit-tested — jsdom does not carry clientX on synthetic
    // PointerEvents (arrives as 0, never hits), and with peaks hidden there is
    // no focusable peak <g> to drive the focus path. The readout is gated on
    // `layers.peaks` in TracePlot; the positive path is covered by D1.
  });

  describe("q-link (cross-panel hover wires)", () => {
    // Outward: hovering a peak (here driven through the focus path, which is the
    // deterministic jsdom mechanism — see D1; the effect keys on internal hoverId
    // so both the frame hit-test and the glyph focus share this site) emits onHoverQ.
    it("Q1: emits onHoverQ(peak.q) on internal hover and onHoverQ(undefined) on leave", () => {
      const onHoverQ = vi.fn();
      const { container } = render(
        <TracePlot
          trace={model}
          width={500}
          height={300}
          interaction={{ onXDomain: () => {} }}
          onHoverQ={onHoverQ}
        />,
      );
      const peakG = getPeakG(container)!;
      // Mount may emit onHoverQ(undefined) once; clear it so we assert the hover delta.
      onHoverQ.mockClear();
      fireEvent.focus(peakG);
      expect(onHoverQ).toHaveBeenCalledWith(0.2);
      onHoverQ.mockClear();
      fireEvent.blur(peakG);
      expect(onHoverQ).toHaveBeenCalledWith(undefined);
    });

    // Incoming: an external hoveredQ at (or very near) a peak's q lights it —
    // full parity with internal hover: the q-readout chip appears showing that q.
    it("Q2: an incoming hoveredQ lights the matching peak (q-readout chip appears)", () => {
      const { container } = render(
        <TracePlot
          trace={model}
          width={500}
          height={300}
          interaction={{ onXDomain: () => {} }}
          hoveredQ={0.2}
        />,
      );
      const readout = container.querySelector('[data-role="q-readout"]');
      expect(readout).toBeTruthy();
      expect(readout!.querySelector("text")!.textContent).toBe("0.200");
    });

    // No feedback loop: an external hoveredQ with no real pointer/focus hover must
    // NOT cause onHoverQ to fire with a defined number (page→hoveredQ→onHoverQ→…).
    it("Q3: incoming hoveredQ does not emit onHoverQ (no feedback loop)", () => {
      const onHoverQ = vi.fn();
      render(
        <TracePlot
          trace={model}
          width={500}
          height={300}
          interaction={{ onXDomain: () => {} }}
          hoveredQ={0.2}
          onHoverQ={onHoverQ}
        />,
      );
      // A mount-time onHoverQ(undefined) is acceptable; a defined-number call is not.
      const definedCalls = onHoverQ.mock.calls.filter(
        ([q]) => typeof q === "number",
      );
      expect(definedCalls).toEqual([]);
    });
  });
});
