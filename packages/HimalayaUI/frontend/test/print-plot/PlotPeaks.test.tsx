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

  it("clamps a peak taller than the ceiling into the plot body (FO-CLIPMARK)", () => {
    // intensity 100000 maps far above the yDomain top (30) → without clamping
    // the glyph apex lands above the plot top (negative y), in the dead margin
    // where pointer clicks are ignored. The clamp keeps the whole glyph at
    // y >= 0 so the marker stays visible and clickable at the top edge.
    const { container } = render(
      <svg>
        <PlotPeaks
          peaks={[{ id: 1, q: 0.2, intensity: 100000, source: "auto" }]}
          projection={proj}
          color="var(--color-accent)"
        />
      </svg>,
    );
    const poly = container.querySelector('[data-role="peak-glyph"] polygon')!;
    const ys = (poly.getAttribute("points") ?? "")
      .trim()
      .split(/\s+/)
      .map((pt) => Number(pt.split(",")[1]));
    expect(Math.min(...ys)).toBeGreaterThanOrEqual(0); // no vertex above the top
    expect(Math.max(...ys)).toBeLessThanOrEqual(200); // within plotHeight
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
    // Paint flows through the animatable color channel (FO-DIM), not a literal.
    expect(qline.getAttribute("stroke")).toBe("currentColor");
  });

  it("a hot peak does not change the mark — no halo/ring, emphasis is the q-line", () => {
    const restColor = "var(--color-success)";
    const { container, rerender } = render(
      <svg>
        <PlotPeaks
          peaks={[{ id: 1, q: 0.2, intensity: 20, source: "auto" }]}
          projection={proj}
          color={restColor}
        />
      </svg>,
    );
    const restingStroke = container
      .querySelector('[data-role="peak-glyph"] [data-shape]')!
      .getAttribute("stroke");
    rerender(
      <svg>
        <PlotPeaks
          peaks={[{ id: 1, q: 0.2, intensity: 20, source: "auto", hot: true }]}
          projection={proj}
          color={restColor}
        />
      </svg>,
    );
    // No separate halo/ring element, and the mark's stroke is unchanged by hot.
    expect(container.querySelector('[data-role="peak-halo"]')).toBeNull();
    const hotStroke = container
      .querySelector('[data-role="peak-glyph"] [data-shape]')!
      .getAttribute("stroke");
    expect(hotStroke).toBe(restingStroke);
    expect(hotStroke).not.toBe("var(--color-accent)");
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
    // The auto peak is filled → the per-peak color (not the layer color) lands
    // on the wrapping <g>'s CSS `color`; the polygon paints via currentColor so
    // the global 120ms ease-out `color` transition animates dim/release.
    const peakG = container.querySelector(
      '[data-role="plot-peaks"] > g',
    ) as SVGGElement;
    expect(peakG.style.color).toBe("var(--color-success)");
    expect(
      container.querySelector('[data-role="peak-glyph"] polygon')?.getAttribute('fill'),
    ).toBe("currentColor");
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
    const peakG = container.querySelector(
      '[data-role="plot-peaks"] > g',
    ) as SVGGElement;
    expect(peakG.style.color).toBe("var(--color-accent)");
    expect(
      container.querySelector('[data-role="peak-glyph"] polygon')?.getAttribute('fill'),
    ).toBe("currentColor");
  });

  it("onPeakFocus alone (read-only marks) does NOT make the <g> focusable or a button", () => {
    // The hover q-link (onPeakFocus) is independent of editability. A mark that
    // exposes role=button while the plate is disarmed would lie: activating it
    // would do nothing. Focusability is gated on onPeakActivate alone.
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
    expect(peakG!.getAttribute("tabindex")).toBeNull();
    expect(peakG!.getAttribute("role")).toBeNull();
    expect(peakG!.getAttribute("aria-label")).toBeNull();
  });

  it("peak <g> has tabIndex/role/aria-label/aria-keyshortcuts when onPeakActivate is provided", () => {
    const { container } = render(
      <svg>
        <PlotPeaks
          peaks={[{ id: 1, q: 0.2, intensity: 20, source: "auto" }]}
          projection={proj}
          color="var(--color-accent)"
          onPeakActivate={vi.fn()}
        />
      </svg>,
    );
    const peakG = container.querySelector('[data-role="plot-peaks"] > g');
    expect(peakG).toBeTruthy();
    expect(peakG!.getAttribute("tabindex")).toBe("0");
    expect(peakG!.getAttribute("role")).toBe("button");
    // FO-RESCORE2 F11: the label names provenance (Auto vs Manual).
    expect(peakG!.getAttribute("aria-label")).toBe("Auto peak at q = 0.2000");
    expect(peakG!.getAttribute("aria-keyshortcuts")).toBe("Enter Space");
  });

  it("names peak provenance and excluded state in the aria-label (FO-RESCORE2 F11)", () => {
    const { container } = render(
      <svg>
        <PlotPeaks
          peaks={[
            { id: 1, q: 0.18, intensity: 20, source: "manual" },
            { id: 2, q: 0.24, intensity: 20, source: "auto", excluded: true },
          ]}
          projection={proj}
          color="var(--color-accent)"
          onPeakActivate={vi.fn()}
        />
      </svg>,
    );
    const gs = container.querySelectorAll('[data-role="plot-peaks"] > g');
    expect(gs[0]!.getAttribute("aria-label")).toBe("Manual peak at q = 0.1800");
    expect(gs[1]!.getAttribute("aria-label")).toBe("Auto peak at q = 0.2400 (excluded)");
  });

  it("a peak OUTSIDE the visible x-domain is not focusable even when armed (FO-ZOOMEDIT)", () => {
    // proj's window is [0.1, 0.3]. A peak at q=0.5 maps off-window: its glyph
    // and focus ring are clipped (invisible), so an armed role=button/tabindex
    // would be a phantom tab stop with no visible target. The <g> still renders
    // (data-peak-id is kept for the focus-re-anchor lookup) but is NOT a control.
    const { container } = render(
      <svg>
        <PlotPeaks
          peaks={[
            { id: 1, q: 0.2, intensity: 20, source: "auto" },
            { id: 2, q: 0.5, intensity: 20, source: "auto" },
          ]}
          projection={proj}
          color="var(--color-accent)"
          onPeakActivate={vi.fn()}
        />
      </svg>,
    );
    const inWindow = container.querySelector('g[data-peak-id="1"]');
    const offWindow = container.querySelector('g[data-peak-id="2"]');
    // In-window peak is a real, focusable control.
    expect(inWindow!.getAttribute("tabindex")).toBe("0");
    expect(inWindow!.getAttribute("role")).toBe("button");
    // Off-window peak still renders (for re-anchor lookup) but is inert.
    expect(offWindow).toBeTruthy();
    expect(offWindow!.getAttribute("tabindex")).toBeNull();
    expect(offWindow!.getAttribute("role")).toBeNull();
    expect(offWindow!.getAttribute("aria-label")).toBeNull();
  });

  it("Enter / Space fire onPeakActivate(id); activation carries no modifier (alt retired)", () => {
    const spy = vi.fn();
    const { container } = render(
      <svg>
        <PlotPeaks
          peaks={[{ id: 7, q: 0.2, intensity: 20, source: "auto" }]}
          projection={proj}
          color="var(--color-accent)"
          onPeakActivate={spy}
        />
      </svg>,
    );
    const peakG = container.querySelector('[data-role="plot-peaks"] > g')!;
    // Enter activates. fireEvent returns false when default was prevented.
    expect(fireEvent.keyDown(peakG, { key: "Enter" })).toBe(false);
    expect(spy).toHaveBeenLastCalledWith(7);
    // Space activates (and must preventDefault so the page does not scroll).
    expect(fireEvent.keyDown(peakG, { key: " " })).toBe(false);
    expect(spy).toHaveBeenLastCalledWith(7);
    // Alt is retired: Alt+Enter is a plain activation, no second argument.
    fireEvent.keyDown(peakG, { key: "Enter", altKey: true });
    expect(spy).toHaveBeenLastCalledWith(7);
    expect(spy).toHaveBeenCalledTimes(3);
    // Unrelated keys are ignored.
    fireEvent.keyDown(peakG, { key: "a" });
    expect(spy).toHaveBeenCalledTimes(3);
    // Key-repeat must not machine-gun mutations.
    fireEvent.keyDown(peakG, { key: "Enter", repeat: true });
    expect(spy).toHaveBeenCalledTimes(3);
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

  describe("focusRequest — keyboard re-anchor after a destructive edit", () => {
    it("the focusable wrapper <g> carries data-peak-id so the survivor can be found", () => {
      const { container } = render(
        <svg>
          <PlotPeaks
            peaks={[{ id: 7, q: 0.2, intensity: 20, source: "auto" }]}
            projection={proj}
            color="var(--color-accent)"
            onPeakActivate={vi.fn()}
          />
        </svg>,
      );
      const wrapper = container.querySelector(
        '[data-role="plot-peaks"] > g[data-peak-id="7"]',
      );
      expect(wrapper).toBeTruthy();
      expect(wrapper!.getAttribute("role")).toBe("button");
    });

    it("calls onFocusFallback once when the requested id has no surviving mark", () => {
      const fallback = vi.fn();
      const { rerender } = render(
        <svg>
          <PlotPeaks
            peaks={[{ id: 1, q: 0.2, intensity: 20, source: "auto" }]}
            projection={proj}
            color="var(--color-accent)"
            onPeakActivate={vi.fn()}
            onFocusFallback={fallback}
          />
        </svg>,
      );
      // A request for an absent id (no survivor) → fallback fires exactly once.
      rerender(
        <svg>
          <PlotPeaks
            peaks={[{ id: 1, q: 0.2, intensity: 20, source: "auto" }]}
            projection={proj}
            color="var(--color-accent)"
            onPeakActivate={vi.fn()}
            onFocusFallback={fallback}
            focusRequest={{ id: null, nonce: 1 }}
          />
        </svg>,
      );
      expect(fallback).toHaveBeenCalledTimes(1);
    });

    it("does not re-fire the same nonce (idempotent per request)", () => {
      const fallback = vi.fn();
      const tree = (nonce: number) => (
        <svg>
          <PlotPeaks
            peaks={[{ id: 1, q: 0.2, intensity: 20, source: "auto" }]}
            projection={proj}
            color="var(--color-accent)"
            onPeakActivate={vi.fn()}
            onFocusFallback={fallback}
            focusRequest={{ id: null, nonce }}
          />
        </svg>
      );
      const { rerender } = render(tree(1));
      expect(fallback).toHaveBeenCalledTimes(1);
      // An unrelated re-render with the SAME nonce must not re-fire (a foreign
      // SSE re-render never steals focus).
      rerender(tree(1));
      expect(fallback).toHaveBeenCalledTimes(1);
      // A NEW nonce fires again.
      rerender(tree(2));
      expect(fallback).toHaveBeenCalledTimes(2);
    });
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

      // Dim colour lands on the wrapping <g>'s CSS `color` (ink-faint); the
      // glyph paints via currentColor so the dim/release eases through the
      // global 120ms ease-out `color` transition instead of snapping.
      expect((g2 as SVGGElement).style.color).toBe("var(--color-ink-faint)");
      expect((g3 as SVGGElement).style.color).toBe("var(--color-ink-faint)");
      expect(
        g2!.querySelector('[data-role="peak-glyph"] polygon')?.getAttribute("fill"),
      ).toBe("currentColor");
      expect(
        g3!.querySelector('[data-role="peak-glyph"] polygon')?.getAttribute("fill"),
      ).toBe("currentColor");
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
