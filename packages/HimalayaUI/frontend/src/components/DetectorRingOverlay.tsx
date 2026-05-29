import { useEffect, useRef, useState, type CSSProperties } from "react";
import { useAppState } from "../state";
import { decideOrient } from "../lib/detectorOrient";
import { qToRadius, nearestRingQ, RING_VIEWBOX } from "../lib/qRing";

interface Props {
  /** The peak q-values to draw as rings. */
  peakQs: number[];
  /** Orientation, when the parent drives it (mirrors DetectorImage). When
   *  omitted the overlay self-observes its wrapper. */
  orient?: "portrait" | "landscape";
  /** Tolerance (in q) for matching hoveredQ to a ring. */
  matchTol?: number;
}

const CENTER = RING_VIEWBOX / 2;

/**
 * DetectorRingOverlay — the q-link ring layer over FocusDetectorPanel's
 * detector image (#180). One ring per peak q, positioned by the normalized
 * `qToRadius` map (no detector geometry exists — rings are presentational).
 *
 * Rotation-aware: in landscape the overlay applies the same `rotate(90deg)`
 * as DetectorImage's canvas, via the shared `decideOrient` rule, so rings
 * stay registered with the (rotated) image. We do NOT touch DetectorImage's
 * own rotation logic.
 *
 * Hovering a ring sets `hoveredQ`; a ring whose q matches `hoveredQ` lights.
 */
export function DetectorRingOverlay({
  peakQs, orient: orientProp, matchTol,
}: Props): JSX.Element {
  const wrapperRef = useRef<HTMLDivElement>(null);
  const [selfOrient, setSelfOrient] =
    useState<"portrait" | "landscape">("portrait");
  const orient = orientProp ?? selfOrient;

  const hoveredQ = useAppState((s) => s.hoveredQ);
  const setHoveredQ = useAppState((s) => s.setHoveredQ);

  // Self-observe only when the parent doesn't drive orient.
  //
  // Defensive / rotation-aware-by-contract: I4.3 (#180) requires a
  // rotation-aware overlay. In the current FocusDetectorPanel composition the
  // detector sits in an `aspect-square` box (containerAspect ~ 1), so
  // `decideOrient` never selects landscape and this stays portrait — and even
  // if it flipped, concentric rings are rotation-invariant so `rotate(90deg)`
  // is a visual no-op. We keep it (rather than hardcoding portrait) so the
  // overlay remains correct if a future layout gives the detector a
  // non-square slot, matching DetectorImage's own rotation contract via the
  // shared `decideOrient` rule.
  useEffect(() => {
    if (orientProp !== undefined) return;
    const wrapper = wrapperRef.current;
    if (!wrapper) return;
    const recompute = (): void => {
      const r = decideOrient({
        containerW: wrapper.clientWidth,
        containerH: wrapper.clientHeight,
        // square overlay → image aspect is 1; rings sit on the image center
        imageW: 1, imageH: 1,
        viewportW: typeof window !== "undefined" ? window.innerWidth : 0,
      });
      setSelfOrient(r.orient);
    };
    recompute();
    const ro = new ResizeObserver(recompute);
    ro.observe(wrapper);
    return () => ro.disconnect();
  }, [orientProp]);

  const qLo = peakQs.length ? Math.min(...peakQs) : 0;
  const qHi = peakQs.length ? Math.max(...peakQs) : 1;
  const tol = matchTol ?? Math.max((qHi - qLo) * 0.02, 1e-6);
  const matched = nearestRingQ(hoveredQ, peakQs, tol);

  const svgStyle: CSSProperties = {
    position: "absolute",
    inset: 0,
    width: "100%",
    height: "100%",
    pointerEvents: "none", // rings re-enable it individually
    ...(orient === "landscape"
      ? { transform: "rotate(90deg)", transformOrigin: "center" }
      : {}),
  };

  return (
    <div ref={wrapperRef} className="pointer-events-none absolute inset-0 z-10">
      <svg
        data-testid="detector-ring-overlay"
        data-orient={orient}
        viewBox={`0 0 ${RING_VIEWBOX} ${RING_VIEWBOX}`}
        style={svgStyle}
        aria-hidden="true"
      >
        {peakQs.map((q) => {
          const r = qToRadius(q, qLo, qHi);
          const hot = matched !== undefined && q === matched;
          return (
            <g key={q}>
              {/* Visible ring. */}
              <circle
                data-testid={`detector-ring-q-${q}`}
                data-hot={hot ? "true" : "false"}
                cx={CENTER}
                cy={CENTER}
                r={r}
                fill="none"
                stroke={hot ? "var(--color-accent)" : "var(--color-ink-faint)"}
                strokeWidth={hot ? 1.6 : 0.8}
                opacity={hot ? 0.95 : 0.5}
                style={{ pointerEvents: "none" }}
              />
              {/* Wide transparent hit-ring (mockup parity) so the thin visible
                  stroke isn't a hard hover target and the pointer lands even
                  over the canvas. */}
              <circle
                data-testid={`detector-ring-hit-${q}`}
                cx={CENTER}
                cy={CENTER}
                r={r}
                fill="none"
                stroke="transparent"
                strokeWidth={5}
                style={{ pointerEvents: "stroke", cursor: "pointer" }}
                onMouseEnter={() => setHoveredQ(q)}
                // Guarded clear: crossing ring A -> B fires enter(B) before
                // leave(A); an unconditional clear here would clobber B's
                // highlight back to undefined. Only clear if we're still the
                // active q (a real leave, not a hand-off to a neighbour).
                onMouseLeave={() =>
                  setHoveredQ(
                    useAppState.getState().hoveredQ === q ? undefined
                      : useAppState.getState().hoveredQ,
                  )}
              />
            </g>
          );
        })}
      </svg>
    </div>
  );
}
