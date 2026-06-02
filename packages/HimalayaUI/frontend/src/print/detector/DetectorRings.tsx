import { type CSSProperties } from "react";
import type { RingPlacement } from "./detectorGeometry";

interface Props {
  /** Beam center in normalized image coords (0-1). The ring origin. */
  beamCenter: { x: number; y: number };
  /** Rings with radii in normalized image coords (fraction of image width). */
  rings: RingPlacement[];
  /** q hovered elsewhere (trace peak / comb). The matching ring lights to the
   *  accent q-link. Props-driven so the renderer is pure; the composite threads
   *  Zustand `hoveredQ`. */
  hoveredQ?: number;
  /** Fired on hit-ring enter (q) / leave (undefined). Omit -> inert rings. */
  onHoverQ?: (q?: number) => void;
  /** Displayed image aspect ratio (width / height). The viewBox is stretched to
   *  the (non-square) image rect by preserveAspectRatio="none", which would
   *  squash a <circle> into an ellipse; we pre-correct the y-radius by this ratio
   *  so rings render as TRUE circles. Defaults to 1 (square). */
  imageAspect?: number;
  /** Orientation when the parent drives it (mirrors DetectorImage). */
  orient?: "portrait" | "landscape";
}

export function DetectorRings({ beamCenter, rings, hoveredQ, onHoverQ, imageAspect, orient }: Props): JSX.Element {
  const svgStyle: CSSProperties = {
    position: "absolute", inset: 0, width: "100%", height: "100%",
    pointerEvents: "none", // hit rings re-enable individually
    ...(orient === "landscape" ? { transform: "rotate(90deg)", transformOrigin: "center" } : {}),
  };
  const TOL = 1e-6;
  // r is a fraction of image WIDTH. After preserveAspectRatio="none" stretches the
  // 0..1 viewBox to a Wp×Hp rect, an x-radius rx renders rx·Wp px and a y-radius ry
  // renders ry·Hp px. To get equal pixel radii (a circle), ry = rx·(Wp/Hp) = r·aspect.
  const aspect = imageAspect ?? 1;

  return (
    <div className="pointer-events-none absolute inset-0 z-10">
      {/* viewBox is the normalized image box; preserveAspectRatio="none" stretches
          it to the displayed image rect so beam center + radii register. The frame's
          overflow:hidden clips the arcs of an off-center beam. */}
      <svg data-testid="detector-rings" data-orient={orient ?? "portrait"}
           viewBox="0 0 1 1" preserveAspectRatio="none" style={svgStyle} aria-hidden="true">
        {rings.map(({ q, r, color, ghost }, i) => {
          const hot = hoveredQ !== undefined && Math.abs(q - hoveredQ) <= TOL;
          // Hover keeps the ring's OWN colour — emphasis is a wider, more opaque
          // stroke, never a terracotta-accent recolour (matches the trace-plot
          // hover rule: the accent is not a hover highlight).
          const stroke = color ?? "var(--color-ink-faint)";
          // Radii/strokes are in the normalized (0..1) space; vector-effect keeps
          // stroke widths crisp despite preserveAspectRatio="none" scaling.
          const sharpW = hot ? 1.8 : ghost ? 0.8 : 1.0;
          const sharpOp = hot ? 0.95 : ghost ? 0.45 : 0.7;
          return (
            <g key={i} data-role="det-ring" data-ring-q={q}
               data-hot={hot ? "true" : undefined} data-ghost={ghost ? "true" : undefined}>
              <ellipse data-role="ring-glow" cx={beamCenter.x} cy={beamCenter.y} rx={r} ry={r * aspect}
                fill="none" stroke={stroke} strokeWidth={2.4} vectorEffect="non-scaling-stroke"
                opacity={hot ? 0.4 : ghost ? 0.06 : 0.18} style={{ pointerEvents: "none" }} />
              <ellipse data-role="ring-sharp" cx={beamCenter.x} cy={beamCenter.y} rx={r} ry={r * aspect}
                fill="none" stroke={stroke} strokeWidth={sharpW} vectorEffect="non-scaling-stroke"
                strokeDasharray={ghost ? "2 2.5" : undefined} opacity={sharpOp}
                style={{ pointerEvents: "none" }} />
              <ellipse data-role="ring-hit" cx={beamCenter.x} cy={beamCenter.y} rx={r} ry={r * aspect}
                fill="none" stroke="transparent" strokeWidth={5} vectorEffect="non-scaling-stroke"
                style={{ pointerEvents: onHoverQ ? "stroke" : "none", cursor: onHoverQ ? "pointer" : "default" }}
                {...(onHoverQ ? { onMouseEnter: () => onHoverQ(q), onMouseLeave: () => onHoverQ(undefined) } : {})} />
            </g>
          );
        })}
      </svg>
    </div>
  );
}
