import type { ImageSize, RingSide } from "../lib/detectorRings";

export interface DetectorRingsView {
  size: ImageSize;
  cx: number;
  cy: number;
  side: RingSide;
  rings: readonly { r: number; color: string }[];
}

interface Props {
  view: DetectorRingsView;
  fit: "cover" | "contain";
}

export function DetectorRings({ view, fit }: Props): JSX.Element {
  const { size, cx, cy, side, rings } = view;
  const arc = (r: number): string => {
    const vertical = side === "left" || side === "right";
    const sweep = side === "right" || side === "top" ? 1 : 0;
    return vertical
      ? `M ${cx} ${cy - r} A ${r} ${r} 0 0 ${sweep} ${cx} ${cy + r}`
      : `M ${cx - r} ${cy} A ${r} ${r} 0 0 ${sweep} ${cx + r} ${cy}`;
  };
  return (
    <svg
      data-testid="detector-rings"
      aria-hidden="true"
      viewBox={`0 0 ${size.width} ${size.height}`}
      preserveAspectRatio={fit === "cover" ? "xMidYMid slice" : "xMidYMid meet"}
      className="pointer-events-none absolute inset-0 h-full w-full"
    >
      {rings.map(({ r, color }, i) => (
        <path
          key={i}
          data-role="det-ring"
          d={arc(r)}
          fill="none"
          stroke={color}
          strokeOpacity={0.85}
          strokeWidth={1.5}
          strokeDasharray="6 4"
          vectorEffect="non-scaling-stroke"
        />
      ))}
    </svg>
  );
}
