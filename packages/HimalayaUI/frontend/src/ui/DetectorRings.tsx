import type { ImageSize } from "../lib/detectorRings";

export interface DetectorRingsView {
  size: ImageSize;
  cx: number;
  cy: number;
  rings: readonly { r: number; color: string }[];
}

interface Props {
  view: DetectorRingsView;
  fit: "cover" | "contain";
}

export function DetectorRings({ view, fit }: Props): JSX.Element {
  const { size, cx, cy, rings } = view;
  return (
    <svg
      data-testid="detector-rings"
      aria-hidden="true"
      viewBox={`0 0 ${size.width} ${size.height}`}
      preserveAspectRatio={fit === "cover" ? "xMidYMid slice" : "xMidYMid meet"}
      className="pointer-events-none absolute inset-0 h-full w-full"
    >
      {rings.map(({ r, color }, i) => (
        <circle
          key={i}
          data-role="det-ring"
          cx={cx}
          cy={cy}
          r={r}
          fill="none"
          stroke={color}
          strokeOpacity={0.85}
          strokeWidth={1.25}
          vectorEffect="non-scaling-stroke"
        />
      ))}
    </svg>
  );
}
