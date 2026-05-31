import { type Axis1D } from "./projection";
import { formatAxis } from "../../lib/plot/formatAxis";

export interface AxisProps {
  axis: Axis1D;
  orientation: "bottom" | "left";
  plotWidth: number;
  plotHeight: number;
  label?: string;
}

const TICK_LABEL_STYLE: React.CSSProperties = {
  fontFamily: "var(--font-mono)",
  fontSize: 10,
  fill: "var(--color-ink-soft)",
};

const AXIS_LABEL_STYLE: React.CSSProperties = {
  fontFamily: "var(--font-sans)",
  fontSize: 11,
  fill: "var(--color-ink-soft)",
};

export function Axis({
  axis,
  orientation,
  plotWidth,
  plotHeight,
  label,
}: AxisProps): JSX.Element {
  const ticks = axis.ticks(6);
  const isBottom = orientation === "bottom";
  return (
    <g data-role={`axis-${orientation}`} aria-hidden="true">
      {ticks.map((t) => {
        const p = axis.to(t);
        if (!Number.isFinite(p)) return null;
        return isBottom ? (
          <g key={t} transform={`translate(${p},${plotHeight})`}>
            <line y2={5} stroke="var(--color-ink-faint)" strokeWidth={1} />
            <text y={16} textAnchor="middle" style={TICK_LABEL_STYLE}>
              {formatAxis(t)}
            </text>
          </g>
        ) : (
          <g key={t} transform={`translate(0,${p})`}>
            <line x2={-5} stroke="var(--color-ink-faint)" strokeWidth={1} />
            <text x={-8} dy="0.32em" textAnchor="end" style={TICK_LABEL_STYLE}>
              {formatAxis(t)}
            </text>
          </g>
        );
      })}
      {label != null &&
        (isBottom ? (
          <text
            x={plotWidth / 2}
            y={36}
            textAnchor="middle"
            style={AXIS_LABEL_STYLE}
          >
            {label}
          </text>
        ) : (
          <text
            transform={`translate(-40,${plotHeight / 2}) rotate(-90)`}
            textAnchor="middle"
            style={AXIS_LABEL_STYLE}
          >
            {label}
          </text>
        ))}
    </g>
  );
}
