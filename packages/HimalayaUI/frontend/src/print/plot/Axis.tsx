import { type Axis1D, axisTicks } from "./projection";
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
  fontSize: 10.5,
  fill: "var(--color-ink-faint)",
};

const AXIS_TITLE_STYLE: React.CSSProperties = {
  fontFamily: "var(--font-sans)",
  fontSize: 11.5,
  fontWeight: 600,
  fill: "var(--color-ink-soft)",
};

export function Axis({
  axis,
  orientation,
  plotWidth,
  plotHeight,
  label,
}: AxisProps): JSX.Element {
  const ticks = axisTicks(axis);
  const isBottom = orientation === "bottom";
  return (
    <g data-role={`axis-${orientation}`} aria-hidden="true">
      {ticks
        .filter((t) => t.kind === "major")
        .map((t) => {
          const p = axis.to(t.value);
          if (!Number.isFinite(p)) return null;
          return isBottom ? (
            <line
              key={`g${t.value}`}
              data-role="gridline"
              x1={p}
              y1={0}
              x2={p}
              y2={plotHeight}
              stroke="var(--color-hair)"
              strokeWidth={1}
            />
          ) : (
            <line
              key={`g${t.value}`}
              data-role="gridline"
              x1={0}
              y1={p}
              x2={plotWidth}
              y2={p}
              stroke="var(--color-hair)"
              strokeWidth={1}
            />
          );
        })}
      {isBottom ? (
        <line
          data-role="axis-spine"
          x1={0}
          y1={plotHeight}
          x2={plotWidth}
          y2={plotHeight}
          stroke="var(--color-hair-strong)"
          strokeWidth={1}
        />
      ) : (
        <line
          data-role="axis-spine"
          x1={0}
          y1={0}
          x2={0}
          y2={plotHeight}
          stroke="var(--color-hair-strong)"
          strokeWidth={1}
        />
      )}
      {ticks.map((t) => {
        const p = axis.to(t.value);
        if (!Number.isFinite(p)) return null;
        const len = t.kind === "major" ? 6 : t.kind === "mid" ? 4 : 3;
        const stroke =
          t.kind === "minor" ? "var(--color-hair)" : "var(--color-hair-strong)";
        const labelled = t.kind !== "minor";
        return isBottom ? (
          <g key={t.value} transform={`translate(${p},${plotHeight})`}>
            <line y2={len} stroke={stroke} strokeWidth={1} />
            {labelled && (
              <text y={len + 11} textAnchor="middle" style={TICK_LABEL_STYLE}>
                {formatAxis(t.value)}
              </text>
            )}
          </g>
        ) : (
          <g key={t.value} transform={`translate(0,${p})`}>
            <line x2={-len} stroke={stroke} strokeWidth={1} />
            {labelled && (
              <text
                x={-(len + 4)}
                dy="0.32em"
                textAnchor="end"
                style={TICK_LABEL_STYLE}
              >
                {formatAxis(t.value)}
              </text>
            )}
          </g>
        );
      })}
      {label != null &&
        (isBottom ? (
          <text
            x={plotWidth / 2}
            y={plotHeight + 40}
            textAnchor="middle"
            style={AXIS_TITLE_STYLE}
          >
            {label}
          </text>
        ) : (
          <text
            transform={`translate(-46,${plotHeight / 2}) rotate(-90)`}
            textAnchor="middle"
            style={AXIS_TITLE_STYLE}
          >
            {label}
          </text>
        ))}
    </g>
  );
}
