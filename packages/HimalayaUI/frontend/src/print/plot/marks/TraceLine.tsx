import { line as d3line, area as d3area } from "d3-shape";
import { type Projection } from "../projection";
import { type Trace } from "../../../api";

export interface TraceLineProps {
  trace: Trace;
  projection: Projection;
  /** Draw the ±σ band behind the line (default true; needs trace.sigma). */
  band?: boolean;
}

export function TraceLine({
  trace,
  projection,
  band = true,
}: TraceLineProps): JSX.Element {
  const { x, y } = projection;
  const n = Math.min(trace.q.length, trace.I.length);
  const idx = Array.from({ length: n }, (_, i) => i);
  const valid = (i: number): boolean =>
    Number.isFinite(trace.q[i]!) &&
    trace.q[i]! > 0 &&
    Number.isFinite(trace.I[i]!) &&
    trace.I[i]! > 0;

  const linePath =
    d3line<number>()
      .defined(valid)
      .x((i) => x.to(trace.q[i]!))
      .y((i) => y.to(trace.I[i]!))(idx) ?? "";

  const sigma = trace.sigma;
  const bandPath =
    band && sigma
      ? d3area<number>()
          .defined((i) => valid(i) && Number.isFinite(sigma[i]!))
          .x((i) => x.to(trace.q[i]!))
          .y0((i) => y.to(Math.max(trace.I[i]! - sigma[i]!, 1e-9)))
          .y1((i) => y.to(trace.I[i]! + sigma[i]!))(idx) ?? ""
      : "";

  return (
    <g data-role="trace-line">
      {bandPath ? (
        <path d={bandPath} fill="var(--color-ink-soft)" opacity={0.12} stroke="none" />
      ) : null}
      <path
        d={linePath}
        fill="none"
        stroke="var(--color-ink-soft)"
        strokeWidth={1.8}
        strokeLinejoin="round"
      />
    </g>
  );
}
