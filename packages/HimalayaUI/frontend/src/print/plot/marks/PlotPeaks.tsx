import { type Projection } from "../projection";
import { PeakGlyph } from "../../ui/PeakGlyph";
import { peakGlyph } from "../../ui/peakMark";

export interface PlotPeak {
  id: number;
  q: number;
  intensity?: number | null;
  source: "auto" | "manual";
  excluded?: boolean;
  predictedAbsent?: boolean;
  hot?: boolean;
}

export interface PlotPeaksProps {
  peaks: PlotPeak[];
  projection: Projection;
  /** Resolved colour for the glyphs (phase colour, usually). */
  color: string;
  /** y data-value for peaks lacking an intensity (anchors them to baseline). */
  baselineI?: number;
  /** Paper colour threaded to PeakGlyph's halo (export-parity). */
  paperColor?: string;
}

export function PlotPeaks({
  peaks,
  projection,
  color,
  baselineI,
  paperColor,
}: PlotPeaksProps): JSX.Element {
  const { x, y } = projection;
  return (
    <g data-role="plot-peaks">
      {peaks.map((p) => {
        if (p.id < 0) return null;
        const px = x.to(p.q);
        if (!Number.isFinite(px)) return null;
        const iVal = p.intensity ?? baselineI;
        const py =
          iVal != null && Number.isFinite(iVal) ? y.to(iVal) : y.range[0];
        const descriptor = peakGlyph({
          source: p.source,
          color,
          ...(p.predictedAbsent ? { predictedAbsent: true } : {}),
          ...(p.excluded ? { excluded: true } : {}),
          ...(p.hot ? { hot: true } : {}),
        });
        return (
          <g key={p.id}>
            {p.hot ? (
              <line
                data-role="peak-qline"
                x1={px}
                y1={y.range[0]}
                x2={px}
                y2={y.range[1]}
                stroke={color}
                strokeWidth={1}
                opacity={0.6}
              />
            ) : null}
            <PeakGlyph
              descriptor={descriptor}
              x={px}
              y={py}
              dataPeakId={p.id}
              {...(paperColor ? { haloStroke: paperColor } : {})}
            />
          </g>
        );
      })}
    </g>
  );
}
