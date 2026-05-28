/**
 * SeriesMiniWaterfall — the read-only mini-waterfall plate on a folio card
 * (R6, finding F-A; mockup `series-folio.html` `figSVG()` + `.card-fig`).
 *
 * A pure SVG render of `buildMiniWaterfall` (see `lib/series/folioFigure.ts`),
 * fed by live member snapshots. NOT the interactive `MultiTracePlot`: a folio
 * is a masonry of many cards, so a heavy Observable-Plot instance per card is
 * the wrong tool — this is the deliberate "read-only geometry variant".
 *
 * Styling follows the mockup's `svg .wf-*` rules: a dashed hair baseline, a
 * faint phase-coloured fill, a phase-coloured line, and short peak ticks.
 */
import type { SeriesMember } from "../api";
import { buildMiniWaterfall } from "../lib/series/folioFigure";

export interface SeriesMiniWaterfallProps {
  members: SeriesMember[];
}

export function SeriesMiniWaterfall({ members }: SeriesMiniWaterfallProps): JSX.Element {
  const wf = buildMiniWaterfall(members);

  return (
    <svg
      data-testid="series-mini-waterfall"
      data-trace-count={members.length}
      viewBox={`0 0 ${wf.width} ${wf.height}`}
      role="img"
      aria-label={`${members.length}-trace waterfall preview`}
      className="block h-auto w-full"
    >
      {wf.traces.map((t, i) => (
        <g key={i}>
          {/* dashed baseline (mockup .wf-base) */}
          <line
            x1={14}
            y1={t.baselineY}
            x2={wf.width - 14}
            y2={t.baselineY}
            stroke="var(--color-hair)"
            strokeWidth={1}
            strokeDasharray="2 3"
          />
          {/* faint fill under the curve (mockup .wf-fill, opacity 0.09) */}
          <path d={t.fillPath} fill={t.color} opacity={0.09} />
          {/* the integration line (mockup .wf-line) */}
          <path
            data-testid="wf-line"
            d={t.linePath}
            fill="none"
            stroke={t.color}
            strokeWidth={1.6}
            strokeLinejoin="round"
          />
          {/* short peak ticks hugging the baseline (mockup .wf-tick) */}
          {t.ticks.map((tk, k) => (
            <line
              key={k}
              x1={tk.x}
              y1={t.baselineY + 1.5}
              x2={tk.x}
              y2={t.baselineY + 1.5 + tk.h}
              stroke={tk.color}
              strokeWidth={1.4}
              strokeLinecap="round"
              opacity={tk.opacity}
            />
          ))}
        </g>
      ))}
    </svg>
  );
}
