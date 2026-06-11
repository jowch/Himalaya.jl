import { RESID_DOMAIN, RESID_BAND } from "../comb";

export type ResidVocabItem = "in-tolerance" | "out-of-tolerance" | "off-scale";

interface ResidualLegendProps {
  /** PLACEMENT-ONLY. */
  className?: string;
}

const ITEMS: ReadonlyArray<ResidVocabItem> = [
  "in-tolerance",
  "out-of-tolerance",
  "off-scale",
];

const LABELS: Record<ResidVocabItem, string> = {
  "in-tolerance": "in tolerance",
  "out-of-tolerance": "out of tolerance",
  "off-scale": "off scale",
};

// Neutral ink, like CombLegend: residual rows are per-phase colored, but the
// vocabulary itself is shape-based (filled dot / hollow dot / chevron).
const NEUTRAL_COLOR = "var(--color-ink-soft)";

// `±2.2%` / `±3%` from the SAME constants ResidualChart draws with, so the note
// cannot drift from the geometry. `+…toFixed(1)` drops a trailing ".0".
function pct(v: number): string {
  return `±${+(v * 100).toFixed(1)}%`;
}

// Glyph geometry mirrors what ResidualChart actually draws: filled circle (in
// tolerance), stroke-only circle (out of tolerance, strokeWidth 1.6), and the
// off-scale edge chevron (half-width 3.2, height 4, strokeWidth 1.3).
function glyph(item: ResidVocabItem): JSX.Element {
  switch (item) {
    case "in-tolerance":
      return <circle data-shape="dot-filled" cx={9} cy={7} r={3} fill={NEUTRAL_COLOR} />;
    case "out-of-tolerance":
      return (
        <circle
          data-shape="dot-hollow"
          cx={9} cy={7} r={3}
          fill="none" stroke={NEUTRAL_COLOR} strokeWidth={1.6}
        />
      );
    case "off-scale":
      return (
        <path
          data-shape="chevron"
          d="M5.8 9 L9 5 L12.2 9"
          fill="none" stroke={NEUTRAL_COLOR} strokeWidth={1.3}
          strokeLinejoin="round" strokeLinecap="round"
        />
      );
  }
}

/** Glyph vocabulary for the residual ("indexing space") view, mirroring CombLegend's
 *  row style, plus the one quiet annotation naming the invisible constants. */
export function ResidualLegend({ className }: ResidualLegendProps): JSX.Element {
  return (
    <div
      data-role="resid-legend"
      className={`flex flex-wrap items-center gap-3 border-t border-hair pt-2${className ? ` ${className}` : ""}`}
    >
      {ITEMS.map((item) => (
        <span key={item} className="inline-flex items-center gap-1.5">
          <svg
            data-role="resid-legend-glyph"
            width="18"
            height="14"
            viewBox="0 0 18 14"
            aria-hidden="true"
          >
            {glyph(item)}
          </svg>
          <span className="font-mono text-caption text-ink-soft">
            {LABELS[item]}
          </span>
        </span>
      ))}
      {/* Names the otherwise-invisible drawing constants once, not per row. */}
      <span data-role="resid-legend-note" className="ml-auto font-mono text-caption text-ink-soft">
        band {pct(RESID_BAND)} · track {pct(RESID_DOMAIN)}
      </span>
    </div>
  );
}
