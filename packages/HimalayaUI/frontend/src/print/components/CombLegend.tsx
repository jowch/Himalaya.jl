import { PeakGlyphFromOpts } from "../ui/peakMark";
import type { PeakMarkOpts } from "../ui/peakMark";

export type GlyphVocabItem = "auto" | "manual" | "predicted-absent" | "excluded";

interface CombLegendProps {
  /** Which vocabulary entries to show, in order. Default: all four. */
  items?: ReadonlyArray<GlyphVocabItem>;
  /** PLACEMENT-ONLY. */
  className?: string;
}

const ALL_ITEMS: ReadonlyArray<GlyphVocabItem> = [
  "auto",
  "manual",
  "predicted-absent",
  "excluded",
];

const LABELS: Record<GlyphVocabItem, string> = {
  auto: "observed",
  manual: "manual",
  "predicted-absent": "predicted, absent",
  excluded: "excluded",
};

const NEUTRAL_COLOR = "var(--color-ink-soft)";

const SWATCH_R = 4;

function glyphOpts(item: GlyphVocabItem): PeakMarkOpts {
  switch (item) {
    case "auto":
      return { source: "auto", color: NEUTRAL_COLOR, r: SWATCH_R };
    case "manual":
      return { source: "manual", color: NEUTRAL_COLOR, r: SWATCH_R };
    case "predicted-absent":
      return { source: "auto", predictedAbsent: true, color: NEUTRAL_COLOR, r: SWATCH_R };
    case "excluded":
      return { source: "auto", excluded: true, color: NEUTRAL_COLOR, r: SWATCH_R };
  }
}

export function CombLegend({ items, className }: CombLegendProps): JSX.Element {
  const shown = items ?? ALL_ITEMS;
  return (
    <div
      className={`flex flex-wrap items-center gap-3 border-t border-hair pt-2${className ? ` ${className}` : ""}`}
    >
      {shown.map((item) => (
        <span key={item} className="inline-flex items-center gap-1.5">
          <svg
            width="18"
            height="14"
            viewBox="0 0 18 14"
            aria-hidden="true"
          >
            <PeakGlyphFromOpts opts={glyphOpts(item)} x={9} y={11} />
          </svg>
          <span className="font-mono text-caption text-ink-faint">
            {LABELS[item]}
          </span>
        </span>
      ))}
    </div>
  );
}
