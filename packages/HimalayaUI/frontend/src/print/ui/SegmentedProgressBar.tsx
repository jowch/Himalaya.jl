import { cx } from "../../lib/cx";
import type { ProgressSegment } from "../../lib/ingestStages";

interface SegmentedProgressBarProps {
  segments: ProgressSegment[];
  /** Accessible name for the whole strip (e.g. "Scan progress"). */
  label?: string;
  className?: string;
}

/** Multi-segment sibling of `ProgressBar`: one equal-width capsule per stage,
 *  each filling against its OWN denominator, so the bar advances left-to-right
 *  through stages instead of resetting a shared scale between them.
 *
 *  Same visual language as `ProgressBar` — 4px `rounded-full bg-hair` rails with
 *  `bg-accent` fills — so the two read as one system. Progress is carried by fill
 *  WIDTH plus segment POSITION, never hue alone, so it survives grayscale (A);
 *  terracotta `accent` stays the single rationed progress mark (B). The completed
 *  stages hold a dimmed `accent` so "done" is distinguishable from "in progress"
 *  without introducing a second hue.
 *
 *  Non-interactive indicator (C/D N/A). Exposes one `role="progressbar"` for the
 *  whole strip, reporting the ACTIVE segment's counts — a screen reader gets
 *  "Analyzing, 92 of 604" rather than an unlabelled aggregate percentage that
 *  would have to invent cross-stage weights. H — flat (rail tonal step + fill). */
export function SegmentedProgressBar({
  segments,
  label,
  className = "",
}: SegmentedProgressBarProps): JSX.Element {
  const active = segments.find((s) => s.active);
  const activeLabel = active ? `${label ? `${label}: ` : ""}${active.label}` : label;

  return (
    <div
      data-testid="segmented-progressbar"
      role="progressbar"
      // ARIA rides the FRACTION, not the raw counts: a zero-total stage renders
      // full (nothing to do) but would report valuemin == valuemax == 0, which is
      // invalid ARIA with undefined AT behavior -- and would contradict the bar.
      // Percent keeps min/max fixed and always distinct.
      aria-valuenow={active ? Math.round(active.fraction * 100) : 0}
      aria-valuemin={0}
      aria-valuemax={100}
      aria-valuetext={
        active
          ? (active.total && active.total > 0
              ? `${active.label}, ${active.processed ?? 0} of ${active.total}`
              : `${active.label}, complete`)
          : undefined
      }
      aria-label={activeLabel}
      className={cx("flex items-center gap-1", className)}
    >
      {segments.map((s) => (
        <span
          key={s.key}
          data-testid={`segment-${s.key}`}
          data-active={s.active ? "true" : undefined}
          data-fraction={s.fraction}
          // flex-1 == equal width; see ingestStages.ts for why not weighted.
          className="block h-1 flex-1 overflow-hidden rounded-full bg-hair"
        >
          <span
            className={cx("block h-full", s.active ? "bg-accent" : "bg-accent/45")}
            style={{ width: `${s.fraction * 100}%` }}
          />
        </span>
      ))}
    </div>
  );
}
