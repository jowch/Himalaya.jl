import type { ReactNode } from "react";
import { IconButton } from "./IconButton";
import { cx } from "../../lib/cx";

const GLYPHS = {
  vertical: ["↑", "↓"],
  horizontal: ["←", "→"],
} as const;

export interface DockStepperProps {
  /** Axis name; drives the faint label and the `Previous/Next {label}` aria. */
  label: string;
  /** "vertical" → ↑/↓ (sample axis); "horizontal" → ←/→ (frame axis). */
  axis: keyof typeof GLYPHS;
  /** Stable id stem: `dock-prev-{base}` / `dock-next-{base}` / `dock-{base}-count`. */
  testIdBase: string;
  onPrev: () => void;
  onNext: () => void;
  prevDisabled?: boolean;
  nextDisabled?: boolean;
  /** Current/total readout (e.g. `3 / 12`). Omit to hide the readout. */
  count?: ReactNode;
  /** Readout min-width (placement). Default fits a two-digit `n / nn`. */
  countWidthClass?: string;
}

/**
 * Labeled prev/next stepper inside the bottom Dock (§7). Six-ish copies across
 * the Focus / Loupe / Series-builder docks collapsed to one primitive; callers
 * vary the label, axis glyphs, disable conditions, and the readout node, while
 * the `data-testid` contract (`dock-{prev,next}-*`, `dock-*-count`) stays fixed.
 */
export function DockStepper({
  label,
  axis,
  testIdBase,
  onPrev,
  onNext,
  prevDisabled,
  nextDisabled,
  count,
  countWidthClass = "min-w-[3.5rem]",
}: DockStepperProps): JSX.Element {
  const [prevGlyph, nextGlyph] = GLYPHS[axis];
  const verb = label.toLowerCase();
  return (
    <div className="flex items-center gap-1">
      <span className="text-meta text-ink-soft">{label}</span>
      <IconButton
        label={`Previous ${verb}`}
        tone="ghost"
        disabled={prevDisabled}
        onClick={onPrev}
        data-testid={`dock-prev-${testIdBase}`}
      >
        {prevGlyph}
      </IconButton>
      {count !== undefined && (
        <span
          className={cx("text-data tabular-nums text-ink text-center", countWidthClass)}
          data-testid={`dock-${testIdBase}-count`}
        >
          {count}
        </span>
      )}
      <IconButton
        label={`Next ${verb}`}
        tone="ghost"
        disabled={nextDisabled}
        onClick={onNext}
        data-testid={`dock-next-${testIdBase}`}
      >
        {nextGlyph}
      </IconButton>
    </div>
  );
}
