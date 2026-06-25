import { cx } from "../../lib/cx";
/**
 * RailBack — floating tab that restores a collapsed editing rail.
 *
 * Presentational contract: NO local state. The consumer decides WHEN to mount
 * it (only when the rail is collapsed); RailBack itself is always visible when
 * rendered. Anchored to the right edge, vertically centered.
 *
 * The label is a vertical-text run (`[writing-mode:vertical-rl]` — an arbitrary
 * LAYOUT property, not appearance) in the uppercase/700/tracked `text-kicker`
 * role, sitting beside a `‹` restore glyph.
 *
 * Radius note: uses `rounded-l-lg` — left-only rounding mirrors the mockup's
 * `8px 0 0 8px` (the tab fuses to the right edge of the viewport).
 */


export interface RailBackProps {
  /** Vertical-text label. Default "Compose". */
  label?: string;
  onClick?: () => void;
  /** PLACEMENT-ONLY, appended last. */
  className?: string;
}

export function RailBack({
  label = "Compose",
  onClick,
  className,
}: RailBackProps): JSX.Element {
  return (
    <button
      type="button"
      data-testid="rail-back"
      onClick={onClick}
      className={cx(
        "fixed right-0 top-1/2 -translate-y-1/2 z-50 flex items-center gap-2",
        "bg-plate border border-hair-strong border-r-0 rounded-l-lg px-2 py-3.5 shadow-lg",
        "text-ink-soft hover:text-ink",
        className,
      )}
    >
      <span aria-hidden="true">&#8249;</span>
      <span className="text-kicker tracking-wider [writing-mode:vertical-rl]">{label}</span>
    </button>
  );
}
