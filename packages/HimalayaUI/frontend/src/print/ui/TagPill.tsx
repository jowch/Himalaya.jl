import type { ReactNode } from "react";

interface TagPillProps {
  children: ReactNode;
  onRemove?: () => void;
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** A single tag pill (mockup `.tag`): neutral chrome — `bg-plate` fill, hairline
 *  border, `rounded-full`, `text-xs font-semibold text-ink-soft`. Carries no
 *  color-coded meaning, so A/B are otherwise N/A for the resting pill.
 *
 *  B — the rationed grease-pencil accent appears ONLY on the optional remove (×)
 *  affordance's hover (→accent) and focus-visible (accent outline), never on the
 *  resting pill. C — the × ships hover + focus-visible; the pill itself is static.
 *  H — flat (`bg-plate` is the allowed surface, no shadow/gradient). */
export function TagPill({
  children,
  onRemove,
  className = "",
}: TagPillProps): JSX.Element {
  return (
    <span
      data-testid="tag-pill"
      className={cx(
        "inline-flex items-center gap-1 text-xs font-semibold text-ink-soft bg-plate border border-hair rounded-full px-2 py-px whitespace-nowrap",
        className,
      )}
    >
      {children}
      {onRemove && (
        <button
          type="button"
          aria-label="Remove tag"
          onClick={onRemove}
          className="text-ink-faint hover:text-accent focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent rounded-sm transition-colors"
        >
          ×
        </button>
      )}
    </span>
  );
}
