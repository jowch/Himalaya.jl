interface FacetChipProps {
  label: string;
  onClick?: () => void;
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** Pill-shaped dropdown trigger: a labeled chip that opens a menu (declared via
 *  `aria-haspopup="menu"`; the menu itself is the consumer's concern). Neutral
 *  chrome — no color-coded meaning (A/B N/A), so the `rounded-full` is the standard
 *  filter/dropdown-pill affordance, not decorative rounding. Distinct from the
 *  (later) FilterChip, which is a toggle rather than a trigger.
 *
 *  C — ships default / hover (`bg-paper-sunk`) / focus-visible (`accent` outline);
 *  active/disabled N/A for a transient menu trigger. D — `px-3 py-1` compact pill;
 *  effective hit area is enlarged by its toolbar/facet-bar container (known-dense
 *  chrome control). */
export function FacetChip({
  label,
  onClick,
  className = "",
}: FacetChipProps): JSX.Element {
  return (
    <button
      type="button"
      data-testid="facet-chip"
      onClick={onClick}
      aria-haspopup="menu"
      className={cx(
        "inline-flex items-center gap-2 border border-hair-strong bg-plate rounded-full px-3 py-1 text-sm font-semibold text-ink transition-colors",
        "hover:bg-paper-sunk focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
        className,
      )}
    >
      {label}
      <span aria-hidden="true" className="text-ink-faint">
        ▾
      </span>
    </button>
  );
}
