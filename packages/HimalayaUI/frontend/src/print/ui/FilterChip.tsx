interface FilterChipProps {
  label: string;
  active: boolean;
  onClick: () => void;
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** Toggle pill (mockup `.chip`): a stateful filter that flips between resting
 *  and selected. ON inverts the surface to ink-on-paper (`bg-ink text-paper
 *  border-ink`) rather than shifting hue — the second channel is the fill/
 *  contrast inversion, so the selected state reads in grayscale (A). Distinct
 *  from FacetChip, which is a dropdown trigger (`aria-haspopup`), not a toggle.
 *
 *  B — no rationed accent on resting/active chrome; terracotta appears only on
 *  the focus-visible ring. C — ships default / hover (border→`ink-faint`) /
 *  focus-visible / active (the inverted ON state); disabled N/A for a filter
 *  toggle. D — `px-3 py-1` compact pill; effective hit area is enlarged by its
 *  facet/filter-bar container (known-dense chrome control). */
export function FilterChip({
  label,
  active,
  onClick,
  className = "",
}: FilterChipProps): JSX.Element {
  return (
    <button
      type="button"
      data-testid="filter-chip"
      aria-pressed={active}
      data-active={active ? "true" : "false"}
      onClick={onClick}
      className={cx(
        "text-sm font-semibold rounded-full px-3 py-1 border transition-colors",
        "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
        active
          ? "bg-ink text-paper border-ink"
          : "bg-plate text-ink-soft border-hair-strong hover:border-ink-faint",
        className,
      )}
    >
      {label}
    </button>
  );
}
