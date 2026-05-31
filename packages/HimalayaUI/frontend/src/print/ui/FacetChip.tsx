import { Chip } from "./Chip";

interface FacetChipProps {
  label: string;
  onClick?: () => void;
  className?: string;
}

/** Pill-shaped dropdown trigger: a labeled chip that opens a menu. Now a thin
 *  wrapper over the base `Chip` `trigger` variant, which owns the appearance plus
 *  the `aria-haspopup="menu"` declaration and the aria-hidden ▾ chevron; FacetChip
 *  keeps its own `data-testid="facet-chip"` contract via `testId`. The menu itself
 *  stays the consumer's concern (Menu-wiring deferred — trigger-only). */
export function FacetChip({
  label,
  onClick,
  className = "",
}: FacetChipProps): JSX.Element {
  return (
    <Chip
      variant="trigger"
      testId="facet-chip"
      {...(onClick ? { onClick } : {})}
      className={className}
    >
      {label}
    </Chip>
  );
}
