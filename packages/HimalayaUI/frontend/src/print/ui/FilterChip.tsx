import { Chip } from "./Chip";
import type { ChipSize } from "./Chip";

interface FilterChipProps {
  label: string;
  active: boolean;
  onClick: () => void;
  /** Size axis, forwarded to the base Chip. Defaults to `"md"` (text-sm). */
  size?: ChipSize;
  /** Optional native tooltip explaining what the filter selects (hover/focus).
   *  Omitted ⇒ no `title` (existing consumers stay byte-identical). */
  title?: string;
  className?: string;
}

/** Toggle pill (mockup `.chip`): a stateful filter that flips between resting
 *  and selected. Now a thin wrapper over the base `Chip` `toggle` variant, which
 *  owns the appearance (ink/paper inversion on active, focus-visible accent ring);
 *  FilterChip keeps its own `data-testid="filter-chip"` contract via `testId`.
 *  Distinct from FacetChip, a dropdown trigger (`aria-haspopup`) not a toggle. */
export function FilterChip({
  label,
  active,
  onClick,
  size = "md",
  title,
  className = "",
}: FilterChipProps): JSX.Element {
  return (
    <Chip
      variant="toggle"
      size={size}
      testId="filter-chip"
      active={active}
      onClick={onClick}
      {...(title !== undefined ? { title } : {})}
      className={className}
    >
      {label}
    </Chip>
  );
}
