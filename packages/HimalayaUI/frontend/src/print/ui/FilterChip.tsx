import { Chip } from "./Chip";

interface FilterChipProps {
  label: string;
  active: boolean;
  onClick: () => void;
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
  className = "",
}: FilterChipProps): JSX.Element {
  return (
    <Chip
      variant="toggle"
      testId="filter-chip"
      active={active}
      onClick={onClick}
      className={className}
    >
      {label}
    </Chip>
  );
}
