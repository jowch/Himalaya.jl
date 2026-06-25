import type { JSX } from "react";

/**
 * SourceChip — a small provenance badge for a geometry field's source
 * (prp / setup / computed / unset / edited). The user-edited state gets the
 * accent wash; every other source is a neutral paper-sunk badge.
 *
 * Shared by the gear Configuration (GeometryLedger) AND the first-run funnel
 * (ConfigurationPage) so the chips look identical in both places. Lives in
 * print/components/ (NOT the design-guard-exempt print/ui/ layer), so it may
 * only use named token utilities — bg-accent-wash is sanctioned because
 * --color-accent-wash is declared in styles.css @theme.
 */
export function SourceChip({
  label,
  emphasized = false,
}: {
  label: string;
  /** The field was user-edited — render the accent wash instead of neutral. */
  emphasized?: boolean;
}): JSX.Element {
  return (
    <span
      className={
        emphasized
          ? "rounded-sm bg-accent-wash px-1.5 py-0.5 text-xs font-bold uppercase text-accent"
          : "rounded-sm bg-paper-sunk px-1.5 py-0.5 text-xs font-bold uppercase text-ink-faint"
      }
    >
      {label}
    </span>
  );
}
