/**
 * SeriesReadingPanel — the rail's derived "phases present" readout (Plan E,
 * Task E-4). Renders the `SeriesReading` shape: one row per indexed phase
 * (phase-coloured name + variable span + lattice trend), plus coexistence and
 * form-factor lines. Pure presentation over the derivation — no data fetching.
 *
 * Phase names/dots carry `phaseColor` via inline `style` (a resolved colour
 * threaded through a style attr, NOT an appearance className) so coexistence
 * rows self-decode — the same idiom as PhaseStrip / MemberMetaRow. Placement
 * lives on the consumer's `className`.
 */
import type { SeriesReading } from "../lib/series/seriesReading";
import { phaseColor } from "../phases";

export interface SeriesReadingPanelProps {
  reading: SeriesReading;
  /** PLACEMENT ONLY: margin / width / grid. */
  className?: string;
}

export function SeriesReadingPanel({ reading, className = "" }: SeriesReadingPanelProps): JSX.Element {
  const { phases, coexistenceAt, formFactorOnlyAt } = reading;
  return (
    <div
      data-testid="series-reading"
      className={`flex flex-col gap-2 rounded-md border border-hair bg-plate px-3 py-3 ${className}`}
    >
      {phases.map((p) => (
        <div key={p.phase} data-testid="reading-phase-row" className="flex flex-col gap-0.5">
          <div className="flex items-center gap-2">
            <span
              aria-hidden="true"
              className="inline-block h-2 w-2 shrink-0 rounded-full"
              style={{ background: phaseColor(p.phase) }}
            />
            <span
              data-testid="reading-phase-name"
              className="text-sm font-bold"
              style={{ color: phaseColor(p.phase) }}
            >
              {p.phase}
            </span>
            <span className="ml-auto font-mono text-xs text-ink-soft tabular-nums">
              {p.spanLabel}
            </span>
          </div>
          <div className="ml-4 font-mono text-xs text-ink-faint tabular-nums">
            {p.latticeTrend}
          </div>
        </div>
      ))}
      {coexistenceAt.length > 0 && (
        <div
          data-testid="reading-coexistence"
          className="border-t border-hair pt-2 font-mono text-xs text-ink-soft"
        >
          coexistence at {coexistenceAt.join(", ")}
        </div>
      )}
      {formFactorOnlyAt.length > 0 && (
        <div
          data-testid="reading-form-factor"
          className="border-t border-hair pt-2 font-mono text-xs text-ink-soft"
        >
          form factor only at {formFactorOnlyAt.join(", ")}
        </div>
      )}
    </div>
  );
}
