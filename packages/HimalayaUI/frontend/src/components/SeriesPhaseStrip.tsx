/**
 * SeriesPhaseStrip — the per-sample phase strip + transition caption on a
 * folio card (R6, finding F-B; mockup `series-folio.html` `.ps` / `.ps-cap`).
 *
 * One cell PER sample, coloured by that member's confirmed phase, so a phase
 * transition reads at a glance across the whole wall. Coexistence renders as a
 * 2-stop gradient (wired, though a single snapshot carries one confirmed index
 * today — see `buildPhaseStrip`). The caption summarises the story:
 * "Pn3m → Lamellar", "Pn3m throughout", or "No clear phase".
 */
import type { SeriesMember } from "../api";
import { buildPhaseStrip } from "../lib/series/folioFigure";
import { phaseColor } from "../phases";

export interface SeriesPhaseStripProps {
  members: SeriesMember[];
}

const UNINDEXED = "var(--color-ink-faint)";

function segBackground(phase: string | null, coexistWith: string | null): string {
  if (phase === null) return UNINDEXED;
  if (coexistWith !== null) {
    return `linear-gradient(100deg, ${phaseColor(phase)} 42%, ${phaseColor(coexistWith)} 58%)`;
  }
  return phaseColor(phase);
}

export function SeriesPhaseStrip({ members }: SeriesPhaseStripProps): JSX.Element {
  const strip = buildPhaseStrip(members);

  return (
    <div data-testid="series-phase-strip">
      <div className="mt-3 flex h-[7px] gap-[2px]">
        {strip.segments.map((seg, i) => (
          <div
            key={i}
            data-testid="ps-seg"
            className="flex-1 rounded-[1.5px]"
            style={{ background: segBackground(seg.phase, seg.coexistWith) }}
          />
        ))}
      </div>
      <div
        data-testid="ps-cap"
        className="mt-1.5 flex items-center gap-1.5 text-[11px] text-ink-soft"
      >
        {strip.kind === "none" && (
          <span className="font-semibold text-ink-faint">No clear phase</span>
        )}
        {strip.kind === "throughout" && strip.first !== null && (
          <>
            <span className="font-semibold" style={{ color: phaseColor(strip.first) }}>
              {strip.first}
            </span>
            <span className="text-ink-faint">throughout</span>
          </>
        )}
        {strip.kind === "transition" && strip.first !== null && strip.last !== null && (
          <>
            <span className="font-semibold" style={{ color: phaseColor(strip.first) }}>
              {strip.first}
            </span>
            <span className="text-ink-faint" aria-hidden="true">
              →
            </span>
            <span className="font-semibold" style={{ color: phaseColor(strip.last) }}>
              {strip.last}
            </span>
          </>
        )}
      </div>
    </div>
  );
}
