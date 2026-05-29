import { phaseColor } from "../phases";
import type { PhaseRead } from "../lib/scoping/dominantPhase";

interface Props {
  reads: PhaseRead[]; // one per member, in display (low→high) order
}

function segBackground(r: PhaseRead): string {
  if (r.dominant === null) return "var(--color-hair)";
  if (r.coexist) {
    return `linear-gradient(100deg, ${phaseColor(r.dominant)} 42%, ${phaseColor(r.coexist)} 58%)`;
  }
  return phaseColor(r.dominant);
}

/**
 * Preview phase strip (series-scoping.html `.preview`): the phase story this
 * order will print — one cell per member (two-phase coexistence = gradient),
 * captioned first→last (or "throughout" when they agree).
 */
export function ScopingPhaseStrip({ reads }: Props): JSX.Element {
  const first = reads.find((r) => r.dominant !== null)?.dominant ?? null;
  const last = [...reads].reverse().find((r) => r.dominant !== null)?.dominant ?? null;
  return (
    <div className="mt-5">
      <div className="mb-1.5 text-[10.5px] font-bold uppercase tracking-wider text-ink-faint">
        Preview: phase across the series
      </div>
      <div className="flex h-2 gap-0.5">
        {reads.map((r, i) => (
          <div
            key={i}
            data-testid="scoping-ps-seg"
            className="flex-1 rounded-[1.5px]"
            style={{ background: segBackground(r) }}
          />
        ))}
      </div>
      <div data-testid="scoping-ps-cap" className="mt-2 flex items-center gap-1.5 text-[11.5px] text-ink-soft">
        {first === null ? (
          <span className="text-ink-faint">Members not yet indexed; phase preview unavailable.</span>
        ) : first === last ? (
          <>
            <span className="font-semibold" style={{ color: phaseColor(first) }}>{first}</span>
            <span className="text-ink-faint">throughout</span>
          </>
        ) : (
          <>
            <span className="font-semibold" style={{ color: phaseColor(first) }}>{first}</span>
            <span className="text-ink-faint">→</span>
            <span className="font-semibold" style={{ color: phaseColor(last!) }}>{last}</span>
          </>
        )}
      </div>
    </div>
  );
}
