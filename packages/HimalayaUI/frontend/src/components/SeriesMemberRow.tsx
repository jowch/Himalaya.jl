/**
 * SeriesMemberRow — a Series rail member row (Plan E, Task E-6).
 *
 * Per member, derived from the snapshot: the lattice (`a` for cubics / `d` for
 * lamellar; BOTH shown under coexistence, e.g. `a 195 · d 60 Å`), the first-peak
 * q₁ = `min(effective_peaks.q)`, and the phase name(s) `phaseColor`-coloured so
 * coexistence rows self-decode. A form-factor / null member shows the
 * "no Bragg peaks · q₁ —" line and a neutral name.
 *
 * Phase colour is threaded through inline `style` (a resolved colour on a style
 * attr, not an appearance className) — the same idiom as PhaseStrip /
 * SeriesReadingPanel.
 */
import type { SeriesMember, MemberSnapshotPhase } from "../api";
import { phaseColor, CUBIC_PHASES } from "../phases";

function phasesOf(m: SeriesMember): MemberSnapshotPhase[] {
  const snap = m.snapshot;
  if (!snap) return [];
  if (snap.assignment_state === "null" || snap.assignment_state === "form_factor") return [];
  if (snap.confirmed_phases && snap.confirmed_phases.length > 0) return snap.confirmed_phases;
  if (snap.confirmed_index) {
    return [{ phase: snap.confirmed_index.phase, lattice_d: snap.confirmed_index.lattice_d }];
  }
  return [];
}

/** q₁ = the smallest observed effective-peak q (the first Bragg reflection). */
function firstPeakQ(m: SeriesMember): number | null {
  const peaks = m.snapshot?.effective_peaks ?? [];
  if (peaks.length === 0) return null;
  let min = Infinity;
  for (const p of peaks) if (p.q < min) min = p.q;
  return Number.isFinite(min) ? min : null;
}

function latticeLabel(phases: MemberSnapshotPhase[]): string {
  const sym = (phase: string) => (CUBIC_PHASES.has(phase) ? "a" : "d");
  const fmt = (v: number | null) => (v == null ? "—" : (Number.isInteger(v) ? String(v) : v.toFixed(1)));
  if (phases.length > 1) {
    // Both lattices under coexistence: `a 195 · d 60 Å`.
    return phases.map((p) => `${sym(p.phase)} ${fmt(p.lattice_d)}`).join(" · ") + " Å";
  }
  const p = phases[0]!;
  return `${sym(p.phase)} = ${fmt(p.lattice_d)} Å`;
}

export interface SeriesMemberRowProps {
  member: SeriesMember;
  /** The member's ordering-variable value (e.g. "1:0.5"). */
  variableLabel: string;
  /** PLACEMENT ONLY. */
  className?: string;
  /** Optional hover/focus callbacks (wired by the page to drive band-hot). */
  onMouseEnter?: () => void;
  onMouseLeave?: () => void;
}

export function SeriesMemberRow({
  member, variableLabel, className = "", onMouseEnter, onMouseLeave,
}: SeriesMemberRowProps): JSX.Element {
  const phases = phasesOf(member);
  const q1 = firstPeakQ(member);
  const isFormFactor = member.snapshot?.assignment_state === "form_factor";
  const isBragg = phases.length > 0;

  return (
    <div
      data-testid="series-member-row"
      data-member-id={String(member.id)}
      className={`flex items-center gap-2 rounded-md px-2 py-1.5 hover:bg-plate ${className}`}
      onMouseEnter={onMouseEnter}
      onMouseLeave={onMouseLeave}
    >
      {/* swatch — phase colour, coexistence gradient, or hollow dashed (no phase) */}
      <span
        aria-hidden="true"
        className={isBragg ? "h-3 w-3 shrink-0 rounded-full" : "h-3 w-3 shrink-0 rounded-full border border-dashed border-hair-strong"}
        style={
          isBragg
            ? {
                background: phases.length > 1
                  ? `linear-gradient(135deg, ${phaseColor(phases[0]!.phase)} 48%, ${phaseColor(phases[1]!.phase)} 52%)`
                  : phaseColor(phases[0]!.phase),
              }
            : undefined
        }
      />
      <div className="flex min-w-0 flex-1 flex-col gap-px">
        <div className="flex items-baseline justify-between gap-2">
          <span data-testid="member-row-name" className="inline-flex items-baseline">
            {isBragg ? (
              phases.map((p, i) => (
                <span key={p.phase} className="inline-flex items-baseline">
                  {i > 0 && <span className="mx-0.5 text-xs text-ink-faint">+</span>}
                  <span className="text-sm font-semibold" style={{ color: phaseColor(p.phase) }}>
                    {p.phase}
                  </span>
                </span>
              ))
            ) : (
              <span className="text-sm font-semibold text-ink-faint">
                {isFormFactor ? "Form factor" : "No phase"}
              </span>
            )}
          </span>
          <span className="font-mono text-xs text-ink-soft tabular-nums">{variableLabel}</span>
        </div>
        <div data-testid="member-row-data" className="font-mono text-xs text-ink-faint tabular-nums">
          {isBragg
            ? `${latticeLabel(phases)} · q₁ ${q1 != null ? q1.toFixed(3) : "—"} Å⁻¹`
            : "no Bragg peaks · q₁ —"}
        </div>
      </div>
    </div>
  );
}
