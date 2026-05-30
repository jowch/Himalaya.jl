/**
 * SeriesPhaseStrip — the Series phase-strip companion (Plan E, Task E-5).
 *
 * Derives a `PhaseSegment[]` from the members' snapshots (in display / variable
 * order) and composes the `ui/PhaseStrip` primitive: one cell per sample,
 * coexistence → a two-phase gradient (`coexistWith`), form-factor → a hollow
 * dashed cell, null → a faint distinct cell. The 3-state comes from the
 * snapshot's `assignment_state` (E-7); the coexisting phase from
 * `confirmed_phases` (E-4). All appearance lives in the primitive — this
 * component is a pure deriver + placement wrapper.
 */
import type { SeriesMember, MemberSnapshotPhase } from "../api";
import { PhaseStrip } from "./ui/PhaseStrip";
import type { PhaseSegment } from "./ui/PhaseStrip";

/** The phases a member's assignment carries, honouring the 3-state (form-factor
 *  / null carry none) and falling back to confirmed_index for legacy snapshots. */
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

/** Derive one `PhaseSegment` per member, in member (variable) order. */
export function segmentsFromMembers(members: SeriesMember[]): PhaseSegment[] {
  return members.map((m) => {
    const state = m.snapshot?.assignment_state;
    if (state === "form_factor") return { phase: null, state: "form_factor" };
    if (state === "null") return { phase: null, state: "null" };
    const phases = phasesOf(m);
    if (phases.length === 0) return { phase: null };
    return {
      phase: phases[0]!.phase,
      coexistWith: phases.length > 1 ? phases[1]!.phase : null,
    };
  });
}

export interface SeriesPhaseStripProps {
  members: SeriesMember[];
  /** PLACEMENT ONLY. */
  className?: string;
}

export function SeriesPhaseStrip({ members, className }: SeriesPhaseStripProps): JSX.Element {
  const segments = segmentsFromMembers(members);
  return (
    <PhaseStrip
      segments={segments}
      emptyLabel="No phase calls yet"
      {...(className !== undefined ? { className } : {})}
    />
  );
}
