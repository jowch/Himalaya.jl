// memberRead.ts — the single shared read of a SeriesMember's rendered
// identity: dominant phase, assignment state, row label, and the annotated
// peak set.
//
// CONSTRAINT (BU-EXPORTDIVERGE): the builder's plate footnote promises
// "what you compose is what you publish", so the plate
// (print/waterfall/waterfallModel) and the export pipeline
// (lib/figure-export + lib/comparison/coloring byPhase) must resolve these
// four reads through the SAME helpers. Re-deriving any fallback chain in one
// layer is a WYSIWYG bug, not a style choice — the coexistence case
// (confirmed_phases[0] disagreeing with confirmed_index.phase) is exactly
// where duplicated chains drift apart.
import type {
  AssignmentState,
  MemberSnapshotPeak,
  SeriesMember,
} from "../../api";

/** Dominant phase: first `confirmed_phases` entry, else
 *  `confirmed_index.phase`, else null. Form-factor / null assignment states
 *  read as phaseless even when a lingering `confirmed_index` is present. */
export function dominantPhase(member: SeriesMember): string | null {
  const snap = member.snapshot;
  if (snap === null) return null;
  if (snap.assignment_state === "form_factor" || snap.assignment_state === "null") return null;
  const cp = snap.confirmed_phases;
  if (cp && cp.length > 0) return cp[0]!.phase;
  return snap.confirmed_index?.phase ?? null;
}

/** Resolve the assignment state, treating a missing snapshot or an undefined
 *  state (pre-E-7 snapshots) as "indexed". */
export function resolveState(member: SeriesMember): AssignmentState {
  const snap = member.snapshot;
  if (snap === null) return "indexed";
  return snap.assignment_state ?? "indexed";
}

/** The plate's row-label register: `label_override`, else
 *  `"exp <exposure_id>"`, else `""`. The export must speak this register too
 *  — never the retired `"Exposure #N"` long form. */
export function memberRowLabel(member: SeriesMember): string {
  if (member.label_override !== null) return member.label_override;
  return member.exposure_id != null ? `exp ${member.exposure_id}` : "";
}

/** The peak set the plate annotates: the confirmed index's anchor peaks only,
 *  ascending q (the plate's 1..n numbering register). Unindexed members and
 *  form-factor / null members have no anchors — NOT all `effective_peaks`. */
export function indexedAnchorPeaks(member: SeriesMember): MemberSnapshotPeak[] {
  const snap = member.snapshot;
  if (snap === null || resolveState(member) !== "indexed") return [];
  const confirmedIndex = snap.confirmed_index;
  if (confirmedIndex === null) return [];
  const peakById = new Map(snap.effective_peaks.map((p) => [p.id, p]));
  // Defensive: legacy snapshot rows can carry a confirmed_index without
  // peak_ids despite the typed contract — read as "no anchors", not a crash.
  return (confirmedIndex.peak_ids ?? [])
    .filter((id) => peakById.has(id))
    .map((id) => peakById.get(id)!)
    .sort((a, b) => a.q - b.q);
}
