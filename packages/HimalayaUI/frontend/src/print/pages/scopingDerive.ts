import type { OrderingRow } from "../../lib/scoping/proposeOrdering";
import type { PhaseRead } from "../../lib/scoping/dominantPhase";
import type { PhaseSegment } from "../ui";

export interface FootState {
  kind: "warn" | "ready";
  text: string;
}

/** The confirm-gate state line. Warn while any value still needs a look; ready
 *  once all are confirmed. `memberCount` is the member total (the denominator
 *  in the "All N values confirmed" copy). */
export function buildFootState(flagCount: number, memberCount: number): FootState {
  if (flagCount > 0) {
    return {
      kind: "warn",
      text: `${flagCount} value${flagCount === 1 ? "" : "s"} to check before you can build`,
    };
  }
  return { kind: "ready", text: `All ${memberCount} values confirmed — ready to build` };
}

/** Build gate (carried contract from the legacy page): an ordering key must
 *  exist, at least one member must be included, and no included member may be
 *  flagged. The batch route 400s on an empty array. */
export function canScopeBuild(rows: OrderingRow[], orderingKey: string | undefined): boolean {
  if (orderingKey === undefined) return false;
  const included = rows.filter((r) => r.include && !r.flagged);
  return included.length > 0 && rows.every((r) => !r.include || !r.flagged);
}

/** The batch payload: only included, non-flagged members are written. */
export function buildScopePayload(rows: OrderingRow[]): { sampleId: number; value: string }[] {
  return rows
    .filter((r) => r.include && !r.flagged)
    .map((r) => ({ sampleId: r.sampleId, value: r.value }));
}

/** Map per-member phase reads (dominant + optional coexist partner) onto the
 *  PhaseStrip preview segments. A null dominant stays a null (unindexed) cell;
 *  a coexist partner becomes the single-element `coexistWith` array (the
 *  greenfield strip takes `string[]`, not the legacy `string | null`). */
export function toPreviewSegments(reads: PhaseRead[]): PhaseSegment[] {
  return reads.map((r) =>
    r.coexist ? { phase: r.dominant, coexistWith: [r.coexist] } : { phase: r.dominant },
  );
}
