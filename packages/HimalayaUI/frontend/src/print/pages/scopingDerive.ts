import type { OrderingRow } from "../../lib/scoping/proposeOrdering";
import type { PhaseRead } from "../../lib/scoping/dominantPhase";
import type { PhaseSegment } from "../ui";
import type { PickerSampleRow } from "../../api";

/** Filter the corpus picker rows to only those whose sample.id is in `seed`.
 *  When `seed` is null the full corpus is returned (no-filter path for a
 *  direct /series/new visit). */
export function filterPickerBySeed(
  rows: PickerSampleRow[],
  seed: number[] | null,
): PickerSampleRow[] {
  if (seed === null) return rows;
  const ids = new Set(seed);
  return rows.filter((r) => ids.has(r.sample.id));
}

export interface FootState {
  kind: "warn" | "ready";
  text: string;
}

/** The confirm-gate state line. `flagged` is reframed as "skip this read from
 *  the write" (the honest Option-A model: every member's parsed read is
 *  committed unless the user skips it). `keptCount` = members whose read will be
 *  committed; `skippedCount` = members the user excluded. Warn only when nothing
 *  is kept (nothing to build); otherwise ready, annotating any skips. */
export function buildFootState(keptCount: number, skippedCount: number): FootState {
  if (keptCount === 0) {
    return { kind: "warn", text: "Keep at least one value to build" };
  }
  const base = `${keptCount} value${keptCount === 1 ? "" : "s"} ready to commit`;
  return { kind: "ready", text: skippedCount > 0 ? `${base} · ${skippedCount} skipped` : base };
}

/** Build gate: an ordering key must exist and at least one KEPT member must
 *  remain (included, not skipped, has a value). A skipped (flagged) member is
 *  excluded from the write, never blocks. */
export function canScopeBuild(rows: OrderingRow[], orderingKey: string | undefined): boolean {
  if (orderingKey === undefined) return false;
  return rows.some((r) => r.include && !r.flagged && r.value !== "");
}

/** The batch payload: only included, non-skipped members WITH a value are
 *  written. The empty-value guard prevents corrupting a sample with `value:""`
 *  (loose matches carry no value and must never reach the write). */
export function buildScopePayload(rows: OrderingRow[]): { sampleId: number; value: string }[] {
  return rows
    .filter((r) => r.include && !r.flagged && r.value !== "")
    .map((r) => ({ sampleId: r.sampleId, value: r.value }));
}

export interface ColdAssignRow {
  sampleId: number;
  sampleName: string;
  value: string;
}

/** Seed the cold-assign worksheet from a list of (id, name) pairs.
 *  All values start empty — the user fills them in. */
export function buildColdAssignRows(
  seed: ReadonlyArray<{ sampleId: number; sampleName: string }>,
): ColdAssignRow[] {
  return seed.map((s) => ({ sampleId: s.sampleId, sampleName: s.sampleName, value: "" }));
}

/** Cold-assign build gate: key must be non-empty and every sample must have a
 *  non-empty value. Never blocks on a partial fill — controls-don't-lie. */
export function canColdBuild(key: string, rows: ColdAssignRow[]): boolean {
  if (key.trim() === "") return false;
  return rows.every((r) => r.value.trim() !== "");
}

/** Cold-assign scope payload: all rows unconditionally (gate is upstream). */
export function buildColdScopePayload(rows: ColdAssignRow[]): { sampleId: number; value: string }[] {
  return rows.map((r) => ({ sampleId: r.sampleId, value: r.value }));
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
