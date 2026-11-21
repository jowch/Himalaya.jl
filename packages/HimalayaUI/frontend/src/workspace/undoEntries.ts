import type { Exposure, IndexEntry, Peak, SampleBridge, SampleLabelState, LabelColor } from "../api";
import type { EditOps } from "../queries";
import { bridgedWith } from "../lib/bridges";
import { sameLabel } from "../lib/labels";
import { peakQTol } from "../lib/queue/peakQTol";
import type { UndoEntry } from "./undo";

// physics: a re-fit after a peak edit moves a lattice's basis slightly; the backend re-attaches calls within 5% (MEMBER_REATTACH_RELTOL).
export const PHASE_BASIS_RELTOL = 0.05;

export interface PhaseIdentity {
  exposureId: number;
  phase: string;
  basis: number;
  kind: IndexEntry["kind"];
}

interface Where {
  label: string;
  sampleId: number;
  exposureId: number | undefined;
}

export function matchIndex(
  indices: readonly IndexEntry[],
  phase: string,
  basis: number,
  keep: (ix: IndexEntry) => boolean = () => true,
): IndexEntry | undefined {
  const tol = Math.max(PHASE_BASIS_RELTOL * Math.abs(basis), 1e-9);
  let best: IndexEntry | undefined;
  for (const ix of indices) {
    if (ix.id <= 0 || ix.phase !== phase || !keep(ix)) continue;
    const d = Math.abs(ix.basis - basis);
    if (d <= tol && (!best || d < Math.abs(best.basis - basis))) best = ix;
  }
  return best;
}

export async function identify(
  ops: EditOps,
  refs: readonly { exposureId: number; indexId: number }[],
): Promise<PhaseIdentity[]> {
  const out: PhaseIdentity[] = [];
  for (const ref of refs) {
    const ix = (await ops.read.knownIndices(ref.exposureId)).find((i) => i.id === ref.indexId);
    if (ix) out.push({ exposureId: ref.exposureId, phase: ix.phase, basis: ix.basis, kind: ix.kind });
  }
  return out;
}

export async function setPhaseTicked(ops: EditOps, id: PhaseIdentity, ticked: boolean): Promise<number> {
  const [indices, assignment] = await Promise.all([ops.read.indices(id.exposureId), ops.read.assignment(id.exposureId)]);
  const members = new Set(assignment.members);
  const hit = matchIndex(indices, id.phase, id.basis, (ix) => ix.kind === id.kind && members.has(ix.id) !== ticked);
  if (!hit) return 0;
  if (ticked) ops.addPhase(id.exposureId, hit.id);
  else ops.removePhase(id.exposureId, hit.id);
  return 1;
}

async function setEach(ops: EditOps, ids: readonly PhaseIdentity[], ticked: boolean): Promise<number> {
  let sent = 0;
  for (const id of ids) sent += await setPhaseTicked(ops, id, ticked);
  return sent;
}

async function restoreFormFactor(ops: EditOps, sampleId: number, exposureId: number): Promise<number> {
  const a = await ops.read.assignment(exposureId);
  if (a.state === "form_factor") return 0;
  ops.setCallState(sampleId, "form_factor");
  return 1;
}

export function tickEntry(ops: EditOps, at: Where & { id: PhaseIdentity; wasFormFactor: boolean }): UndoEntry {
  return {
    ...at,
    undo: async () => {
      const sent = await setPhaseTicked(ops, at.id, false);
      return sent + (at.wasFormFactor ? await restoreFormFactor(ops, at.sampleId, at.id.exposureId) : 0);
    },
    redo: () => setPhaseTicked(ops, at.id, true),
  };
}

export function untickEntry(ops: EditOps, at: Where & { ids: Promise<PhaseIdentity[]> }): UndoEntry {
  return {
    ...at,
    undo: async () => setEach(ops, await at.ids, true),
    redo: async () => setEach(ops, await at.ids, false),
  };
}

export function noPhaseEntry(
  ops: EditOps,
  at: Where & { exposureId: number; on: boolean; ticked: Promise<PhaseIdentity[]> },
): UndoEntry {
  const set = async (formFactor: boolean): Promise<number> => {
    const a = await ops.read.assignment(at.exposureId);
    if ((a.state === "form_factor") === formFactor) return 0;
    ops.setCallState(at.sampleId, formFactor ? "form_factor" : "null");
    return 1;
  };
  return {
    ...at,
    undo: async () => {
      const sent = await set(!at.on);
      return at.on && sent > 0 ? sent + (await setEach(ops, await at.ticked, true)) : sent;
    },
    redo: () => set(at.on),
  };
}

export function dropEntry(
  ops: EditOps,
  at: Where & { exposureId: number; prev: Exposure["status"]; next: Exposure["status"] },
): UndoEntry {
  return {
    ...at,
    undo: async () => {
      ops.setStatus(at.sampleId, at.exposureId, at.prev);
      return 1;
    },
    redo: async () => {
      ops.setStatus(at.sampleId, at.exposureId, at.next);
      return 1;
    },
  };
}

function near(peaks: readonly Peak[], q: number, auto: boolean): Peak | undefined {
  return peaks.find((p) => p.id > 0 && (p.source === "auto") === auto && Math.abs(p.q - q) < peakQTol(q));
}

async function setManualPeak(ops: EditOps, exposureId: number, q: number, present: boolean): Promise<number> {
  const hit = near(await ops.read.peaks(exposureId), q, false);
  if (present === Boolean(hit)) return 0;
  if (hit) ops.removePeak(exposureId, hit.id);
  else ops.addPeak(exposureId, q);
  return 1;
}

export function peakEntry(ops: EditOps, at: Where & { exposureId: number; q: number; added: boolean }): UndoEntry {
  return {
    ...at,
    undo: () => setManualPeak(ops, at.exposureId, at.q, !at.added),
    redo: () => setManualPeak(ops, at.exposureId, at.q, at.added),
  };
}

export function excludeEntry(ops: EditOps, at: Where & { exposureId: number; q: number; excluded: boolean }): UndoEntry {
  const set = async (excluded: boolean): Promise<number> => {
    const hit = near(await ops.read.peaks(at.exposureId), at.q, true);
    if (!hit || hit.excluded === excluded) return 0;
    ops.setPeakExcluded(at.exposureId, hit.id, excluded);
    return 1;
  };
  return { ...at, undo: () => set(!at.excluded), redo: () => set(at.excluded) };
}

export function customIndexEntry(
  ops: EditOps,
  at: Where & { exposureId: number; phase: string; basis: number; ratios: number[]; wasFormFactor: boolean },
): UndoEntry {
  const find = async () =>
    matchIndex(await ops.read.indices(at.exposureId), at.phase, at.basis, (ix) => ix.kind === "speculative");
  return {
    ...at,
    undo: async () => {
      const hit = await find();
      if (!hit) return 0;
      ops.deleteIndex(at.exposureId, hit.id);
      return 1 + (at.wasFormFactor ? await restoreFormFactor(ops, at.sampleId, at.exposureId) : 0);
    },
    redo: async () => {
      if (await find()) return 0;
      ops.commitCustom(at.exposureId, at.phase, at.basis, at.ratios);
      return 1;
    },
  };
}

export function labelEntry(
  ops: EditOps,
  at: Where & { tag: string; added: boolean; prev: SampleLabelState | undefined },
): UndoEntry {
  const current = async () =>
    (await ops.read.labels()).find((l) => l.sample_id === at.sampleId && sameLabel(l.label, at.tag))?.state;
  const set = async (present: boolean): Promise<number> => {
    if ((await current()) === (present ? "confirmed" : undefined)) return 0;
    if (present) ops.addLabel(at.sampleId, at.tag);
    else ops.removeLabel(at.sampleId, at.tag);
    return 1;
  };
  const resuggest = async (): Promise<number> => {
    if ((await current()) === "suggested") return 0;
    ops.restoreLabelSuggestion(at.sampleId, at.tag);
    return 1;
  };
  return {
    ...at,
    undo: () => (at.prev === "suggested" ? resuggest() : set(!at.added)),
    redo: () => set(at.added),
  };
}

export function combinedEntry(at: Where, parts: readonly UndoEntry[]): UndoEntry {
  return {
    ...at,
    undo: async () => {
      let sent = 0;
      for (const p of [...parts].reverse()) sent += await p.undo();
      return sent;
    },
    redo: async () => {
      let sent = 0;
      for (const p of parts) sent += await p.redo();
      return sent;
    },
  };
}

export function labelColorEntry(
  ops: EditOps,
  at: Where & { tag: string; prev: LabelColor | null; next: LabelColor | null },
): UndoEntry {
  return {
    ...at,
    undo: async () => {
      ops.setLabelColor(at.tag, at.prev);
      return 1;
    },
    redo: async () => {
      ops.setLabelColor(at.tag, at.next);
      return 1;
    },
  };
}

function together(list: readonly SampleBridge[], a: number, b: number): boolean {
  return bridgedWith(list, a).includes(b);
}

export function bridgeEntry(
  ops: EditOps,
  at: Where & { from: number; to: number; bridges: readonly SampleBridge[] },
): UndoEntry {
  const sideA = [at.from, ...bridgedWith(at.bridges, at.from)];
  const sideB = [at.to, ...bridgedWith(at.bridges, at.to)];
  const split = sideA.length <= sideB.length ? sideA : sideB;
  return {
    ...at,
    undo: async () => {
      if (!together(await ops.read.bridges(), at.from, at.to)) return 0;
      for (const s of split) ops.unbridge(s);
      for (const s of split.slice(1)) ops.bridge(split[0]!, s);
      return split.length * 2 - 1;
    },
    redo: async () => {
      if (together(await ops.read.bridges(), at.from, at.to)) return 0;
      ops.bridge(at.from, at.to);
      return 1;
    },
  };
}

export function unbridgeEntry(
  ops: EditOps,
  at: Where & { target: number; bridges: readonly SampleBridge[] },
): UndoEntry {
  const others = bridgedWith(at.bridges, at.target);
  return {
    ...at,
    undo: async () => {
      const partner = others[0];
      if (partner === undefined || together(await ops.read.bridges(), at.target, partner)) return 0;
      ops.bridge(at.target, partner);
      return 1;
    },
    redo: async () => {
      if (bridgedWith(await ops.read.bridges(), at.target).length === 0) return 0;
      ops.unbridge(at.target);
      return 1;
    },
  };
}
