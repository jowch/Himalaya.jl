import { selectedFrames } from "../lib/exposure";
import type { AssignmentState, Exposure, IndexEntry, Peak, SampleBridge, SampleCallMember, SampleLabel, SampleLabelState, LabelColor } from "../api";
import type { KeepCall } from "../lib/bridgeCall";
import type { EditOps } from "../queries";
import type { SetupPatch } from "../lib/queue/mutators/setup";
import type { ExperimentPatch } from "../lib/queue/mutators/experiment";
import { bridgedWith } from "../lib/bridges";
import { comparesOf } from "../lib/compares";
import { labelKey, sameLabel, sharedAcrossBridge } from "../lib/labels";
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
  at: Where & { exposureId: number; owner?: number; prev: Exposure["status"]; next: Exposure["status"]; selected?: boolean },
): UndoEntry {
  const owner = at.owner ?? at.sampleId;
  const set = async (from: Exposure["status"], to: Exposure["status"], selected?: boolean): Promise<number> => {
    const hit = (await ops.read.exposures(owner)).find((e) => e.id === at.exposureId);
    if (!hit || (hit.status ?? null) !== (from ?? null)) return 0;
    if (selected === undefined) ops.setStatus(owner, at.exposureId, to);
    else ops.setStatus(owner, at.exposureId, to, selected);
    return 1;
  };
  return {
    ...at,
    undo: () => set(at.next, at.prev, at.next === "rejected" ? at.selected : undefined),
    redo: () => set(at.prev, at.next, at.prev === "rejected" ? at.selected : undefined),
  };
}

export function selectionEntry(
  ops: EditOps,
  at: Where & { exposureId: number | undefined; owner?: number; prev: readonly number[]; next: readonly number[] },
): UndoEntry {
  const owner = at.owner ?? at.sampleId;
  const key = (ids: readonly number[]): string => [...ids].sort((a, b) => a - b).join(",");
  const set = async (from: readonly number[], to: readonly number[]): Promise<number> => {
    const now = selectedFrames(await ops.read.exposures(owner)).map((e) => e.id);
    if (key(now) !== key(from)) return 0;
    ops.setSelection(owner, [...to]);
    return 1;
  };
  return {
    ...at,
    undo: () => set(at.next, at.prev),
    redo: () => set(at.prev, at.next),
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

async function stillTicked(ops: EditOps, id: PhaseIdentity): Promise<boolean> {
  const [indices, assignment] = await Promise.all([ops.read.indices(id.exposureId), ops.read.assignment(id.exposureId)]);
  const members = new Set(assignment.members);
  return matchIndex(indices, id.phase, id.basis, (ix) => ix.kind === id.kind && members.has(ix.id)) !== undefined;
}

export async function droppedByEdit(ops: EditOps, before: readonly PhaseIdentity[]): Promise<PhaseIdentity[]> {
  if (before.length === 0) return [];
  try {
    await ops.whenSaved();
    const out: PhaseIdentity[] = [];
    for (const id of before) if (!(await stillTicked(ops, id))) out.push(id);
    return out;
  } catch {
    return [];
  }
}

async function thenRetick(ops: EditOps, sent: number, dropped: Promise<readonly PhaseIdentity[]> | undefined): Promise<number> {
  if (sent === 0 || !dropped) return sent;
  const lost = await dropped;
  if (lost.length === 0) return sent;
  await ops.whenSaved();
  return sent + (await setEach(ops, lost, true));
}

type PeakEdit = Where & { exposureId: number; q: number; dropped?: Promise<readonly PhaseIdentity[]> };

export function peakEntry(ops: EditOps, at: PeakEdit & { added: boolean }): UndoEntry {
  return {
    ...at,
    undo: async () => thenRetick(ops, await setManualPeak(ops, at.exposureId, at.q, !at.added), at.dropped),
    redo: () => setManualPeak(ops, at.exposureId, at.q, at.added),
  };
}

export function excludeEntry(ops: EditOps, at: PeakEdit & { excluded: boolean }): UndoEntry {
  const set = async (excluded: boolean): Promise<number> => {
    const hit = near(await ops.read.peaks(at.exposureId), at.q, true);
    if (!hit || hit.excluded === excluded) return 0;
    ops.setPeakExcluded(at.exposureId, hit.id, excluded);
    return 1;
  };
  return { ...at, undo: async () => thenRetick(ops, await set(!at.excluded), at.dropped), redo: () => set(at.excluded) };
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
  const set = async (from: LabelColor | null, to: LabelColor | null): Promise<number> => {
    const stored = (await ops.read.labelColors()).find((e) => sameLabel(e.label, at.tag))?.color ?? null;
    if (stored !== from) return 0;
    ops.setLabelColor(at.tag, to);
    return 1;
  };
  return {
    ...at,
    undo: () => set(at.next, at.prev),
    redo: () => set(at.prev, at.next),
  };
}

const experimentField = (p: { name?: string | null; description?: string | null }, field: keyof ExperimentPatch): string | null => {
  const v = p[field];
  return v === undefined || v === null || v === "" ? null : v;
};

export function experimentEntry(ops: EditOps, at: Where & { prev: ExperimentPatch; next: ExperimentPatch }): UndoEntry {
  const apply = async (from: ExperimentPatch, to: ExperimentPatch): Promise<number> => {
    const server = await ops.read.experiment();
    const patch: ExperimentPatch = {};
    if (to.name !== undefined && experimentField(server, "name") === experimentField(from, "name")) patch.name = to.name;
    if (to.description !== undefined && experimentField(server, "description") === experimentField(from, "description")) {
      patch.description = to.description;
    }
    const sent = Object.keys(patch).length;
    if (sent > 0) ops.patchExperiment(patch);
    return sent;
  };
  return {
    ...at,
    undo: () => apply(at.next, at.prev),
    redo: () => apply(at.prev, at.next),
  };
}

const setupName = (v: string | null | undefined, folder: string): string | null => {
  const t = (v ?? "").trim();
  return t === "" || t === folder ? null : t;
};

export function setupEntry(ops: EditOps, at: Where & { setupId: number; prev: SetupPatch; next: SetupPatch }): UndoEntry {
  const apply = async (from: SetupPatch, to: SetupPatch): Promise<number> => {
    const server = (await ops.read.experiment()).setups?.find((s) => s.id === at.setupId);
    if (!server) return 0;
    const patch: SetupPatch = {};
    if (to.name !== undefined && setupName(server.display_name, server.folder) === setupName(from.name, server.folder)) patch.name = to.name;
    if (to.included !== undefined && server.included === from.included) patch.included = to.included;
    const sent = Object.keys(patch).length;
    if (sent > 0) ops.patchSetup(at.setupId, patch);
    return sent;
  };
  return {
    ...at,
    undo: () => apply(at.next, at.prev),
    redo: () => apply(at.prev, at.next),
  };
}

function together(list: readonly SampleBridge[], a: number, b: number): boolean {
  return bridgedWith(list, a).includes(b);
}

function bridgeGains(members: readonly number[], labels: readonly SampleLabel[]): Map<number, Map<string, SampleLabelState | undefined>> {
  const mine = (sid: number): SampleLabel[] => labels.filter((l) => l.sample_id === sid && sharedAcrossBridge(l.label));
  const union = new Map<string, string>();
  for (const sid of members) for (const l of mine(sid)) if (l.state === "confirmed" && !union.has(labelKey(l.label))) union.set(labelKey(l.label), l.label);
  return new Map(members.map((sid) => {
    const had = new Map(mine(sid).map((l) => [labelKey(l.label), l.state] as const));
    return [sid, new Map([...union.keys()].filter((k) => had.get(k) !== "confirmed").map((k) => [k, had.get(k)] as const))];
  }));
}

export function bridgeEntry(
  ops: EditOps,
  at: Where & {
    from: number;
    to: number;
    bridges: readonly SampleBridge[];
    labels: readonly SampleLabel[];
    keep: KeepCall | undefined;
    states: { from: AssignmentState; to: AssignmentState };
    lost: readonly SampleCallMember[];
  },
): UndoEntry {
  const sideA = [at.from, ...bridgedWith(at.bridges, at.from)];
  const sideB = [at.to, ...bridgedWith(at.bridges, at.to)];
  const split = sideA.length <= sideB.length ? sideA : sideB;
  const gains = bridgeGains([...sideA, ...sideB], at.labels);
  return {
    ...at,
    undo: async () => {
      if (!together(await ops.read.bridges(), at.from, at.to)) return 0;
      const now = await ops.read.labels();
      const call = await ops.read.call(at.from);
      for (const s of split) ops.unbridge(s);
      for (const s of split.slice(1)) ops.bridge(split[0]!, s, "both");
      let sent = split.length * 2 - 1;
      for (const [side, prev] of [[sideA, at.states.from], [sideB, at.states.to]] as const) {
        if (prev === call.state || (prev === "null" && call.members.length > 0)) continue;
        ops.setCallState(side[0]!, prev);
        sent += 1;
      }
      for (const m of at.lost) {
        if (call.members.some((x) => x.exposure_id === m.exposure_id && x.index_id === m.index_id)) continue;
        ops.addPhase(m.exposure_id, m.index_id);
        sent += 1;
      }
      const undone = new Set<string>();
      for (const side of [sideA, sideB]) {
        for (const sid of side) {
          for (const [key, prev] of gains.get(sid) ?? []) {
            const hit = now.find((l) => l.sample_id === sid && labelKey(l.label) === key && l.state === "confirmed");
            if (!hit || undone.has(`${side[0]}:${key}`)) continue;
            undone.add(`${side[0]}:${key}`);
            if (prev === "suggested") ops.restoreLabelSuggestion(sid, hit.label);
            else ops.removeLabel(sid, hit.label);
            sent += 1;
          }
        }
      }
      return sent;
    },
    redo: async () => {
      if (together(await ops.read.bridges(), at.from, at.to)) return 0;
      ops.bridge(at.from, at.to, at.keep);
      return 1;
    },
  };
}

export function settleCallEntry(
  ops: EditOps,
  at: Where & { from: number; to: number; keep: number; lost: readonly SampleCallMember[] },
): UndoEntry {
  const has = (members: readonly SampleCallMember[], m: SampleCallMember): boolean =>
    members.some((x) => x.exposure_id === m.exposure_id && x.index_id === m.index_id);
  return {
    ...at,
    undo: async () => {
      const call = await ops.read.call(at.from);
      let sent = 0;
      for (const m of at.lost) {
        if (has(call.members, m)) continue;
        ops.addPhase(m.exposure_id, m.index_id);
        sent += 1;
      }
      return sent;
    },
    redo: async () => {
      ops.bridge(at.from, at.to, at.keep);
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
      ops.bridge(at.target, partner, "both");
      return 1;
    },
    redo: async () => {
      if (bridgedWith(await ops.read.bridges(), at.target).length === 0) return 0;
      ops.unbridge(at.target);
      return 1;
    },
  };
}

export type CompareEdit =
  | { kind: "add"; ref: number }
  | { kind: "remove"; ref: number; subtracted: boolean }
  | { kind: "subtract"; from: number | undefined; to: number | undefined };

export function compareEntry(ops: EditOps, at: Where & { edit: CompareEdit }): UndoEntry {
  const now = async () => comparesOf(await ops.read.compares(), at.sampleId);
  const e = at.edit;
  const send = (change: Parameters<EditOps["patchCompares"]>[1]): number => {
    ops.patchCompares(at.sampleId, change);
    return 1;
  };
  if (e.kind === "subtract") {
    const flip = async (expect: number | undefined, next: number | undefined): Promise<number> =>
      (await now()).subtraction === expect ? send({ subtraction: next ?? null }) : 0;
    return { ...at, undo: () => flip(e.to, e.from), redo: () => flip(e.from, e.to) };
  }
  const put = async (): Promise<number> => {
    const s = await now();
    if (s.refs.includes(e.ref)) return 0;
    return send(e.kind === "remove" && e.subtracted && s.subtraction === undefined ? { add: [e.ref], subtraction: e.ref } : { add: [e.ref] });
  };
  const take = async (): Promise<number> => ((await now()).refs.includes(e.ref) ? send({ remove: [e.ref] }) : 0);
  return e.kind === "add" ? { ...at, undo: take, redo: put } : { ...at, undo: put, redo: take };
}
