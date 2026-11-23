import type { CorpusSample, Exposure, IndexEntry } from "../api";
import { samePeak, type PeakRef } from "../lib/peakRef";

export type SampleStatus =
  | { kind: "unreviewed" }
  | { kind: "form_factor" }
  | { kind: "phase"; phases: string[] };

export function sampleStatus(s: Pick<CorpusSample, "assignment_state" | "phase" | "call_phases">): SampleStatus {
  if (s.assignment_state === "form_factor") return { kind: "form_factor" };
  if (s.assignment_state !== "indexed" || !s.phase) return { kind: "unreviewed" };
  return { kind: "phase", phases: [s.phase, ...(s.call_phases ?? []).filter((p) => p !== s.phase)] };
}

export function firstUnreviewed<T extends Pick<CorpusSample, "assignment_state" | "phase">>(samples: readonly T[]): T | undefined {
  return samples.find((s) => sampleStatus(s).kind === "unreviewed") ?? samples[0];
}

export function nextUnreviewed(
  samples: readonly CorpusSample[],
  currentId: number | undefined,
  dir: 1 | -1,
): number | undefined {
  const n = samples.length;
  if (n === 0) return undefined;
  const start = Math.max(0, samples.findIndex((s) => s.id === currentId));
  for (let k = 1; k <= n; k++) {
    const s = samples[(((start + dir * k) % n) + n) % n]!;
    if (s.id !== currentId && sampleStatus(s).kind === "unreviewed") return s.id;
  }
  return undefined;
}

export function rankCandidates(indices: readonly IndexEntry[] | undefined): IndexEntry[] {
  return (indices ?? [])
    .filter((ix) => ix.kind === "speculative" || ix.score !== null)
    .sort((a, b) => {
      const custom = Number(b.kind === "speculative") - Number(a.kind === "speculative");
      if (custom !== 0) return custom;
      return a.kind === "speculative" ? b.id - a.id : (b.score ?? 0) - (a.score ?? 0);
    });
}

export function stepExposure(
  exposures: readonly Pick<Exposure, "id">[],
  currentId: number | undefined,
  dir: 1 | -1,
): number | undefined {
  const i = exposures.findIndex((e) => e.id === currentId);
  if (i < 0) return exposures[0]?.id;
  return exposures[Math.min(exposures.length - 1, Math.max(0, i + dir))]!.id;
}

export function nearestKeptExposure(
  exposures: readonly Pick<Exposure, "id" | "status">[],
  droppedId: number,
): number | undefined {
  const at = exposures.findIndex((e) => e.id === droppedId);
  for (let d = 1; d < exposures.length; d++) {
    for (const i of [at + d, at - d]) {
      const e = exposures[i];
      if (e && e.status !== "rejected") return e.id;
    }
  }
  return undefined;
}

export function imageUrl(e: Pick<Exposure, "id" | "image_path" | "image_version">, thumb: boolean): string | null {
  if (e.image_path === null) return null;
  const params = new URLSearchParams();
  if (thumb) params.set("thumb", "1");
  if (e.image_version) params.set("v", e.image_version);
  const qs = params.toString();
  return `/api/exposures/${e.id}/image${qs ? `?${qs}` : ""}`;
}

export type PeakClick = { kind: "exclude"; excluded: boolean } | { kind: "remove" };

export function peakClick(
  peaks: readonly (PeakRef & { excluded: boolean })[],
  ref: PeakRef,
): PeakClick | null {
  const p = peaks.find((pk) => samePeak(pk, ref));
  if (!p) return null;
  return p.source === "auto" ? { kind: "exclude", excluded: !p.excluded } : { kind: "remove" };
}

// physics: a _tot.dat holds raw detector counts per pixel. A beam-off frame records
// physics: only the dark floor, so its median is ~0, while even a dim real frame keeps
// physics: ~10 counts. The detector decides; the ion chamber can't, because a frame
// physics: with a negative I₀ can still hold real scattering.
export const DARK_MEDIAN_COUNTS = 1;

export function traceIsDark(trace: { I: readonly number[] } | undefined): boolean {
  if (!trace || trace.I.length === 0) return false;
  const sorted = [...trace.I].sort((a, b) => a - b);
  return sorted[Math.floor(sorted.length / 2)]! < DARK_MEDIAN_COUNTS;
}

// physics: with the shutter closed the ion chamber's electrometer reads a negative
// physics: offset, so I₀ at or below zero supports a beam-off reading but never
// physics: decides it alone.
export function beamOffWarning(e: Exposure): string | undefined {
  return e.beam_off && e.status !== "rejected" ? "Beam likely off" : undefined;
}

export function i0BelowZero(frame: { incident_monitor: number | null } | undefined): boolean {
  const i0 = frame?.incident_monitor;
  return i0 != null && i0 <= 0;
}
