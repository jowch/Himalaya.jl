import type { CorpusSample, Exposure, IndexEntry } from "../api";
import { samePeak, type PeakRef } from "../lib/peakRef";

export type SampleStatus =
  | { kind: "unreviewed" }
  | { kind: "form_factor" }
  | { kind: "phase"; phase: string };

export const MAX_CANDIDATES = 9;

export function sampleStatus(s: Pick<CorpusSample, "assignment_state" | "phase">): SampleStatus {
  if (s.assignment_state === "form_factor") return { kind: "form_factor" };
  if (s.assignment_state === "indexed" && s.phase) return { kind: "phase", phase: s.phase };
  return { kind: "unreviewed" };
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
    })
    .slice(0, MAX_CANDIDATES);
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

// physics: I₀ is the upstream ion chamber. With the shutter closed its electrometer
// physics: reads a negative offset, about -0.07 to -0.12 V across ~100k real SSRL
// physics: frames, never 0. A beam at low chamber gain reads only ~0.08 V, so the
// physics: magnitude can't separate the two; the sign can. No frame reads between
// physics: -0.07 and +0.08 V, and the signed cutoff sits in that gap with margin.
export const BEAM_OFF_I0_VOLTS = 0.04;

export function beamLikelyOff(frame: { incident_monitor: number | null } | undefined): boolean {
  const i0 = frame?.incident_monitor;
  return i0 != null && i0 < BEAM_OFF_I0_VOLTS;
}
