import type { CorpusSample, Exposure, IndexEntry } from "../api";

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
    .filter((ix) => ix.score !== null)
    .sort((a, b) => (b.score ?? 0) - (a.score ?? 0))
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
