import type { CorpusSample, Exposure } from "../api";

export interface FrameGroup {
  sampleId: number;
  title: string | undefined;
  exposures: readonly Exposure[];
}

export interface FrameSource {
  sample: CorpusSample;
  frames: readonly Exposure[];
}

export function frameGroups(
  own: FrameSource,
  mates: readonly FrameSource[],
  rank: (s: CorpusSample) => number,
  title: (s: CorpusSample) => string,
): FrameGroup[] {
  const others = mates
    .filter((m) => m.frames.length > 0)
    .map((m, i) => ({ m, i }))
    .sort((a, b) => rank(a.m.sample) - rank(b.m.sample) || a.i - b.i)
    .map(({ m }) => m);
  const named = others.length > 0;
  return [own, ...others].map((src) => ({
    sampleId: src.sample.id,
    title: named ? title(src.sample) : undefined,
    exposures: src.frames,
  }));
}

export function allFrames(groups: readonly FrameGroup[]): Exposure[] {
  return groups.flatMap((g) => g.exposures);
}

export function frameOf(
  groups: readonly FrameGroup[],
  exposureId: number | undefined,
): { group: FrameGroup; index: number } | undefined {
  for (const group of groups) {
    const index = group.exposures.findIndex((e) => e.id === exposureId);
    if (index >= 0) return { group, index };
  }
  return undefined;
}
