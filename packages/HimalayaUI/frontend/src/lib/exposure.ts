import type { Exposure } from "../api";

type Frame = Pick<Exposure, "id" | "selected" | "status">;

export function isKept(exposure: Pick<Exposure, "status">): boolean {
  return exposure.status !== "rejected";
}

export function keptExposures<T extends Pick<Exposure, "status">>(
  exposures: readonly T[] | undefined,
): T[] {
  return (exposures ?? []).filter(isKept);
}

export function representativeExposure<T extends Frame>(
  exposures: readonly T[] | undefined,
): T | undefined {
  const newestFirst = keptExposures(exposures).sort((a, b) => b.id - a.id);
  return newestFirst.find((e) => e.selected) ?? newestFirst[0];
}
