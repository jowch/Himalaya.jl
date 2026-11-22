import type { Exposure } from "../api";

type Frame = Pick<Exposure, "id" | "status">;

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
  return keptExposures(exposures).reduce<T | undefined>(
    (best, e) => (best === undefined || e.id > best.id ? e : best), undefined);
}
