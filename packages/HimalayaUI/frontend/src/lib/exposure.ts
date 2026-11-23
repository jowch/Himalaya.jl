import type { Exposure } from "../api";

type Frame = Pick<Exposure, "id" | "status">;
type Kinded = Frame & Pick<Exposure, "kind">;

export const MIN_COMBINED = 2;

export function isKept(exposure: Pick<Exposure, "status">): boolean {
  return exposure.status !== "rejected";
}

export function isCombined(exposure: Pick<Exposure, "kind">): boolean {
  return exposure.kind === "combined";
}

export function isIncluded(exposure: Pick<Exposure, "status" | "included">): boolean {
  return isKept(exposure) && exposure.included !== false;
}

export function keptExposures<T extends Pick<Exposure, "status">>(
  exposures: readonly T[] | undefined,
): T[] {
  return (exposures ?? []).filter(isKept);
}

export function imageFrames<T extends Pick<Exposure, "kind">>(exposures: readonly T[] | undefined): T[] {
  return (exposures ?? []).filter((e) => !isCombined(e));
}

export function newestKeptFrame<T extends Kinded>(exposures: readonly T[] | undefined): T | undefined {
  return keptExposures(imageFrames(exposures)).reduce<T | undefined>(
    (best, e) => (best === undefined || e.id > best.id ? e : best), undefined);
}

export interface Combination<T> {
  combined: T;
  members: T[];
}

export function combinationOf<T extends Kinded & Pick<Exposure, "included">>(
  exposures: readonly T[] | undefined,
): Combination<T> | undefined {
  const combined = (exposures ?? []).find(isCombined);
  if (!combined) return undefined;
  const members = imageFrames(exposures).filter(isIncluded);
  return members.length >= MIN_COMBINED ? { combined, members } : undefined;
}

export function representativeExposure<T extends Kinded & Pick<Exposure, "included">>(
  exposures: readonly T[] | undefined,
): T | undefined {
  return combinationOf(exposures)?.combined ?? newestKeptFrame(exposures);
}
