import type { Exposure } from "../api";

type Frame = Pick<Exposure, "id" | "status">;
type Kinded = Frame & Pick<Exposure, "kind">;
type Selectable = Kinded & Pick<Exposure, "selected">;
type Sourced = Selectable & { sources?: readonly ({ source_exposure_id: number } | number)[] };

export const MIN_COMBINED = 2;

export function isKept(exposure: Pick<Exposure, "status">): boolean {
  return exposure.status !== "rejected";
}

export function isCombined(exposure: Pick<Exposure, "kind">): boolean {
  return exposure.kind === "combined";
}

export function isSelected(exposure: Selectable): boolean {
  return !isCombined(exposure) && isKept(exposure) && exposure.selected !== false;
}

export function keptExposures<T extends Pick<Exposure, "status">>(
  exposures: readonly T[] | undefined,
): T[] {
  return (exposures ?? []).filter(isKept);
}

export function imageFrames<T extends Pick<Exposure, "kind">>(exposures: readonly T[] | undefined): T[] {
  return (exposures ?? []).filter((e) => !isCombined(e));
}

export function selectedFrames<T extends Selectable>(exposures: readonly T[] | undefined): T[] {
  return (exposures ?? []).filter(isSelected).sort((a, b) => a.id - b.id);
}

export function newestKeptFrame<T extends Kinded>(exposures: readonly T[] | undefined): T | undefined {
  return keptExposures(imageFrames(exposures)).reduce<T | undefined>(
    (best, e) => (best === undefined || e.id > best.id ? e : best), undefined);
}

export function sourceIds(exposure: Sourced): number[] {
  return (exposure.sources ?? [])
    .map((s) => (typeof s === "number" ? s : s.source_exposure_id))
    .sort((a, b) => a - b);
}

export interface Analysis<T> {
  selected: T[];
  analyzed: T | undefined;
}

export function analysisOf<T extends Sourced>(exposures: readonly T[] | undefined): Analysis<T> {
  const selected = selectedFrames(exposures);
  if (selected.length < MIN_COMBINED) return { selected, analyzed: selected[0] };
  const key = selected.map((e) => e.id).join(",");
  const analyzed = (exposures ?? []).find((e) => isCombined(e) && sourceIds(e).join(",") === key);
  return { selected, analyzed };
}

export function analyzedExposure<T extends Sourced>(exposures: readonly T[] | undefined): T | undefined {
  return analysisOf(exposures).analyzed;
}
