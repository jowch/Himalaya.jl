import type { Peak } from "../api";

export type PeakRef = Pick<Peak, "id" | "source">;

export const peakKey = (p: PeakRef): string => `${p.source}:${p.id}`;

export const samePeak = (a: PeakRef, b: PeakRef): boolean => a.id === b.id && a.source === b.source;
