import type { Peak } from "../../src/api";

type Spec = number | { q: number; excluded?: boolean };

const specQ = (s: Spec): number => (typeof s === "number" ? s : s.q);

export function peakList(exposureId: number, sources: { auto?: Spec[]; manual?: Spec[] }): Peak[] {
  const mint = (source: Peak["source"], specs: Spec[]): Peak[] =>
    specs.map((s, i) => ({
      id: i + 1,
      exposure_id: exposureId,
      q: specQ(s),
      intensity: null,
      prominence: null,
      sharpness: null,
      source,
      excluded: typeof s === "number" ? false : Boolean(s.excluded),
    }));
  return [...mint("auto", sources.auto ?? []), ...mint("manual", sources.manual ?? [])].sort((a, b) => a.q - b.q);
}
