import type { SampleLabel } from "../api";
import { buffersFirst } from "../lib/labels";
import { naturalSortBy } from "../lib/naturalSort";
import { byRecent } from "../lib/recents";

interface Named {
  id: number;
  name: string;
}

const idKey = (s: Named): string => String(s.id);

export function compareOrder<T extends Named>(samples: readonly T[], labels: readonly SampleLabel[], recent: readonly string[]): T[] {
  return byRecent(buffersFirst(naturalSortBy(samples, (s) => s.name), labels), idKey, recent);
}

export interface BackgroundOrderInput<T extends Named> {
  samples: readonly T[];
  open: number;
  current: number | null;
  runBuffers: readonly number[];
  setupBuffers: readonly number[];
  recent: readonly string[];
}

export function backgroundOrder<T extends Named>(p: BackgroundOrderInput<T>): T[] {
  const byId = new Map(p.samples.map((s) => [s.id, s]));
  const taken = new Set<number>(p.current === null ? [p.open] : [p.open, p.current]);
  const out: T[] = [];
  const add = (id: number): void => {
    const s = byId.get(id);
    if (!s || taken.has(id)) return;
    taken.add(id);
    out.push(s);
  };
  p.runBuffers.forEach(add);
  p.setupBuffers.forEach(add);
  naturalSortBy(p.samples, (s) => s.name).forEach((s) => add(s.id));
  return byRecent(out, idKey, p.recent);
}
