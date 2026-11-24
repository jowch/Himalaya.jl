import type { SampleBackground, Trace } from "../api";
import { compareScale, subtractTrace, type CompareScale } from "./compare";
import type { TraceStatus } from "./useCompareTraces";

export type SkipReason = "buffer" | "none" | "no-default" | "no-module" | "hidden" | "loading" | "missing" | "no-trace";

export const SKIP_NOTES: Record<SkipReason, string> = {
  buffer: "Not subtracted: this is a buffer",
  none: "Not subtracted: its background is set to none",
  "no-default": "Not subtracted: its run has no buffer to subtract",
  "no-module": "Not subtracted: only autosampler samples have a background",
  hidden: "Not subtracted: the setup of its background is hidden",
  loading: "Loading the background trace",
  missing: "Not subtracted: no trace for its background",
  "no-trace": "Not subtracted: no trace for this sample",
};

export type MinusPlan =
  | { kind: "drawn"; id: number; background: number; q: number[]; d: number[]; scale: CompareScale }
  | { kind: "skipped"; id: number; reason: SkipReason };

interface PlanInput {
  chips: readonly number[];
  backgrounds: readonly SampleBackground[];
  traces: ReadonlyMap<number, Trace>;
  statuses: ReadonlyMap<number, TraceStatus>;
  hidden: ReadonlySet<number>;
  setupOf: (sampleId: number) => number | null;
  peakQs: readonly number[];
}

type Background = { ok: true; ref: number } | { ok: false; reason: SkipReason };

function backgroundOf(entry: SampleBackground | undefined, hidden: ReadonlySet<number>): Background {
  if (!entry) return { ok: false, reason: "no-module" };
  if (entry.is_buffer) return { ok: false, reason: "buffer" };
  if (entry.ref_sample_id === null) return { ok: false, reason: entry.source === "user" ? "none" : "no-default" };
  if (hidden.has(entry.ref_sample_id)) return { ok: false, reason: "hidden" };
  return { ok: true, ref: entry.ref_sample_id };
}

export function backgroundsToFetch(chips: readonly number[], backgrounds: readonly SampleBackground[], hidden: ReadonlySet<number>): number[] {
  const byId = new Map(backgrounds.map((b) => [b.sample_id, b]));
  const refs = chips.flatMap((id) => {
    const b = backgroundOf(byId.get(id), hidden);
    return b.ok ? [b.ref] : [];
  });
  return [...new Set(refs)];
}

export function planMinusCompares({ chips, backgrounds, traces, statuses, hidden, setupOf, peakQs }: PlanInput): MinusPlan[] {
  const byId = new Map(backgrounds.map((b) => [b.sample_id, b]));
  return chips.map((id): MinusPlan => {
    const background = backgroundOf(byId.get(id), hidden);
    if (!background.ok) return { kind: "skipped", id, reason: background.reason };
    const own = traces.get(id);
    if (!own) return { kind: "skipped", id, reason: statuses.get(id) === "loading" ? "loading" : "no-trace" };
    const ref = traces.get(background.ref);
    if (!ref) return { kind: "skipped", id, reason: statuses.get(background.ref) === "loading" ? "loading" : "missing" };
    const ownSetup = setupOf(id);
    const scale = compareScale(own, ref, peakQs, ownSetup !== null && ownSetup === setupOf(background.ref));
    const { q, d } = subtractTrace(own, ref, scale.scale);
    return { kind: "drawn", id, background: background.ref, q, d, scale };
  });
}
