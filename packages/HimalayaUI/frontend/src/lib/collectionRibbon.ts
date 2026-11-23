import type { CollectionTimeline } from "../api";

export const JOIN_GAP_S = 600;

export interface RibbonSpan { from: number; to: number; reviewed: boolean }
export interface RibbonGap { from: number; to: number }

export type RibbonPoint =
  | { kind: "sample"; from: number; to: number; reviewed: boolean; beamOff: number }
  | { kind: "gap"; from: number; to: number };

export function ribbonLayout(tl: CollectionTimeline): { spans: RibbonSpan[]; gaps: RibbonGap[]; end: number } {
  const spans: RibbonSpan[] = [];
  const gaps: RibbonGap[] = [];
  let end = 0;
  tl.samples.forEach(([from, to, reviewed], i) => {
    end = Math.max(end, to);
    const next = tl.samples[i + 1];
    const joined = next !== undefined && next[0] - end <= JOIN_GAP_S ? Math.max(to, next[0]) : to;
    if (next !== undefined && next[0] - end > JOIN_GAP_S) gaps.push({ from: end, to: next[0] });
    spans.push({ from, to: joined, reviewed: reviewed === 1 });
  });
  return { spans, gaps, end };
}

export function ribbonAt(tl: CollectionTimeline, t: number): RibbonPoint | null {
  const { spans, gaps } = ribbonLayout(tl);
  const i = spans.findIndex((s) => t >= s.from && t <= s.to);
  if (i >= 0) {
    const [from, to, reviewed] = tl.samples[i]!;
    const beamOff = tl.beam_off.filter((b) => b >= from && b <= to).length;
    return { kind: "sample", from, to, reviewed: reviewed === 1, beamOff };
  }
  const gap = gaps.find((g) => t > g.from && t < g.to);
  return gap ? { kind: "gap", ...gap } : null;
}

const DAYS = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"];
const MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
const pad = (n: number): string => String(n).padStart(2, "0");

function wallClock(start: string, offsetS: number): Date {
  const [d, t] = start.split("T");
  const [y, mo, da] = (d ?? "").split("-").map(Number);
  const [h, mi, s] = (t ?? "").split(":").map(Number);
  return new Date(Date.UTC(y!, mo! - 1, da!, h!, mi!, s!) + offsetS * 1000);
}

export function clockTime(start: string, offsetS: number): string {
  const d = wallClock(start, offsetS);
  return `${pad(d.getUTCHours())}:${pad(d.getUTCMinutes())}`;
}

export function dayAndTime(start: string, offsetS: number): string {
  const d = wallClock(start, offsetS);
  return `${DAYS[d.getUTCDay()]} ${d.getUTCDate()} ${MONTHS[d.getUTCMonth()]} ${clockTime(start, offsetS)}`;
}

export function duration(seconds: number): string {
  const m = Math.round(seconds / 60);
  if (m < 60) return `${m} min`;
  return m % 60 === 0 ? `${m / 60} h` : `${Math.floor(m / 60)} h ${m % 60} min`;
}

export function describeRibbon(tl: CollectionTimeline): { started: string; length: string; left: string } {
  const left = tl.samples.filter(([, , reviewed]) => reviewed !== 1).length;
  return {
    started: dayAndTime(tl.start, 0),
    length: duration(ribbonLayout(tl).end),
    left: left === 0 ? "All reviewed" : `${left} to review`,
  };
}
