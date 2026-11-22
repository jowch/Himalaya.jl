import { LABEL_COLORS, type SampleLabel, type LabelColor, type LabelColorEntry } from "../api";

export const STANDARD_LABELS = ["Buffer", "Calibrant", "Empty capillary"] as const;
export const MAX_LABEL_LENGTH = 64;

const fold = (s: string): string => s.trim().toLowerCase();
export const labelKey = (s: string): string => s.trim().replace(/[A-Z]+/g, (m) => m.toLowerCase());

const DISTANCE_LABEL = /^\s*\d+(\.\d+)?\s*m\s*$/;

export function sharedAcrossBridge(label: string): boolean {
  return !DISTANCE_LABEL.test(label);
}

export function sameLabel(a: string, b: string): boolean {
  return labelKey(a) === labelKey(b);
}

export function standardLabel(label: string): (typeof STANDARD_LABELS)[number] | undefined {
  return STANDARD_LABELS.find((t) => sameLabel(t, label));
}

export function displayLabel(label: string): string {
  return standardLabel(label) ?? label;
}

export function cleanLabel(raw: string): string | undefined {
  const t = raw.trim();
  if (t.length === 0 || t.length > MAX_LABEL_LENGTH || /[\r\n]/.test(t)) return undefined;
  return t;
}

export function labelsFor(list: readonly SampleLabel[], sampleId: number): SampleLabel[] {
  return list.filter((l) => l.sample_id === sampleId);
}

export function bufferRole(list: readonly SampleLabel[], sampleId: number): SampleLabel["state"] | undefined {
  const hit = list.filter((l) => l.sample_id === sampleId && sameLabel(l.label, "Buffer"));
  if (hit.some((l) => l.state === "confirmed")) return "confirmed";
  return hit.length > 0 ? "suggested" : undefined;
}

export function buffersFirst<T extends { id: number }>(samples: readonly T[], list: readonly SampleLabel[]): T[] {
  const rank = (id: number): number => {
    const role = bufferRole(list, id);
    return role === "confirmed" ? 0 : role === "suggested" ? 1 : 2;
  };
  return samples.map((s) => ({ s, r: rank(s.id) })).sort((a, b) => a.r - b.r).map(({ s }) => s);
}

export function applyAddLabel(list: readonly SampleLabel[], sampleId: number, label: string): SampleLabel[] {
  const hit = list.find((l) => l.sample_id === sampleId && sameLabel(l.label, label));
  if (hit) return list.map((l) => (l === hit ? { ...l, state: "confirmed" } : l));
  const fresh = Math.min(0, ...list.map((l) => l.id)) - 1;
  return [...list, { id: fresh, sample_id: sampleId, label: label.trim(), state: "confirmed", source: "user" }];
}

export function applyRemoveLabel(list: readonly SampleLabel[], sampleId: number, label: string): SampleLabel[] {
  return list.filter((l) => !(l.sample_id === sampleId && sameLabel(l.label, label)));
}

export function applyRestoreLabelSuggestion(list: readonly SampleLabel[], sampleId: number, label: string): SampleLabel[] {
  const hit = list.find((l) => l.sample_id === sampleId && sameLabel(l.label, label));
  if (hit) return list.map((l) => (l === hit ? { ...l, state: "suggested", source: "ingest" } : l));
  const fresh = Math.min(0, ...list.map((l) => l.id)) - 1;
  return [...list, { id: fresh, sample_id: sampleId, label: label.trim(), state: "suggested", source: "ingest" }];
}

export function labelVocabulary(list: readonly SampleLabel[]): string[] {
  const counts = new Map<string, { label: string; n: number }>();
  for (const l of list) {
    if (l.state !== "confirmed") continue;
    const k = labelKey(l.label);
    const cur = counts.get(k);
    counts.set(k, { label: cur?.label ?? displayLabel(l.label), n: (cur?.n ?? 0) + 1 });
  }
  const used = [...counts.values()].sort((a, b) => b.n - a.n || a.label.localeCompare(b.label)).map((c) => c.label);
  return [...used, ...STANDARD_LABELS.filter((t) => !counts.has(labelKey(t)))];
}

function editDistance(a: string, b: string): number {
  let prev = Array.from({ length: b.length + 1 }, (_, j) => j);
  for (let i = 1; i <= a.length; i++) {
    const cur = [i];
    for (let j = 1; j <= b.length; j++) {
      cur[j] = Math.min(prev[j]! + 1, cur[j - 1]! + 1, prev[j - 1]! + (a[i - 1] === b[j - 1] ? 0 : 1));
    }
    prev = cur;
  }
  return prev[b.length]!;
}

function closeness(label: string, query: string): number | undefined {
  const t = fold(label);
  if (t.startsWith(query)) return 0;
  if (t.includes(query)) return 1;
  if (query.length < 3) return undefined;
  const allowed = query.length >= 6 ? 2 : 1;
  const head = t.slice(0, query.length);
  return Math.min(editDistance(head, query), editDistance(t, query)) <= allowed ? 2 : undefined;
}

export interface LabelOptions {
  matches: string[];
  create: string | undefined;
}

export function labelOptions(vocabulary: readonly string[], query: string, onSample: readonly string[]): LabelOptions {
  const taken = (t: string): boolean => onSample.some((s) => sameLabel(s, t));
  const free = vocabulary.filter((t) => !taken(t));
  const q = fold(query);
  if (q.length === 0) return { matches: free, create: undefined };
  const ranked = free
    .map((t, i) => ({ t, i, c: closeness(t, q) }))
    .filter((r): r is { t: string; i: number; c: number } => r.c !== undefined)
    .sort((a, b) => a.c - b.c || a.i - b.i)
    .map((r) => r.t);
  const typed = cleanLabel(query);
  const known = typed !== undefined && (vocabulary.some((t) => sameLabel(t, typed)) || taken(typed));
  return { matches: ranked, create: typed !== undefined && !known ? typed : undefined };
}

export function sampleMatchesFilter(name: string, labels: readonly SampleLabel[], filter: string): boolean {
  const f = fold(filter);
  if (f.length === 0) return true;
  return name.toLowerCase().includes(f) || labels.some((l) => l.label.toLowerCase().includes(f));
}

export function applyLabelColor(list: readonly LabelColorEntry[], label: string, color: LabelColor | null): LabelColorEntry[] {
  const rest = list.filter((e) => !sameLabel(e.label, label));
  return color === null ? rest : [...rest, { label: label.trim(), color }];
}

export const STANDARD_LABEL_COLORS: Readonly<Record<(typeof STANDARD_LABELS)[number], LabelColor>> = {
  Buffer: "blue",
  Calibrant: "ochre",
  "Empty capillary": "stone",
};

const SYSTEM_COLORS: readonly LabelColor[] = Object.values(STANDARD_LABEL_COLORS);
export const USER_LABEL_COLORS: readonly LabelColor[] = LABEL_COLORS.filter((c) => !SYSTEM_COLORS.includes(c));

export function labelColorIsFixed(label: string): boolean {
  return standardLabel(label) !== undefined;
}

export function defaultLabelColor(label: string): LabelColor {
  const standard = standardLabel(label);
  return standard ? STANDARD_LABEL_COLORS[standard] : "violet";
}

export function labelColor(colors: readonly LabelColorEntry[], label: string): LabelColor {
  if (labelColorIsFixed(label)) return defaultLabelColor(label);
  const stored = colors.find((e) => sameLabel(e.label, label))?.color;
  return stored !== undefined && USER_LABEL_COLORS.includes(stored) ? stored : defaultLabelColor(label);
}

export function labelColorIsDefault(colors: readonly LabelColorEntry[], label: string): boolean {
  return labelColorIsFixed(label) || !colors.some((e) => sameLabel(e.label, label) && USER_LABEL_COLORS.includes(e.color));
}

export function firstUnusedColor(colors: readonly LabelColorEntry[], labels: readonly string[]): LabelColor {
  const uses = new Map<LabelColor, number>(USER_LABEL_COLORS.map((c) => [c, 0]));
  for (const l of labels) {
    if (labelColorIsFixed(l)) continue;
    const c = labelColor(colors, l);
    uses.set(c, (uses.get(c) ?? 0) + 1);
  }
  const least = Math.min(...uses.values());
  return USER_LABEL_COLORS.find((c) => uses.get(c) === least)!;
}
