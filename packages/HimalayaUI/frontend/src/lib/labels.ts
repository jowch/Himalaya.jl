import type { SampleLabel, TagColor, TagColorEntry } from "../api";

export const STANDARD_TAGS = ["Buffer", "Calibrant", "Empty capillary"] as const;
export const MAX_TAG_LENGTH = 64;

const fold = (s: string): string => s.trim().toLowerCase();
const labelKey = (s: string): string => s.trim().replace(/[A-Z]+/g, (m) => m.toLowerCase());

export function sameLabel(a: string, b: string): boolean {
  return labelKey(a) === labelKey(b);
}

export function standardTag(label: string): string | undefined {
  return STANDARD_TAGS.find((t) => sameLabel(t, label));
}

export function displayTag(label: string): string {
  return standardTag(label) ?? label;
}

export function cleanTag(raw: string): string | undefined {
  const t = raw.trim();
  if (t.length === 0 || t.length > MAX_TAG_LENGTH || /[\r\n]/.test(t)) return undefined;
  return t;
}

export function labelsFor(list: readonly SampleLabel[], sampleId: number): SampleLabel[] {
  return list.filter((l) => l.sample_id === sampleId);
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

export function tagVocabulary(list: readonly SampleLabel[]): string[] {
  const counts = new Map<string, { label: string; n: number }>();
  for (const l of list) {
    if (l.state !== "confirmed") continue;
    const k = labelKey(l.label);
    const cur = counts.get(k);
    counts.set(k, { label: cur?.label ?? displayTag(l.label), n: (cur?.n ?? 0) + 1 });
  }
  const used = [...counts.values()].sort((a, b) => b.n - a.n || a.label.localeCompare(b.label)).map((c) => c.label);
  return [...used, ...STANDARD_TAGS.filter((t) => !counts.has(labelKey(t)))];
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

function closeness(tag: string, query: string): number | undefined {
  const t = fold(tag);
  if (t.startsWith(query)) return 0;
  if (t.includes(query)) return 1;
  if (query.length < 3) return undefined;
  const allowed = query.length >= 6 ? 2 : 1;
  const head = t.slice(0, query.length);
  return Math.min(editDistance(head, query), editDistance(t, query)) <= allowed ? 2 : undefined;
}

export interface TagOptions {
  matches: string[];
  create: string | undefined;
}

export function tagOptions(vocabulary: readonly string[], query: string, onSample: readonly string[]): TagOptions {
  const taken = (t: string): boolean => onSample.some((s) => sameLabel(s, t));
  const free = vocabulary.filter((t) => !taken(t));
  const q = fold(query);
  if (q.length === 0) return { matches: free, create: undefined };
  const ranked = free
    .map((t, i) => ({ t, i, c: closeness(t, q) }))
    .filter((r): r is { t: string; i: number; c: number } => r.c !== undefined)
    .sort((a, b) => a.c - b.c || a.i - b.i)
    .map((r) => r.t);
  const typed = cleanTag(query);
  const known = typed !== undefined && (vocabulary.some((t) => sameLabel(t, typed)) || taken(typed));
  return { matches: ranked, create: typed !== undefined && !known ? typed : undefined };
}

export function sampleMatchesFilter(name: string, labels: readonly SampleLabel[], filter: string): boolean {
  const f = fold(filter);
  if (f.length === 0) return true;
  return name.toLowerCase().includes(f) || labels.some((l) => l.label.toLowerCase().includes(f));
}

export function applyTagColor(list: readonly TagColorEntry[], label: string, color: TagColor | null): TagColorEntry[] {
  const rest = list.filter((e) => !sameLabel(e.label, label));
  return color === null ? rest : [...rest, { label: label.trim(), color }];
}

export function defaultTagColor(label: string): TagColor {
  return standardTag(label) ? "blue" : "violet";
}

export function tagColor(colors: readonly TagColorEntry[], label: string): TagColor {
  return colors.find((e) => sameLabel(e.label, label))?.color ?? defaultTagColor(label);
}

export function tagColorIsDefault(colors: readonly TagColorEntry[], label: string): boolean {
  return !colors.some((e) => sameLabel(e.label, label));
}
