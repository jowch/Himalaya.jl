export interface PaletteItem {
  id: string;
  label: string;
  name?: string;
  value?: string | undefined;
  hint?: string | undefined;
  shortcut?: string | undefined;
  page?: string;
  form?: () => PaletteForm;
  run?: () => void;
  swatch?: string | undefined;
  danger?: boolean;
  typeOnly?: boolean;
}

export interface FieldOption {
  value: string;
  label: string;
  swatch?: string | undefined;
}

export interface FormField {
  id: string;
  label: string;
  kind: "text" | "choice";
  value: string;
  placeholder?: string;
  options?: readonly FieldOption[];
  check?: (value: string) => string | undefined;
}

export interface PaletteForm {
  title: string;
  verb: string;
  fields: readonly FormField[];
  submit: (values: Readonly<Record<string, string>>) => void;
}

export interface PaletteGroup {
  title: string;
  items: readonly PaletteItem[];
  searchOnly?: boolean;
}

export interface PalettePage {
  title: string;
  arg?: string;
  placeholder?: string;
  items: readonly PaletteItem[];
  deep?: boolean;
  empty?: string;
  extra?: (query: string) => readonly PaletteItem[];
}

export type PaletteRow =
  | { kind: "head"; title: string }
  | { kind: "item"; item: PaletteItem; path?: string };

export const DEEP_TITLE = "Inside this sample's pages";
export const DEEP_LIMIT = 6;
export const SEARCH_LIMIT = 20;

export function matches(text: string, query: string): boolean {
  const hay = text.toLowerCase();
  return query.toLowerCase().split(/\s+/).filter(Boolean).every((w) => hay.includes(w));
}

export function rootRows(groups: readonly PaletteGroup[], pages: Readonly<Record<string, PalettePage>>, query: string): PaletteRow[] {
  const q = query.trim();
  const out: PaletteRow[] = [];
  const pushGroup = (title: string, rows: PaletteRow[]): void => {
    if (rows.length === 0) return;
    out.push({ kind: "head", title }, ...rows);
  };
  const searchGroup = (g: PaletteGroup): void => {
    const hits = g.items.filter((i) => matches(i.label, q));
    const exact = hits.filter((i) => i.name !== undefined && i.name.toLowerCase() === q.toLowerCase());
    const ordered = [...exact, ...hits.filter((i) => !exact.includes(i))].slice(0, SEARCH_LIMIT);
    pushGroup(g.title, ordered.map((item) => ({ kind: "item", item })));
  };
  const exactGroups = q ? groups.filter((g) => g.searchOnly && g.items.some((i) => i.name !== undefined && i.name.toLowerCase() === q.toLowerCase())) : [];
  for (const g of exactGroups) searchGroup(g);
  for (const g of groups) {
    if (g.searchOnly) continue;
    const items = g.items.filter((i) => (q ? matches(i.label, q) : !i.typeOnly));
    pushGroup(g.title, items.map((item) => ({ kind: "item", item })));
  }
  if (!q) return out;
  const deep: PaletteRow[] = [];
  for (const page of Object.values(pages)) {
    if (!page.deep) continue;
    for (const item of page.items) if (matches(`${page.title} ${item.label}`, q)) deep.push({ kind: "item", item, path: page.title });
  }
  pushGroup(DEEP_TITLE, deep.slice(0, DEEP_LIMIT));
  for (const g of groups) if (g.searchOnly && !exactGroups.includes(g)) searchGroup(g);
  return out;
}

export function pageRows(page: PalettePage, query: string): PaletteRow[] {
  const q = query.trim();
  const items = [...page.items.filter((i) => !q || matches(i.label, q)), ...(page.extra?.(q) ?? [])];
  return items.map((item) => ({ kind: "item", item }));
}

export function fieldErrors(fields: readonly FormField[], values: Readonly<Record<string, string>>): Record<string, string> {
  const out: Record<string, string> = {};
  for (const f of fields) {
    const why = f.check?.(values[f.id] ?? "");
    if (why) out[f.id] = why;
  }
  return out;
}

export function cycleOption(options: readonly FieldOption[], value: string, dir: 1 | -1): string {
  if (options.length === 0) return value;
  const i = options.findIndex((o) => o.value === value);
  if (i < 0) return options[dir === 1 ? 0 : options.length - 1]!.value;
  return options[(i + dir + options.length) % options.length]!.value;
}

export function pickable(rows: readonly PaletteRow[]): { item: PaletteItem; path?: string }[] {
  return rows.flatMap((r) => (r.kind === "item" ? [r.path === undefined ? { item: r.item } : { item: r.item, path: r.path }] : []));
}

export interface PaletteSection {
  title: string | undefined;
  rows: { item: PaletteItem; path?: string; index: number }[];
}

export function sections(rows: readonly PaletteRow[]): PaletteSection[] {
  const out: PaletteSection[] = [];
  let index = 0;
  for (const r of rows) {
    if (r.kind === "head") {
      out.push({ title: r.title, rows: [] });
      continue;
    }
    if (out.length === 0) out.push({ title: undefined, rows: [] });
    out[out.length - 1]!.rows.push(r.path === undefined ? { item: r.item, index: index++ } : { item: r.item, path: r.path, index: index++ });
  }
  return out;
}
