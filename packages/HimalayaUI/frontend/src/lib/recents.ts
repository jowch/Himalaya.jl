const PREFIX = "himalaya-ui:recent:";
export const MAX_RECENTS = 20;

function storageKey(kind: string): string {
  return PREFIX + kind;
}

export function recentRank(kind: string): string[] {
  try {
    const raw = localStorage.getItem(storageKey(kind));
    const list: unknown = raw ? JSON.parse(raw) : [];
    return Array.isArray(list) ? list.filter((k): k is string => typeof k === "string") : [];
  } catch {
    return [];
  }
}

export function recordRecent(kind: string, key: string): void {
  const next = [key, ...recentRank(kind).filter((k) => k !== key)].slice(0, MAX_RECENTS);
  try {
    localStorage.setItem(storageKey(kind), JSON.stringify(next));
  } catch {
    return;
  }
}

export function byRecent<T>(items: readonly T[], keyOf: (item: T) => string, recent: readonly string[]): T[] {
  const rank = new Map(recent.map((k, i) => [k, i]));
  return items
    .map((item, i) => ({ item, i, r: rank.get(keyOf(item)) }))
    .sort((a, b) => (a.r ?? Infinity) - (b.r ?? Infinity) || a.i - b.i)
    .map((x) => x.item);
}
